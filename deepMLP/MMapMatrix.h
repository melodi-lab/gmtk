/*
 * MMapMatrix.h - a Matrix class that mmap()s a file
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (C) 2013 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

#pragma once

#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>

#include <assert.h>
#include <string.h>

#include "error.h"
#include "Matrix.h"
#include "FileBackedMatrix.h"

class MMapMatrix : public FileBackedMatrix {
  friend class MutableVector;

  int     fd;            // backing file descriptor
  long    pagesize;      // VM system page size

public:

  MMapMatrix() : FileBackedMatrix() { }

  // Just alias mat, not file-backed...
  MMapMatrix(const Matrix & mat) : FileBackedMatrix(mat) { }

  // Create a new file-backed matrix. Put in values with PutCols
  MMapMatrix(int numR, int numC, int ld, bool trans=false) 
    : FileBackedMatrix(numR, numC, ld, trans)
  {
    // O_EXCL so it fails if we didn't create the file
    fd = open(fileName, O_RDWR|O_CREAT|O_EXCL, S_IRWXU);
    if (fd == -1) {
      perror(fileName);
      error("ERROR: unable to create temporary file '%s'\n", fileName);
    }
    // memory mappings must be a multiple of the VM system page size
#if defined(_SC_PAGESIZE)
    pagesize = sysconf(_SC_PAGESIZE);
#elif defined(PAGE_SIZE)
    pagesize = sysconf(PAGE_SIZE);
#else
    error("ERROR: unknown pagesize\n");
#endif
    file_length = ( ((off_t)numC * (off_t)ld * sizeof(double) - 1) / (off_t)pagesize + 1 ) * (off_t)pagesize;
    // write to end of file to make it file_length bytes long so we can
    // memory map that many bytes
    if (lseek(fd, (off_t)file_length, SEEK_SET) == (off_t) -1) {
      perror(fileName);
      error("ERROR: unable to seek to offset %u in temporary file '%s'\n", file_length, fileName);
    }
    if (write(fd, &numR, sizeof(numR)) != sizeof(numR)) {
      perror(fileName);
      error("ERROR: error writing to temporary file '%s'\n", fileName);
    }
    _start = (const double *)mmap(NULL, file_length, PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0);
    if (!_start || _start == (void *) -1) {
      perror(fileName);
      error("ERROR: unable to allocate %ld bytes via mmap()\n", file_length);
    }
    ref_count[string(fileName)] = 1;
    backed = true;
  }


  // deep copy ctor - create a new backing file and copy that's contents into it 
  MMapMatrix(MMapMatrix const &that) 
    : FileBackedMatrix(NULL, that._numR, that._numC, that._ld, that._trans), pagesize(that.pagesize)
  { 
    file_length = that.file_length;
    if (that.backed) {
      fileName = tempname(dmlpTempDir, "gmtk_XXXXX");
      if (!fileName) {
	error("ERROR: unable to generate temporary file name\n");
      }
      fd = open(fileName, O_RDWR|O_CREAT|O_EXCL, S_IRWXU);
      if (fd == -1) {
	perror(fileName);
	error("ERROR: unable to create temporary file '%s'\n", fileName);
      }
      // appearantly there's an upper limit on how much you can write in one call...
      ssize_t bwritten = 0, result;
      do {
	char *buf = ((char *)that._start)+bwritten;
	result = write(fd, buf, file_length-bwritten);
	if (result < 0) {
	  perror(fileName);
	  error("ERROR: error duplicating temporary file '%s' to '%s'\n", that.fileName, fileName);
	}
	bwritten += result;
      } while (bwritten < file_length);
      _start = (const double *)mmap(NULL, file_length, PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0);
      if (!_start || _start == (void *) -1) {
	perror(fileName);
	error("ERROR: unable to allocate %ld bytes via mmap()\n", file_length);
      }
      string fn(fileName);
      ref_count[fn] = 1;
      backed = true;
    } else {
      // just a shallow copy in this case
      _start = that._start;
      backed = false;
      fileName = NULL;
    }
  }


  // shallow copy - just alias the backing file
  MMapMatrix & operator= (const MMapMatrix &rhs) {
    _numR  = rhs._numR;
    _numC  = rhs._numC;
    _ld    = rhs._ld;
    _trans = rhs._trans;

    if (backed && _start) {
      // free old this
      string fn(fileName);
      if (ref_count[fn] == 1) {
	// this was the last user of the backing file
	munmap((void *)_start, file_length);
	if (close(fd)) {
	  perror(fileName);
	  error("ERROR: closing temporary file '%s'\n", fileName);
	}
	if (unlink(fileName)) {
	  perror(fileName);
	  error("ERROR: deleting temporary file '%s'\n", fileName);
	}
	ref_count.erase(fn);
	backed = false;
      } else {
	// other instances are still using the backing file
	ref_count[fn] = ref_count[fn] - 1;
      }
      if (fileName) free(fileName);
      fileName = NULL;
    }
    if (rhs.backed) {
      // alias the rhs file
      file_length = rhs.file_length;
      pagesize = rhs.pagesize;
      fileName = strdup(rhs.fileName);
      if (!fileName) {
	error("ERROR: unable to duplicate temporary file name\n");
      }
      fd = rhs.fd;
      string fn(fileName);
      ref_count[fn] = ref_count[fn] + 1;
      backed = true;
    } else {
      // no backing file to alias ...
      backed = false;
      fileName = NULL;
    }
    _start = rhs._start;
    return *this;
  }


  ~MMapMatrix() {
    if (backed && _start) {
      string fn(fileName);
      if (ref_count[fn] == 1) {
	// this was the last user of the backing file
	munmap((void *)_start, file_length);
	close(fd);
	unlink(fileName);
	ref_count.erase(fn);
      } else {
	// other instances are still using the backing file
	ref_count[fn] = ref_count[fn] - 1;
      }
    }
    if (fileName) free(fileName);
  }


  // Copy Matrix(src, numRows, numCols, srcLd) into this, starting at column destCol
  void PutCols(double *src, int numCols, int numRows, int srcLd, int destCol) {
    assert(_start && _start != (void *)-1);
    assert(destCol + numCols <= _numC);
    double *dest = const_cast<double *>(_start) + (ssize_t)destCol * (ssize_t)_ld;
    double *end = src + numCols * srcLd;
    do {
      memcpy((void *)dest, (void *)src, numRows * sizeof(double));
      dest += _ld;
      src += srcLd;
    } while (src != end);
  }

}; 
