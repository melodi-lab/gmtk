/*
 * MMapMatrix.h
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
#include <map>

#include "error.h"
#include "GMTK_TrainingSchedule.h"
#include "Matrix.h"


#define MAX_TEMP_FILENAME_LENGTH 1024

// A subclass of Matrix (see Matrix.h) that is backed by a file.
// You can only do O(1) space accesses via GetCols, Submatrix, and the like.
// Trying to access the entire matrix (e.g., via Start or Vec) will fail.
// You can do O(1) space writes via PutCols.

class MMapMatrix : public Matrix {
  friend class MutableVector;

  bool    mapped;        // is this instance mmapped?
  char   *fileName;      // backing file name
  int     fd;            // backing file descriptor
  long    pagesize;      // VM system page size
  ssize_t map_length;    // size of mmapped buffer
  static  map<string, unsigned> ref_count;  // Reference counting garbage collector.
                                            // The assignment operator causes the MMapMatrix
                                            // instances to share the same backing file since
                                            // the Matrix class is supposed to be like a double*.
                                            // Thus we can't unmap the file until the last user
                                            // is destroyed.

  static  unsigned fileNumber;              // counter for generating temporary file names


  // Generate backing file names. POSIX has a number of functions that do this,
  // but they all warn of issues in the man pages or during compilation. I don't
  // think the issues are very applicable to our use case, but it's not hard to
  // roll our own...
  char *tempname(char const *dir, char const *dummy) {
    char  tempname_buf[MAX_TEMP_FILENAME_LENGTH];
    char const *tempdir = getenv("GMTKTMPDIR");
    if (!tempdir) tempdir = dir;
    assert(tempdir);
    size_t len = strlen(tempdir);
    if (len > MAX_TEMP_FILENAME_LENGTH - 32) {
      error("ERROR: TMPDIR '%s' is too long, must be less than %u characters\n", tempdir, MAX_TEMP_FILENAME_LENGTH-32);
    }
    sprintf(tempname_buf, "%s/gmtk_%08X.%08X", tempdir, fileNumber++,(unsigned)getpid());
    return strdup(tempname_buf);
  }

public:

  static char const *dmlpTempDir; // temp dir from command line

  MMapMatrix() : Matrix(), mapped(false), fileName(NULL) { }

  MMapMatrix(const Matrix & mat) : Matrix(mat), mapped(false), fileName(NULL) { }

#if 0
  MMapMatrix(double *start, int numR, int numC, int ld, bool trans=false) 
    : Matrix(start, numR, numC, ld, trans), mapped(false), fileName(NULL) 
  { }
#endif

  MMapMatrix(int numR, int numC, int ld, bool trans=false) 
    : Matrix(NULL, numR, numC, ld, trans), mapped(false) 
  {
    assert(numC > 0 && ld > 0 && numR > 0);
    // create backing file
    fileName = tempname(dmlpTempDir, "gmtk_XXXXX");
    if (!fileName) {
      error("ERROR: unable to generate temporary file name\n");
    }
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
    map_length = ( ((off_t)numC * (off_t)ld * sizeof(double) - 1) / (off_t)pagesize + 1 ) * (off_t)pagesize;
    // write to end of file to make it map_length bytes long so we can
    // memory map that many bytes
    if (lseek(fd, (off_t)map_length, SEEK_SET) == (off_t) -1) {
      perror(fileName);
      error("ERROR: unable to seek to offset %u in temporary file '%s'\n", map_length, fileName);
    }
    if (write(fd, &numR, sizeof(numR)) != sizeof(numR)) {
      perror(fileName);
      error("ERROR: error writing to temporary file '%s'\n", fileName);
    }
    _start = (const double *)mmap(NULL, map_length, PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0);
    if (!_start || _start == (void *) -1) {
      perror(fileName);
      error("ERROR: unable to allocate %ld bytes via mmap()\n", map_length);
    }
    ref_count[string(fileName)] = 1;
    mapped = true;
  }


  // Copy ctor - create a new backing file 
  MMapMatrix(MMapMatrix const &that) 
    : Matrix(NULL, that._numR, that._numC, that._ld, that._trans), 
      pagesize(that.pagesize), map_length(that.map_length)
  { 
    if (that.mapped) {
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
	result = write(fd, buf, map_length-bwritten);
	if (result < 0) {
	  perror(fileName);
	  error("ERROR: error duplicating temporary file '%s' to '%s'\n", that.fileName, fileName);
	}
	bwritten += result;
      } while (bwritten < map_length);
      _start = (const double *)mmap(NULL, map_length, PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0);
      if (!_start || _start == (void *) -1) {
	perror(fileName);
	error("ERROR: unable to allocate %ld bytes via mmap()\n", map_length);
      }
      string fn(fileName);
      ref_count[fn] = 1;
      mapped = true;
    } else {
      // just a shallow copy in this case
      _start = that._start;
      mapped = false;
      fileName = NULL;
    }
  }


  MMapMatrix & operator= (const MMapMatrix &rhs) {
    _numR  = rhs._numR;
    _numC  = rhs._numC;
    _ld    = rhs._ld;
    _trans = rhs._trans;

    if (mapped && _start) {
      // free old this
      string fn(fileName);
      if (ref_count[fn] == 1) {
	// this was the last user of the backing file
	munmap((void *)_start, map_length);
	close(fd);
	unlink(fileName);
	ref_count.erase(fn);
	fileName = NULL;
	mapped = false;
      } else {
	// other instances are still using the backing file
	ref_count[fn] = ref_count[fn] - 1;
      }
      if (fileName) free(fileName);
    }
    if (rhs.mapped) {
      // alias the rhs file
      map_length = rhs.map_length;
      pagesize = rhs.pagesize;
      fileName = strdup(rhs.fileName);
      if (!fileName) {
	error("ERROR: unable to duplicate temporary file name\n");
      }
      fd = rhs.fd;
      string fn(fileName);
      ref_count[fn] = ref_count[fn] + 1;
      mapped = true;
    } else {
      // no backing file to alias ...
      mapped = false;
      fileName = NULL;
    }
    _start = rhs._start;
    return *this;
  }


  ~MMapMatrix() {
    if (mapped && _start) {
      string fn(fileName);
      if (ref_count[fn] == 1) {
	// this was the last user of the backing file
	munmap((void *)_start, map_length);
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


  // delete all MMapMatrix backing files
  static void GarbageCollect() {
    for (map<string,unsigned>::iterator it = ref_count.begin(); it != ref_count.end(); ++it) {
      unlink(it->first.c_str());
    }
  }


  char *GetFileName() {
    return mapped ? fileName : NULL;
  }


  // Copy Matrix m into this, starting at column destCol
  void PutCols(Matrix const &m, int destCol) {
    assert(_start && _start != (void *)-1);
    int numCols = m.NumC();
    assert(destCol + numCols <= _numC);
    double *dest = const_cast<double *>(_start) + (ssize_t)destCol * (ssize_t)_ld;
    double const *src  = m.Start();
    double const *end = m.End();
    int numToCopy = m.NumR(), mLd = m.Ld();
    do {
      memcpy((void *)dest, (void *)src, numToCopy * sizeof(double));
      dest += _ld;
      src  += mLd;
    } while (src != end);
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

  
  // non-constant space accessors not supported
  const double *Start() const { assert(0); return NULL; }
  const double *End() const { assert(0); return NULL; }
  Vector Vec() const { assert(0); return Vector(); }
  Vector GetRow(int r) const { assert(0); return Vector(); }
  Matrix GetRows(int beginRow, int endRow) const { assert(0); return Matrix(); }

}; 
