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
#include <fcntl.h>

#include <assert.h>
#include <string.h>
#include <map>

#include "error.h"
#include "GMTK_TrainingSchedule.h"
#include "Matrix.h"


#define MAX_TEMP_FILENAME_LENGTH 1024

class MMapMatrix : public Matrix {
  friend class MutableVector;

  bool    mapped;
  char   *fileName;
  int     fd;
  long    pagesize;
  ssize_t map_length;
  static  map<string, unsigned> ref_count;
  static  unsigned fileNumber;


  char *tempname(char const *dir, char const *dummy) {
    char  tempname_buf[MAX_TEMP_FILENAME_LENGTH];
    char const *tempdir = getenv("GMTKTMPDIR");
    if (!tempdir) tempdir = dir;
    assert(tempdir);
    size_t len = strlen(tempdir);
    if (len > MAX_TEMP_FILENAME_LENGTH - 16) {
      error("ERROR: TMPDIR '%s' is too long, must be less than %u characters\n", tempdir, MAX_TEMP_FILENAME_LENGTH-16);
    }
    sprintf(tempname_buf, "%s/gmtk_%08X", tempdir, fileNumber++);
    return strdup(tempname_buf);
  }

public:

  static char const *dmlpTempDir;

  MMapMatrix() : Matrix(), mapped(false), fileName(NULL) { }


  MMapMatrix(const Matrix & mat) : Matrix(mat), mapped(false), fileName(NULL) { }


  MMapMatrix(double *start, int numR, int numC, int ld, bool trans=false) 
    : Matrix(start, numR, numC, ld, trans), mapped(false), fileName(NULL) 
  { }


  MMapMatrix(int numR, int numC, int ld, bool trans=false) 
    : Matrix(NULL, numR, numC, ld, trans), mapped(false) 
  {
    assert(numC > 0 && ld > 0 && numR > 0);
    fileName = tempname(dmlpTempDir, "gmtk_XXXXX");
    if (!fileName) {
      error("ERROR: unable to generate temporary file name\n");
    }
    fd = open(fileName, O_RDWR|O_CREAT|O_EXCL, S_IRWXU);
    if (fd == -1) {
      perror(fileName);
      error("ERROR: unable to create temporary file '%s'\n", fileName);
    }
#if defined(_SC_PAGESIZE)
    pagesize = sysconf(_SC_PAGESIZE);
#elif defined(PAGE_SIZE)
    pagesize = sysconf(PAGE_SIZE);
#else
    error("ERROR: unknown pagesize\n");
#endif
    map_length = ( (numC * ld * sizeof(double) - 1) / pagesize + 1 ) * pagesize;
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
      string fn(fileName);
      if (ref_count[fn] == 1) {
	munmap((void *)_start, map_length);
	close(fd);
	unlink(fileName);
	ref_count.erase(fn);
	fileName = NULL;
	mapped = false;
      } else {
	ref_count[fn] = ref_count[fn] - 1;
      }
      if (fileName) free(fileName);
    }
    if (rhs.mapped) {
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
	munmap((void *)_start, map_length);
	close(fd);
	unlink(fileName);
	ref_count.erase(fn);
      } else {
	ref_count[fn] = ref_count[fn] - 1;
      }
    }
    if (fileName) free(fileName);
  }

  static void GarbageCollect() {
    for (map<string,unsigned>::iterator it = ref_count.begin(); it != ref_count.end(); ++it) {
      unlink(it->first.c_str());
    }
  }


#if 0
void debug() {
  printf("%p %d x %d + %d @ %p mapped %c trans %c   %s\n",
	 this, _numR, _numC, _ld, _start, mapped?'T':'F',_trans?'T':'F', mapped?fileName:"");
#if 0
  int prec=3;
  int len = 4 + prec;
  std::string s = "%-l.pf";
  s[2] = '0' + len;
  s[4] = '0' + prec;
  for (int r = 0; r < 5; ++r) {
    for (int c = 0; c < 5; ++c) {
      printf(s.c_str(), At(r,c));
    }
    printf("\n");
  }
#endif
}
#endif


  char *GetFileName() {
    return mapped ? fileName : NULL;
  }


  void PutCols(Matrix const &m, int destCol) {
    assert(_start && _start != (void *)-1);
    int numCols = m.NumC();
    assert(destCol + numCols <= _numC);
    double *dest = const_cast<double *>(_start) + destCol * _ld;
    double const *src  = m.Start();
    double const *end = m.End();
    int numToCopy = m.NumR(), mLd = m.Ld();
    do {
      memcpy((void *)dest, (void *)src, numToCopy * sizeof(double));
      dest += _ld;
      src  += mLd;
    } while (src != end);
  }


  void PutCols(double *src, int numCols, int numRows, int srcLd, int destCol) {
    assert(_start && _start != (void *)-1);
    assert(destCol + numCols <= _numC);
    double *dest = const_cast<double *>(_start) + destCol * _ld;
    double *end = src + numCols * srcLd;
    do {
      memcpy((void *)dest, (void *)src, numRows * sizeof(double));
      dest += _ld;
      src += srcLd;
    } while (src != end);
  }

  
  const double *Start() const { assert(0); return NULL; }
  const double *End() const { assert(0); return NULL; }

  
  Vector Vec() const { assert(0); return Vector(); }

 
  const double & At(int r, int c) const {
    if (_trans) std::swap(r, c);
    if (r < 0) r += _numR;
    if (c < 0) c += _numC;
    assert (r >= 0 && r < _numR);
    assert (c >= 0 && c < _numC);
    return *(_start + c * _ld + r);
  }


  Vector GetCol(int c) const {
    if (_trans) return Trans().GetRow(c);

    if (c < 0) c += _numC;
    assert (0 <= c && c < _numC);

    return Vector(Start() + c * _ld, _numR, 1, _ld);
  }


  Vector GetRow(int r) const {
    assert(0);
    return Vector();
  }


  Matrix GetCols(int beginCol, int endCol) const {
    return SubMatrix(0, -1, beginCol, endCol);
  }


  Matrix SubMatrix(int beginRow, int endRow, int beginCol, int endCol) const {
#if 0
    assert(0);
    return Matrix();
#else
    if (_trans) {
      std::swap(beginRow, beginCol);
      std::swap(endRow, endCol);
    }

    if (beginRow < 0) beginRow = _numR - ~beginRow;
    if (endRow < 0) endRow = _numR - ~endRow;
    assert (0 <= beginRow && beginRow <= endRow && endRow <= _numR);

    if (beginCol < 0) beginCol = _numC - ~beginCol;
    if (endCol < 0) endCol = _numC - ~endCol;
    assert (0 <= beginCol && beginCol <= endCol && endCol <= _numC);

    return Matrix(_start + beginCol * _ld + beginRow, endRow - beginRow, endCol - beginCol, _ld, _trans);
#endif
  }
}; 
