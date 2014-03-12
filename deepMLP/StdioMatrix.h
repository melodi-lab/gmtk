/*
 * StdioMatrix.h - a Matrix class backed by a stdio file
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (C) 2014 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

#pragma once

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif
#if HAVE_FSEEKO
#  define FM_SEEK fseeko
#else
#  define FM_SEEK fseek
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <vector>

#include "error.h"
#include "Matrix.h"
#include "FileBackedMatrix.h"

// The copy ctor copys the copyee's backing file in chunks of this size (# of floats / chunk)
#define FILE_COPY_BUFFER_SIZE ((1<<20)/sizeof(float))

class StdioMatrix : public FileBackedMatrix {
  friend class MutableVector;

  FILE   *f;                               // backing file descriptor

  static float buf[FILE_COPY_BUFFER_SIZE]; // for copying files in copy ctor
  static double at;                        // To hold ::At() result

  mutable vector<double> doubleBuf; // The backing file stores floats, so we need
  mutable vector<float>  floatBuf;  // space to convert them to doubles. And vice versa.

public:

  StdioMatrix() : FileBackedMatrix(), doubleBuf(0), floatBuf(0) { }

  // Just alias mat, not file-backed...
  StdioMatrix(const Matrix & mat) : FileBackedMatrix(mat), doubleBuf(0), floatBuf(0) { }

  // Create a new file-backed matrix. Put in values with PutCols
  StdioMatrix(int numR, int numC, int ld, bool trans=false) 
    : FileBackedMatrix(numR, numC, ld, trans), doubleBuf(0), floatBuf(0)
  {
    // O_EXCL so it fails if we didn't create the file
    int fd = open(fileName, O_RDWR|O_CREAT|O_EXCL, S_IRWXU);
    if (fd == -1) {
      perror(fileName);
      error("ERROR: unable to create temporary file '%s'\n", fileName);
    }
#if HAVE_FDOPEN && !defined(__CYGWIN__)
    f = fdopen(fd, "r+b");
    if (!f) {
      perror(fileName);
      error("ERROR: unable to open temporary file '%s'\n", fileName);
    }
#else
    // Race condition -- file may disappear between close() and fopen()
    if (close(fd)) {
      perror(fileName);
      error("ERROR: unable to close temporary file '%s'\n", fileName);
    }
    f = fopen(fileName, "r+b");
    if (!f) {
      perror(fileName);
      error("ERROR: unable to open temporary file '%s'\n", fileName);
    }
#endif
    file_length = (off_t)numC * (off_t)ld * sizeof(float);
    ref_count[string(fileName)] = 1;
    backed = true;
  }


  // deep copy ctor - create a new backing file and copy that's contents into it
  StdioMatrix(StdioMatrix const &that) 
    : FileBackedMatrix(NULL, that._numR, that._numC, that._ld, that._trans), 
      doubleBuf(0), floatBuf(0)
  { 
    file_length = that.file_length;
    if (that.backed) {
      fileName = tempname(dmlpTempDir, "gmtk_XXXXX");
      if (!fileName) {
	error("ERROR: unable to generate temporary file name\n");
      }
      int fd = open(fileName, O_RDWR|O_CREAT|O_EXCL, S_IRWXU);
      if (fd == -1) {
	perror(fileName);
	error("ERROR: unable to create temporary file '%s'\n", fileName);
      }
#if HAVE_FDOPEN && !defined(__CYGWIN__)
      f = fdopen(fd, "bw+");
      if (!f) {
	perror(fileName);
	error("ERROR: unable to create temporary file '%s'\n", fileName);
      }
#else
      // Race condition -- file may disappear between close() and fopen()
      if (close(fd)) {
	perror(fileName);
	error("ERROR: unable to close temporary file '%s'\n", fileName);
      }
      f = fopen(fileName, "r+b");
      if (!f) {
	perror(fileName);
	error("ERROR: unable to open temporary file '%s'\n", fileName);
      }
#endif
      size_t bwritten = 0, read_result, write_result;
      assert(sizeof(size_t) == sizeof(off_t));
      if (FM_SEEK(that.f, 0, SEEK_SET)) {
	perror(that.fileName);
	error("ERROR: seek failed in temporary file '%s'\n", that.fileName);
      }

      do { // copy that's backing file
	read_result = fread(buf, sizeof(float), FILE_COPY_BUFFER_SIZE, that.f);
	if (ferror(that.f)) {
	  perror(that.fileName);
	  error("ERROR: error duplicating temporary file '%s' to '%s'\n", that.fileName, fileName);
	}
	write_result = fwrite(buf, sizeof(float), read_result, f);
	if (write_result != read_result) {
	  perror(fileName);
	  error("ERROR: error duplicating temporary file '%s' to '%s'\n", that.fileName, fileName);
	}
	bwritten += write_result;
      } while (bwritten < (size_t)file_length);
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
  StdioMatrix & operator= (const StdioMatrix &rhs) {
    _numR  = rhs._numR;
    _numC  = rhs._numC;
    _ld    = rhs._ld;
    _trans = rhs._trans;

    if (backed) {
      // free old this
      string fn(fileName);
      if (ref_count[fn] == 1) {
	// this was the last user of the backing file
	if (fclose(f)) {
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
      fileName = strdup(rhs.fileName);
      if (!fileName) {
	error("ERROR: unable to duplicate temporary file name\n");
      }
      f = rhs.f;
      string fn(fileName);
      ref_count[fn] = ref_count[fn] + 1;
      backed = true;
    } else {
      // no backing file to alias ...
      backed = false;
      fileName = NULL;
    }
    return *this;
  }


  ~StdioMatrix() {
    if (backed) {
      string fn(fileName);
      if (ref_count[fn] == 1) {
	// this was the last user of the backing file
	if (fclose(f)) {
	  perror(fileName);
	  error("ERROR: closing temporary file '%s'\n", fileName);
	}
	if (unlink(fileName)) {
	  perror(fileName);
	  error("ERROR: deleting temporary file '%s'\n", fileName);
	}
	ref_count.erase(fn);
      } else {
	// other instances are still using the backing file
	ref_count[fn] = ref_count[fn] - 1;
      }
    }
    if (fileName) free(fileName);
    fileName = NULL;
  }

#if 0
  // this is now in base class FileBackedMatrix

  // Copy Matrix m into this, starting at column destCol
  void PutCols(Matrix const &m, int destCol) {
    PutCols(const_cast<double *>(m.Start()), m.NumC(), m.NumR(), m.Ld(), destCol);
  }
#endif

  // Copy Matrix(src, numRows, numCols, srcLd) into this, starting at column destCol
  void PutCols(double *src, int numCols, int numRows, int srcLd, int destCol) {
    assert(backed);
    assert(destCol + numCols <= _numC);
    assert(numRows == _ld);
    if ((unsigned) numRows > floatBuf.size()) {
      floatBuf.resize(numRows);
    }
    // seek to destCol
    off_t dest = (off_t)destCol * (off_t)_ld * sizeof(float);
    if (FM_SEEK(f, dest, SEEK_SET)) {
      perror(fileName);
      error("ERROR: failed to seek in temporary file '%s'\n", fileName);
    }
    // write input matrix column by column (converting to float)
    double *end = src + numCols * srcLd;
    do {
      for (int i=0; i < numRows; i+=1)
	floatBuf[i] = (float)src[i];
      if (fwrite(&floatBuf[0], sizeof(float), numRows, f) != (size_t)numRows) {
	perror(fileName);
	error("ERROR: writing to temporary file '%s'\n", fileName);
      }
      src  += srcLd;
    } while (src != end);
    fflush(f);
  }
  
  static double dAt; // gotta convert float to double
  virtual const double & At(int r, int c) const {
    if (!backed) return Matrix::At(r,c);
    if (_trans) std::swap(r, c);
    if (r < 0) r += _numR;
    if (c < 0) c += _numC;
    assert (r >= 0 && r < _numR);
    assert (c >= 0 && c < _numC);

    off_t dest = ((off_t)c * (off_t)_ld + (off_t) r) * sizeof(float);
    if (FM_SEEK(f, dest, SEEK_SET)) {
      perror(fileName);
      error("ERROR: failed to seek in temporary file '%s'\n", fileName);
    }
    float element;
    if (fread(&element, sizeof(float), 1, f) != 1) {
      perror(fileName);
      error("ERROR: failed to read element from temporary file '%s'\n", fileName);
    }
    dAt = (double) element;
    return dAt;
  }

  // seek to requested column, convert it to double
  Vector GetCol(int c) {
    if (!backed) return Matrix::GetCol(c);
    if (_trans) return Trans().GetRow(c);  // fails if file backed since not O(1) space

    if (c < 0) c += _numC;
    assert (0 <= c && c < _numC);

    if (floatBuf.size() < (size_t)_numR) floatBuf.resize((size_t)_numR);
    off_t dest = (off_t)c * (off_t)_ld * sizeof(float);
    if (FM_SEEK(f, dest, SEEK_SET)) {
      perror(fileName);
      error("ERROR: failed to seek in temporary file '%s'\n", fileName);
    }
    float *fp = const_cast<float *>(&floatBuf[0]);
    if (fread(fp, sizeof(float), _numR, f) != (size_t)_numR) {
      perror(fileName);
      error("ERROR: StdioMatrix::GetCol() failed to read from temporary file '%s'\n", fileName);
    }
    if (doubleBuf.size() < (size_t)_numR) doubleBuf.resize(_numR);
    for (int i=0; i < _numR; i+=1) doubleBuf[i] = (double)floatBuf[i];
    return Vector(&doubleBuf[0], _numR);
  }

  // seek to requested submatrix, convert to double
  virtual Matrix SubMatrix(int beginRow, int endRow, int beginCol, int endCol) const {
    if (!backed) return Matrix::SubMatrix(beginRow, endRow, beginCol, endCol);

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

    size_t numRows = endRow - beginRow;
    size_t numCols = endCol - beginCol;
    size_t numInputs = numCols * (size_t)_ld;
    size_t numOutputs= numCols * numRows;
    if (floatBuf.size()  < numInputs) floatBuf.resize(numInputs);
    if (doubleBuf.size() < numOutputs) doubleBuf.resize(numOutputs);

    off_t dest = ((off_t)beginCol * (off_t)_ld + (off_t)beginRow) * sizeof(float);
    if (FM_SEEK(f, dest, SEEK_SET)) {
      perror(fileName);
      error("ERROR: failed to seek in temporary file '%s'\n", fileName);
    }
    float *fp = &floatBuf[0];
    size_t nread;
    if ((nread=fread(fp, sizeof(float), numInputs, f)) != numInputs) {
      perror(fileName);
      error("ERROR: StdioMatrix::SubMatrix() read %u of %u from temporary file '%s'\n", nread, numInputs, fileName);
    }
    unsigned i = 0;
    for (size_t c=0; c < numCols; c+=1) {
      for (size_t r=0; r < numRows; r+=1) {
	doubleBuf[i] = (double) (fp[r]);
	i += 1;
      }
      fp += _ld;
    }
    return Matrix(&doubleBuf[0], numRows, numCols, numRows, _trans);
  }

}; 
  
