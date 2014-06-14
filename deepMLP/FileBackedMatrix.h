/*
 * FileBackedMatrix.h - Base class for file-backed matrices
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
#include "Matrix.h"


#define MAX_TEMP_FILENAME_LENGTH 1024

// A subclass of Matrix (see Matrix.h) that is backed by a file.
// You can only do O(1) space accesses via GetCols, Submatrix, and the like.
// Trying to access the entire matrix (e.g., via Start or Vec) will fail.
// You can do O(1) space writes via PutCols.

// The StdioMatrix subclass uses stdio operations to access the files, and 
// so requires "large file support" to handle large temporary files. It also 
// stores the data as floats rather than doubles to save space.

// The MMapMatrix subclass memory maps the double array, and thus requires
// 64-bit addresses to support large temporary files.

class FileBackedMatrix : public Matrix {
  friend class MutableVector;

 protected:

  bool    backed;        // is this instance actually backed by a file?
  char   *fileName;      // backing file name
  ssize_t file_length;    // size of backing file (in bytes)

  static  map<string, unsigned> ref_count;  // Reference counting garbage collector.
                                            // The assignment operator causes the FileBackedMatrix
                                            // instances to share the same backing file since
                                            // the Matrix class is supposed to act like a double*.
                                            // Thus we can't delete the file until the last user
                                            // is destroyed.

  static  unsigned fileNumber;              // counter for generating temporary file names


  // Generate backing file names. POSIX has a number of functions that do this,
  // but they all warn of issues in the man pages or during compilation. I don't
  // think the issues are very applicable to our use case, but it's not hard to
  // roll our own...
  char *tempname(char const *tempdir, char const *dummy) {
    char  tempname_buf[MAX_TEMP_FILENAME_LENGTH];
    if (!tempdir) {
      char const *envdir = getenv("GMTKTMPDIR");
      if (!envdir) {
	error("ERROR: you must specify a directory for temporary files using either the -tempDir command line option or the GMTKTMPDIR environment variable\n");
      }
      tempdir = envdir;
    }
    assert(tempdir);
    size_t len = strlen(tempdir);
    if (len > MAX_TEMP_FILENAME_LENGTH - 32) {
      error("ERROR: temporary file directory path '%s' is too long, must be less than %u characters\n", tempdir, MAX_TEMP_FILENAME_LENGTH-32);
    }
    sprintf(tempname_buf, "%s/gmtk_%08X.%08X", tempdir, fileNumber++,(unsigned)getpid());
    return strdup(tempname_buf);
  }

 public:

  static char const *dmlpTempDir; // temp dir from command line

  FileBackedMatrix() : Matrix(), backed(false), fileName(NULL) { }

  // Just alias mat, not file-backed...
  FileBackedMatrix(const Matrix & mat) : Matrix(mat), backed(false), fileName(NULL) { }


  // Just alias Matrix(start, numR, numC, ld, trans), not file-backed...
  FileBackedMatrix(double *start, int numR, int numC, int ld, bool trans=false) 
    : Matrix(start, numR, numC, ld, trans), backed(false) // subclass ctor will set to true
  { }

  // Create a new file-backed matrix. Put in values with PutCols
  FileBackedMatrix(int numR, int numC, int ld, bool trans=false) 
    : Matrix(NULL, numR, numC, ld, trans), backed(false) 
  {
    assert(numC > 0 && ld > 0 && numR > 0);
    // create backing file
    fileName = tempname(dmlpTempDir, "gmtk_XXXXX");
    if (!fileName) {
      error("ERROR: unable to generate temporary file name\n");
    }
  }

  // Delete all backing files - the subclass = operators create alias (shallow copies),
  // and the instances might not have all be destroyed when we're finished using them...
  static void GarbageCollect() {
    for (auto it = ref_count.begin(); it != ref_count.end(); ++it) {
      if (unlink(it->first.c_str())) {
	perror(it->first.c_str());
	error("ERROR: deleting temporary file '%s'\n", it->first.c_str());
      }
    }
  }

  char *GetFileName() {
    return backed ? fileName : NULL;
  }


  // Copy Matrix m into this, starting at column destCol
  virtual void PutColsM(Matrix const &m, int destCol) {
    PutCols(const_cast<double *>(m.Start()), m.NumC(), m.NumR(), m.Ld(), destCol);
  }

  // Copy Matrix(src, numRows, numCols, srcLd) into this, starting at column destCol
  virtual void PutCols(double *src, int numCols, int numRows, int srcLd, int destCol) = 0;
  
  // non-constant space accessors not supported
  virtual const double *Start() const { assert(0); return NULL; }
  virtual const double *End() const { assert(0); return NULL; }
  virtual Vector Vec() const { assert(0); return Vector(); }
  virtual Vector GetRow(int r) const { assert(0); return Vector(); }
  virtual Matrix GetRows(int beginRow, int endRow) const { assert(0); return Matrix(); }

}; 
