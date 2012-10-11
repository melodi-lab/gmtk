
/*
 * GMTK_AffineFilter.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2012, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */

#ifndef GMTK_AFFINEFILTER_H
#define GMTK_AFFINEFILTER_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "machine-dependent.h"
#include "GMTK_Filter.h"

// Apply an affine transform to the continuous features
// y = Bx + c

// Actually done with matrix-matrix multiply to process
// multiple frames at once for better performance:
// Y = XB' + [c c ... c]

// If there are n continuous features and m frames, X is 
// stored linearly in memory as
// x_1[t] ... x_n[t] x_1[t+1] ... x_n[t+1] ... x_1[t+m-1] ... x_n[t+m-1]
// as an m x n matrix with stride n. 

// FIXME - would be nice to use stride = numFeatures() to 
// skip over integer features, but that would need a
// float x float or float x double mul_mdmd_md

// FIXME - would be nice to use the PHiPAC A*B+C, which 
// might work with Cstride=0. The problem is that the
// result overwrites C.

// So B should have n rows. If B has k columns, the output of
// the filter will have k continuous features. Store B in the 
// file thusly:

/*

n k
b11 ... b1k
b21 ... b2k
 :       :
bn1 ... bnk
c_1 ... c_k

*/


class AffineFilter: public Filter {

  unsigned  rows;  // must match input # continuous features
  unsigned  cols;  // # output continuous features

  double   *B;     // transform coefficient matrix  rows x cols, stride cols
  double   *c;     // constant term  rows x 1

  Data32  *buffer;   // result buffer
  unsigned buffSize; // in Data32's

  double  *workBuffer;   // holds XB'
  unsigned workBuffSize; // in doubles

  double  *xBuffer;      // holds double version of X
  unsigned xBuffSize;    // in doubles

 public:
  
  AffineFilter() 
    : B(NULL), c(NULL),
      buffer(NULL), buffSize(0),
      workBuffer(NULL), workBuffSize(0),
      xBuffer(NULL), xBuffSize(0)
  {}

  // Read B and c from fileName
  AffineFilter(char const *fileName, Filter *nextFilter);

  // B or c may optionally be NULL to do y = B x  or  y = x + c
  AffineFilter(unsigned rows, unsigned cols, double *B, double *c, Filter *nextFilter=NULL) 
    : Filter(nextFilter),
      rows(rows), cols(cols), 
      B(B), c(c), 
      buffer(NULL), buffSize(0),
      workBuffer(NULL), workBuffSize(0),
      xBuffer(NULL), xBuffSize(0)
  {}

  ~AffineFilter() {
    if (B) delete[] B;
    if (c) delete[] c;
    if (buffer) free(buffer);
    if (workBuffer) free (workBuffer);
    if (xBuffer) free (xBuffer);
  }

  void numOutputFeatures(unsigned inputContinuous, unsigned inputDiscrete,
			 unsigned &outputContinuous, unsigned &outputDiscrete) 
  {
    outputContinuous = cols;
    outputDiscrete = inputDiscrete;
  }

  subMatrixDescriptor *
  getRequiredInput(unsigned first, unsigned count, 
		   unsigned inputContinous, unsigned inputDiscrete,
		   unsigned inputTotalFrames);


  subMatrixDescriptor describeLocalOutput(subMatrixDescriptor const &inputDescription);
  

  Data32 const *localTransform(Data32 const *inputSubMatrix,
			       subMatrixDescriptor const &inputDescription,
			       subMatrixDescriptor *outputDescription=NULL); 
};

#endif
