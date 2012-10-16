
/*
 * GMTK_FIRFilter.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2012, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */

#ifndef GMTK_FIRFILTER_H
#define GMTK_FIRFILTER_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "machine-dependent.h"
#include "GMTK_SubmatrixDescriptor.h"
#include "GMTK_Filter.h"

// Apply a Finite Impulse Response filter to the continuous features

class FIRFilter: public Filter {

  unsigned  order; // order N requires frames [t-N,t] to produce output for frame t
  unsigned  numFeatures;
  float    *B;     // filter coefficient matrix
  float    *c;     // constant term

  // For a single (scalar) feature x, an order N FIR is 
  // y[t] = \sum_{j=0}^N b_i x[t-j]
  // This implementation allows a constant term:
  // y[t] = \sum_{j=0}^N b_i x[t-j] + c
  // and treats B and X as a matrices, and and y and c as vectors 
  // y[t] = B' X[t,t-N] + c

  // B is N+1 rows by numFeatures columns stored linearly as <row 0> <row 1> ... <row N>
  //   where b_{ij} is the coefficient for x_j [t-i]  ie, the jth continuous feature at
  //   time lag i

  Data32  *buffer;
  unsigned buffSize; // in Data32's

 public:
  
  FIRFilter() {buffer=NULL; buffSize=0;}

  // Read B and c from fileName
  FIRFilter(char const *fileName, Filter *nextFilter);

  // B or c may optionally be NULL to do y[t] = B x[t,t-N]'  or  y[t] = t[t] + c
  FIRFilter(unsigned order, unsigned numFeatures, float *B, float *c, Filter *nextFilter=NULL) 
    : Filter(nextFilter), order(order), numFeatures(numFeatures), B(B), c(c), buffer(NULL), buffSize(0)
  {}

  ~FIRFilter() {
    if (B) delete[] B;
    if (c) delete[] c;
    if (buffer) free(buffer);
  }

  subMatrixDescriptor *
  getRequiredInput(unsigned first, unsigned count, 
		   unsigned inputContinous, unsigned inputDiscrete,
		   unsigned inputTotalFrames);


  subMatrixDescriptor describeLocalOutput(subMatrixDescriptor const &inputDescription);
  
  virtual void getNextFrameInfo(unsigned &numNewIn, unsigned &dropOldIn, unsigned &numNewOut,
				unsigned inputContinuous, unsigned inputDiscrete,
				subMatrixDescriptor &input);

  Data32 const *localTransform(Data32 const *inputSubMatrix,
			       subMatrixDescriptor const &inputDescription,
			       subMatrixDescriptor *outputDescription=NULL); 
};

#endif
