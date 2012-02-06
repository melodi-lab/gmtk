
/*
 * GMTK_Filter.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2012, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */

#ifndef GMTK_FILTER_H
#define GMTK_FILTER_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "machine-dependent.h"

// Packages up the "shape" of a section of an observation
// matrix.
class subMatrixDescriptor {
 public:
  unsigned firstFrame;            // frame sub-matrix starts at
  unsigned numFrames;             // # frames in sub-matrix
  unsigned numContinuous;
  unsigned numDiscrete;
  unsigned fullMatrixFrameCount;  // # frames in full matrix

  subMatrixDescriptor(unsigned firstFrame, unsigned numFrames,
		      unsigned numContinuous, unsigned numDiscrete,
		      unsigned fullMatrixFrameCount)
    : firstFrame(firstFrame), 
      numFrames(numFrames),
      numContinuous(numContinuous),
      numDiscrete(numDiscrete),
      fullMatrixFrameCount(fullMAtrixFrameCount)
  {
  }
} 


//     input submatrix -->  Filter --> output submatrix

class Filter {

  Filter *nextFilter;

 public:
  
  // The filter's client (eg, inference) needs the portion
  // of the filter's output described by requestedSubMatrix.
  // getRequiredInput() returns the portion of the fitler's 
  // input necessary to produce the requested output portion.
  virtual subMatrixDescriptor 
    getRequiredInput(subMatrixDescriptor const &requestedSubMatrix) 
  {
    // transform requested to what I need as input
    subMatrixDescriptor myInputSubMatrix = requestedSubMatrix;
    if (nextFilter) {
      return nextFilter->getRequiredInput(myInputSubMatrix);
    } else {
      return myInputSubMatrix;
    }
  }

  // What will the output look like if the input described by
  // inputDescription is feed into the filter?
  virtual subMatrixDescriptor
    describeOutput(subMatrixDescriptor const &inputDescription) 
  {
    subMatrixDescriptor myInputSubMatrix;
    if (nextFilter) {
      myInputSubmatrix = nextFilter->describeOutput(inputDescription);
    } else {
      myInputSubmatrix = inputDescription;
    }
    // describe my output given myInputSubMatrix
    return myInputSubMatrix;
  }

  // Returns the filter's output given the inputSubMatrix
  // described by inputDescription (presumably created by 
  // getRequiredInput()).  If outputDescription is non-NULL, 
  // it will describe the portion of filter's output returned 
  // (which should hopefully match the requestedSubMatrix passed
  // to getRequiredInput()).
  virtual Data32 const *
    transform(Data32 const *inputSubMatrix,
	      subMatrixDescriptor const &inputDescription,
	      subMatrixDescriptor *outputDescription=NULL) 
  {
    Data32 const *inputData;
    Data32 const *outputData;
    subMatrixDescriptor myInputDescriptor, myOutputDescriptor;

    if (nextFilter) {
      inputData = nextFilter->transform(inputSubMatrix,
					inputDescription,
					&myInputDescriptor);
    } else {
      myInputDescriptor = inputDescription;
      inputData = inputSubMatrix;
    }
    // transform inputData
    outputData = inputData;
    if (outputDescription)
      *outputDescription = inputDescription;
    return outputData;
  }
}

#endif
