
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

#include <stdlib.h>
#include <assert.h>

#include "machine-dependent.h"

// Packages up the "shape" of a section of an observation matrix.
class subMatrixDescriptor {
 public:
  unsigned firstFrame;            // frame sub-matrix starts at
  unsigned numFrames;             // # frames in sub-matrix
  unsigned historyFrames;         // # frames just there for history
  unsigned futureFrames;          // # frames just there for future
  unsigned numContinuous;
  unsigned numDiscrete;
  unsigned fullMatrixFrameCount;  // # frames in full matrix
  subMatrixDescriptor *next;

  subMatrixDescriptor() {next = NULL;}

  subMatrixDescriptor(unsigned firstFrame, unsigned numFrames,
		      unsigned historyFrames, unsigned futureFrames,
		      unsigned numContinuous, unsigned numDiscrete,
		      unsigned fullMatrixFrameCount, 
		      subMatrixDescriptor *next = NULL)
    : firstFrame(firstFrame), numFrames(numFrames),
      historyFrames(historyFrames), futureFrames(futureFrames),
      numContinuous(numContinuous), numDiscrete(numDiscrete),
      fullMatrixFrameCount(fullMatrixFrameCount),
      next(next)
  {}

  void deallocate() {
    if (next) {
      next->deallocate();
      delete next;
    }
  }
}; 


// This is the base class for observation filtering.
// It implements the identity filter as an example.

class Filter {

//   input submatrix -->  Filter ( --> nextFilter)*  --> output submatrix
// (ObservationSource)                                     (Inference)
  
 protected:
  Filter *nextFilter;
  
 public:
  
  Filter(Filter *nextFilter = NULL) : nextFilter(nextFilter) {}

  
  virtual ~Filter() { if (nextFilter) delete nextFilter; }

  
  // The filter's client (eg, inference) needs the 
  // [first,first+count)frames of the filter's output.
  // getRequiredInput() returns the portion of the fitler's 
  // input (eg, ObservationFile) necessary to produce the 
  // requested output frames.
  //
  // Note that first and count describe the desired OUTPUT frames.
  // inputContinuous, inputDiscrete and inputTotalFrames are the number
  // of continuous & discrete features and frames in the Filter's INPUT.

  virtual subMatrixDescriptor *
    getRequiredInput(unsigned first, unsigned count, 
		     unsigned inputContinuous, unsigned inputDiscrete,
		     unsigned inputTotalFrames)
  {
    subMatrixDescriptor *nextFilterInput = NULL;
    if (nextFilter) {
      // the client's getting the data from the Filter(s) AFTER me, so
      // ask it what frames I need to provide it as input (it gets as
      // input the # of continuous and discrete features and frames
      // I produce as output).
      unsigned outputContinuous = inputContinuous;  // # floats I produce
      unsigned outputDiscrete   = inputDiscrete;    // # ints I produce
      unsigned outputFrames     = inputTotalFrames; // total # frames I produce
      nextFilterInput = nextFilter->getRequiredInput(first, count, 
						     outputContinuous, outputDiscrete,
						     outputFrames);
      assert(nextFilterInput);
      first = nextFilterInput->firstFrame;
      count = nextFilterInput->numFrames;
    } 
#if 0
    required.firstFrame    = first; // the first frame I need as input
    required.numFrames     = count; // the number of frames I need as input
    required.historyFrames = 0;
    required.futureFrames  = 0;     // # of frames that are just for context
    required.numContinuous = inputContinuous;
    required.numDiscrete   = inputDiscrete;
    required.fullMatrixFrameCount = inputTotalFrames;
#endif
    return new subMatrixDescriptor(first, count, 0, 0, inputContinuous, 
				   inputDiscrete, inputTotalFrames, nextFilterInput);
  }
 

  // What will the output look like if the input described by
  // inputDescription is feed into this filter?
  virtual subMatrixDescriptor
    describeLocalOutput(subMatrixDescriptor const &inputDescription)
  {
    // describe my output given the specified input
    subMatrixDescriptor myOutput = inputDescription;
    return myOutput;
  }


  // What will the output look like if the input described by
  // inputDescription is feed through the entire filter stack?
  subMatrixDescriptor
    describeOutput(subMatrixDescriptor const &inputDescription)
  {
    subMatrixDescriptor myOutput = describeLocalOutput(inputDescription);
    if (nextFilter) {
      // the next filter gets my output as its input
      myOutput = nextFilter->describeOutput(myOutput);
    } 
    return myOutput;
  }
  

  // Returns the filter's output given the inputSubMatrix
  // described by inputDescription (presumably created by 
  // getRequiredInput()).  
  //
  // If outputDescription is non-NULL, 
  // it will describe the portion of filter's output returned 

  virtual Data32 const *
    localTransform(Data32 const *inputSubMatrix,
		   subMatrixDescriptor const &inputDescription,
		   subMatrixDescriptor *outputDescription=NULL) 
  {
    Data32 const *outputData = inputSubMatrix; // perform my transformation
    subMatrixDescriptor myOutput = describeLocalOutput(inputDescription);
    if (outputDescription) *outputDescription = myOutput;
    return outputData;
  }


  Data32 const *
    transform(Data32 const *inputSubMatrix,
	      subMatrixDescriptor const &inputDescription,
	      subMatrixDescriptor *outputDescription=NULL) 
  {
    subMatrixDescriptor myOutput;
    Data32 const *outputData = localTransform(inputSubMatrix,
					      inputDescription, &myOutput);
    if (nextFilter) {
      // FIXME - check that myOutput == *(inputDescription->next)
      outputData = 
	nextFilter->transform(outputData, *(inputDescription.next), outputDescription);
    } else {
      if (outputDescription) *outputDescription = myOutput;
    }
    return outputData;
  }

};

#endif
