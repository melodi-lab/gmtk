
/*
 * GMTK_ARMAFilter.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2012, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */

#ifndef GMTK_ARMAFILTER_H
#define GMTK_ARMAFILTER_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <assert.h>

#include "machine-dependent.h"
#include "GMTK_Filter.h"

#include "GMTK_SubmatrixDescriptor.h"


class ARMAFilter: public Filter {

  Data32   *buffer;     // result buffer
  unsigned  buffSize;   // in Data32's
  
  unsigned  order;
  unsigned  numFloats;
  unsigned  numFrames;  // total # frames in entire input stream if known, else 0

  // have to remember the previous order outputs
  unsigned  firstRememberedFrame;  // first frame # of remembered outputs
  float    *output;                // last order outputs
  unsigned  outputIndex;           // index (in output) of oldest remembered frame 

  float    *outputSum;             // sum of the remembered outputs
  unsigned  numRemembered;         // how many frames' worth of outputs are remembered
  bool      initialized;           // have I remembered anything yet?

 public:
  
  ARMAFilter(unsigned order, unsigned numFloats, Filter *nextFilter = NULL) 
    : buffer(NULL), buffSize(0), order(order), numFloats(numFloats), numFrames(0),
      firstRememberedFrame(0), output(0), numRemembered(0), initialized(false)
  {
    output = new float[order * numFloats];
    outputSum = new float[numFloats];
    this->nextFilter = nextFilter;
  }
  
  unsigned getTotalOutputFrames(unsigned totalInputFrames) {
    // reset for next segment
    numFrames = totalInputFrames;
    if (order > numFrames/2) 
      error("ERROR: In applying tranforms: ARMA filter order (%d) has to be less than "
	    "half the number of frames (%d).", order, numFrames);

    outputIndex = 0;
    numRemembered = 0;
    initialized = false;
    for (unsigned f=0; f < order; f+=1) {
      for (unsigned i=0; i < numFloats; i+=1) {
	output[numFloats * f + i] = 0.0;
      }
    }
    for (unsigned i=0; i < numFloats; i+=1) {
      outputSum[i] = 0.0;
    }

    if (nextFilter) {
      return nextFilter->getTotalOutputFrames(totalInputFrames);
    } else {
      return totalInputFrames;
    }
  }

  virtual ~ARMAFilter() {
    if (buffer)     free(buffer);
    if (output)     delete[] output;
    if (outputSum)  delete[] outputSum;
  }

  // The filter's client (e.g. inference) needs the 
  // [first,first+count)frames of the filter's output.
  // getRequiredInput() returns the portion of the fitler's 
  // input (e.g. ObservationFile) necessary to produce the 
  // requested output frames.
  //
  // Note that first and count describe the desired OUTPUT frames.
  // inputContinuous, inputDiscrete and inputTotalFrames are the number
  // of continuous & discrete features and frames in the Filter's INPUT.

  virtual subMatrixDescriptor *
    getRequiredInput(unsigned first, unsigned count, 
		     unsigned inputContinuous, unsigned inputDiscrete,
		     unsigned inputTotalFrames);

  // What will the output look like if the input described by
  // inputDescription is feed into this Filter?
  virtual subMatrixDescriptor
    describeLocalOutput(subMatrixDescriptor const &inputDescription);




  void getNextFrameInfo(unsigned &numNewIn, unsigned &dropOldIn, unsigned &numNewOut,
			unsigned inputContinuous, unsigned inputDiscrete,
			subMatrixDescriptor &input)
  {
    if (frameNum == 0) {
      numNewIn = order + 1;
      dropOldIn = 0;
      frameNum = 1;
    } else {
      numNewIn = 1;
      dropOldIn = 1;
    }
    numNewOut = 1;
    input.firstFrame = 0;
    input.numFrames = order + 1;
    input.historyFrames = 0;
    input.futureFrames = order;
    input.numContinuous = inputContinuous;
    input.numDiscrete = inputDiscrete;
    input.fullMatrixFrameCount = order+1;
    input.requestedFirst = 0;
    input.requestedCount = 1;
    input.next = NULL;
  }


  void getEOSFrameInfo(int numFramesShort, unsigned &numNewOut, subMatrixDescriptor &input) {
    assert(0 < numFramesShort && numFramesShort <= (int)(order+1));
    assert(numFramesShort == 1);
    numNewOut = 0;
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
		   subMatrixDescriptor *outputDescription=NULL);

};

#endif
