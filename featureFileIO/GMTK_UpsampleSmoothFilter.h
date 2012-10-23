
/*
 * GMTK_UpsampleSmoothFilter.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2012, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */

#ifndef GMTK_UPSAMPLESMOOTHFILTER_H
#define GMTK_UPSAMPLESMOOTHFILTER_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <assert.h>

#include "machine-dependent.h"
#include "GMTK_Filter.h"

#include "GMTK_SubmatrixDescriptor.h"


class UpsampleSmoothFilter: public Filter {

  Data32   *buffer;     // result buffer
  unsigned  buffSize;   // in Data32's

  unsigned upsample;

 public:

  UpsampleSmoothFilter(unsigned upsample, Filter *nextFilter = NULL) 
    : buffer(NULL), buffSize(0), upsample(upsample+1)
  {
    this->nextFilter = nextFilter;
  }
  
  
  virtual ~UpsampleSmoothFilter() { 
    if (buffer)     free(buffer);
  }


  virtual unsigned getTotalOutputFrames(unsigned totalInputFrames) {
    unsigned frames = 1 + (totalInputFrames - 1) * upsample;
    if (nextFilter) {
      return nextFilter->getTotalOutputFrames(frames);
    } else {
      return frames;
    }
  }

  void getNextFrameInfo(unsigned &numNewIn, unsigned &dropOldIn, unsigned &numNewOut,
			unsigned inputContinuous, unsigned inputDiscrete,
			subMatrixDescriptor &input)
  {
    if (frameNum == 0) {
      numNewIn = 2;
      dropOldIn = 0;
      frameNum = 1;
    } else {
      numNewIn = 1;
      dropOldIn = 1;
    }
    numNewOut = upsample;
    input.firstFrame = 0;
    input.numFrames = 2;
    input.historyFrames = 0;
    input.futureFrames = 1;
    input.numContinuous = inputContinuous;
    input.numDiscrete = inputDiscrete;
    input.fullMatrixFrameCount = 2;
    input.requestedFirst = 0;
    input.requestedCount = upsample;
    input.next = NULL;
  }


  void getEOSFrameInfo(int numFramesShort, unsigned &numNewOut, subMatrixDescriptor &input) {
    if (numFramesShort == 1) {
      numNewOut = 1;
      input.numFrames = 1;
      input.fullMatrixFrameCount = 1;
      input.requestedCount = 1;
    } else {
      numNewOut = 0;
    }
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
