
/*
 * GMTK_ARMAFilter.cc
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2012, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "error.h"
#include "general.h"

#include "file_utils.h"
#include "GMTK_ARMAFilter.h"

  
  // The filter's client (e.g. inference) needs the 
  // [first,first+count)frames of the filter's output.
  // getRequiredInput() returns the portion of the fitler's 
  // input (e.g. ObservationFile) necessary to produce the 
  // requested output frames.
  //
  // Note that first and count describe the desired OUTPUT frames.
  // inputContinuous, inputDiscrete and inputTotalFrames are the number
  // of continuous & discrete features and total frames (not just 
  // [first,first+count), but all available input frames)in the Filter's INPUT.

subMatrixDescriptor *
ARMAFilter::getRequiredInput(unsigned first, unsigned count, 
			     unsigned inputContinuous, unsigned inputDiscrete,
			     unsigned inputTotalFrames)
{

  subMatrixDescriptor *nextFilterInput = NULL;
  if (nextFilter) {
    // the client's getting the data from the Filter(s) AFTER me, so
    // ask it what frames I need to provide it as input (it gets as
    // input the # of continuous and discrete features and frames
    // I produce as output).
    unsigned outputContinuous = inputContinuous;
    unsigned outputDiscrete   = inputDiscrete;
    unsigned outputFrames     = inputTotalFrames;
    nextFilterInput = nextFilter->getRequiredInput(first, count, 
                                                   outputContinuous, outputDiscrete,
                                                   outputFrames);
    assert(nextFilterInput);
    first = nextFilterInput->requestedFirst;
    count = nextFilterInput->requestedCount;
  } 
  unsigned requiredFirst = first;
  unsigned requiredLast  = first + count - 1 + order;
  if (numFrames > 0 && requiredLast > numFrames-1) 
    requiredLast = numFrames-1;
  unsigned requiredCount = requiredLast - requiredFirst + 1;

  return subMatrixDescriptor::getSMD(requiredFirst, requiredCount, 0, order,
				     inputContinuous, inputDiscrete, 
				     inputTotalFrames, first, count, 
				     nextFilterInput);
}

// What will the output look like if the input described by
// inputDescription is feed into this Filter?
subMatrixDescriptor
ARMAFilter::describeLocalOutput(subMatrixDescriptor const &inputDescription) {
  subMatrixDescriptor myOutput = inputDescription;
  myOutput.firstFrame = inputDescription.requestedFirst;
  myOutput.numFrames  = inputDescription.requestedCount;
  myOutput.historyFrames = 0;
  myOutput.futureFrames = 0;
  myOutput.next = NULL;
  return myOutput;
}

// Returns the filter's output given the inputSubMatrix
// described by inputDescription (presumably created by 
// getRequiredInput()).  
//
// If outputDescription is non-NULL, 
// it will describe the portion of filter's output returned 
Data32 const *
ARMAFilter::localTransform(Data32 const *inputSubMatrix,
			   subMatrixDescriptor const &inputDescription,
			   subMatrixDescriptor *outputDescription)
{
  subMatrixDescriptor myOutput = describeLocalOutput(inputDescription);
  unsigned numContinuous = inputDescription.numContinuous;
  unsigned numDiscrete  = inputDescription.numDiscrete;
  unsigned stride = numContinuous + numDiscrete;
  unsigned needed = stride * inputDescription.requestedCount;
  if (buffSize < needed) {
    buffer = (Data32 *) realloc(buffer, needed * sizeof(Data32));
    assert(buffer);
    buffSize = needed;
  }

  assert(numContinuous == numFloats);
  float *floatOut = (float *)buffer;
  float *floatIn  = (float *)inputSubMatrix;

  if (initialized && 
    // FIXME - [[ inputDescription.requestedFirst != 0 && ]]  below is a hack to support FilterStream,
    // which lies to the filter repeated asking for frame 0. The input frames are managed
    // such that the filter will produce the desired output. This error should still be triggered
    // by a backward pass or island algorithm, which are not supported unless non-constant 
    // memory observation input is used
      inputDescription.requestedFirst != 0 &&
      inputDescription.requestedFirst < numFrames - order &&
      inputDescription.requestedFirst != firstRememberedFrame + order)
  {
    error("ERROR: ARMA filter doesn't support random access (island or constantSpace options)"
	  " requested frame %u but expected frame %u\n", inputDescription.firstFrame,
	  firstRememberedFrame + order);
  }

  unsigned cnt = 0;
  if (!initialized) {

#if 0
    printf("-----------------------------------------------------\n");
    for (unsigned i=0; i < numContinuous; i+=1) {
      printf("%f ", outputSum[i]);
    }
    printf("\n");
    printf("=====================================================\n");
#endif

    // FIXME - [[ && inputDescription.firstFrame != 0 ]] below is a hack to support FilterStream,
    // which lies to the filter repeated asking for frame 0. The input frames are managed
    // such that the filter will produce the desired output. This error should still be triggered
    // by a backward pass or island algorithm, which are not supported unless non-constant 
    // memory observation input is used
    if (inputDescription.firstFrame != numRemembered && inputDescription.firstFrame != 0) {
      error("ERROR: initializing ARMA filter requires frame %u, but got frame %u\n", 
	    numRemembered, inputDescription.firstFrame);
    }
    for (unsigned f=0; f < inputDescription.requestedCount && numRemembered < order; f+=1, numRemembered+=1, cnt+=1) {
      for (unsigned i=0; i < numContinuous; i+=1) {
	float x = floatIn[stride * f + i];
	floatOut[stride * f + i] = x;
	output[numContinuous * numRemembered + i] = x;
	outputSum[i] = outputSum[i] + x;
      }
      for (unsigned i = numContinuous; i < stride; i+=1) {
	buffer[stride * f + i] = inputSubMatrix[stride * f + i];
      }	  
    }
    initialized = numRemembered == order;
#if 0
    if (initialized) {
      for (unsigned f=0; f < order; f+=1) {
	for (unsigned i=0; i < numContinuous; i+=1) {
	  printf("%f ", output[numContinuous * f + i]);
	}
	printf("\n");
      }
      printf("-----------------------------------------------------\n");
      for (unsigned i=0; i < numContinuous; i+=1) {
	printf("%f ", outputSum[i]);
      }
      printf("\n\n");
      printf("-----------------------------------------------------\n");
    }
#endif
  }
  
  if (cnt < inputDescription.requestedCount && !initialized) {
    error("ERROR: ARMA filter requested frame %u, but only the first %u frames have been processed\n", 
	  inputDescription.firstFrame + cnt, numRemembered);
  }

  for ( ; cnt < inputDescription.requestedCount; cnt += 1) {
    if (numFrames == 0 || inputDescription.firstFrame + cnt < numFrames - order) {
      for (unsigned i = 0; i < numContinuous; i+=1) {
	float x = outputSum[i] + floatIn[stride * cnt + i];
	for (unsigned tau=1; tau <= order; tau+=1) {
	  x += floatIn[stride * (cnt + tau) + i];
	}
	x /= (2*order + 1);
	floatOut[stride * cnt + i] = x;
	outputSum[i] = outputSum[i] - output[numFloats * outputIndex + i] + x;
	output[numFloats * outputIndex + i] = x;
      }
      firstRememberedFrame += 1;
      outputIndex = (outputIndex + 1) % order;
      for (unsigned i = numContinuous; i < stride; i+=1) {
	buffer[stride * cnt + i] = inputSubMatrix[stride * cnt + i];
      }	  
    } else { // just copy the last order frames
      for (unsigned i=0; i < stride; i += 1) {
	buffer[stride * cnt + i] = inputSubMatrix[stride * cnt + i];
      }
    }
  }

  if (outputDescription) *outputDescription = myOutput;
  return buffer;
}
