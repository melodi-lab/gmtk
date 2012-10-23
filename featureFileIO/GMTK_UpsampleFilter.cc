
/*
 * GMTK_UpsampleFilter.cc
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
#include "GMTK_UpsampleFilter.h"

  
  // The filter's client (e.g. inference) needs the 
  // [first,first+count)frames of the filter's output.
  // getRequiredInput() returns the portion of the fitler's 
  // input (e.g. ObservationFile) necessary to produce the 
  // requested output frames.
  //
  // Note that first and count describe the desired OUTPUT frames.
  // inputContinuous, inputDiscrete and inputTotalFrames are the number
  // of continuous & discrete features and frames in the Filter's INPUT.

subMatrixDescriptor *
UpsampleFilter::getRequiredInput(unsigned first, unsigned count, 
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
    unsigned outputFrames     = inputTotalFrames * upsample;
    nextFilterInput = nextFilter->getRequiredInput(first, count, 
                                                   outputContinuous, outputDiscrete,
                                                   outputFrames);
    assert(nextFilterInput);
    first = nextFilterInput->requestedFirst;
    count = nextFilterInput->requestedCount;
  } 
  unsigned reqFirst = first / upsample;
  unsigned reqLast  = (first + count - 1) / upsample;
  unsigned reqCount = reqLast - reqFirst + 1;

  return subMatrixDescriptor::getSMD(reqFirst, reqCount, 0, 0, 
				     inputContinuous, inputDiscrete, 
				     inputTotalFrames, first, count, 
				     nextFilterInput);
}

// What will the output look like if the input described by
// inputDescription is feed into this Filter?
subMatrixDescriptor
UpsampleFilter::describeLocalOutput(subMatrixDescriptor const &inputDescription) {
  subMatrixDescriptor myOutput = inputDescription;
  myOutput.firstFrame = inputDescription.requestedFirst;
  myOutput.numFrames  = inputDescription.requestedCount;
  myOutput.historyFrames = 0;
  myOutput.futureFrames = 0;
  myOutput.fullMatrixFrameCount *= upsample;
  myOutput.next = NULL;

#if 0
printf("UH dLO: ( [%u,%u)  -%u +%u  %u %u  %u  %u %u )\n", 
       inputDescription.firstFrame,
       inputDescription.firstFrame + inputDescription.numFrames,
       inputDescription.historyFrames, inputDescription.futureFrames,
       inputDescription.numContinuous, inputDescription.numDiscrete,
       inputDescription.fullMatrixFrameCount,
       inputDescription.requestedFirst, inputDescription.requestedCount);
printf("UH dLO: ( [%u,%u)  -%u +%u  %u %u  %u  %u %u )\n", 
       myOutput.firstFrame,
       myOutput.firstFrame + myOutput.numFrames,
       myOutput.historyFrames, myOutput.futureFrames,
       myOutput.numContinuous, myOutput.numDiscrete,
       myOutput.fullMatrixFrameCount,
       myOutput.requestedFirst, myOutput.requestedCount);
#endif
  return myOutput;
}

// Returns the filter's output given the inputSubMatrix
// described by inputDescription (presumably created by 
// getRequiredInput()).  
//
// If outputDescription is non-NULL, 
// it will describe the portion of filter's output returned 
Data32 const *
UpsampleFilter::localTransform(Data32 const *inputSubMatrix,
			       subMatrixDescriptor const &inputDescription,
			       subMatrixDescriptor *outputDescription)
{
  subMatrixDescriptor myOutput = describeLocalOutput(inputDescription);
  unsigned stride = inputDescription.numContinuous + inputDescription.numDiscrete;
  unsigned needed = stride * inputDescription.requestedCount;
  if (buffSize < needed) {
    buffer = (Data32 *) realloc(buffer, needed * sizeof(Data32));
    assert(buffer);
    buffSize = needed;
  }
  unsigned cnt = 0;
  unsigned firstFrameReps = upsample - (inputDescription.requestedFirst % upsample);

  for (unsigned i=0; i < firstFrameReps && i < inputDescription.requestedCount; i+=1, cnt+=1) {
    for (unsigned j=0; j < stride; j+=1) {
      buffer[cnt*stride+j]=inputSubMatrix[j];
    }
  }

  for (unsigned i=1; 
       i < inputDescription.numFrames && cnt < inputDescription.requestedCount; 
       i+=1)
  {
    for (unsigned k=0; k < upsample && cnt < inputDescription.requestedCount; k+=1, cnt+=1) {
      for (unsigned j=0; j < stride; j+=1) {
	buffer[cnt*stride+j]=inputSubMatrix[i*stride+j];
      }
    }
  }
  if (outputDescription) *outputDescription = myOutput;
  return buffer;
}
