
/*
 * GMTK_UpsampleSmoothFilter.cc
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
#include "GMTK_UpsampleSmoothFilter.h"

  
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
UpsampleSmoothFilter::getRequiredInput(unsigned first, unsigned count, 
				 unsigned inputContinuous, unsigned inputDiscrete,
				 unsigned inputTotalFrames)
{

  // requires 1 future frame if available...

  subMatrixDescriptor *nextFilterInput = NULL;
  if (nextFilter) {
    // the client's getting the data from the Filter(s) AFTER me, so
    // ask it what frames I need to provide it as input (it gets as
    // input the # of continuous and discrete features and frames
    // I produce as output).
    unsigned outputContinuous = inputContinuous;
    unsigned outputDiscrete   = inputDiscrete;
    unsigned outputFrames     = 1 + (inputTotalFrames-1) * upsample;
    nextFilterInput = nextFilter->getRequiredInput(first, count, 
                                                   outputContinuous, outputDiscrete,
                                                   outputFrames);
    assert(nextFilterInput);
    first = nextFilterInput->requestedFirst;
    count = nextFilterInput->requestedCount;
  } 
  unsigned requiredFirst = first / upsample;
  unsigned requiredLast  = (first + count - 1) / upsample;
  unsigned futureFrames = (requiredLast < inputTotalFrames - 1) ? 1 : 0; // can't smooth last input frame
  unsigned requiredCount = requiredLast - requiredFirst + 1 + futureFrames;

  return subMatrixDescriptor::getSMD(requiredFirst, requiredCount, 0, futureFrames,
				     inputContinuous, inputDiscrete, 
				     inputTotalFrames, first, count, 
				     nextFilterInput);
}

// What will the output look like if the input described by
// inputDescription is feed into this Filter?
subMatrixDescriptor
UpsampleSmoothFilter::describeLocalOutput(subMatrixDescriptor const &inputDescription) {
  subMatrixDescriptor myOutput = inputDescription;
  myOutput.firstFrame = inputDescription.requestedFirst;
  myOutput.numFrames  = inputDescription.requestedCount;
  myOutput.historyFrames = 0;
  myOutput.futureFrames = 0;

  myOutput.fullMatrixFrameCount = 1 + (inputDescription.fullMatrixFrameCount - 1) * upsample;
  myOutput.next = NULL;

#if 0
printf("US dLO: ( [%u,%u)  -%u +%u  %u %u  %u  %u %u )\n", 
       inputDescription.firstFrame,
       inputDescription.firstFrame + inputDescription.numFrames,
       inputDescription.historyFrames, inputDescription.futureFrames,
       inputDescription.numContinuous, inputDescription.numDiscrete,
       inputDescription.fullMatrixFrameCount,
       inputDescription.requestedFirst, inputDescription.requestedCount);
printf("US dLO: ( [%u,%u)  -%u +%u  %u %u  %u  %u %u )\n", 
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
UpsampleSmoothFilter::localTransform(Data32 const *inputSubMatrix,
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
  unsigned cnt = 0;

  // Each input frame appears upsample times in the complete filter
  // ouput, but the requested portion may start in the middle of a
  // sequence of repeated frames. firstFrameSkip is the # of repeats
  // of the first input frame that are not requested in the output. 
  // Likewise, firstFrameReps is the # of repeats of the first input
  // frame that are included in the output.
  unsigned firstFrameSkip = inputDescription.requestedFirst % upsample;
  unsigned firstFrameReps = upsample - firstFrameSkip > inputDescription.requestedCount ?
    inputDescription.requestedCount : upsample - firstFrameSkip;

  double inc;
  double additive_const;

  float *input  = (float *)inputSubMatrix;
  float *output = (float *)buffer;
  
  // This should work OK if the unrepeated final input frame (of the full
  // matrix) is the first requested - in that case the requested count must be 1

  for (unsigned i=0; i < firstFrameReps; i+=1, cnt+=1) {
    for (unsigned j=0; j < numContinuous; j+=1) {
      additive_const =  input[stride + j] - input[j];
      additive_const /= upsample;
      inc = (firstFrameSkip + i) * additive_const;
      output[cnt*stride+j] = input[j] + (float)inc;
    }
    for (unsigned j=0; j < numDiscrete; j+=1) {
      buffer[cnt*stride + numContinuous + j] = inputSubMatrix[numContinuous + j];
    }
  }

  unsigned repeatingFrames = inputDescription.numFrames - inputDescription.futureFrames;
  bool     doFinalFrame;
  if (inputDescription.requestedFirst + inputDescription.requestedCount == 
      //inputDescription.fullMatrixFrameCount
      1 + (inputDescription.fullMatrixFrameCount - 1) * upsample)
  {
    repeatingFrames -= 1;      // don't repeat last frame
    doFinalFrame     = true;   // just append it once
  } else {
    doFinalFrame     = false;  // repeat all frames normally
  }

  for (unsigned i=1; i < repeatingFrames && cnt < inputDescription.requestedCount; i+=1)   {
    for (unsigned k=0; k < upsample && cnt < inputDescription.requestedCount; k+=1, cnt+=1) {
      for (unsigned j=0; j < stride; j+=1) {
	additive_const =  input[(i+1)*stride + j] - input[i*stride + j];
	additive_const /= upsample;
	inc = k * additive_const;
	output[cnt*stride+j] = input[i*stride+j] + (float)inc;
      }
      for (unsigned j=0; j < numDiscrete; j+=1) {
	buffer[cnt*stride + numContinuous + j] = inputSubMatrix[i*stride + numContinuous + j];
      }
    }
  }
  
  if (doFinalFrame && cnt < inputDescription.requestedCount) { // tack on 1 copy of last frame if needed
    for (unsigned j = 0; j < stride; j += 1) {
      buffer[cnt*stride + j] = inputSubMatrix[(inputDescription.numFrames-1) * stride + j];
    }
  }
  assert(cnt == inputDescription.requestedCount);
  if (outputDescription) *outputDescription = myOutput;
  return buffer;
}
