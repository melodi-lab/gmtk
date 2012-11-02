
/*
 * GMTK_FIRFilter.cc
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2012, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */

#if HAVE_CONFIG_H
#include <config.h>
#endif
#if HAVE_HG_H
#include "hgstamp.h"
#endif

#include "error.h"
#include "general.h"
VCID(HGID)

#include <string.h>
#include <stdio.h>

#include "GMTK_FIRFilter.h"

FIRFilter::FIRFilter(char const *fileName, Filter *nextFilter) 
  : Filter(nextFilter)
{
  FILE *f = fopen(fileName, "r");
  if (!f) {
    error("FIRFilter: cannot open '%s' for reading\n", fileName);
  }
  if (fscanf(f, "%u %u", &order, &numFeatures) != 2) {
    error("FIRFilter: error reading FIR filter size from '%s'\n", fileName);
  }

  c = new float[numFeatures];
  if (!c) {
    error("FIRFilter: failed to allocate constant vector\n");
  }
  unsigned Bsize = (order+1) * numFeatures;
  B = new float[Bsize];
  if (!B) {
    error("FIRFilter: failed to allocate coefficient matrix\n");
  }
  for (unsigned i=0; i < Bsize; i+=1) {
    if (fscanf(f,"%e", B+i) != 1) {
      error("FIRFilter: error reading coefficient matrix\n");
    }
  }
  for (unsigned i=0; i < numFeatures; i+=1) {
    if (fscanf(f,"%e", c+i) != 1) {
      error("FIRFilter: error reading constant vector\n");
    }
  }
  fclose(f);
  buffer=NULL; buffSize=0;
}


subMatrixDescriptor *
FIRFilter::getRequiredInput(unsigned first, unsigned count, 
			    unsigned inputContinuous, unsigned inputDiscrete,
			    unsigned inputTotalFrames)
{
  subMatrixDescriptor *nextFilterInput = NULL;
  if (nextFilter) {
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
  unsigned history;
  unsigned requiredFirst = first;
  unsigned requiredCount = count;
  if (requiredFirst >= order) {
    requiredFirst -= order;
    requiredCount += order;
    history = order;
  } else {
    requiredCount += first;
    history = first;
    requiredFirst = 0;
  }
  return subMatrixDescriptor::getSMD(requiredFirst, requiredCount, history, 0, inputContinuous, 
				     inputDiscrete, inputTotalFrames, first, count,
				     nextFilterInput);
}


void
FIRFilter::getNextFrameInfo(unsigned &numNewIn, unsigned &dropOldIn, unsigned &numNewOut,
			    unsigned inputContinuous, unsigned inputDiscrete,
			    subMatrixDescriptor &input)
{
  numNewIn = 1;
  dropOldIn = frameNum == 0 ? 0 : 1; 
  numNewOut= 1;
  subMatrixDescriptor *d = getRequiredInput(frameNum, 1, inputContinuous, inputDiscrete, order+1);
  input = *d;
  subMatrixDescriptor::freeSMD(d);
  frameNum += (frameNum > order) ? 0 : 1; // only keep track up to order		     
}


subMatrixDescriptor
FIRFilter::describeLocalOutput(subMatrixDescriptor const &inputDescription) {
  subMatrixDescriptor myOutput = inputDescription;
  myOutput.firstFrame += inputDescription.historyFrames;
  myOutput.numFrames  -= inputDescription.historyFrames;
  myOutput.historyFrames = 0;
  myOutput.futureFrames = 0;
  myOutput.fullMatrixFrameCount -= inputDescription.historyFrames;
  myOutput.next = NULL;
  return myOutput;
}


Data32 const *
FIRFilter::localTransform(Data32 const *inputSubMatrix, 
			  subMatrixDescriptor const &inputDescription,
			  subMatrixDescriptor *outputDescription) 
{
  if (numFeatures != inputDescription.numContinuous) {
    error("FIRFilter: filter expects %u features, but input has %u\n", 
	  numFeatures, inputDescription.numContinuous);
  }
  unsigned stride = inputDescription.numContinuous + inputDescription.numDiscrete;

  subMatrixDescriptor myOutput = describeLocalOutput(inputDescription);
  assert(myOutput.numContinuous == numFeatures);
  assert(myOutput.numDiscrete == inputDescription.numDiscrete);

  unsigned needed = stride * myOutput.numFrames;

  //  printf("xfrm in  [%u, %u, %u, %u]  :  stride %u  :  out [%u, %u]\n", inputDescription.firstFrame, inputDescription.historyFrames, inputDescription.numFrames, inputDescription.futureFrames, stride, myOutput.firstFrame, myOutput.numFrames);
  if (needed > buffSize) {
    buffer = (Data32 *) realloc(buffer, needed * sizeof(Data32));
    assert(buffer);
    buffSize = needed;
  }

  for (unsigned i=0; i < myOutput.numFrames; i+=1) {
    float *outputCont = (float *)(buffer + i * stride);
    if (c) {
      memcpy(outputCont, c, myOutput.numContinuous * sizeof(Data32));
    } else {
      memset(outputCont, 0, myOutput.numContinuous * sizeof(Data32));
    }
    Uint32 *outputDisc = (Uint32 *)(buffer + i * stride + myOutput.numContinuous);
#if 0
    Uint32 *inputDisc  = (Uint32 *)(inputSubMatrix + 
				    (i + inputDescription.historyFrames) * stride + 
				    inputDescription.numDiscrete);
#else
    Uint32 *inputDisc  = (Uint32 *)(inputSubMatrix + 
				    (i + inputDescription.historyFrames) * stride + 
				    inputDescription.numContinuous);
#endif
    memcpy(outputDisc, inputDisc, myOutput.numDiscrete * sizeof(Data32));
  }

  if (B==NULL) {
    B = new float[myOutput.numContinuous];
    assert(B);
    for (unsigned i=0; i < myOutput.numContinuous; i+=1)
      B[i] = 1.0f;
  }
  unsigned availableHistory = inputDescription.historyFrames;
  for (unsigned outr=0; outr < myOutput.numFrames; outr+=1) {
    float *outputRow = (float *)(buffer + outr * stride);
    for (unsigned h=0; h <= availableHistory; h+=1) {
      float *inputRow  = (float *)(inputSubMatrix + (outr + inputDescription.historyFrames - h) * stride);
      for (unsigned j=0; j < myOutput.numContinuous; j+=1) {
	outputRow[j] += inputRow[j] * (B + h * myOutput.numContinuous)[j];
      }
    }
  }
  availableHistory = (availableHistory < order) ? availableHistory + 1  :  order;

  if (outputDescription) *outputDescription = myOutput;
  return buffer;
}
