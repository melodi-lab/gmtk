
/*
 * GMTK_FilterStream.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2012, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */

#ifndef GMTK_FILTERSTREAM_H
#define GMTK_FILTERSTREAM_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include <stdio.h>
using namespace std;

#include "machine-dependent.h"
#include "debug.h"
#include "error.h"

#include "GMTK_Filter.h"
#include "GMTK_ObservationStream.h"


class FilterStream: public ObservationStream {
  ObservationStream *incomingStream;
  Filter *filter;

  Data32 *inputBuffer;
  unsigned numInputBufferFrames; // # frames currently in buffer
  unsigned inputBufferCapacity; // in frames

  Data32 *outputBuffer;
  unsigned numOutputBufferFrames;  // # frames currently in buffer
  unsigned outputBufferCapacity; // in frames

  unsigned outputFrame;
  unsigned inStride, outStride;

  bool eos;
  bool lastFrame;                // true iff we cooked the last input frame in the current segment

 public:

  FilterStream() 
    : incomingStream(NULL), filter(NULL),
      inputBuffer(NULL), numInputBufferFrames(0), inputBufferCapacity(0), 
      outputBuffer(NULL), numOutputBufferFrames(0), outputBufferCapacity(0),
      outputFrame(0), inStride(0), outStride(0), eos(true), lastFrame(true)
  {}


  FilterStream(ObservationStream *stream, Filter *filter)
    : incomingStream(stream), filter(filter),
      inputBuffer(NULL), numInputBufferFrames(0), inputBufferCapacity(0),
      outputBuffer(NULL), numOutputBufferFrames(0), outputBufferCapacity(0),
      outputFrame(0), eos(false), lastFrame(false)
    {
      inStride = stream->numLogicalContinuous() + stream->numLogicalDiscrete();
      subMatrixDescriptor wholeSegment(0U, 1U, 0U, 0U,
				       stream->numLogicalContinuous(),
				       stream->numLogicalDiscrete(), 1U, 0U, 0U);
      subMatrixDescriptor output = filter->describeOutput(wholeSegment);
      nFloat = output.numContinuous;
      nInt   = output.numDiscrete;
      outStride = nFloat + nInt;
      logicalFrameData = new Data32[outStride];
    }

  ~FilterStream() {
    if (filter) delete filter;
    if (incomingStream) delete incomingStream;
    if (inputBuffer) free(inputBuffer);
    if (outputBuffer) free(outputBuffer);
  }

  bool EOS() {return eos;}

  Data32 const *getNextFrame() {
    if (outputFrame < numOutputBufferFrames) {
      // needed frame is already in the output buffer
      infoMsg(IM::ObsStream,IM::High-2,"FilterStream: cache hit %u / %u\n", outputFrame, numOutputBufferFrames);
      frameData = outputBuffer + outStride * outputFrame++;
      return frameData;
    }

    if (lastFrame) { // all input frames in this segment have been cooked
      // reset for next segment 
      numInputBufferFrames = 0;
      numOutputBufferFrames = 0;
      filter->nextStreamSegment();
      lastFrame = false;
      return NULL;
    }

    // need to cook more raw frames
    unsigned numNewIn, dropOldIn, numNewOut;
    subMatrixDescriptor inputDesc, outputDesc;
    
    filter->getNextFrameInfo(numNewIn, dropOldIn, numNewOut, 
			     incomingStream->numLogicalContinuous(), 
			     incomingStream->numLogicalDiscrete(), inputDesc);

    assert(dropOldIn <= numInputBufferFrames);
    
    infoMsg(IM::ObsStream,IM::High-2,"FilterStream: queue %u - %u + %u\n", numInputBufferFrames, dropOldIn, numNewIn);
    unsigned newInputSize = numInputBufferFrames - dropOldIn + numNewIn;
    if (newInputSize > inputBufferCapacity) {
      infoMsg(IM::ObsStream,IM::High-2,"FilterStream: resize queue %u -> %u\n", inputBufferCapacity, newInputSize);
      inputBuffer = (Data32 *)realloc(inputBuffer, newInputSize * inStride * sizeof(Data32));
      if (!inputBuffer) {
	error("FilterStream::getNextFrame failed to resize input buffer\n");
	inputBufferCapacity = newInputSize;
      }
      inputBufferCapacity = newInputSize;
    }

    Data32 *newFrameDest;
    if (dropOldIn == numInputBufferFrames) {
      // dropping all old frames, so just over-write them
      infoMsg(IM::ObsStream,IM::High-2,"FilterStream: overwriting old input buffer\n");
      newFrameDest = inputBuffer;
    } else {
      // need to shift some frames to start of input buffer
      unsigned numPreserved = numInputBufferFrames - dropOldIn;
      if (numPreserved <= dropOldIn) {
	// no overlap, so shift with memcpy
	infoMsg(IM::ObsStream,IM::High-2,"FilterStream: memcpy %u frames\n", numPreserved);
	memcpy(inputBuffer, inputBuffer + dropOldIn * inStride, numPreserved * inStride * sizeof(Data32));
      } else {
	// overlap, so shift with memmove
	infoMsg(IM::ObsStream,IM::High-2,"FilterStream: memmove %u frames\n", numPreserved);
	memmove(inputBuffer, inputBuffer + dropOldIn * inStride, numPreserved * inStride * sizeof(Data32));
      }
      newFrameDest = inputBuffer + numPreserved * inStride;
    }
    numInputBufferFrames -= dropOldIn;

    // enqueue new frames from incomingStream
    unsigned numEnqueuedFrames = 0;
    for (unsigned i=0; i < numNewIn; i+=1) {
      Data32 const *newFrame = incomingStream->getNextLogicalFrame();
      if (!newFrame) {
	eos = incomingStream->EOS();
	lastFrame = true;
	break;
      }
      infoMsg(IM::ObsStream,IM::High-2,"FilterStream: enqueue %u / %u frames @ %u\n", 
	      i, numNewIn, numInputBufferFrames);
      memcpy(newFrameDest, newFrame, inStride * sizeof(Data32));
      newFrameDest += inStride;
      numInputBufferFrames += 1;
      numEnqueuedFrames += 1;
    }



    if (lastFrame) {  // segment had less than numNewIn frames left
      int numFramesShort = (int)numNewIn - (int)numEnqueuedFrames;
      filter->getEOSFrameInfo(numFramesShort, numNewOut, inputDesc);
      if (numNewOut == 0) {
	outputFrame = 1;   // ensure outputBuffer miss
	numOutputBufferFrames = 0;
	return NULL;
      }
    }


    if (numNewOut > outputBufferCapacity) {
      outputBuffer = (Data32 *)realloc(outputBuffer, numNewOut * outStride * sizeof(Data32));
      if (!outputBuffer) {
	error("FilterStream::getNextFrame failed to resize output buffer\n");
      }
      outputBufferCapacity = numNewOut;
    }
#if 0
    if (IM::messageGlb(IM::ObsStream, IM::High-2)) {
      for (unsigned i=0; i < numInputBufferFrames; i+=1) {
	printf("in  %u:", i);
	for (unsigned j=0; j < incomingStream->numLogicalContinuous(); j+=1) {
	  printf(" %f", ((float *)inputBuffer)[i * inStride + j]);
	}
	printf("\n");
      }
      for (unsigned i=0; i < numNewOut; i+=1) {
	printf("out %u:", i);
	for (unsigned j=0; j < nFloat; j+=1) {
	  printf(" %f", ((float *)outputBuffer)[i * outStride + j]);
	}
	printf("\n");
      }
      printf("---------------------------------------------------------\n");
    }
#endif


    // use outputDesc.numFrames
    // may be 0 -> return NULL

    Data32 const *filterOutput = filter->transform(inputBuffer, inputDesc, &outputDesc);
    memcpy(outputBuffer, filterOutput, outputDesc.numFrames * outStride * sizeof(Data32));

#if 0
    if (IM::messageGlb(IM::ObsStream, IM::High-2)) {
      for (unsigned i=0; i < numInputBufferFrames; i+=1) {
	printf("in  %u:", i);
	for (unsigned j=0; j < incomingStream->numLogicalContinuous(); j+=1) {
	  printf(" %f", ((float *)inputBuffer)[i * inStride + j]);
	}
	printf("\n");
      }
      for (unsigned i=0; i < numNewOut; i+=1) {
	printf("out %u:", i);
	for (unsigned j=0; j < nFloat; j+=1) {
	  printf(" %f", ((float *)outputBuffer)[i * outStride + j]);
	}
	printf("\n");
      }
      printf("=========================================================\n");
    }

    infoMsg(IM::ObsStream,IM::High-2, "FilterStream: produced %u (%u) frames\n", numNewOut, outputDesc.numFrames);
#endif
    outputFrame = 1; // returning first new output frame now
    numOutputBufferFrames = numNewOut;

    frameData = outputBuffer;
    return frameData;
  }

};

#endif

