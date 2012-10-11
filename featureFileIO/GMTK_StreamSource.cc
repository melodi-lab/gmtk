/*
 * GMTK_StreamSource.cc
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2012, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */

#include <string.h>

#include "debug.h"
#include "error.h"

#include "GMTK_StreamSource.h"

StreamSource::StreamSource() {
  cookedBuffer = NULL;
  cookedBuffSize = 0;
  maxCookedFrames = 0;
  currentCookedFrames = 0;
  numFramesInSegment = 0;
  segmentNum = -1;
  _startSkip = 0;
  nFloat = 0;
  nInt = 0;
  nFeatures = 0;
}

StreamSource::StreamSource(unsigned nStreams,
			   ObservationStream *stream[], 
			   unsigned queueLength,
			   char *filterStr,
			   unsigned startSkip)
  : cookedBuffSize(queueLength),
    segmentNum(-1),
    _startSkip(startSkip)
{
  if (nStreams > 1) {
    this->stream = new MergeStream(stream, nStreams);
  } else {
    this->stream = stream[0];
  }
  
  if (filterStr) {
    Filter *filters = instantiateFilters(filterStr, this->stream->numLogicalContinuous(),this->stream->numLogicalDiscrete());
    while (filters) {
      // unlike the FilterFile, each filter needs its own independent FilterStream
      Filter *next = filters->nextFilter;
      filters->nextFilter = NULL;
      this->stream = new FilterStream(this->stream, filters);
      filters = next;
    }
  }
  initialize(queueLength, startSkip);
}

void
StreamSource::initialize(unsigned queueLength, unsigned startSkip)
{
  assert(stream);
  cookedBuffSize = queueLength;
  nFloat = stream->numLogicalContinuous();
  nInt   = stream->numLogicalDiscrete();
  nFeatures = nFloat + nInt;
  segmentNum = -1;
  _startSkip = startSkip;

  cookedBuffer = new Data32[cookedBuffSize];
  maxCookedFrames = cookedBuffSize / nFeatures;
  currentCookedFrames = 0;
  firstCookedFrameNum = 0;
  firstCookedFrameIndex = 0;
  numFramesInSegment = 0;
}


void
StreamSource::preloadFrames(unsigned nFrames) {
//fprintf(stdout, "preloadFrames(%u)\n", nFrames);
  if (nFrames >  maxCookedFrames / 2) {
      error("ERROR: StreamSource::enqueueFrames -streamBufferSize must be at least %u MB\n",
	    1 + nFrames * 2 / (1024*1024) );
  }
  currentCookedFrames = 0;
  firstCookedFrameNum = 0;
  firstCookedFrameIndex = 0;
  numFramesInSegment = 0;
  segmentNum += 1;
  
  Data32 *newFrameDest = cookedBuffer;
  unsigned numEnqueued;
  for (numEnqueued=0; numEnqueued < nFrames; numEnqueued+=1) {
    Data32 const *newFrame = stream->getNextLogicalFrame();
    if (!newFrame) {
      // discovered segment length
      numFramesInSegment = currentCookedFrames;
      break;
    }
    memcpy(newFrameDest, newFrame, nFeatures*sizeof(Data32));
    newFrameDest += nFeatures;
    currentCookedFrames += 1;
    infoMsg(IM::ObsStream, IM::High, "preload: [%u,%u] @ [%u,%u / %u]\n", 
	    firstCookedFrameNum, firstCookedFrameNum+currentCookedFrames-1,
	    firstCookedFrameIndex, firstCookedFrameIndex+currentCookedFrames-1,
	    maxCookedFrames);
  }
}


Data32 const *
StreamSource::loadFrames(unsigned first, unsigned count) {
//fprintf(stdout, "loadFrames(%u,%u)\n", first, count);  
  first += _startSkip;
  if (count > maxCookedFrames) {
    error("ERROR: StreamSource::loadFrames: requested %u frames, but the frame queue can only hold %u. Increase -streamBufferSize\n", count, maxCookedFrames);
  }

  if (numFramesInSegment > 0 && first + count > numFramesInSegment) {
    error("ERROR: StreamSource::loadFrames: requested frame %u, but the segment only has %u frames\n", first+count, numFramesInSegment);
  }

  if (first < firstCookedFrameNum) {
    error("ERROR: StreamSource::loadFrames: requested frame %u which is no longer available; the earliest available frame is %u\n", first, firstCookedFrameNum);
  }

  if (first > firstCookedFrameNum + currentCookedFrames) {
//fprintf(stdout, "  first %u  num %u\n", firstCookedFrameNum, currentCookedFrames);
    error("ERROR: StreamSource::loadFrames: requested frames [%u,%u] would require skipping %u frames past [%u,%u]\n", 
	  first, first + count -1, 
	  first - firstCookedFrameNum - currentCookedFrames,
	  firstCookedFrameNum, firstCookedFrameNum+currentCookedFrames-1);
  }
  
  if (first + count <= firstCookedFrameNum + currentCookedFrames) {
    // all requested frames are already in the queue
//fprintf(stdout, "StreamSource > [%u,%u)  all frames in queue\n", first,first+count);
    infoMsg(IM::ObsStream, IM::High, "loadFrames: loading [%u,%u] from [%u,%u]+%u @ %u\n",
	    first, first+count-1,
	    firstCookedFrameNum, firstCookedFrameNum+currentCookedFrames-1,
	    first - firstCookedFrameNum,
	    firstCookedFrameIndex);
    return cookedBuffer + nFeatures * (firstCookedFrameIndex + first - firstCookedFrameNum);
  }

  error("StreamSource::loadFrames requested frames [%u,%u], but only [%u,%u] are available\n",
	first, first+count-1, firstCookedFrameNum, firstCookedFrameNum+currentCookedFrames-1);
  return NULL;
}



unsigned 
StreamSource::enqueueFrames(unsigned nFrames) {
//fprintf(stdout, "enqueueFrames(%u)\n", nFrames);
  if (numFramesInSegment != 0)
    return 0; // current segment's done - call preloadFrames() to reset for the next

  if (nFrames > currentCookedFrames) {
    error("ERROR: StreamSource::enqueueFrames doesn't support enqueuing more than %u frames\n",
	  currentCookedFrames);
  }

  if (firstCookedFrameIndex + currentCookedFrames + nFrames > maxCookedFrames) {
    // need to move active frames back to start of buffer

    infoMsg(IM::ObsStream, IM::High, "enqueueFrames: moving [%u,%u] @ %u to start of buffer\n",
	    firstCookedFrameNum, firstCookedFrameNum+currentCookedFrames-1, firstCookedFrameIndex);

    // verify it's safe to memcpy (no overlap)
    if (maxCookedFrames / 2 < currentCookedFrames) {
      error("ERROR: StreamSource::enqueueFrames -streamBufferSize must be at least %u MB\n",
	    1 + currentCookedFrames * 2 / (1024*1024) );
    }
    Data32 *newFirstFrame = cookedBuffer + firstCookedFrameIndex * nFeatures;

    unsigned bytesToMove = currentCookedFrames * nFeatures * sizeof(Data32);

    memcpy(cookedBuffer, newFirstFrame, bytesToMove);
    
    firstCookedFrameIndex = 0;
  }
  Data32 *newFrameDest = cookedBuffer + (firstCookedFrameIndex + currentCookedFrames) * nFeatures;
  unsigned numEnqueued;
  for (numEnqueued=0; numEnqueued < nFrames; numEnqueued+=1) {
    Data32 const *newFrame = stream->getNextLogicalFrame();
    if (!newFrame) {
      // discovered segment length
      numFramesInSegment = firstCookedFrameNum + currentCookedFrames;
      break;
    }
    if (IM::messageGlb(IM::ObsStream,IM::High)) {
      printf("enqueue frame %u:", firstCookedFrameNum + currentCookedFrames);
      for (unsigned i=0; i < stream->numContinuous(); i+=1) {
	printf(" %f", ((float *)newFrame)[i]);
      }
      for (unsigned i=0; i < stream->numDiscrete(); i+=1) {
	printf(" %u", ((unsigned *)newFrame)[i]);
      }
      printf("\n");
    }

    memcpy(newFrameDest, newFrame, nFeatures*sizeof(Data32));
    newFrameDest += nFeatures;
    firstCookedFrameIndex += 1;
    firstCookedFrameNum += 1;
    infoMsg(IM::ObsStream, IM::High, "enqueue: [%u,%u] @ [%u,%u / %u]\n", 
	    firstCookedFrameNum, firstCookedFrameNum+currentCookedFrames-1,
	    firstCookedFrameIndex, firstCookedFrameIndex+currentCookedFrames-1,
	    maxCookedFrames);
  }
  return numEnqueued;
}

