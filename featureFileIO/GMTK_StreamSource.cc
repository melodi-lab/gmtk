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

#include "error.h"

#include "GMTK_StreamSource.h"

StreamSource::StreamSource() {
  cookedBuffer = NULL;
  cookedBuffSize = 0;
  maxCookedFrames = 0;
  currentCookedFrames = 0;
  numFramesInSegment = 0;
  rawBuffer = NULL;
  rawBuffSize = 0;
  maxRawFrames = 0;
  currentRawFrames = 0;
  nFloat = 0;
  nInt = 0;
  nFeatures = 0;
  nStreams = 0;
  floatStart = NULL;
  intStart = NULL;
  filter = NULL;
}

StreamSource::StreamSource(unsigned nStreams,
			   ObservationStream *stream[], 
			   unsigned queueLength,
			   Filter *filter,
			   unsigned startSkip)
  : cookedBuffSize(queueLength),
    nStreams(nStreams),
    stream(stream),
    filter(filter),
    _startSkip(startSkip)
    
{
  initialize(nStreams, stream, queueLength, filter, startSkip);
}

void
StreamSource::initialize(unsigned nStreams, ObservationStream *stream[], 
			 unsigned queueLength, Filter *filter, unsigned startSkip)
{
  assert(stream);
  cookedBuffSize = queueLength;
  this->nStreams = nStreams;
  this->stream = new ObservationStream *[nStreams];
  unsigned floatCount = 0, intCount = 0;
  for (unsigned i=0; i < nStreams; i+=1) {
    assert(stream[i]);
    this->stream[i] = stream[i];
    floatCount += stream[i]->numLogicalContinuous();
    intCount   += stream[i]->numLogicalDiscrete();
  }
  rawFloats = floatCount;
  rawInts   = intCount;
  numRawFeatures = floatCount + intCount;
  this->filter = filter;
  if (filter) {
    subMatrixDescriptor wholeSegment(0U, 0U, 0U, 0U, floatCount, intCount, 0u);
    subMatrixDescriptor output = filter->describeOutput(wholeSegment);
    nFloat = output.numContinuous;
    nInt   = output.numDiscrete;
  } else {
    nFloat = floatCount;
    nInt   = intCount;
  }
  nFeatures = nFloat + nInt;
  _startSkip = startSkip;

  cookedBuffer = new Data32[cookedBuffSize];
  assert(cookedBuffer);
  maxCookedFrames = cookedBuffSize / nFeatures;
  currentCookedFrames = 0;
  firstCookedFrameNum = 0;

  // FIXME - pass or compute raw buff size
#define RAW_BUF_SIZE (4*1024*1024)
  rawBuffSize = RAW_BUF_SIZE;
  rawBuffer = new Data32[rawBuffSize];
  assert(rawBuffer);
  maxRawFrames = rawBuffSize / numRawFeatures;
  currentRawFrames = 0;
  firstRawFrameNum = 0;
  floatStart   = new Data32 *[nStreams];
  intStart     = new Data32 *[nStreams];
  assert(floatStart && intStart);
  unsigned offset = 0;
  for (unsigned i=0; i < nStreams; i+=1) {
    floatStart[i] = rawBuffer + offset;
    offset += stream[i]->numLogicalContinuous();
  }
  for (unsigned i=0; i < nStreams; i+=1) {
    intStart[i] = rawBuffer + offset;
    offset += stream[i]->numLogicalDiscrete();
  }

  numFramesInSegment = 0;
}


void
StreamSource::preloadFrames(unsigned nFrames) {
//fprintf(stdout, "preloadFrames(%u)\n", nFrames);
  assert(nFrames < maxCookedFrames);
  currentRawFrames = 0;
  firstRawFrameNum = 0;
  currentCookedFrames = 0;
  firstCookedFrameNum = 0;
  numFramesInSegment = 0;
  assert(enqueueFrames(nFrames) == nFrames);
}


Data32 const *
StreamSource::loadFrames(unsigned first, unsigned count) {
//fprintf(stdout, "loadFrames(%u,%u)\n", first, count);  
  if (count > maxCookedFrames) {
    error("ERROR: StreamSource::loadFrames: requested %u frames, but the frame queue can only hold %u\n", count, maxCookedFrames);
  }

  if (numFramesInSegment > 0 && first + count > numFramesInSegment) {
    error("ERROR: StreamSource::loadFrames: requested frame %u, but the segment only has %u frames\n", first+count, numFramesInSegment);
  }

  if (first < firstCookedFrameNum) {
    error("ERROR: StreamSource::loadFrames: requested frame %u which is no longer available; the earliest available frame is %u\n", first, firstCookedFrameNum);
  }

  if (first > firstCookedFrameNum + currentCookedFrames) {
//fprintf(stdout, "  first %u  num %u\n", firstCookedFrameNum, currentCookedFrames);
    error("ERROR: StreamSource::loadFrames: requested frame %u would require skipping %u frames\n", first - firstCookedFrameNum - currentCookedFrames);
  }
  
  if (first + count <= firstCookedFrameNum + currentCookedFrames) {
    // all requested frames are already in the queue
//fprintf(stdout, "StreamSource > [%u,%u)  all frames in queue\n", first,first+count);
    return cookedBuffer + nFeatures * (first - firstCookedFrameNum);
  }

  unsigned numNewFramesNeeded = 
    (first + count) - (firstCookedFrameNum + currentCookedFrames);

  // enqueue the missing frames - may change firstCookedFrameNum etc.
#if 0
	firstCookedFrameNum + currentCookedFrames, 
	firstCookedFrameNum + currentCookedFrames + numNewFramesNeeded, 
	first, first+count);
#endif
  unsigned numAdded = enqueueFrames(numNewFramesNeeded);
  if (numAdded > 0) {
    return cookedBuffer + nFeatures * (first - firstCookedFrameNum);
  } else {
    return NULL;
  }
}



unsigned 
StreamSource::enqueueFrames(unsigned nFrames) {
//fprintf(stdout, "enqueueFrames(%u)\n", nFrames);
  if (numFramesInSegment != 0 && firstCookedFrameNum + currentCookedFrames >= numFramesInSegment)
    return 0; // current segment's done - call preloadFrames() to reset for the next
  if (nFrames > maxCookedFrames) {
    error("ERROR: StreamSource::enqueueFrames: requested %u frames, but the frame queue can only hold %u\n", nFrames, maxCookedFrames);
  }
  if (currentCookedFrames + nFrames > maxCookedFrames) {
    // need to make room
#if 0
    FIXME - check that any required frames are preserved;
#endif
    // FIXME - better strategy for deciding how much to dump?
    unsigned numFramesToDrop = 
      currentCookedFrames + nFrames - maxCookedFrames;

    Data32 *newFirstFrame = cookedBuffer + numFramesToDrop * nFeatures;
    unsigned bytesToMove = 
      (currentCookedFrames - numFramesToDrop) * nFeatures * sizeof(Data32);

//fprintf(stdout, "  moving %u frames, %u -> %u first frame\n", numFramesToDrop, firstCookedFrameNum, firstCookedFrameNum + numFramesToDrop);
    // FIXME - memcpy doesn't permit overlap. Could memcpy numFramesToDrop
    // frames at a time, but would the overhead be worse than a single
    // memmove? Or, just code up an efficient move ...
    memmove(cookedBuffer, newFirstFrame, bytesToMove);
 
    currentCookedFrames -= numFramesToDrop;
    firstCookedFrameNum += numFramesToDrop;
  }
  Data32 *newFrameDest = cookedBuffer + currentCookedFrames * nFeatures;
  unsigned result = cookFrames(newFrameDest, firstCookedFrameNum + currentCookedFrames, nFrames);
  return result;
}


unsigned
StreamSource::cookFrames(Data32 *destination, unsigned first, unsigned count) {
  subMatrixDescriptor *requiredRaw;
  if (filter) {
    // which raw frames should we cook?
    requiredRaw = filter->getRequiredInput(first, count, nFloat, nInt, 0);
    assert(requiredRaw);
    first = requiredRaw->firstFrame;
    count = requiredRaw->numFrames;
  }
  bool eos = false;
  Data32 const *rawFrames = loadRawFrames(first, count, eos);
  Data32 const *cookedFrames;
  if (count > 0) {
    if (filter) {
//fprintf(stdout, "  filtering\n");
      requiredRaw->numFrames = count;
      subMatrixDescriptor output;
      cookedFrames = filter->transform(rawFrames, *requiredRaw, &output);
      first = output.firstFrame;
      count = output.numFrames;
    } else {
//fprintf(stdout,"  no filter\n");
      cookedFrames = rawFrames;
    }
    memcpy(destination, cookedFrames, count * nFeatures * sizeof(Data32));
    currentCookedFrames += count;
  }
  if (eos) {
    numFramesInSegment = first + count;
//fprintf(stdout, "cookFrames: set segment length to %u\n", numFramesInSegment);
  }
  return count;
}


Data32 const *
StreamSource::loadRawFrames(unsigned first, unsigned &count, bool &eos) {
  
  if (count > maxRawFrames) {
    error("ERROR: StreamSource::loadRawFrames: requested %u frames, but the frame queue can only hold %u\n", count, maxRawFrames);
  }

  if (first < firstRawFrameNum) {
    error("ERROR: StreamSource::loadRawFrames: requested frame %u which is no longer available; the earliest available frame is %u\n", first, firstRawFrameNum);
  }

  if (first > firstRawFrameNum + currentRawFrames) {
    error("ERROR: StreamSource::loadRawFrames: requested frame %u would require skipping %u frames\n", first - firstRawFrameNum - currentRawFrames);
  }
  
  if (first + count <= firstRawFrameNum + currentRawFrames) {
    // all requested frames are already in the queue
    return rawBuffer + numRawFeatures * (first - firstRawFrameNum);
  }

  unsigned numNewFramesNeeded = 
    (first + count) - (firstRawFrameNum + currentRawFrames);
  assert(numNewFramesNeeded > 0);
  // enqueue the missing frames - may change firstRawFrameNum etc.
  count = enqueueRawFrames(numNewFramesNeeded, eos);
//if (eos) fprintf(stdout,"loadRawFrames: detected eos\n");
  return rawBuffer + numRawFeatures * (first - firstRawFrameNum);
}


unsigned 
StreamSource::enqueueRawFrames(unsigned nFrames, bool &eos) {
  assert(nFrames > 0);
  if (nFrames > maxRawFrames) {
    error("ERROR: StreamSource::enqueueRawFrames: requested %u frames, but the frame queue can only hold %u\n", nFrames, maxRawFrames);
  }
  if (currentRawFrames + nFrames > maxRawFrames) {
    // need to make room
#if 0
    check that any required frames are preserved;
#endif
    // FIXME - better strategy for deciding how much to dump?
    unsigned numFramesToDrop = 
      currentRawFrames + nFrames - maxRawFrames;

    Data32 *newFirstFrame = rawBuffer + numFramesToDrop * numRawFeatures;
    unsigned bytesToMove = 
      (currentRawFrames - numFramesToDrop) * numRawFeatures * sizeof(Data32);

    // FIXME - memcpy doesn't permit overlap. Could memcpy numFramesToDrop
    // frames at a time, but would the overhead be worse than a single
    // memmove? Or, just code up an efficient move ...  Or alternate between
    // two queue buffers
    memmove(rawBuffer, newFirstFrame, bytesToMove);
//fprintf(stdout,"enqueueRawFrames: dropping %u frames - [%u,%u) ", numFramesToDrop, firstRawFrameNum, firstRawFrameNum + currentRawFrames);
    currentRawFrames -= numFramesToDrop;
    firstRawFrameNum += numFramesToDrop;
//fprintf(stdout,"-> [%u,%u)\n", firstRawFrameNum, firstRawFrameNum + currentRawFrames);
  }

//fprintf(stdout, "enqueueRawFrames : requested %u frames\n", nFrames);
  unsigned rawBuffOffset = currentRawFrames * numRawFeatures;
  eos = false;
  for (unsigned i=0; i < nFrames && !eos; i+=1) {
    for (unsigned j=0; j < nStreams; j+=1) {
      Data32 *newFrameDest = floatStart[j] + rawBuffOffset;
      Data32 const *newFrame = stream[j]->getNextLogicalFrame();
      if (!newFrame && !eos) {
	nFrames = i;
//fprintf(stdout, "enqueueRawFrames: detected end of segement @ %u\n", nFrames);
	eos = true;
	if (j != 0) {
	  error("ERROR: StreamSource::enqueueRawFrames: streams disagree on length of current segment");
	}
      } else if (newFrame && eos) {
	error("ERROR: StreamSource::enqueueRawFrames: streams disagree on length of current segment");
      }
      if (newFrame) {
#if 0
fprintf(stdout, "of%u:", j);
float *fp = (float *)newFrame;
for (unsigned k=0; k < stream[j]->numLogicalContinuous(); k+=1) fprintf(stdout," %f", fp[k]);
unsigned *up = (unsigned *)(newFrame + stream[j]->numLogicalContinuous());
for (unsigned k=0; k < stream[j]->numLogicalDiscrete(); k+=1) fprintf(stdout," %u", up[k]);
fprintf(stdout,"\n");
#endif
	memcpy(newFrameDest, newFrame, stream[j]->numLogicalContinuous() * sizeof(Data32));
	newFrameDest = intStart[j] + rawBuffOffset;
	memcpy(newFrameDest, newFrame+stream[j]->numLogicalContinuous(), stream[j]->numLogicalDiscrete() * sizeof(Data32));
      }
    }
    currentRawFrames += 1;
    rawBuffOffset += numRawFeatures;
  }
  return nFrames;
}

