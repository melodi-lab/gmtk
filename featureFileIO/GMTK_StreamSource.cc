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

StreamSource::StreamSource(ObservationStream *stream, 
			   unsigned queueLength,
			   Filter *filter,
			   unsigned startSkip)
  : cookedBuffSize(queueLength),
    stream(stream),
    filter(filter),
    _startSkip(startSkip)
    
{
  assert(stream);
  nFloat = stream->numContinuous();
  nInt   = stream->numDiscrete();
  nFeatures = nFloat + nInt;

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
  maxRawFrames = rawBuffSize / nFeatures;
  currentRawFrames = 0;
  firstRawFrameNum = 0;
  numFramesInSegment = 0;
}


void
StreamSource::preloadFrames(unsigned nFrames) {
//fprintf(stderr, "preloadFrames(%u)\n", nFrames);
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
//fprintf(stderr, "loadFrames(%u,%u)\n", first, count);  
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
//fprintf(stderr, "  first %u  num %u\n", firstCookedFrameNum, currentCookedFrames);
    error("ERROR: StreamSource::loadFrames: requested frame %u would require skipping %u frames\n", first - firstCookedFrameNum - currentCookedFrames);
  }
  
  if (first + count <= firstCookedFrameNum + currentCookedFrames) {
    // all requested frames are already in the queue
//fprintf(stderr, "  all frames in queue\n");
    return cookedBuffer + nFeatures * (first - firstCookedFrameNum);
  }

  unsigned numNewFramesNeeded = 
    (first + count) - (firstCookedFrameNum + currentCookedFrames);

  // enqueue the missing frames - may change firstCookedFrameNum etc.
  unsigned numAdded = enqueueFrames(numNewFramesNeeded);
  if (numAdded > 0) {
    return cookedBuffer + nFeatures * (first - firstCookedFrameNum);
  } else {
    return NULL;
  }
}



unsigned 
StreamSource::enqueueFrames(unsigned nFrames) {
//fprintf(stderr, "enqueueFrames(%u)\n", nFrames);
  if (nFrames > maxCookedFrames) {
    error("ERROR: StreamSource::enqueueFrames: requested %u frames, but the frame queue can only hold %u\n", nFrames, maxCookedFrames);
  }
  if (currentCookedFrames + nFrames > maxCookedFrames) {
    // need to make room
#if 0
    check that any required frames are preserved;
#endif
    // FIXME - better strategy for deciding how much to dump?
    unsigned numFramesToDrop = 
      currentCookedFrames + nFrames - maxCookedFrames;

    Data32 *newFirstFrame = cookedBuffer + numFramesToDrop * nFeatures;
    unsigned bytesToMove = 
      (currentCookedFrames - numFramesToDrop) * nFeatures * sizeof(Data32);

//fprintf(stderr, "  moving %u frames, %u -> %u first frame\n", numFramesToDrop, firstCookedFrameNum, firstCookedFrameNum + numFramesToDrop);
    // FIXME - memcpy doesn't permit overlap. Could memcpy numFramesToDrop
    // frames at a time, but would the overhead be worse than a single
    // memmove? Or, just code up an efficient move ...
    memmove(cookedBuffer, newFirstFrame, bytesToMove);
 
    currentCookedFrames -= numFramesToDrop;
    firstCookedFrameNum += numFramesToDrop;
  }
  Data32 *newFrameDest = cookedBuffer + currentCookedFrames * nFeatures;
  return cookFrames(newFrameDest, firstCookedFrameNum + currentCookedFrames, nFrames);
}


unsigned
StreamSource::cookFrames(Data32 *destination, unsigned first, unsigned count) {
//fprintf(stderr,"cookFrames(%u,%u)\n", first, count);
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
  if (filter) {
//fprintf(stderr, "  filtering\n");
    subMatrixDescriptor output;
    cookedFrames = filter->transform(rawFrames, *requiredRaw, &output);
    first = output.firstFrame;
    count = output.numFrames;
  } else {
//fprintf(stderr,"  no filter\n");
    cookedFrames = rawFrames;
  }
  if (eos) numFramesInSegment = first + count;
  memcpy(destination, cookedFrames, count * nFeatures * sizeof(Data32));
  currentCookedFrames += count;
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
    return rawBuffer + stream->numLogicalFeatures() * (first - firstRawFrameNum);
  }

  unsigned numNewFramesNeeded = 
    (first + count) - (firstRawFrameNum + currentRawFrames);

  // enqueue the missing frames - may change firstRawFrameNum etc.
  count = enqueueRawFrames(numNewFramesNeeded, eos);
  return rawBuffer + nFeatures * (first - firstRawFrameNum);
}


unsigned 
StreamSource::enqueueRawFrames(unsigned nFrames, bool &eos) {
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

    Data32 *newFirstFrame = rawBuffer + numFramesToDrop * stream->numLogicalFeatures();
    unsigned bytesToMove = 
      (currentCookedFrames - numFramesToDrop) * stream->numLogicalFeatures() * sizeof(Data32);

    // FIXME - memcpy doesn't permit overlap. Could memcpy numFramesToDrop
    // frames at a time, but would the overhead be worse than a single
    // memmove? Or, just code up an efficient move ...
    memmove(rawBuffer, newFirstFrame, bytesToMove);
    currentRawFrames -= numFramesToDrop;
    firstRawFrameNum += numFramesToDrop;
  }
  Data32 *newFrameDest = 
    rawBuffer + currentRawFrames * stream->numLogicalFeatures();
  for (unsigned i=0; i < nFrames; i+=1, newFrameDest += stream->numLogicalFeatures()) {
    Data32 const *newFrame = stream->getNextLogicalFrame();
    if (!newFrame) {
      nFrames = i;
      eos = true;
      break;
    }
    memcpy(newFrameDest, newFrame, stream->numLogicalFeatures() * sizeof(Data32));
    currentRawFrames += 1;
  }
  return nFrames;
}





#if 0
  
  
  if (EOS()) {
    error("ERROR: StreamSource::loadFrames: end-of-stream reached\n");
  }


  for ( ; curFrame < first; curFrame+=1) 
    if (! stream->getNextFrame() ) {
      error("ERROR: StreamSource::loadFrames: requested frames %u to %u, but only $u frames appear to be available in the segment\n", first, first+count-1, curFrame);
    }
  unsigned needed = count * numFeatures();
  if (buffSize < needed) {
    cookedBuffer = (Data32 *) realloc(cookedBuffer, needed * sizeof(Data32));
    assert(cookedBuffer);
    buffSize = needed;
  }

  Data32 *dst = cookedBuffer;
  for ( ; curFrame < first+count; curFrame+=1) {
    Data32 const *src = stream->getNextFrame();
    if (! src) {
      // FIXME - probably shouldn't warn here
      warning("StreamSource::loadFrames: requested frames %u to %u, but only $u frames appear to be available in the segment\n", first, first+count-1, curFrame);
      curFrame = 0; // new segment
      if (_startSkip > 0) {
	loadFrames(0, _startSkip);
      }     
    }
    memcpy(dst, src, needed * sizeof(Data32));
    dst += numFeatures();
  }
  return cookedBuffer;
    // The current design loops over observation segments,
    // loading them into the ObservationMatrix, then inference
    // iterates over the modified partitions of the current
    // segment. The plan is to instead have the inference code
    // call loadFrames() with the frame range needed for the
    // current modified partition (with prefetching, caching,
    // data transformations, etc. happening behind the scenes
    // in StreamSource). [ Note that for online inference, we
    // can just save the inference output directly into the
    // Viterbi unpacking buffers and printing the completed
    // original frames without actually packing & unpacking ]

    // But there may be some frame overlap between consecutive
    // modified partitions. In the new Viterbi printing case, 
    // a single RV instance was used to store the shared data
    // for any modified partitions that contained it. That was
    // possible because the Viterbi printing algorithm worked
    // with sets of *RV, so the pointers in each partition's set 
    // could point to the shared RV instance. For inference, 
    // adding another layer of indirection around the observed
    // Data32s in order to share the Data32 instances between
    // modified partitions would not be good for performance.

    // Fortunately, that only poses a problem when the cookedBuffer
    // fills - the last several frames may still be needed for
    // upcoming modified partitions. So we can copy the needed 
    // frames to the beginning of the cookedBuffer and then load
    // in the new needed frames following them.

    
    // if @ end of cookedBuffer
    //   copy any need frames to beginning of cookedBuffer
    //   adjust cookedBuffer destination
    //   adjust (first,count)to account for overlap -> (first',count')
    // until count' frames read, EOF, or timeout:
    //   getNextFrame() from each stream into cookedBuffer
    // return &cookedBuffer + offset

#endif
