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
			   unsigned startSkip)
  : frameQueueLength(queueLength),
    stream(stream),
    _startSkip(startSkip)
    
{
  assert(stream);
  cookedBuffer = NULL;
  buffSize = 0;
  curFrame = 0;
  if (_startSkip > 0) {
    loadFrames(0, _startSkip);
  }
}


Data32 const *
StreamSource::loadFrames(unsigned first, unsigned count) {
  
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
}

