/*
 * GMTK_StreamSource.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2011, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */

#ifndef GMTK_STREAMSOURCE_H
#define GMTK_STREAMSOURCE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "machine-dependent.h"
#include "GMTK_ObservationSource.h"
#include "GMTK_ObservationStream.h"

// The StreamSource subclasses (see GMTK_StreamSource.h) handle
// non-random access data sources like pipes and sockets.
// Online inference assumes a StreamSource.

class StreamSource : ObservationSource {

 private:

  // buffer to hold transfored frames - enough to do inference
  // for the current modified partition (possilby includes 
  // some "pseudo-future")
  Data32 cookedBuffer[];

  // the streams assembled to form the observations
  ObservationStream stream[];

 public:
  
  // each stream provides its data for the next frame,
  // transformed and "read to eat." So the StreamSource
  // just manages assembling them to satisfy the loadFrames()
  // calls
  StreamSource(ObservationStream stream[]); 

  Data32 *loadFrames(unsigned first, unsigned count) {
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

    return NULL;
    
    // if @ end of cookedBuffer
    //   copy any need frames to beginning of cookedBuffer
    //   adjust cookedBuffer destination
    //   decrement count by amount of overlap
    // until count frames read, EOF, or timeout:
    //   getNextFrame() from each stream into cookedBuffer
    // return &cookedBuffer + offset

  }
}

#endif
