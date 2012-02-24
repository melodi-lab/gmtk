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

#include <assert.h>

#include "machine-dependent.h"
#include "GMTK_ObservationSource.h"
#include "GMTK_ObservationStream.h"

// The StreamSource subclasses (see GMTK_StreamSource.h) handle
// non-random access data sources like pipes and sockets.
// Online inference assumes a StreamSource.

class StreamSource : public ObservationSource {

 private:

  // buffer to hold transfored frames - enough to do inference
  // for the current modified partition (possilby includes 
  // some "pseudo-future")
  Data32 *cookedBuffer;
  unsigned buffSize;
  unsigned frameQueueLength;

  // the streams assembled to form the observations
  ObservationStream *stream;

  unsigned curFrame;

  unsigned _startSkip;

 public:
  
  StreamSource(ObservationStream *stream, unsigned queueLength, unsigned startSkip=0); 


  // After discussing with Jeff, we decided online inference
  // should support multiple segments, so we need a way to
  // indicate the end of a segment - prehaps by returning 
  // NULL or adding a bool& EOS parameter. Also, there may
  // be (M,S) incompatibility between the model triangulation
  // and the amount of data supplied by the stream  - should
  // that be an error, and if so, how to indicate it? Perhaps
  // make count an unsigned& returning the # of missing frames
  // (should be 0 for success)
  Data32 const *loadFrames(unsigned first, unsigned count);
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


  bool EOS() { 
    assert(stream); 
    return stream->EOS();
  }

  
  // The number of continuous, discrete, total features

  unsigned numContinuous() {
    assert(stream);
    return stream->numContinuous();
  }


  unsigned numDiscrete() {
    assert(stream);
    return stream->numDiscrete();
  }
  

  unsigned numFeatures() {
    assert(stream);
    return stream->numFeatures();
  }


  unsigned stride() {
    assert(stream);
    return numFeatures();
  }


  // number of frames to skip at the beginning
  unsigned startSkip() {return _startSkip;};

  float *const floatVecAtFrame(unsigned f) {return NULL;}

  float *const floatVecAtFrame(unsigned f, const unsigned startFeature) {
    return NULL;
  }

  unsigned *const unsignedVecAtFrame(unsigned f) {
    return NULL;
  }

  unsigned &unsignedAtFrame(const unsigned frame, const unsigned feature) {
    return *(new unsigned);
  }

  Data32 const * const baseAtFrame(unsigned f) {
    return NULL;
  }

  bool elementIsDiscrete(unsigned el) {
    return numContinuous() <= el && el < numFeatures();
  }

};

#endif
