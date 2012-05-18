/*
 * GMTK_ObservationSource.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2011, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */

#ifndef GMTK_OBSERVATIONSOURCE_H
#define GMTK_OBSERVATIONSOURCE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "machine-dependent.h"

// This is the base class for the island-ified and online
// ObservationMatrix. The inference code should only use
// this interface to get to the observation data. It manages
// all the loading, transformation, and caching of the data.

// The StreamSource subclasses (see GMTK_StreamSource.h) handle
// non-random access data sources like pipes and sockets.
// Online inference assumes a StreamSource.

// The FileSource subclasses (see GMTK_FileSource.h) handle
// random access data sources like the various binary file
// formats.
class ObservationSource {

 public:
  // Load count frames of observed data, starting from first.
  // For a StreamSource, it may take arbitrarily long to
  // fullfill the request since the data may be coming from
  // an asynchronous pipe or socket. For finite data sets of
  // sufficiently small size, count may be 0 to request loading
  // the entire data set (frames [first, numFrames)).
  //
  // Note that StreamSources (see GMTK_StreamSource.h) do not
  // support random access, so:
  //   the first first must be >= 0
  //   subsequent firsts must be >= previous(first+count)
  //   perhaps require no gaps ?
  // 
  // error if first+count >= numFrames(), or can it return short counts?
  virtual Data32 const *loadFrames(unsigned first, unsigned count) = 0;
  
  // The number of continuous, discrete, total features
  virtual unsigned numContinuous() = 0;
  virtual unsigned numDiscrete() = 0;
  virtual unsigned numFeatures() = 0;

  // The number of Data32's between each frame
  virtual unsigned stride() = 0;

  // If frames [i,j] are loaded with loadFrames() or floatVecAtFrame() 
  // (in the latter case, i=j), ensure that [i-minPastFrames, j+minFutureFrames]
  // are present in the frame cache
  virtual unsigned minPastFrames() = 0;
  virtual unsigned minFutureFrames() = 0;
  virtual void setMinPastFrames(unsigned n) = 0;
  virtual void setMinFutureFrames(unsigned n) = 0;

  // number of frames to skip at the beginning
  virtual unsigned startSkip() = 0;

  // can't use endSkip() in streams since length is unkown?
  virtual unsigned endSkip() = 0;

  virtual float *const floatVecAtFrame(unsigned f) = 0;

  virtual float *const floatVecAtFrame(unsigned f, const unsigned startFeature) = 0;

  virtual unsigned *const unsignedVecAtFrame(unsigned f) = 0;

  virtual unsigned &unsignedAtFrame(const unsigned frame, const unsigned feature) = 0;

  virtual Data32 const * const baseAtFrame(unsigned f) = 0;

  virtual bool elementIsDiscrete(unsigned el) = 0;

  virtual bool active() = 0;

  virtual bool openSegment(unsigned seg) = 0;

  virtual unsigned segmentNumber() = 0;

  virtual unsigned numFrames() = 0;
};

extern ObservationSource *globalObservationMatrix;

#endif
