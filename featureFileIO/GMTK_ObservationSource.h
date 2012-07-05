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

// This is the base class for the new observation implementation.
// The inference code should use only this interface to get to the 
// observation data. It manages all the loading, transformation, 
// caching, etc. of the observation data.

// The StreamSource subclass (see GMTK_StreamSource.h) handles
// non-random access data sources like pipes and sockets.
// Online inference assumes a StreamSource.

// The FileSource subclass (see GMTK_FileSource.h) handles
// random access data sources like the various binary file
// formats.

class ObservationSource {

 public:


  // The number of continuous, discrete, total features in the observed data
  virtual unsigned numContinuous() = 0;
  virtual unsigned numDiscrete() = 0;
  virtual unsigned numFeatures() = 0;

  // The number of Data32's between each frame (equals numFeatures())
  virtual unsigned stride() = 0;

  // If frames [i,j] are loaded with loadFrames() or floatVecAtFrame() 
  // (in the latter case, i=j), ensure that [i-minPastFrames, j+minFutureFrames]
  // are present in the frame cache. This is necessary because dlinks,
  // VE CPTs, or filters may require some preceding or following frames
  // in order to compute the requested frames.
  virtual unsigned minPastFrames() = 0;
  virtual unsigned minFutureFrames() = 0;
  virtual void setMinPastFrames(unsigned n) = 0;
  virtual void setMinFutureFrames(unsigned n) = 0;

  // number of frames to skip at the beginning/end of a segement
  virtual unsigned startSkip() = 0;

  // Note that you can't use endSkip() with StreamSource since the segment
  // length may be unkown
  virtual unsigned endSkip() = 0;

  // Returns true iff a segment is open
  virtual bool active() = 0;


  // Specify the segment to use. Note that this call is unnecessary for 
  // StreamSource, which will error() if openSegment() is called with a
  // segment number other than that of the currently active segment - see
  // GMTK_StreamSource.h for more details. For FileSource, a segment must
  // be open before the following methods will work. They all work on
  // the segment specified in the most recent openSegment() call. Returns
  // true iff the segment is successfully openned, generally aborts with
  // error() if something fails.
  virtual bool openSegment(unsigned seg) = 0;

  // Returns the number of the currently active segment (most recent
  // openSegment() argument for FileStream).
  virtual unsigned segmentNumber() = 0;

  // Returns the number of frames in the currently active segment. Note
  // that numFrames() returns 0 for StreamSource until the length of the
  // current segment becomes known.
  virtual unsigned numFrames() = 0;

  // Load count frames of observed data, starting from first.
  // For a StreamSource, it may take arbitrarily long to
  // fullfill the request since the data may be coming from
  // an asynchronous pipe or socket. 
  //
  // Note that StreamSources (see GMTK_StreamSource.h) do not
  // support random access, so no skipping or re-reading any frames.
  //
  virtual Data32 const *loadFrames(unsigned first, unsigned count) = 0;
  

  // Return the vector/feature @ the specified frame. Note that these use
  // loadFrames() internally, and so are subject to the same restrictions
  // discussed above.
  virtual float *const floatVecAtFrame(unsigned f) = 0;
  virtual float *const floatVecAtFrame(unsigned f, const unsigned startFeature) = 0;

  virtual unsigned *const unsignedVecAtFrame(unsigned f) = 0;
  virtual unsigned &unsignedAtFrame(const unsigned frame, const unsigned feature) = 0;

  virtual Data32 const * const baseAtFrame(unsigned f) = 0;

  // Returns true iff the specified element is discrete
  virtual bool elementIsDiscrete(unsigned el) {
    return numContinuous() <= el && el < numFeatures();
  }

};

// Global variable that the inference code uses to access the observation
// data. It will be a pointer to a FileSource instance except for gmtkOnline
// which uses a StreamSource.
extern ObservationSource *globalObservationMatrix;

#endif
