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
#include "error.h"
#include "GMTK_ObservationSource.h"
#include "GMTK_ObservationStream.h"
#include "GMTK_Filter.h"

// The StreamSource subclasses (see GMTK_StreamSource.h) handle
// non-random access data sources like pipes and sockets.
// Online inference assumes a StreamSource.

class StreamSource : public ObservationSource {

 private:

  // buffer to hold transfored frames - enough to do inference
  // for the current modified partition (possilby includes 
  // some "pseudo-future")
  Data32 *cookedBuffer;
  unsigned cookedBuffSize;
  
  unsigned maxCookedFrames;        // how many frames fit in the queue?
  unsigned currentCookedFrames;    // how many frames are in it now?
  unsigned firstCookedFrameNum;    // frame # of first frame in queue

  unsigned numFramesInSegment;     // # of frames in current segment
                                   // 0 until we know what it is
  Data32 *rawBuffer;
  unsigned rawBuffSize;

  unsigned maxRawFrames;
  unsigned rawFloats;              // # of continuous features in raw frames
  unsigned rawInts;                // # of discrete features in raw frames
  unsigned numRawFeatures;         // sum of the above
  unsigned currentRawFrames;
  unsigned firstRawFrameNum;

  unsigned nFloat;
  unsigned nInt;
  unsigned nFeatures;

  // low-level stream driver (ASCII, binary)
  unsigned nStreams; 
  ObservationStream **stream;

  Data32 **floatStart;
  Data32 **intStart;

  // transform stack to cook the raw frames with
  Filter *filter;

  unsigned curFrame;

  unsigned _startSkip;
  // FIXME - add endSkip, probably to base class?

  // Try to append count cooked frames starting from cooked frame
  // # first to destination. 
  // Returns # of frames actually appended. May be short if 
  // the segment or stream ends.
  unsigned cookFrames(Data32 *destination, unsigned first, unsigned count);

  // Get count raw frames from the buffer, starting at first,
  // calling stream->getNextFrame() as needed. 
  // count will contain the actual # of frames retrieved, and
  // eos is true iff the last frame returned is the last in the stream
  Data32 const *loadRawFrames(unsigned first, unsigned &count, bool &eos);

  unsigned enqueueRawFrames(unsigned nFrames, bool &eos);

 public:

  StreamSource();

  StreamSource(unsigned nStreams, ObservationStream *stream[], unsigned queueLength, 
	       Filter *filter = NULL, unsigned startSkip=0); 

  ~StreamSource() {
    if (cookedBuffer) delete[] cookedBuffer;
    if (rawBuffer) delete[] rawBuffer;
    if (stream) delete stream;
  }

  void initialize(unsigned nStreams, ObservationStream *stream[], unsigned queueLength, 
		  Filter *filter = NULL, unsigned startSkip=0); 

  // Resets queue state for starting a new segment & preloads
  // the requested # of frames

  // side effect: may set numFramesInSegment

  // FIXME - error conditions?
  void preloadFrames(unsigned nFrames);

  // Add the requested # of frames to the queue, possibly 
  // triggering a flush of the older frames in the queue to
  // make room. Returns number of frames actually added.
  
  // side effect: may set numFramesInSegment
  unsigned enqueueFrames(unsigned nFrames);

  // After discussing with Jeff, we decided online inference
  // should support multiple segments, so we need a way to
  // indicate the end of a segment - prehaps by returning 
  // NULL or adding a bool& EOS parameter. Also, there may
  // be (M,S) incompatibility between the model triangulation
  // and the amount of data supplied by the stream  - should
  // that be an error, and if so, how to indicate it? Perhaps
  // make count an unsigned& returning the # of missing frames
  // (should be 0 for success)

  // side effect: may set numFramesInSegment
  Data32 const *loadFrames(unsigned first, unsigned count);


  bool EOS() { 
    assert(stream);
    bool eos = stream[0]->EOS();
    for (unsigned i=1; i < nStreams; i+=1) {
      if (eos != stream[i]->EOS()) {
	// FIXME - maybe just warn and take min # segments?
        error("ERROR: StreamSource: streams disagree on number of segments");
      }
    }
    return eos;
  }

  
  // returns 0 until the length of the current segment is known
  unsigned segmentLength() { return numFramesInSegment; }


  // The number of continuous, discrete, total features

  unsigned numContinuous() {
    assert(stream);
    return nFloat;
  }


  unsigned numDiscrete() {
    assert(stream);
    return nInt;
  }
  

  unsigned numFeatures() {
    assert(stream);
    return nFeatures;
  }


  unsigned stride() {
    assert(stream);
    return nFeatures;
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
