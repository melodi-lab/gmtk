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
#include "GMTK_FilterStream.h"
#include "GMTK_MergeStream.h"

// The StreamSource subclasses (see GMTK_StreamSource.h) handle
// sequential access data sources like pipes and sockets.
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
  unsigned firstCookedFrameIndex;  // starting position (frame #) of first queued frame in cookedBuffer
  unsigned numFramesInSegment;     // # of frames in current segment
                                   // 0 until we know what it is
  
  int      segmentNum;             // current segment
  
  unsigned nFloat;
  unsigned nInt;
  unsigned nFeatures;

  // low-level stream driver (ASCII, binary, merge, filter)
  ObservationStream *stream;

  unsigned curFrame;

  unsigned _startSkip;

  unsigned _minPastFrames;
  unsigned _minFutureFrames;

 public:

  StreamSource();

  StreamSource(unsigned nStreams, ObservationStream *stream[], unsigned queueLength, 
	       char *filterStr = NULL, unsigned startSkip=0); 

  ~StreamSource() {
    if (cookedBuffer) delete[] cookedBuffer;
    if (stream) delete stream;
  }

  void initialize(unsigned queueLength, unsigned startSkip=0); 

  // Resets queue state for starting a new segment & preloads
  // the requested # of frames

  // side effect: may set numFramesInSegment

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


  void resetFrameNumbers(unsigned firstFrameNumber) {
    infoMsg(IM::ObsStream, IM::Default, "StreamSource::reset [%u,%u] -> [%u,%u]\n",
	    firstCookedFrameNum, firstCookedFrameNum + currentCookedFrames - 1,
	    firstFrameNumber, firstFrameNumber + currentCookedFrames - 1);
    firstCookedFrameNum = firstFrameNumber;
  }

  bool EOS() { 
    assert(stream);
    return stream->EOS();
  }

  
  // returns 0 until the length of the current segment is known
  unsigned numFrames() { return numFramesInSegment; }


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

  // streams don't suppot endSkip
  unsigned endSkip() {return 0;}

  unsigned minPastFrames() {return _minPastFrames;}
  unsigned minFutureFrames() {return _minFutureFrames;}
  void setMinPastFrames(unsigned n) {_minPastFrames = n;}
  void setMinFutureFrames(unsigned n) {_minFutureFrames = n;}

  float *const floatVecAtFrame(unsigned f) {return (float *)loadFrames(f,1);}

  float *const floatVecAtFrame(unsigned f, const unsigned startFeature) {
    assert(startFeature < numContinuous());
    return floatVecAtFrame(f)+startFeature;
  }

  unsigned *const unsignedVecAtFrame(unsigned f) {
    return (unsigned *)loadFrames(f,1) + numContinuous();
  }

  unsigned &unsignedAtFrame(const unsigned frame, const unsigned feature) {
    assert(numContinuous() <= feature && feature < numFeatures()); 
    return ((unsigned *)baseAtFrame(frame))[feature];
  }

  Data32 const * const baseAtFrame(unsigned f) {
    return loadFrames(f,1);
  }

  bool elementIsDiscrete(unsigned el) {
    return numContinuous() <= el && el < numFeatures();
  }

  bool active() { return true; }

  unsigned segmentNumber() { return (unsigned)segmentNum; }

  bool openSegment(unsigned seg) {
    if (seg != (unsigned) segmentNum) {
      error("ERROR: tried to open segment %u while segment %d is active", seg, segmentNum);
    }
    return true;
  }
};

#endif
