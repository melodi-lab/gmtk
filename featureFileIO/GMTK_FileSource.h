/*
 * GMTK_FileSource.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2011, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */

#ifndef GMTK_FILESOURCE_H
#define GMTK_FILESOURCE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "machine-dependent.h"
#include "GMTK_FilterFile.h"
#include "GMTK_ObservationSource.h"
#include "GMTK_Filter.h"

// FileSource handles
// random access data sources like the various binary file
// formats.

class FileSource: public ObservationSource {

 private:

  // buffer to hold transformed cahced/prefected observations
  Data32 *cookedBuffer;
  unsigned bufStride;  // increment between frames in cookedBuffer


  unsigned firstBufferedFrame;        // frame # of first buffered frame
  unsigned firstBufferedFrameIndex;   // index in cooked buffer of where the first buffered frame starts
  unsigned numBufferedFrames;         // # of frames in cookedBuffer
  unsigned bufferFrames;              // size of cookedBuffer in frames
  unsigned bufferSize;                // size of cookedBuffer in Data32

  unsigned window;                    // # of frames to load at once
  unsigned delta;                     // load more frames when within delta frames of the edge of the window

  ObservationFile *file;

  // Justification
  int justificationMode;
  unsigned justificationOffset;

  // number of frames to skip at the beginning
  unsigned _startSkip;
  unsigned _endSkip;
  
  int ftrcombo;

  unsigned _minPastFrames;
  unsigned _minFutureFrames;

  unsigned _numCacheableFrames;   // after considering -fdiffactX  & -gpr only
                                  // This is the # of frames that can be cached

  unsigned _numFrames;            // after -fdiffactX, -gpr, -startSkip, -endSkip, -justification
                                  // This is the # of frames actually accessible to clients

  // _numFrames <= _numCacheableFrames.  _numFrames will be less when -startSkip
  // or -endSkip reserve some frames at the start or end of each segment. These 
  // reserved frames (for dlinks, VECPTs) must be present in the cache when the first
  // or last several frames are accessed, but are not directly accessable to
  // clients via loadFrames(), floatVecAtFrame(), etc.

  int      segment;               // currently open segment; -1 if none yet

  bool     constantSpace;         // if true, load only O(1) frames at a time
                                  // if false, load the entire segment, resizing the buffer 
                                  //           if necessary

  // Load requested frames into cooked buffer, starting at the specified index.
  // index is in frames, so the frames will start at 
  // cookedBuffer[bufferIndex * buffStride]
  Data32 const *loadFrames(unsigned bufferIndex, unsigned first, unsigned count);

 public:

#define DEFAULT_BUFFER_SIZE       (16 * 1024 * 1024)
#ifndef DEFAULT_FILE_WINDOW_SIZE
#define DEFAULT_FILE_WINDOW_SIZE  (4)
#define DEFAULT_FILE_WINDOW_DELTA (100)
#endif
#define DEFAULT_FILE_WINDOW_BYTES (DEFAULT_FILE_WINDOW_SIZE * 1024 * 1024)

  // each file provides its data for the requested submatrix,
  // transformed and "ready to eat." So the FileSource
  // just manages assembling them to satisfy the loadFrames()
  // calls, prefetching/caching. For archipelagos, each thread
  // gets its own FileSource (all aimed at the same files, of course).
  FileSource(ObservationFile *file,
	     unsigned windowBytes = DEFAULT_FILE_WINDOW_BYTES, 
	     unsigned deltaFrames = DEFAULT_FILE_WINDOW_DELTA,
	     unsigned bufferSize  = DEFAULT_BUFFER_SIZE, 
	     unsigned startSkip=0, unsigned endSkip=0,
	     int justificationMode=0, bool constantSpace = false); 

  FileSource() {
    cookedBuffer = NULL;
    bufferSize = 0;
    bufStride = 0;
    window = 0;
    delta = 0;
    numBufferedFrames = 0;
    file = NULL;
    _startSkip = 0;
    _endSkip = 0;
    justificationMode = 0;
    justificationOffset = 0;
    _minPastFrames = 0;
    _minFutureFrames = 0;
    segment = -1;
    constantSpace = false;
  }

  virtual ~FileSource() {
    if (cookedBuffer) delete [] cookedBuffer;
  }

  void initialize(ObservationFile *file,
		  unsigned windowBytes = DEFAULT_FILE_WINDOW_BYTES, 
		  unsigned deltaFrames = DEFAULT_FILE_WINDOW_DELTA,
		  unsigned bufferSize  = DEFAULT_BUFFER_SIZE,
		  unsigned startSkip=0, unsigned endSkip=0,
		  int justificationMode = 0, bool constantSpace = false);

  // The number of available segments.
  unsigned numSegments() { 
    assert(file);
    return file->numSegments(); 
  }

  // The current segment 
  unsigned segmentNumber() {
    assert(segment >= 0);
    return (unsigned) segment;
  }

  // Begin sourcing data from the requested segment.
  // Must be called before any other operations are performed on a segment.
  bool openSegment(unsigned seg);

  // Apply left, center, or right justification to the usable frames
  // in the current segment. numUsableFrames must be <= _numFrames
  // prior to calling justifySegment. After the call, _numFrames = numUsableFrames
  void justifySegment(unsigned numUsableFrames);

  // The number of frames available in the currently open segment.
  unsigned numFrames() {
    assert(segment >= 0);
    return _numFrames;
  }

  Data32 const *loadFrames(unsigned first, unsigned count);
    // if requested frames are already in cookedBuffer
    //   return &cookedBuffer[index of first requested frame]
    // adjust (first,count) for prefetching -> (first',count')
    // adjust cookedBuffer destination
    // getFrames(first', count') from each file into cookedBuffer
    // return &cookedBuffer[destination of first]


  // The number of continuous, discrete, total features
  unsigned numContinuous();
  unsigned numDiscrete();
  unsigned numFeatures();

  // The number of Data32's between each frame
  unsigned stride();

  // number of frames to skip at the beginning
  unsigned startSkip() {return _startSkip;};
  // number of frames to skip at the end
  unsigned endSkip() {return _endSkip;} 

  unsigned minPastFrames() {return _minPastFrames;}
  unsigned minFutureFrames() {return _minFutureFrames;}
  void setMinPastFrames(unsigned n) {_minPastFrames = n;}
  void setMinFutureFrames(unsigned n) {_minFutureFrames = n;}

  float *const floatVecAtFrame(unsigned f);

  float *const floatVecAtFrame(unsigned f, const unsigned startFeature);

  unsigned *const unsignedVecAtFrame(unsigned f);

  unsigned &unsignedAtFrame(const unsigned frame, const unsigned feature);

  Data32 const * const baseAtFrame(unsigned f);

  bool elementIsDiscrete(unsigned el);

  bool active();

};

#endif
