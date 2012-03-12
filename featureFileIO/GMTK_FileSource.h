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

  Filter *filter;      // -posttrans filters

  // buffer to hold transformed cahced/prefected observations
  Data32 *cookedBuffer;

  Data32 **floatStart; // the ith file's continuous features start here
  Data32 **intStart;   // the ith file's discrete features start here
  unsigned bufStride;  // increment between frames in cookedBuffer

  // the files assembled to form the observations
  unsigned nFiles;
  ObservationFile **file;

  // Global Per-segment final frame Range
  char const *globalFrameRangeStr; // -gpr
  Range *globalFrameRange;

  // number of frames to skip at the beginning
  unsigned _startSkip;
  unsigned _endSkip;

  unsigned const *sdiffact;   // how to adjust for files w/ different # of segs
  unsigned const *fdiffact;   //  ... frames
  unsigned _numSegments;      // after considering -sdiffactX
  int segment;

  unsigned _numInputFrames;   // after considering -fdiffactX only
  unsigned _numFrames;        // after -fdiffactX and -gpr

  int _numContinuous;         // # of continuous & discrete features after -posttrans
  int _numDiscrete;           //   these are computed once if -1 then cached

  unsigned adjustForSdiffact(unsigned fileNum, unsigned seg);

 public:
  
  // each file provides its data for the requested submatrix,
  // transformed and "ready to eat." So the FileSource
  // just manages assembling them to satisfy the loadFrames()
  // calls, prefetching/caching. For archipelagos, each thread
  // gets its own FileSource (all aimed at the same files, of course).
  FileSource(unsigned _nFiles, ObservationFile *file[], 
	     char const *_globalFrameRangeStr = NULL, 
	     unsigned const *sdiffact = NULL, 
	     unsigned const *fdiffact = NULL,
	     unsigned startSkip=0, unsigned endSkip=0,
	     Filter *posttrans = NULL); 

  FileSource() {
    cookedBuffer = NULL;
    floatStart = NULL;
    intStart = NULL;
    file = NULL;
    globalFrameRangeStr = NULL;
    globalFrameRange = NULL;
    _startSkip = 0;
    _endSkip = 0;
    segment = -1;
    filter = NULL;
  }

  // FIXME - dtor

  void initialize(unsigned nFiles, ObservationFile *file[], 
		  unsigned const *sdiffact = NULL, 
		  unsigned const *fdiffact = NULL,
		  char const *globalFrameRangeStr = NULL, 
		  unsigned startSkip=0, unsigned endSkip=0,
		  Filter *posttrans = NULL);

  // The number of available segments.
  unsigned numSegments() { return _numSegments; }

  // The current segment 
  unsigned segmentNumber() {
    assert(segment >= 0);
    return (unsigned) segment;
  }

  // Begin sourcing data from the requested segment.
  // Must be called before any other operations are performed on a segment.
  bool openSegment(unsigned seg);

  // The number of frames in the currently open segment.
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

  float *const floatVecAtFrame(unsigned f);

  float *const floatVecAtFrame(unsigned f, const unsigned startFeature);

  unsigned *const unsignedVecAtFrame(unsigned f);

  unsigned &unsignedAtFrame(const unsigned frame, const unsigned feature);

  Data32 const * const baseAtFrame(unsigned f);

  bool elementIsDiscrete(unsigned el);

  bool active();

};

extern FileSource globalObservationMatrix;

#endif
