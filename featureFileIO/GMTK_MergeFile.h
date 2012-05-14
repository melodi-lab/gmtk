
/*
 * GMTK_MergeFile.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2012, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */

#ifndef GMTK_MERGEFILE_H
#define GMTK_MERGEFILE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "machine-dependent.h"
#include "range.h"

#include "GMTK_ObservationFile.h"
#include "GMTK_Filter.h"

class MergeFile: public ObservationFile {

  // the files assembled to form the observations
  unsigned nFiles;
  ObservationFile **file;

  unsigned *floatStart; // the ith file's continuous features start here
  unsigned *intStart;   // the ith file's discrete features start here
  
  int       ftrcombo;   // how should the continuous features from the multiple files be combined

  unsigned const *sdiffact;   // how to adjust for files w/ different # of segs
  unsigned const *fdiffact;   //  ... frames
  unsigned _numSegments;      // after considering -sdiffactX
  unsigned _numFrames;        // after considering -fdiffactX
  unsigned _numContinuous;
  unsigned _numDiscrete;

  int      segment;           // currently open segment, -1 if none yet

  Data32  *buffer;            // data for current segment
  unsigned buffSize;          // in Data32's
  unsigned bufStride;         // increment between frames in buffer

  // map requested global segment # to logical segment # in specified file
  unsigned adjustForSdiffact(unsigned fileNum, unsigned seg);

  // map requested merged frame range to individual files' frame range
  void adjustForFdiffact(unsigned first, unsigned count, unsigned fileNum,
			 unsigned &adjFirst, unsigned &adjCount, unsigned &deltaT);
 public:

  MergeFile(unsigned nFiles, ObservationFile *file[], 
	    unsigned const *sdiffact = NULL, 
	    unsigned const *fdiffact = NULL,
	    int ftrcombo=FTROP_NONE);

  virtual ~MergeFile() {
    if (file) delete [] file;
    if (floatStart) delete [] floatStart;
    if (intStart) delete [] intStart;
    if (buffer) free(buffer);
  }

  // The number of segments after -sdiffact
  unsigned numSegments() { return _numSegments; }

  bool openSegment(unsigned seg);

  // The number of frames in the currently open segment after -fdiffact
  unsigned numFrames() {
    assert(segment >= 0);
    return _numFrames;
  }

  // Load count frames of observed data, starting from first (physical),
  // in the current segment. count may be 0 to request loading
  // the entire data segment (frames [first, numFrames)).
  Data32 const *getFrames(unsigned first, unsigned count);

  unsigned numContinuous() { return _numContinuous; }

  unsigned numDiscrete() { return _numDiscrete; }

};

#endif
