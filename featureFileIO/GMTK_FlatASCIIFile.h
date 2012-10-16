
/*
 * GMTK_FlatASCIIFile.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2012, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */

#ifndef GMTK_FLATASCIIFILE_H
#define GMTK_FLATASCIIFILE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <string>
using namespace std;
#include <vector>

#include "machine-dependent.h"
#include "error.h"
#include "general.h"

#include "GMTK_ObservationFile.h"


// Reads in flat ASCII files, which consist of a single file
// with lines in the format:
//
// S F f_0 ... f_n i_0 ... i_m
//
// where S is the segment number, F is the frame number,
// f_i are the continuous features and i_j are the 
// discrete features. The entire file is read into
// memory.

class FlatASCIIFile: public ObservationFile {

  vector<unsigned>  nFrames;       // nFrames[i] is the # of frames in the ith segment
  unsigned          nSegments;
  int               currSegment;

  bool        cppIfAscii;
  char const *cppCommandOptions;

  Data32     *buffer;              // data for current segment
  Data32    **segment;             // the frames for the ith segment start at segment[i]

 public:

  FlatASCIIFile(const char *name, unsigned nfloats, unsigned nints, unsigned num, 
                bool cppIfAscii, char const *cppCommandOptions=NULL,
		char const *contFeatureRangeStr_=NULL, 
		char const *discFeatureRangeStr_=NULL, 
		char const *preFrameRangeStr_=NULL, 
		char const *segRangeStr_=NULL);


  ~FlatASCIIFile() {
    if (buffer) delete [] buffer;
    if (segment) delete [] segment;
  }
 
  // The number of available segments.
  unsigned numSegments() {return nSegments;}

  
  // Set the frame range for the segment
  bool openSegment(unsigned seg) {
    assert(seg < nSegments);
    currSegment = seg;
    if (preFrameRange)
      delete preFrameRange;
    if (preFrameRangeStr) {
      preFrameRange = new Range(preFrameRangeStr,0,nFrames[seg]);
      assert(preFrameRange);
    }
    return true;
  }

  // The number of frames in the currently open segment.
  unsigned numFrames()  {
    assert(currSegment >= 0);
    return nFrames[currSegment];
  }


  // get the starting frame of the current segment and then
  // add the offset for the first requested frame
  Data32 const *getFrames(unsigned first, unsigned count) {
    assert(currSegment > -1); 
    assert(first < nFrames[currSegment]);
    assert(first + count <= nFrames[currSegment]);
    return segment[currSegment] + first * _numFeatures;
  }

  // Number of continuous/discrete/total features in the file
  // after applying -frX and -irX
  unsigned numLogicalContinuous() { return _numLogicalContinuousFeatures; }
  unsigned numLogicalDiscrete()   { return _numLogicalDiscreteFeatures; }
  unsigned numLogicalFeatures()   { return _numLogicalFeatures; }

};

#endif
