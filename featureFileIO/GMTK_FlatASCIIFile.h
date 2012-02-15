
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


// The ObservationFile provides a simple API to wrap around
// random access data file formats.
// Just subclass ObservationFile and implement getFrames
// for the new type of file (and add command line options 
// to instantiate it - aspect-oriented programming?) and the
// new file type is supported by GMTK.
//
// Planned subtypes:
//   ASCIIFile    -   ASCII files (read entirely into memory)
//   PFileFile    -   indexed PFiles (non-indexed read entirely into memory)
//   HDF5File     
//   HTKFile      
//   BinaryFile   
//   FilterFile   -   ObservationFile wrapper for IIR, ARMA, etc transforms

class FlatASCIIFile: public ObservationFile {

  unsigned          nFloats;
  unsigned          nInts;
  vector<unsigned>  nFrames;
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

  // Begin sourcing data from the requested segment.
  // Must be called before any other operations are performed on a segment.
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

  // Load count frames of observed data, starting from first,
  // in the current segment. count may be 0 to request loading
  // the entire data segment (frames [first, numFrames)).
  Data32 const *getFrames(unsigned first, unsigned count) {
    assert(currSegment > -1); 
    assert(first < nFrames[currSegment]);
    assert(first + count <= nFrames[currSegment]);
    return segment[currSegment] + first * (nFloats + nInts);
  }

  unsigned numContinuous() {return nFloats;}
  unsigned numDiscrete() {return nInts;}
  unsigned numFeatures() {return nFloats + nInts;}

};

#endif
