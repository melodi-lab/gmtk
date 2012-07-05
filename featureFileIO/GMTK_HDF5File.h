
/*
 * GMTK_HDF5File.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2012, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */

#ifndef GMTK_HDF5FILE_H
#define GMTK_HDF5FILE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <string>
using namespace std;

#include "machine-dependent.h"
#include "error.h"
#include "general.h"

#include "GMTK_ObservationFile.h"

#if HAVE_LIBHDF5_CPP

// If Autoconf couldn't find the library, disable the code 
// that requires it.


// The argument to -ofX is the name of a file that contains
// a list of:
//
// filename:groupname:x_start,y_start;x_stride,y_stride;x_count,y_count
//
// Each entry in the file represents a segment. 
// See https://j.ee.washington.edu/trac/gmtk/ticket/21 for 
// additional details and discussion

#include "H5Cpp.h"
#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

class HDF5File: public ObservationFile {

  bool         cppIfAscii;  // for the FoF  
  char const  *cppCommandOptions;

  unsigned     numSlabs;   // # hyperslabs specified in FoF
  unsigned     curSegment; // currently open segment

  unsigned     nFrames;    // # of frames in current segment

  Data32      *buffer;     
  unsigned     bufSize;    // in Data32's

  H5File      *curFile;    // file for current segment

  char const **fileName;
  char const **groupName;  // groupName[i] is the group to read in fileName[i]
  unsigned    *xStart;
  unsigned    *yStart;     // hyperslab in fileName[i]:groupName[i] starts at (xStart[i],yStart[1])
  unsigned    *xStride;
  unsigned    *yStride;    // stride in ith hyperslab
  unsigned    *xCount;
  unsigned    *yCount;     // size of ith hyperslab

  // ith segment is the 2D hyperslab 
  // fileName[i]:groupName[i]:xStart[i],yStart[i];xStride[i],yStride[i];xCount[i],yCount[i]

 public:

  HDF5File(const char *name, unsigned num, 
	   bool cppIfAscii, char const* cppCommandOptions=NULL,
	   char const *contFeatureRangeStr_=NULL, 
	   char const *discFeatureRangeStr_=NULL, 
	   char const *preFrameRangeStr_=NULL, 
	   char const *segRangeStr_=NULL);

  ~HDF5File();
 
  // The number of available segments.
  unsigned numSegments() { return numSlabs; }

  bool openSegment(unsigned seg);

  // The number of frames in the currently open segment.
  unsigned numFrames() {
    assert(curFile); // segment open
    return nFrames; 
  }

  Data32 const *getFrames(unsigned first, unsigned count);

  // Number of continuous/discrete/total features in the file
  // after applying -frX and -irX
  unsigned numLogicalContinuous() { return _numLogicalContinuousFeatures; }
  unsigned numLogicalDiscrete()   { return _numLogicalDiscreteFeatures; }
  unsigned numLogicalFeatures()   { return _numLogicalFeatures; }
};

#else

// Fake do-nothing implementation in case the
// HDF5 library is unavailable

class HDF5File: public ObservationFile {

 public:

  HDF5File(const char *name, unsigned num,
	   bool cppIfAscii, char const* cppCommandOptions=NULL,
	   char const *contFeatureRangeStr_=NULL, 
	   char const *discFeatureRangeStr_=NULL, 
	   char const *preFrameRangeStr_=NULL, 
	   char const *segRangeStr_=NULL)
  {
    error("This GMTK build does not support HDF5 files");
  }

  ~HDF5File() {}
 
  // The number of available segments.
  unsigned numSegments() { return 0; }

  // Begin sourcing data from the requested segment.
  // Must be called before any other operations are performed on a segment.
  bool openSegment(unsigned seg) { return false; }

  // The number of frames in the currently open segment.
  unsigned numFrames()  { return 0; }

  // Load count frames of observed data, starting from first,
  // in the current segment. count may be 0 to request loading
  // the entire data segment (frames [first, numFrames)).
  Data32 const *getFrames(unsigned first, unsigned count) { return NULL; }

  unsigned numContinuous() { return 0; }

  unsigned numDiscrete() { return 0; }

  unsigned numFeatures() {return 0; }

};

#endif
#endif
