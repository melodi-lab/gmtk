
/*
 * GMTK_ASCIIFile.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2011, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */

#ifndef GMTK_ASCIIFILE_H
#define GMTK_ASCIIFILE_H

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


class ASCIIFile: public ObservationFile {

  bool        cppIfAscii;
  char const *cppCommandOptions;
  unsigned    numFileNames;             // number of filenames

  char      *fofName;
  FILE      *fofFile;                   // this file (list of file names)

  unsigned   nFrames;                   // # of frames in current segment

  char     **dataNames;  // pointers to individual filenames

  Data32    *buffer;     // data for current segment
  unsigned   bufferSize; // in Data32s

 public:

  ASCIIFile(const char *name, unsigned nfloats, unsigned nints, 
	    unsigned num, bool cppIfAscii, char const* cppCommandOptions=NULL,
	    char const *contFeatureRangeStr_=NULL, 
	    char const *discFeatureRangeStr_=NULL, 
	    char const *preFrameRangeStr_=NULL, 
	    char const *segRangeStr_=NULL);


  ~ASCIIFile() {
    if (fofName) delete fofName;
    if (dataNames) {
      for (unsigned i = 0; i < numFileNames; i++)
	if (dataNames[i])
	  delete [] dataNames[i];
      delete [] dataNames;
    }
    if (buffer) free(buffer);
  }
 
  // The number of available segments.
  unsigned numSegments() {return numFileNames;}

  // Begin sourcing data from the requested segment.
  // Must be called before any other operations are performed on a segment.
  bool openSegment(unsigned seg);

  // The number of frames in the current segment is set by the openSegment() call.
  unsigned numFrames()  {
    assert(buffer);
    return nFrames;
  }

  // The data for the segment is loaded during the openSegment() call.
  Data32 const *getFrames(unsigned first, unsigned count) {
    assert(buffer);
    assert(first < nFrames);
    assert(first + count <= nFrames);
    return buffer + first * _numFeatures;
  }

  // Number of continuous/discrete/total features in the file
  // after applying -frX and -irX
  unsigned numLogicalContinuous() { return _numLogicalContinuousFeatures; }
  unsigned numLogicalDiscrete()   { return _numLogicalDiscreteFeatures; }
  unsigned numLogicalFeatures()   { return _numLogicalFeatures; }
};

#endif
