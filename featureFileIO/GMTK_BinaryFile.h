
/*
 * GMTK_BinaryFile.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2012, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */

#ifndef GMTK_BINARYFILE_H
#define GMTK_BINARYFILE_H

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

#ifdef HAVE_FSEEKO
#  define bin_fseek(a,b,c) fseeko(a,b,c)
#  define bin_ftell(a)     ftello(a)
   typedef off_t bin_off_t;
#else
#  define bin_fseek(a,b,c) fseek(a,b,c)
#  define bin_ftell(a)     ftell(a)
   typedef long bin_off_t;
#endif

class BinaryFile: public ObservationFile {

  bool        swap;
  bool        cppIfAscii;
  char const *cppCommandOptions;
  unsigned    numFileNames;             // number of filenames

  char      *fofName;
  FILE      *fofFile;                   // this file (list of file names)
  FILE      *curDataFile;               // currently open segment

  unsigned   nFrames;                   // # frames in current segment before -preprX

  unsigned   curSegment;                // current segment #

  char     **dataNames;  // pointers to individual filenames

  Data32    *buffer;     // data for current segment
  unsigned   buffSize;   // in Data32's
 public:

  BinaryFile(const char *name, unsigned nfloats, unsigned nints, 
	     unsigned num, bool swap, bool cppIfAscii, 
	     char const* cppCommandOptions=NULL,
	     char const *contFeatureRangeStr_=NULL, 
	     char const *discFeatureRangeStr_=NULL, 
	     char const *preFrameRangeStr_=NULL, 
	     char const *segRangeStr_=NULL);


  ~BinaryFile() {
    if (fofName) delete fofName;
    if (dataNames) {
      for (unsigned i = 0; i < numFileNames; i++)
	if (dataNames[i])
	  delete[] dataNames[i];
      delete[] dataNames;
    }
    if (buffer) delete [] buffer;
    if (curDataFile)
      if (fclose(curDataFile))
	warning("BinaryFile: Error closing data file\n");
  }
 
  // The number of available segments.
  unsigned numSegments() {return numFileNames;}

  // Begin sourcing data from the requested segment.
  // Must be called before any other operations are performed on a segment.
  bool openSegment(unsigned seg);

  // The number of frames in the currently open segment.
  unsigned numFrames()  {
    assert(curDataFile);
    return nFrames;
  }

  Data32 const *getFrames(unsigned first, unsigned count);

  // Number of continuous/discrete/total features in the file
  // after applying -frX and -irX
  unsigned numLogicalContinuous() { return _numLogicalContinuousFeatures; }
  unsigned numLogicalDiscrete()   { return _numLogicalDiscreteFeatures; }
  unsigned numLogicalFeatures()   { return _numLogicalFeatures; }
};

#endif
