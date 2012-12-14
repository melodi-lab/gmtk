
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

  // for writable files
  FILE       *writeFile;      // current segment file
  FILE       *listFile;       // list of files
  char const *outputFileName;
  char const *outputNameSeparatorStr;
  bool        oswap;

  unsigned    currSegment;
  unsigned    currFrame;
  unsigned    currFeature;

 public:

  BinaryFile(const char *name, unsigned nfloats, unsigned nints, 
	     unsigned num, bool swap, bool cppIfAscii, 
	     char const* cppCommandOptions=NULL,
	     char const *contFeatureRangeStr_=NULL, 
	     char const *discFeatureRangeStr_=NULL, 
	     char const *preFrameRangeStr_=NULL, 
	     char const *segRangeStr_=NULL);

  
  BinaryFile(char const *listFileName, char const *outputFileName, 
	     char const *outputNameSeparatorStr, bool swap, unsigned nfloats, unsigned nints) 
    : outputFileName(outputFileName), outputNameSeparatorStr(outputNameSeparatorStr), oswap(swap)
  {
    fileName = listFileName;
    if(fileName != NULL) {
      if ((listFile = fopen(fileName, "w")) == NULL) {
       error("Couldn't open output list (%s) for writing.\n", fileName);
      }
    } else {
      error("ERROR: null ASCII list file name\n");
    }
    _numContinuousFeatures = nfloats;
    _numDiscreteFeatures = nints;
    _numFeatures = nfloats + nints;
    
    writeFile = NULL;
    currSegment = 0;
    currFrame = 0;
    currFeature = 0;
  }


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
    if (listFile) {
      if (fclose(listFile)) {
       error("Error closing output file '%s'\n", fileName);
      }
    }
    if (writeFile) {
      if (fclose(writeFile)) {
       error("Error closing output file\n");
      }
    }
  }
 

  // Write frame to the file (call endOfSegment after last frame of a segment)
  void writeFrame(Data32 const *frame);

  // Write frame to the file (call endOfSegment after last frame of a segment)
  void writeFeature(Data32 x);

  // Call after last writeFrame of a segment
  void endOfSegment();


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
