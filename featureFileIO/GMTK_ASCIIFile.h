
/*
 * GMTK_ASCIIFile.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (C) 2011 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * 
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

  ASCIIFile(const char *name, unsigned nfloats, unsigned nints, 
	    unsigned num, bool cppIfAscii, char const* cppCommandOptions=NULL,
	    char const *contFeatureRangeStr_=NULL, 
	    char const *discFeatureRangeStr_=NULL, 
	    char const *preFrameRangeStr_=NULL, 
	    char const *segRangeStr_=NULL, 
	    unsigned leftPad=0, unsigned rightPad=0);

  
  ASCIIFile(char const *listFileName, char const *outputFileName, 
	    char const *outputNameSeparatorStr, unsigned nfloats, unsigned nints) 
    : outputFileName(outputFileName), outputNameSeparatorStr(outputNameSeparatorStr)
  {
    cppIfAscii = false;
    cppCommandOptions = NULL;
    fofName = NULL;
    fofFile = NULL;
    dataNames = NULL;
    buffer = NULL;
    
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


  ~ASCIIFile() {
    if (fofName) delete fofName;
    if (dataNames) {
      for (unsigned i = 0; i < numFileNames; i++)
	if (dataNames[i])
	  delete [] dataNames[i];
      delete [] dataNames;
    }
    if (buffer) free(buffer);
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
 
  // Write segment to the file (no need to call endOfSegment)
  void writeSegment(Data32 const *segment, unsigned nFrames);

  // Set frame # to write within current segemnt (partially supported)
  void setFrame(unsigned frame) {
    if (frame == currFrame) return;
    assert(currFeature == 0);
    if (frame > currFrame) {
      for (unsigned i=0; i < frame - currFrame; i+=1) {
	for (unsigned j=0; j < _numFeatures; j+=1) {
	  writeFeature(0);
	}
      }
    } else { // frame < currFrame - can't go backwards
      error("ERROR: ASCII files do not support random access\n");
    }
  }

  // Write frame to the file (call endOfSegment after last frame of a segment)
  void writeFrame(Data32 const *frame);

  // Write the next feature in the current frame (call endOfSegment after last frame of a segment)
  void writeFeature(Data32 x);

  // Call after last writeFrame of a segment
  void endOfSegment();

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
