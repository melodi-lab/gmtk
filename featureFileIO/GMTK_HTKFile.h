
/*
 * GMTK_HTKFile.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2011, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */

#ifndef GMTK_HTKFILE_H
#define GMTK_HTKFILE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <string>
using namespace std;

#include "machine-dependent.h"
#include "error.h"
#include "general.h"

#include "GMTK_Stream.h"

#include "GMTK_ObservationFile.h"


// Most of the implementation is recycled from the
// previous implementation

class HTKFile: public ObservationFile {

  StreamInfo *info;
  unsigned    nLogicalFrames; // in current segment

  Data32     *buffer;
  unsigned    bufferSize;     // in Data32's

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

  HTKFile(char const *name, unsigned nfloats, unsigned nints, unsigned num, 
	  bool swap, bool cppIfAscii, char const* cppCommandOptions=NULL,
	  char const *contFeatureRangeStr_=NULL, 
	  char const *discFeatureRangeStr_=NULL, 
	  char const *preFrameRangeStr_=NULL, 
	  char const *segRangeStr_=NULL)  
    : ObservationFile(contFeatureRangeStr_, 
		      discFeatureRangeStr_, 
		      preFrameRangeStr_,
		      segRangeStr_)
  {
    listFile = NULL;
    writeFile = NULL;
    fileName = name;
    unsigned format = HTK;
    info = new StreamInfo(name, 
			  contFeatureRangeStr_,
			  discFeatureRangeStr_,
			  &nfloats, &nints,
			  &format, 
			  swap, 
			  num, 
			  cppIfAscii, (char *)cppCommandOptions,
			  segRangeStr_);
    _numContinuousFeatures = nfloats;
    _numDiscreteFeatures   = nints;
    _numFeatures           = nfloats + nints;

    if (contFeatureRangeStr) {
      contFeatureRange = new Range(contFeatureRangeStr_,0,numContinuous());
      assert(contFeatureRange);
    }
    _numLogicalContinuousFeatures = info->getNumFloatsUsed();

    if (discFeatureRangeStr) {
      discFeatureRange = new Range(discFeatureRangeStr_,0,numDiscrete());
      assert(discFeatureRange);
    }
    _numLogicalDiscreteFeatures = info->getNumIntsUsed();

    _numLogicalFeatures = _numLogicalContinuousFeatures + _numLogicalDiscreteFeatures;

    buffer = NULL;
    bufferSize = 0;
  }
  
  HTKFile(char const *listFileName, char const *outputFileName, 
	  char const *outputNameSeparatorStr, bool swap,
	  unsigned nfloats, unsigned nints) 
    : outputFileName(outputFileName), outputNameSeparatorStr(outputNameSeparatorStr), oswap(swap)
  {
    fileName = listFileName;
    if(fileName != NULL) {
      if ((listFile = fopen(fileName, "w")) == NULL) {
	error("Couldn't open output list (%s) for writing.\n", fileName);
      }
    } else {
      error("ERROR: null HTK list file name\n");
    }
    _numContinuousFeatures = nfloats;
    _numDiscreteFeatures = nints;
    _numFeatures = nfloats + nints;
    
    writeFile = NULL;
    currSegment = 0;
    currFrame = 0;
    currFeature = 0;
  }

  ~HTKFile() {
    if (info)   delete info;
    if (buffer) free(buffer);
    if (listFile) {
      if (fclose(listFile)) {
	error("ERROR: failed to close output list file '%s'\n", fileName);
      }
    }
    if (writeFile) {
      if (fclose(writeFile)) {
	error("ERROR: failed to close output file\n");
      }
    }
  }
 
  // Write segment to the file
  void writeSegment(Data32 const *segment, unsigned nFrames);

  // Write frame to the file (call endOfSegment after last frame of a segment)
  void writeFrame(Data32 const *frame);

  // Write the next feature in the current frame (call endOfSegment after last frame of a segment)
  void writeFeature(Data32 x);

  // Call after last writeFrame of a segment
  void endOfSegment();

  // The number of available physical segments.
  unsigned numSegments() {
    assert(info);
    return (unsigned) info->getFullFofSize();
  }

  // The number of available logical segments.
  unsigned numLogicalSegments() {
    assert(info);
    return (unsigned) info->getFofSize();
  }

  // Begin sourcing data from the requested physical segment.
  // Must be called before any other operations are performed on a segment.
  bool openSegment(unsigned seg);

  bool openLogicalSegment(unsigned seg) {
    assert(info);
    return openSegment( info->mapToValueInRange(seg) );
  }

  // The number of frames in the currently open physical segment.
  unsigned numFrames() {
    assert(info && info->curHTKFileInfo);
    return info->getCurNumFrames();
  }

  unsigned numLogicalFrames() {
    assert(info && info->curHTKFileInfo);
    return nLogicalFrames;
  }

  Data32 const *getFrames(unsigned first, unsigned count);

  // Number of continuous/discrete/total features in the file
  // after applying -frX and -irX
  unsigned numLogicalContinuous() { return _numLogicalContinuousFeatures; }
  unsigned numLogicalDiscrete()   { return _numLogicalDiscreteFeatures; }
  unsigned numLogicalFeatures()   { return _numLogicalFeatures; }

};

#endif
