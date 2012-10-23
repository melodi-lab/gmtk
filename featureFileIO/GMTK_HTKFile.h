
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
  
  ~HTKFile() {
    if (info)   delete info;
    if (buffer) free(buffer);
  }
 
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
