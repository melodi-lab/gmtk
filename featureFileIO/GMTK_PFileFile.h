
/*
 * GMTK_PFileFile.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2012, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */

#ifndef GMTK_PFILEFILE_H
#define GMTK_PFILEFILE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <string>
using namespace std;

#include "machine-dependent.h"
#include "error.h"
#include "general.h"
#include "pfile.h"

#include "GMTK_ObservationFile.h"

class PFileFile: public ObservationFile {

  Data32    *buffer;     // data for current segment
  unsigned   bufferSize; // in Data32's

  // These are the buffers that the low-level pfile code
  // reads the data into. These buffers are then copied
  // into buffer with the continuous and discrete features
  // merged together as GMTK requires.
  float     *contBuf;
  UInt32    *discBuf;

  unsigned   currentSegment;
  char      *fileName;
  FILE      *dataFile;
  InFtrLabStream_PFile *pfile;

  unsigned  _numFrames; // # physical frames in current segment

 public:

  PFileFile(const char *name, unsigned nfloats, unsigned nints,
	    unsigned num, bool bswap, 
	    char const *contFeatureRangeStr_=NULL, 
	    char const *discFeatureRangeStr_=NULL, 
	    char const *preFrameRangeStr_=NULL, 
	    char const *segRangeStr_=NULL);

  ~PFileFile() {
    if (pfile)    delete pfile;
    if (dataFile) fclose(dataFile);
    if (buffer)   free(buffer);
    if (contBuf)  free(contBuf);
    if (discBuf)  free(discBuf);
    if (fileName) free(fileName);
  }
 
  // The number of available segments.
  unsigned numSegments() {
    assert(pfile);
    return pfile->num_segs();
  }

  // Begin sourcing data from the requested segment.
  // Must be called before any other operations are performed on a segment.
  bool openSegment(unsigned seg);

  // The number of frames in the currently open segment.
  unsigned numFrames()  {
    assert(pfile);
    return _numFrames;
  }

  // Load count frames of observed data, starting from first (physical),
  // in the current segment. 
  Data32 const *getFrames(unsigned first, unsigned count);

  // Number of continuous/discrete/total features in the file
  // after applying -frX and -irX
  unsigned numLogicalContinuous() { return _numLogicalContinuousFeatures; }
  unsigned numLogicalDiscrete()   { return _numLogicalDiscreteFeatures; }
  unsigned numLogicalFeatures()   { return _numLogicalFeatures; }
};

#endif
