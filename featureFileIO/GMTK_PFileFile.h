
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
  unsigned   currentSegment;
  FILE      *dataFile;
  InFtrLabStream_PFile *pfile;

 public:

  PFileFile(const char *name, unsigned nfloats, unsigned nints,
	    unsigned num, bool bswap, 
	    char const *_contFeatureRangeStr=NULL, 
	    char const *_discFeatureRangeStr=NULL, 
	    char const *_preFrameRangeStr=NULL, 
	    char const *_segRangeStr=NULL);

  ~PFileFile() {
    if (pfile) delete pfile;
    if (dataFile) fclose(dataFile);
    if (buffer) free(buffer);
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
    return pfile->num_frames(currentSegment);
  }

  // Load count frames of observed data, starting from first,
  // in the current segment. count may be 0 to request loading
  // the entire data segment (frames [first, numFrames)).
  Data32 const *getFrames(unsigned first, unsigned count);

  unsigned numContinuous() {
    assert(pfile);
    return pfile->num_ftrs();
  }

  unsigned numDiscrete() {
    assert(pfile);
    return pfile->num_labs();
  }

  unsigned numFeatures() {return numContinuous() + numDiscrete();}

};

#endif
