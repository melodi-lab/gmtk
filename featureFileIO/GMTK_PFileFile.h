
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

  // for writable files
  FILE                  *out_fp;
  OutFtrLabStream_PFile *out_stream;
  unsigned               currentFeature;

 public:

  PFileFile(const char *name, unsigned nfloats, unsigned nints,
	    unsigned num, bool bswap, 
	    char const *contFeatureRangeStr_=NULL, 
	    char const *discFeatureRangeStr_=NULL, 
	    char const *preFrameRangeStr_=NULL, 
	    char const *segRangeStr_=NULL);

  PFileFile(const char *name, unsigned nfloats, unsigned nints, bool swap, int debug_level=0) 
    : fileName((char *)name)
  {
    buffer = NULL;
    contBuf = NULL;
    discBuf = NULL;
    fileName = NULL;
    dataFile = NULL;
    pfile = NULL;

    if (name == NULL)     
      error("PFileFile: output file name is NULL\n");    
    if ((out_fp = fopen(name,"wb")) == NULL)
      error("PFileFile: Can't open '%s' for output", name);
    out_stream = new OutFtrLabStream_PFile(debug_level,"",out_fp,nfloats,nints,1,swap ? 1 : 0);
    _numContinuousFeatures = nfloats;
    _numDiscreteFeatures   = nints;
    _numFeatures           = nfloats + nints;
    currentSegment = 0;
    currentFeature = 0;
  }

  ~PFileFile() {
    if (pfile)    delete pfile;
    if (dataFile) fclose(dataFile);
    if (buffer)   free(buffer);
    if (contBuf)  free(contBuf);
    if (discBuf)  free(discBuf);
    if (fileName) free(fileName);
    if (out_stream) {
      delete out_stream;
    }
    if (out_fp) {
      if (fclose(out_fp)) {
	error("Error closing output file '%s'\n", fileName);
      }
    }
  }
 

  // Write segment to the file (no need to call endOfSegment)
  void writeSegment(Data32 const *segment, unsigned nFrames);

  // returns true iff file supports random access writes via setFrame()
  bool seekable() { return true; }

  // Set frame # to write within current segemnt
  void setFrame(unsigned frame);

  // Write frame to the file (call endOfSegment after last frame of a segment)
  void writeFrame(Data32 const *frame);

  // Write frame to the file (call endOfSegment after last frame of a segment)
  void writeFeature(Data32 x);

  // Call after last writeFrame of a segment
  void endOfSegment();


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
