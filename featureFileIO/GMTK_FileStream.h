
/*
 * GMTK_FileStream.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2012, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */

#ifndef GMTK_FILESTREAM_H
#define GMTK_FILESTREAM_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <assert.h>
#include <stdio.h>
#include <string.h>

using namespace std;

#include "machine-dependent.h"
#include "error.h"

#include "GMTK_ObservationFile.h"
#include "GMTK_ObservationStream.h"

// adapt an ObservationFile into an ObservationStream

class FileStream: public ObservationStream {

 protected:
  
  ObservationFile *obsFile;

  unsigned curSegment;
  unsigned numSegments;
  unsigned curFrame;  
  unsigned numFrames;

 public:

  FileStream(ObservationFile *file)
    :ObservationStream(), obsFile(file)
  {
    assert(file);
    nFloat      = file->numLogicalContinuous();
    nInt        = file->numLogicalDiscrete();
    numSegments = file->numLogicalSegments();
    frameData = new Data32[nFloat+nInt]; assert(frameData);
    logicalFrameData = new Data32[numLogicalFeatures()]; assert(logicalFrameData);

    curSegment = 0;
    curFrame = 0;
    if (! file->openLogicalSegment(0)) {
      error("ERROR: FileStream: failed to open segment 0\n");
    }
    numFrames = file->numLogicalFrames();
  }

  ~FileStream() { 
    if (obsFile) delete obsFile; 
    if (frameData) {
      delete [] frameData;
      frameData = NULL;
    }
    if (logicalFrameData) {
      delete [] logicalFrameData;
      logicalFrameData = NULL;
    }
  }

  Data32 const *getNextFrame() {
    if (curFrame < numFrames) {
      Data32 const *frame = obsFile->getLogicalFrames(curFrame++, 1);
#if 0
printf("FileStream %3u %3u >", curSegment, curFrame);
float *fp = (float *)frame;
for (unsigned i=0; i < nFloat; i+=1) {
  printf(" %f", *(fp++));
}
unsigned *ip = (unsigned *) fp;
for (unsigned i=0; i < nInt; i+=1) {
  printf(" %u", *(ip++));
}
printf("\n");
#endif
      memcpy(frameData, frame, (nFloat + nInt) * sizeof(Data32));
      return frameData;
    }
    if (curSegment < numSegments) curSegment += 1;
    if (curSegment < numSegments) {
      if (! obsFile->openLogicalSegment(curSegment)) {
	error("ERROR: FileStream::getNextFrame failed to open segment %u\n", curSegment);
      }
      numFrames = obsFile->numLogicalFrames();
      curFrame  = 0;
    }
    return NULL;
  }

  bool EOS() {return curSegment >= numSegments;}
};

#endif
