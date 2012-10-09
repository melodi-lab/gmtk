
/*
 * GMTK_FileSrcStream.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2012, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */

#ifndef GMTK_FILESRCSTREAM_H
#define GMTK_FILESRCSTREAM_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <assert.h>
#include <stdio.h>
using namespace std;

#include "machine-dependent.h"
#include "error.h"

#include "GMTK_FileSource.h"
#include "GMTK_ObservationStream.h"

// adapt a FileSource into an ObservationStream

class FileSrcStream: public ObservationStream {

 protected:
  
  FileSource *obsFile;

  unsigned curSegment;
  unsigned numSegments;
  unsigned curFrame;  
  unsigned numFrames;

 public:

  FileSrcStream(FileSource *file)
    :ObservationStream(), obsFile(file)
  {
    if (!file) {
      error("ERROR: no observation input source specified\n");
    }
    nFloat      = file->numContinuous();
    nInt        = file->numDiscrete();
    numSegments = file->numSegments();
    frameData = new Data32[nFloat+nInt]; assert(frameData);
    logicalFrameData = new Data32[numLogicalFeatures()]; assert(logicalFrameData);

    curSegment = 0;
    curFrame = 0;
    if (! file->openSegment(0)) {
      error("ERROR: FileStream: failed to open segment 0\n");
    }
    numFrames = file->numFrames();
  }

  ~FileSrcStream() {
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
      Data32 const *frame = obsFile->loadFrames(curFrame++, 1);
      memcpy(frameData, frame, (nFloat + nInt) * sizeof(Data32));
      return frameData;
    }
    if (curSegment < numSegments) curSegment += 1;
    if (curSegment < numSegments) {
      if (! obsFile->openSegment(curSegment)) {
	error("ERROR: FileStream::getNextFrame failed to open segment %u\n", curSegment);
      }
      numFrames = obsFile->numFrames();
      curFrame  = 0;
    }
    return NULL;
  }

  bool EOS() {return curSegment >= numSegments;}
};

#endif
