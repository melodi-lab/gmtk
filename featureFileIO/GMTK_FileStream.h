
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
using namespace std;

#include "machine-dependent.h"
#include "error.h"

#include "GMTK_FileSource.h"
#include "GMTK_ObservationStream.h"

// adapt a FileSource into an ObservationStream

class FileStream: public ObservationStream {

 protected:
  
  FileSource *obsFile;

  unsigned curSegment;
  unsigned numSegments;
  unsigned curFrame;  
  unsigned numFrames;

 public:

  FileStream(FileSource *file)
    :ObservationStream(), obsFile(file)
  {
    assert(file);
    nFloat      = file->numContinuous();
    nInt        = file->numDiscrete();
    numSegments = file->numSegments();

    curSegment = 0;
    curFrame = 0;
    if (! file->openSegment(0)) {
      error("ERROR: FileStream: failed to open segment 0\n");
    }
    numFrames = file->numFrames();
  }

  ~FileStream() { if (obsFile) delete obsFile; }

  Data32 const *getNextFrame() {
    if (curFrame < numFrames) {
      return obsFile->loadFrames(curFrame++, 1);
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
