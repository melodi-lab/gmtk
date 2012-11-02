
/*
 * GMTK_FileSourceNoCache.cc
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2012, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <string.h>

#include "error.h"
#include "general.h"
#include "machine-dependent.h"
#include "debug.h"

#include "GMTK_FileSourceNoCache.h"

FileSourceNoCache::
FileSourceNoCache(ObservationFile *file, 
		  unsigned windowBytes, unsigned deltaFrames, unsigned bufferSize, 
		  unsigned startSkip, unsigned endSkip, int justificationMode)
{
  initialize(file, windowBytes, deltaFrames, bufferSize, startSkip, endSkip, justificationMode, false);
}


// This method loads the requested frames (along with any necessary
// preceding or following frames) into the cookedBuffer if they are
// not already present and returns a pointer to the first requested
// frame. It tries to prefetch a full window of frames at once when
// the inference code gets close to the end of the buffered frames.
Data32 const *
FileSourceNoCache::loadFrames(unsigned first, unsigned count) {

  first += _startSkip + justificationOffset; // adjust for -startSkip and -justification

  if (_minPastFrames <= first && first + count + _minFutureFrames <= numBufferedFrames)
  { 
    return cookedBuffer + first * bufStride;
  } else {
    error("ERROR: FileSourceNoCache::loadFrames: requested frames [%u,%u), but only "
	  "frames [0,%u) are available", first, first+count+_minFutureFrames, numBufferedFrames);
    return NULL; // this statement isn't reached, just here to silence warning
  }
}

