/*
 * GMTK_FileSourceNoCache.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2012, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */

#ifndef GMTK_FILESOURCENOCACHE_H
#define GMTK_FILESOURCENOCACHE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "GMTK_FileSource.h"

class FileSourceNoCache: public FileSource {

 public:

#define DEFAULT_BUFFER_SIZE       (16 * 1024 * 1024)
#ifndef DEFAULT_FILE_WINDOW_SIZE
#define DEFAULT_FILE_WINDOW_SIZE  (4)
#define DEFAULT_FILE_WINDOW_DELTA (100)
#endif
#define DEFAULT_FILE_WINDOW_BYTES (DEFAULT_FILE_WINDOW_SIZE * 1024 * 1024)

  FileSourceNoCache(ObservationFile *file,
		    unsigned windowBytes = DEFAULT_FILE_WINDOW_BYTES, 
		    unsigned deltaFrames = DEFAULT_FILE_WINDOW_DELTA,
		    unsigned bufferSize  = DEFAULT_BUFFER_SIZE, 
		    unsigned startSkip=0, unsigned endSkip=0,
		    int justificationMode=0); 

  // Construct an invalid empty FileSource. Turn it into a valid one
  // with the initialize() method.
  FileSourceNoCache() : FileSource() { }

  // Returns a pointer to the requested frames in the cookedBuffer.
  // Any required frames before/after (see _minPastFrames and
  // _minFutureFrames above) are guaranteed to be in the cookedBuffer
  // as well. 
  Data32 const *loadFrames(unsigned first, unsigned count);

};

#endif
