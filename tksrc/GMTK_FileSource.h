/*
 * GMTK_FileSource.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2011, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */

#ifndef GMTK_FILESOURCE_H
#define GMTK_FILESOURCE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "machine-dependent.h"
#include "GMTK_ObservationSource.h"
#include "GMTK_ObservationFile.h"

// FileSource handles
// random access data sources like the various binary file
// formats.

class FileSource : ObservationSource {

 private:

  // buffer to hold transformed cahced/prefected observations
  Data32 cookedBuffer[];

  // the files assembled to form the observations
  ObservationFile file[];

 public:
  
  // each file provides its data for the requested submatrix,
  // transformed and "read to eat." So the FileSource
  // just manages assembling them to satisfy the loadFrames()
  // calls, prefetching/caching. For archipelagos, each thread
  // gets its own FileSource (all aimed at the same files, of course).
  FileSource(ObservationFile file[]); 


  // The number of available segments.
  unsigned numSegments();

  // Begin sourcing data from the requested segment.
  // Must be called before any other operations are performed on a segment.
  void openSegment(unsigned seg);

  // The number of frames in the currently open segment.
  unsigned numFrames();

  Data32 *loadFrames(unsigned first, unsigned count) {
    // if requested frames are already in cookedBuffer
    //   return &cookedBuffer[index of first requested frame]
    // adjust (first,count) for prefetching -> (first',count')
    // adjust cookedBuffer destination
    // getFrames(first', count') from each file into cookedBuffer
    // return &cookedBuffer[destination of first]
    return NULL;
  }
}

#endif
