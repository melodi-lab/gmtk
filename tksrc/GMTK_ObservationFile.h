
/*
 * GMTK_ObservationFile.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2011, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */

#ifndef GMTK_OBSERVATIONFILE_H
#define GMTK_OBSERVATIONFILE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "machine-dependent.h"

// The ObservationFile provides a simple API to wrap around
// random access data file formats.
// Just subclass ObservationFile and implement getFrames
// for the new type of file (and add command line options 
// to instantiate it - aspect-oriented programming?) and the
// new file type is supported by GMTK.
//
// Planned subtypes:
//   ASCIIFile    -   ASCII files (read entirely into memory)
//   PFileFile    -   indexed PFiles (non-indexed read entirely into memory)
//   HDF5File     
//   HTKFile      
//   BinaryFile   
//   FilterFile   -   ObservationFile wrapper for IIR, ARMA, etc transforms

class ObservationFile {

 private:
  // -frX
  // -irX

 public:

  // The number of available segments.
  virtual unsigned numSegments();

  // Begin sourcing data from the requested segment.
  // Must be called before any other operations are performed on a segment.
  virtual void openSegment(unsigned seg);

  // The number of frames in the currently open segment.
  virtual unsigned numFrames();

  // Load count frames of observed data, starting from first,
  // in the current segment. count may be 0 to request loading
  // the entire data segment (frames [first, numFrames)).
  virtual Data32 *getFrames(unsigned first, unsigned count);
}

#endif
