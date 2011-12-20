/*
 * GMTK_ObservationStream.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2011, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */

#ifndef GMTK_OBSERVATIONSTREAM_H
#define GMTK_OBSERVATIONSTREAM_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "machine-dependent.h"

// The ObservationStream provides a simple API to wrap around
// non-random access data sources like pipes and sockets.
// Just subclass ObservationStream and implement getNextFrame
// for the new type of stream (and add command line options 
// to instantiate it - aspect-oriented programming?) and the
// new stream type is supported by GMTK.
//
// Planned subtypes:
//   ASCIIStream  -   ASCII files
//   PFileStream  -   non-indexed PFiles
//   SocketStream -   Data streamed over IP (in what format?)
//   FilterStream -   ObservationStream wrapper for FIR, ARMA, etc transforms
//   FileStream   -   ObservationStream wrapper for seekable files (ObservationFile)

class ObservationStream {

 private:
  // -frX
  // -irX

 public:

  // Return the next frame from the stream.
  // May wait for data, may need a timeout.

  virtual Data32 *getNextFrame();
}

#endif
