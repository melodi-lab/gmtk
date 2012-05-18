
/*
 * GMTK_ASCIIFile.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2011, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */

#ifndef GMTK_ASCIIFILE_H
#define GMTK_ASCIIFILE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <string>
using namespace std;

#include "machine-dependent.h"
#include "error.h"
#include "general.h"

#include "GMTK_ObservationFile.h"


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

class ASCIIFile: public ObservationFile {

  bool        cppIfAscii;
  char const *cppCommandOptions;
  unsigned    numFileNames;             // number of filenames

  char      *fofName;
  FILE      *fofFile;                   // this file (list of file names)
  unsigned   nFloats;
  unsigned   nInts;
  unsigned   nFrames;                   // in current segment

  char     **dataNames;  // pointers to individual filenames (into fofBuf)

  Data32    *buffer;     // data for current segment

 public:

  ASCIIFile(const char *name, unsigned nfloats, unsigned nints, 
	    unsigned num, bool cppIfAscii, char const* cppCommandOptions=NULL,
	    char const *contFeatureRangeStr_=NULL, 
	    char const *discFeatureRangeStr_=NULL, 
	    char const *preFrameRangeStr_=NULL, 
	    char const *segRangeStr_=NULL);


  ~ASCIIFile() {
    if (fofName) delete fofName;
    if (dataNames) {
      for (unsigned i = 0; i < numFileNames; i++)
	if (dataNames[i])
	  delete [] dataNames[i];
      delete [] dataNames;
    }
    if (buffer) delete [] buffer;
  }
 
  // The number of available segments.
  unsigned numSegments() {return numFileNames;}

  // Begin sourcing data from the requested segment.
  // Must be called before any other operations are performed on a segment.
  bool openSegment(unsigned seg);

  // The number of frames in the currently open segment.
  unsigned numFrames()  {
    assert(buffer);
    return nFrames;
  }

  // Load count frames of observed data, starting from first,
  // in the current segment. count may be 0 to request loading
  // the entire data segment (frames [first, numFrames)).
  Data32 const *getFrames(unsigned first, unsigned count) {
    assert(buffer);
    assert(first < nFrames);
    assert(first + count <= nFrames);
    return buffer + first * (nFloats + nInts);
  }

  unsigned numContinuous() {return nFloats;}
  unsigned numDiscrete() {return nInts;}
  unsigned numFeatures() {return nFloats + nInts;}

};

#endif
