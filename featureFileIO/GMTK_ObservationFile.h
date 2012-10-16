
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

#include <stdlib.h>

#include "error.h"
#include "machine-dependent.h"
#include "range.h"

// The ObservationFile provides a simple API to wrap around
// random access data file formats. To support a new file
// format, just subclass ObservationFile and implement:
//
//   numSegments()
//   openSegment()
//   numFrames()
//   getFrames()
//   numContinuous()
//   numDiscrete()
//
// for the new type of file (and add command line options 
// to instantiate it) and the new file type is supported by 
// GMTK. You should add support for the new file type to the
// instantiateFile() and formatStrToNumber functions to 
// facilitate command line argument parsing.
//
// The methods above all work with the file's true physical
// layout. Client code should access the file through the
// *Logical* methods, which apply the feature, frame, and segment
// range selections.
//
// Subclasses should initialize the *RangeStr members to support
// the -frX -irX -preprX -srX arguments. The 
// _num{Continuous,Discrete,}Features should be set in
// the constructors. Initializing them is sufficient to get them 
// working; the *Logical*() methods can be over-ridden if the file
// format supports a better implementation than the simple-minded
// one ObservationFile provides.
// 
// Subclasses:
//   ASCIIFile     -   ASCII files (reads segment entirely into memory)
//   FlatAsciiFile -   seg #   frame #   f_0 ... f_n i_0 ... i_m
//   PFileFile     -   indexed PFiles (indexing specified at file open)
//   HDF5File     
//   HTKFile      
//   BinaryFile   
//   FilterFile    -   ObservationFile wrapper for IIR, FIR, etc transforms
//   MergeFile     -   Combines multiple ObservationFiles into a single 
//                     logical file

class ObservationFile {

 protected:

  char const *contFeatureRangeStr;  // -frX
  Range      *contFeatureRange;
  char const *discFeatureRangeStr;  // -irX
  Range      *discFeatureRange;
  char const *preFrameRangeStr;     // -preprX
  Range      *preFrameRange;
  char const *segRangeStr;          // -srX
  Range      *segRange;

  Data32 *logicalObservationBuffer;
  unsigned logicalObsBufSize;       // (in Data32's)

  // # of physical features
  unsigned _numContinuousFeatures;
  unsigned _numDiscreteFeatures;
  unsigned _numFeatures;

  // # of logical features (after -frX and -irX)
  unsigned _numLogicalContinuousFeatures;
  unsigned _numLogicalDiscreteFeatures;
  unsigned _numLogicalFeatures;

 public:

  ObservationFile(char const *contFeatureRangeStr_=NULL, 
		  char const *discFeatureRangeStr_=NULL,
		  char const *preFrameRangeStr_=NULL, 
		  char const *segRangeStr_=NULL)
    : contFeatureRangeStr(contFeatureRangeStr_), contFeatureRange(NULL),
      discFeatureRangeStr(discFeatureRangeStr_), discFeatureRange(NULL),
      preFrameRangeStr(preFrameRangeStr_), preFrameRange(NULL),
      segRangeStr(segRangeStr_), segRange(NULL), 
      logicalObservationBuffer(NULL), logicalObsBufSize(0)
    {
    }

  virtual ~ObservationFile() {
    if (contFeatureRange) delete contFeatureRange;
    if (discFeatureRange) delete discFeatureRange;
    if (preFrameRange) delete preFrameRange;
    if (segRange) delete segRange;
    if (logicalObservationBuffer) {
      free(logicalObservationBuffer);
      logicalObservationBuffer = NULL;
    }
  }

  // The number of physical (before -srX) segments in the file.
  virtual unsigned numSegments() = 0;

  // The number of ObservationFiles combined into the observation matrix.
  virtual unsigned numFiles() { return 1; };

  // Begin sourcing data from the requested physical (before -srX) segment.
  // Must be called before any other operations are performed on a segment.
  virtual bool openSegment(unsigned seg) = 0;

  // The number of physical frames (before -preprX) in the currently open segment.
  virtual unsigned numFrames() = 0;

  // Load count frames of observed data, starting from first (physical),
  // in the current segment. 
  virtual Data32 const *getFrames(unsigned first, unsigned count) = 0;

  // The number of continuous, discrete, total features (all physical).
  // Note that they may be called before openSegment()

  virtual unsigned numContinuous() { return _numContinuousFeatures; }
  virtual unsigned numDiscrete()   { return _numDiscreteFeatures; }
  virtual unsigned numFeatures()   { return _numFeatures; }

  // logical -> physical translation
  
  // simple-minded default implementations provided, but more
  // sophisticated file formats (e.g., HDF5) might benefit from
  // implementing them directly with their native API

  virtual unsigned numLogicalSegments();  // after -srX
 
  virtual bool openLogicalSegment(unsigned seg); // after -srX

  virtual unsigned numLogicalFrames();   // after -preprX

  virtual Data32 const *getLogicalFrames(unsigned first, unsigned count); // after -preprX

  virtual unsigned numLogicalContinuous();   // after -frX

  virtual unsigned numLogicalDiscrete();     // after -irX

  virtual unsigned numLogicalFeatures() {
    // be aware of -frX and -irX if you over-ride this
    return _numLogicalFeatures;
  }

};


// Instantiate the appropriate ObservationFile subclass according 
// to the command line arguments. number is the file number (the
// X in the -argX arguments).
ObservationFile *
instantiateFile(unsigned ifmt, char *ofs, unsigned nfs, unsigned nis,
		unsigned number, bool iswp, bool Cpp_If_Ascii, 
		char *cppCommandOptions, char const *frs, char const *irs, 
		char const *prepr, char const *sr);

// Converts the command line -fmtX string to an integer (enum)
int
formatStrToNumber(char const *fmt);

#endif
