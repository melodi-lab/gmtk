
/*
 * GMTK_MergeFile.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2012, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */

#ifndef GMTK_MERGEFILE_H
#define GMTK_MERGEFILE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "machine-dependent.h"
#include "range.h"

#include "GMTK_ObservationFile.h"
#include "GMTK_Filter.h"

// Combine multiple ObservationFile instances into a single virtual file.

class MergeFile: public ObservationFile {

  // the files assembled to form the observations
  unsigned nFiles;
  ObservationFile **file;

  unsigned *floatStart; // The ith file's continuous features start here.
                        //    This is the offset from the start of a complete 
                        //    frame of the ith file's first continuous feature.
  unsigned *intStart;   // Likewise for discrete features
  
  int       ftrcombo;   // How should the continuous features from the multiple files be combined?
                        //   FTROP_NONE   files are concatenated (floats first, then ints)
                        //   FTROP_{ADD,SUB,MUL,DIV} continuous features are combined with the
                        //     the indicated operator (ints still concatenated)

  unsigned const *sdiffact;   // how to adjust for files with different # of segmentss
  unsigned const *fdiffact;   // how to adjust for segments with different # of frames

  unsigned _numSegments;      // after considering -sdiffact
  unsigned _numFrames;        // after considering -fdiffact

  int      segment;           // currently open segment, -1 if none yet

  Data32  *buffer;            // data for current segment
  unsigned buffSize;          // in Data32's
  unsigned bufStride;         // increment between frames in buffer


  // Map requested global segment # to logical segment # in specified file
  unsigned adjustForSdiffact(unsigned fileNum, unsigned seg);


  // Map requested merged frame range [first,first+count) to individual file 
  // fileNum's pre-fdiffactX frame range [adjFirst,adjFirst+adjCount). Also 
  // returns deltaT, the difference in frames between the merged segment 
  // length and the individual file's segment length.
  void adjustForFdiffact(unsigned first, unsigned count, unsigned fileNum,
			 unsigned &adjFirst, unsigned &adjCount, unsigned &deltaT);
 public:

  MergeFile(unsigned nFiles, ObservationFile *file[], 
	    unsigned const *sdiffact = NULL, 
	    unsigned const *fdiffact = NULL,
	    int ftrcombo=FTROP_NONE);

  virtual ~MergeFile() {
    if (file) {
      for (unsigned i=0; i < nFiles; i+=1) {
	if (file[i]) delete file[i];
      }
      delete [] file;
    }
    if (floatStart) delete [] floatStart;
    if (intStart) delete [] intStart;
    if (buffer) free(buffer);
  }

  // The number of segments after -sdiffact
  unsigned numSegments() { return _numSegments; }

  // The number of ObservationFiles combined into the observation matrix.
  unsigned numFiles() { 
    unsigned sum = 0;
    for (unsigned i=0; i < nFiles; i+=1) {
      assert(file);
      sum += file[i]->numFiles();
    }
    return sum;
  };

  bool openSegment(unsigned seg);

  // The number of frames in the currently open segment after -fdiffact
  unsigned numFrames() {
    assert(segment >= 0);
    return _numFrames;
  }

  Data32 const *getFrames(unsigned first, unsigned count);

  // Number of continuous/discrete/total features in the file
  // after applying -frX and -irX
  unsigned numLogicalContinuous() { return _numLogicalContinuousFeatures; }
  unsigned numLogicalDiscrete()   { return _numLogicalDiscreteFeatures; }
  unsigned numLogicalFeatures()   { return _numLogicalFeatures; }

};

#endif
