
/*
 * GMTK_FilterFile.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2012, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */

#ifndef GMTK_FILTERFILE_H
#define GMTK_FILTERFILE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "machine-dependent.h"
#include "range.h"
#include "GMTK_ObservationFile.h"
#include "GMTK_Filter.h"

class FilterFile: public ObservationFile {

  Filter          *filter;
  ObservationFile *file;

 public:
  
  // Apply filter to the ObservationFile and present the
  // result as an ObservationFile

  // FilterFile handles -postprX
  Filterfile(Filter *filter, ObservationFile *file, Range *postpr)
    :filter(filter), file(file), preFrameRange(postpr) {
  }

  // We have to use the ObservationFile's logical methods 
  // so that it can handle -srX, -frX, -irX, -preprX

  // filtering shouldn't change the number of segments...
  // -sdiffact might, but that should be handled by the
  // FileSource
  unsigned numSegments() {
    return file->numLogicalSegments();
  }

  bool openSegment(unsigned seg) {
    return file->openLogicalSegment(seg);
  }

  // The number of frames in the filter's output - ask the
  // filter how many frames it will produce for the whole
  // input file
  unsigned numFrames() {
    subMatrixDescriptor wholeSegment(0, file->numLogicalFrames(),
				     file->numLogicalContinuous(),
				     file->numLogicalDiscrete());
    subMatrixDescriptor output = filter->describeOutput(wholeSegment);
    return output.numFrames;
  }

  // Get the transformed data
  Data32 const *getFrames(unsigned first, unsigned count) {
    return filter.transform(file, first, count);
  }

  // Number of continuous, discrete, total features in filter's output

  // Ask the filter how many continuous features it will produce 
  unsigned numContinuous() {
    // ask for the full segment output, since the filter might 
    // have a minimum number of input frames?
    subMatrixDescriptor wholeSegment(0, file->numLogicalFrames(),
				     file->numLogicalContinuous(),
				     file->numLogicalDiscrete());
    subMatrixDescriptor output = filter->describeOutput(wholeSegment);
    return output.numContinous;
  }

  unsigned numDiscrete() {
    subMatrixDescriptor wholeSegment(0, file->numLogicalFrames(),
				     file->numLogicalContinuous(),
				     file->numLogicalDiscrete());
    subMatrixDescriptor output = filter->describeOutput(wholeSegment);
    return output.numDiscrete;
  }

  unsigned numFeatures() {
    return numContinuous() + numDiscrete();
  }

}

#endif
