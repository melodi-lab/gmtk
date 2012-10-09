
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


// This is an adaptor class that applies a linked list of 
// Filter objects (see GMTK_Filter.h) to an ObservationFile
// instance.
//
// FilterFile handles -postprX -frX -irX -transX for the
// individual ObservationFile instances. Another instance
// of FilterFile on top of MergeFile (see GMTK_MergeFile.h)
// handles -gpr and -posttrans. We could easily implement
// a -gsr option to select the global segment range by 
// passing in a segment range string. At present, the
// segment range selection is handled by the GMTK application
// programs themselves.

class FilterFile: public ObservationFile {

  Filter          *filter;
  ObservationFile *file;

  unsigned        _numFrames;

 public:
  
  FilterFile(Filter *filter, ObservationFile *file, 
	     char const *contFeatureRangeStr = NULL,
	     char const *discFeatureRangeStr = NULL,
	     char const *postpr = NULL)
    : ObservationFile(contFeatureRangeStr, discFeatureRangeStr, postpr),
      filter(filter), file(file)
  {
    subMatrixDescriptor wholeSegment(0U, 1U, 0U, 0U,
				     file->numLogicalContinuous(),
				     file->numLogicalDiscrete(), 1U, 0U, 0U);
    subMatrixDescriptor output = filter->describeOutput(wholeSegment);
    _numContinuousFeatures = output.numContinuous;
    _numDiscreteFeatures   = output.numDiscrete;
    _numFeatures           = _numContinuousFeatures + _numDiscreteFeatures;

    if (contFeatureRangeStr) {
      contFeatureRange = new Range(contFeatureRangeStr, 0, _numContinuousFeatures);
      assert(contFeatureRange);
      _numLogicalContinuousFeatures = contFeatureRange->length();
    } else
      _numLogicalContinuousFeatures = _numContinuousFeatures;
    if (discFeatureRangeStr) {
      discFeatureRange = new Range(discFeatureRangeStr, 0, _numDiscreteFeatures);
      assert(discFeatureRange);
      _numLogicalDiscreteFeatures = discFeatureRange->length();
    } else
      _numLogicalDiscreteFeatures = _numDiscreteFeatures;
    _numLogicalFeatures = _numLogicalContinuousFeatures + _numLogicalDiscreteFeatures;
  }


  ~FilterFile() {
    if (filter) delete filter;
    if (file) delete file;
  }


  // We have to use the ObservationFile's logical methods 
  // so that it can handle -srX, -frX, -irX, -preprX

  // filtering shouldn't change the number of segments...
  // -sdiffact might, but that is handled by MergeFile
  unsigned numSegments() {
    return file->numLogicalSegments();
  }


  // The number of ObservationFiles combined into the observation matrix.
  unsigned numFiles() { 
    assert(file);
    return file->numFiles(); 
  }

  bool openSegment(unsigned seg) {
    bool result = file->openLogicalSegment(seg);
    _numFrames = filter->getTotalOutputFrames(file->numLogicalFrames());
    return result;
  }


  // The number of frames in the filter's output - ask the
  // filter how many frames it will produce for the whole
  // input segment
  unsigned numFrames() { return _numFrames; }


  // Get the transformed data
  Data32 const *getFrames(unsigned first, unsigned count) {
    // find out which frames we need to feed the Filter to compute
    // the requested frames
    subMatrixDescriptor *inputDesc =
      filter->getRequiredInput(first, count, file->numLogicalContinuous(),
			       file->numLogicalDiscrete(),
			       file->numLogicalFrames());
    Data32 const *inputData =
      file->getLogicalFrames(inputDesc->firstFrame, inputDesc->numFrames);
    Data32 const *result = filter->transform(inputData, *inputDesc);
    subMatrixDescriptor::freeSMD(inputDesc);
    return result;
  }


  // Number of continuous, discrete, total features in filter's output
  // after applying -frX and -irX
  unsigned numLogicalContinuous() { return _numLogicalContinuousFeatures; }
  unsigned numLogicalDiscrete()   { return _numLogicalDiscreteFeatures; }
  unsigned numLogicalFeatures()   { return _numLogicalFeatures; }

};

#endif
