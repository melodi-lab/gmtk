
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
  
  // Apply filter stack to the ObservationFile and present the
  // result as an ObservationFile

  // FilterFile handles -postprX
  FilterFile(Filter *filter, ObservationFile *file, 
	     char const *contFeatureRangeStr = NULL,
	     char const *discFeatureRangeStr = NULL,
	     char const *postpr = NULL)
    : ObservationFile(contFeatureRangeStr, discFeatureRangeStr, postpr),
      filter(filter), file(file)
  { }


  ~FilterFile() {
    if (filter) delete filter;
    if (file) delete file;
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
    subMatrixDescriptor wholeSegment(0U, file->numLogicalFrames(), 0U, 0U,
				     file->numLogicalContinuous(),
				     file->numLogicalDiscrete(),
				     file->numLogicalFrames());
    subMatrixDescriptor output = filter->describeOutput(wholeSegment);
    return output.numFrames;
  }


#if 0
  // needed for -fdiffactX -- it applies to the number of frames in a
  // segment before -transX is performed 
  // FIXME - this might be somewhat silly; it looks like the fdiffact
  // actions are actually applied to the transform output although
  // the checking & setup are done using the # of input frames...
  unsigned preFilterFrameCount() { return file->numLogicalFrames(); }
#endif

  // Get the transformed data
  Data32 const *getFrames(unsigned first, unsigned count) {
//printf(" requesting [%u,%u) : ", first, first+count);
    subMatrixDescriptor *inputDesc =
      filter->getRequiredInput(first, count, file->numLogicalContinuous(),
			       file->numLogicalDiscrete(),
			       file->numLogicalFrames());
#if 0
    printf("input [%u + %u, %u - %u) \n", inputDesc->firstFrame, inputDesc->historyFrames,
                                          inputDesc->firstFrame + inputDesc->numFrames,
                                          inputDesc->futureFrames);
#endif
#if 1
    Data32 const *inputData =
      file->getLogicalFrames(inputDesc->firstFrame, inputDesc->numFrames);
    return filter->transform(inputData, *inputDesc);
#else
    Data32 const *inputData =
      file->getLogicalFrames(inputDesc->firstFrame, inputDesc->numFrames);
    subMatrixDescriptor outputDesc;
    Data32 const *outputData = filter->transform(inputData, *inputDesc, &outputDesc);
    printf("output [%u + %u, %u - %u) : ", outputDesc.firstFrame, outputDesc.historyFrames,
                                          outputDesc.firstFrame + outputDesc.numFrames,
                                          outputDesc.futureFrames);
    return outputData;
#endif
  }


  // Number of continuous, discrete, total features in filter's output

  // Ask the filter how many continuous features it will produce 
  unsigned numContinuous() {
    subMatrixDescriptor wholeSegment(0U, 1U, 0U, 0U,
				     file->numLogicalContinuous(),
				     file->numLogicalDiscrete(), 1U);
    subMatrixDescriptor output = filter->describeOutput(wholeSegment);
    return output.numContinuous;
  }


  unsigned numDiscrete() {
    subMatrixDescriptor wholeSegment(0U, 1U, 0U, 0U,
				     file->numLogicalContinuous(),
				     file->numLogicalDiscrete(), 1U);
    subMatrixDescriptor output = filter->describeOutput(wholeSegment);
    return output.numDiscrete;
  }


  unsigned numFeatures() {
    return numContinuous() + numDiscrete();
  }

};

#endif
