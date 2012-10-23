
/*
 * GMTK_Filter.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2012, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */

#ifndef GMTK_FILTER_H
#define GMTK_FILTER_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <assert.h>

#include "machine-dependent.h"

#include "GMTK_SubmatrixDescriptor.h"


// This is the base class for observation filtering/transformation.
// It implements the identity filter as an example. Subclasses
// implement FIR, IIR, ARMA filters, etc. The Filter instances
// can be stacked in a linked list to apply multiple transformations
// to the observations. The FilterFile (see GMTK_FilterFile.h) class
// is responsible for applying the Filters to ObservationFile instances.
// FilterStream applies a single Filter to an ObservationStream.
//
// Subclasses should implement 
//   getReqiredInput()
//   describeLocalOutput()
//   localTransform()
//   getNextFrameInfo()
//   getEOSFrameInfo()
//   getTotalOutputFrames()
//   numOutputFeatures()

// You should also update parseTransform() and instantiateFilters()
// to facilitate command line argument parsing.

class StreamSource;

class Filter {
  friend class StreamSource;

// Transformed frame data flows this way (note that the arrows run
// in the opposite direction of the nextFilter links)
//
//   input submatrix -->  Filter ( --> nextFilter)*  --> output submatrix
// (ObservationSource)                                     (Inference)
  
 protected:
  Filter *nextFilter;

  unsigned frameNum; // state for streaming mode

 public:

  Filter(Filter *nextFilter = NULL) : nextFilter(nextFilter), frameNum(0) {}

  
  virtual ~Filter() { 
    if (nextFilter) delete nextFilter; 
  }

  void appendFilter(Filter *filter) {
    if (nextFilter == NULL) {
      nextFilter = filter;
      return;
    }
    nextFilter->appendFilter(filter);
  }


  // This is called by FilterStream after the last frame in a segment
  // is read. Any state can be reset in preparation for the next segment.
  virtual void nextStreamSegment() { frameNum = 0; }

  // FileStream has a random access backing store so that the frames 
  // needed to produce the requested filter output can always be
  // accessed. In contrast, FilterStream can only maintain a small
  // queue of the most recent input frames to which the filter is applied
  // to produce the next frame(s) of filter output. getNextFrameInfo()
  // tells FilterStream how to manage its input frame queue - drop dropOldIn
  // frames for the queue, enqueue numNewIn new frames; then applying the
  // filter will produce numNewOut output frames. This method requires state
  // since the first and last frames may need to be handled differently from
  // the middle frames within a segment. inputContinuous and inputDiscrete
  // are the number of real and discrete features in the input frames. The
  // input subMatrixDescriptor is setup to tell the filter how to interpret 
  // the input frames.

  virtual void getNextFrameInfo(unsigned &numNewIn, unsigned &dropOldIn, unsigned &numNewOut,
				unsigned inputContinuous, unsigned inputDiscrete,
				subMatrixDescriptor &input)
  {
    numNewIn = 1;
    dropOldIn = frameNum == 0 ? 0 : 1; frameNum = 1; // only need to know first frame or not
    numNewOut= 1;
    input.firstFrame = 0; 
    input.numFrames = 1;
    input.historyFrames = 0; 
    input.futureFrames = 0;
    input.numContinuous = inputContinuous; 
    input.numDiscrete = inputDiscrete;
    input.fullMatrixFrameCount = 1;
    input.requestedFirst = 0; 
    input.requestedCount = 1;
  }

  // A stream segment may end before the FilterStream can read the numNewIn 
  // frames requested by getNExtFrameInfo(). In that case, the FitlerStream
  // will call getEOSFrameInfo, passing in the number of requested frames that
  // were unavailable in numFramesShort. It returns numNewOut, the number of 
  // output frames that will be produced from these last frames in the segment.
  // The input subMatrixDescriptor again describes how the input frames should
  // be interpreted by the Fitler.

  virtual void getEOSFrameInfo(int numFramesShort, unsigned &numNewOut, subMatrixDescriptor &input) {
    numNewOut = 0;
  }

  
  // The filter's client (e.g. inference) needs the 
  // [first,first+count)frames of the filter's output.
  // getRequiredInput() returns the portion of the fitler's 
  // input (e.g. ObservationFile) necessary to produce the 
  // requested output frames.
  //
  // Note that first and count describe the desired OUTPUT frames.
  // inputContinuous, inputDiscrete and inputTotalFrames are the number
  // of continuous & discrete features and frames in the Filter's INPUT.

  virtual subMatrixDescriptor *
    getRequiredInput(unsigned first, unsigned count, 
		     unsigned inputContinuous, unsigned inputDiscrete,
		     unsigned inputTotalFrames)
  {
    subMatrixDescriptor *nextFilterInput = NULL;
    if (nextFilter) {
      // the client's getting the data from the Filter(s) AFTER me, so
      // ask it what frames I need to provide it as input (it gets as
      // input the # of continuous and discrete features and frames
      // I produce as output).
      unsigned outputContinuous = inputContinuous;  // # floats I produce
      unsigned outputDiscrete   = inputDiscrete;    // # ints I produce
      unsigned outputFrames     = inputTotalFrames; // total # frames I produce
      nextFilterInput = nextFilter->getRequiredInput(first, count, 
						     outputContinuous, outputDiscrete,
						     outputFrames);
      assert(nextFilterInput);
      first = nextFilterInput->requestedFirst;
      count = nextFilterInput->requestedCount;
    }
    return subMatrixDescriptor::getSMD(first, count, 0, 0, inputContinuous, 
	     inputDiscrete, inputTotalFrames, first, count, nextFilterInput);
  }


  // How many total frames will the Filter produce for a segment of
  // length totalInputFrames?
  virtual unsigned getTotalOutputFrames(unsigned totalInputFrames) {
    if (nextFilter) {
      return nextFilter->getTotalOutputFrames(totalInputFrames);
    } else {
      return totalInputFrames;
    }
  }


  // How many continuous/discrete features in the output given
  // the number of input continous/discrete features.

  virtual void numOutputFeatures(unsigned inputContinuous, unsigned inputDiscrete,
				 unsigned &outputContinuous, unsigned &outputDiscrete) 
  {
    outputContinuous = inputContinuous;
    outputDiscrete = inputDiscrete;
  }

  // What will the output look like if the input described by
  // inputDescription is feed into this Filter?
  virtual subMatrixDescriptor
    describeLocalOutput(subMatrixDescriptor const &inputDescription)
  {
    // describe my output given the specified input
    subMatrixDescriptor myOutput = inputDescription;
    return myOutput;
  }


  // What will the output look like if the input described by
  // inputDescription is feed through the entire filter stack?
  subMatrixDescriptor
    describeOutput(subMatrixDescriptor const &inputDescription)
  {
    subMatrixDescriptor myOutput = describeLocalOutput(inputDescription);
    if (nextFilter) {
      // the next filter gets my output as its input
      myOutput = nextFilter->describeOutput(myOutput);
    } 
    return myOutput;
  }
  

  // Returns the filter's output given the inputSubMatrix
  // described by inputDescription (presumably created by 
  // getRequiredInput()).  
  //
  // If outputDescription is non-NULL, 
  // it will describe the portion of filter's output returned 

  virtual Data32 const *
    localTransform(Data32 const *inputSubMatrix,
		   subMatrixDescriptor const &inputDescription,
		   subMatrixDescriptor *outputDescription=NULL) 
  {
    Data32 const *outputData = inputSubMatrix; // perform my transformation
    subMatrixDescriptor myOutput = describeLocalOutput(inputDescription);
    if (outputDescription) *outputDescription = myOutput;
    return outputData;
  }


  // apply the local transform, then pass its output down the Filter stack
#if 0
  Data32 const *
    transform(Data32 const *inputSubMatrix,
	      subMatrixDescriptor const &inputDescription,
	      subMatrixDescriptor *outputDescription=NULL) 
  {
    error("Deprecated transform");
    return transform(inputSubMatrix, inputDescription, 
		     inputDescription.firstFrame, 
		     inputDescription.numFrames,
		     outputDescription);
  }
#endif
  Data32 const *
    transform(Data32 const *inputSubMatrix,
	      subMatrixDescriptor const &inputDescription,
	      subMatrixDescriptor *outputDescription=NULL) 
  {
    subMatrixDescriptor myOutput;
    Data32 const *outputData = localTransform(inputSubMatrix,
					      inputDescription, 
					      &myOutput);
    if (nextFilter) {
      outputData = 
	nextFilter->transform(outputData, *(inputDescription.next), outputDescription);
    } else {
      if (outputDescription) *outputDescription = myOutput;
    }
    return outputData;
  }

};



/**
 * parses a transformation string and returns the next transformation
 * to perform.  Advances the string pointer to the next
 * transformation.
 *
 * side effects: string is modified to point to the next transformation substring.
 *
 * returns END_STR (-1) if the end of the string is reached,
 * UNRECOGNIZED_TRANSFORM (-2) if the transformation substring is not
 * recognized.  Otherwise returns the unsigned that corresponds to the
 * transformation substring.
 * */
int 
parseTransform(char*& trans_str, int& magic_int, double& magic_double, char *&filterFileName);


Filter *
instantiateFilters(char *filterStr, unsigned numContinuous, unsigned numDiscrete);

#define NONE_LETTER 'X'
#define TRANS_AFFINE_LETTER 'A'
#define TRANS_NORMALIZATION_LETTER 'N'
#define TRANS_MEAN_SUB_LETTER 'E'
#define TRANS_UPSAMPLING_LETTER 'U'
#define TRANS_HOLD_LETTER 'H'
#define TRANS_SMOOTH_LETTER 'S'
#define TRANS_MULTIPLICATION_LETTER 'M'
#define TRANS_ARMA_LETTER 'R'
#define TRANS_OFFSET_LETTER 'O'
#define FILTER_LETTER 'F'
#define END_STR (-1)
#define UNRECOGNIZED_TRANSFORM (-2)
#define SEPARATOR "_"

#define MAX_NUM_STORED_FILTERS 10
#define MAX_FILTER_LEN 100
#define MIN_FILTER_LEN 1

#define ALLOW_VARIABLE_DIM_COMBINED_STREAMS 1


// transformations
enum {
  NONE,
  UPSAMPLE_HOLD,
  UPSAMPLE_SMOOTH,
  DOWNSAMPLE,
  DELTAS,
  DOUBLE_DELTAS,
  MULTIPLY,
  NORMALIZE,
  MEAN_SUB,
  ARMA,
  AFFINE,
  FILTER,
  OFFSET
};

enum {
  FRAMEMATCH_ERROR,
  FRAMEMATCH_REPEAT_LAST,
  FRAMEMATCH_REPEAT_FIRST,
  FRAMEMATCH_EXPAND_SEGMENTALLY,
  FRAMEMATCH_TRUNCATE_FROM_END,
  FRAMEMATCH_TRUNCATE_FROM_START,
  FRAMEMATCH_REPEAT_RANGE_AT_START, // not used yet
  FRAMEMATCH_REPEAT_RANGE_AT_END   // not used yet
};

enum {
  SEGMATCH_ERROR,
  SEGMATCH_TRUNCATE_FROM_END,
  SEGMATCH_REPEAT_LAST,
  SEGMATCH_WRAP_AROUND
};

// ftrcombo operation flags
enum {
  FTROP_NONE = 0,
  FTROP_ADD,
  FTROP_SUB,
  FTROP_MUL,
  FTROP_DIV
};


enum {
  FRAMEJUSTIFICATION_LEFT,
  FRAMEJUSTIFICATION_CENTER,
  FRAMEJUSTIFICATION_RIGHT
};

#endif
