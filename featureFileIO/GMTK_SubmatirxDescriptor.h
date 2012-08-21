
/*
 * GMTK_SubmatirxDescriptor.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2012, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */

#ifndef GMTK_SUBMATRIXDESCRIPTOR_H
#define GMTK_SUBMATRIXDESCRIPTOR_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <assert.h>

// Packages up the "shape" of a subsection of an observation matrix.
// This is used to ask a Filter which frames it needs to produce
// the requested output frames and to describe the Filter input and
// output frame ranges.

class subMatrixDescriptor {

 public:

  unsigned firstFrame;            // frame sub-matrix starts at
  unsigned numFrames;             // # frames in sub-matrix
  unsigned historyFrames;         // # frames just there for history, i.e., not part of the
                                  //   desired output but necessary to compute it
  unsigned futureFrames;          // # frames just there for future
  unsigned numContinuous;
  unsigned numDiscrete;
  unsigned fullMatrixFrameCount;  // # frames in full matrix

  unsigned requestedFirst;        // requested filter output
  unsigned requestedCount;

  // Maintain a list of pre-allocated instances to minimize new and ctor calls.

  subMatrixDescriptor *next;
  static subMatrixDescriptor *freeList;

  static subMatrixDescriptor *getSMD() {
    subMatrixDescriptor *result;
    if (freeList) {
      result = freeList;
      freeList = freeList->next;
      result->next = NULL;
    } else {
      result = new subMatrixDescriptor();
    }
    assert(result);
    return result;
  }

  static subMatrixDescriptor *
  getSMD(unsigned firstFrame, unsigned numFrames,
	 unsigned historyFrames, unsigned futureFrames,
	 unsigned numContinuous, unsigned numDiscrete,
	 unsigned fullMatrixFrameCount, 
	 unsigned requestedFirst,
	 unsigned requestedCount,
	 subMatrixDescriptor *next = NULL)
  {
    subMatrixDescriptor *result  = getSMD();
    result->firstFrame           = firstFrame;
    result->numFrames            = numFrames;
    result->historyFrames        = historyFrames;
    result->futureFrames         = futureFrames;
    result->numContinuous        = numContinuous;
    result->numDiscrete          = numDiscrete;
    result->fullMatrixFrameCount = fullMatrixFrameCount;
    result->requestedFirst       = requestedFirst;
    result->requestedCount       = requestedCount;
    result->next                 = next;
    return result;
  }

  static void freeSMD(subMatrixDescriptor *&p) {
    if (p->next) freeSMD(p->next);
    p->next  = freeList;
    freeList = p;
    p = NULL;
  }


  subMatrixDescriptor() {next = NULL;}

  subMatrixDescriptor(unsigned firstFrame, unsigned numFrames,
		      unsigned historyFrames, unsigned futureFrames,
		      unsigned numContinuous, unsigned numDiscrete,
		      unsigned fullMatrixFrameCount, 
		      unsigned requestedFirst,
		      unsigned requestedCount,
		      subMatrixDescriptor *next = NULL)
    : firstFrame(firstFrame), numFrames(numFrames),
      historyFrames(historyFrames), futureFrames(futureFrames),
      numContinuous(numContinuous), numDiscrete(numDiscrete),
      fullMatrixFrameCount(fullMatrixFrameCount),
      requestedFirst(requestedFirst), requestedCount(requestedCount),
      next(next)
  {}

  void deallocate() {
    if (next) {
      next->deallocate();
      delete next;
    }
  }

  ~subMatrixDescriptor() {
    deallocate();
  }

}; 


#endif
