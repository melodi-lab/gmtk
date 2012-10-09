
/*
 * GMTK_MergeStream.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2011, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */

#ifndef GMTK_MERGESTREAM_H
#define GMTK_MERGESTREAM_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include <stdio.h>
using namespace std;

#include "machine-dependent.h"
#include "error.h"

#include "GMTK_Filter.h"
#include "GMTK_ObservationStream.h"


class MergeStream: public ObservationStream {
  ObservationStream **stream;
  unsigned nStreams;

  unsigned totalFeatures;

  unsigned *floatStart; // The i^th stream's float features start at mergedFrame[i]
  unsigned *intStart;   //   likewise for int features

  bool eos;

 public:

  MergeStream() : stream(NULL), nStreams(0), totalFeatures(0), 
                  floatStart(NULL), intStart(NULL) {}
  
  MergeStream(ObservationStream *stream[], unsigned nStreams)
    : nStreams(nStreams), eos(false)
  {
    if (!stream || nStreams == 0) {
      error("MergeStream: no input streams to merge!=\n");
    }
    this->stream = new ObservationStream *[nStreams];
    for (unsigned i=0; i < nStreams; i+=1) {
      assert(stream[i]);
      this->stream[i] = stream[i];
    }
    floatStart = new unsigned[nStreams];
    intStart   = new unsigned[nStreams];
    totalFeatures = 0;
    for (unsigned i=0; i < nStreams; i+=1) {
      floatStart[i] = totalFeatures;
      nFloat += stream[i]->numLogicalContinuous();
      totalFeatures += stream[i]->numLogicalContinuous();
    }
    for (unsigned i=0; i < nStreams; i+=1) {
      intStart[i] = totalFeatures;
      nInt += stream[i]->numLogicalDiscrete();
      totalFeatures += stream[i]->numLogicalDiscrete();
    }
    frameData = new Data32[totalFeatures];
  }

  ~MergeStream() {
    
    if (stream) {
      for (unsigned i=0; i < nStreams; i+=1)
	if (stream[i]) delete stream[i];
     delete[] stream;
    }
    if (floatStart) delete[] floatStart;
    if (intStart)   delete[] intStart;
  }

  bool EOS() {return eos;}

  Data32 const *getNextLogicalFrame() { return getNextFrame(); }

  Data32 const *getNextFrame() {
    for (unsigned i=0; i < nStreams; i+=1) {
      Data32 const *partialFrame = stream[i]->getNextLogicalFrame();
      if (!partialFrame) {
	if (i != 0) {
	  error("MergeStream: input stream %u ended before the other streams\n");
	}
	eos = stream[i]->EOS();
	for (i+=1; i < nStreams; i+=1) {
	  if (stream[i]->EOS() != eos) {
	    error("MergeStream: input streams disagree on number of segements\n");
	  }
	}
	return NULL;
      }
      unsigned nC = stream[i]->numLogicalContinuous();
      unsigned nD = stream[i]->numLogicalDiscrete();
      memcpy(frameData + floatStart[i], partialFrame     , nC * sizeof(Data32));
      memcpy(frameData + intStart[i],   partialFrame + nC, nD * sizeof(Data32));
    }
    return frameData;
  }

};

#endif

