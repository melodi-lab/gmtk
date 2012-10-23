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

#include <assert.h>

#include "machine-dependent.h"
#include "range.h"

// The ObservationStream provides a simple API to wrap around
// non-random access data sources like pipes and sockets.
// Just subclass ObservationStream and implement getNextFrame
// for the new type of stream (and add command line options 
// to instantiate it - aspect-oriented programming?) and the
// new stream type is supported by GMTK.
//
// Planned subtypes:
//   ASCIIStream  -   ASCII protocol over a FILE (eg pipe or stdin)
//   SocketStream -   Data streamed over IP (in what format?)
//   FilterStream -   ObservationStream wrapper for IIR, ARMA, etc transforms
//   FileStream   -   ObservationStream wrapper for seekable files (ObservationFile)

class ObservationStream {

 protected:

  unsigned nFloat;
  unsigned nInt;

  char const *contFeatureRangeStr;  // -frX
  Range      *contFeatureRange;
  char const *discFeatureRangeStr;  // -irX
  Range      *discFeatureRange;

  Data32 *frameData;                // data as read from source
  Data32 *logicalFrameData;         // -frX and -irX applied

#if 0
  // don't make sense for streams?
  char const *preFrameRangeStr;     // -preprX
  Range      *preFrameRange;
  char const *segRangeStr;          // -srX
  Range      *segRange;
#endif

 public:

  ObservationStream() {
    contFeatureRangeStr = NULL; contFeatureRange = NULL;
    discFeatureRangeStr = NULL; discFeatureRange = NULL;
    nFloat = 0; nInt = 0;
    frameData = NULL; logicalFrameData = NULL;
  }

  ObservationStream(unsigned nFloat, unsigned nInt, char const *contFeatureRangeStr_=NULL, char const *discFeatureRangeStr_=NULL)
    : nFloat(nFloat), nInt(nInt),
      contFeatureRangeStr(contFeatureRangeStr_), contFeatureRange(NULL),
      discFeatureRangeStr(discFeatureRangeStr_), discFeatureRange(NULL)
  {
    frameData = new Data32[nFloat+nInt];
    logicalFrameData = new Data32[numLogicalFeatures()];
  }

  virtual ~ObservationStream() {
    if (contFeatureRange) delete contFeatureRange;
    if (discFeatureRange) delete discFeatureRange;
    if (frameData) delete [] frameData;
    if (logicalFrameData) delete [] logicalFrameData;
  }

  virtual unsigned numContinuous() {return nFloat;}
  virtual unsigned numDiscrete()   {return nInt;}
  virtual unsigned numFeatures()   {return numContinuous() + numDiscrete();}

  virtual unsigned numLogicalContinuous() {
    if (!contFeatureRange && contFeatureRangeStr) {
      contFeatureRange = new Range(contFeatureRangeStr, 0, numContinuous());
      assert(contFeatureRange);
    }
    if (contFeatureRange) {
      return contFeatureRange->length();
    } else {
      return numContinuous();
    }
  }
  
  virtual unsigned numLogicalDiscrete()  {
    if (!discFeatureRange && discFeatureRangeStr) {
      discFeatureRange = new Range(discFeatureRangeStr, 0, numDiscrete());
      assert(discFeatureRange);
    }
    if (discFeatureRange) {
      return discFeatureRange->length();
    } else {
      return numDiscrete();
    }
  }

  virtual unsigned numLogicalFeatures() {return numLogicalContinuous() + numLogicalDiscrete();}

  // Return the next frame from the stream.
  // May wait for data, may need a timeout.

  virtual bool EOS() {return true;} // true iff no more data available in the stream

  // Returns the next frame, or NULL to indicate the end of a segment.
  // There may be more segments, so call EOS() after getting a NULL from getNextFrame()

  virtual Data32 const *getNextFrame() {return NULL;}

  virtual Data32 const *getNextLogicalFrame() { 
    if (getNextFrame()) {
      assert(frameData && logicalFrameData);
      unsigned nlc = numLogicalContinuous();
      unsigned nld = numLogicalDiscrete();
#if 0
printf("ObsStr GLF:");
for (unsigned i=0; i<nlc; i+=1) printf(" %f", ((float *)frameData)[i]);
for (unsigned i=0; i<nld; i+=1) printf(" %u", ((unsigned *)frameData)[nlc+i]);
printf("\n");
//fprintf(stderr, "nlc %u   nlf %u\n", nlc, nlf);
#endif
      for (unsigned i=0; i < nlc; i+=1) {
	unsigned srcIdx = contFeatureRange ? contFeatureRange->index(i) : i;
//printf("%u -> %u %f\n", i, srcIdx,((float *)frameData)[srcIdx]);
	logicalFrameData[i] = frameData[srcIdx];
      }
      for (unsigned i = 0; i < nld; i+=1) {
	unsigned srcIdx = discFeatureRange ? discFeatureRange->index(i) : i;
//printf("%u -> %u %u\n", nlc+i, srcIdx, ((unsigned *)frameData)[nlc+i]);
//fprintf(stderr,"gNLF input %u -> output %u\n", srcIdx, i);  
        logicalFrameData[nlc+i] = frameData[numContinuous() + srcIdx];
      }
#if 0
printf("ObsStr GLF:");
for (unsigned i=0; i<nlc; i+=1) printf(" %f", ((float *)logicalFrameData)[i]);
for (unsigned i=0; i<nld; i+=1) printf(" %u", ((unsigned *)logicalFrameData)[nlc+i]);
printf("\n");
#endif
      return logicalFrameData;
    }
    return NULL;
  }

};

#endif
