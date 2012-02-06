
/*
 * GMTK_FileSource.cc
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2011, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <string.h>

#include "error.h"
#include "general.h"
#include "machine-dependent.h"

#include "GMTK_ObservationFile.h"
#include "GMTK_FileSource.h"

#define BUFFER_SIZE (1024*1024)
  // each file provides its data for the requested submatrix,
  // transformed and "ready to eat." So the FileSource
  // just manages assembling them to satisfy the loadFrames()
  // calls, prefetching/caching. For archipelagos, each thread
  // gets its own FileSource (all aimed at the same files, of course).
FileSource::FileSource(unsigned _nFiles, ObservationFile *file[], 
		       char const *_globalFrameRangeStr, unsigned startSkip,
		       unsigned endSkip) 
{
  initialize(_nFiles, file, _globalFrameRangeStr, startSkip, endSkip);
}

void 
FileSource::initialize(unsigned nFiles, ObservationFile *file[], char const *globalFrameRangeStr, unsigned startSkip, unsigned endSkip) {
  this->nFiles = nFiles;
  this->globalFrameRangeStr = globalFrameRangeStr;
  this->globalFrameRange = NULL;
  this->_startSkip = startSkip;
  this->_endSkip = endSkip;
  this->file = new ObservationFile *[nFiles];
  for (unsigned i=0; i < nFiles; i+=1) {
    this->file[i] = file[i];
  }
  cookedBuffer = new Data32[BUFFER_SIZE];
  floatStart   = new Data32 *[nFiles];
  intStart     = new Data32 *[nFiles];
  this->segment = -1;
  // FIXME - FileSource should know about feature subranges unless the file format can do it better
  //  no it shouldn't - ObservationFile base class does naive feature subranges, subclasses
  //  can over-ride the *Logical* methods if they can do better
  unsigned offset = 0;
  for (unsigned i=0; i < nFiles; i+=1) {
    floatStart[i] = cookedBuffer + offset;
    offset += file[i]->numLogicalContinuous();
  }
  for (unsigned i=0; i < nFiles; i+=1) {
    intStart[i] = cookedBuffer + offset;
    offset += file[i]->numLogicalDiscrete();
  }
  bufStride = offset;
}


  // The number of available segments.
unsigned 
FileSource::numSegments() {
// FIXME - FileSource should know about segment subranges unless the file format can do it better
  return file[0]->numSegments();
}

  // Begin sourcing data from the requested segment.
  // Must be called before any other operations are performed on a segment.
bool 
FileSource::openSegment(unsigned seg) {
  // FIXME - FileSource should know about segment subranges unless the file format can do it better
  //  no it shouldn't - ObservationFile base class does naive feature subranges, subclasses
  //  can over-ride the *Logical* methods if they can do better
  bool success = true;

  // FIXME - error handling
  for (unsigned i=0; i < nFiles; i+=1) {
    success = success && file[i]->openLogicalSegment(seg);
  }
  this->segment = seg;
  if (globalFrameRange) delete globalFrameRange;
  // FIXME - max should be max logical frames over files
  globalFrameRange = new Range(globalFrameRangeStr, 0, file[0]->numLogicalFrames());
  return success;
}

  // The number of frames in the currently open segment.
unsigned 
FileSource::numFrames() {
  // FIXME - FileSource should know about frame subranges unless the file format can do it better
  //  no it shouldn't - ObservationFile base class does naive feature subranges, subclasses
  //  can over-ride the *Logical* methods if they can do better
  return globalFrameRange->length();
}

Data32 const*
FileSource::loadFrames(unsigned first, unsigned count) {
  // FIXME - startSkip and endSkip
  assert(count > 0);
  assert(count * numFeatures() < BUFFER_SIZE);
  for (unsigned int i=0; i < nFiles; i+=1) {
    Data32 const *fileBuf = file[i]->getLogicalFrames(first,count);
    assert(fileBuf);
    Data32 *dst = floatStart[i];
    Data32 const *src = fileBuf;
    unsigned srcStride = file[i]->numLogicalFeatures();
    for (unsigned f=0; f < count; f += 1, src += srcStride, dst += bufStride) 
    {
      memcpy((void *)dst, (const void *)src, file[i]->numLogicalContinuous() * sizeof(Data32));
    }
    src = fileBuf + file[i]->numLogicalContinuous();
    dst=intStart[i];
    for (unsigned f=0; f < count; f += 1, src += srcStride, dst += bufStride) 
    {
      memcpy((void *)dst, (const void *)src, file[i]->numLogicalDiscrete() * sizeof(Data32));
    }
  }
  // if requested frames are already in cookedBuffer
  //   return &cookedBuffer[index of first requested frame]
  // adjust (first,count) for prefetching -> (first',count')
  // adjust cookedBuffer destination
  // getFrames(first', count') from each file into cookedBuffer
  // return &cookedBuffer[destination of first]
  return cookedBuffer;
}



 // FIXME - FileSource should know about feature subranges unless the file format can do it better
//  no it shouldn't - ObservationFile base class does naive feature subranges, subclasses
//  can over-ride the *Logical* methods if they can do better

// The number of continuous, discrete, total features
unsigned 
FileSource::numContinuous() {
  unsigned sum = 0;
  for (unsigned i=0; i < nFiles; i+=1)
    sum += file[i]->numLogicalContinuous();
  return sum;
}

unsigned 
FileSource::numDiscrete() {
  unsigned sum = 0;
  for (unsigned i=0; i < nFiles; i+=1)
    sum += file[i]->numLogicalDiscrete();
  return sum;
}

unsigned 
FileSource::numFeatures() {
  return numContinuous() + numDiscrete();
}

// The number of Data32's between each frame
unsigned 
FileSource::stride() {
  return numFeatures();
}



float *const 
FileSource::floatVecAtFrame(unsigned f) {
  assert(0 <= f && f < numFrames());
  return (float *)loadFrames(f, 1);
}

unsigned *const 
FileSource::unsignedVecAtFrame(unsigned f) {
  assert(0 <=f && f < numFrames());
  return (unsigned *)loadFrames(f,1);
}

unsigned &
FileSource::unsignedAtFrame(const unsigned frame, const unsigned feature) {
  assert (0 <= frame && frame < numFrames());
  assert (feature >= numContinuous()
	  &&
	  feature <  numFeatures());
  return *(unsigned*)(loadFrames(frame,1)+feature);
}


float *const 
FileSource::floatVecAtFrame(unsigned f, const unsigned startFeature) {
  assert (0 <= f && f < numFrames());
  return (float*)(loadFrames(f,1) + startFeature);
}

Data32 const * const
FileSource::baseAtFrame(unsigned f) {
  assert(0 <= f && f < numFrames());
  return loadFrames(f,1);
}

bool 
FileSource::elementIsDiscrete(unsigned el) {
  return numContinuous() <= el && el < numFeatures();
}

bool
FileSource::active() {
  return this->segment > -1;
}
