
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
#include "debug.h"

#include "GMTK_FilterFile.h"
#include "GMTK_FileSource.h"
#include "GMTK_Filter.h"
#include "GMTK_MergeFile.h"

#define BUFFER_SIZE (1024*1024)
  // each file provides its data for the requested submatrix,
  // transformed and "ready to eat." So the FileSource
  // just manages assembling them to satisfy the loadFrames()
  // calls, prefetching/caching. For archipelagos, each thread
  // gets its own FileSource (all aimed at the same files, of course).
FileSource::FileSource(ObservationFile *file, 
		       unsigned bufferSize,
		       unsigned startSkip, unsigned endSkip,
		       int justificationMode) 
{
  initialize(file, bufferSize, startSkip, endSkip, justificationMode);
}


#if 0
static
void
dumpFrames(Data32 *buf, unsigned nFloat, unsigned nInt,
	   unsigned nFrames)
{
  for (unsigned j=0; j < nFrames; j+=1) {
    float *fp = (float *)buf;
    printf("%3u: ", j);
    for (unsigned i=0; i < nFloat; i+=1) {
      printf(" %f", *(fp++));
    }
    unsigned *ip = (unsigned *)fp;
    for (unsigned i=0; i < nInt; i+=1) {
      printf(" %u", *(ip++));
    }
    printf("\n");
  }
}
#endif
    

void 
FileSource::initialize(ObservationFile *file, 
		       unsigned bufferSize,
		       unsigned startSkip, unsigned endSkip,
		       int justificationMode) 
{
  assert(file);
  assert( 0 <= justificationMode && justificationMode <= FRAMEJUSTIFICATION_RIGHT );
  this->_startSkip = startSkip;
  this->_endSkip = endSkip;
  this->justificationMode = justificationMode;
  justificationOffset = 0;
  _minPastFrames = 0;
  _minFutureFrames = 0;
  this->file = file;
  cookedBuffer = new Data32[bufferSize];
  if (!cookedBuffer) {
    error("ERROR: FileSource::intialize: failed to allocate frame buffer");
  }
  this->bufferSize = bufferSize;
  numBufferedFrames = 0;
  segment = -1;       // no openSegment() call yet
  bufStride = file->numFeatures();
  bufferFrames = bufferSize / bufStride;
}


bool 
FileSource::openSegment(unsigned seg) {
  assert(file);
  if (seg >= numSegments()) {
    error("ERROR: FileSource::openSegment: requested segment %u, but only up to %u are available\n", 
	  seg, numSegments()-1);
  }

  this->segment = seg;
  bool success = file->openSegment(seg);
  _numCacheableFrames = file->numLogicalFrames();
  if (_numCacheableFrames < _startSkip + _endSkip) {
    error("ERROR: segment %u has only %u frames, but -startSkip %u and -endSkip %u requires at least %u frames", seg, _numCacheableFrames, _startSkip, _endSkip, _startSkip + _endSkip + 1);
  }
  _numFrames = _numCacheableFrames;
  _numFrames -= _startSkip; // reserve frames at start of segment
  _numFrames -= _endSkip;   // reserve frames at end of segment

#ifndef JEFFS_STRICT_DEBUG_OUTPUT_TEST
  infoMsg(IM::ObsFile,IM::Low,"%u / %u input usable/cacheable frames in segment %d (unjustified)\n",
	  _numFrames, _numCacheableFrames, seg);
#endif
  
  justificationOffset = 0;
  numBufferedFrames = 0;
  return success;
}


void 
FileSource::justifySegment(unsigned numUsableFrames) {
  if (numUsableFrames > _numFrames) {
    error("ERROR: FileSource::justifySegment: numUsableFrames (%u) must not be larger than the number of available frames (%u)",
	  numUsableFrames, _numFrames);
  }
  assert(segment >= 0);
  switch (justificationMode) {
  case FRAMEJUSTIFICATION_LEFT:
    justificationOffset = 0;
    break;
  case FRAMEJUSTIFICATION_CENTER:
    justificationOffset = (_numFrames - numUsableFrames) / 2;
    break;
  case FRAMEJUSTIFICATION_RIGHT:
    justificationOffset = _numFrames - numUsableFrames;
    break;
  default:
    error("ERROR: FileSource::justifySegment: unknown justification mode %d", justificationMode);
  }
  infoMsg(IM::ObsFile,IM::Low,"justification mode %u offset = %u\n", justificationMode, justificationOffset);
  _numFrames = numUsableFrames;
}


Data32 const*
FileSource::loadFrames(unsigned bufferIndex, unsigned first, unsigned count) {
  assert(0 <= bufferIndex && bufferIndex < bufferFrames);
  assert(0 < count && count < bufferFrames-bufferIndex);
  
  unsigned buffOffset = bufferIndex * bufStride;

  if (first > _numCacheableFrames || first + count > _numCacheableFrames) {
    error("ERROR: FileSource::loadFrames: requested frames [%u,%u), but %u is the last available frame\n", first, first+count, _numCacheableFrames);
  }
    
  Data32 const *fileBuf = file->getLogicalFrames(first,count);
  assert(fileBuf);
  Data32 *dst = cookedBuffer + buffOffset;

  memcpy((void *)dst, (const void *)fileBuf, file->numLogicalFeatures() * count * sizeof(Data32));
#if 0
  warning("fetching [%u,%u) -> %u", first, first+count, buffOffset);
for (unsigned x=0; x < count; x+=1) {
  printf("%3u:", x+first);
  for (unsigned c=0; c < numContinuous(); c+=1) printf(" %f", *((float    *)(cookedBuffer+buffOffset+c+x*bufStride)));
  for (unsigned d=0; d < numDiscrete()  ; d+=1) printf(" %u", *((unsigned *)(cookedBuffer+buffOffset+d+x*bufStride+numContinuous())));
  printf("  @ %u\n",buffOffset+x*bufStride );
}
#endif
  return cookedBuffer + buffOffset;
}


#define window 20
#define delta 5

// FIXME - have to ensure startSkip frames before and endSkip frames after
// every request ?
Data32 const *
FileSource::loadFrames(unsigned first, unsigned count) {
  unsigned preFirst;  // first frame # to prefetch
  unsigned preCount;  // # of frames in prefetch request

#if 0
unsigned requestedFirst = first; 
unsigned requestedCount = count;
#endif

  // FIXME - adjust first+count checking for endSkip & justification
  if (first + count > _numFrames) {
    error("ERROR: FileSource::loadFrames: requested frames [%u,%u), but only %u frames are available", first, first+count, _numFrames);
  }


  first += _startSkip + justificationOffset;

  if (first >= _minPastFrames) {
#if 0
    // I think the _min{Past,Future}Frames adjustment needs to apply to pre{First,Count}
    first -= _minPastFrames;
    count += _minPastFrames;
#endif
  } else {
    error("ERROR: FileSource::loadFrames: requested frames [%u,%u), but there must be %u frames before the first requested frame",
	  first, first+count, _minPastFrames);
  }
  if (first + count - 1 + _minFutureFrames <= _numCacheableFrames) {
#if 0
    // I think the _min{Past,Future}Frames adjustment needs to apply to pre{First,Count}
    count += _minFutureFrames;
#endif
  } else {
    error("ERROR: FileSource::loadFrames: requested frames [%u,%u), but there must be %u frames after the last requested frame",
	  first, first+count, _minFutureFrames);
  }

#if 0
printf("loadFrames [%u,%u) -> [%u,%u)  requires [%u,%u)\n", 
       requestedFirst, requestedFirst + requestedCount,
       first, first+count,
       first - _minPastFrames, first + count + _minPastFrames + _minFutureFrames);
#endif
  if (numBufferedFrames == 0) {
    // cache empty

    // FIXME - move this above if
    // also need to check against count + _minFutureFrames + _minPastFrames
    if (bufferFrames < count) {
      error("ERROR: FileSource::loadFrames: requested %u frames, but buffer can only hold %u", count, bufferFrames);
    }
    // FIXME - first + count + _minFutureFrames ?
    if (first + count > _numCacheableFrames) {
      error("ERROR: FileSource::loadFrames: requested frames %u to %u, but there are only %u frames available", first, first+count-1, _numFrames);
    }

    if (count + _minPastFrames + _minFutureFrames + 2 * window < bufferFrames) {
      preFirst = ( first > window + _minPastFrames ) ? first - window - _minPastFrames : 0;
      preCount = ( preFirst + count + _minPastFrames + _minFutureFrames + 2 * window < _numCacheableFrames ) ?
        count + _minPastFrames + _minFutureFrames + 2 * window : _numCacheableFrames - preFirst;
    } else {
      assert(first >= _minPastFrames);
      preFirst = first - _minPastFrames;
      preCount = count + _minPastFrames + _minFutureFrames;
      assert(preFirst + preCount <= _numCacheableFrames);
    }
    assert(preCount > 0);
    firstBufferedFrameIndex = (bufferFrames - preCount) / 2;
    firstBufferedFrame = preFirst;
    numBufferedFrames = preCount;
    if (firstBufferedFrameIndex + numBufferedFrames > bufferFrames) {
      error("ERROR: FileSource:loadFrames:  attempted to load %u frames at index %u, which overflows the frame buffer", numBufferedFrames, firstBufferedFrameIndex*bufStride);
    }
    Data32 const *frames = loadFrames(firstBufferedFrameIndex, preFirst, preCount);
    assert(frames == cookedBuffer + firstBufferedFrameIndex * bufStride);
#if 0
warning("CACHE EMPTY [%u,%u) -> [%u,%u)  cached [%u,%u) @ %u",
	requestedFirst, requestedFirst + count,
	first, first+count, 
	firstBufferedFrame, firstBufferedFrame+numBufferedFrames,
	(firstBufferedFrameIndex + (first - firstBufferedFrame))*bufStride);

for (unsigned x=preFirst; x < preFirst+preCount; x+=1) {
  unsigned buffOffset = firstBufferedFrameIndex * bufStride;
  printf("%3u:", x);
  for (unsigned c=0; c < numContinuous(); c+=1) printf(" %f", *((float    *)(cookedBuffer+buffOffset+c+x*bufStride)));
  for (unsigned d=0; d < numDiscrete()  ; d+=1) printf(" %u", *((unsigned *)(cookedBuffer+buffOffset+d+x*bufStride+numContinuous())));
  printf("  @ %u\n",buffOffset+x*bufStride );
}
#endif
#if 0
for (unsigned frm=0; frm < preCount; frm+=1) {
  printf("   C %u @ %u :", preFirst+frm, firstBufferedFrameIndex+frm);
  float *ppp = (float *)frames;
  for (unsigned j=0; j < 2; j+=1) {
    printf(" %f", *(ppp++));
  }
  printf("\n");
}
#endif
    return frames + (first - preFirst) * bufStride;
  }
  // FIXME - i think the -1 should be gone
  if (firstBufferedFrame + _minPastFrames <= first && 
      first + count - 1 + _minFutureFrames <= firstBufferedFrame + numBufferedFrames)
  { 
    // cache hit
    
    if (first - _minPastFrames - firstBufferedFrame <= delta && firstBufferedFrame > 0) {  // prefetch backward?
      preFirst = firstBufferedFrame > window ? firstBufferedFrame - window : 0;
      preCount = firstBufferedFrame - preFirst;
      if (firstBufferedFrameIndex + numBufferedFrames + preCount < bufferFrames) { // do prefetch
	firstBufferedFrame = preFirst;
	firstBufferedFrameIndex -= preCount;
#if 0
warning("CACHE HIT< [%u,%u)  cached [%u,%u)", first,first+count, firstBufferedFrame, firstBufferedFrame+numBufferedFrames);
#endif
#if 1
        (void)loadFrames(firstBufferedFrameIndex, preFirst, preCount);
#else
	Data32 const *frames = loadFrames(firstBufferedFrameIndex, preFirst, preCount);
#endif
	numBufferedFrames += preCount;
      } 
    } else if (firstBufferedFrame + numBufferedFrames - (first + count + _minFutureFrames) <= delta && 
	firstBufferedFrame + numBufferedFrames < _numCacheableFrames) 
    {  // prefetch forward?
      preFirst = firstBufferedFrame + numBufferedFrames;
      preCount = firstBufferedFrame + numBufferedFrames + window < _numCacheableFrames ? window : _numCacheableFrames - preFirst;
      if (firstBufferedFrameIndex + numBufferedFrames + preCount < bufferFrames) { // do prefetch
	// FIXME - does the above need a -1 ?
#if 0
warning("PREFETCH > [%u,%u)", preFirst, preFirst+preCount);
#endif
#if 1
        (void) loadFrames(firstBufferedFrameIndex+numBufferedFrames, preFirst, preCount);
#else
	Data32 const *frames = loadFrames(firstBufferedFrameIndex+numBufferedFrames, preFirst, preCount);
#endif
	numBufferedFrames += preCount;
#if 0
warning("CACHE HIT> [%u,%u)  cached [%u,%u)", first,first+count, firstBufferedFrame, firstBufferedFrame+numBufferedFrames);
#endif
	return cookedBuffer + (firstBufferedFrameIndex + (first - firstBufferedFrame)) * bufStride;
      } 
    } else {
    // no prefetch
#if 0
warning("CACHE HIT [%u,%u)  cached [%u,%u) @ %u", 
	first,first+count, 
	firstBufferedFrame, firstBufferedFrame+numBufferedFrames, 
	(firstBufferedFrameIndex + first - firstBufferedFrame)*bufStride);
#endif
    }
    assert(first - firstBufferedFrame >= _minPastFrames);
    return cookedBuffer + (firstBufferedFrameIndex + (first - firstBufferedFrame)) * bufStride;
  }

  // cache miss
#if 0
warning("CACHE MISS [%u,%u)  cached [%u,%u)", first,first+count, firstBufferedFrame, firstBufferedFrame+numBufferedFrames);
#endif
  if (count + _minPastFrames + _minFutureFrames + 2 * window < bufferFrames) {
    preFirst = ( first > window + _minPastFrames ) ? first - window - _minPastFrames : 0;
    preCount = ( preFirst + count + _minPastFrames + _minFutureFrames + 2 * window < _numCacheableFrames ) ?
      count + _minPastFrames + _minFutureFrames + 2 * window : _numCacheableFrames - preFirst;
  } else {
    assert(first >= _minPastFrames);
    preFirst = first - _minPastFrames;
    preCount = count + _minPastFrames + _minFutureFrames;
    assert(preFirst + preCount <= _numCacheableFrames);
  }
  firstBufferedFrameIndex = (bufferFrames - preCount) / 2;
  firstBufferedFrame = preFirst;
  numBufferedFrames = preCount;
#if 0
warning("CACHE FLUSH [%u,%u)  cached [%u,%u)", first, first+count, firstBufferedFrame, firstBufferedFrame+numBufferedFrames);
#endif
  Data32 const *frames = loadFrames(firstBufferedFrameIndex, preFirst, preCount);
  return frames + (first - preFirst) * bufStride;

#if 0
  if (bufferFrames < count) {
    error("ERROR: FileSource::loadFrames: requested %u frames, but buffer can only hold %u", count, bufferFrames);
  }
  firstBufferedFrameIndex = (bufferFrames - count) / 2;
  firstBufferedFrame = first;
  numBufferedFrames = count;
  return loadFrames(firstBufferedFrameIndex, first, count);
#endif
}


unsigned 
FileSource::numContinuous() {
  assert(file);
  return file->numLogicalContinuous();
}


unsigned 
FileSource::numDiscrete() {
  assert(file);
  return file->numLogicalDiscrete();
}


unsigned 
FileSource::numFeatures() {
  return numContinuous() + numDiscrete();
}


unsigned 
FileSource::stride() {
  return numFeatures();
}


// loadFrames() handles -justification and -startSkip

float *const 
FileSource::floatVecAtFrame(unsigned f) {
  assert(0 <= f && f < numFrames());
  Data32 const * buf = loadFrames(f, 1);
#if 0
  printf("%3u:", f);
  for (unsigned i=0; i < numContinuous(); i+=1) printf(" %f", ((float *)buf)[i]);
  for (unsigned i=numContinuous(); i < numFeatures(); i+=1) printf(" %u", ((unsigned *)buf)[i]);
#endif
  return(float *const) buf;
}


unsigned *const 
FileSource::unsignedVecAtFrame(unsigned f) {
  assert(0 <=f && f < numFrames());
#if 0
unsigned *up = (unsigned *)(loadFrames(f,1)+_numContinuous);
printf("uVec(%3u):", f);
 for (unsigned i=0; i < (unsigned)_numDiscrete; i+=1) printf(" %u", up[i]);
printf("\n"); 
#endif
  return (unsigned *)(loadFrames(f,1) + numContinuous());
}


unsigned &
FileSource::unsignedAtFrame(const unsigned frame, const unsigned feature) {
  assert (0 <= frame && frame < numFrames());
  assert (feature >= numContinuous()
	  &&
	  feature <  numFeatures());
#if 0
unsigned *up = (unsigned *)(loadFrames(frame,1)+_numContinuous);
printf("uVec(%3u,%u):", frame, feature);
 for (unsigned i=0; i < (unsigned)_numDiscrete; i+=1) printf(" %u", up[i]);
printf("\n"); 
#endif
  return *(unsigned*)(loadFrames(frame,1)+feature);
}


float *const 
FileSource::floatVecAtFrame(unsigned f, const unsigned startFeature) {
  assert (0 <= f && f < numFrames());
  float *result = (float*)(loadFrames(f,1) + startFeature);
  return result;
}

Data32 const * const
FileSource::baseAtFrame(unsigned f) {
  assert(0 <= f && f < numFrames());
  Data32 const * featuresBase = loadFrames(f,1);
#if 0
printf("baseAtFrame(%3u / %3u):", f, f + _startSkip);
float *fp = (float *)(featuresBase);
for(unsigned i=0; i < (unsigned)_numContinuous; i+=1) printf(" %f", fp[i]);
unsigned *up = (unsigned *)(featuresBase);
for(unsigned i=_numContinuous; i < stride(); i+=1) printf(" %u", up[i]);
printf("\n");
#endif
  return featuresBase;
}

bool 
FileSource::elementIsDiscrete(unsigned el) {
  return numContinuous() <= el && el < numFeatures();
}

bool
FileSource::active() {
  return this->segment > -1;
}

