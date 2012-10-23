
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


FileSource::FileSource(ObservationFile *file, 
		       unsigned windowBytes, unsigned deltaFrames, unsigned bufferSize, 
		       unsigned startSkip, unsigned endSkip,
		       int justificationMode, bool constantSpace)
{
  initialize(file, windowBytes, deltaFrames, bufferSize, startSkip, endSkip, justificationMode, constantSpace);
}


#if 0
// helpful for debugging

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
    
#define MEBIBYTE (1048576)

void 
FileSource::initialize(ObservationFile *file, 
		       unsigned windowBytes, unsigned deltaFrames, unsigned bufferSize, 
		       unsigned startSkip, unsigned endSkip,
		       int justificationMode, bool constantSpace) 
{
  assert(file);
  assert( 0 <= justificationMode && justificationMode <= FRAMEJUSTIFICATION_RIGHT );
  if (windowBytes > bufferSize * sizeof(Data32)) {
    error("ERROR: fileWindowSize (%u MB) must be less than fileBufferSize (%u MB)",
	  windowBytes/MEBIBYTE, bufferSize * sizeof(Data32) / MEBIBYTE);
  }
  // window size in frames
  this->window = windowBytes / (file->numLogicalFeatures() * sizeof(Data32));
  if (this->window < 2 * deltaFrames + 1) {
    error("ERROR: fileWindowSize (%u MB) must be at least %u MB to hold %u frame",
	  windowBytes / MEBIBYTE, file->numLogicalFeatures() * sizeof(Data32) / MEBIBYTE,
	  2 * deltaFrames + 1);
  }

  // Note: if the number of buffered frames less than 2 delta, it's ambiguous whether
  //       we should prefetch forward or backward. That's a goofy case thoug...
  this->delta = deltaFrames;
  this->_startSkip = startSkip;
  this->_endSkip = endSkip;
  this->constantSpace = constantSpace;
  this->justificationMode = justificationMode;
  justificationOffset = 0;
  _minPastFrames = 0;
  _minFutureFrames = 0;
  this->file = file;
  if (bufferSize > 0) {
    cookedBuffer = new Data32[bufferSize];
    if (!cookedBuffer) {
      error("ERROR: FileSource::intialize: failed to allocate frame buffer");
    }
  } else {
    cookedBuffer = NULL;
  }
  this->bufferSize = bufferSize;
  bufStride = file->numLogicalFeatures();
  bufferFrames = bufferSize / bufStride;
  if (bufferFrames < 1 && bufferSize > 0) {
    unsigned minSize = bufStride * sizeof(Data32) / MEBIBYTE;
    minSize = (minSize < 1) ? 1 : minSize;
    error("ERROR: file buffer size must be at least %u MB", minSize);
  }

  infoMsg(IM::ObsFile, IM::Low, "FileSource window = %u frames, %u MiB\n"
	  "  cookedBuffer = %u frames, %u MiB\n", this->window, windowBytes/MEBIBYTE,
	  bufferFrames, bufferSize * sizeof(Data32) / MEBIBYTE);
  numBufferedFrames = 0;
  segment = -1; // no openSegment() call yet

  numContinuousFeatures = file->numLogicalContinuous();
  numDiscreteFeatures = file->numLogicalDiscrete();
  _numFeatures = numContinuousFeatures + numDiscreteFeatures;
  if (_numFeatures == 0) {
    error("ERROR: No features (continuous or discrete) were selected.  Check the feature ranges.");
  }
}


bool 
FileSource::openSegment(unsigned seg) {
  assert(file);
  if (seg >= numSegments()) {
    error("ERROR: FileSource::openSegment: requested segment %u, but only up to %u are available\n", 
	  seg, numSegments()-1);
  }

  this->segment = seg;
  bool success = file->openLogicalSegment(seg);

  _numCacheableFrames = file->numLogicalFrames();  // the file handles -gpr, so this is what's left after that
  if (_numCacheableFrames < _startSkip + _endSkip) {
    error("ERROR: segment %u has only %u frames, but -startSkip %u and -endSkip %u requires at least %u frames", 
	  seg, _numCacheableFrames, _startSkip, _endSkip, _startSkip + _endSkip + 1);
  }
  _numFrames  = _numCacheableFrames;
  _numFrames -= _startSkip; // reserve frames at start of segment
  _numFrames -= _endSkip;   // reserve frames at end of segment

  // The modular debugging tests require the output to match the output 
  // from before the O(1) observation code, so infoMsg() calls in this
  // code will cause the tests to fail. So defining JEFFS_STRICT_DEBUG_OUTPUT_TEST
  // turns off the infoMsg() to allow the tests to pass.
#ifndef JEFFS_STRICT_DEBUG_OUTPUT_TEST
  infoMsg(IM::ObsFile,IM::Low,"%u / %u input accessable/cacheable frames in segment %d (unjustified)\n",
	  _numFrames, _numCacheableFrames, seg);
#endif
  
  justificationOffset = 0;  // default to left justification until justifySegment() is called
  numBufferedFrames = 0;    // empty the cache for the new segment

  if (!constantSpace) {     // load the entire segment
    if (_numCacheableFrames > bufferFrames) { // need to enlarge the buffer
      if (cookedBuffer) delete [] cookedBuffer;
      bufferFrames = _numCacheableFrames;
      bufferSize = numFeatures() * _numCacheableFrames;
      cookedBuffer = new Data32[bufferSize];
      if (!cookedBuffer) {
	error("ERROR: FileSource::openSegment: failed to allocate %u frame buffer for segment %u",
	      bufferSize, seg);
      }
    }
    unsigned bytesPerFrame = numFeatures() * sizeof(Data32);
    unsigned framesPerGulp;
    if (bytesPerFrame > DEFAULT_BUFFER_SIZE) {
      // The frames are extremely large. Try to read them all in in one go
      framesPerGulp = _numCacheableFrames;
    } else {
      // The frames are reasonably sized - read them incrementally
      framesPerGulp = DEFAULT_BUFFER_SIZE / bytesPerFrame;
    }
    // load all the frames
    unsigned remainder = _numCacheableFrames % framesPerGulp;
    if (remainder > 0)
      (void) loadFrames(0, 0, remainder);
    firstBufferedFrame = 0;
    firstBufferedFrameIndex = 0;
    numBufferedFrames = remainder;
    for (unsigned frame=remainder; frame < _numCacheableFrames; frame+=framesPerGulp) {
      (void) loadFrames(frame, frame, framesPerGulp);
      numBufferedFrames += framesPerGulp;
    }
  }
  return success;
}


void 
FileSource::justifySegment(unsigned numUsableFrames) {
  if (numUsableFrames > _numFrames) {
    error("ERROR: FileSource::justifySegment: numUsableFrames (%u) must not be "
	  "larger than the number of available frames (%u)",
	  numUsableFrames, _numFrames);
  }
  assert(segment >= 0);
  // set justificationOffset to the # of frames required to achieve
  // the specified justificaton
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
    error("ERROR: FileSource::justifySegment: unknown justification mode %d", 
	  justificationMode);
  }
  infoMsg(IM::ObsFile,IM::Low,"justification mode %u offset = %u\n", 
	  justificationMode, justificationOffset);
  _numFrames = numUsableFrames - _startSkip - _endSkip;
}


#if 0
static unsigned cachemiss = 0;
static unsigned cachehit  = 0;
#endif


// FIXME - check for and handle NaNs ?

// This loadFrames() method just loads the frames from the ObservationFile
// into the indicated position within the cookedBuffer. The other loadFrames()
// method below is the main workhorse. 
Data32 const*
FileSource::loadFrames(unsigned bufferIndex, unsigned first, unsigned count) {
  assert(0 <= bufferIndex && bufferIndex < bufferFrames);
  assert(0 < count && count <= bufferFrames-bufferIndex);
  
  unsigned buffOffset = bufferIndex * bufStride;

  assert(!(first > _numCacheableFrames || first + count > _numCacheableFrames));
    
  Data32 const *fileBuf = file->getLogicalFrames(first,count);
  assert(fileBuf);
  Data32 *dst = cookedBuffer + buffOffset;

  memcpy((void *)dst, (const void *)fileBuf, file->numLogicalFeatures() * count * sizeof(Data32));
#ifdef WARNING_ON_NAN
  if (!ObservationsAllowNan) {
    for (float *fp = (float *)dst, unsigned i=0; i < count; i+=1, fp += file->numLogicalFeatures()) {
      for (unsigned j=0; j < file->numLogicalContinuousFeatures(); j+=1) {
	if (isnan(fp[j])) {
	  error("ERROR: Found NaN or +/-INF at %u'th float in frame %u, segment %u\n",
                j, first+i, segmentNumber());
	}
      }
    }
  }
#endif
  return cookedBuffer + buffOffset;
}


// This method loads the requested frames (along with any necessary
// preceding or following frames) into the cookedBuffer if they are
// not already present and returns a pointer to the first requested
// frame. It tries to prefetch a full window of frames at once when
// the inference code gets close to the end of the buffered frames.
Data32 const *
FileSource::loadFrames(unsigned first, unsigned count) {
  unsigned preFirst;  // first frame # to prefetch
  unsigned preCount;  // # of frames in prefetch request

  if (first + count > _numFrames) {
    error("ERROR: FileSource::loadFrames: requested frames [%u,%u), but "
	  "only frames [0,%u) are available", 
	  first, first+count, _numFrames);
  }
  
  if (bufferFrames < count + _minFutureFrames + _minPastFrames) {
    error("ERROR: FileSource::loadFrames: requested %u frames, but buffer "
	  "can only hold %u", count + _minFutureFrames + _minPastFrames, 
	  bufferFrames);
  }
  first += _startSkip + justificationOffset; // adjust for -startSkip and -justification

  if (numBufferedFrames == 0) {
    // the cookedBuffer is empty - load the first window of frames into the middle of the buffer
    if (count + _minPastFrames + _minFutureFrames + window < bufferFrames) {
#if 0
      infoMsg(IM::ObsFile, IM::Giga, "%u+%u+%u frames + %u window fits in %u buffer capacity\n", 
	      _minPastFrames,count,_minFutureFrames, window, bufferFrames);
#endif
      // if there's enough room, prefetch a window's worth of frames before 
      // and after the requested range

      preFirst = ( first > window/2 + _minPastFrames ) ? first - window/2 - _minPastFrames : 0;
      preCount = ( preFirst + count + _minPastFrames + _minFutureFrames + window < _numCacheableFrames ) ?
        count + _minPastFrames + _minFutureFrames + window : _numCacheableFrames - preFirst;
    } else {
#if 0
      infoMsg(IM::ObsFile, IM::Giga, "%u+%u+%u frames + %u window > %u buffer capacity\n", 
	      _minPastFrames,count,_minFutureFrames, window, bufferFrames);
#endif
      // otherwise, just prefetch the requested range and the minimum
      // number of preceding and subsequent frames
      assert(first >= _minPastFrames);
      preFirst = first - _minPastFrames;
      preCount = count + _minPastFrames + _minFutureFrames;
      assert(preFirst + preCount <= _numCacheableFrames);
    }
    assert(preCount > 0);
    // put the frames in the middle of the cookedBuffer so it can grow
    // in either direction, since we don't know if inference is doing a
    // forward or backward pass
    firstBufferedFrameIndex = (bufferFrames - preCount) / 2;  
    firstBufferedFrame = preFirst;
    numBufferedFrames = preCount;
#if 0
    infoMsg(IM::ObsFile, IM::Giga, "empty cache, fetching [%u,%u)@%u for [%u,%u)\n",
	    preFirst, preFirst+preCount, firstBufferedFrameIndex, first, first+count);
#endif
    if (firstBufferedFrameIndex + numBufferedFrames > bufferFrames) {
      error("ERROR: FileSource:loadFrames:  attempted to load %u frames at index %u, "
	    "which overflows the frame buffer", numBufferedFrames, firstBufferedFrameIndex*bufStride);
    }
    Data32 const *frames = loadFrames(firstBufferedFrameIndex, preFirst, preCount);
    assert(frames == cookedBuffer + firstBufferedFrameIndex * bufStride);
    // we may have loaded more than was asked for, so just return the 
    // requested range
    return frames + (first - preFirst) * bufStride;
  } // END of empty cache case


  if (firstBufferedFrame + _minPastFrames <= first && 
      first + count + _minFutureFrames <= firstBufferedFrame + numBufferedFrames)
  { 
    // cache hit

    assert(firstBufferedFrame <= first); // A
    //assert(first < firstBufferedFrame + numBufferedFrames); // redundant w/ B

    //assert(firstBufferedFrame < first + count);  // redundant w/ A
    assert(first + count <= firstBufferedFrame + numBufferedFrames); // B
    
    assert(count > 0); // redundancies require count > 0

#if 0
cachehit+=1;
if (cachehit % 80000 == 0) 
infoMsg(IM::ObsFile, IM::Giga, "frames [%7u,%7u) cache hit %7u  cache miss %7u\n",
	first, first+count, cachehit, cachemiss);
#endif
    if (first - _minPastFrames - firstBufferedFrame < delta && firstBufferedFrame > 0) { 
      // prefetch backward if within delta frames of the first buffered frame
      preFirst = firstBufferedFrame > window ? firstBufferedFrame - window : 0;
      preCount = firstBufferedFrame - preFirst;
      if (preCount <= firstBufferedFrameIndex) { // do prefetch
#if 0
	infoMsg(IM::ObsFile, IM::Giga, "prefetch <  [%u,%u) + [%u,%u) for [%u,%u)\n",
		firstBufferedFrame, firstBufferedFrame + numBufferedFrames,
		preFirst, preFirst + preCount,
		first, first + count);
#endif
	firstBufferedFrame = preFirst;
	firstBufferedFrameIndex -= preCount;
        (void)loadFrames(firstBufferedFrameIndex, preFirst, preCount);
	numBufferedFrames += preCount;
      } 
    } else if (firstBufferedFrame + numBufferedFrames - (first + count + _minFutureFrames) < delta && 
	firstBufferedFrame + numBufferedFrames < _numCacheableFrames) 
    { 
      // prefetch forward if within delta frames of the last buffered frame
      preFirst = firstBufferedFrame + numBufferedFrames;
      preCount = firstBufferedFrame + numBufferedFrames + window < _numCacheableFrames ? 
	window : _numCacheableFrames - preFirst;
      if (firstBufferedFrameIndex + numBufferedFrames + preCount < bufferFrames) { // do prefetch
#if 0
	infoMsg(IM::ObsFile, IM::Giga, "prefetch >  [%u,%u) + [%u,%u)@%u for [%u,%u)\n",
		firstBufferedFrame, firstBufferedFrame + numBufferedFrames,
		preFirst, preFirst + preCount, firstBufferedFrameIndex+numBufferedFrames,
		first, first + count);
#endif
        (void) loadFrames(firstBufferedFrameIndex+numBufferedFrames, preFirst, preCount);
	numBufferedFrames += preCount;
      }
#if 0 
      else
	infoMsg(IM::ObsFile, IM::Giga, "no fetch >  [%u,%u) + [%u,%u)@%u for [%u,%u) would overflow %u\n",
		firstBufferedFrame, firstBufferedFrame + numBufferedFrames,
		preFirst, preFirst + preCount, firstBufferedFrameIndex+numBufferedFrames,
		first, first + count, bufferFrames);
#endif
    } else {
      // no prefetch
    }
    assert(first - firstBufferedFrame >= _minPastFrames);
    return cookedBuffer + (firstBufferedFrameIndex + (first - firstBufferedFrame)) * bufStride;
  } // END of cache hit case

  // cache miss - drop the current frames and load in the new
#if 0
infoMsg(IM::ObsFile, IM::Giga, "cache miss %u\n", ++cachemiss);
#endif
  if (count + _minPastFrames + _minFutureFrames + window < bufferFrames) {
    preFirst = ( first > window/2 + _minPastFrames ) ? first - window/2 - _minPastFrames : 0;
    preCount = ( preFirst + count + _minPastFrames + _minFutureFrames + window < _numCacheableFrames ) ?
      count + _minPastFrames + _minFutureFrames + window : _numCacheableFrames - preFirst;
#if 0
    infoMsg(IM::ObsFile, IM::Giga, "cache flush: %u+%u+%u frames + %u window fits in %u buffer capacity\n", 
	    _minPastFrames,count,_minFutureFrames, window, bufferFrames);
#endif
  } else {
    assert(first >= _minPastFrames);
    preFirst = first - _minPastFrames;
    preCount = count + _minPastFrames + _minFutureFrames;
    assert(preFirst + preCount <= _numCacheableFrames);
#if 0
    infoMsg(IM::ObsFile, IM::Giga, "cache flush: %u+%u+%u frames + %u window > %u buffer capacity\n", 
	    _minPastFrames,count,_minFutureFrames, window, bufferFrames);
#endif
  }
#if 0
  infoMsg(IM::ObsFile, IM::Giga, "cache miss on [%u,%u)/%u for [%u,%u)   loading [%u,%u)@%u\n",
	  firstBufferedFrame, firstBufferedFrame+numBufferedFrames, delta,
	  first, first+count, preFirst, preFirst+preCount, (bufferFrames - preCount) / 2);
#endif
  firstBufferedFrameIndex = (bufferFrames - preCount) / 2;
  firstBufferedFrame = preFirst;
  numBufferedFrames = preCount;
  Data32 const *frames = loadFrames(firstBufferedFrameIndex, preFirst, preCount);
  return frames + (first - preFirst) * bufStride;
}


bool
FileSource::active() {
  return this->segment > -1;
}


// The following all just use loadFrames() to ensure the requested frame
// is present in the cookedBuffer.

float *const 
FileSource::floatVecAtFrame(unsigned f) {
  assert(0 <= f && f < numFrames());
  Data32 const * buf = loadFrames(f, 1);
  return(float *const) buf;
}


unsigned *const 
FileSource::unsignedVecAtFrame(unsigned f) {
  assert(0 <=f && f < numFrames());
  return (unsigned *)(loadFrames(f,1) + numContinuous());
}


unsigned &
FileSource::unsignedAtFrame(const unsigned frame, const unsigned feature) {
  assert (0 <= frame && frame < numFrames());
  assert (feature >= numContinuous() && feature <  numFeatures());
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
  return featuresBase;
}

