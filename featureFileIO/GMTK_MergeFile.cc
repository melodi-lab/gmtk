
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <string.h>

#include "error.h"
#include "general.h"
#include "machine-dependent.h"
#include "debug.h"

#include "GMTK_MergeFile.h"
#include "GMTK_Filter.h"

// Adjust segment number seg according to -sdiffactX.
// Returns the physical segment number corresponding
// to the logical segment number seg.
unsigned 
MergeFile::adjustForSdiffact(unsigned fileNum, unsigned seg) {
  if (!sdiffact || sdiffact[fileNum] == SEGMATCH_TRUNCATE_FROM_END)
    return seg;
  if (sdiffact[fileNum] == SEGMATCH_REPEAT_LAST) {
    if (seg < file[fileNum]->numLogicalSegments())
      return seg;
    else
      return file[fileNum]->numLogicalSegments() - 1;
  } else if (sdiffact[fileNum] == SEGMATCH_WRAP_AROUND) {
    return seg % file[fileNum]->numLogicalSegments();
  } else if (sdiffact[fileNum] == SEGMATCH_ERROR) {
    return seg;
  } else 
    error("ERROR: MergeFile::adjustForSdiffact: unknown -sdiffact%u %d",
	  fileNum, sdiffact[fileNum]);
  return 0; // impossilbe to get here, but compiler warns about no return
}


// Returns the number of frames in the current segment
// after the -fdiffactX options are applied
static
unsigned
checkNumFrames(unsigned nFiles, ObservationFile *file[],
	       unsigned const *fdiffact)
{
  unsigned min_len=file[0]->numLogicalFrames();
  unsigned max_len=min_len;
  if(max_len == 0) 
    warning("WARNING: MergeFile::openSegment:checkNumFrames:  segment 0 is empty");

  bool got_error    = false;  // an input file insists all files have the same segment length
  bool got_truncate = false;  // an input file will shorten itself to match the others
  bool got_expand   = false;  // an input file will lengthen itself to match the others

  for(unsigned file_no=0; file_no < nFiles; file_no += 1) {
    unsigned len=file[file_no]->numLogicalFrames();
    if(len < min_len) min_len=len;
    if(len > max_len) max_len=len;

    if (fdiffact != NULL) {
      switch (fdiffact[file_no]) {
      case FRAMEMATCH_TRUNCATE_FROM_START:
      case FRAMEMATCH_TRUNCATE_FROM_END:
	got_truncate = true;
	break;
      case FRAMEMATCH_REPEAT_FIRST:
      case FRAMEMATCH_REPEAT_LAST:
      case FRAMEMATCH_EXPAND_SEGMENTALLY:
	got_expand = true;
	break;
      case FRAMEMATCH_ERROR:
	got_error = true;
	break;
      }
    } else {
      got_error = true; // FRAMEMATCH_ERROR is the default
    }
    if(got_truncate && got_expand)
      error("ERROR: MergeFile::openSegment:checkNumFrames: Cannot specify "
	    "both truncate and expand actions when segments have different lengths."
            "  Check the supplied and default -fdiff options.\n");
  }
  if(max_len == 0)
    error("ERROR: MergeFile::openSegment:checkNumFrames:  all segments have zero length");
  
#define EXPANSIVE(action,file_no) ( action && (action[file_no] == FRAMEMATCH_REPEAT_FIRST || \
                                               action[file_no] == FRAMEMATCH_REPEAT_LAST  || \
                                               action[file_no] == FRAMEMATCH_EXPAND_SEGMENTALLY) )

#define CONTRACTIVE(action,file_no) ( action && (action[file_no] == FRAMEMATCH_TRUNCATE_FROM_START || \
                                                 action[file_no] == FRAMEMATCH_TRUNCATE_FROM_END) )

  for(unsigned file_no=0; file_no < nFiles; ++file_no) {
    unsigned len=file[file_no]->numLogicalFrames();
    if( got_expand && len < max_len && !EXPANSIVE(fdiffact,file_no) ) {
      error("ERROR: observation file %u needs an -fdiffact%u that expands", file_no+1, file_no+1);
    }
    if ( got_truncate && len > min_len && !CONTRACTIVE(fdiffact,file_no) ) {
      error("ERROR: observation file %u needs an -fdiffact%u that truncates", file_no+1, file_no+1);
    }
    if (!fdiffact || fdiffact[file_no] == FRAMEMATCH_ERROR) {
      if (got_truncate && len != min_len) {
	error("ERROR: observation file %u needs an -fdiffact%u that truncates", file_no+1, file_no+1);
      } else if (len != max_len) {
	error("ERROR: observation file %u needs an -fdiffact%u that expands", file_no+1, file_no+1);
      }
    }
  }
  
  if(got_truncate) {
    // Adjust length of segments in the truncate case
    // We don't need to do that in the expand case
    if(min_len == 0)
      error("ERROR: MergeFile:openSegment:checkNumFrames: minimum segment length is zero");
    return min_len;
  } else {
    return max_len;  // if there is no expand, it means min_len == max_len
  }
}


// Returns the number of segments after adjusting for -sdiffactX
static
unsigned
checkNumSegments(unsigned nFiles, ObservationFile *file[],
		 unsigned const *sdiffact)
{
  unsigned min_len=file[0]->numLogicalSegments();
  unsigned max_len=min_len;
  if(max_len == 0) 
    warning("WARNING: MergeFile::checkNumSegments:  file 0 is empty");

  bool got_error    = false;
  bool got_truncate = false;
  bool got_expand   = false;

  for(unsigned file_no=0; file_no < nFiles; ++file_no) {
    unsigned len=file[file_no]->numLogicalSegments();
    if(len < min_len) min_len=len;
    if(len > max_len) max_len=len;

    if(	(sdiffact != NULL && sdiffact[file_no] == SEGMATCH_ERROR) ) {
      got_error = true;
    // the default is always to truncate (i.e. NULL case below)
    } else if(sdiffact == NULL || sdiffact[file_no] == SEGMATCH_TRUNCATE_FROM_END) {
      got_truncate =true;
    } else if( sdiffact != NULL &&
	     (sdiffact[file_no] == SEGMATCH_REPEAT_LAST ||
	      sdiffact[file_no] == SEGMATCH_WRAP_AROUND)
	     )
    {
      got_expand = true;
    }
  }
  if(max_len == 0)
    error("ERROR: MergeFile::checkNumSegments:  All files have zero length");

  // error checking
  if(min_len != max_len) {
    if(got_error)
      error("ERROR: MergeFile::checkNumSegments: Files have different # of segments (min=%d, max=%d)",min_len,max_len);

    if(got_truncate && got_expand)
      error("ERROR: MergeFile::checkNumSegments: Cannot specify both truncate and expand actions when using observation files of different lengths");
  }

  if(got_truncate) {
    // Adjust length of streams in the truncate case
    // We don't need to do that in the expand case
    if(min_len == 0)
      error("ERROR: MergeFile::checkNumSegments: minimum file length is zero");
    return min_len;
  } else {
    return max_len;  // if there is no expand, it means min_len=max_len
  }
}


MergeFile::MergeFile(unsigned nFiles, ObservationFile *file[], 
		     unsigned const *sdiffact, 
		     unsigned const *fdiffact,
		     int ftrcombo)
  : ObservationFile(NULL,NULL,NULL,NULL),
    nFiles(nFiles), ftrcombo(ftrcombo)
{
  this->sdiffact = sdiffact;
  this->fdiffact = fdiffact;
  this->file = new ObservationFile *[nFiles];
  for (unsigned i=0; i < nFiles; i+=1) {
    assert(file[i]);
    this->file[i] = file[i];
  }
  this->buffer         = NULL;
  this->buffSize       = 0;
  this->floatStart     = new unsigned[nFiles];
  this->intStart       = new unsigned[nFiles];
  this->segment        = -1;  // no openSegment() call yet

  // setup floatStart - append files for FTROP_NONE, otherwise 0 so they all overlap
  unsigned offset = 0;
  unsigned maxFloats = 0;
  for (unsigned i=0; i < nFiles; i+=1) {
    this->floatStart[i] = offset;
    unsigned nlc = file[i]->numLogicalContinuous();
    if (ftrcombo == FTROP_NONE) offset += nlc;
    if (maxFloats < nlc) maxFloats = nlc;
  }
  // discrete features always append
  if (ftrcombo != FTROP_NONE) {
    offset = maxFloats;
  }
  this->_numContinuousFeatures = offset;
  this->_numDiscreteFeatures = 0;
  for (unsigned i=0; i < nFiles; i+=1) {
    this->intStart[i] = offset;
    unsigned nld = file[i]->numLogicalDiscrete();
    offset += nld;
    this->_numDiscreteFeatures += nld;

#if !(ALLOW_VARIABLE_DIM_COMBINED_STREAMS)
    if(ftrcombo!=FTROP_NONE && i > 0) {
      if(file[i]->numLogicalContinuous() != file[i-1]->numLogicalContinuous()) {
	error("ERROR: MergeFile: When doing feature combination, the number of floats across files has to be the same.\n");
      }
    }
#endif
  }
  this->_numFeatures = _numContinuousFeatures + _numDiscreteFeatures;
  this->_numLogicalContinuousFeatures = _numContinuousFeatures;
  this->_numLogicalDiscreteFeatures   = _numDiscreteFeatures;
  this->_numLogicalFeatures = _numFeatures;

  this->bufStride = offset;
  this->_numSegments = checkNumSegments(nFiles, file, sdiffact);
}


bool
MergeFile::openSegment(unsigned seg) {
  if (seg >= _numSegments) {
    error("ERROR: MergeFile::openSegment: requested segment %u, but only up to %u are available", 
	  seg, _numSegments-1);
  }
  for (unsigned i=0; i < nFiles; i+=1) {
    if (!file[i]->openLogicalSegment(adjustForSdiffact(i,seg)))
      return false;
  }
  this->segment = seg;

  unsigned numInputFrames = checkNumFrames(nFiles, file, fdiffact);
  this->_numFrames = numInputFrames;
  // The modular debugging tests require the output to match the output 
  // from before the O(1) observation code, so infoMsg() calls in this
  // code will cause the tests to fail. So defining JEFFS_STRICT_DEBUG_OUTPUT_TEST
  // turns off the infoMsg() to allow the tests to pass.
#ifndef JEFFS_STRICT_DEBUG_OUTPUT_TEST
  infoMsg(IM::ObsFile,IM::Low,"%u input frames in segment %d\n", numInputFrames, seg);
#endif
  return true;
}


// repCount() returns the number of times a merged frame should be repeated in the
// merged getFrames() output to achieve any expansion requested by -fdiffactX.
//
// frameNum      pre-expansion frame # we need repetition count for
// deltaT        difference in frames between target segment length and actual 
//                 pre-expansion length
// firstFrame    first frame # desired (post-expansion)
// frameCount    # of desired frames (post-expansion) - i.e., trying to get merged 
//                 frames [firstFrame,firstFrame+frameCount)
// segLength     total # of frames in pre-expansion segment

static unsigned 
repCount(unsigned frameNum, unsigned const *fdiffact, unsigned fileNum, 
	 unsigned deltaT, unsigned firstFrame, unsigned frameCount, unsigned segLength)
{
  if (!fdiffact) return 1; // no expansion needed

  if (fdiffact[fileNum] == FRAMEMATCH_REPEAT_FIRST && frameNum == 0) {
    // repeat first frame to make the file match the target frameCount
    return frameCount > deltaT + 1 - firstFrame  ?  
      deltaT + 1 - firstFrame  :  frameCount;
  } else if (fdiffact[fileNum] == FRAMEMATCH_REPEAT_LAST && frameNum == segLength-1) {
    // repeat last frame to make the file match the target frameCount
    return firstFrame >= segLength  ?
      frameCount  :  frameCount - (segLength - firstFrame) + 1;
  } else if (fdiffact[fileNum] == FRAMEMATCH_EXPAND_SEGMENTALLY) {

    // most frames are repeated frameReps times, but the first remainder are
    // repeated 1 extra time
    unsigned frameReps = (segLength + deltaT) / segLength;
    unsigned remainder = (segLength + deltaT) % segLength;

    // map the post-expansion [firstFrame,firstFrame+frameCount) frame
    // range to the corresponding pre-expansion [preExpFirst,preExpLast]
    unsigned preExpFirst = (firstFrame < remainder * (frameReps+1)) ? 
      (firstFrame / (frameReps+1)) : 
      (remainder + (firstFrame - remainder * (frameReps + 1)) / frameReps);

    unsigned postExpLast = firstFrame + frameCount - 1;
    unsigned preExpLast = (postExpLast < remainder * (frameReps+1)) ? 
      (postExpLast / (frameReps+1)) : 
      (remainder + (postExpLast - remainder * (frameReps + 1)) / frameReps);

    assert(preExpFirst <= frameNum && frameNum <= preExpLast);

    if (preExpFirst == preExpLast) return frameCount;

    unsigned runLength = (frameNum < remainder) ?  frameReps + 1  :  frameReps;

    // The requested frame is repeated runLength times in the merged
    // segment, but we need to determin the overlap between that and
    // the requested frame range.
    if (frameNum == preExpFirst) {
      if (firstFrame < remainder * (frameReps + 1)) 
	firstFrame -=  remainder * (frameReps + 1);
      unsigned positionInRun = firstFrame % runLength;
      return runLength - positionInRun;
    }

    if (frameNum == preExpLast) {
      if (postExpLast < remainder * (frameReps + 1)) 
	postExpLast -=  remainder * (frameReps + 1);
      unsigned positionInRun = postExpLast % runLength;
      return 1 + positionInRun;
    }

    return runLength;
  }
  return 1;
}

  // Map requested merged frame range [first,first+count) to individual file 
  // fileNum's pre-fdiffactX frame range [adjFirst,adjFirst+adjCount). Also 
  // returns deltaT, the difference in frames between the merged segment 
  // length and the individual file's segment length.

  // precondition: already verified that fdiffact != NULL and 
  //               fdiffact[i] != FRAMEMATCH_ERROR

void 
MergeFile::adjustForFdiffact(unsigned first, unsigned count, unsigned i,
			     unsigned &adjFirst, unsigned &adjCount, unsigned &deltaT)
{
  if (fdiffact[i] == FRAMEMATCH_REPEAT_FIRST) {
    deltaT = _numFrames - file[i]->numLogicalFrames();
    assert(_numFrames >= file[i]->numLogicalFrames());
    adjFirst =  first >= deltaT  ?  first - deltaT  :  0;
    adjCount = count;
  } else if (fdiffact[i] == FRAMEMATCH_REPEAT_LAST) {
    deltaT = _numFrames - file[i]->numLogicalFrames();
    assert(_numFrames >= file[i]->numLogicalFrames());
    adjFirst = first < file[i]->numLogicalFrames() ?
      first  :  file[i]->numLogicalFrames() - 1;
    adjCount =  first + count <= file[i]->numLogicalFrames()  ?
      count  :  file[i]->numLogicalFrames() - adjFirst;
  } else if (fdiffact[i] == FRAMEMATCH_EXPAND_SEGMENTALLY) {
    deltaT = _numFrames - file[i]->numLogicalFrames();
    assert(_numFrames >= file[i]->numLogicalFrames());
    
    unsigned frameReps = _numFrames / file[i]->numLogicalFrames();
    unsigned remainder =  _numFrames % file[i]->numLogicalFrames();
    unsigned last = first + count - 1;
    unsigned adjLast;
    if (first < remainder * (frameReps+1)) {
      adjFirst = first / (frameReps + 1);
    } else {
      adjFirst = remainder + 
	(first - remainder * (frameReps + 1)) / frameReps;
    }
    if (last < remainder * (frameReps+1)) {
      adjLast = last / (frameReps + 1);
    } else {
      adjLast = remainder + 
	(last - remainder * (frameReps + 1)) / frameReps;
    }
    assert(adjLast >= adjFirst);
    adjCount = adjLast - adjFirst + 1;
    assert(adjCount >= 1);
    
  } else if (fdiffact[i] == FRAMEMATCH_TRUNCATE_FROM_END) {
    // don't have to do anything for this case - the error checking
    // prevents accessing any frames after numFrames()
    deltaT = file[i]->numLogicalFrames() - _numFrames;
    assert(file[i]->numLogicalFrames() >= _numFrames);
    adjFirst = first;
    adjCount = count;
  } else if (fdiffact[i] == FRAMEMATCH_TRUNCATE_FROM_START) {
    deltaT = file[i]->numLogicalFrames() - _numFrames;
    assert(file[i]->numLogicalFrames() >= _numFrames);
    adjFirst = first + deltaT;
    adjCount = count;
  }
}


Data32 const *
MergeFile::getFrames(unsigned first, unsigned count) {
  assert(first < _numFrames);
  assert(first + count <= _numFrames);
  unsigned needed = count * numFeatures();
  if (needed > buffSize) {
    buffer = (Data32 *) realloc(buffer, needed * sizeof(Data32));
    assert(buffer);
    buffSize = needed;
  }
  memset(buffer, 0, sizeof(Data32) * needed);
  
  for (unsigned int i=0; i < nFiles; i+=1) {
    
    // in case we need to adjust for fdiffact
    unsigned adjFirst, adjCount, deltaT; 

    if (!fdiffact || fdiffact[i] == FRAMEMATCH_ERROR) {
      deltaT   = 0;
      adjFirst = first;
      adjCount = count;
    } else {
      adjustForFdiffact(first, count, i, adjFirst, adjCount, deltaT);
    }

    // if -fdiffacti expands this segment, we may need to
    // fiddle with (first,count):
    //   rf  numFrames() - file[i]->numLogicalFrames() extra first frames
    //   rl  numFrames() - file[i]->numLogicalFrames() extra last  frames
    //   se  approx. numFrames() / file[i]->numLogicalFrames extra of each frame
    Data32 const *fileBuf = file[i]->getLogicalFrames(adjFirst,adjCount);
    assert(fileBuf);
    Data32 *dst = buffer + floatStart[i];
    Data32 const *src = fileBuf;
    unsigned srcStride = file[i]->numLogicalFeatures();

#if 0
    // This original version has a case statement in the inner loop, which
    // requires less code but is expensive at run-time. 

    for (unsigned f=0; f < adjCount; f += 1) {
      unsigned frameReps;
      if (fdiffact[i] == FRAMEMATCH_ERROR)
	frameReps = 1;
      else
	frameReps = repCount(adjFirst + f,fdiffact, i, deltaT, first, count, file[i]->numLogicalFrames());
      assert(0 < frameReps && frameReps <= count);
      for (unsigned j=0; j < frameReps; j+=1) {
	float *fdst = (float *) dst;
	float *fsrc = (float *) src;
	switch(ftrcombo) {
	case FTROP_NONE:
	  memcpy((void *)dst, (const void *)src, file[i]->numLogicalContinuous() * sizeof(Data32));
	  break;
	case FTROP_ADD:
	  for (unsigned k=0; k < file[i]->numLogicalContinuous(); k+=1) {
	    fdst[k] += fsrc[k];
	  }
	  break;
	case FTROP_SUB:
	  if (i == 0) {
	    memcpy((void *)dst, (const void *)src, file[i]->numLogicalContinuous() * sizeof(Data32));
	  } else {
	    for (unsigned k=0; k < file[i]->numLogicalContinuous(); k+=1) {
	      fdst[k] -= fsrc[k];
	    }
	  }
	  break;
	case FTROP_MUL:
	  if (i == 0) {
	    memcpy((void *)dst, (const void *)src, file[i]->numLogicalContinuous() * sizeof(Data32));
	  } else {
	    for (unsigned k=0; k < file[i]->numLogicalContinuous(); k+=1) {
	      fdst[k] *= fsrc[k];
	    }
	  }
	  break;
	case FTROP_DIV:
	  if (i == 0) {
	    memcpy((void *)dst, (const void *)src, file[i]->numLogicalContinuous() * sizeof(Data32));
	  } else {
	    for (unsigned k=0; k < file[i]->numLogicalContinuous(); k+=1) {
	      if (fsrc[k] != 0) {
		fdst[k] /= fsrc[k];
	      }
	    }
	  }
	  break;
	default:
	  error("ERROR: unknown -comb option value %d", ftrcombo);
	}
	dst += bufStride;
      }
      src += srcStride;
    }

#else
    
    // This version exchanges the case statement with the loop over the frames
    // (relative to the implementation above), which results in longer but more efficient code
    
    switch(ftrcombo) {
    case FTROP_NONE:
      
      for (unsigned f=0; f < adjCount; f += 1) {
	unsigned frameReps = 1;	
	if (fdiffact && fdiffact[i] != FRAMEMATCH_ERROR)
	  frameReps = repCount(adjFirst + f,fdiffact, i, deltaT, first, count, file[i]->numLogicalFrames());
#if 0
	assert(0 < frameReps && frameReps <= count);
#endif
	do {
	  memcpy((void *)dst, (const void *)src, file[i]->numLogicalContinuous() * sizeof(Data32));
	  dst += bufStride;
	} while (--frameReps);
	src += srcStride;
      }
      break;
      

    case FTROP_ADD:
      
      for (unsigned f=0; f < adjCount; f += 1) {
	unsigned frameReps = 1;	
	if (fdiffact && fdiffact[i] != FRAMEMATCH_ERROR)
	  frameReps = repCount(adjFirst + f,fdiffact, i, deltaT, first, count, file[i]->numLogicalFrames());
#if 0
	assert(0 < frameReps && frameReps <= count);
#endif
	float *fdst = (float *) dst;
	float *fsrc = (float *) src;
	do {
	  for (unsigned k=0; k < file[i]->numLogicalContinuous(); k+=1) {
	    fdst[k] += fsrc[k];
	  }
	  dst += bufStride;
	} while (--frameReps);
	src += srcStride;
      }
      break;

    case FTROP_SUB:

      for (unsigned f=0; f < adjCount; f += 1) {
	unsigned frameReps = 1;		
	if (fdiffact && fdiffact[i] != FRAMEMATCH_ERROR)
	  frameReps = repCount(adjFirst + f,fdiffact, i, deltaT, first, count, file[i]->numLogicalFrames());
#if 0
	assert(0 < frameReps && frameReps <= count);
#endif
	float *fdst = (float *) dst;
	float *fsrc = (float *) src;
	do {
	  if (i == 0) {
	    memcpy((void *)dst, (const void *)src, file[i]->numLogicalContinuous() * sizeof(Data32));
	  } else {
	    for (unsigned k=0; k < file[i]->numLogicalContinuous(); k+=1) {
	      fdst[k] -= fsrc[k];
	    }
	  }
	  dst += bufStride;
	} while (--frameReps);
	src += srcStride;
      }
      break;

    case FTROP_MUL:

      for (unsigned f=0; f < adjCount; f += 1) {
	unsigned frameReps = 1;		
	if (fdiffact && fdiffact[i] != FRAMEMATCH_ERROR)
	  frameReps = repCount(adjFirst + f,fdiffact, i, deltaT, first, count, file[i]->numLogicalFrames());
#if 0
	assert(0 < frameReps && frameReps <= count);
#endif
	float *fdst = (float *) dst;
	float *fsrc = (float *) src;
	do {
	  if (i == 0) {
	    memcpy((void *)dst, (const void *)src, file[i]->numLogicalContinuous() * sizeof(Data32));
	  } else {
	    for (unsigned k=0; k < file[i]->numLogicalContinuous(); k+=1) {
	      fdst[k] *= fsrc[k];
	    }
	  }
	  dst += bufStride;
	} while (--frameReps);
	src += srcStride;
      }
      break;

    case FTROP_DIV:

      for (unsigned f=0; f < adjCount; f += 1) {
	unsigned frameReps = 1;		
	if (fdiffact && fdiffact[i] != FRAMEMATCH_ERROR)
	  frameReps = repCount(adjFirst + f,fdiffact, i, deltaT, first, count, file[i]->numLogicalFrames());
#if 0
	assert(0 < frameReps && frameReps <= count);
#endif
	float *fdst = (float *) dst;
	float *fsrc = (float *) src;
	do {
	  if (i == 0) {
	    memcpy((void *)dst, (const void *)src, file[i]->numLogicalContinuous() * sizeof(Data32));
	  } else {
	    for (unsigned k=0; k < file[i]->numLogicalContinuous(); k+=1) {
	      if (fsrc[k] != 0) {
		fdst[k] /= fsrc[k];
	      }
	    }
	  }
	  dst += bufStride;
	} while (--frameReps);
	src += srcStride;
      }
      break;

    default:
      error("ERROR: unknown -comb option value %d", ftrcombo);
    }
#endif

    src = fileBuf + file[i]->numLogicalContinuous();
    dst = buffer + intStart[i];
    for (unsigned f=0; f < adjCount; f += 1) {
      unsigned frameReps = 1;
      if (fdiffact && fdiffact[i] != FRAMEMATCH_ERROR)
	frameReps = repCount(adjFirst + f,fdiffact, i, deltaT, first, count, file[i]->numLogicalFrames());
#if 0
      assert(0 < frameReps && frameReps <= count);
#endif
      do {
        memcpy((void *)dst, (const void *)src, file[i]->numLogicalDiscrete() * sizeof(Data32));
	dst += bufStride;
      } while (--frameReps);
      src += srcStride;
    }
  }
  return buffer;
}
