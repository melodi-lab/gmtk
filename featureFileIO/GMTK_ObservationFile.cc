
/*
 * GMTK_ObservationFile.cc
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2012, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <assert.h>

#include "GMTK_ObservationFile.h"

unsigned 
ObservationFile::numLogicalSegments() {
  // Trying to simplify life for subclass authors by hiding the
  // range handling, so (redundantly) check if the segment range
  // needs to be instantiated here...
  if (!segRange && segRangeStr) {
    segRange = new Range(segRangeStr, 0, numSegments());
    assert(segRange);
  }
  if (segRange) {
    return segRange->length();
  } else {
    return numSegments();
  }
}  

bool 
ObservationFile::openLogicalSegment(unsigned seg) {
  // Trying to simplify life for subclass authors by hiding the
  // range handling, so (redundantly) check if the segment range
  // needs to be instantiated here...
  if (!segRange && segRangeStr) {
    segRange = new Range(segRangeStr, 0, numSegments());
    assert(segRange);
  }
  bool result;
  if (segRange) {
    result = openSegment(segRange->index((int)seg));
  } else {
    result = openSegment(seg);
  }
  if (!result) return false;

  // Can't determine the number of physical frames in a segment until
  // it's open, so the preFrameRange has to be handled here
  if (preFrameRangeStr) {
    if (preFrameRange) delete preFrameRange;
    preFrameRange = new Range(preFrameRangeStr, 0, numFrames());
    assert(preFrameRange);
  }
  return result;
}

unsigned 
ObservationFile::numLogicalFrames() {
  // openSegment() must be called first, so preFrameRange is already handled
  if (preFrameRange) {
    return preFrameRange->length();
  } else {
    return numFrames();
  }
}

Data32 const *
ObservationFile::getLogicalFrames(unsigned first, unsigned count) {
  if (!preFrameRange && !contFeatureRange && !discFeatureRange) {
    return getFrames(first, count);
  }

  unsigned needed = numLogicalFeatures() * count;
  assert(count > 0);
  if (needed > logicalObsBufSize) {
    logicalObservationBuffer = (Data32 *) 
      realloc(logicalObservationBuffer, needed * sizeof(Data32));
    assert(logicalObservationBuffer);
    logicalObsBufSize = needed;
  }
    Data32 *dest = logicalObservationBuffer;
  for (unsigned f = first; f < first+count; f+=1) {
    unsigned physFrameIdx = preFrameRange ? preFrameRange->index(f) : f;
    Data32 const *physicalFrame = getFrames(physFrameIdx, 1);
    assert(physicalFrame);
    for (unsigned i=0; i < numLogicalContinuous(); i+=1) {
      *(dest++) = physicalFrame[ contFeatureRange->index(i) ];
    }
    for (unsigned i=0; i < numLogicalDiscrete(); i+=1) {
      *(dest++) = physicalFrame[ numLogicalContinuous() +
				 discFeatureRange->index(i) ];
    }
  }
  return logicalObservationBuffer;
}

unsigned 
ObservationFile::numLogicalContinuous() {
  // Trying to simplify life for subclass authors by hiding the
  // range handling, so (redundantly) check if the feature range
  // needs to be instantiated here...
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

unsigned 
ObservationFile::numLogicalDiscrete() {
  // Trying to simplify life for subclass authors by hiding the
  // range handling, so (redundantly) check if the feature range
  // needs to be instantiated here...
  if (!discFeatureRange && discFeatureRangeStr) {
    discFeatureRange = new Range(discFeatureRangeStr, 0, numDiscrete());
    assert(contFeatureRange);
  }
  if (discFeatureRange) {
    return discFeatureRange->length();
  } else {
    return numDiscrete();
  }
}

