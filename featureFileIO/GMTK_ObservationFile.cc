
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
#include <string.h>

#include "GMTK_ObservationFile.h"

#include "GMTK_ASCIIFile.h"
#include "GMTK_FlatASCIIFile.h"
#include "GMTK_PFileFile.h"
#include "GMTK_HTKFile.h"
#include "GMTK_HDF5File.h"
#include "GMTK_BinaryFile.h"

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
  if (preFrameRange) {
    return preFrameRange->length();
  } else {
    return numFrames();
  }
}

Data32 const *
ObservationFile::getLogicalFrames(unsigned first, unsigned count) {
  // if no range selection, logical = physical
  if (!preFrameRange && !contFeatureRange && !discFeatureRange) {
    return getFrames(first, count);
  }

  unsigned needed = _numLogicalFeatures * count;
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
    for (unsigned i=0; i < _numLogicalContinuousFeatures; i+=1) {
      unsigned srcIdx = contFeatureRange ? contFeatureRange->index(i) : i;
      *(dest++) = physicalFrame[srcIdx];
    }
    unsigned discOffset = _numContinuousFeatures;
    for (unsigned i=0; i < _numLogicalDiscreteFeatures; i+=1) {
      unsigned srcIdx = discFeatureRange ? discFeatureRange->index(i) : i;
      *(dest++) = physicalFrame[discOffset+srcIdx];
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
    contFeatureRange = new Range(contFeatureRangeStr, 0, _numContinuousFeatures);
    assert(contFeatureRange);
    _numLogicalContinuousFeatures = contFeatureRange->length();
  }
  return _numLogicalContinuousFeatures;
}

unsigned 
ObservationFile::numLogicalDiscrete() {
  // Trying to simplify life for subclass authors by hiding the
  // range handling, so (redundantly) check if the feature range
  // needs to be instantiated here...
  if (!discFeatureRange && discFeatureRangeStr) {
    discFeatureRange = new Range(discFeatureRangeStr, 0, _numDiscreteFeatures);
    assert(discFeatureRange);
    _numLogicalDiscreteFeatures = discFeatureRange->length();
  }
  return _numLogicalDiscreteFeatures;
}


ObservationFile *
instantiateFile(unsigned ifmt, char *ofs, unsigned nfs, unsigned nis,
		unsigned number, bool iswp, bool Cpp_If_Ascii, 
		char *cppCommandOptions, char const *frs, char const *irs, 
		char const *prepr, char const *sr)
{
  ObservationFile *obsFile = NULL;
  switch (ifmt) {
  case RAWASC:
    obsFile = new ASCIIFile(ofs, nfs, nis, number,
			       Cpp_If_Ascii, cppCommandOptions,
			       frs, irs, prepr, sr);
    break;
  case PFILE:
    obsFile = new PFileFile(ofs, nfs, nis, number, iswp, frs, irs, prepr, sr);
    break;
  case HTK:
    obsFile = new HTKFile(ofs, nfs, nis, number, iswp, Cpp_If_Ascii, cppCommandOptions,
			     frs, irs, prepr, sr);
    break;
  case HDF5:
    obsFile = new HDF5File(ofs, number, Cpp_If_Ascii, cppCommandOptions,
			      frs, irs, prepr, sr);
    break;
  case FLATASC:
    obsFile = new FlatASCIIFile(ofs, nfs, nis, number, Cpp_If_Ascii, cppCommandOptions,
				   frs, irs, prepr, sr);
    break;
  case RAWBIN:
    obsFile = new BinaryFile(ofs, nfs, nis, number, iswp, Cpp_If_Ascii, cppCommandOptions,
				frs, irs, prepr, sr);
    break;
  default:
    error("ERROR: Unknown observation file format type: '%s'\n", ifmt);
  }
  return obsFile;
}


int
formatStrToNumber(char const *fmt) {
  if      (strcmp(fmt,"htk") == 0)
    return HTK;
  else if (strcmp(fmt,"binary") == 0)
    return RAWBIN;
  else if (strcmp(fmt,"ascii") == 0)
    return RAWASC;
  else if (strcmp(fmt,"flatascii") == 0)
    return FLATASC;
  else if (strcmp(fmt,"hdf5") == 0) 
    return HDF5;
  else if (strcmp(fmt,"pfile") == 0)
    return PFILE;
  return -1;
}

