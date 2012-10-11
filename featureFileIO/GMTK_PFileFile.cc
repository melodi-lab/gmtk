/*
 * GMTK_PFileFile.cc
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

#include <stdio.h>
#include <string.h>

#include "error.h"
#include "general.h"
#include "debug.h"
#include "pfile.h"

#include "GMTK_PFileFile.h"

using namespace std;

PFileFile::PFileFile(const char *name, unsigned nfloats, unsigned nints,
		     unsigned num, bool bswap, 
		     char const *contFeatureRangeStr_, 
		     char const *discFeatureRangeStr_, 
		     char const *preFrameRangeStr_, 
		     char const *segRangeStr_)
  : ObservationFile(contFeatureRangeStr_, 
		    discFeatureRangeStr_, 
		    preFrameRangeStr_,
		    segRangeStr_)
{
  assert(name);
  fileName = strdup(name);
  buffer = NULL;
  bufferSize = 0;
  contBuf = NULL;
  discBuf = NULL;
  currentSegment = 0;

  if (name == NULL) 	
    error("PFileFile: File name is NULL for stream %i",num);	
  if ((dataFile = fopen(name,"rb")) == NULL)
       error("PFileFile: Can't open '%s' for input", name);
  pfile = new InFtrLabStream_PFile(0,name,dataFile,1,bswap);
  assert(pfile);
  if (pfile->num_ftrs() != nfloats) 
    error("PFileFile: File %s has %i floats, expected %i",
	     name,
	     pfile->num_ftrs(),
	     nfloats);
     
  if (pfile->num_labs() != nints)
    error("PFileFile: File %s has %i ints, expected %i",
	     name,
	     pfile->num_labs(),
	     nints);

  _numContinuousFeatures = pfile->num_ftrs();
  _numDiscreteFeatures = pfile->num_labs();
  _numFeatures = _numContinuousFeatures + _numDiscreteFeatures;

  if (contFeatureRangeStr) {
    contFeatureRange = new Range(contFeatureRangeStr, 0, _numContinuousFeatures);
    assert(contFeatureRange);
infoMsg(IM::ObsFile, IM::Low, "-prefr%u '%s'\n", num, contFeatureRange);
    _numLogicalContinuousFeatures = contFeatureRange->length();
  } else
    _numLogicalContinuousFeatures = nfloats;
  if (discFeatureRangeStr) {
    discFeatureRange = new Range(discFeatureRangeStr, 0, _numDiscreteFeatures);
infoMsg(IM::ObsFile, IM::Low, "-preir%u '%s'\n", num, discFeatureRange);
    assert(discFeatureRange);
    _numLogicalDiscreteFeatures = discFeatureRange->length();
  } else
    _numLogicalDiscreteFeatures = nints;
  _numLogicalFeatures = _numLogicalContinuousFeatures + _numLogicalDiscreteFeatures;
}


bool
PFileFile::openSegment(unsigned seg) {
  assert(pfile);
  assert(seg < pfile->num_segs());
  long segId = pfile->set_pos(seg, 0);
  assert(segId != SEGID_BAD);
  currentSegment = seg;
  _numFrames = pfile->num_frames(currentSegment);
  return true;
}


Data32 const *
PFileFile::getFrames(unsigned first, unsigned count) {
  assert(pfile);
  assert(first < _numFrames  &&  first + count <= _numFrames);
  unsigned needed = _numFeatures * count;
  if (!buffer || needed > bufferSize) {
    infoMsg(IM::ObsFile, IM::Low, "PFileFile buffer resize %u - > %u B for [%u,%u)\n",
	    bufferSize * sizeof(Data32), needed * sizeof(Data32), first, first+count);
    buffer = (Data32 *) realloc(buffer, needed * sizeof(Data32));
    assert(buffer);
    bufferSize = needed;

    contBuf = (float *) realloc(contBuf, count * _numContinuousFeatures * sizeof(Data32));
    discBuf = (UInt32*) realloc(discBuf, count * _numDiscreteFeatures   * sizeof(Data32));
    assert(_numContinuousFeatures == 0 || contBuf != NULL);
    assert(_numDiscreteFeatures == 0   || discBuf != NULL);
  }
  assert(buffer);
  float  *contSrc = contBuf;
  UInt32 *discSrc = discBuf;
  if (pfile->set_pos(currentSegment, first) == SEGID_BAD) {
    error("ERROR: PFileFile: unable to seek in PFile %s", fileName);
  }
  unsigned framesRead = pfile->read_ftrslabs(count, contBuf, discBuf);
  assert(framesRead == count);
  Data32 *floatDest = buffer;
  Data32 *intDest   = buffer + _numContinuousFeatures;
  // FIXME - this wants a strided memcpy
  for (unsigned i=0; i < count; i+=1) {
    memcpy(floatDest, contSrc, _numContinuousFeatures * sizeof(Data32));
    memcpy(intDest, discSrc, _numDiscreteFeatures * sizeof(Data32));
    floatDest += _numFeatures;
    intDest   += _numFeatures;
    contSrc   += _numContinuousFeatures;
    discSrc   += _numDiscreteFeatures;
  }
  return buffer;
}

