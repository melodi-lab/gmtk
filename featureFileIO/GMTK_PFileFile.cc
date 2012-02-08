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
  buffer = NULL;
  bufferSize = 0;
  currentSegment = 0;

  if (name == NULL) 	
    error("PFileFile: File name is NULL for stream %i\n",num);	
  if ((dataFile = fopen(name,"rb")) == NULL)
       error("PFileFile: Can't open '%s' for input\n", name);
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
}


// Begin sourcing data from the requested segment.
// Must be called before any other operations are performed on a segment.
bool
PFileFile::openSegment(unsigned seg) {
  assert(pfile);
  assert(seg < pfile->num_segs());
  long segId = pfile->set_pos(seg, 0);
  assert(segId != SEGID_BAD);
  currentSegment = seg;
  return true;
}


Data32 const *
PFileFile::getFrames(unsigned first, unsigned count) {
  assert(pfile);
  assert(first < numFrames()  &&  first + count <= numFrames());
  unsigned needed = numFeatures() * count;
  if (!buffer || needed > bufferSize) {
    buffer = (Data32 *) realloc(buffer, needed * sizeof(Data32));
    bufferSize = needed;
  }
  assert(buffer);
  // FIXME - move the new -> ctor and realloc iff too small; free in dtor
  float  *contBuf = new float[numContinuous() * needed];
  float  *contSrc = contBuf;
  UInt32 *discBuf = new UInt32[numDiscrete() * needed];
  UInt32 *discSrc = discBuf;
  assert(contBuf && discBuf);
  if (pfile->set_pos(currentSegment, first) == SEGID_BAD) {
    // FIXME - remember file name for error reporting
    error("ERROR: PFileFile: unable to seek in PFile");
  }
  unsigned framesRead = pfile->read_ftrslabs(count, contBuf, discBuf);
  assert(framesRead == count);
  Data32 *floatDest = buffer;
  Data32 *intDest   = buffer + numContinuous();
  for (unsigned i=0; i < count; i+=1) {
    memcpy(floatDest, contSrc, numContinuous() * sizeof(Data32));
    memcpy(intDest, discSrc, numDiscrete() * sizeof(Data32));
    floatDest += numFeatures();
    intDest += numFeatures();
    contSrc += numContinuous();
    discSrc += numDiscrete();
  }
  delete discBuf;
  delete contBuf;
  return buffer;
}

