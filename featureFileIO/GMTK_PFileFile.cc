
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
		     char const *_contFeatureRangeStr, 
		     char const *_discFeatureRangeStr, 
		     char const *_preFrameRangeStr, 
		     char const *_segRangeStr)
  : ObservationFile(_contFeatureRangeStr, 
		    _discFeatureRangeStr, 
		    _preFrameRangeStr,
		    _segRangeStr)
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
fprintf(stderr,"PFile ctor '%s'\n", name);
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
fprintf(stderr,"PFile open seg %u\n", seg);
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
fprintf(stderr,"PFile: reading [%u,%u) in seg %u\n", first, first+count, currentSegment);
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

