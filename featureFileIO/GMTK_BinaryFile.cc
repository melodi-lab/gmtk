
/*
 * GMTK_BinaryFile.cc
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
#include "vbyteswapping.h"

#include "file_utils.h"
#include "GMTK_BinaryFile.h"

using namespace std;

BinaryFile::BinaryFile(const char *name, unsigned nfloats, unsigned nints, 
		       unsigned num, bool swap, bool cppIfAscii, 
		       char const* cppCommandOptions,
		       char const *contFeatureRangeStr_, 
		       char const *discFeatureRangeStr_, 
		       char const *preFrameRangeStr_, 
		       char const *segRangeStr_)
  : ObservationFile(contFeatureRangeStr_, 
		    discFeatureRangeStr_, 
		    preFrameRangeStr_,
		    segRangeStr_),
    swap(swap),
    cppIfAscii(cppIfAscii),
    cppCommandOptions(cppCommandOptions)
{
  buffer = NULL;
  buffSize = 0;
  curDataFile = NULL;
  if (name == NULL) 	
    error("BinaryFile: File name is NULL for stream %i\n",num);	
  if (nfloats == 0 && nints == 0)
    error("BinaryFile: number of float and int features cannot both be zero");

  // local copy of file name
  fofName = new char[strlen(name)+1];
  strcpy(fofName,name);
  nFloats = nfloats;
  nInts = nints;

  if (contFeatureRangeStr) {
    contFeatureRange = new Range(contFeatureRangeStr, 0, nFloats);
    assert(contFeatureRange);
  }
  if (discFeatureRangeStr) {
    discFeatureRange = new Range(discFeatureRangeStr, 0, nInts);
    assert(discFeatureRange);
  }
  fofFile = openCPPableFile(fofName, cppIfAscii, cppCommandOptions);
  if (!fofFile)
    error("BinaryFile: couldn't open '%s' for reading\n", fofName);

  // for some reason this asserts while the following doesn't
  // assert(numFileNames == readFof(fofFile));
  
  unsigned rfof = readFof(fofFile, fofName, numFileNames, dataNames,
			  cppIfAscii, cppCommandOptions);
  assert(numFileNames == rfof);
  closeCPPableFile(fofFile, cppIfAscii);
  if (segRangeStr) {
    segRange = new Range(segRangeStr,0,numFileNames);
    assert(segRange);
  }
}


// Begin sourcing data from the requested segment.
// Must be called before any other operations are performed on a segment.
bool
BinaryFile::openSegment(unsigned seg) {
  assert(seg < numFileNames);
  char *fname = dataNames[seg];

  if (fname == NULL) {
    warning("BinaryFile::openSegment: Filename is NULL for segment %u\n",seg);
    return false;
  }

  if (curDataFile) {
    if (fclose(curDataFile)) {
      // FIXME - track file name for error messages
      warning("BinaryFile::openSegment: failed to close data file\n");
      return false;
    }
    curDataFile = NULL;
  }
  curDataFile = fopen(dataNames[seg], "rb");
  if (curDataFile == NULL) {
    warning("BinaryFile::openSegment: Can't open '%s' for input\n",fname);
    return false;
  }

  if (fseek(curDataFile,0L,SEEK_END) == -1) {
    warning("BinaryFile::openSegment:: Can't skip to end of file %s",
	      dataNames[seg]);
  }
  unsigned long fsize = (unsigned long)ftell(curDataFile);

  if ((fsize % numFeatures()) > 0)
      error("BinaryFile::openSegment: odd number of bytes in file %s\n",dataNames[seg]);

  nFrames = fsize / sizeof(Data32) / numFeatures();

  if (preFrameRange) 
    delete preFrameRange;
  if (preFrameRangeStr) {
    preFrameRange = new Range(preFrameRangeStr, 0, nFrames);
    assert(preFrameRange);
  }
  return true;
}


Data32 const *
BinaryFile::getFrames(unsigned first, unsigned count) {
  assert(curDataFile);
  assert(first < nFrames);
  assert(first + count <= nFrames);
  unsigned needed = count * numFeatures();
  if (needed > buffSize) {
    buffer = (Data32 *) realloc(buffer, needed * sizeof(Data32));
    assert(buffer);
    buffSize = needed;
  }
  if (fseek(curDataFile,first * sizeof(Data32) * numFeatures(),SEEK_SET) == -1) {
    warning("BinaryFile::getFrames: Can't seek to frame %u", first);
    return NULL;
  }
  float* float_buffer_ptr = (float *) buffer;
  Int32* int_buffer_ptr   = (Int32 *)(buffer) + numContinuous();
  unsigned n_read = 0;
  for(unsigned s=0; s < count; s+=1) {
    n_read += fread((float*)float_buffer_ptr,sizeof(float), numContinuous(), curDataFile);
    // swap if needed.
    if(swap) {
      float tmp_float[1];
      for (unsigned i=0; i<numContinuous(); i+=1) {
	swapb_vf32_vf32(1,(float_buffer_ptr+i),tmp_float);
	float_buffer_ptr[i]=tmp_float[0];
      }
    }
    n_read += fread((Int32*)int_buffer_ptr,  sizeof(Int32), numDiscrete(), curDataFile);
    // swap if needed
    if(swap) {
      for (unsigned i=0; i<numDiscrete(); i+=1) {
	int_buffer_ptr[i]=swapb_i32_i32(int_buffer_ptr[i]);
      }
    }
    float_buffer_ptr += numFeatures();
    int_buffer_ptr   += numFeatures();
  }
  if (n_read != needed) {
    warning("BinaryFile::getFrames: read %i items, expected %i",n_read,needed);
    return NULL;
  }
  return buffer;
}
