
/*
 * GMTK_HDF5File.cc
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
#include <stdio.h>
#include <string.h>
#include <string>
using namespace std;

#include "machine-dependent.h"
#include "error.h"
#include "general.h"
#include "file_utils.h"

#include "GMTK_HDF5File.h"

// FIXME - configure should check for 1.8, won't build with 1.6
#if HAVE_LIBHDF5_CPP

// split up 
// filename:groupname:x_start,y_start;x_stride,y_stride;x_count,y_count
// into filename, groupname, etc.

void
parseHDF5slab(char *slabSpecifier, 
	      char const *&fileName, char const *&groupName,
	      unsigned &xStart,  unsigned &yStart,
	      unsigned &xStride, unsigned &yStride,
	      unsigned &xCount,  unsigned &yCount)
{
  // FIXME - make parts optional, add error checking
  char *t = strchr(slabSpecifier, ':');
  assert(t);
  *t = 0;
  fileName = strdup(slabSpecifier);

  char *s = t+1;
  t = strchr(s, ':');
  assert(t);
  *t = 0;
  groupName = strdup(s);

  xStart = (unsigned) strtoul(t+1, &t, 10);
  assert(*t==',');
  yStart = (unsigned) strtoul(t+1, &t, 10);
  assert(*t==';');

  xStride = (unsigned) strtoul(t+1, &t, 10);
  assert(*t==',');
  yStride = (unsigned) strtoul(t+1, &t, 10);
  assert(*t==';');

  xCount = (unsigned) strtoul(t+1, &t, 10);
  assert(*t==',');
  yCount = (unsigned) strtoul(t+1, &t, 10);
  assert(*t==0);
}

#define ALLOC_NUM_SEGS(type,name) \
  name = new type[numSlabs]; assert(name);

HDF5File::HDF5File(const char *name, unsigned num,
		   bool cppIfAscii, char const* cppCommandOptions,
		   char const *contFeatureRangeStr_, 
		   char const *discFeatureRangeStr_, 
		   char const *preFrameRangeStr_, 
		   char const *segRangeStr_)
  : ObservationFile(contFeatureRangeStr_, 
		    discFeatureRangeStr_, 
		    preFrameRangeStr_,
		    segRangeStr_),
    cppIfAscii(cppIfAscii),
    cppCommandOptions(cppCommandOptions),
    buffer(NULL), bufSize(0),
    curFile(NULL)
{
  FILE *fofFile = openCPPableFile(name, cppIfAscii, cppCommandOptions);
  if (!fofFile) 
    error("HDF5File: couldn't open '%s' fore reading\n", name);
  char **hdf5Name;
  (void) readFof(fofFile, name, numSlabs, hdf5Name, 
		 cppIfAscii, cppCommandOptions);
  closeCPPableFile(fofFile, cppIfAscii);
  
  ALLOC_NUM_SEGS(char const*,fileName);
  ALLOC_NUM_SEGS(char const*,groupName);
  ALLOC_NUM_SEGS(unsigned, xStart);
  ALLOC_NUM_SEGS(unsigned, yStart);
  ALLOC_NUM_SEGS(unsigned, xStride);
  ALLOC_NUM_SEGS(unsigned, yStride);
  ALLOC_NUM_SEGS(unsigned, xCount);
  ALLOC_NUM_SEGS(unsigned, yCount);
  
  for (unsigned i=0; i < numSlabs; i+=1) {
    parseHDF5slab(hdf5Name[i], fileName[i], groupName[i],
		  xStart[i],  yStart[i],
		  xStride[i], yStride[i],
		  xCount[i],  yCount[i]);
    free(hdf5Name[i]);
  }
  delete [] hdf5Name;

#if 0
  for (unsigned i=0; i < numSlabs; i+=1) {
    printf("%03u: %s : %s : %u , %u ; %u , %u ; %u , %u\n",
	   i, fileName[i], groupName[i], xStart[i], yStart[i],
	   xStride[i], yStride[i], xCount[i], yCount[i]);
  }
#endif
  assert(numSlabs);

  // determine # floats & # ints by looking into segment 0.
  // require all subsequent segments to have the same shape.

  H5File      tmpFile(fileName[0], H5F_ACC_RDONLY);
  DataSet     dataset   = tmpFile.openDataSet(groupName[0]);
  H5T_class_t dataType  = dataset.getTypeClass();
  DataSpace   dataspace = dataset.getSpace();
  hsize_t start[2]  = {xStart[0],  yStart[0]};
  hsize_t stride[2] = {xStride[0], yStride[0]};
  hsize_t count[2]  = {xCount[0],  yCount[0]};
  hsize_t end[2];

  dataspace.selectHyperslab(H5S_SELECT_SET, count, start, stride);
  dataspace.getSelectBounds(start,end);
  switch (dataType) {
  case H5T_INTEGER:
    nContinuous = 0;
    nDiscrete = end[1] - start[1] + 1;
    break;
  case H5T_FLOAT:
    nDiscrete = 0;
    nContinuous = end[1] - start[1] + 1;
    break;
  default:
    error("HDF5File: dataset '%s:%s' is not of integer or float type\n",
	  fileName[0],groupName[0]);
  }

  if (segRangeStr) segRange = new Range(segRangeStr, 0, numSlabs);
  if (contFeatureRangeStr) contFeatureRange = new Range(contFeatureRangeStr, 0, nContinuous);
  if (discFeatureRangeStr) discFeatureRange = new Range(discFeatureRangeStr, 0, nDiscrete);
}


HDF5File::~HDF5File() {
  if (fileName) {
    for (unsigned i=0; i < numSlabs; i+=1) 
      free((void*)fileName[i]);
    delete [] fileName;
  }
  if (groupName) {
    for (unsigned i=0; i < numSlabs; i+=1)
      free((void *)groupName[i]);
    delete [] groupName;
  }
  if (xStart)  delete xStart;
  if (yStart)  delete yStart;
  if (xStride) delete xStride;
  if (yStride) delete yStride;
  if (xCount)  delete xCount;
  if (yCount)  delete yCount;
  if (curFile) {
    curFile->close();
    delete curFile;
  }
  if (buffer) free(buffer);
}  


bool 
HDF5File::openSegment(unsigned seg) {
  assert(0 <= seg && seg < numSlabs);

  // throw in curFile != NULL to ensure the segment has in fact
  // been openned, rather than an unfortunate unintialized curSegment
  // that happens to equal seg.
  if (curFile && seg == curSegment) return true;

  if (curFile) {
    curFile->close();
    delete curFile;
    curFile = NULL;
  }
  curSegment = seg;
  curFile = new H5File(fileName[seg], H5F_ACC_RDONLY);
  DataSet     dataset   = curFile->openDataSet(groupName[seg]);
  H5T_class_t dataType  = dataset.getTypeClass();
  DataSpace   dataspace = dataset.getSpace();
  hsize_t start[2]  = {xStart[seg],  yStart[seg]};
  hsize_t stride[2] = {xStride[seg], yStride[seg]};
  hsize_t count[2]  = {xCount[seg],  yCount[seg]};
  hsize_t end[2];

  dataspace.selectHyperslab(H5S_SELECT_SET, count, start, stride);
  dataspace.getSelectBounds(start,end);
  unsigned nFeatures = end[1] - start[1] + 1;

  switch (dataType) {
  case H5T_INTEGER:
    if (nContinuous != 0) 
      error("HDF5File: dataset '%s:%s' must be of integer type\n",
	    fileName[seg], groupName[seg]);
    if (nDiscrete != nFeatures) 
      // FIXME - need to add :x,y;dx,dy;nx,ny for specificity
      error("HDF5File: dataset '%s:%s' has %u features, expected %u\n",
	    fileName[seg], groupName[seg], nFeatures, nDiscrete);
    break;
  case H5T_FLOAT:
    if (nDiscrete != 0) 
      error("HDF5File: dataset '%s:%s' must be of float type\n",
	    fileName[seg], groupName[seg]);
    if (nContinuous != nFeatures) 
      // FIXME - need to add :x,y;dx,dy;nx,ny for specificity
      error("HDF5File: dataset '%s:%s' has %u features, expected %u\n",
	    fileName[seg], groupName[seg], nFeatures, nContinuous);
    break;
  default:
    error("HDF5File: dataset '%s:%s' is not of integer or float type\n",
	  fileName[seg],groupName[seg]);
  }
  nFrames = end[0] - start[0] + 1;

  if (preFrameRange) delete preFrameRange;
  if (preFrameRangeStr) {
    preFrameRange = new Range(preFrameRangeStr, 0, nFrames);
  }
  return true;
}


Data32 const *
HDF5File::getFrames(unsigned first, unsigned count) {
  assert(curFile); // segment open
  assert(first+count <= nFrames);
    
  DataSet     dataset   = curFile->openDataSet(groupName[curSegment]);
  H5T_class_t dataType  = dataset.getTypeClass();
  DataSpace   dataspace = dataset.getSpace();
  hsize_t start[2]  = {xStart[curSegment] + first * xStride[curSegment],
		       yStart[curSegment]};
  hsize_t stride[2] = {xStride[curSegment], yStride[curSegment]};
  hsize_t length[2]  = {count,  yCount[curSegment]};
  hsize_t end[2];

  dataspace.selectHyperslab(H5S_SELECT_SET, length, start, stride);
  dataspace.getSelectBounds(start,end);
  unsigned nFeatures = end[1] - start[1] + 1;

  unsigned needed = count * nFeatures;
  if (needed > bufSize) {
    buffer = (Data32 *) realloc(buffer, needed * sizeof(Data32));
    assert(buffer);
    bufSize = needed;
  }
  hsize_t memdim[1] = {needed};
  DataSpace memspace(1, memdim);

  switch (dataType) {
  case H5T_INTEGER:
    dataset.read(buffer, PredType::NATIVE_UINT32, memspace, dataspace);
    break;
  case H5T_FLOAT:
    dataset.read(buffer, PredType::NATIVE_FLOAT, memspace, dataspace);
    break;
  default:
    error("HDF5File: dataset '%s:%s' is not of integer or float type\n",
	  fileName[curSegment],groupName[curSegment]);
  }
  return buffer;
}

#endif
