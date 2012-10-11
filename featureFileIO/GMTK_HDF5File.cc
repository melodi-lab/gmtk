
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

#if HAVE_LIBHDF5_CPP
// Only use HDF5 library if Autoconf found it


void
skipSpace(char *&s) {
  while(s && isspace(*s)) s+=1;
}


// split up 
// filename:groupname:x_start,y_start;x_stride,y_stride;x_count,y_count
// into filename, groupname, etc.

void
parseHDF5slab(char *slabSpecifier, 
	      char const *&fileName, char const *&groupName,
	      unsigned &xStart,  unsigned &yStart,
	      unsigned &xStride, unsigned &yStride,
	      unsigned &xCount,  unsigned &yCount,
	      const char *fofName)
{

  // filename:groupname: WS* int WS* , WS* int WS* ; WS* int WS* , WS* int WS* ; WS* int WS* , WS* int WS* \0

  char *slabDup = strdup(slabSpecifier); // for error messages
  assert(slabDup);

  char *t = strchr(slabSpecifier, ':');
  if (!t) {
    error("Expected ':' after filename in slab specifier '%s' in %s", 
	  slabSpecifier, fofName);
  }
  *t = 0;
  fileName = strdup(slabSpecifier);

  char *s = t+1;
  t = strchr(s, ':');
  if (!t) {
    error("Expected ':' after groupname in slab specifier '%s' in %s", 
	  slabDup, fofName);
  }
  *t = 0;
  groupName = strdup(s);

  const char *nptr = t+1;
  xStart = (unsigned) strtoul(nptr, &t, 10);  // strtoul() will eat the WS* before the int
  if (t == nptr) {
    error("Expected integer for x start in slab specifier '%s' in %s", 
	  slabDup, fofName);
  }
  skipSpace(t);
  if (!t || *t != ',') {
    error("Expected ',' after x start in slab specifier '%s' in %s", 
	  slabDup, fofName);
  }

  nptr = t+1;
  yStart = (unsigned) strtoul(nptr, &t, 10);
  if (t == nptr) {
    error("Expected integer for y start in slab specifier '%s' in %s", 
	  slabDup, fofName);
  }
  skipSpace(t);
  if (!t || *t != ';') {
    error("Expected ';' after y start in slab specifier '%s' in %s", 
	  slabDup, fofName);
  }

  nptr = t+1;
  xStride = (unsigned) strtoul(nptr, &t, 10);
  if (t == nptr) {
    error("Expected integer for x stride in slab specifier '%s' in %s", 
	  slabDup, fofName);
  }
  skipSpace(t);
  if (!t || *t != ',') {
    error("Expected ',' after x stride in slab specifier '%s' in %s", 
	  slabDup, fofName);
  }

  nptr = t+1;
  yStride = (unsigned) strtoul(nptr, &t, 10);
  if (t == nptr) {
    error("Expected integer for y stride in slab specifier '%s' in %s", 
	  slabDup, fofName);
  }
  skipSpace(t);
  if (!t || *t != ';') {
    error("Expected ';' after y stride in slab specifier '%s' in %s", 
	  slabDup, fofName);
  }
  
  nptr = t+1;
  xCount = (unsigned) strtoul(nptr, &t, 10);
  if (t == nptr) {
    error("Expected integer for x count in slab specifier '%s' in %s", 
	  slabDup, fofName);
  }
  skipSpace(t);
  if (!t || *t != ',') {
    error("Expected ',' after x count in slab specifier '%s' in %s", 
	  slabDup, fofName);
  }

  nptr = t+1;
  yCount = (unsigned) strtoul(nptr, &t, 10);
  if (t == nptr) {
    error("Expected integer for y count in slab specifier '%s' in %s", 
	  slabDup, fofName);
  }
  skipSpace(t);
  if (!t || *t != 0) {
    error("Expected ':' after groupname in slab specifier '%s' in %s", 
	  slabDup, fofName);
  }
  free(slabDup);
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
    error("HDF5File: couldn't open '%s' fore reading", name);
  char **hdf5Name;
  (void) readFof(fofFile, name, numSlabs, hdf5Name, 
		 cppIfAscii, cppCommandOptions);
  closeCPPableFile(fofFile, cppIfAscii);

  // parse the segment descriptors into arrays  
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
		  xCount[i],  yCount[i], name);
    free(hdf5Name[i]);
  }
  delete [] hdf5Name;
  if (numSlabs == 0) {
    error("HDF5 file '%s' contained no segments", name);
  }

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
    _numContinuousFeatures = 0;
    _numDiscreteFeatures = end[1] - start[1] + 1;
    break;
  case H5T_FLOAT:
    _numDiscreteFeatures = 0;
    _numContinuousFeatures = end[1] - start[1] + 1;
    break;
  default:
    error("HDF5File: dataset '%s:%s' is not of integer or float type",
	  fileName[0],groupName[0]);
  }
  _numFeatures = _numContinuousFeatures + _numDiscreteFeatures;

  if (segRangeStr) segRange = new Range(segRangeStr, 0, numSlabs);

  if (contFeatureRangeStr) {
    contFeatureRange = new Range(contFeatureRangeStr, 0, _numContinuousFeatures);
    assert(contFeatureRange);
    _numLogicalContinuousFeatures = contFeatureRange->length();
  } else
    _numLogicalContinuousFeatures = _numContinuousFeatures;

  if (discFeatureRangeStr) {
    discFeatureRange = new Range(discFeatureRangeStr, 0, _numDiscreteFeatures);
    assert(discFeatureRange);
    _numLogicalDiscreteFeatures = discFeatureRange->length();
  } else 
    _numLogicalDiscreteFeatures = _numDiscreteFeatures;

  _numLogicalFeatures = _numLogicalContinuousFeatures + _numLogicalDiscreteFeatures;
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
    if (_numContinuousFeatures != 0) 
      error("HDF5File: dataset '%s:%s:%u,%u;%u,%u;%u,%u' must be of integer type",
	    fileName[seg], groupName[seg], 
	    xStart[seg], yStart[seg], xStride[seg],yStride[seg],xCount[seg],yCount[seg]);
    if (_numDiscreteFeatures != nFeatures) 
      error("HDF5File: dataset '%s:%s:%u,%u;%u,%u;%u,%u' has %u features, expected %u",
	    fileName[seg], groupName[seg],
	    xStart[seg], yStart[seg], xStride[seg],yStride[seg],xCount[seg],yCount[seg],
	    nFeatures, _numDiscreteFeatures);
    break;
  case H5T_FLOAT:
    if (_numDiscreteFeatures != 0) 
      error("HDF5File: dataset '%s:%s:%u,%u;%u,%u;%u,%u' must be of float type\n",
	    fileName[seg], groupName[seg],
	    xStart[seg], yStart[seg], xStride[seg],yStride[seg],xCount[seg],yCount[seg]);
    if (_numContinuousFeatures != nFeatures) 
      error("HDF5File: dataset '%s:%s:%u,%u;%u,%u;%u,%u' has %u features, expected %u",
	    fileName[seg], groupName[seg], 
	    xStart[seg], yStart[seg], xStride[seg],yStride[seg],xCount[seg],yCount[seg],
	    nFeatures, _numContinuousFeatures);
    break;
  default:
    error("HDF5File: dataset '%s:%s:%u,%u;%u,%u;%u,%u' is not of integer or float type",
	  fileName[seg],groupName[seg],
	  xStart[seg], yStart[seg], xStride[seg],yStride[seg],xCount[seg],yCount[seg]);
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
