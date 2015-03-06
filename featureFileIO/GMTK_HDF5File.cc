
/*
 * GMTK_HDF5File.cc
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (C) 2012 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 * 
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
    error("Expected ':' after filename in slab specifier '%s' in %s\n", 
	  slabSpecifier, fofName);
  }
  *t = 0;
  fileName = strdup(slabSpecifier);

  char *s = t+1;
  t = strchr(s, ':');
  if (!t) {
    error("Expected ':' after groupname in slab specifier '%s' in %s\n", 
	  slabDup, fofName);
  }
  *t = 0;
  groupName = strdup(s);

  const char *nptr = t+1;
  xStart = (unsigned) strtoul(nptr, &t, 10);  // strtoul() will eat the WS* before the int
  if (t == nptr) {
    error("Expected integer for x start in slab specifier '%s' in %s\n", 
	  slabDup, fofName);
  }
  skipSpace(t);
  if (!t || *t != ',') {
    error("Expected ',' after x start in slab specifier '%s' in %s\n", 
	  slabDup, fofName);
  }

  nptr = t+1;
  yStart = (unsigned) strtoul(nptr, &t, 10);
  if (t == nptr) {
    error("Expected integer for y start in slab specifier '%s' in %s\n", 
	  slabDup, fofName);
  }
  skipSpace(t);
  if (!t || *t != ';') {
    error("Expected ';' after y start in slab specifier '%s' in %s\n", 
	  slabDup, fofName);
  }

  nptr = t+1;
  xStride = (unsigned) strtoul(nptr, &t, 10);
  if (t == nptr) {
    error("Expected integer for x stride in slab specifier '%s' in %s\n", 
	  slabDup, fofName);
  }
  skipSpace(t);
  if (!t || *t != ',') {
    error("Expected ',' after x stride in slab specifier '%s' in %s\n", 
	  slabDup, fofName);
  }

  nptr = t+1;
  yStride = (unsigned) strtoul(nptr, &t, 10);
  if (t == nptr) {
    error("Expected integer for y stride in slab specifier '%s' in %s\n", 
	  slabDup, fofName);
  }
  skipSpace(t);
  if (!t || *t != ';') {
    error("Expected ';' after y stride in slab specifier '%s' in %s\n", 
	  slabDup, fofName);
  }
  
  nptr = t+1;
  xCount = (unsigned) strtoul(nptr, &t, 10);
  if (t == nptr) {
    error("Expected integer for x count in slab specifier '%s' in %s\n", 
	  slabDup, fofName);
  }
  skipSpace(t);
  if (!t || *t != ',') {
    error("Expected ',' after x count in slab specifier '%s' in %s\n", 
	  slabDup, fofName);
  }

  nptr = t+1;
  yCount = (unsigned) strtoul(nptr, &t, 10);
  if (t == nptr) {
    error("Expected integer for y count in slab specifier '%s' in %s\n", 
	  slabDup, fofName);
  }
  skipSpace(t);
  if (!t || *t != 0) {
    error("Expected ':' after groupname in slab specifier '%s' in %s\n", 
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
		   char const *segRangeStr_,
		   unsigned leftPad, unsigned rightPad)
  : ObservationFile(name, num,
		    contFeatureRangeStr_, 
		    discFeatureRangeStr_, 
		    preFrameRangeStr_,
		    segRangeStr_,
		    leftPad, rightPad),
    cppIfAscii(cppIfAscii),
    cppCommandOptions(cppCommandOptions),
    buffer(NULL), bufSize(0),
    curFile(NULL), 
    listFile(NULL), fofName(NULL), outputFile(NULL), frameBuffer(NULL),
    contDataspace(NULL), discDataspace(NULL)
{
  FILE *fofFile = openCPPableFile(name, cppIfAscii, cppCommandOptions);
  if (!fofFile) 
    error("HDF5File: couldn't open '%s' for reading\n", name);
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
    error("HDF5 file '%s' contained no segments\n", name);
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
    error("HDF5File: dataset '%s:%s' is not of integer or float type in '%s'\n",
	  fileName[0],groupName[0], name);
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


// write ctor
HDF5File::HDF5File(char const *listFileName, char const *outputFileName, unsigned nfloats, unsigned nints)
  : fofName(listFileName), hdf5Name(outputFileName)
{
  if (fofName) {
    if ((listFile = fopen(fofName, "w")) == NULL) {
      error("ERROR: Couldn't open output list (%s) for writing.\n", fofName);
    }
  } else {
    error("ERROR: you must specify an HDF5 list file name\n");
  }
  _numContinuousFeatures = nfloats;
  _numDiscreteFeatures = nints;
  _numFeatures = nfloats + nints;

  frameBuffer = new Data32[_numFeatures];
  
  curSegment = 0;
  curFrame = 0;
  curFeature = 0;
  frameCount = 0;
  segStart = 0;

  try {
    // Let GMTK do the error logging
    Exception::dontPrint();

    // create the HDF5 file
    outputFile = new H5File(outputFileName, H5F_ACC_TRUNC);

    hsize_t dims[2];
    hsize_t maxdims[2];
    hsize_t chunk_dims[2];

    // Create continuous & discrete dataspaces & datasets

    if (nfloats > 0) {
      dims[0] = 1;
      dims[1] = nfloats;
      maxdims[0] = H5S_UNLIMITED;
      maxdims[1] = nfloats;
      chunk_dims[0] = 1;
      chunk_dims[1] = nfloats;
      contDataspace = new DataSpace(2, dims, maxdims);
      
      // Modify dataset creation property to enable chunking
      DSetCreatPropList prop;
      prop.setChunk(2, chunk_dims);

      const H5std_string datasetName("continuous");
      contDataset = outputFile->createDataSet(datasetName, PredType::NATIVE_FLOAT, *contDataspace, prop);
    }
    if (nints > 0) {
      dims[0] = 1;
      dims[1] = nints;
      maxdims[0] = H5S_UNLIMITED;
      maxdims[1] = nints;
      chunk_dims[0] = 1;
      chunk_dims[1] = nints;
      discDataspace = new DataSpace(2, dims, maxdims);

      // Modify dataset creation property to enable chunking
      DSetCreatPropList prop;
      prop.setChunk(2, chunk_dims);

      const H5std_string datasetName("discrete");
      discDataset = outputFile->createDataSet(datasetName, PredType::NATIVE_UINT32, *discDataspace, prop);
    }
  } catch(Exception err) {
err.printError(stderr);
fprintf(stderr, "%s\n", err.getCFuncName());
fflush(stderr);
    string errmsg("ERROR: creating HDF5 output file '");
    errmsg = (errmsg + outputFileName) + "' failed:\n";
    errmsg = errmsg + err.getCDetailMsg();
    error(errmsg.c_str() );
  }
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
  if (listFile) {
    if (fclose(listFile)) {
      error("ERROR: failed to close output list file '%s'\n", fofName);
    }
  }
  if (outputFile) {
    delete outputFile;
  }
  if (buffer) free(buffer);
  if (frameBuffer) delete[] frameBuffer;
  if (contDataspace) delete contDataspace;
  if (discDataspace) delete discDataspace;
}  


// Write frame to the file (call endOfSegment after last frame of a segment)
void 
HDF5File::writeFeature(Data32 x) {
  assert(frameBuffer);
  assert(curFeature < _numFeatures);
  frameBuffer[curFeature++] = x;
  if (curFeature == _numFeatures) {
    try {
      Exception::dontPrint();
      if (_numContinuousFeatures > 0) {
	hsize_t dims[2] = {1, _numContinuousFeatures};
	DataSpace memspace(2, dims, NULL); // shape of the vector (in memory) that we're going to write
	if (curFrame > 0) { // have to extend the data set for new frame
	  hsize_t size[2] = {curFrame+1, _numContinuousFeatures};
	  contDataset.extend(size);
	}
	hsize_t offset[2] = {curFrame, 0}; // offset for new frame
	DataSpace filespace(contDataset.getSpace());
	filespace.selectHyperslab(H5S_SELECT_SET, dims, offset);
	contDataset.write(frameBuffer, PredType::NATIVE_FLOAT, memspace, filespace);
      }
      if (_numDiscreteFeatures > 0) {
	hsize_t dims[2] = {1, _numDiscreteFeatures};
	DataSpace memspace(2, dims, NULL); // shape of the vector (in memory) that we're going to write

	if (curFrame > 0) { // have to extend the data set for new frame
	  hsize_t size[2] = {curFrame+1, _numDiscreteFeatures};
	  discDataset.extend(size);
	}
	hsize_t offset[2] = {curFrame, 0}; // offset for new frame
	DataSpace filespace(discDataset.getSpace());
	filespace.selectHyperslab(H5S_SELECT_SET, dims, offset);
	discDataset.write(frameBuffer+_numContinuousFeatures, PredType::NATIVE_UINT32, memspace, filespace);
      }
    } catch(Exception err) {
err.printError(stderr);
fprintf(stderr, "%s\n", err.getCFuncName());
fflush(stderr);
      string errmsg("ERROR: writing to HDF5 output file '");
      errmsg = (errmsg + hdf5Name) + "' failed:\n";
      errmsg = errmsg + err.getCDetailMsg();
      error(errmsg.c_str() );
    }
    curFeature = 0;
    curFrame += 1;
  }
}


void 
HDF5File::endOfSegment() {
// filename:groupname:x_start,y_start;x_stride,y_stride;x_count,y_count
  assert(listFile);
  assert(curFeature == 0);
  if (_numContinuousFeatures > 0) {
    fprintf(listFile,"%s:/continuous:%u,0;1,1;%u,%u\n", hdf5Name, segStart, curFrame-segStart, _numContinuousFeatures);
  }
  if (_numDiscreteFeatures > 0) {
    fprintf(listFile,"%s:/discrete:%u,0;1,1;%u,%u\n", hdf5Name, segStart, curFrame-segStart, _numDiscreteFeatures);
  }
  segStart = curFrame;
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
      error("HDF5File: dataset '%s:%s:%u,%u;%u,%u;%u,%u' must be of integer type in observation file '%s'\n",
	    fileName[seg], groupName[seg], 
	    xStart[seg], yStart[seg], xStride[seg],yStride[seg],xCount[seg],yCount[seg], observationFileName);
    if (_numDiscreteFeatures != nFeatures) 
      error("HDF5File: dataset '%s:%s:%u,%u;%u,%u;%u,%u' has %u features, expected %u in observation file '%s'\n",
	    fileName[seg], groupName[seg],
	    xStart[seg], yStart[seg], xStride[seg],yStride[seg],xCount[seg],yCount[seg],
	    nFeatures, _numDiscreteFeatures, observationFileName);
    break;
  case H5T_FLOAT:
    if (_numDiscreteFeatures != 0) 
      error("HDF5File: dataset '%s:%s:%u,%u;%u,%u;%u,%u' must be of float type in observation file '%s'\n",
	    fileName[seg], groupName[seg],
	    xStart[seg], yStart[seg], xStride[seg],yStride[seg],xCount[seg],yCount[seg], observationFileName);
    if (_numContinuousFeatures != nFeatures) 
      error("HDF5File: dataset '%s:%s:%u,%u;%u,%u;%u,%u' has %u features, expected %u in observation file '%s'\n",
	    fileName[seg], groupName[seg], 
	    xStart[seg], yStart[seg], xStride[seg],yStride[seg],xCount[seg],yCount[seg],
	    nFeatures, _numContinuousFeatures, observationFileName);
    break;
  default:
    error("HDF5File: dataset '%s:%s:%u,%u;%u,%u;%u,%u' is not of integer or float type in observation file '%s'\n",
	  fileName[seg],groupName[seg],
	  xStart[seg], yStart[seg], xStride[seg],yStride[seg],xCount[seg],yCount[seg], observationFileName);
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
    error("HDF5File: dataset '%s:%s' is not of integer or float type in observation file '%s'\n",
	  fileName[curSegment],groupName[curSegment], observationFileName);
  }
  return buffer;
}

#else

// ticket 90: get rid of a warning from the OS X ranlib
char make_osx_ranlib_shutup_about_no_symbols_HDF5;

#endif
