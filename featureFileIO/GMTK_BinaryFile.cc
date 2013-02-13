
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
  _numContinuousFeatures = nfloats;
  _numDiscreteFeatures   = nints;
  _numFeatures           = nfloats + nints;

  if (contFeatureRangeStr) {
    contFeatureRange = new Range(contFeatureRangeStr, 0, _numContinuousFeatures);
    assert(contFeatureRange);
    _numLogicalContinuousFeatures = contFeatureRange->length();
  } else
    _numLogicalContinuousFeatures = nfloats;
  if (discFeatureRangeStr) {
    discFeatureRange = new Range(discFeatureRangeStr, 0, _numDiscreteFeatures);
    assert(discFeatureRange);
    _numLogicalDiscreteFeatures = discFeatureRange->length();
  } else
    _numLogicalDiscreteFeatures = nints;
  _numLogicalFeatures = _numLogicalContinuousFeatures + _numLogicalDiscreteFeatures;

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


// Set frame # to write within current segemnt

void 
BinaryFile::setFrame(unsigned frame) {
  assert(currFeature == 0);
  if (frame == currFrame) return;

  if (!writeFile) { // previous EoS (or ctor) should have closed it
    assert(listFile);
    char* current_output_fname = new char[strlen(outputFileName)+strlen(outputNameSeparatorStr)+50];
    sprintf(current_output_fname,"%s%s%d",outputFileName,outputNameSeparatorStr,currSegment);
    if ((writeFile = fopen(current_output_fname, "wb")) == NULL) {
      error("Couldn't open output file '%s' for writing.",current_output_fname);
    }
    fprintf(listFile,"%s\n",current_output_fname);
    delete [] current_output_fname;
  }
  
  if (gmtk_fseek(writeFile, (gmtk_off_t)(frame * _numFeatures * sizeof(Data32)), SEEK_SET) == -1) {
    error("ERROR: failed to seek to frame %u in file '%s' segment %u\n", frame, fofName, curSegment);
  }
  currFrame = frame;
}


void 
BinaryFile::writeFrame(Data32 const *frame) {
  void (*copy_swap_func_ptr)(size_t, const intv_int32_t*, intv_int32_t*)=NULL;
  if(oswap) {
    copy_swap_func_ptr=&swapb_vi32_vi32;
  } else {
    copy_swap_func_ptr=&copy_vi32_vi32;
  }
  assert(currFeature == 0);
  if (!writeFile) { // previous EoS (or ctor) should have closed it
    assert(listFile);
    char* current_output_fname = new char[strlen(outputFileName)+strlen(outputNameSeparatorStr)+50];
    sprintf(current_output_fname,"%s%s%d",outputFileName,outputNameSeparatorStr,currSegment);
    if ((writeFile = fopen(current_output_fname, "wb")) == NULL) {
      error("Couldn't open output file '%s' for writing.",current_output_fname);
    }
    fprintf(listFile,"%s\n",current_output_fname);
    delete [] current_output_fname;
  }
  float  *cont_buf_p = new float[_numContinuousFeatures];
  copy_swap_func_ptr(_numContinuousFeatures, (const intv_int32_t *) frame, (intv_int32_t *)cont_buf_p);
  UInt32 *disc_buf_p = new UInt32[_numDiscreteFeatures];
  copy_swap_func_ptr(_numDiscreteFeatures, (const intv_int32_t *) frame + _numContinuousFeatures, (intv_int32_t *)disc_buf_p);

  /// Print continuous part of frame /////////////////////////////
  if (_numContinuousFeatures) {
    if (fwrite(cont_buf_p, sizeof(Data32), _numContinuousFeatures, writeFile) != _numContinuousFeatures) {
      error("Error writing to output file\n");
    }
  }
  /// Print discrete part of the frame ///////////////////////////
  if (_numDiscreteFeatures) {
    if (fwrite(disc_buf_p, sizeof(Data32), _numDiscreteFeatures, writeFile) != _numDiscreteFeatures) {
      error("Error writing to output file\n");
    }
  }
  currFrame += 1;
  currFeature = 0;
  delete []cont_buf_p;
  delete []disc_buf_p;
}


void 
BinaryFile::writeFeature(Data32 x) {
  void (*copy_swap_func_ptr)(size_t, const intv_int32_t*, intv_int32_t*)=NULL;
  if(oswap) {
    copy_swap_func_ptr=&swapb_vi32_vi32;
  } else {
    copy_swap_func_ptr=&copy_vi32_vi32;
  }
  if (!writeFile) { // previous EoS (or ctor) should have closed it
    assert(listFile);
    assert(currFrame == 0 && currFeature == 0); 
    char* current_output_fname = new char[strlen(outputFileName)+strlen(outputNameSeparatorStr)+50];
    sprintf(current_output_fname,"%s%s%d",outputFileName,outputNameSeparatorStr,currSegment);
    if ((writeFile = fopen(current_output_fname, "wb")) == NULL) {
      error("Couldn't open output file '%s' for writing.",current_output_fname);
    }
    fprintf(listFile,"%s\n",current_output_fname);
    delete [] current_output_fname;
  }
  copy_swap_func_ptr(1, &x, &x);
  if (fwrite(&x, sizeof(Data32), 1, writeFile) != 1) {
    error("Error writing to output file\n");
  }
  currFeature += 1;
  if (currFeature == _numFeatures) {
    currFrame += 1;
    currFeature = 0;
  }
}


void 
BinaryFile::endOfSegment() {
  assert(currFeature == 0);
#if 0
  if (currFrame == 0) {
    assert(listFile);
    assert(!writeFile); // previous EoS (or ctor) should have closed it
    char* current_output_fname = new char[strlen(outputFileName)+strlen(outputNameSeparatorStr)+50];
    sprintf(current_output_fname,"%s%s%d",outputFileName,outputNameSeparatorStr,currSegment);
    if ((writeFile = fopen(current_output_fname, "wb")) == NULL) {
      error("Couldn't open output file '%s' for writing.",current_output_fname);
    }
    fprintf(listFile,"%s\n",current_output_fname);
    delete [] current_output_fname;
  } else {
    assert(writeFile); // if not, what have we been writing to?
  }
#else
  assert(writeFile);
#endif
  if (fclose(writeFile)) {
    char* current_output_fname = new char[strlen(outputFileName)+strlen(outputNameSeparatorStr)+50];
    sprintf(current_output_fname,"%s%s%d",outputFileName,outputNameSeparatorStr,currSegment);
    error("ERROR closing output file '%s'\n", current_output_fname);
  }
  writeFile = NULL;
  currSegment += 1;
  currFrame = 0;
  currFeature = 0;
}


// Begin sourcing data from the requested segment.
// Must be called before any other operations are performed on a segment.
bool
BinaryFile::openSegment(unsigned seg) {
  assert(seg < numFileNames);
  unsigned prevSegment = curSegment;
  curSegment = seg;
  char *fname = dataNames[seg];

  if (fname == NULL) {
    warning("BinaryFile::openSegment: Filename is NULL for segment %u\n",seg);
    return false;
  }

  if (curDataFile) {
    if (fclose(curDataFile)) {
      warning("BinaryFile::openSegment: failed to close data file %s\n", dataNames[prevSegment]);
      return false;
    }
    curDataFile = NULL;
  }
  curDataFile = fopen(dataNames[seg], "rb");
  if (curDataFile == NULL) {
    warning("BinaryFile::openSegment: Can't open '%s' for input\n",fname);
    return false;
  }

  // get the file length to determine the # of frames in the segment
  if (gmtk_fseek(curDataFile,(gmtk_off_t)0,SEEK_END) == -1) {
    warning("BinaryFile::openSegment:: Can't skip to end of file %s",
	      dataNames[seg]);
  }
  gmtk_off_t fsize = gmtk_ftell(curDataFile);

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
  if (gmtk_fseek(curDataFile,(gmtk_off_t)(first * sizeof(Data32) * numFeatures()),SEEK_SET) == -1) {
    warning("BinaryFile::getFrames: Can't seek to frame %u in %s\n", first, dataNames[curSegment]);
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
    warning("BinaryFile::getFrames: read %i items, expected %i in %s\n",n_read,needed, dataNames[curSegment]);
    return NULL;
  }
  return buffer;
}
