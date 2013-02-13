
/*
 * GMTK_HTKFile.cc
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2011, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#if HAVE_INTTYPES_H
   // The ISO C99 standard specifies that the macros in inttypes.h must
   //  only be defined if explicitly requested. 
#  define __STDC_FORMAT_MACROS 1
#  include <inttypes.h>
#endif
#if HAVE_STDINT_H
#  include <stdint.h>
#endif

#include <stdio.h>
#include <string.h>

#include "error.h"
#include "general.h"
#include "vbyteswapping.h"

#include "file_utils.h"
#include "GMTK_HTKFile.h"

using namespace std;




void 
parseSentenceSpec(const string& sentLoc, int* startFrame, int* endFrame, 
		  string& fnameStr, char *fofName) 
{
  size_t fNameLen;
  *startFrame=0, *endFrame=-1; //these are the right values if the frame range is not specified
  if(sentLoc[sentLoc.length()-1]==']'){
    //have a subrange spec
    fNameLen=sentLoc.find_last_of('[');
    if (fNameLen==string::npos){
      error("ERROR: parseSentenceSpec: '%s' is an invalid sentence location in observation file '%s'.  "
	    "Must be of the form 'filename[startFrame:endFrame]'\n",sentLoc.c_str(), fofName);
    }

    string range= sentLoc.substr(fNameLen+1,sentLoc.length()-2-fNameLen);
    if (sscanf(range.c_str(),"%d:%d",startFrame,endFrame) != 2)
      error("ERROR: parseSentenceSpec: '%s' is an invalid sentence location in observation file '%s'.  "
	    "Must be of the form 'filename[startFrame:endFrame]'\n",sentLoc.c_str(), fofName);		  
    
    if(*endFrame<*startFrame)
      error("ERROR: parseSentenceSpec: %s has the last frame smaller than first frame.\n",sentLoc.c_str());	  
  }
  else{
    fNameLen = sentLoc.length();
  }
  fnameStr = sentLoc.substr(0,fNameLen);
  
}


/** Actually open the file designated by fname, seek to the beginning of the data,(past the header).
 * The HTK file is presumed to be in the following format (taken directly
 * from the HTK book:
 *
 *
 * 5.10.1 HTK Format Parameter Files
 *
 * HTK format files consist of a contiguous sequence of samples preceded
 * by a header. Each sample is a vector of either 2-byte integers or
 * 4-byte floats. 2-byte integers are used for compressed forms as
 * described below and for vector quantised data as described later in
 * section 5.14. HTK format data files can also be used to store speech
 * waveforms as described in section 5.11.  The HTK file format header is
 * 12 bytes long and contains the following data:
 *
 * nSamples - number of samples in file (4-byte integer)
 * sampPeriod -  sample period in 100ns units (4-byte integer)
 * sampSize -  number of bytes per sample (2-byte integer)
 * parmKind -  a code indicating the sample kind (2-byte integer)
 *
 * The parameter kind consists of a 6 bit code representing the basic
 * parameter kind plus additional bits for each of the possible
 * qualifiers. The basic parameter kind codes are
 *
 * 0 WAVEFORM sampled waveform
 * 1 LPC linear prediction filter coefficients
 * 2 LPREFC linear prediction reflection coefficients
 * 3 LPCEPSTRA LPC cepstral coefficients
 * 4 LPDELCEP LPC cepstra plus delta coefficients
 * 5 IREFC LPC reflection coef in 16 bit integer format
 * 6 MFCC mel-frequency cepstral coefficients
 * 7 FBANK log mel-filter bank channel outputs
 * 8 MELSPEC linear mel-filter bank channel outputs
 * 9 USER user defined sample kind
 * 10 DISCRETE vector quantised data
 *
 *  @pre f->curHTKFileInfo and f->curDataFile must be NULL
 *  @post f->curHTKFileInfo and f->curDataFile and opened and filled in appropriately
 *  @param fname the name of file to open
 *  @param f associated stream
 *  @return the total number of samples in the opened file
 */
unsigned
openHTKFile2(const string& fname, StreamInfo *f) {

  // structure of HTK header
  Int32 n_samples;
  Int32 samp_period;
  short samp_size;
  short parm_kind;

  Int32 tmp1,tmp2;

  short stmp1,stmp2;

  bool bswap = f->swap();
  int nints = f->nInts;
  int nfloats = f->nFloats;

  //printf ("call openHTKFile2(%s)\n",fname.c_str());

  if ((f->curDataFile = fopen(fname.c_str(),"rb")) == NULL) {
    error("ERROR: openHTKFile: Can't open observation file '%s' for input\n", fname.c_str());
  }

  if (fread(&tmp1,sizeof(Int32),1,f->curDataFile) != 1) {
    error("ERROR: openHTKFile: Can't read number of samples in observation file '%s'\n", fname.c_str());
  }


  if (fread((Int32 *)&tmp2,sizeof(Int32),1,f->curDataFile) != 1) {
    error("ERROR: openHTKFile: Can't read sample period from observation file '%s'\n", fname.c_str());
  }

  if (fread((short *)&stmp1,sizeof(short),1,f->curDataFile) != 1) {
    error("ERROR: openHTKFile: Can't read sample size from observation file '%s'\n", fname.c_str());
    return 0;
  }

  if (fread(&stmp2,sizeof(short),1,f->curDataFile) != 1) {
    error("ERROR: openHTKFile: Can't read parm kind from observation file '%s'\n", fname.c_str());
  }

  if (bswap) {
    n_samples = swapb_i32_i32(tmp1);
    samp_period = swapb_i32_i32(tmp2);
    samp_size = swapb_short_short(stmp1);
    parm_kind = swapb_short_short(stmp2);
  }
  else {
    n_samples = tmp1;
    samp_period = tmp2;
    samp_size = stmp1;
    parm_kind = stmp2;
  }

  if (n_samples <= 0) {
    error("ERROR: openHTKFile: number of samples %d should be > 0 in observation file '%s'\n",
	  n_samples, fname.c_str());
  }

  if (samp_period <= 0 || samp_period > 1000000) {
    warning("WARNING: openHTKFile: sample period is %d in observation file '%s' - "
	    "must be between 0 and 1000000\n", samp_period, fname.c_str());
  }

  if (samp_size <= 0 || samp_size > 5000) {
    warning("WARNING: openHTKFile: sample size %d in observation file '%s' - "
	    "must be between 0 and 5000\n", samp_size, fname.c_str());
  }

  short pk = parm_kind & BASEMASK;
  bool isCompressed = parm_kind & IS_COMPRESSED;

  if (pk < WAVEFORM || pk > ANON) {
    warning("WARNING: openHTKFile: Undefined parameter kind %d for HTK feature file '%s'. "
	    "Will assume float features.\n",pk, fname.c_str());
  }

  // For now we don't support the WAVEFORM and IREFC parameter kind.  It uses
  // shorts instead of floats and that requires special treatment.
  if (pk == WAVEFORM) {
    warning("WARNING: openHTKFile: HTK WAVEFORM parameter kind %d in not supported in observation file '%s'\n",
	    pk, fname.c_str());
  }
  else if (pk == IREFC) {
    warning("WARNING: openHTKFile: HTK IREFC parameter kind %d not supported in observation file '%s'\n",
	    pk, fname.c_str());
  }

  int n_fea;

  // parameter kind DISCRETE = all discrete features

  if (nfloats == 0) {
    // we divide by sizeof(short) which should be presumably size of a 2-byte int
    // (but this is MACHINE DEPENDENT).
    n_fea = samp_size / sizeof(short) ;
    if (n_fea != nints) {
      error("ERROR: openHTKFile:  Number of features (%i) does not match number of ints "
	    "specified (%i) in observation file '%s'\n", n_fea,nints, fname.c_str());
    }
    if(parm_kind != DISCRETE) {
      warning("WARNING: openHTKFile: Number of floats specified is 0 but the HTK parameter "
	      "kind is not DISCRETE for observation file '%s'.\n", fname.c_str());
    }
  }
  // otherwise all continuous features
  else {
    if(isCompressed)
      n_fea = samp_size / sizeof(short);
    else
      n_fea = samp_size / sizeof(float);
		
    if (n_fea != nfloats) {
      error("ERROR: openHTKFile:  Number of features (%i) does not match number of floats "
	    "specified (%i) in observation file '%s'\n", n_fea,nfloats, fname.c_str());
    }
    if(parm_kind == DISCRETE) {
      warning("WARNING: openHTKFile:  Number of floats specified (%i) is not 0 but the HTK "
	      "parameter kind is DISCRETE for observation file '%s'.\n",nfloats, fname.c_str());
    }
  }
	  
  float* scale =NULL;	  
  float* offset =NULL;
	  
  if (isCompressed){
    //the scale and offset floats take 8 bytes in total, 
    //which is the same as 4 compressed samples
    scale = new float[n_fea]; //A in htk book	  
    offset = new float[n_fea]; //B in htk book	  
    float* tmp = new float[n_fea];
    if (fread(tmp,sizeof(float),n_fea,f->curDataFile) != (unsigned short)n_fea) {
      error("ERROR: openHTKFile: Can't read scales for decompressing '%s'.\n", fname.c_str());
    }
    if(bswap)
      swapb_vf32_vf32(n_fea,tmp,scale);
    else
      copy_vf32_vf32(n_fea,tmp,scale);

    if (fread(tmp,sizeof(float),n_fea,f->curDataFile) != (unsigned short)n_fea) {
      error("ERROR: openHTKFile: Can't read offsets for decompressing '%s'.\n", fname.c_str());
    }
    if(bswap)
      swapb_vf32_vf32(n_fea,tmp,offset);
    else
      copy_vf32_vf32(n_fea,tmp,offset);

    n_samples -= 4;

  }

  f->curHTKFileInfo= new HTKFileInfo(samp_size,n_samples,gmtk_ftell(f->curDataFile),isCompressed,scale,offset);
  f->curDataFilename=fname;
	  
  return n_samples;
	
}

/**  open an HTK file for segment 'sentno' if necessary and seek to the beginning of data.
 *
 * HTK sentence locators (lines in the FoF file) can be filesnames 
 * or they can be of format filename[sentStart:sentEnd], where sentStart and sentEnd
 * are 0-based, indexes into the file.  So now you can specify an utterance as a 
 * consecutive subset of frames of the htk file, as described in the htk book, section
 * "4.2 Script Files", starting probably with HTK version 3.3
 *
 * 
 *  If another HTK file is open for a different segment is already open, it is first closed.
 *  If the sentno is contained within the previously opened HTK file, only an fseek is done
 *  instead of reopening the file. 
 *
 */
size_t 
openHTKFile(StreamInfo *f, size_t sentno) {

  DBGFPRINTF((stderr,"In openHTKFile, sentno %d\n",sentno));
  unsigned long htkfile_size = f->getFullFofSize();
  if(sentno < 0 || sentno >= htkfile_size) {
    error("ERROR: openHTKFile: Requested segment no %li of observation file '%s' "
	  "but the max num of segments in list of HTK files is %li\n",
	  sentno, f->fofName, htkfile_size);
  }

  //  assert(sentno >= 0 && sentno < _numSegments);
  if (f->dataNames[sentno] == NULL) {
    error("ERROR: openHTKFile: Filename is NULL for segment %li in observation file '%s'\n",
	  f->dataNames[sentno], f->fofName);
  }

  
  int startFrame, endFrame;
  string fnameStr;
  parseSentenceSpec(f->dataNames[sentno], &startFrame, &endFrame, fnameStr, f->fofName);

  if(f->curDataFilename != fnameStr && f->curDataFile){
    DBGFPRINTF((stderr,"In openHTKFile, f->curDataFilename  %s fnameStr %s f->curDataFile  %d\n", 
		f->curDataFilename.c_str(), fnameStr.c_str(), f->curDataFile));
    //the wrong file is open
    fclose(f->curDataFile);
    f->curDataFile = NULL;
    delete f->curHTKFileInfo;
    f->curHTKFileInfo= NULL;
  }
  
  if(!f->curDataFile)   
    openHTKFile2(fnameStr,f);
  
  const HTKFileInfo* htkInfo = f->curHTKFileInfo;
   
  if(endFrame<0) //endFrame, and therefore the range has not been specified: use the whole file
    endFrame= htkInfo->n_samples-1;
  
  if(endFrame>=htkInfo->n_samples)
    error("ERROR: openHTKFile: '%s' has the last frame at %d, beyond %d, which is the number of frames in file.\n",
	  f->dataNames[sentno], endFrame, htkInfo->n_samples);	  
  
  f->curNumFrames=endFrame-startFrame+1;
  
  
  DBGFPRINTF((stderr,"In openHTKFile, curNumFrames %d\n",f->curNumFrames));
  
  //now we seek to the start frame 
  if (gmtk_fseek(f->curDataFile, (gmtk_off_t)(htkInfo->startOfData+startFrame*htkInfo->samp_size), SEEK_SET)) {
    error("ERROR: openHTKFile: fseek() failed for '%s'\n", f->dataNames[sentno]);
  }
  return f->curNumFrames;
}


bool 
HTKFile::openSegment(unsigned seg) {
  assert(info);
  unsigned numPhysicalFrames = openHTKFile(info, seg);
  if (preFrameRange) {
    delete preFrameRange;
  }
  if (preFrameRangeStr) {
    preFrameRange = new Range(preFrameRangeStr, 0, numPhysicalFrames);
    assert(preFrameRange);
    nLogicalFrames = preFrameRange->length();
  } else {
    nLogicalFrames = numPhysicalFrames;
  }
  return true;
}


Data32 const *
HTKFile::getFrames(unsigned first, unsigned count) {
  assert(info && info->curHTKFileInfo && info->curDataFile);
  assert(count > 0);
  unsigned needed = numFeatures() * count;
  if (needed > bufferSize) {
    buffer = (Data32 *) realloc(buffer, needed * sizeof(Data32));
    assert(buffer);
    bufferSize = needed;
  }
  const HTKFileInfo *htkInfo = info->curHTKFileInfo;
  if (gmtk_fseek(info->curDataFile, (gmtk_off_t)(htkInfo->startOfData + first * htkInfo->samp_size), SEEK_SET)) {
    error("HTKFile: fseek() failed for '%s'\n", info->fofName);
  }

  // HTK files are either all discrete or all continuous
  unsigned featuresPerFrame;
  if (numContinuous()) {
    assert(numDiscrete() == 0);
    featuresPerFrame = numContinuous();
  } else {
    assert(numContinuous());
    featuresPerFrame = numDiscrete();
  }
  unsigned totalFeatures = featuresPerFrame * count;
  size_t nread;

  if (numDiscrete() > 0 || htkInfo->isCompressed) {
    Int16 *tmpBuf = new Int16[totalFeatures];
    assert(tmpBuf);
    nread = fread(tmpBuf, sizeof(Int16), totalFeatures, info->curDataFile);
    if (nread != totalFeatures) {
      error("HTKFile: read %u items, expected %u in observation file '%s'\n",
	    nread, totalFeatures, info->fofName);
    }
    if (info->swap()) {
      for (unsigned i=0; i<totalFeatures; i+=1) {
	tmpBuf[i]=swapb_short_short(tmpBuf[i]);
      }
    }
    if (numDiscrete() > 0) { // 16-bit integers
      for (unsigned i=0; i<totalFeatures; i+=1) {
	buffer[i]=(Data32)tmpBuf[i];
      }
    } else { // compressed 16-bit fixed point
      float *floatBuf = (float *)buffer;
      for (unsigned i=0; i<totalFeatures; i+=1) {
	floatBuf[i]=(float)tmpBuf[i];
      }
    }
    if (numContinuous()) { // uncompress
      for (unsigned i=0; i < count; i+=1) {
	float* curSampPtr=(float *)buffer+i*featuresPerFrame;
	copy_add_vf32_vf32(featuresPerFrame,htkInfo->offset,curSampPtr);
	copy_div_vf32_vf32(featuresPerFrame,htkInfo->scale, curSampPtr);	
      }
    }
    delete [] tmpBuf;
    return buffer;
  }
  if (numContinuous() > 0 && !htkInfo->isCompressed) {
    float *tmpBuf = (float *)buffer;
    nread = fread(tmpBuf,sizeof(float),totalFeatures, info->curDataFile);
    if (nread != totalFeatures) {
      error("HTKFile: read %i items, expected %u from observation file '%s'\n", 
	    nread, totalFeatures, info->fofName);
    }
    // swap if needed.
    if(info->swap()) {
      float tmp_float[1];
      for (unsigned i=0; i<totalFeatures; i+=1) {
	swapb_vf32_vf32(1,(tmpBuf+i),tmp_float);
	tmpBuf[i]=tmp_float[0];
      }
    }
    return buffer;
  }
  
  assert(false); // should never get here!
  return NULL;
}


static
void 
printHTKHeader(FILE* ofp, bool oswap, int numInts, int numFloats, int numSamples) {

  DBGFPRINTF((stderr,"obsPrint: Printing HTK header.\n"));
  // structure of HTK header
  //  Int32 numSamples;
  Int32 samplePeriod=1;
  short parameterKind;
  short sampleSize;
  if(numFloats > 0) {
    sampleSize=numFloats*sizeof(float);
    parameterKind=USER;  // enum in GMTK_Stream.h
  } else {
    sampleSize=numInts*sizeof(short);
    parameterKind=DISCRETE;  // enum in GMTK_Stream.h
  }

  if (oswap) {
    numSamples = swapb_i32_i32(numSamples);
    samplePeriod = swapb_i32_i32(samplePeriod);
    sampleSize = swapb_short_short(sampleSize);
    parameterKind = swapb_short_short(parameterKind);
  }

  if (fwrite(&numSamples,sizeof(Int32),1,ofp) != 1) {
    error("Cannot write HTK number of samples\n");
  }
  if (fwrite((Int32 *)&samplePeriod,sizeof(Int32),1,ofp) != 1) {
    error("Cannot write HTK sample period\n");
  }

  if (fwrite((short *)&sampleSize,sizeof(short),1,ofp) != 1) {
    error("Cannot write HTK sample size\n");
  }

  if (fwrite(&parameterKind,sizeof(short),1,ofp) != 1) {
    error("Cannot write HTK parm kind\n");
  }

  DBGFPRINTF((stderr,"obsPrint: Finished printing HTK header.\n"));
}


void 
HTKFile::writeSegment(Data32 const *segment, unsigned nFrames) {
  assert(currFeature == 0);
  if (writeFile) {
    if (fclose(writeFile)) {
      error("ERROR closing output file\n");
    }
    writeFile = NULL;
  }
  char* current_output_fname = new char[strlen(outputFileName)+strlen(outputNameSeparatorStr)+50];
  sprintf(current_output_fname,"%s%s%d",outputFileName,outputNameSeparatorStr,currSegment);
  if ((writeFile = fopen(current_output_fname, "wb")) == NULL) {
    error("Couldn't open output file '%s' for writing.",current_output_fname);
  }
  fprintf(listFile,"%s\n",current_output_fname);

  void (*copy_swap_func_ptr)(size_t, const intv_int32_t*, intv_int32_t*)=NULL;
  if(oswap) {
    copy_swap_func_ptr=&swapb_vi32_vi32;
  } else {
    copy_swap_func_ptr=&copy_vi32_vi32;
  }
  
  printHTKHeader(writeFile, oswap, _numDiscreteFeatures, _numContinuousFeatures, nFrames);

  //////////  Print the frames ////////////////////////////////////
  float  *cont_buf_p = (float  *) segment;
  UInt32 *disc_buf_p = (UInt32 *) segment + _numContinuousFeatures;
  size_t fwrite_result;

  for (unsigned frame_no=0; frame_no < nFrames ; frame_no+=1) {

    /// Print continuous part of frame /////////////////////////////
    for (unsigned frit=0; frit < _numContinuousFeatures; frit+=1) {
      DBGFPRINTF((stderr,"obsPrint: Printing HTK float %f.\n", cont_buf_p[frit]));
      copy_swap_func_ptr(1, (int*)&cont_buf_p[frit],(int*)&cont_buf_p[frit]);
      fwrite_result = fwrite(&cont_buf_p[frit], sizeof(cont_buf_p[frit]), 1, writeFile);
      if (fwrite_result != 1) {
	error("Error writing to output file '%s'\n", current_output_fname);
      }
    }
    ///////////////////////////////////////////////////////////////
    
    /// Print discrete part of the frame ///////////////////////////
    for (unsigned lrit=0; lrit < _numDiscreteFeatures; lrit+=1) {
      if (_numContinuousFeatures > 0) {
	copy_swap_func_ptr(1, (int*)&disc_buf_p[lrit],(int*)&disc_buf_p[lrit]);
	fwrite_result = fwrite(&disc_buf_p[lrit], sizeof(disc_buf_p[lrit]), 1, writeFile);
	if (fwrite_result != 1) {
	  error("Error writing to output file '%s'\n", current_output_fname);
	}
      } else if (_numContinuousFeatures == 0) { // in the HTK format we
	// cannot mix floats with discrete data; that's why if there is at
	// least one float component everyting is written out as a float.
	short short_lab_buf_p=(short)disc_buf_p[0];
	DBGFPRINTF((stderr,"obsPrint: Printing HTK short %d.\n", short_lab_buf_p));
	if (oswap) {
	  short_lab_buf_p = swapb_short_short(short_lab_buf_p);
	}
	fwrite_result = fwrite(&short_lab_buf_p, sizeof(short_lab_buf_p), 1, writeFile);
	if (fwrite_result != 1) {
	  error("Error writing to output file '%s'\n", current_output_fname);
	}
      }
    } // end of for (unsigned lrit=0;lrit<num_discrete; ++lrit)

    cont_buf_p += _numFeatures;
    disc_buf_p += _numFeatures;
  }  // end of for (unsigned frame_no=0; frame_no < num_frames ; ++frame_no)
  //////////////////////////////////////////////////////////////////////////// 

  if (fclose(writeFile)) {
    error("ERROR closing output file '%s'\n", current_output_fname);
  }
  writeFile = NULL;
  delete []  current_output_fname;

  currSegment += 1;
  currFrame = 0;
  currFeature = 0;
}


void 
HTKFile::setFrame(unsigned frame) {
}


void 
HTKFile::writeFrame(Data32 const *frame) {
  assert(currFeature == 0);
  void (*copy_swap_func_ptr)(size_t, const intv_int32_t*, intv_int32_t*)=NULL;
  if(oswap) {
    copy_swap_func_ptr=&swapb_vi32_vi32;
  } else {
    copy_swap_func_ptr=&copy_vi32_vi32;
  }
  
  if (currFrame == 0) {
    assert(!writeFile); // previous EoS (or ctor) should have closed it
    char* current_output_fname = new char[strlen(outputFileName)+strlen(outputNameSeparatorStr)+50];
    sprintf(current_output_fname,"%s%s%d",outputFileName,outputNameSeparatorStr,currSegment);
    if ((writeFile = fopen(current_output_fname, "wb")) == NULL) {
      error("Couldn't open output file '%s' for writing.",current_output_fname);
    }
    fprintf(listFile,"%s\n",current_output_fname);
    delete [] current_output_fname;
    // write fake header - EoS will seek back & wirte correct header
    printHTKHeader(writeFile, oswap, _numDiscreteFeatures, _numContinuousFeatures, 1);
  }
  float  *cont_buf_p = (float  *) frame;
  UInt32 *disc_buf_p = (UInt32 *) frame + _numContinuousFeatures;
  size_t fwrite_result;

  /// Print continuous part of frame /////////////////////////////
  for (unsigned frit=0; frit < _numContinuousFeatures; frit+=1) {
    DBGFPRINTF((stderr,"obsPrint: Printing HTK float %f.\n", cont_buf_p[frit]));
    copy_swap_func_ptr(1, (int*)&cont_buf_p[frit],(int*)&cont_buf_p[frit]);
    fwrite_result = fwrite(&cont_buf_p[frit], sizeof(cont_buf_p[frit]), 1, writeFile);
    if (fwrite_result != 1) {
      error("Error writing to output file\n");
    }
  }
  ///////////////////////////////////////////////////////////////
  
  /// Print discrete part of the frame ///////////////////////////
  for (unsigned lrit=0; lrit < _numDiscreteFeatures; lrit+=1) {
    if (_numContinuousFeatures > 0) {
      copy_swap_func_ptr(1, (int*)&disc_buf_p[lrit],(int*)&disc_buf_p[lrit]);
      fwrite_result = fwrite(&disc_buf_p[lrit], sizeof(disc_buf_p[lrit]), 1, writeFile);
      if (fwrite_result != 1) {
	error("Error writing to output file\n");
      }
    } else if (_numContinuousFeatures == 0) { // in the HTK format we
      // cannot mix floats with discrete data; that's why if there is at
      // least one float component everyting is written out as a float.
      short short_lab_buf_p=(short)disc_buf_p[0];
      DBGFPRINTF((stderr,"obsPrint: Printing HTK short %d.\n", short_lab_buf_p));
      if (oswap) {
	short_lab_buf_p = swapb_short_short(short_lab_buf_p);
      }
      fwrite_result = fwrite(&short_lab_buf_p, sizeof(short_lab_buf_p), 1, writeFile);
      if (fwrite_result != 1) {
	error("Error writing to output file\n");
      }
    }
  } // end of for (unsigned lrit=0;lrit<num_discrete; ++lrit) 
  currFrame += 1;
  currFeature = 0;
}


void 
HTKFile::writeFeature(Data32 x) {
  void (*copy_swap_func_ptr)(size_t, const intv_int32_t*, intv_int32_t*)=NULL;
  if(oswap) {
    copy_swap_func_ptr=&swapb_vi32_vi32;
  } else {
    copy_swap_func_ptr=&copy_vi32_vi32;
  }
  
  if (currFrame == 0 && currFeature == 0) {
    assert(!writeFile); // previous EoS (or ctor) should have closed it
    char* current_output_fname = new char[strlen(outputFileName)+strlen(outputNameSeparatorStr)+50];
    sprintf(current_output_fname,"%s%s%d",outputFileName,outputNameSeparatorStr,currSegment);
    if ((writeFile = fopen(current_output_fname, "wb")) == NULL) {
      error("Couldn't open output file '%s' for writing.",current_output_fname);
    }
    fprintf(listFile,"%s\n",current_output_fname);
    delete [] current_output_fname;
    // write fake header - EoS will seek back & wirte correct header
    printHTKHeader(writeFile, oswap, _numDiscreteFeatures, _numContinuousFeatures, 1);
  }
  size_t fwrite_result;

  /// Print continuous part of frame /////////////////////////////
  if (currFeature < _numContinuousFeatures) {
    copy_swap_func_ptr(1, (intv_int32_t*)&x,(intv_int32_t*)&x);
    fwrite_result = fwrite(&x, sizeof(x), 1, writeFile);
    if (fwrite_result != 1) {
      error("Error writing to output file\n");
    }
  } else {
    /// Print discrete part of the frame ///////////////////////////
    if (_numContinuousFeatures > 0) {
      copy_swap_func_ptr(1, (intv_int32_t*)&x,(intv_int32_t*)&x);
      fwrite_result = fwrite(&x, sizeof(x), 1, writeFile);
      if (fwrite_result != 1) {
	error("Error writing to output file\n");
      }
    } else if (_numContinuousFeatures == 0) { // in the HTK format we
      // cannot mix floats with discrete data; that's why if there is at
      // least one float component everyting is written out as a float.
      short short_lab_buf_p=(short)x;
      if (oswap) {
	short_lab_buf_p = swapb_short_short(short_lab_buf_p);
      }
      fwrite_result = fwrite(&short_lab_buf_p, sizeof(short_lab_buf_p), 1, writeFile);
      if (fwrite_result != 1) {
	error("Error writing to output file\n");
      }
    }
  } 

  currFeature += 1;
  if (currFeature == _numFeatures) {
    currFrame   += 1;
    currFeature  = 0;
  }
}


void 
HTKFile::endOfSegment() {
  assert(currFeature == 0);
  if (currFrame == 0) {
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
    // Now that we know the number of frames, re-do the header correctly
    if (fseek(writeFile, 0, SEEK_SET)) {
      error("ERROR seeking in output file\n");
    }
  }
  printHTKHeader(writeFile, oswap, _numDiscreteFeatures, _numContinuousFeatures, currFrame);

  if (fclose(writeFile)) {
    error("ERROR closing output file\n");
  }
  writeFile = NULL;
  currSegment += 1;
  currFrame = 0;
  currFeature = 0;
}

