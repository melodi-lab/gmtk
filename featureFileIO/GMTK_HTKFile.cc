
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
		  string& fnameStr) 
{
  size_t fNameLen;
  *startFrame=0, *endFrame=-1; //these are the right values if the frame range is not specified
  if(sentLoc[sentLoc.length()-1]==']'){
    //have a subrange spec
    fNameLen=sentLoc.find_last_of('[');
    if (fNameLen==string::npos){
      error("ERROR: parseSentenceSpec: '%s' is an invalid sentence location.  Must be of the form 'filename[startFrame:endFrame]'",sentLoc.c_str());		  
    }

    string range= sentLoc.substr(fNameLen+1,sentLoc.length()-2-fNameLen);
    if (sscanf(range.c_str(),"%d:%d",startFrame,endFrame) != 2)
      error("ERROR: parseSentenceSpec: '%s' is an invalid sentence location.  Must be of the form 'filename[startFrame:endFrame]'",sentLoc.c_str());		  
    
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
    error("ERROR: openHTKFile: Can't open '%s' for input\n",fname.c_str());
  }

  if (fread(&tmp1,sizeof(Int32),1,f->curDataFile) != 1) {
    error("ERROR: openHTKFile: Can't read number of samples\n");
  }


  if (fread((Int32 *)&tmp2,sizeof(Int32),1,f->curDataFile) != 1) {
    error("ERROR: openHTKFile: Can't read sample period\n");
  }

  if (fread((short *)&stmp1,sizeof(short),1,f->curDataFile) != 1) {
    error("ERROR: openHTKFile: Can't read sample size\n");
    return 0;
  }

  if (fread(&stmp2,sizeof(short),1,f->curDataFile) != 1) {
    error("ERROR: openHTKFile: Can't read parm kind\n");
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
    error("ERROR: openHTKFile: number of samples is %i\n",n_samples);
  }

  if (samp_period <= 0 || samp_period > 1000000) {
    warning("WARNING: openHTKFile: sample period is %i - must be between 0 and 1000000\n", samp_period);
  }

  if (samp_size <= 0 || samp_size > 5000) {
    warning("WARNING: openHTKFile: sample size is %i - must be between 0 and 5000\n",samp_size);
  }

  short pk = parm_kind & BASEMASK;
  bool isCompressed = parm_kind & IS_COMPRESSED;

  if (pk < WAVEFORM || pk > ANON) {
    warning("WARNING: openHTKFile: Undefined parameter kind for HTK feature file: %i. Will assume float features.\n",pk);
  }

  // For now we don't support the WAVEFORM and IREFC parameter kind.  It uses
  // shorts instead of floats and that requires special treatment.
  if (pk == WAVEFORM) {
    warning("WARNING: openHTKFile: HTK WAVEFORM parameter kind not supported: %i\n",pk);
  }
  else if (pk == IREFC) {
    warning("WARNING: openHTKFile: HTK IREFC parameter kind not supported: %i\n",pk);
  }

  int n_fea;

  // parameter kind DISCRETE = all discrete features

  if (nfloats == 0) {
    // we divide by sizeof(short) which should be presumably size of a 2-byte int
    // (but this is MACHINE DEPENDENT).
    n_fea = samp_size / sizeof(short) ;
    if (n_fea != nints) {
      error("ERROR: openHTKFile:  Number of features in file (%i) does not match number of ints specified (%i)\n", n_fea,nints);
    }
    if(parm_kind != DISCRETE) {
      warning("WARNING: openHTKFile: Number of floats specified is 0 but the HTK parameter kind is not DISCRETE.\n");
    }
  }
  // otherwise all continuous features
  else {
    if(isCompressed)
      n_fea = samp_size / sizeof(short);
    else
      n_fea = samp_size / sizeof(float);
		
    if (n_fea != nfloats) {
      error("ERROR: openHTKFile:  Number of features in file (%i) does not match number of floats specified (%i)\n", n_fea,nfloats);
    }
    if(parm_kind == DISCRETE) {
      warning("WARNING: openHTKFile:  Number of floats specified (%i) is not 0 but the HTK parameter kind is DISCRETE.\n",nfloats);
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
      error("ERROR: openHTKFile: Can't read scales for decompressing.\n");
    }
    if(bswap)
      swapb_vf32_vf32(n_fea,tmp,scale);
    else
      copy_vf32_vf32(n_fea,tmp,scale);

    if (fread(tmp,sizeof(float),n_fea,f->curDataFile) != (unsigned short)n_fea) {
      error("ERROR: openHTKFile: Can't read offsets for decompressing.\n");
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
    error("ERROR: openHTKFile: Requested segment no %li of observation file '%s' but the max num of segments in list of HTK files is %li",sentno,f->fofName,htkfile_size);
  }

  //  assert(sentno >= 0 && sentno < _numSegments);
  if (f->dataNames[sentno] == NULL) {
    error("ERROR: openHTKFile: Filename is NULL for segment %li\n",f->dataNames[sentno]);
  }

  
  int startFrame, endFrame;
  string fnameStr;
  parseSentenceSpec(f->dataNames[sentno], &startFrame, &endFrame, fnameStr);

  if(f->curDataFilename != fnameStr && f->curDataFile){
    DBGFPRINTF((stderr,"In openHTKFile, f->curDataFilename  %s fnameStr %s f->curDataFile  %d \n", f->curDataFilename.c_str(), fnameStr.c_str(), f->curDataFile));
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
    error("ERROR: openHTKFile: %s has the last frame at %d, beyond %d, which is the number of frames in file.\n",f->dataNames[sentno], endFrame, htkInfo->n_samples);	  
  
  f->curNumFrames=endFrame-startFrame+1;
  
  
  DBGFPRINTF((stderr,"In openHTKFile, curNumFrames %d\n",f->curNumFrames));
  
  //now we seek to the start frame 
  if (gmtk_fseek(f->curDataFile, (gmtk_off_t)(htkInfo->startOfData+startFrame*htkInfo->samp_size), SEEK_SET)) {
    error("ERROR: openHTKFile: fseek() failed for '%s'", f->dataNames[sentno]);
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
    error("HTKFile: fseek() failed for '%s'", info->fofName);
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
      error("HTKFile: read %u items, expected %u", nread, totalFeatures);
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
      error("HTKFile: read %i items, expected %i", nread,totalFeatures);
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

