
/*
 * GMTK_Filter.cc
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
#include <assert.h>
#include <string.h>

#include "error.h"
#include "general.h"

#include "file_utils.h"
#include "GMTK_Filter.h"

#include "GMTK_FilterFile.h"
#include "GMTK_FIRFilter.h"
#include "GMTK_ARMAFilter.h"
#include "GMTK_AffineFilter.h"
#include "GMTK_UpsampleFilter.h"
#include "GMTK_UpsampleSmoothFilter.h"

subMatrixDescriptor *subMatrixDescriptor::freeList = NULL;

static long parse_long(const char*const s) {
    size_t len = strlen(s);
    char *ptr;
    long val;

    val = strtol(s, &ptr, 0);

    if (ptr != (s+len))
        error("Not an integer argument.");

    return val;
}


static double parse_float(const char*const s) {
    size_t len = strlen(s);
    char *ptr;
    double val;
    val = strtod(s, &ptr);
    if (ptr != (s+len))
        error("Not an floating point argument.");
    return val;
}


static double conv2double(char* str, unsigned& len, char delimiter,bool conv2int=false) {
#define NUMBER_STRING_MAX_LEN 20
  char* str_ptr=str;
  char new_str[NUMBER_STRING_MAX_LEN];
  double val;

  DBGFPRINTF((stderr,"str=%s\n",str));

  int i=0;
  while(*str_ptr!='_' && *str_ptr!='\0') {
    new_str[i]=*str_ptr; str_ptr++; i++;
    if(i>=NUMBER_STRING_MAX_LEN)
      error("ERROR: While reading observation files, a number string was too long while converting to float\n");
  }

  new_str[i]='\0';
  if(conv2int)
    val=parse_long(new_str);
  else
    val=parse_float(new_str);

  len=i;

  return (val);

}


int 
parseTransform(char*& trans_str, int& magic_int, double& magic_double, char *&filterFileName) {
  if(*trans_str=='\0')
    return (END_STR);

  unsigned len;
  unsigned return_val;
  // remove leading spaces
  trans_str += strspn(trans_str, " ");

  char c=*trans_str;
  DBGFPRINTF((stderr,"In parseTransform: c=%c\n",c));
  switch(c) {
    //see GMTK_ObservationMatrix.h for definitions of LETTERs below
  case TRANS_NORMALIZATION_LETTER:   // 'N'
    ++trans_str;
    if(*trans_str=='_') ++trans_str;  // get rid of separator
    return NORMALIZE;
  case TRANS_MEAN_SUB_LETTER:  // 'E'
    ++trans_str;
    if(*trans_str=='_') ++trans_str;  // get rid of separator
    return MEAN_SUB;
  case TRANS_MULTIPLICATION_LETTER: // 'M'
    ++trans_str;
    // get multiplier
    //DBGFPRINTF((stderr,"trans_str=%s\n",trans_str));
    magic_double=conv2double(trans_str,len,'_');
    if(len==0)
      error("ERROR: In parsing tranforms: Need to supply multiplicative factor with the 'M' transformation.\n");
    trans_str+=len;
    if(*trans_str=='_') ++trans_str;  // get rid of separator
    return MULTIPLY;
  case TRANS_OFFSET_LETTER: // 'O'
    ++trans_str;
    // get offset
    magic_double=conv2double(trans_str,len,'_');
    if(len==0)
      error("ERROR: In parsing tranforms: Need to supply offset with 'O' transformation.\n");
    trans_str+=len;
    if(*trans_str=='_') ++trans_str;  // get rid of separator
    return OFFSET;
  case TRANS_UPSAMPLING_LETTER: // 'U'
    // an upsampling style should follow
    ++trans_str;
    if(*(trans_str)==TRANS_SMOOTH_LETTER)  // 'S'
      return_val=UPSAMPLE_SMOOTH;
    else if(*(trans_str)==TRANS_HOLD_LETTER)  // 'H'
      return_val=UPSAMPLE_HOLD;
    else
      return UNRECOGNIZED_TRANSFORM;
    // a number to upsample by should follow
    ++trans_str;
    magic_double=conv2double(trans_str,len,'_');
    if(len==0)
      error("ERROR: In parsing tranforms: Need to supply upsampling factor with 'UH' or 'US' transformations.\n");
    trans_str+=len;
    if(*trans_str=='_') ++trans_str;  // get rid of separator
    return return_val;
  case TRANS_ARMA_LETTER: // 'R'
    ++trans_str;
    //DBGFPRINTF((stderr,"trans_str=%s\n",trans_str));
    // get order of ARMA filter
    magic_int=(int)conv2double(trans_str,len,'_',true); // conv2int is true
    if(len==0)
      error("ERROR: In parsing tranforms: Need to supply order of arma filter with 'R' transformation.\n");
    trans_str+=len;
    if(*trans_str=='_') ++trans_str;
    return ARMA;
  case TRANS_AFFINE_LETTER: // 'A'
  case FILTER_LETTER:       // 'F'
    ++trans_str;
    if(*(trans_str)=='@') {
      magic_int=-1;  //we'll read the filter from a file
#define MAX_TMP_STRING_LEN 200
      char tmpString[MAX_TMP_STRING_LEN];
      ++trans_str;
      unsigned i=0;
      while(*trans_str != '\0' && *trans_str != '_') {
	tmpString[i++]=*trans_str;
	++trans_str;
      }
      tmpString[i]='\0';
      DBGFPRINTF((stderr,"Read filter file name '%s'\n",tmpString));
      unsigned tmpStringLen=strlen(tmpString);
      if(tmpStringLen==0) error("ERROR: In parsing tranforms: no filter file name specified\n");
      filterFileName=new char[tmpStringLen];
      strcpy(filterFileName,tmpString);
      if(*trans_str=='_') ++trans_str;
    }
    else if(*trans_str=='_') ++trans_str;

    return c == TRANS_AFFINE_LETTER ? AFFINE : FILTER;
  case NONE_LETTER:
    ++trans_str;
    if(*trans_str=='_') ++trans_str;
    return NONE;
  default:
    error("ERROR: In parsing tranforms: Unrecognized transformation substring (%s)\n",trans_str);
    //DBGFPRINTF((stderr,"In parseTransform: Unrecognized transform @ trans_str=%s\n",trans_str));
    //return UNRECOGNIZED_TRANSFORM;

  }
  return UNRECOGNIZED_TRANSFORM;
}


Filter *
instantiateFilters(char *filterStr, unsigned numContinuous, unsigned numDiscrete) {
  int magicInt;
  double magicDouble;
  char * filterFileName;
  int xform;
  float  *B;
  double *cc;
  
  Filter *xformer = NULL;
  Filter *firstFilter = NULL;
  if (filterStr) {
    while((xform=parseTransform(filterStr, magicInt, magicDouble, filterFileName)) != END_STR) {
      switch (xform) {
      case NONE: 
	xformer = new Filter(NULL); // identity filter
	break;
      case UPSAMPLE_HOLD:
	xformer = new UpsampleFilter((unsigned)magicDouble);
	break;
      case  UPSAMPLE_SMOOTH:
	xformer = new UpsampleSmoothFilter((unsigned)magicDouble);
	break;
#if 0
      // unsupported
      case DOWNSAMPLE:
	printf("downsample stream %u\n", i);
	break;
      case DELTAS:
	printf("deltas stream %u\n", i);
	break;
      case DOUBLE_DELTAS:
	printf("double deltas stream %u\n", i);
	break;
#endif
      case MULTIPLY:
	//	printf("multiply by %f (%u)\n", magicDouble, numContinuous);
	B = new float[numContinuous];
	for (unsigned j=0; j < numContinuous; j+=1)
	  B[j] = (float)magicDouble;
	xformer = new FIRFilter(0, numContinuous, B, NULL, NULL);
	break;
      case NORMALIZE:
	error("normalize stream - now handled by affine\n");
	break;
      case MEAN_SUB:
	error("mean sub stream - now handled by affine\n");
	break;
      case ARMA:
	xformer = new ARMAFilter(magicInt, numContinuous);
	break;
      case AFFINE:
	//printf("affine with %s\n", filterFileName);
	xformer = new AffineFilter(filterFileName, NULL);
	xformer->numOutputFeatures(numContinuous, numDiscrete, numContinuous, numDiscrete);
	break;
      case FILTER:
	//printf("filter stream with %s\n", filterFileName);
	xformer = new FIRFilter(filterFileName, NULL);
	break;
      case OFFSET:
	//printf("offset stream by %f\n", magicDouble);
	cc = new double[numContinuous];
	for (unsigned j=0; j < numContinuous; j+=1)
	  cc[j] = magicDouble;
	xformer = new AffineFilter(numContinuous, numContinuous, NULL, cc, NULL);
	break;
      default:
	error("ERROR: unknown transform '%s'\n", filterStr);
      }
      if (!firstFilter) {
	firstFilter = xformer;
      } else {
	firstFilter->appendFilter(xformer);
      }
    }
  }
  if (firstFilter)
    return firstFilter;
  return new Filter(NULL);
}
