
/*
 * GMTK_ASCIIFile.cc
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

#include <stdio.h>
#include <string.h>

#include "error.h"
#include "general.h"

#include "file_utils.h"
#include "GMTK_ASCIIFile.h"

using namespace std;

ASCIIFile::ASCIIFile(const char *name, unsigned nfloats, unsigned nints, 
		     unsigned num, bool cppIfAscii, char const* cppCommandOptions,
		     char const *contFeatureRangeStr_, 
		     char const *discFeatureRangeStr_, 
		     char const *preFrameRangeStr_, 
		     char const *segRangeStr_)
  : ObservationFile(contFeatureRangeStr_, 
		    discFeatureRangeStr_, 
		    preFrameRangeStr_,
		    segRangeStr_),
    cppIfAscii(cppIfAscii),
    cppCommandOptions(cppCommandOptions)
{
  buffer = NULL;
  bufferSize = 0;
  if (name == NULL) 	
    error("ASCIIFile: File name is NULL for stream %i\n",num);	
  if (nfloats == 0 && nints == 0)
    error("ASCIIFile: number of float and int features cannot both be zero");

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
    error("ASCIIFile: couldn't open '%s' for reading\n", fofName);

  // for some reason this asserts while the following doesn't
  // assert(numFileNames == readFof(fofFile));
  
  unsigned rfof = readFof(fofFile, fofName, numFileNames, dataNames,
			  cppIfAscii, cppCommandOptions);
  assert(numFileNames == rfof);

  closeCPPableFile(fofFile, cppIfAscii);
}


// Begin sourcing data from the requested segment.
// Must be called before any other operations are performed on a segment.
// Loads the segment's data into the buffer and set the # of frames.
bool
ASCIIFile::openSegment(unsigned seg) {
  assert(seg < numFileNames);
  bool found_new_line=true;
  char *fname = dataNames[seg];
  unsigned n_samples = 0;
  char ch;
  FILE *curDataFile = NULL;

  if (fname == NULL) {
    warning("ASCIIFile::openSegment: Filename is NULL for segment %u\n",seg);
    return false;
  }

  curDataFile = openCPPableFile(fname, cppIfAscii, cppCommandOptions);
  if (curDataFile == NULL) {
    warning("ASCIIFile::openSegment: Can't open '%s' for input\n",fname);
    return false;
  }
  
  // Read through the segment once to determine its length, then again 
  // to actually load the data.

  // For ASCII, newline is record delimiter - additional or missing nl's
  // will cause error messages.
  int tmp;
  while ((tmp = fgetc(curDataFile)) != EOF) {
    ch=tmp;
    if(ch=='\n') { found_new_line=true; continue;}
    if(ch==CPP_DIRECTIVE_CHAR) { found_new_line=false; continue; }
    if(found_new_line) {
      found_new_line=false;
      n_samples++;
    }
  }
  closeCPPableFile(curDataFile, cppIfAscii);
  curDataFile = openCPPableFile(fname, cppIfAscii, cppCommandOptions);
  if (curDataFile == NULL) {
    warning("ASCIIFile::openSegment: Can't open '%s' for input\n",fname);
    return false;
  }

  nFrames = n_samples;
  if (n_samples > bufferSize) {
    buffer = (Data32 *) realloc(buffer, n_samples * _numFeatures * sizeof(Data32));
    bufferSize = n_samples;
  }
  Data32 *dest = buffer;

  int lineNum=0;

  // consume CPP special directives if any
#ifdef PIPE_ASCII_FILES_THROUGH_CPP
  if(cppIfAscii) {
    while((tmp=fgetc(curDataFile))==CPP_DIRECTIVE_CHAR) {
      while((tmp=fgetc(curDataFile))!='\n');
      lineNum++;
    }
    ungetc(tmp,curDataFile);
  }
#endif

  // could be made a bit more efficient since we check whether
  // num_floats and num_ints > 0 for each frame.
  for(unsigned s=0; s < n_samples; ++s) {
    lineNum++;
    for (unsigned n = 0; n < _numContinuousFeatures; n+=1) {
      if (fscanf(curDataFile,"%e", (float *)(dest++)) != 1) {
	error("ERROR: ASCIIFile::openSegment: couldn't read %u'th item in frame %u\n",n,s);
      }
    }
    for (unsigned n = 0; n < _numDiscreteFeatures; n+=1) {
      if (fscanf(curDataFile,"%d", (Int32 *)(dest++)) != 1) {
	error("ERROR: ASCIIFile::openSegment: couldn't read %u'th item in frame %u\n",n,s);
      }
    }
  }
  closeCPPableFile(curDataFile, cppIfAscii);
  return true;
}

