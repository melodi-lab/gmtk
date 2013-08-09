
/*
 * GMTK_ASCIIFile.cc
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (C) 2011 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * 
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
		     char const *segRangeStr_,
		     unsigned leftPad, unsigned rightPad)
  : ObservationFile(name, num,
		    contFeatureRangeStr_, 
		    discFeatureRangeStr_, 
		    preFrameRangeStr_,
		    segRangeStr_, 
		    leftPad, rightPad),
    cppIfAscii(cppIfAscii),
    cppCommandOptions(cppCommandOptions)
{
  buffer = NULL;
  bufferSize = 0;
  if (name == NULL) 	
    error("ASCIIFile: File name is NULL for observation file %u\n",num);	
  if (nfloats == 0 && nints == 0)
    error("ASCIIFile: observation file %u '%s': number of float and int features cannot both be zero\n", num, name);

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
    error("ASCIIFile: observation file %u: couldn't open '%s' for reading\n", num, fofName);

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
    warning("ASCIIFile::openSegment: observation file %u '%s': filename is NULL for segment %u\n",
	    observationFileNum, observationFileName, seg);
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
    for (unsigned n = 0; n < _numContinuousFeatures; n+=1, dest+=1) {
      if (fscanf(curDataFile,"%e", (float *)(dest)) != 1) {
	error("ERROR: ASCIIFile::openSegment: observation file %u '%s' segment %u: couldn't read %u'th item in frame %u\n",
	      observationFileNum, observationFileName, seg, n, s);
      }
    }
    for (unsigned n = 0; n < _numDiscreteFeatures; n+=1, dest+=1) {
      if (fscanf(curDataFile,"%d", (Int32 *)(dest)) != 1) {
	error("ERROR: ASCIIFile::openSegment: observation file %u '%s' segment %u: couldn't read %u'th item in frame %u\n",
	      observationFileNum, observationFileName, seg, _numContinuousFeatures+n, s);
      }
    }
  }
  closeCPPableFile(curDataFile, cppIfAscii);
  return true;
}


void 
ASCIIFile::writeSegment(Data32 const *segment, unsigned nFrames) {
  assert(currFeature == 0);
  if (writeFile) {
    if (fclose(writeFile)) {
      error("ERROR closing output file\n");
    }
    writeFile = NULL;
  }
  char* current_output_fname = new char[strlen(outputFileName)+strlen(outputNameSeparatorStr)+50];
  sprintf(current_output_fname,"%s%s%d",outputFileName,outputNameSeparatorStr,currSegment);
  if ((writeFile = fopen(current_output_fname, "w")) == NULL) {
    error("Couldn't open output file '%s' for writing.",current_output_fname);
  }
  assert(listFile);
  fprintf(listFile,"%s\n",current_output_fname);

  //////////  Print the frames ////////////////////////////////////
  float  *cont_buf_p = (float  *) segment;
  UInt32 *disc_buf_p = (UInt32 *) segment + _numContinuousFeatures;

  for (unsigned frame_no=0; frame_no < nFrames ; frame_no+=1) {

    /// Print continuous part of frame /////////////////////////////
    for (unsigned frit=0; frit < _numContinuousFeatures; frit+=1) {
      fprintf(writeFile, " %.8e", cont_buf_p[frit]);
    }
    ///////////////////////////////////////////////////////////////
    
    /// Print discrete part of the frame ///////////////////////////
    for (unsigned lrit=0; lrit < _numDiscreteFeatures; lrit+=1) {
      fprintf(writeFile, " %u", disc_buf_p[lrit]);
    }
    fprintf(writeFile, "\n");

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
ASCIIFile::writeFrame(Data32 const *frame) {
  assert(currFeature == 0);
  if (currFrame == 0) {
    assert(listFile);
    assert(!writeFile); // previous EoS (or ctor) should have closed it
    char* current_output_fname = new char[strlen(outputFileName)+strlen(outputNameSeparatorStr)+50];
    sprintf(current_output_fname,"%s%s%d",outputFileName,outputNameSeparatorStr,currSegment);
    if ((writeFile = fopen(current_output_fname, "w")) == NULL) {
      error("Couldn't open output file '%s' for writing.",current_output_fname);
    }
    fprintf(listFile,"%s\n",current_output_fname);
    delete [] current_output_fname;
  }
  float  *cont_buf_p = (float  *) frame;
  UInt32 *disc_buf_p = (UInt32 *) frame + _numContinuousFeatures;

  /// Print continuous part of frame /////////////////////////////
  for (unsigned frit=0; frit < _numContinuousFeatures; frit+=1) {
    fprintf(writeFile, " %.8e", cont_buf_p[frit]);
  }
  ///////////////////////////////////////////////////////////////
  
  /// Print discrete part of the frame ///////////////////////////
  for (unsigned lrit=0; lrit < _numDiscreteFeatures; lrit+=1) {
    fprintf(writeFile, " %u", disc_buf_p[lrit]);
  }
  fprintf(writeFile, "\n");
  currFrame += 1;
  currFeature = 0;
}


void 
ASCIIFile::writeFeature(Data32 x) {
  if (currFrame == 0 && currFeature == 0) {
    assert(listFile);
    assert(!writeFile); // previous EoS (or ctor) should have closed it
    char* current_output_fname = new char[strlen(outputFileName)+strlen(outputNameSeparatorStr)+50];
    sprintf(current_output_fname,"%s%s%d",outputFileName,outputNameSeparatorStr,currSegment);
    if ((writeFile = fopen(current_output_fname, "w")) == NULL) {
      error("Couldn't open output file '%s' for writing.",current_output_fname);
    }
    fprintf(listFile,"%s\n",current_output_fname);
    delete [] current_output_fname;
  }

  union {
    float  f;
    UInt32 i;
  } fea;
  fea.i = (UInt32) x;
  /// Print continuous part of frame /////////////////////////////
  if (currFeature < _numContinuousFeatures) {
    fprintf(writeFile, " %.8e", fea.f);
  } else {
    /// Print discrete part of the frame ///////////////////////////
    fprintf(writeFile, " %u", fea.i);
  }
  
  currFeature += 1;
  if (currFeature == _numFeatures) {
    fprintf(writeFile, "\n");
    currFrame += 1;
    currFeature = 0;
  }
}


void 
ASCIIFile::endOfSegment() {
  assert(currFeature == 0);
  if (currFrame == 0) {
    assert(listFile);
    assert(!writeFile); // previous EoS (or ctor) should have closed it
    char* current_output_fname = new char[strlen(outputFileName)+strlen(outputNameSeparatorStr)+50];
    sprintf(current_output_fname,"%s%s%d",outputFileName,outputNameSeparatorStr,currSegment);
    if ((writeFile = fopen(current_output_fname, "w")) == NULL) {
      error("Couldn't open output file '%s' for writing.",current_output_fname);
    }
    fprintf(listFile,"%s\n",current_output_fname);
    delete [] current_output_fname;
  } else {
    assert(writeFile); // if not, what have we been writing to?
  }
  if (fclose(writeFile)) {
    error("ERROR closing output file\n");
  }
  writeFile = NULL;
  currSegment += 1;
  currFrame = 0;
  currFeature = 0;
}

