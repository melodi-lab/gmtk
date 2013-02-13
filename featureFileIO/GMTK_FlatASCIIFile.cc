
/*
 * GMTK_ASCIIFile.cc
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

#include "file_utils.h"
#include "GMTK_FlatASCIIFile.h"

using namespace std;

FlatASCIIFile::FlatASCIIFile(const char *name, unsigned nfloats, unsigned nints, 
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
    cppCommandOptions(cppCommandOptions),
    close(true)
{
  buffer = NULL;
  fileName = name;
  writeFile = NULL;
  if (name == NULL) 	
    error("FlatASCIIFile: File name is NULL for stream %i",num);	
  if (nfloats == 0 && nints == 0)
    error("FlatASCIIFile: number of float and int features cannot both be zero");

  _numContinuousFeatures = nfloats; 
  _numDiscreteFeatures   = nints;
  _numFeatures           = nfloats + nints;

  FILE *f = openCPPableFile(name, cppIfAscii, cppCommandOptions);
  if (!f)
    error("FlatASCIIFile: couldn't open '%s' for reading", name);

  unsigned lineNum = 1;
  unsigned totalFrames = 0;
  int prevSegment = -1, currSegment;
  int prevFrame=-1, currFrame;
  int tmp;
  nSegments = 0;

  // consume CPP special directives if any
#ifdef PIPE_ASCII_FILES_THROUGH_CPP
  if(cppIfAscii) {
    while((tmp=fgetc(f))==CPP_DIRECTIVE_CHAR) {
      while((tmp=fgetc(f))!='\n');
      lineNum++;
    }
    ungetc(tmp,f);
  }
#endif

  // read entire file once to count segments, frames per segment
  while (!feof(f)) {
    // read a frame
    int nread = fscanf(f, "%d %d", &currSegment, &currFrame);
    if (nread == EOF) break;
    if (nread != 2) {
      error("ERROR: FlatASCIIFile: error reading '%s'", name);
    }
    if (currSegment != prevSegment) {
      if (currSegment != prevSegment+1)
	error("ERROR: FlatASCIIFile: expected segment %d, but got %d at line %u in '%s'",
	      prevSegment+1, currSegment, lineNum, name);
      if (currFrame != 0) 
	error("ERROR: FlatASCIIFile: expected frame 0, but got %d at line %u in '%s'",
	      currFrame, lineNum, name);
      if (prevSegment > -1) {
	nFrames.push_back(prevFrame+1); // count of frames in previous segment
	nSegments += 1;
      }
      prevSegment = currSegment;
    } else {
      if (currFrame != prevFrame+1) 
	error("ERROR: FlatASCIIFile: expected frame %d, but got %d at line %u in '%s'",
	      prevFrame+1, currFrame, lineNum, name);
    }
    prevFrame = currFrame;
    // skip rest of line
    while((tmp=fgetc(f)) != '\n')
      ;
    lineNum += 1;
    totalFrames += 1;
  }
  nSegments += 1;
  nFrames.push_back(prevFrame+1);
  closeCPPableFile(f, cppIfAscii);

  buffer = new Data32[totalFrames * _numFeatures];
  if (!buffer) 
    error("ERROR: FlatASCIIFile: unable to allocate memory for %u frames", totalFrames);
  segment = new Data32 *[nSegments];
  assert(segment);

  // read file again, loading data into memory

  f = openCPPableFile(name, cppIfAscii, cppCommandOptions);
  if (!f)
    error("FlatASCIIFile: couldn't open '%s' for reading", name);
  
  lineNum = 1;
  // consume CPP special directives if any
#ifdef PIPE_ASCII_FILES_THROUGH_CPP
  if(cppIfAscii) {
    while((tmp=fgetc(f))==CPP_DIRECTIVE_CHAR) {
      while((tmp=fgetc(f))!='\n');
      lineNum++;
    }
    ungetc(tmp,f);
  }
#endif

  float *fDest;
  Int32 *iDest = (Int32 *)buffer;
  for (unsigned seg=0; seg < nSegments; seg+=1) {
    fDest = (float *) iDest;
    segment[seg] = (Data32 *) fDest;
    //printf("seg %u  has %u frames\n",seg,nFrames[seg]);
    for (unsigned frame=0; frame < nFrames[seg]; frame+=1) {
      if (fscanf(f, "%d %d", &currSegment, &currFrame) != 2) {
	error("ERROR: FlatASCIIFile: error reading '%s'", name);
      }
      assert(currSegment == (int)seg);
      assert(currFrame == (int)frame);
      for (unsigned n = 0; n < _numContinuousFeatures; n+=1) {
	if (fscanf(f,"%e", fDest++) != 1) {
	  error("ERROR: FlatASCIIFile: couldn't read %u'th float in segment %u, frame %u",
		n, seg, frame);
	}
      }
      iDest = (Int32 *)fDest;
      for (unsigned n = 0; n < _numDiscreteFeatures; n+=1) {
	if (fscanf(f,"%d", iDest++) != 1) {
	error("ERROR: FlatASCIIFile: couldn't read %u'th int in segment %u, frame %u",
	      n, seg, frame);
	}
      }
      lineNum += 1;
    }
  }
  closeCPPableFile(f, cppIfAscii);

  if (contFeatureRangeStr) {
    contFeatureRange = new Range(contFeatureRangeStr_,0,_numContinuousFeatures);
    assert(contFeatureRange);
    _numLogicalContinuousFeatures = contFeatureRange->length();
  } else
    _numLogicalContinuousFeatures = nfloats;

  if (discFeatureRangeStr) {
    discFeatureRange = new Range(discFeatureRangeStr_,0,_numDiscreteFeatures);
    assert(discFeatureRange);
    _numLogicalDiscreteFeatures = discFeatureRange->length();
  } else
    _numLogicalDiscreteFeatures = nints;
  
  _numLogicalFeatures = _numLogicalContinuousFeatures + _numLogicalDiscreteFeatures;

  if (segRangeStr_)
    segRange = new Range(segRangeStr_,0,nSegments);
}


FlatASCIIFile::FlatASCIIFile(const char *name, unsigned nfloats, unsigned nints) 
  : ObservationFile(NULL, NULL, NULL, NULL),
    cppIfAscii(false), cppCommandOptions(NULL)
{
  _numContinuousFeatures = _numLogicalContinuousFeatures = nfloats;
  _numDiscreteFeatures   = _numLogicalDiscreteFeatures   = nints;
  _numFeatures           = _numLogicalFeatures           = nfloats + nints;

  logicalObservationBuffer = NULL;
  logicalObsBufSize = 0;
  buffer = NULL;
  segment = NULL;

  nSegments = 0;
  currSegment = 0;
  currFrame = 0;
  currFeature = 0;

  if (strcmp(name,"-")) {
    writeFile = fopen(name, "w+");
    if (!writeFile) {
      error("ERROR: failed to open output file '%s'\n", name);
    }
    close = true;
  } else {
    writeFile = stdout;
    close = false;
  }
  if (!writeFile) {
    error("FlatASCIIFile: couldn't open '%s' for reading\n", name);
  }
}


void 
FlatASCIIFile::writeSegment(Data32 const *segment, unsigned nFrames) {
  assert(currFeature == 0);
  for (unsigned f=0; f < nFrames; f+=1) {
    fprintf(writeFile, "%u %u", currSegment, f);
    for (unsigned i=0; i < _numContinuousFeatures; i+=1) {
      fprintf(writeFile, " %.8e", *((float *)segment++));
    }
    for (unsigned i=0; i < _numDiscreteFeatures; i+=1) {
      fprintf(writeFile, " %u", *((unsigned *)segment++));
    }
    fprintf(writeFile, "\n");
  }
  currSegment += 1;
  nSegments += 1;
  currFrame = 0;
  currFeature = 0;
}


void 
FlatASCIIFile::writeFrame(Data32 const *frame) {
  assert(currFeature == 0);
  fprintf(writeFile, "%u %u", currSegment, currFrame++);
  for (unsigned i=0; i < _numContinuousFeatures; i+=1) {
    fprintf(writeFile, " %.8e", *((float *)frame++));
  }
  for (unsigned i=0; i < _numDiscreteFeatures; i+=1) {
    fprintf(writeFile, " %u", *((unsigned *)frame++));
  }
  fprintf(writeFile, "\n");
  currFeature = 0;
}


void 
FlatASCIIFile::writeFeature(Data32 x) {
  if (currFeature == 0) {
    fprintf(writeFile, "%u %u", currSegment, currFrame);
  }
  if (currFeature < _numContinuousFeatures) {
    fprintf(writeFile, " %.8e", *(float *)&x);
  } else {
    fprintf(writeFile, " %u", *(unsigned *)&x);
  }
  currFeature += 1;
  if (currFeature == _numFeatures) {
    fprintf(writeFile, "\n");
    currFrame += 1;
    currFeature = 0;
  }
}


void 
FlatASCIIFile::endOfSegment() {
  assert(currFeature == 0);
  currSegment += 1;
  currFrame = 0;
  currFeature = 0;
  nSegments += 1;
}

