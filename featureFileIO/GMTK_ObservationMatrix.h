/*
 * GMTK_ObservationMatrix: buffer for input data
 *
 * Written by Katrin Kirchhoff <katrin@ee.washington.edu> and Jeff Bilmes
 * <bilmes@ee.washington.edu> 
 *
 * Copyright (c) 2001, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about the suitability of this software
 * for any purpose. It is provided "as is" without express or implied warranty.
 *
 */

#ifndef GMTK_OBSERVATIONMATRIX_H
#define GMTK_OBSERVATIONMATRIX_H

#include "logp.h"
#include "sArray.h"
#include "machine-dependent.h"
#include "pfile.h"
#include "bp_range.h"
#include "GMTK_Stream.h"

#define MAXBUFSIZE 1000 


/* ObservationMatrix: basic data structure for input feature buffer 
 * contains one or more input streams, where each stream is a list of filenames 
 * or a pfile. Must be created using constructor and initialized using
 * 'openFiles'. Reads specified number of features ( continuous and/or discrete)
 * from each input stream and creates one global feature buffer, where all
 * continuous features come first (per frame), followed by all discrete features */



class ObservationMatrix
{


  int _numStreams; // number of input streams 

  unsigned _totalContinuous;  // total num of continuous features
  unsigned _totalDiscrete;    // total num of discrete features 

  unsigned _maxContinuous;  // max number of cont. features in any stream
  unsigned _maxDiscrete;    // max number of disc. features in any stream

  // total number of segments in input streams (identical for all streams)

  unsigned _numSegments; 

 // temporary feature buffers for single frame 

  sArray<float> _contFea; 
  sArray<Int32> _discFea;

  Data32 *_cont_p; // pointers into continuous/discrete feature block
  Data32 *_disc_p; 

  StreamInfo **_inStreams; // input streams

  size_t _bufSize; // maximum number of frames in buffer

  // read pfile features
  bool readPFloats(InFtrLabStream_PFile *, BP_Range *);	
  bool readPInts(InFtrLabStream_PFile *, BP_Range *);

  // read binary features
  bool readBinFloats(unsigned, FILE *, BP_Range *, bool);
  bool readBinInts(unsigned,FILE *,BP_Range *, bool);
  
  // read ascii features
  bool readAscFloats(unsigned, FILE *,BP_Range *);
  bool readAscInts(unsigned, FILE *, BP_Range *);


  void reset();	 // resets pointers to beginning of obs matrix

  void resize(size_t); // resize obs matrix


  // get pointer to individual features in current frame

  float *getContFea(unsigned short n);

  Int32 *getDiscFea(unsigned short n);

  // read single frame into observation matrix

  void readFrame(size_t);	

  // file opening routines

  size_t openBinaryFile(StreamInfo *,size_t);
  size_t openAsciiFile(StreamInfo *,size_t);
  size_t openHTKFile(StreamInfo *,size_t);
  size_t openPFile(StreamInfo *,size_t);

  void closeDataFiles();

 public:

  // the segment number
  size_t segmentNumber;

  // the number of frames in this segment
  unsigned numFrames;

  // number of continuous features in this segment
  unsigned numContinuous;

  // number of discrete features
  unsigned numDiscrete;

  // sum of the above two
  unsigned numFeatures;

  // stride, which is the number of Data32's between
  // each frame, as it might be different than numFeatures

  unsigned stride;

  sArray< Data32 > features; // matrix of features

  ObservationMatrix();
  ~ObservationMatrix();

  // skip to beginning of frame f	 
  
  Data32*const baseAtFrame(unsigned f) {
    assert (f >= 0 && f < numFrames);
    return features.ptr + stride*f;
  }

  float*const floatVecAtFrame(unsigned f) {
    assert (f >= 0 && f < numFrames);
    return (float*)(features.ptr + stride*f);
  }
  
  float*const floatAtFrame(unsigned f) {
    assert (f >= 0 && f < numFrames);
    return (float*)(features.ptr + stride*f);
  }

  float*const floatVecAtFrame(unsigned f, 
			      const unsigned startFeature,
			      const unsigned len) {
    assert (f >= 0 && f < numFrames);
    assert (startFeature >= 0 && startFeature + len <= numContinuous);
    return (float*)(features.ptr + stride*f + startFeature);
  }


  float*const floatVecAtFrame(unsigned f, 
			      const unsigned startFeature) {
    assert (f >= 0 && f < numFrames);
    return (float*)(features.ptr + stride*f + startFeature);
  }


  unsigned*const unsignedAtFrame(unsigned f) {
    assert (f >= 0 && f < numFrames);
    return (unsigned*)(features.ptr + stride*f + numContinuous);
  }

  unsigned& unsignedAtFrame(const unsigned frame, const unsigned feature) {
    assert (frame >= 0 && frame < numFrames);
    assert (feature >= numContinuous
	    &&
	    feature < numFeatures);
    return *(unsigned*)(features.ptr+stride*frame+feature);
  }

  bool elementIsDiscrete(unsigned el) {
    return (el >= numContinuous && el < numFeatures);
  }

  bool elementIsContinuous(unsigned el) {
    return (el >= 0 && el < numContinuous);
  }

  // initialize input streams and allocate obs matrix
  
  void openFiles(int n_files,  
		 const char **fof_names,
		 const char **cont_range_str,
		 const char **disc_range_str,
		 unsigned *n_floats,
		 unsigned *n_ints,
		 unsigned *formats,
		 bool *swapflags);

  // load data for single segment (utterance)

  void loadSegment(const unsigned seg);

  unsigned numSegments() { return _numSegments ; }


  bool active() { return (_numStreams > 0); }

  void printSegmentInfo();

  // print frame

  void printFrame(FILE *, size_t no);


};


////////////////////////////////////////////////
// The global matrix object, must be
// actually defined near where main() is defined.

extern ObservationMatrix globalObservationMatrix;

#endif


