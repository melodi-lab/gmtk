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
#include "GMTK_FileDescription.h"


/* ObservationMatrix: basic data structure for input feature buffer */


class ObservationMatrix
{


  // the segment number
  size_t _segmentNumber;

  // the number of frames in this segment
  unsigned _numFrames;

  // number of continuous features in this segment
  unsigned _numContinuous;

  // number of discrete features
  unsigned _numDiscrete;

  // sum of the above two
  unsigned _numFeatures;

  // stride, which is the number of Data32's between
  // each frame, as it might be different than numFeatures

  unsigned _stride;

  // max number of discrete/continuous features in input streams

  unsigned _maxContinuous;

  unsigned _maxDiscrete;

  Data32 *_cont_p; // pointers into continuous or discrete features
  Data32 *_disc_p; 

 // temporary feature buffers for single frame

  sArray<float> cont_fea; 
  sArray<Int32> disc_fea;

  size_t _bufSize; // maximum number of frames in buffer

  // read pfile features
  bool readPFloats(InFtrLabStream_PFile *, BP_Range *);	

  bool readPInts(InFtrLabStream_PFile *, BP_Range *);

  // read binary features

  bool readBinFloats(unsigned, FILE *, BP_Range *, bool);

  bool readCRCFloats(unsigned, FILE *, BP_Range *, bool);

  bool readBinInts(unsigned,FILE *,BP_Range *, bool);
  
  // read ascii features

  bool readAscFloats(unsigned, FILE *,BP_Range *);

  bool readAscInts(unsigned, FILE *, BP_Range *);

  // advance pointers for reading in next frame

  void next();

public:

  sArray< Data32 > features; // matrix of features

  ObservationMatrix(size_t,unsigned,unsigned,unsigned,unsigned);
  ~ObservationMatrix();

  unsigned getSegmentNumber() { return _segmentNumber ;}
  unsigned getNumFrames() { return _numFrames; }
  unsigned getNumContinuous() { return _numContinuous; }
  unsigned getNumDiscrete() { return _numDiscrete; }
  unsigned getNumFeatures() { return _numFeatures; }
  unsigned getStride() {return _stride; }
  size_t getBufSize() { return _bufSize; }

  void setSegNo(size_t n) { _segmentNumber = n; }

  void setNumFrames(size_t n) { _numFrames = n; }

  void reset();	 // resets pointers to beginning of obs matrix, 
                 // call after every utterance

  void resize(size_t); // resize obs matrix

  // set internal pointers to beginning of frame f

  void gotoFrame(size_t f);

  // get pointer to individual features in current frame

  float *getContFea(unsigned short n);

  Int32 *getDiscFea(unsigned short n);

  // skip to beginning of frame f	 
  
  Data32*const baseAtFrame(unsigned f) {
    assert (f >= 0 && f < _numFrames);
    return features.ptr + _stride*f;
  }

  float*const floatVecAtFrame(unsigned f) {
    assert (f >= 0 && f < _numFrames);
    return (float*)(features.ptr + _stride*f);
  }
  
  float*const floatAtFrame(unsigned f) {
    assert (f >= 0 && f < _numFrames);
    return (float*)(features.ptr + _stride*f);
  }

  float*const floatVecAtFrame(unsigned f, 
			      const unsigned startFeature,
			      const unsigned len) {
    assert (f >= 0 && f < _numFrames);
    assert (startFeature >= 0 && startFeature + len <= _numContinuous);
    return (float*)(features.ptr + _stride*f + startFeature);
  }


  float*const floatVecAtFrame(unsigned f, 
			      const unsigned startFeature) {
    assert (f >= 0 && f < _numFrames);
    return (float*)(features.ptr + _stride*f + startFeature);
  }


  unsigned*const unsignedAtFrame(unsigned f) {
    assert (f >= 0 && f < _numFrames);
    return (unsigned*)(features.ptr + _stride*f + _numContinuous);
  }

  unsigned& unsignedAtFrame(const unsigned frame, const unsigned feature) {
    assert (frame >= 0 && frame < _numFrames);
    assert (feature >= _numContinuous
	    &&
	    feature < _numFeatures);
    return *(unsigned*)(features.ptr+stride*frame+feature);
  }

  bool elementIsDiscrete(unsigned el) {
    return (el >= _numContinuous && el < _numFeatures);
  }

  bool elementIsContinuous(unsigned el) {
    return (el >= 0 && el < _numContinuous);
  }

  // read single frame into observation matrix

  void readFrame(size_t,FileDescription**,int);	

  // print frame

  void printFrame(FILE *stream, size_t no);


};


////////////////////////////////////////////////
// The global matrix object, must be
// actually defined near where main() is defined.

extern ObservationMatrix globalObservationMatrix;

#endif


