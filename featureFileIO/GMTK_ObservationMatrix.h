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
#include "bp_range.h"
#include "GMTK_Stream.h"




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
  bool readPFrame(InFtrLabStream_PFile *, BP_Range *, BP_Range *);

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


  // the segment number
  size_t _segmentNumber;

  // the number of frames in this segment
  unsigned _numFrames;

  // the number of frames in this segment that are not skipped
  unsigned _numNonSkippedFrames;

  // number of continuous features in this segment
  unsigned _numContinuous;

  // number of discrete features
  unsigned _numDiscrete;

  // sum of the above two
  unsigned _numFeatures;

  // stride, which is the number of Data32's between
  // each frame, as it might be different than numFeatures

  unsigned _stride;

  // number of frames to skip at the beginning
  unsigned _startSkip;

  // number of frames to skip a the end.
  unsigned _endSkip;

  /////////////////////////////////////////////
  unsigned _totalSkip;

 public:
  
  /////////////////////////////////////////////////
  // constructor just makes an inactive object.
  ObservationMatrix();
  ~ObservationMatrix();

  /////////////////////////////////////////////////////////
  // the true constructor, in that it initializes the
  // input streams and allocate obs matrix.
  void openFiles(int n_files,  
		 const char **fof_names,
		 const char **cont_range_str,
		 const char **disc_range_str,
		 unsigned *n_floats,
		 unsigned *n_ints,
		 unsigned *formats,
		 bool *swapflags,
		 const unsigned _startSkip = 0,
		 const unsigned _endSkip = 0);

  
  ///////////////////////////////////////////////////////////
  // returns true if the current observation matrix
  // is "active" in the sence that there are open
  // files, and data can be read from them. 
  bool active() { return (_numStreams > 0); }

  // the segment number
  size_t segmentNumber() { return _segmentNumber; }

  // the number of frames not including the skip in this segment
  unsigned numFrames() { return _numFrames ; }

  // the number of "real" frames in this segment
  unsigned numNonSkippedFrames() { return _numNonSkippedFrames; }

  // number of continuous features in this segment
  unsigned numContinuous() { return _numContinuous; }

  // number of discrete features
  unsigned numDiscrete() { return _numDiscrete; }

  // sum of the above two
  unsigned numFeatures() { return _numFeatures; }

  // stride, which is the number of Data32's between
  // each frame, as it might be different than numFeatures

  unsigned stride() { return _stride; }

  // number of frames to skip at the beginning
  unsigned startSkip() { return _startSkip; }

  // number of frames to skip a the end.
  unsigned endSkip() { return _endSkip; }

  unsigned numSegments() { return _numSegments ; }

  // The actual matrix of features, which may be used directly
  // if so desired. This is a matrix of Data32s, some of which
  // might refer to single precision floating point numbers,
  // and some of which might refer to 32 bit unsigned integers.
  // The access routines below will index into this for convenient
  // user access. 
  sArray< Data32 > features; 

  /////////////////////////////////////////
  // A pointer to the starting base of: 
  Data32 *featuresBase;


  /////////////////////////////////////////////
  // these access routines respect the start frame and end frame.
  Data32*const baseAtFrame(unsigned f) {
    assert (f >= 0 && f < _numFrames);
    return featuresBase + _stride*f;
  }

  float*const floatVecAtFrame(unsigned f) {
    assert (f >= 0 && f < _numFrames);
    return (float*)(featuresBase + _stride*f);
  }
  
  float*const floatAtFrame(unsigned f) {
    assert (f >= 0 && f < _numFrames);
    return (float*)(featuresBase + _stride*f);
  }

  float*const floatVecAtFrame(unsigned f, 
			      const unsigned startFeature,
			      const unsigned len) {
    assert (f >= 0 && f < _numFrames);
    assert (startFeature >= 0 && startFeature + len <= _numContinuous);
    return (float*)(featuresBase + _stride*f + startFeature);
  }


  float*const floatVecAtFrame(unsigned f, 
			      const unsigned startFeature) {
    assert (f >= 0 && f < _numFrames);
    return (float*)(featuresBase + _stride*f + startFeature);
  }


  unsigned*const unsignedAtFrame(unsigned f) {
    assert (f >= 0 && f < _numFrames);
    return (unsigned*)(featuresBase + _stride*f + _numContinuous);
  }

  unsigned& unsignedAtFrame(const unsigned frame, const unsigned feature) {
    assert (frame >= 0 && frame < _numFrames);
    assert (feature >= _numContinuous
	    &&
	    feature < _numFeatures);
    return *(unsigned*)(featuresBase+_stride*frame+feature);
  }

  bool elementIsDiscrete(unsigned el) {
    return (el >= _numContinuous && el < _numFeatures);
  }

  bool elementIsContinuous(unsigned el) {
    return (el >= 0 && el < _numContinuous);
  }


  // load data for single segment (utterance)

  void loadSegment(const unsigned seg);
  void storeSegment(const unsigned seg) { error("not implemented\n"); }

  void printSegmentInfo();

  // print frame

  void printFrame(FILE *, size_t no);


};


////////////////////////////////////////////////
// The global matrix object, must be
// actually defined near where main() is defined.

extern ObservationMatrix globalObservationMatrix;

#endif


