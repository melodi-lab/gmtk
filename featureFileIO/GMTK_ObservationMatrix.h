/*
 * GMTK_ObservationMatrix: buffer for input data
 *
 * Written by Katrin Kirchhoff <katrin@ee.washington.edu> and Jeff Bilmes
 * <bilmes@ee.washington.edu> 
 *
 *  Modified by Karim Filali <karim@cs.washington.edu> on 08sep2003:
 *
 *   - ASCII files are piped through CPP
 *
 *   - Changed the way data is read into memory.  It used to be read
 *   frame by frame, iterating over all streams. Now it is read as a
 *   whole sentence for each stream, then the sentences are glued
 *   together in the final buffer.  This should be more efficient
 *   since disk reads are consolidated.  TODO: maybe a add a pre-load
 *   function that reads a bunch of sentences in memory upfront to
 *   avoid the need to go through the whole loading process over and
 *   over again for that subset of sentences.
 *
 *   - A per-sentence range can be applied to each stream individually
 *
 *   - The length of sentences across streams can be adjusted so that
 *   they match.  The following options are supported when there is a length mismatch:
 *
 *       a) Report an error
 *
 *       b) Repeat the last frame up to the maximum number of frames
 *
 *       c) Repeat the first frame
 *
 *       d) Repeat all frames "segmentally" i.e. in a segmental k-means initialization fashion
 *
 *       e) Truncate from the end
 *
 *       f) Truncate from the start
 *
 *       g) TODO: repeat a range of frames, for example , the last three, all frames...
 *
 *   - The number of sentences per stream can be adjusted so that they
 *   match across streams.  The following options are available:
 *
 *       a) Report an error 
 *
 *       b) Truncate 
 *
 *       c) Repeat the last sentence
 *
 *       d) Wrap around:  repeat setences from the beginning after the end is reached. 
 *
 *   - The following transformations can be applied to the data:
 *
 *        a) Mean and variance normalization
 *       
 *        b) Multiplication by some constant
 *
 *        c) Upsampling (either by simple repetition of frames or interpolation)
 *
 * $Header$
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
#include "GMTK_Stream.h"
#include "ieeeFPsetup.h"

#define INITIAL_NUM_FRAMES_IN_BUFFER 1000

#define NONE_LETTER 'X'
#define TRANS_NORMALIZATION_LETTER 'N'
#define TRANS_MEAN_SUB_LETTER 'E'
#define TRANS_UPSAMPLING_LETTER 'U'
#define TRANS_HOLD_LETTER 'H'
#define TRANS_SMOOTH_LETTER 'S'
#define TRANS_MULTIPLICATION_LETTER 'M'
#define TRANS_ARMA_LETTER 'R'
#define TRANS_OFFSET_LETTER 'O'
#define FILTER_LETTER 'F'
#define END_STR (-1)
#define UNRECOGNIZED_TRANSFORM (-2)
#define SEPARATOR "_"

#define MAX_NUM_STORED_FILTERS 10
#define MAX_FILTER_LEN 100
#define MIN_FILTER_LEN 1

#define ALLOW_VARIABLE_DIM_COMBINED_STREAMS 1

#ifdef DEBUG
#define DBGFPRINTF(_x_) fprintf _x_
#else
#define DBGFPRINTF(_x_)
#endif


// transformations
enum {
  NONE,
  UPSAMPLE_HOLD,
  UPSAMPLE_SMOOTH,
  DOWNSAMPLE,
  DELTAS,
  DOUBLE_DELTAS,
  MULTIPLY,
  NORMALIZE,
  MEAN_SUB,
  ARMA,
  FILTER,
  OFFSET
};

enum {
  FRAMEMATCH_ERROR,
  FRAMEMATCH_REPEAT_LAST,
  FRAMEMATCH_REPEAT_FIRST,
  FRAMEMATCH_EXPAND_SEGMENTALLY,
  FRAMEMATCH_TRUNCATE_FROM_END,
  FRAMEMATCH_TRUNCATE_FROM_START,
  FRAMEMATCH_REPEAT_RANGE_AT_START, // not used yet
  FRAMEMATCH_REPEAT_RANGE_AT_END   // not used yet
};

enum {
  SEGMATCH_ERROR,
  SEGMATCH_TRUNCATE_FROM_END,
  SEGMATCH_REPEAT_LAST,
  SEGMATCH_WRAP_AROUND
};

// ftrcombo operation flags
enum {
  FTROP_NONE = 0,
  FTROP_ADD,
  FTROP_SUB,
  FTROP_MUL,
  FTROP_DIV
  };

/* ObservationMatrix: basic data structure for input feature buffer 
 * contains one or more input streams, where each stream is a list of filenames 
 * or a pfile. Must be created using constructor and initialized using
 * 'openFiles'. Reads specified number of features ( continuous and/or discrete)
 * from each input stream and creates one global feature buffer, where all
 * continuous features come first (per frame), followed by all discrete features */



class ObservationMatrix {


 // Actions when the number of sentences is different across streams
 

  ////////////  Global info independent of the segment ///////////////
  unsigned      _numStreams;           // number of input streams 
  unsigned      _numContinuous;        // number of continuous features (summed over all streams)
  unsigned      _numDiscrete;          // number of discrete feature    (summed over all streams)
  unsigned      _numFeatures;          // sum of the above two (_numContinuous and _numDiscrete)
  // stride, which is the number of Data32's between
  // each frame, as it might be different than numFeatures
  unsigned      _stride;
  unsigned      _maxContinuous;        // max number of cont. features in any stream
  unsigned      _maxDiscrete;          // max number of disc. features in any stream
  unsigned      _numSegments;          // total number of segments in input streams (identical for all streams)
  size_t        _bufSize;              // maximum number of frames in buffer;  is dynamically increased when needed
  unsigned      _startSkip;            // number of frames to skip at the beginning
  unsigned      _endSkip;              // number of frames to skip a the end.
  unsigned      _totalSkip;            // _startSkip + _endSkip
  StreamInfo ** _inStreams;            // input streams
  //////////////////////////////////////////////////////////////

  ////////////  Current segment info  //////////////////////////
  size_t        _segmentNumber;        // the segment number
  unsigned      _numFrames;            // the number of frames in this segment
  unsigned      _numNonSkippedFrames;  // the number of frames in this segment that are not skipped
  ///////////////////////////////////////////////////////////////

  ////////////  Data structures for applying various transformations to the data //////
  //                              karim - 29aug2003
  bool          _cppIfAscii;
  char*         _cppCommandOptions;
  unsigned*     _actionIfDiffNumFrames;
  unsigned*     _actionIfDiffNumSents;
  const char**  _prrngStr;
  const char**  _preTransFrameRangeStr;
  const char*   _finalFrameRangeStr;
  char**        _perStreamPreTransforms;
  char*         _postTransforms;
  unsigned      _ftrcombo;
  char*         _filterFileName;
  float         _filterCoeffs[MAX_FILTER_LEN];
  float         _preStoredFilterCoeffs[MAX_NUM_STORED_FILTERS][MAX_FILTER_LEN];
  unsigned      _startFilterSkip,_endFilterSkip;
  sArray<int>   _repeat;               // array used when adjusting the length of sentences
  /////////////////////////////////////////////////////////////////////////////////////////


  ////////////  Temporary buffers for a sentence  //////////////
  sArray<float> _tmpFloatSenBuffer;
  sArray<Int32> _tmpIntSenBuffer;
  //////////////////////////////////////////////////////////////


  // read a whole sentence into memory
  bool readBinSentence(float* float_buffer, unsigned n_floats, Int32* int_buffer, unsigned n_ints,StreamInfo* s);
  bool readPfileSentence(const unsigned segno, float* float_buffer,  Int32* int_buffer, InFtrLabStream_PFile *f);
  bool readAsciiSentence(float* float_buffer, unsigned num_floats, Int32* int_buffer, unsigned num_ints,unsigned n_samples, FILE *f);
  

  void copyToFinalBuffer(unsigned stream_no,float* ,Int32*,Range* ,Range*,Range*);
  void copyAndAdjustLengthToFinalBuffer(unsigned stream_no,float* float_buf,Int32* int_buf,Range* float_rng,Range* int_rng,Range* pr_rng,unsigned prrng_n_samps);
  
  
  /////////////////         data transformation routines      //////////////////  
  
  template<class T> unsigned applyPreTransformFrameRange(sArray<T>* tmp_sen_buffer, unsigned vec_size, unsigned stride, unsigned num_frames, const char* preTransFrameRangeStr);

  unsigned applyFinalFrameRange(unsigned num_frames,const char* final_frame_range_str);

  int  parseTransform(char*& trans_str, int& magic_int, double& magic_double);
  void applyTransforms(char* trans_str, unsigned num_floats, unsigned num_ints, unsigned num_frames);
  void applyPostTransforms(char* trans_str, unsigned num_floats, unsigned num_ints, unsigned num_frames);
  
  void inPlaceMeanSubVarNorm(float* x, int vec_size, int stride, int num_frames);
  void inPlaceMeanSub(float* x, int vec_size, int stride, int num_frames);
  void multiply(float* x, unsigned vec_size, unsigned stride, unsigned num_frames, double multiplier);
  void addOffset(float* x, unsigned vec_size, unsigned stride, unsigned num_frames, double offset);
  void arma(float* x, unsigned vec_size, unsigned stride, unsigned num_frames,unsigned order);
  template<class T> void upsampleHold(sArray<T>* x, unsigned vec_size, unsigned stride, unsigned num_frames,unsigned upsample);
  template<class T> void upsampleSmooth(sArray<T>* tmp_sen_buffer, unsigned vec_size, unsigned stride, unsigned num_frames, unsigned upsample);
  void upsampleIntNoSmooth(sArray<int>* tmp_sen_buffer, unsigned vec_size, unsigned stride, unsigned num_frames, unsigned upsample);

  void upsampleHold(unsigned num_frames, unsigned upsample);
  void upsampleSmooth(unsigned num_frames, unsigned upsample);

  void filter(float* x, unsigned vec_size, unsigned stride, unsigned num_frames,float* filter_coeffs, unsigned filter_len);
  
  void readFilterFromFile(char* filter_file_name,float* filter_coeffs, unsigned& filter_len);
  
  ///////////////////////////////////////////////////////////////////////
  
 
  void   reset();	       // resets pointers to beginning of obs matrix
  void   resize(size_t);  // resize obs matrix
  
  // get pointer to individual features in current frame
  float *getContFea(unsigned short n);
  Int32 *getDiscFea(unsigned short n);
  
  ///////////////        file opening routines    //////////////////
  
  size_t openBinaryFile(StreamInfo *,size_t);
  size_t openAsciiFile (StreamInfo *,size_t);
  size_t openHTKFile   (StreamInfo *,size_t);
  size_t openPFile     (StreamInfo *,size_t);
  
  void   closeDataFiles();
  
  /////////////////////////////////////////////////////////////////
  

   // Auxilliary functions
  unsigned checkNumSegments(StreamInfo* streams[], unsigned n_streams, unsigned* action_if_diff_seg_len);
  void     checkNumFeatures(StreamInfo* streams[], unsigned n_streams, unsigned* n_floats, unsigned* n_ints, unsigned* max_n_floats, unsigned* max_n_ints, unsigned ftrcombo);
  

 public:
  
  /////////////////////////////////////////////////
  // constructor just makes an inactive object.
  ObservationMatrix();
  ~ObservationMatrix();

  /////////////////////////////////////////////////////////
  // the true constructor, in that it initializes the
  // input streams and allocate obs matrix.
  void openFiles(int            n_files,  
		 const char **  fof_names,
		 const char **  cont_range_str,
		 const char **  disc_range_str,
		 unsigned *     n_floats,
		 unsigned *     n_ints,
		 unsigned *     formats,
		 bool *         swapflags,
		 const unsigned _startSkip = 0,
		 const unsigned _endSkip   = 0,
		 bool           cppIfAscii = true,
		 char *         cppCommandOptions      = NULL,
		 const char **  pr_range_str           = NULL,
		 unsigned*      actionIfDiffNumFrames  = NULL,
		 unsigned*      actionIfDiffNumSents   = NULL,
		 char**         perStreamPreTransforms = NULL,
		 char*          postTransforms         = NULL,
		 unsigned       ftrcombo               = FTROP_NONE,
		 const char**   sr_range_str           = NULL,
		 const char**   preTransFrameRangeStr  = NULL,
		 const char*    finalFrameRangeStr     = NULL
);

  /////////////////////////////////////////////////////////////////////////
  // A version of the above that just opens one file.
  void openFile(const char *   f_name,
		const char *   cont_range_str,
		const char *   disc_range_str,
		unsigned       n_floats,
		unsigned       n_ints,
		unsigned       formats,
		bool           swapflags,
		const unsigned _startSkip = 0,
		const unsigned _endSkip   = 0,
		bool           cppIfAscii = true,
		char *         cppCommandOptions      = NULL,
		const char **  pr_range_str           = NULL,
		unsigned*      actionIfDiffNumFrames  = NULL,
		unsigned*      actionIfDiffNumSents   = NULL,
		char**         perStreamPreTransforms = NULL,
		char*          postTransforms         = NULL,
		unsigned       ftrcombo               = FTROP_NONE,
		const char**   sr_range_str           = NULL,
		const char**   preTransFramerangeStr  = NULL,
		const char*    finalFrameRangeStr     = NULL
		) {
    openFiles(1,&f_name,&cont_range_str,&disc_range_str,&n_floats,&n_ints,&formats,&swapflags,
	      _startSkip,_endSkip,cppIfAscii,cppCommandOptions,pr_range_str,
	      actionIfDiffNumFrames,actionIfDiffNumSents,
	      perStreamPreTransforms,postTransforms,ftrcombo,sr_range_str,preTransFramerangeStr,finalFrameRangeStr);
  }

  unsigned formatStrToNumber(const char * fmt) {
    if (strcmp(fmt,"htk") == 0)
      return HTK;
    else if (strcmp(fmt,"binary") == 0)
      return RAWBIN;
    else if (strcmp(fmt,"ascii") == 0)
      return RAWASC;
    else if (strcmp(fmt,"pfile") == 0)
      return  PFILE;
    else {
      error("ERROR: Unknown observation file format type: '%s'\n",fmt);
      return 0;
    }
  }


  
  ///////////////////////////////////////////////////////////
  // returns true if the current observation matrix
  // is "active" in the sense that there are open
  // files, and data can be read from them. 
  bool     active()        { return (_numStreams > 0); }
  unsigned numStreams()    { return _numStreams ;  }
  unsigned numSegments()   { return _numSegments ;  }
  size_t   segmentNumber() { return _segmentNumber; }  // the segment number
  unsigned numFrames()     { return _numFrames ;    }  // the number of frames not including the skip in this segment
  unsigned numFrames(unsigned segno);  // The difference between this method the numFrames() is that numFrames() returns the number of frames in teh CURRENT segment (the one that has been loaded last using loadSegment(...); while numFrames(segno) loads sentence number segno.
  unsigned numNonSkippedFrames() { return _numNonSkippedFrames; }    // the number of "real" frames in this segment
  unsigned numContinuous() { return _numContinuous; }  // number of continuous features in this segment
  unsigned numDiscrete()   { return _numDiscrete;   }  // number of discrete features
  unsigned numFeatures()   { return _numFeatures;   }  // sum of the above two

  // stride, which is the number of Data32's between
  // each frame, as it might be different than numFeatures
  unsigned stride()        { return _stride;        }


  unsigned startSkip()     { return _startSkip;     }  // number of frames to skip at the beginning
  unsigned endSkip()       { return _endSkip;       }  // number of frames to skip a the end.


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

  // advances the featuresBase pointer by n frames
  // Side effect:  decreases the number of frames by n
  void mvFeatBaseByNumFrames(unsigned n) {
    assert (n >= 0 && n < _numFrames);
    featuresBase += (_stride*n);
    _numFrames -= n;
  }

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

  float floatAtFrame(unsigned f, 
		     const unsigned startFeature) {
    assert (f >= 0 && f < _numFrames);
    return *(float*)(featuresBase + _stride*f + startFeature);
  }

  // TODO: change name to unsignedVecAtFrame
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
  /////////////////////////////////////////////////////////////

  bool   checkIfSameNumSamples(unsigned segno, unsigned&,unsigned&);

  void   loadSegment(unsigned seg);   // load data for single segment (utterance)
  void   storeSegment() { error("not implemented\n"); }

  void   printSegmentInfo();
  void   printFrame(FILE *, size_t no);

  float* getPhysicalStartOfFloatFeaturesBuffer() { return (float*) features.ptr;                  }
  Int32* getPhysicalStartOfIntFeaturesBuffer()   { return (Int32*) features.ptr + _numContinuous; }

  StreamInfo* getStream(unsigned stream_no)      { 
    if(stream_no >= _numStreams)
      return NULL;
    else
      return _inStreams[stream_no];
  }

  //////// old stuff -- not used anymore ////////////////////////////
  void   framewiseLoadSegment(const unsigned segno); // old version that
  // used to load the sentence frame by frame, which is both
  // inefficient and inflexible when it comes to manipulationg the
  // frame set
  ///////////////////////////////////////////////////////////////////
};



template<class T>
void ObservationMatrix::upsampleHold(sArray<T>* tmp_sen_buffer, unsigned vec_size, unsigned stride, unsigned num_frames,unsigned upsample) {

  if(vec_size==0 || stride ==0)
    return;

  T* tmp_buf=new T[num_frames*vec_size*(upsample+1)];
  T* x=tmp_sen_buffer->ptr;

  // copy buffer into a temporary one
  for(unsigned i=0;i<num_frames;++i)  
    for(unsigned j=0;j<vec_size;++j) {
	tmp_buf[i*vec_size+j]=x[i*stride+j];
    }

  tmp_sen_buffer->resize(num_frames*stride*(upsample+1)*2);
  resize(num_frames*(upsample+1)*2);
   
  x=tmp_sen_buffer->ptr;

  unsigned cnt=0;
  for(unsigned i=0;i<num_frames;++i)  
    for(unsigned k=0;k<upsample+1;++k) {
      for(unsigned j=0;j<vec_size;++j) {
	x[cnt*stride+j]=tmp_buf[i*vec_size+j];
      }
      cnt++;
    }

  delete [] tmp_buf;

}

template<class T> void ObservationMatrix::upsampleSmooth(sArray<T>* tmp_sen_buffer, unsigned vec_size, unsigned stride, unsigned num_frames, unsigned upsample) {

  if(vec_size==0 || stride ==0)
    return;
  
  T* tmp_buf=new T[num_frames*vec_size*(upsample+1)];
  T* x=tmp_sen_buffer->ptr;

  // copy buffer into a temporary one
  for(unsigned i=0;i<num_frames;++i)  
    for(unsigned j=0;j<vec_size;++j) {
	tmp_buf[i*vec_size+j]=x[i*stride+j];
    }

  unsigned after_transform_num_frames = ((num_frames-1)*(upsample+1)+1)*2;

  tmp_sen_buffer->resize(after_transform_num_frames*stride);
  resize(after_transform_num_frames);
   
  x=tmp_sen_buffer->ptr;


  double inc;
  double additive_const;
  unsigned cnt=0;
  for(unsigned i=0;i<num_frames-1;++i) {
    inc=0;
    for(unsigned k=0;k<upsample+1;++k) {
      for(unsigned j=0;j<vec_size;++j) {
	additive_const =  tmp_buf[(i+1)*vec_size+j] - tmp_buf[i*vec_size+j];
	additive_const /= (upsample+1);
	inc = k*additive_const;
	x[cnt*stride+j]=tmp_buf[i*vec_size+j] + (T)inc;
      }
      cnt++;
    }
  }

    for(unsigned j=0;j<vec_size;++j) {
      x[cnt*stride+j]=tmp_buf[(num_frames-1)*vec_size+j];
    }

    delete [] tmp_buf;

}


/**
 *
 * side effects: - updates the number of frames num_frames
*/

template<class T>
unsigned ObservationMatrix::applyPreTransformFrameRange(sArray<T>* tmp_sen_buffer, unsigned vec_size, unsigned stride, unsigned num_frames, const char* preTransFrameRangeStr) {

  assert(vec_size!=0 && stride != 0);

  DBGFPRINTF((stderr,"In ObservationMatrix::applyPreTransformFrameRange: creating preTransFrameRange preTransFrameRangeStr= %s, num_frames=%d.\n",preTransFrameRangeStr,num_frames));
  Range* preTransFrameRange = new Range(preTransFrameRangeStr==NULL?NULL:preTransFrameRangeStr,0,num_frames);
 assert(preTransFrameRange != NULL);
  DBGFPRINTF((stderr,"In ObservationMatrix::applyPreTransformFrameRange: after creating preTransFrameRange preTransFrameRangeStr= %s, num_frames=%d.\n",preTransFrameRangeStr,num_frames));

  unsigned new_num_frames= preTransFrameRange->length();

  if(preTransFrameRange->full()) return new_num_frames;

  T* tmp_buf=new T[num_frames*vec_size];
  T* x=tmp_sen_buffer->ptr;

  // copy buffer into a temporary one
  for(unsigned i=0;i<num_frames;++i)  
    for(unsigned j=0;j<vec_size;++j) {
	tmp_buf[i*vec_size+j]=x[i*stride+j];
    }

  // No need to resize tmp_sen_buffer because it is already assigned the proper size 

  unsigned cnt=0;
  for(Range::iterator frame_it = preTransFrameRange->begin(); !frame_it.at_end();++frame_it,++cnt) {
      for(unsigned j=0;j<vec_size;++j) { 
	x[cnt*stride+j]=tmp_buf[*frame_it*vec_size+j];
      }
  }  

  delete [] tmp_buf;
  delete preTransFrameRange;

  return new_num_frames;

}



////////////////////////////////////////////////
// The global matrix object, must be
// actually defined near where main() is defined.

extern ObservationMatrix globalObservationMatrix;

#endif


