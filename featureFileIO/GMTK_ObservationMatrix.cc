/*
 * Written by Katrin Kirchhoff <katrin@ee.washington.edu>
 *
 * Modified by Karim Filali <karim@cs.washington.edu> on 08sep2003:
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
 *   (to achieve downsampling for example or remove parts of the data)
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
 *        d) ARMA filtering
 *
 *        e) General filtering.  The filter is either read from a file or from pre-stored filters (each of which corresponds to some downsampling factor)
 *
 *        f) Add an offset
 *
 *        g) Mean substraction only
 *
 *      a, b, d, e, f and g can also be applied as a last global transformation after prrng and frame adjustement.  
 *
 *   - Streams can be combined instead of concatenated.  The following operations are supported:
 *      
 *      addition, substraction, multiplication, and division.       
 *
 *   TODO: 
 *    1- Add deltas and deltas
 *    2- Add flatascii and flatbin formats
 *   (could treat them as single sentences or have each frame prefixed
 *   by the sentence number and maybe the frame number)
 *    3- Add support for the WAVEFORM HTK paramter kind 
 *   
 *            
 *   Dec 02, 2003 -- Fixed discrete data reading bug when using the
 *   HTK format.  HTK uses shorts not ints.
 *   
 * $Header$
 *
 * Copyright (c) 2001
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about the suitability of this software
 * for any purpose. It is provided "as is" without express or implied warranty.
 *
 * for any purpose. It is provided "as is" without express or implied warranty.
 * */

#include <stdio.h>
#include "GMTK_ObservationMatrix.h"

#define INITIAL_BUF_NUM_FRAMES 1000 


// karim - 29aug2003
#include "fileParser.h"
#ifdef PIPE_ASCII_FILES_THROUGH_CPP
#ifndef DECLARE_POPEN_FUNCTIONS_EXTERN_C
// extern "C" {
//   FILE     *popen(const char *, const char *) __THROW;
//   int pclose(FILE *stream) __THROW;
// };
#endif
#endif
#ifdef PIPE_ASCII_FILES_THROUGH_CPP
#define CPP_DIRECTIVE_CHAR '#'
#endif

#define WARNING_ON_NAN 1

//#define DEBUG

#ifdef DEBUG
#define DBGFPRINTF(_x_) fprintf _x_
#else
#define DBGFPRINTF(_x_)
#endif






// create ObservationMatrix object
ObservationMatrix::ObservationMatrix() {

  _numFrames      = _numNonSkippedFrames = 0;

  _numContinuous  = 0;
  _numDiscrete    = 0;
  _numFeatures    = 0;
  _maxContinuous  = 0;
  _maxDiscrete    = 0;

  _stride         = 0;

  _segmentNumber  = 0;
  _startSkip      = _endSkip = _totalSkip = 0;

  _bufSize        = 0;

  _numSegments    = 0;

  _numStreams     = 0;
  
  _inStreams      = NULL;

  _filterFileName = NULL;
}

// initializes input streams, allocates feature buffer 

void ObservationMatrix::openFiles(int n_files,
			     const char **fof_names,
			     const char **cont_range_str,
			     const char **disc_range_str,
			     unsigned *n_floats,
			     unsigned *n_ints,
			     unsigned *formats,
			     bool *swapflag,
			     const unsigned _startSkipa,
			     const unsigned _endSkipa,
			     bool cppIfAscii,
			     char* cppCommandOptions,
			     const char **pr_range_str,
			     unsigned* actionIfDiffNumFrames,
			     unsigned* actionIfDiffNumSents,
			     char** perStreamPreTransforms,
			     char* postTransforms,
			     unsigned ftrcombo) 
{
  
  assert (n_files > 0);

  unsigned max_num_segments=0;

  _prrngStr=pr_range_str;
  _actionIfDiffNumFrames=actionIfDiffNumFrames;
  _actionIfDiffNumSents=actionIfDiffNumSents;
  _perStreamPreTransforms=perStreamPreTransforms;
  _postTransforms=postTransforms;

  _numStreams = n_files;

  _ftrcombo=ftrcombo;

  // obligatory info

  if (fof_names == NULL)
    error("ObservationMatrix::openFiles: list of file names is NULL\n");

  if (formats == NULL)
    error("ObservationMatrix::openFiles: list of file formats is NULL\n");

  if (n_floats == NULL)
    error("ObservationMatrix::openFiles: list of number of floats is NULL\n");

  if (n_ints == NULL)
    error("ObservationMatrix::openFiles: list of number of ints is NULL\n");

  _inStreams = new StreamInfo*[_numStreams];
  for (unsigned i = 0; i < _numStreams; i++) {
    _inStreams[i] = NULL;
  }

  // create stream info for each stream

  for (unsigned i = 0; i < _numStreams; i++) {

    if (fof_names[i] == NULL)
      error("ObservationMatrix::openFiles: file list for stream %i is NULL\n",i);

    // when range string is missing, assume default = all features

    const char *crng, *drng;
    
    if (cont_range_str == NULL || cont_range_str[i] == NULL )
      crng = "all"; 
    else
      crng = cont_range_str[i];

    if (disc_range_str == NULL || disc_range_str[i] == NULL)
      drng = "all";
    else
      drng = disc_range_str[i];


    // assume default = no swapping

    bool sflag;

    if (swapflag == NULL || swapflag[i] == false)
      sflag = false; 
    else
      sflag = swapflag[i];

    // If we want to add deltas, something needs to be done about the
    // cont range: it wouldn't be valid to have a range extend into
    // the deltas for example.  Therefore we need to know about
    // delats at this point already

    _inStreams[i] = new StreamInfo(fof_names[i],crng,drng,&n_floats[i],&n_ints[i],&formats[i],sflag,i,cppIfAscii,cppCommandOptions);

    // check stream sizes and take action according to the value of _actionIfDiffNumSents
    if(_actionIfDiffNumSents==NULL || _actionIfDiffNumSents[i]==ERROR  || _actionIfDiffNumSents[i]==TRUNCATE_FROM_END) {
      if (i > 0 ) {
	size_t a = _inStreams[i-1]->fofSize;
	size_t b = _inStreams[i]->fofSize;
	if (a != b) {
	  if(_actionIfDiffNumSents!=NULL && _actionIfDiffNumSents[i]==ERROR) {
	    error("ERROR: ObservationMatrix: different number of files in '%s' (%li) and '%s' (%li)\n",fof_names[i-1],a,fof_names[i], b);
	  }
	  else {
	    warning("WARNING ObservationMatrix: different number of files in '%s' (%li) and '%s' (%li) - will only read minimum number\n",fof_names[i-1],a,fof_names[i], b);
	  
	    a < b ? _inStreams[i]->fofSize = a : _inStreams[i-1]->fofSize = b;
	  }
	}
      }
      max_num_segments =  _inStreams[0]->fofSize; // same for all streams;
    }
    else {
      //find maximumm number of sentences/segments
      if(_inStreams[i]->fofSize > max_num_segments)
	max_num_segments=_inStreams[i]->fofSize;
    }

    if (n_floats[i] > _maxContinuous)
      _maxContinuous = n_floats[i];

    if (n_ints[i] > _maxDiscrete)
      _maxDiscrete = n_ints[i];

    // add up features used for observation matrix
    _numContinuous += _inStreams[i]->getNumFloatsUsed();
    _numDiscrete += _inStreams[i]->getNumIntsUsed();
#if !(ALLOW_VARIABLE_DIM_COMBINED_STREAMS)
    if(ftrcombo!=FTROP_NONE && i >0) {
      if(n_floats[i] != n_floats[i-1]) {
	error("ERROR:  When doing feature combination, the number of floats across streams has to be the same.\n");
      }
    }
#endif
  }  // end of for (unsigned i = 0; i < _numStreams; i++)

  if(ftrcombo!=FTROP_NONE) {
    // we'll combine features streams so the overall number of floats will be the max
    // labels are not affected by the feature combination.  Could be changed in the future.
    _numContinuous=_maxContinuous;
  }

  _numSegments = max_num_segments;

  _numFeatures = _numContinuous + _numDiscrete;

  if(_numFeatures == 0) { // possible this can happen if we all the ranges are -1 (i.e. special case in which no feature is selected)
    error("ERROR: No features (continuous or discrete) were selected.  Check the feature ranges.");
  }

  _stride = _numFeatures; 

  _startSkip = _startSkipa;
  _endSkip = _endSkipa;
  _totalSkip = _startSkip + _endSkip;

  // initialize feature buffer 
  _bufSize = INITIAL_BUF_NUM_FRAMES;

#if ALLOW_VARIABLE_DIM_COMBINED_STREAMS
  features.resizeAndZero(_bufSize * _stride); 
#else
  features.resize(_bufSize * _stride); 
#endif
  featuresBase = features.ptr + _stride*_startSkip;

  DBGFPRINTF((stderr,"Allocating %dx%d=%d entries each temp float buffer\n",_bufSize,_maxContinuous,_bufSize * _maxContinuous));
  _tmpFloatSenBuffer.resize(_bufSize * _maxContinuous);
  _tmpIntSenBuffer.resize(_bufSize * _maxDiscrete);

  _repeat.resize(_bufSize);

  //  _contFea.resize(_maxContinuous); // temporary buffers for 1 frame of input
  //_discFea.resize(_maxDiscrete);   
  
  //_cont_p = features.ptr;  // pointer to continuous block 
  //_disc_p = features.ptr + _numContinuous; // pointer to discrete block

  // karim - 29aug2003
  _cppIfAscii=cppIfAscii;
  _cppCommandOptions=cppCommandOptions;
}

ObservationMatrix::~ObservationMatrix() {
  for (unsigned i = 0; i < _numStreams; i++)
    delete _inStreams[i];
  delete [] _inStreams;
  
  if(  _filterFileName != NULL ) 
    delete [] _filterFileName;
}

/**
 *  returns true if true if all streams have the same number of frames
 *  at segno.  Else, returns false (if the action when the number of
 *  frames is different is ERROR, we exit with an error) 
 *
 * max_num_samples get the maximum number of frames over all streams.
 *
 * prrng_num_samples gets either the maximum or the minimum number of
 * frames over all rangers AFTER the per-sentence range has been
 * applied.  
 * */

bool ObservationMatrix::checkIfSameNumSamples(unsigned segno, unsigned& max_num_samples, unsigned& prrng_num_samples){

  bool sameNumSamples=true;
  bool gotAnExpand=false;
  bool gotATruncate=false;
  unsigned i;
  size_t prrng_max_num_samples=0, prrng_min_num_samples=0;
  size_t cur_n_samps = 0;
  size_t cur_prrng_n_samps = 0;
  StreamInfo *s;
  char *fname = "Unknown file";
  char *sname;

  max_num_samples=cur_n_samps;  // i.e. 0
  prrng_max_num_samples=cur_prrng_n_samps;  // i.e. 0

  if (segno < 0 || segno >= _numSegments)  
    error("ObservationMatrix::loadSegment: segment number (%li) outside range of 0 - %li\n",segno,_numSegments);

  for (i = 0; i < _numStreams; i++) {
    s = _inStreams[i];
    if (s == NULL)
      error("ObservationMatrix::loadSegment: stream %i is NULL\n",i);
    
    if (s->fofName == NULL)
      sname = "(unset)";
    else
      sname = s->fofName;

    if(_actionIfDiffNumSents!=NULL &&  segno > s->fofSize-1) {
      if(_actionIfDiffNumSents[i]==REPEAT_LAST) {
	segno=s->fofSize-1;
      }
      else if(_actionIfDiffNumSents[i]==WRAP_AROUND) {
	segno= segno % s->fofSize;
      }
    }

    if (s->dataFormat != PFILE) {
       if(segno >= s->getNumFileNames()) {  // segno starts a 0, hence the >= rather than >
	 error("ObservationMatrix::loadSegment: segment number requested is %d, but list of filenames has only  %d  entries.\n",segno,s->getNumFileNames());
      }

      if (s->dataNames == NULL)
        error("ObservationMatrix::loadSegment: List of file names for stream %i (%s) is NULL\n",i,sname);

      fname = s->dataNames[segno];
      if (fname == NULL)
        error("ObservationMatrix::loadSegment: Filename %li is NULL in stream %i (%s)\n",segno,i,sname);
    }
    else {
      if (s->pfile_istr == NULL)
        error("ObservationMatrix::loadSegment: pfile stream %i (%s) is NULL\n",
              i,sname);
    }

    switch(s->dataFormat) {
    case RAWBIN:
      cur_n_samps = openBinaryFile(s,segno);
      break;
    case RAWASC:
      cur_n_samps = openAsciiFile(s,segno);
      break;
    case HTK:
      cur_n_samps = openHTKFile(s,segno);
      break;
    case PFILE:
      cur_n_samps = openPFile(s,segno);
      break;
    default:
      error("ObservationMatrix::checkIfSameNumSamples: Invalid file format specified for stream %i\n",i);
    }

    DBGFPRINTF((stderr,"In ObservationMatrix::checkIfSameNumSamples, stream %d, segment number %d has %d frames\n",i, segno,cur_n_samps));

    if (cur_n_samps == 0)
      error("ObservationMatrix::CheckNumSamples: failure to read segment %i: no samples found.",i);

    DBGFPRINTF((stderr,"In ObservationMatrix::checkIfSameNumSamples, before cheking per stream transforms.\n"));    

    // check whether there is a transformation that changes the number of frames, e.g: upsampling
    // Post transformations are not allowed to change the number of frames
    if(_perStreamPreTransforms!=NULL && _perStreamPreTransforms[i] != NULL) {
      char* trans_str_ptr =  _perStreamPreTransforms[i];
      int trans;
      int magic_int;
      double magic_double;
      while ( (trans=parseTransform(trans_str_ptr,magic_int,magic_double) ) != END_STR) {
	switch(trans) {
	case UPSAMPLE_HOLD:
	  cur_n_samps *= ( (unsigned) magic_double + 1);
	  break;
	case UPSAMPLE_SMOOTH:
	   cur_n_samps = ( (cur_n_samps-1)*((unsigned)magic_double+1) + 1 );
	  break;
	case DELTAS:
	  // change _numContinuous and _numFeatures and s->nFloats
	  break;
	}
      }
    }

    DBGFPRINTF((stderr,"In ObservationMatrix::checkIfSameNumSamples, before calling s->setAfterTransformCurNumFrames(cur_n_samps);.\n"));    

    s->setAfterTransformCurNumFrames(cur_n_samps);

    DBGFPRINTF((stderr,"In ObservationMatrix::checkIfSameNumSamples, before setting up prrng.\n"));    

    BP_Range* prrng=new BP_Range(_prrngStr==NULL?NULL:_prrngStr[i],0,cur_n_samps);
    DBGFPRINTF((stderr,"In ObservationMatrix::checkIfSameNumSamples, setting up prrng 1.\n"));    
    cur_prrng_n_samps=prrng->length();
    DBGFPRINTF((stderr,"In ObservationMatrix::checkIfSameNumSamples, setting up prrng 2.\n"));    
    s->setPrrngCurNumFrames(cur_prrng_n_samps);
    DBGFPRINTF((stderr,"In ObservationMatrix::checkIfSameNumSamples, setting up prrng 3.\n"));    
    delete prrng;

    DBGFPRINTF((stderr,"In ObservationMatrix::checkIfSameNumSamples, checking if different num frames.\n"));

    if( _actionIfDiffNumFrames != NULL ) {
      if(_actionIfDiffNumFrames[i]==TRUNCATE_FROM_START || _actionIfDiffNumFrames[i]==TRUNCATE_FROM_END) 
	gotATruncate=true;
      if(_actionIfDiffNumFrames[i]==REPEAT_FIRST || _actionIfDiffNumFrames[i]==REPEAT_LAST || _actionIfDiffNumFrames[i]==EXPAND_SEGMENTALLY) 
	gotAnExpand=true;
    }

    // unless we adjust the number of frames, report an error
    if (i > 0 && cur_prrng_n_samps != _inStreams[i-1]->getPrrngCurNumFrames()) {
      if(_actionIfDiffNumFrames==NULL || (_actionIfDiffNumFrames[i]==ERROR && _actionIfDiffNumFrames[i-1]==ERROR)) {
	error("ObservationMatrix::checkIfSameNumSamples: Number of samples for sentence %i don't match for streams %s and %s (%li vs. %li)\n",segno,_inStreams[i-1]->fofName,s->fofName,_inStreams[i-1]->getPrrngCurNumFrames(),cur_prrng_n_samps);
      }
      else {
	if(gotATruncate && gotAnExpand)
	  error("ObservationMatrix: Cannot specify both truncate AND expand actions for when the number of frames is different\n");
	sameNumSamples=false;
      }
    }

    if(cur_n_samps > max_num_samples)
      max_num_samples=cur_n_samps;

    if(cur_prrng_n_samps > prrng_max_num_samples)
      prrng_max_num_samples=cur_prrng_n_samps;

    if(i==0)
      prrng_min_num_samples=cur_prrng_n_samps;
    else if(cur_prrng_n_samps < prrng_min_num_samples)
      prrng_min_num_samples=cur_prrng_n_samps;
    
        DBGFPRINTF((stderr,"In ObservationMatrix::checkIfSameNumSamples, finished processing stream %d, segment number %d\n",i, segno));

  }  // end for (i = 0; i < _numStreams; i++)

  DBGFPRINTF((stderr,"At end of  ObservationMatrix::checkIfSameNumSamples, segment number %d has %d frames\n",segno,cur_n_samps));

  for (unsigned i = 0; i < _numStreams; i++) {
    s = _inStreams[i];
    if(gotAnExpand && s->getPrrngCurNumFrames() < prrng_max_num_samples && ( _actionIfDiffNumFrames==NULL || (_actionIfDiffNumFrames[i]!=REPEAT_FIRST && _actionIfDiffNumFrames[i]!=REPEAT_LAST && _actionIfDiffNumFrames[i]!=EXPAND_SEGMENTALLY))  ) {
      error("ERROR: Stream no %d does not have expand action associated with it and its number of frames, %d, is less than the maximum, %d, across streams",i,s->getPrrngCurNumFrames(),prrng_max_num_samples);
    }
    if(gotATruncate && s->getPrrngCurNumFrames() > prrng_min_num_samples && ( _actionIfDiffNumFrames==NULL || (_actionIfDiffNumFrames[i]!=TRUNCATE_FROM_START && _actionIfDiffNumFrames[i]!=TRUNCATE_FROM_END))) { 
      error("ERROR: Stream no %d does not have truncate action associated with it and its number of frames, %d, is larger than the minimum, %d, across streams",i,s->getPrrngCurNumFrames(),prrng_min_num_samples);
    }
  }

  // prrng_num_samples gets the effective number of samples
  if(gotAnExpand)
    prrng_num_samples=prrng_max_num_samples;
  else
    prrng_num_samples=prrng_min_num_samples;

  return sameNumSamples;
}


// Overloaded method that returns the number of frames in a given segment
unsigned ObservationMatrix::numFrames(unsigned segno) {
  unsigned num_frames=0;


  // max_n_samps and prrng_n_samps are initialized in the next fct call
  unsigned max_n_samps,prrng_n_samps;  
  checkIfSameNumSamples(segno,max_n_samps,prrng_n_samps);

  // maybe move this to the end after potential upsampling takes place
  if (_totalSkip >= prrng_n_samps)
    error("ERROR: number of real frames (%d) for segment %d of input observation files is less than or equal to startSkip + endSkip = %d.",prrng_n_samps,segno,_totalSkip);
  
  num_frames = prrng_n_samps - _totalSkip;

  return num_frames;
}


// Modified to read whole sentences at once instead of frame by frame
// That's also needed because we need to manipulate buffers (to repeat
// frames for example)
//   - karim (karim@cs) 02Sep2003

/* read input for single segment 'segno' */

void ObservationMatrix::loadSegment(unsigned segno) {
  DBGFPRINTF((stderr,"ObservationMatrix::loadSegment: Loading segment no %d\n",segno));

  BP_Range* prrng;
  StreamInfo *s;
  // max_n_samps and prrng_n_samps are initialized in the next fct call
  unsigned max_n_samps,prrng_n_samps;  
  bool same_num_samples = checkIfSameNumSamples(segno,max_n_samps,prrng_n_samps);

  //reset();  // reset the observation buffer to the beginning.
   


  // maybe move this to the end after potential upsampling takes place
  if (_totalSkip >= prrng_n_samps)
    error("ERROR: number of real frames (%d) for segment %d of input observation files is less than or equal to startSkip + endSkip = %d.",prrng_n_samps,segno,_totalSkip);
  

  /////////////////////////////////////////////
  // resize buffer if necessary
  if (max_n_samps > _bufSize) {
    DBGFPRINTF((stderr,"ObservationMatrix::loadSegment: Resizing buffer for segment no %d\n",segno));
      DBGFPRINTF((stderr,"Allocating %dx%dx%d=%d entries for temp float buffer\n",max_n_samps,_maxContinuous,2,max_n_samps * _maxContinuous*2));
    resize(max_n_samps*2);
    _tmpFloatSenBuffer.resize(max_n_samps*_maxContinuous*2);
    _tmpIntSenBuffer.resize(max_n_samps*_maxDiscrete*2);
    _repeat.resize(max_n_samps*2);
  }

  _numNonSkippedFrames = prrng_n_samps;
  
  for (unsigned i = 0; i < _numStreams; i++) {
    s = _inStreams[i];
    unsigned num_floats = s->getNumFloats();
    unsigned num_ints   = s->getNumInts();
    // all checks were done in the function checkNumSamples
    
    DBGFPRINTF((stderr,"In loadSegment(): s->getAfterTransformCurNumFrames() = %d\n",s->getAfterTransformCurNumFrames()));
    prrng=new BP_Range(_prrngStr==NULL?NULL:_prrngStr[i],0,s->getAfterTransformCurNumFrames());
    assert(prrng != NULL);

    //delete prrng;

    if(_actionIfDiffNumSents!=NULL &&  segno > s->fofSize-1) {
      if(_actionIfDiffNumSents[i]==REPEAT_LAST) {
	segno=s->fofSize-1;
      }
      else if(_actionIfDiffNumSents[i]==WRAP_AROUND) {
	segno= segno % s->fofSize;
      }
    }

    // Used to before the above section of code that modifies segno.  It seems obvious it should follow it instead, but the coment is here just in case there was a reason for that I don't see now. 
    _segmentNumber = segno;

    DBGFPRINTF((stderr,"In loadSegment(): Reading sentence.\n"));

    switch(s->dataFormat) {
    case RAWBIN:
    case HTK:
      //readBinSentence(_tmpFloatSenBuffer.ptr,num_floats,_tmpIntSenBuffer.ptr,num_ints,s->curNumFrames,s->curDataFile,s->dataFormat,s->swap());
      readBinSentence(_tmpFloatSenBuffer.ptr,num_floats,_tmpIntSenBuffer.ptr,num_ints,s);
      break;
    case RAWASC:
      readAsciiSentence(_tmpFloatSenBuffer.ptr,num_floats,_tmpIntSenBuffer.ptr,num_ints,s->curNumFrames,s->curDataFile);
      break;
    case PFILE:
      readPfileSentence(segno,_tmpFloatSenBuffer.ptr,_tmpIntSenBuffer.ptr,s->pfile_istr);
      break;
    default:
      error("ObservationMatrix::loadSegment: Invalid file format specified for stream %i\n",i);
    }


#ifdef DEBUG
    DBGFPRINTF((stderr,"In loadSegment(): num_floats = %d and _tmpFloatSenBuffer.ptr = \n",num_floats));
    for(unsigned fr_no=0;fr_no<_numNonSkippedFrames;++fr_no) 
      for(unsigned j=0; j<num_floats;j++)   
	//	DBGFPRINTF((stderr,"%f\n", _tmpFloatSenBuffer.ptr[fr_no*num_floats+j]));
	DBGFPRINTF((stderr,"%f ", _tmpFloatSenBuffer.ptr[fr_no*num_floats+j]));
    DBGFPRINTF((stderr,"\n"));

    DBGFPRINTF((stderr,"In loadSegment(): num_ints = %d and _tmpintSenBuffer.ptr = \n",num_ints));
    for(unsigned fr_no=0;fr_no<_numNonSkippedFrames;++fr_no) 
      for(unsigned j=0; j<num_ints;j++)    
	//	DBGFPRINTF((stderr,"%d", _tmpIntSenBuffer.ptr[fr_no*num_ints+j]));
	DBGFPRINTF((stderr,"%d ", _tmpIntSenBuffer.ptr[fr_no*num_ints+j]));
    DBGFPRINTF((stderr,"\n"));
#endif


    DBGFPRINTF((stderr,"In loadSegment(): Before applying transforms.\n"));

    // Apply transformations.  For now we'll do in place
    // transformations one at a time even if it would be more
    // efficient to combine all transformations as it is done in the
    // copy...() functions below
    if(_perStreamPreTransforms!=NULL && _perStreamPreTransforms[i] != NULL) {
      applyTransforms(_perStreamPreTransforms[i],num_floats, num_ints, s->curNumFrames);
    }

    DBGFPRINTF((stderr,"In loadSegment(): preparing to copy data to final buffer.\n"));

    if(same_num_samples)
      copyToFinalBuffer(i,_tmpFloatSenBuffer.ptr,_tmpIntSenBuffer.ptr,s->cont_rng,s->disc_rng,prrng);
    else
      copyAndAdjustLengthToFinalBuffer(i,_tmpFloatSenBuffer.ptr,_tmpIntSenBuffer.ptr,s->cont_rng,s->disc_rng,prrng,prrng_n_samps);

    DBGFPRINTF((stderr,"In loadSegment(): Copied data in stream %d.  Deleting prrng...\n",i));

    delete prrng;

    DBGFPRINTF((stderr,"In loadSegment(): Deleted prrng.\n"));
  }  //end for(int i = 0; i < _numStreams; i++) loop

  _numFrames = _numNonSkippedFrames - _totalSkip;

  DBGFPRINTF((stderr,"In loadSegment(): Applying post transforms.\n"));  
  if(_postTransforms!=NULL) {
    applyPostTransforms(_postTransforms,numContinuous(), numDiscrete(), _numNonSkippedFrames);
  }

  DBGFPRINTF((stderr,"In loadSegment(): Closing data files.\n"));
  closeDataFiles();
}



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

double conv2double(char* str, unsigned& len, char delimiter,bool conv2int=false) {
#define NUMBER_STRING_MAX_LEN 20
  char* str_ptr=str;
  char new_str[NUMBER_STRING_MAX_LEN];
  double val;

  DBGFPRINTF((stderr,"str=%s\n",str));

  int i=0;
  while(*str_ptr!='_' && *str_ptr!='\0') { 
    new_str[i]=*str_ptr; str_ptr++; i++;
    if(i>=NUMBER_STRING_MAX_LEN) 
      error("ERROR: Number string too long while converting to float\n");
  }
  
  new_str[i]='\0';
  if(conv2int)
    val=parse_long(new_str);
  else
    val=parse_float(new_str);
  
  len=i;

  return (val);
  
}
 
/**
 * parses a transformation string and returns the next transformation
 * to perform.  Advances the string pointer to the next
 * transformation.  
 *
 * side effects: string is modified to point to the next transformation substring.
 *
 * returns END_STR (-1) if the end of the string is reached,
 * UNRECOGNIZED_TRANSFORM (-2) if the transformation substring is not
 * recognized.  Otherwise returns the unsigned that corresponds to the
 * transformation substring.
 * */
int ObservationMatrix::parseTransform(char*& trans_str, int& magic_int, double& magic_double) {

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
    trans_str+=len;
    if(*trans_str=='_') ++trans_str;  // get rid of separator
    return MULTIPLY;
  case TRANS_OFFSET_LETTER: // 'O'
    ++trans_str; 
    // get offset
    magic_double=conv2double(trans_str,len,'_');
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
    trans_str+=len;
    if(*trans_str=='_') ++trans_str;  // get rid of separator
    return return_val;
  case TRANS_ARMA_LETTER:
    ++trans_str; 
    //DBGFPRINTF((stderr,"trans_str=%s\n",trans_str));
    // get order of ARMA filter
    magic_int=(int)conv2double(trans_str,len,'_',true); // conv2int is true
    trans_str+=len;
    if(*trans_str=='_') ++trans_str; 
    return ARMA;
  case FILTER_LETTER:
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
      unsigned tmpStringLen=strlen(tmpString);
      if(tmpStringLen==0) error("ERROR: parseTransform: no filter file name specified\n");
      _filterFileName=new char[tmpStringLen];
      strcpy(_filterFileName,tmpString);
    }
    else if(*trans_str=='_') ++trans_str;
    
    return FILTER;
  case NONE_LETTER:
    ++trans_str; 
    if(*trans_str=='_') ++trans_str; 
    return NONE;
  default:
    error("ERROR: parseTransform: Unrecognized transformation substring (%s)\n",trans_str);
    //DBGFPRINTF((stderr,"In parseTransform: Unrecognized transform @ trans_str=%s\n",trans_str));
    //return UNRECOGNIZED_TRANSFORM;
    
  }

  return UNRECOGNIZED_TRANSFORM;

}

/**
 * apply transformations to the buffers float_buf and int_buf
 *
 *
 *
 */
void ObservationMatrix::applyTransforms(char* trans_str, unsigned num_floats, unsigned num_ints, unsigned num_frames) {
  
  if(trans_str==NULL)
    return;  // nothing to do 

  DBGFPRINTF((stderr,"At start of ObservationMatrix::applyTransforms: trans_str=%s\n",trans_str));

  char* trans_str_ptr=trans_str;
  int trans;
  int magic_int; // a number that holds a value useful for different transformations below
  double magic_double;
  unsigned filter_len=0;

  // parse the transformation string and apply the transformation in order
  while ( (trans=parseTransform(trans_str_ptr,magic_int,magic_double) ) != END_STR) {
    DBGFPRINTF((stderr,"In ObservationMatrix::applyTransforms: trans no = %d (trans_str_ptr=%s)\n",trans,trans_str_ptr));
    switch(trans) {
    case UPSAMPLE_HOLD:
      upsampleHold(&_tmpFloatSenBuffer,num_floats,num_floats,num_frames,(unsigned)magic_double);
      upsampleHold(&_tmpIntSenBuffer,num_ints,num_ints,num_frames,(unsigned)magic_double);
      num_frames *= ( (unsigned) magic_double + 1);  // needed in case another upsample op follows
      break;
    case UPSAMPLE_SMOOTH:
      upsampleSmooth(&_tmpFloatSenBuffer, num_floats,num_floats,num_frames,(unsigned)magic_double);
      upsampleSmooth(&_tmpIntSenBuffer, num_ints,num_ints,num_frames,(unsigned)magic_double);
      num_frames = ( (num_frames-1)*((unsigned)magic_double+1) + 1 );  // needed in case another upsample op follows
      break;
    case NORMALIZE:
      inPlaceMeanSubVarNorm(_tmpFloatSenBuffer.ptr,num_floats,num_floats,num_frames);
      break;
    case MEAN_SUB:
      inPlaceMeanSub(_tmpFloatSenBuffer.ptr,num_floats,num_floats,num_frames);
      break;
    case MULTIPLY:
      multiply(_tmpFloatSenBuffer.ptr,num_floats,num_floats,num_frames,magic_double);
      break;
    case ARMA:
      arma(_tmpFloatSenBuffer.ptr,num_floats,num_floats,num_frames,magic_int);
      break;
    case FILTER:
      if(magic_int==-1)
	readFilterFromFile(_filterFileName, _filterCoeffs, filter_len);
      filter(_tmpFloatSenBuffer.ptr,num_floats,num_floats,num_frames,_filterCoeffs,filter_len);
      break;
    case OFFSET:
      addOffset(_tmpFloatSenBuffer.ptr,num_floats,num_floats,num_frames,magic_double);
      break;
    case UNRECOGNIZED_TRANSFORM:
      error("ERROR: ObservationMatrix::applyTransforms: Unrecognized transformation substring (%s)\n",trans_str_ptr);
    default:
      error("ERROR: ObservationMatrix::applyTransforms: Unrecognized transformation substring (%s)\n",trans_str_ptr);
    }
  }

}



/**
 * apply transformations to the buffer float_buf as a last step after
 * initial transformations, frame selection and adjustement.
 *
 *
 * */
void ObservationMatrix::applyPostTransforms(char* trans_str, unsigned num_floats, unsigned num_ints, unsigned num_frames) {
  
  if(trans_str==NULL)
    return;  // nothing to do 

  DBGFPRINTF((stderr,"At start of ObservationMatrix::applyTransforms: trans_str=%s\n",trans_str));

  char* trans_str_ptr=trans_str;
  int trans;
  int magic_int; // a number that holds a value useful for different transformations below
  double magic_double;
  unsigned filter_len=0;

  // parse the transformation string and apply the transformation in order
  while ( (trans=parseTransform(trans_str_ptr,magic_int,magic_double) ) != END_STR) {
    DBGFPRINTF((stderr,"In ObservationMatrix::applyTransforms: trans no = %d (trans_str_ptr=%s)\n",trans,trans_str_ptr));
    switch(trans) {
    case NORMALIZE:
      inPlaceMeanSubVarNorm((float*)features.ptr,numContinuous(),stride(),num_frames);
      break;
    case MEAN_SUB:
      inPlaceMeanSub((float*)features.ptr,numContinuous(),stride(),num_frames);
      break;
    case ARMA:
      arma((float*)features.ptr,numContinuous(),stride(),num_frames,magic_int);
      break;
    case FILTER:
      if(magic_int==-1)
	readFilterFromFile(_filterFileName, _filterCoeffs, filter_len);
      filter((float*)features.ptr,numContinuous(),stride(),num_frames,_filterCoeffs,filter_len);
      break;
    case OFFSET:
      addOffset((float*)features.ptr,numContinuous(),stride(),num_frames,magic_double);
      break;
    case MULTIPLY:
      multiply((float*)features.ptr,numContinuous(),stride(),num_frames,magic_double);
      break;
    case UNRECOGNIZED_TRANSFORM:
      error("ERROR: ObservationMatrix::applyPostTransforms: Unrecognized transformation substring (%s)\n",trans_str_ptr);
    default:
      error("ERROR: ObservationMatrix::applyPostTransforms: Unrecognized transformation substring (%s)\n",trans_str_ptr);
    }
  }

}


/**
 * applies the continuous features (float_rng), the discrete features
 * (int_rng), and the per-sentence (pr_rng) ranges to the current
 * segment of stream no (stream_no).  Copies the result directly to
 * the final destination buffer.
 *
 * Side effects: none
 *
 * */
void ObservationMatrix::copyToFinalBuffer(unsigned stream_no,float* float_buf,Int32* int_buf,BP_Range* float_rng,BP_Range* int_rng,BP_Range* pr_rng) {

  DBGFPRINTF((stderr,"In ObservationMatrix::copyToFinalBuffer.\n"));

  StreamInfo* s            = _inStreams[stream_no];
  unsigned num_floats      = s->getNumFloats();
  unsigned num_ints        = s->getNumInts();
  unsigned stride          = numFeatures(); 
  unsigned start_float_pos = 0;
  unsigned start_int_pos   = 0;
  bool swap                = s->swap();
  unsigned fmt             = s->getDataFormat();
  if(fmt == PFILE || fmt == RAWASC) swap=false;  // byte swapping is already taken care of for pfiles and is not needed for ascii 

  DBGFPRINTF((stderr,"In ObservationMatrix::copyToFinalBuffer.  After init.\n"));

  void (*copy_swap_func_ptr)(size_t, const intv_int32_t*, intv_int32_t*)=NULL;

  if(stream_no==0) {
    if(swap) copy_swap_func_ptr=&swapb_vi32_vi32;
    else copy_swap_func_ptr=&copy_vi32_vi32;
  }
  else {
    if(swap) {
      switch(_ftrcombo) {
      case FTROP_NONE:
	copy_swap_func_ptr=&swapb_vi32_vi32; break;
      case FTROP_ADD:
	copy_swap_func_ptr=&swapb_add_vi32_vi32; break;
      case FTROP_MUL:
	copy_swap_func_ptr=&swapb_mul_vi32_vi32; break;
      case FTROP_SUB:
	copy_swap_func_ptr=&swapb_sub_vi32_vi32; break;
      case FTROP_DIV:
	copy_swap_func_ptr=&swapb_div_vi32_vi32; break;
      default:
	error("ERROR: Unrecognized feature combination type (%d)",_ftrcombo);
      }
    }
    else {
      switch(_ftrcombo) {
      case FTROP_NONE:
	copy_swap_func_ptr=&copy_vi32_vi32; break;
      case FTROP_ADD:
	copy_swap_func_ptr=&copy_add_vi32_vi32; break;
      case FTROP_MUL:
	copy_swap_func_ptr=&copy_mul_vi32_vi32; break;
      case FTROP_SUB:
	copy_swap_func_ptr=&copy_sub_vi32_vi32; break;
      case FTROP_DIV:
	copy_swap_func_ptr=&copy_div_vi32_vi32; break;
      default:
	error("ERROR: Unrecognized feature combination type (%d)",_ftrcombo);
      }
    }
  }

  if(_ftrcombo==FTROP_NONE) {
    for (unsigned i=0; i < stream_no; ++i) {
      start_float_pos += _inStreams[i]->getNumFloatsUsed();
      start_int_pos += _inStreams[i]->getNumIntsUsed();
    }
  }
  else {
    // we combine float features, so the start is always 0
    start_float_pos = 0;
    // labels are unaffected...for now.
    for (unsigned i=0; i < stream_no; ++i) {
      start_int_pos += _inStreams[i]->getNumIntsUsed();
    }
  }

  DBGFPRINTF((stderr,"In ObservationMatrix::copyToFinalBuffer.  Before copy.\n"));

  float* start_float_buf = getPhysicalStartOfFloatFeaturesBuffer() + start_float_pos;
  float* float_buf_ptr;
  Int32* start_int_buf   = getPhysicalStartOfIntFeaturesBuffer()   + start_int_pos;
  Int32* int_buf_ptr;

  DBGFPRINTF((stderr,"In ObservationMatrix::copyToFinalBuffer.  Copying floats (num_floats = %d).\n",num_floats));  

  if(num_floats>0) {
    for (BP_Range::iterator pr_it = pr_rng->begin(); !pr_it.isAtEnd(); pr_it++,start_float_buf+=stride) {
      float_buf_ptr=start_float_buf;
      for(BP_Range::iterator it = float_rng->begin(); !it.isAtEnd(); it++) {
	copy_swap_func_ptr(1,(const int *) &float_buf[*it+(*pr_it)*num_floats], (int *)float_buf_ptr++);
	if(!finite(*(float_buf_ptr-1))) {
#ifdef WARNING_ON_NAN
	  warning("WARNING: ObservationMatrix::copyToFinalBuffer: found NaN or +/-INF at %i'th float in frame %i, sentence %i  and stream %i\n",*it,*pr_it,_segmentNumber,stream_no);
#else
	  error("ERROR: ObservationMatrix::copyToFinalBuffer: found NaN or +/-INF at %i'th float in frame %i, sentence %i and stream %i\n",*it,*pr_it,_segmentNumber,stream_no);
#endif
	}
      }
    }
  }

  DBGFPRINTF((stderr,"In ObservationMatrix::copyToFinalBuffer.  Copying ints (num_ints = %d).\n",num_ints));

  if(num_ints>0) {
    for (BP_Range::iterator pr_it = pr_rng->begin(); !pr_it.isAtEnd(); pr_it++,start_int_buf+=stride) {
      int_buf_ptr=start_int_buf;
      for(BP_Range::iterator it = int_rng->begin(); !it.isAtEnd(); it++) {
	copy_swap_func_ptr(1,(const int *) &int_buf[*it+(*pr_it)*num_ints], (int *)int_buf_ptr++);
      }
    }
  }

    DBGFPRINTF((stderr,"In ObservationMatrix::copyToFinalBuffer.  END.\n"));

}

/**
 * applies the continuous features (float_rng), the discrete features
 * (int_rng), and the per-sentence (pr_rng) ranges and adjusts the
 * length of the current segment of stream no (stream_no).  Copies the
 * result directly to the final destination buffer
 *
 * Side effects: none
 *
 * */
void ObservationMatrix::copyAndAdjustLengthToFinalBuffer(unsigned stream_no,float* float_buf,Int32* int_buf,BP_Range* float_rng,BP_Range* int_rng,BP_Range* pr_rng,unsigned prrng_n_samps) {

  StreamInfo* s            = _inStreams[stream_no];
  unsigned num_floats      = s->getNumFloats();
  unsigned num_ints        = s->getNumInts();
  unsigned stride          = numFeatures(); 
  unsigned start_float_pos = 0;
  unsigned start_int_pos   = 0;
  bool swap                = s->swap();
  unsigned fmt             = s->getDataFormat();
  if(fmt == PFILE || fmt == RAWASC) swap=false;  // byte swapping is already taken care of for pfiles and is not needed for ascii 

  unsigned cur_stream_prrng_num_frames=s->getPrrngCurNumFrames();
  unsigned diff_in_num_frames=prrng_n_samps-cur_stream_prrng_num_frames;
  int num_per_frame =  prrng_n_samps / cur_stream_prrng_num_frames;
  int rem =   prrng_n_samps % cur_stream_prrng_num_frames;
  // find out how to repeat/truncate frames
  if( _actionIfDiffNumFrames!=NULL ) {
     switch(_actionIfDiffNumFrames[stream_no]) {
     case ERROR:
       _repeat[0]=cur_stream_prrng_num_frames;
       _repeat[1]=1;
       _repeat[2]=-1;
       break;
     case REPEAT_LAST:
       _repeat[0]=cur_stream_prrng_num_frames-1;
       _repeat[1]=1; // print all but the last frame once
       _repeat[2]=1;
       _repeat[3]=diff_in_num_frames+1; // last frame is printed diff_in_num_frames+1 times
       _repeat[4]=-1;  //we stop here
       break;
     case REPEAT_FIRST:
       _repeat[0]=1;
       _repeat[1]=diff_in_num_frames+1; // last frame is printed diff_in_num_frames+1 times
       _repeat[2]=cur_stream_prrng_num_frames-1;
       _repeat[3]=1; // print all but the last frame once
       _repeat[4]=-1;  //we stop here
       break;
     case EXPAND_SEGMENTALLY:
       _repeat[0]=rem;
       _repeat[1]=num_per_frame+1;
       _repeat[2]=cur_stream_prrng_num_frames-rem;
       _repeat[3]=num_per_frame; // last frame is printed diff_in_num_frames+1 times
       _repeat[4]=-1;  //we stop here
       break;
     case TRUNCATE_FROM_END:
       _repeat[0]=prrng_n_samps;
       _repeat[1]=1;
       _repeat[2]=-diff_in_num_frames;
       _repeat[3]=0;
       _repeat[4]=-1;  //we stop here
       break;
     case TRUNCATE_FROM_START:
       _repeat[0]=-diff_in_num_frames;
       _repeat[1]=0; // delete extra frames
       _repeat[2]=prrng_n_samps;
       _repeat[3]=1;
       _repeat[4]=-1;  //we stop here
       break;
     default:
       error("ObservationMatrix:::copyAndAdjustLengthToFinalBuffer: Invalid action (%d) for when the number of frames in sentences is different across streams",_actionIfDiffNumFrames[stream_no]);
     }
  }
  else {
    //same as the ERROR case above
    _repeat[0]=cur_stream_prrng_num_frames;
    _repeat[1]=1;
    _repeat[2]=-1;
  }

  void (*copy_swap_func_ptr)(size_t, const intv_int32_t*, intv_int32_t*)=NULL;

  if(stream_no==0) {
    copy_swap_func_ptr=swap?&swapb_vi32_vi32:&copy_vi32_vi32;
  }
  else {
    if(swap) {
      switch(_ftrcombo) {
      case FTROP_NONE:
	copy_swap_func_ptr=&swapb_vi32_vi32; break;
      case FTROP_ADD:
	copy_swap_func_ptr=&swapb_add_vi32_vi32; break;
      case FTROP_MUL:
	copy_swap_func_ptr=&swapb_mul_vi32_vi32; break;
      case FTROP_SUB:
	copy_swap_func_ptr=&swapb_sub_vi32_vi32; break;
      case FTROP_DIV:
	copy_swap_func_ptr=&swapb_div_vi32_vi32; break;
      default:
      error("ERROR: Unrecognized feature combination type (%d)",_ftrcombo);
      }
      //copy_swap_func_ptr=&swapb_vi32_vi32;
    }
    else {
      switch(_ftrcombo) {
      case FTROP_NONE:
      copy_swap_func_ptr=&copy_vi32_vi32; break;
      case FTROP_ADD:
      copy_swap_func_ptr=&copy_add_vi32_vi32; break;
      case FTROP_MUL:
      copy_swap_func_ptr=&copy_mul_vi32_vi32; break;
      case FTROP_SUB:
	copy_swap_func_ptr=&copy_sub_vi32_vi32; break;
      case FTROP_DIV:
      copy_swap_func_ptr=&copy_div_vi32_vi32; break;
      default:
	error("ERROR: Unrecognized feature combination type (%d)",_ftrcombo);
      }
      //copy_swap_func_ptr=&copy_vi32_vi32;
    }
  }

 if(_ftrcombo==FTROP_NONE) {
    for (unsigned i=0; i < stream_no; ++i) {
      start_float_pos += _inStreams[i]->getNumFloatsUsed();
      start_int_pos += _inStreams[i]->getNumIntsUsed();
    }
  }
  else {
    // we combine float features, so the start is always 0
    start_float_pos = 0;
    // labels are unaffected...for now.
    for (unsigned i=0; i < stream_no; ++i) {
      start_int_pos += _inStreams[i]->getNumIntsUsed();
    }
  }

  float* start_float_buf = getPhysicalStartOfFloatFeaturesBuffer() + start_float_pos;
  float* float_buf_ptr;
  Int32* start_int_buf   = getPhysicalStartOfIntFeaturesBuffer()   + start_int_pos;
  Int32* int_buf_ptr;

  if(num_floats>0) {
    int repeat, cnt=0;
    int* rep_ptr=_repeat.ptr;
    for (BP_Range::iterator pr_it = pr_rng->begin(); !pr_it.isAtEnd(); pr_it++,cnt++) {
      if(cnt>=*rep_ptr) {
	cnt=0;
	rep_ptr+=2;
      }
      repeat=*(rep_ptr+1);
      assert(*rep_ptr != -1);  // sanity check
      for(int i=0; i<repeat;++i) {
	float_buf_ptr=start_float_buf;
	for(BP_Range::iterator it = float_rng->begin(); !it.isAtEnd(); it++) {
	  copy_swap_func_ptr(1,(const int *) &float_buf[*it+(*pr_it)*num_floats], (int *)float_buf_ptr++);
	  if(!finite(*(float_buf_ptr-1))) {
#ifdef WARNING_ON_NAN
	    warning("WARNING: ObservationMatrix::copyToFinalBuffer: found NaN or +/-INF at %i'th float in frame %i and stream %i\n",*it,*pr_it,stream_no);
#else
	    error("ERROR: ObservationMatrix::copyToFinalBuffer: found NaN or +/-INF at %i'th float in frame %i and stream %i\n",*it,*pr_it,stream_no);
#endif
	  }
	}
	start_float_buf+=stride;
      }
    }
  }
  if(num_ints>0) {
    int repeat, cnt=0;
    int* rep_ptr=_repeat.ptr;
    for (BP_Range::iterator pr_it = pr_rng->begin(); !pr_it.isAtEnd(); pr_it++,cnt++) {
      if(cnt>=*rep_ptr) {
	cnt=0;
	rep_ptr+=2;
      }
      repeat=*(rep_ptr+1);
      assert(*rep_ptr != -1);  // sanity check
      for(int i=0; i<repeat;++i) {
	int_buf_ptr=start_int_buf;
	for(BP_Range::iterator it = int_rng->begin(); !it.isAtEnd(); it++) {
	  copy_swap_func_ptr(1,(const int *) &int_buf[*it+(*pr_it)*num_ints], (int *)int_buf_ptr++);
	}
	start_int_buf+=stride;
      }
    }
  }
}



/* reset pointers to beginning of data buffer */

void ObservationMatrix::reset() {

  //_cont_p = features.ptr;
  //_disc_p = features.ptr + _numContinuous;
}


/* resize data buffer
 * n_frames: new number of frames
 */

void ObservationMatrix::resize(size_t n_frames) 
{
  features.growIfNeeded(n_frames * _stride);
  featuresBase = features.ptr + _stride*_startSkip;
  _bufSize = n_frames;
}

/* print frame from data buffer
 * stream: output stream
 * frameno: frame number
 */

void ObservationMatrix::printFrame(FILE *stream, size_t absoluteFrameno) {
  
  unsigned f;

  assert(absoluteFrameno < _numNonSkippedFrames && absoluteFrameno >= 0);

  unsigned offset = _stride*absoluteFrameno;

  Data32 *p = features.ptr + offset;

  for (f = 0; f < _numContinuous; f++,p++) {
    float *fp = (float *)p;
    fprintf(stream,"%e ",*fp);
  }
  for (f = 0; f < _numDiscrete; f++,p++) {
    Int32 *ip = (Int32 *)p;
    fprintf(stream,"%i ",*ip);
  }
  fprintf(stream,"\n");

}

/* get info for segment 'sentno' from pfile stream 'f' */

size_t ObservationMatrix::openPFile(StreamInfo *f, size_t sentno) {

  assert(sentno >= 0 && sentno < _numSegments);

  if (f->pfile_istr == NULL) {
    error("ObservationMatrix::openPFile: stream is NULL");
    return 0;
  }

  if (f->pfile_istr->set_pos(sentno,0) == SEGID_BAD) {
    warning("ObservationMatrix::openPFile: Can't skip to sent %li frame 0",
	    sentno);
    return 0;
  }

  f->curNumFrames = f->pfile_istr->num_frames(sentno);

  return f->curNumFrames;
}

/* open binary file for segment 'sentno' */


size_t ObservationMatrix::openBinaryFile(StreamInfo *f, size_t sentno) {

  assert(sentno >= 0 && sentno < _numSegments);

  char *fname = f->dataNames[sentno];
  int nfloats = f->nFloats;
  int nints = f->nInts;
  size_t fsize;

  if (fname == NULL) {
    warning("ObservationMatrix::openBinaryFile: Filename is NULL for segment %li\n",sentno);
    return 0;
  }
  
  if ((f->curDataFile = fopen(fname,"rb")) == NULL) {
    error("ObservationMatrix::openBinaryFile: Can't open '%s' for input\n",
	  fname);
	  
  }
    
    //sanity check on number of bytes
    
    if (fseek(f->curDataFile,0L,SEEK_END) == -1) {
      warning("ObservationMatrix::openBinaryFile: Can't skip to end of file %s",
	      fname);
    }
    
    fsize = ftell(f->curDataFile);
    
    rewind(f->curDataFile);
    
    int rec_size = nfloats * sizeof(float) + nints * sizeof(int);
    
    if ((fsize % rec_size) > 0)
      error("ObservationMatrix::openBinaryFile: odd number of bytes in file %s\n",fname);
    
    int n_samples = fsize / rec_size;
    
    f->curNumFrames = n_samples;
    
    return n_samples;
}

/* open ascii file for segment 'sentno' */

size_t ObservationMatrix::openAsciiFile(StreamInfo *f,size_t sentno) {
  bool found_new_line=true;
  char *fname = f->dataNames[sentno];
  size_t n_samples = 0;
  char ch;

  if (fname == NULL) {
    warning("ObservationMatrix::openAsciiFile: Filename is NULL for segment %li\n",sentno);
    return 0;
  }

#ifdef PIPE_ASCII_FILES_THROUGH_CPP     
     if(_cppIfAscii) { 
       string cppCommand = string("cpp");
       if (_cppCommandOptions != NULL) {
	 cppCommand = cppCommand + string(" ") + string(_cppCommandOptions);
       }
       // make sure the file  exists first.
       if ((f->curDataFile = ::fopen(fname,"r")) == NULL) {
	 error("ERROR: unable to open file (%s) for reading",fname);
       }
       fclose(f->curDataFile);
       cppCommand = cppCommand + string(" ") + string(fname);
       f->curDataFile = ::popen(cppCommand.c_str(),"r");    
       if (f->curDataFile == NULL)
	 error("ERROR, can't open file (%s)",f->curDataFile);

       DBGFPRINTF((stderr,"In ObservationMatrix::openAsciiFile(), the cpp processed input ascii file is\n"));
       int tmp;
       while ((tmp = fgetc(f->curDataFile)) != EOF) {
	 DBGFPRINTF((stderr,"%c",tmp));
	 ch=tmp;
	 if(ch=='\n') { found_new_line=true; continue;}
	 if(ch==CPP_DIRECTIVE_CHAR) { found_new_line=false; continue; }
	 if(found_new_line) {
	   found_new_line=false;
	   n_samples++;
	 }
       }
       	 DBGFPRINTF((stderr,"\n"));
       fclose(f->curDataFile);
       f->curDataFile = ::popen(cppCommand.c_str(),"r");    
     }
     else {
       if ((f->curDataFile = fopen(fname,"r")) == NULL) {
	 warning("ObservationMatrix::openAsciiFile: Can't open '%s' for input\n",fname);
	 return 0;
       }
        /* for ascii, newline is record delimiter - additional or missing nl's will cause error messages.  Fixed this -- karim 07sep2003*/
       int tmp;
       while ((tmp = fgetc(f->curDataFile)) != EOF) {
	 ch = tmp;
	 //if (ch == '\n')
	 //n_samples++;
	 if(ch=='\n') { found_new_line=true; continue;}
	 if(ch==CPP_DIRECTIVE_CHAR) { found_new_line=false; continue; }
	 if(found_new_line) {
	   found_new_line=false;
	   n_samples++;
	 }
       }
       rewind(f->curDataFile);
     }
#else
  if ((f->curDataFile = fopen(fname,"r")) == NULL) {
    warning("ObservationMatrix::openAsciiFile: Can't open '%s' for input\n",fname);
    return 0;
  }
  /* for ascii, newline is record delimiter - additional or missing nl's will cause error messages */
  int tmp;
  while ((tmp = fgetc(f->curDataFile)) != EOF) {
    ch = tmp;
    if(ch=='\n') { found_new_line=true; continue;}
    if(ch==CPP_DIRECTIVE_CHAR) { found_new_line=false; continue; }
    if(found_new_line) {
      found_new_line=false;
      n_samples++;
    }
    //if (ch == '\n')
    //n_samples++;
  }
 rewind(f->curDataFile);
#endif

 f->curNumFrames = n_samples;

 DBGFPRINTF((stderr,"In  ObservationMatrix::openAsciiFile: n_samples = %d\n",n_samples));

 return n_samples;
}

/* open htk file for segment 'sentno' */

size_t ObservationMatrix::openHTKFile(StreamInfo *f, size_t sentno) {

  assert(sentno >= 0 && sentno < _numSegments);

  char *fname = f->dataNames[sentno];

  // structure of HTK header
  Int32 n_samples;
  Int32 samp_period;
  short samp_size;
  short parm_kind;

  Int32 tmp1,tmp2;

  short stmp1,stmp2;

  bool bswap = f->bswap;
  int nints = f->nInts;
  int nfloats = f->nFloats;

  if (fname == NULL) {
    warning("ObservationMatrix::openHTKFile: Filename is NULL for segment %li\n",sentno);
    return 0;

  }

  if ((f->curDataFile = fopen(fname,"rb")) == NULL) {
    warning("ObservationMatrix::openHTKFile: Can't open '%s' for input\n",fname);
    return 0;
  }

  if (fread(&tmp1,sizeof(Int32),1,f->curDataFile) != 1) {
    warning("ObservationMatrix::openHTKFile: Can't read number of samples\n");
    return 0;
  }

  
  if (fread((short *)&tmp2,sizeof(Int32),1,f->curDataFile) != 1) {
    warning("ObservationMatrix::openHTKFile: Can't read sample period\n");
    return 0;
  }

  if (fread((short *)&stmp1,sizeof(short),1,f->curDataFile) != 1) {
    warning("ObservationMatrix::openHTKFile: Can't read sample size\n");
    return 0;
  }

  if (fread(&stmp2,sizeof(short),1,f->curDataFile) != 1) {
    warning("ObservationMatrix::openHTKFile: Can't read parm kind\n");
    return 0;
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
    warning("ObservationMatrix::openHTKFile: number of samples is %i\n",n_samples);
    return 0;
  }

  if (samp_period <= 0 || samp_period > 1000000) {
    warning("ObservationMatrix::openHTKFile: sample period is %i - must be between 0 and 1000000\n", samp_period);
    return 0;
  }

  if (samp_size <= 0 || samp_size > 5000) {
    warning("ObservationMatrix::openHTKFile: sample size is %i - must be between 0 and 5000\n",samp_size);
    return 0;
  }

  short pk = parm_kind & BASEMASK;

  if (pk < WAVEFORM || pk > ANON) {
    warning("Undefined parameter kind for HTK feature file: %i\n",pk);
    return 0;
  }

  // For now we don't support the WAVEFORM parameter kind.  It uses
  // shorts instead of floats and that requires special treatment.
  if (pk == WAVEFORM) {
    warning("HTK WAVEFORM parameter kind not supported: %i\n",pk);
    return 0;
  }


  int n_fea;

  // parameter kind DISCRETE = all discrete features

  if (parm_kind == DISCRETE) {
    //    n_fea = samp_size / sizeof(int) ;
    n_fea = samp_size / sizeof(short) ;
    if (n_fea != nints) {
      warning("ObservationMatrix::openHTKFile:  Number of features in file (%i) does not match number of ints specified (%i)\n", n_fea,nints);
      return 0;
    }
  }
  
  // otherwise all continuous features

  else {
    n_fea = samp_size / sizeof(float);
    if (n_fea != nfloats) {
      warning("ObservationMatrix::openHTKFile:  Number of features in file (%i) does not match number of floats specified (%i)\n", n_fea,nfloats);
      return 0;
    }
  }

  f->curNumFrames = n_samples;

  return n_samples;
}

/* close individual data files, except for pfile. */
void ObservationMatrix::closeDataFiles() {

  for (unsigned i = 0; i < _numStreams; i++)
    if (_inStreams[i]->dataFormat != PFILE)
      fclose(_inStreams[i]->curDataFile);
}



void ObservationMatrix::printSegmentInfo() {

  printf("Processing segment # %d. Number of frames = %d, number non skipped frames = %d.\n",
	 (unsigned)_segmentNumber,_numFrames,_numNonSkippedFrames);
}








/* read a binary sentence into data buffer
 * n_floats: number of floats to be read
 * n_samples: number of frames to be read
 * f: stream to read from
 * returns true if successful, else false
 *
 * Pre-conditions: - either the number of floats or the number of ints has to be positive.
 *                 - the buffers corresponding a non-zero datatypes
 *                 need to be non-null and enough space needs to be
 *                 allocated for them
 */

bool ObservationMatrix::readBinSentence(float* float_buffer, unsigned num_floats, Int32* int_buffer, unsigned num_ints,StreamInfo* s) {
  
  assert(num_floats > 0 || num_ints > 0);

  unsigned n_samples= s->curNumFrames;
  FILE* f=s->curDataFile;
  unsigned data_format=s->dataFormat;
  unsigned swap = s->swap();

  unsigned n_read=0;	
  unsigned total_num_floats=num_floats*n_samples;
  unsigned total_num_ints=num_ints*n_samples;

  // if we have only one type read complete sentence
  if(num_ints==0) {
    assert(float_buffer !=NULL);
    n_read = fread((float *)float_buffer,sizeof(float),total_num_floats,f);
    if (n_read != total_num_floats) {
      warning("ObservationMatrix::readBinFloats: read %i items, expected %i",
	      n_read,total_num_floats);
      return false;
    }
  } 
  else if(num_floats==0) {
    assert(int_buffer != NULL);
    if(data_format==HTK) {
      short* tmp_short_buffer = new short[total_num_ints];
      n_read = fread((short*)tmp_short_buffer,sizeof(short),total_num_ints,f);
      for (unsigned i=0; i<total_num_ints; ++i) {
	if(swap) {
	  tmp_short_buffer[i]=swapb_short_short(tmp_short_buffer[i]); 
	  s->setSwap(false);  // we did the swap here. We don't want to repeat it later on.
	}
	int_buffer[i]=(int)tmp_short_buffer[i];
      }
      delete [] tmp_short_buffer;
    }
    else {
      n_read = fread((Int32*)int_buffer,sizeof(Int32),total_num_ints,f);
    }
    if (n_read != total_num_ints) {
      warning("ObservationMatrix::readBinFloats: read %i items, expected %i",n_read,total_num_ints);
      return false;
    }
  }
  else { // mixed case
    assert(float_buffer != NULL && int_buffer != NULL);
    float* float_buffer_ptr = float_buffer;
    Int32* int_buffer_ptr   = int_buffer;
    for(unsigned s=0; s < n_samples; ++s) {
      n_read += fread((float*)float_buffer_ptr,sizeof(float), num_floats, f);
      n_read += fread((Int32*)int_buffer_ptr,  sizeof(Int32), num_ints,   f);
      float_buffer_ptr += num_floats;
      int_buffer_ptr   += num_ints;
    }
    if (n_read != (total_num_ints + total_num_floats)) {
      warning("ObservationMatrix::readBinFloats: read %i items, expected %i",n_read,total_num_ints+total_num_floats);
      return false;
    }
  }

  // we could do byte swapping here but since we are gonna be copying
  // the buffer anyway later, we can do it at the same time

  return true;
}

/* read an ascii sentence into data buffer
 * n_floats: number of floats to be read
 * n_samples: number of frames to be read
 * f: stream to read from
 * returns true if successful, else false
 *
 * Pre-conditions: - either the number of floats or the number of ints has to be positive.
 *                 - the buffers corresponding a non-zero datatypes
 *                 need to be non-null and enough space needs to be
 *                 allocated for them
 */

bool ObservationMatrix::readAsciiSentence(float* float_buffer, unsigned num_floats, Int32* int_buffer, unsigned num_ints,unsigned n_samples, FILE *f) {

  assert(num_floats > 0 || num_ints > 0);

  int tmp,lineNum=0;
  float* float_buffer_ptr = float_buffer;
  Int32* int_buffer_ptr   = int_buffer;

// consume CPP special directives if any
#ifdef PIPE_ASCII_FILES_THROUGH_CPP     
  if(_cppIfAscii) {
    while((tmp=fgetc(f))==CPP_DIRECTIVE_CHAR) {
      while((tmp=fgetc(f))!='\n');
      lineNum++;
    }
    ungetc(tmp,f);
  }
#endif

  // could be made a bit more efficient since we check whether
  // num_floats and num_ints > 0 for each frame.
  DBGFPRINTF((stderr,"Reading ascii sentence...\n"));
  for(unsigned s=0; s < n_samples; ++s) { 
    DBGFPRINTF((stderr,"%d:  ",lineNum));
    lineNum++;
    if(num_floats > 0) {
      for (unsigned n = 0; n < num_floats; ++n,++float_buffer_ptr) {
	if (fscanf(f,"%e",float_buffer_ptr) != 1) {
	  error("ERROR: ObservationMatrix::readAsciiSentence: couldn't read %i'th item in frame %i\n",n,s);
	}
	DBGFPRINTF((stderr,"%f ",*float_buffer_ptr));
      }
    }
    if(num_ints > 0) {
      for (unsigned n = 0; n < num_ints; ++n,++int_buffer_ptr) {
	if (fscanf(f,"%d",int_buffer_ptr) != 1) {
	    error("ERROR: ObservationMatrix::readAsciiSentence: couldn't read %i'th item in frame %i\n",n,s);
	}
	DBGFPRINTF((stderr,"%d ",*int_buffer_ptr));
      }
    }
    DBGFPRINTF((stderr,"\n"));
  }

  return true;
}

/* read a pfile sentence into data buffer
 * f: stream to read from
 * returns true if successful, else false
 *
 * Pre-conditions: - either the number of floats or the number of ints has to be positive.
 *                 - the buffers corresponding a non-zero datatypes
 *                 need to be non-null and enough space needs to be
 *                 allocated for them
 */

bool ObservationMatrix::readPfileSentence(const unsigned segno, float* float_buffer,  Int32* int_buffer, InFtrLabStream_PFile *f) {
  
  unsigned num_ints    = f->num_labs();
  unsigned num_floats  = f->num_ftrs();
  unsigned num_frames  = f->num_frames(segno); 

  DBGFPRINTF((stderr,"In ObservationMatrix::readPfileSentence, num_ints=%d, num_floats=%d, num_frames=%d\n", num_ints,num_floats,num_frames));

  if(num_ints > 0 and num_floats > 0) {
    if (f->read_ftrslabs((size_t)num_frames,(float *)float_buffer,(UInt32 *)int_buffer) != num_frames) {
      warning("ObservationMatrix::readPfileSentence: read ftrslabs failed\n");
      return 0;
    }
  }
  else {
    if(num_floats>0) {
      if (f->read_ftrs((size_t)num_frames,(float *)float_buffer) != num_frames) {
	warning("ObservationMatrix::readPfileSentence: ftr read failed\n");
	return 0;
      }
    }
    else if(num_ints>0) {
      if (f->read_labs((size_t)num_frames,(UInt32 *)int_buffer) != num_frames) {
	warning("ObservationMatrix::readPfileSentence: lab read failed\n");
	return 0;
      }
    }
  }

    DBGFPRINTF((stderr,"In ObservationMatrix::readPfileSentence, the first 2 frames are:\n"));
#ifdef DEBUG
    for(unsigned f=0; f< 2;++f) {
      for(unsigned i=0; i< num_floats;++i)
	DBGFPRINTF((stderr,"%f ",float_buffer[i+f*num_floats]));
      DBGFPRINTF((stderr,"\n"));
    }
#endif

  return true;
}
