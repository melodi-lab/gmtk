/*
 * Written by Katrin Kirchhoff <katrin@ee.washington.edu>
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
 *
 */

#include <stdio.h>
#include "GMTK_ObservationMatrix.h"

// create ObservationMatrix object

ObservationMatrix::ObservationMatrix() {

  numFrames = 0;
  numContinuous = 0;
  numDiscrete = 0;
  numFeatures = 0;
  stride = 0;
  segmentNumber = 0;

  _totalContinuous = 0;
  _totalDiscrete = 0;
  _maxContinuous = 0;
  _maxDiscrete = 0;
  _bufSize = 0;
  _numSegments = 0;
  _numStreams = 0;
  
  _cont_p = NULL;
  _disc_p = NULL;
  _inStreams = NULL;

}

// initializes input streams, allocates feature buffer 

void
ObservationMatrix::openFiles(int n_files,
			     const char **fof_names,
			     const char **cont_range_str,
			     const char **disc_range_str,
			     unsigned *n_floats,
			     unsigned *n_ints,
			     unsigned *formats,
			     bool *swapflag) {
  
  assert (n_files > 0);

  _numStreams = n_files;

  // obligatory info

  if (fof_names == NULL)
    error("ObservationMatrix::openFiles: list of file names is NULL\n");

  if (formats == NULL)
    error("ERROR DataInStream: list of file formats is NULL\n");

  if (n_floats == NULL)
    error("DataInStream: list of number of floats is NULL\n");

  if (n_ints == NULL)
    error("DataInStream: list of number of ints is NULL\n");

  _inStreams = new StreamInfo*[_numStreams];

  // create stream info for each stream

  for (int i = 0; i < _numStreams; i++) {

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

    if (swapflag == NULL || swapflag[i] == NULL)
      sflag = 0; 
    else
      sflag = swapflag[i];

    _inStreams[i] = new StreamInfo(fof_names[i],
				  crng,
				  drng,
				  &n_floats[i],
				  &n_ints[i],
				  &formats[i],
				  sflag,i);

    // check stream sizes and set to smallest common value

    if (i > 0 ) {

      size_t a = _inStreams[i-1]->fofSize;
      size_t b = _inStreams[i]->fofSize;
      
      if (a != b) {
        warning("WARNING: DataInStream: different number of files in '%s' (%li) an\
d '%s' (%li) - will only read minimum number\n",fof_names[i-1],a,
                fof_names[i], b);
	
        a < b ? _inStreams[i]->fofSize = a : _inStreams[i-1]->fofSize = b;
      }
    }

    // add up features used for observation matrix

    numContinuous += _inStreams[i]->nFloatsUsed;
    numDiscrete += _inStreams[i]->nIntsUsed;

    if (n_floats[i] > _maxContinuous)
      _maxContinuous = n_floats[i];

    if (n_ints[i] > _maxDiscrete)
      _maxDiscrete = n_ints[i];
  }


  _numSegments = _inStreams[0]->fofSize; // same for all streams

  numFeatures = numContinuous + numDiscrete;

  stride = numFeatures; 

  // initialize feature buffer 

  _bufSize = MAXBUFSIZE;
  
  features.resize(_bufSize * stride); 

  _contFea.resize(_maxContinuous); // temporary buffers for 1 frame of input
  _discFea.resize(_maxDiscrete);   

  
  _cont_p = features.ptr;  // pointer to continuous block 
  _disc_p = features.ptr + numContinuous; // pointer to discrete block

}

ObservationMatrix::~ObservationMatrix() {

  for (int i = 0; i < _numStreams; i++)
    delete _inStreams[i];
  delete [] _inStreams;

}

/* read input for single segment 'segno' */

void
ObservationMatrix::loadSegment(const unsigned segno) {

  int i;
  size_t n_samps;
  StreamInfo *s;
  char *fname;


  if (segno < 0 || segno > _numSegments)
    error("ObservationMatrix::loadSegment: segment number (%li) outside range of 0 - %li\n",segno,_numSegments);


  reset();

  for (i = 0; i < _numStreams; i++) {

    s = _inStreams[i];

    char *sname = s->fofName;

    if (s == NULL)
      error("ObservationMatrix::loadSegment: stream %s is NULL\n",sname);

    if (s->dataFormat != PFILE) {

      if (s->dataNames == NULL)
        error("ObservationMatrix::loadSegment: List of file names for stream %s is NULL\n",sname);

      fname = s->dataNames[segno];

      if (fname == NULL)
        error("ObservationMatrix::loadSegment: Filename %li is NULL in stream %s\n",segno,sname);
    }
    else {
      if (s->pfile_istr == NULL)
        error("ObservationMatrix::loadSegment: pfile stream %s is NULL\n",
              s->fofName);
    }
  
    switch(s->dataFormat) {

    case RAWBIN:
      n_samps = openBinaryFile(s,segno);
      break;
    case RAWASC:
      n_samps = openAsciiFile(s,segno);
      break;
    case HTK:
      n_samps = openHTKFile(s,segno);
      break;
    case PFILE:
      n_samps = openPFile(s,segno);
      break;
    default:
      error("ObservationMatrix::loadSegment: Invalid file format specified for stream %i\n",i);
    }

    if (n_samps == 0)
      error("ObservationMatrix::loadSegment: failure to read file '%s'\n",fname);

    if (i > 0 && n_samps != _inStreams[i-1]->curNumFrames) {
      error("ObservationMatrix::loadSegment: Number of samples for sentence %i don't match for streams %s and %s (%li vs. %li)\n",
	    segno,
	    _inStreams[i-1]->fofName,
	    s->fofName,
	    _inStreams[i-1]->curNumFrames,
	    n_samps);
    }
  }
  numFrames = n_samps;

  // resize buffer if necessary

  if (numFrames > _bufSize) 
    resize(numFrames*2);

  for (size_t n = 0; n < numFrames; n++) 
    readFrame(n);
  
  segmentNumber = segno;

  closeDataFiles();
}
  

/* read binary floats into data buffer
 * n_floats: number of floats to be read
 * f: stream to read from
 * cont_rng: range of features to actually use (from total frame read)
 * bswap: bytes will be swapped if true
 * returns 1 if successful, else 0
 */

bool
ObservationMatrix::readBinFloats(unsigned n_floats, FILE *f, 
                                 BP_Range *cont_rng,bool bswap) {

  unsigned n_read;	
  
  
  if (_cont_p == NULL) {
    warning("ObservationMatrix::readBinFloats: Data buffer is NULL\n");
    return 0;
  }

  if (_contFea.ptr == NULL) {
    warning("ObservationMatrix::readBinFloats: Temporary buffer is NULL\n");
    return 0;
  }

  // read complete frame
  
  n_read = fread((float *)_contFea.ptr,sizeof(float),n_floats,f);
  
  if (n_read != n_floats) {
    warning("ObservationMatrix::readBinFloats: read %i items, expected %i",
	    n_read,n_floats);
    return 0;
  }

  // select features
  
  for (BP_Range::iterator it = cont_rng->begin(); it <= cont_rng->end(); it++,_cont_p++) {

    int i = (int) *it;
    float *fp = (float *)_cont_p;

    if (fp == NULL) {
      warning("ObservationMatrix::readBinFloats: memory error reading %i'th feature",i);
      return 0;
    }

    if (bswap) 
      swapb_vi32_vi32(1,(const int *)&_contFea.ptr[i],(int *)fp);
    else
      *fp = _contFea.ptr[i];
    
  }
  return 1;
}

/* read floats from pfile stream
 * str: pfile stream 
 * cont_rng: range of feature to use
 * returns 1 if successful, else 0
 */

bool
ObservationMatrix::readPFloats(InFtrLabStream_PFile *str, BP_Range *cont_rng) {


  unsigned n_read;

  if (_cont_p == NULL) {
    warning("ObservationMatrix::readPFLoats: Data buffer is NULL\n");
    return 0;
  }
  
  if (_contFea.ptr == NULL) {
    warning("ObservationMatrix::readPFLoats: Temporary buffer is NULL\n");
    return 0;
  }
  
  if ((n_read = str->read_ftrslabs((size_t)1,(float *)_contFea.ptr,NULL)) != 1) {
    warning("ObservationMatrix::readPFLoats: read failed\n");
    return 0;
  }

  for (BP_Range::iterator it = cont_rng->begin(); it <= cont_rng->end(); it++,_cont_p++) {
    int i = (int) *it;
    float *fp = (float *)_cont_p;

    if (fp == NULL) {
      warning("ObservationMatrix::readPFloats: memory error trying to read %i'th feature",i);
      return 0;
    }

    *fp = _contFea.ptr[i];
  }

  return 1;
}

/* read ascii floats
 * n_floats: number of floats to read
 * f: file to read from
 * cont_rng: range of features to use
 * returns 1 if successful, else 0
 */

bool
ObservationMatrix::readAscFloats(unsigned n_floats, 
				 FILE *f, 
				 BP_Range *cont_rng) {
  
  if (_cont_p == NULL) {
    warning("ObservationMatrix::readAscFloats: Data buffer is NULL\n");
    return 0;
  }

  if (_contFea.ptr == NULL) {
    warning("ObservationMatrix::readAscFloats: Temporary buffer is NULL\n");
    return 0;
  }
	
  float *fea_p = _contFea.ptr;
  
  for (unsigned n = 0; n < n_floats; n++,fea_p++)
    if (fscanf(f,"%e",fea_p) != 1) {
      warning("ObservationMatrix::readAscFloats: couldn't read %i'th item in frame\n",n);
      return 0;
    }
  
  for (BP_Range::iterator it = cont_rng->begin(); it <= cont_rng->end(); it++,_cont_p++) {
    int i = (int) *it;
    float *fp = (float *)_cont_p;

    if (fp == NULL) {
      warning("ObservationMatrix::readAscFloats: memory error trying to read %i'th feature",i);
      return 0;
    }

    *fp = _contFea.ptr[i];
  }
  return 1;
}

/* read integers from pfile stream
 * str: pfile stream
 * disc_rng: range of features to use
 * returns 1 if successful, else 0
 */

bool
ObservationMatrix::readPInts(InFtrLabStream_PFile *str, BP_Range *disc_rng) {


  if (_disc_p == NULL) {
    warning("ObservationMatrix::readPInts: Data buffer is NULL");
    return 0;
  }
  
  if (_discFea.ptr == NULL) {
    warning("ObservationMatrix::readPInts: Temporary buffer is NULL");
    return 0;
  }
  
  if (str->read_labs((size_t)1,(UInt32 *)_discFea.ptr) != 1) {
    warning("ObservationMatrix::readPInts: read failed\n");
    return 0;
  }
  for (BP_Range::iterator it = disc_rng->begin(); it <= disc_rng->end(); it++,_disc_p++) {
    int i = (int)*it;

    if ( _disc_p == NULL) {
      warning("ObservationMatrix::readPInts: memory error trying to read %i'the feature",i);
      return 0;
    }

    (Int32 *) _disc_p = _discFea.ptr[i];
  }
  return 1;
}

/* read binary integers 
 * n_ints: number of features to read
 * f: stream to read from
 * disc_rng: range of features to use
 * bswap: bytes will be swapped if true
 * returns 1 if successful, else 0
 */

bool
ObservationMatrix::readBinInts(unsigned n_ints, FILE *f, BP_Range *disc_rng, 
			       bool bswap) {

  unsigned n_read;  

  if (_disc_p == NULL) {
    warning("ObservationMatrix::readBinInts: Data buffer is NULL");
    return 0;
  }
  if (_discFea.ptr == NULL) {
    warning("ObservationMatrix::readBinInts: Temporary buffer is NULL");
    return 0;
  }
  
  n_read = fread((Int32 *)_discFea.ptr,sizeof(Int32),n_ints,f);

  if (n_read != n_ints) {
    warning("ObservationMatrix::readBinInts: only read %i ints, expected %i\n",
	    n_read,n_ints);
    return 0;
  }

  
  for (BP_Range::iterator it = disc_rng->begin(); it <= disc_rng->end(); it++,_disc_p++) {
  
  int i = (int)*it;

    Int32 *ip = (Int32 *)_disc_p;

    if (ip == NULL) {
      warning("ObservationMatrix:::readBinInts: memory error trying to read %i'th feature",i);
      return 0;
    }

    if (bswap) 
      *ip  = swapb_i32_i32(_discFea.ptr[i]);
    else 
      *ip = _discFea.ptr[i];
  }
  return 1;
}

/* read ascii integers 
 * n_ints: number of features to read
 * f: stream to read from
 * disc_rng: range of features to use
 * returns 1 if successful, else 0
 */

bool
ObservationMatrix::readAscInts(unsigned n_ints, FILE *f, BP_Range *disc_rng) {


  if (_disc_p == NULL) {
    warning("ObservationMatrix::readAscInts: Data buffer is NULL");
    return 0;
  }

   if (_discFea.ptr == NULL) {
    warning("ObservationMatrix::readAscInts: Temporary buffer is NULL");
    return 0;
  }

  Int32 *dp = _discFea.ptr;
 
  for (unsigned a = 0; a < n_ints; a++,dp++) {
    if (fscanf(f,"%i",dp) != 1) {
      warning("ObservationMatrix::readAscInts: couldn't read %i'th item in frame\n",a);
      return 0;
    }
  }
  
  for (BP_Range::iterator it = disc_rng->begin(); it <= disc_rng->end(); it++,_disc_p++) {

    int i = (int)*it;

    Int32 *ip = (Int32 *)_disc_p;

    if (ip == NULL) {
      warning("ObservationMatrix:::readAscInts: memory error trying to read %i'th feature",i);
      return 0;
    } 
    *ip = _discFea.ptr[i];
  }
  return 1;
}


/* read single input frame  
 * frameno: number of frame to be read
 */

void
ObservationMatrix::readFrame(size_t frameno) {

  
  StreamInfo *f;
  unsigned num_floats,num_ints;
  int i;

  
  for (i = 0; i < _numStreams; i++) {
    
    f = _inStreams[i];	
    
    if (f == NULL)
      error("ObservationMatrix::readFrame: Cannot access info for stream %i\n",i);
    
    num_floats = f->nFloats;
    
    assert(num_floats <= _maxContinuous);
    
    num_ints = f->nInts;

    assert(num_ints <= _maxDiscrete);
	
    switch(f->dataFormat) {

    case PFILE:
      if (num_floats > 0)
	if (!(readPFloats(f->pfile_istr,f->cont_rng)))
	  error("ObservationMatrix::readFrame: can't read floats for frame %li from stream %i\n",frameno,i);
      
      if (num_ints > 0)
	if (!(readPInts(f->pfile_istr,f->disc_rng)))
	  error("ObservationMatrix::readFrame: can't read ints for frame %li from stream %i\n",frameno,i);
      break;

    case RAWASC:
      if (num_floats > 0) 
	if (!(readAscFloats(num_floats,f->curDataFile,f->cont_rng)))
	  error("ObservationMatrix::readFrame: can't read floats for frame %li from stream %i\n",frameno,i);
      
      if (num_ints > 0)
	 if (!(readAscInts(num_ints,f->curDataFile,f->disc_rng)))
	   error("ObservationMatrix::readFrame: can't read ints for frame %li from stream %i\n",frameno,i);
      break;

    case RAWBIN:
    case HTK:
      if (num_floats > 0) 
	if (!(readBinFloats(num_floats,f->curDataFile,f->cont_rng,f->bswap)))
	  error("ObservationMatrix::readFrame: can't read floats for frame %li from stream %i \n",frameno,i);
      
      if (num_ints > 0) 
	if (!(readBinInts(num_ints,f->curDataFile,f->disc_rng,f->bswap)))
	  error("ObservationMatrix::readFrame: can't read ints for frame %li from stream %i \n",frameno,i);
      break;
    default:
      error("Unsupported data file format\n");	
    }
  }
}


/* reset pointers to beginning of data buffer */

void
ObservationMatrix::reset() {

 _cont_p = features.ptr;
 _disc_p = features.ptr + numContinuous;
}


/* resize data buffer
 * n_frames: new number of frames
 */

void
ObservationMatrix::resize(size_t n_frames) {

   features.growIfNeeded(n_frames * stride);
   _bufSize = n_frames;
}

/* print frame from data buffer
 * stream: output stream
 * frameno: frame number
 */

void
ObservationMatrix::printFrame(FILE *stream, size_t frameno) {
  
  unsigned f;

  Data32 *p = features.ptr + stride*frameno;

  for (f = 0; f < numContinuous; f++,p++) {
	float *fp = (float *)p;
	fprintf(stream,"%e ",*fp);
  }
  for (f = 0; f < numDiscrete; f++,p++)
     fprintf(stream,"%i ",(Int32)*p);

  fprintf(stream,"\n");

}

// get cont. feature number 'n' in current frame
	
float *
ObservationMatrix::getContFea(unsigned short n) {


  assert(n < numFeatures);

  if (n > (numContinuous-1))
    warning("ObservationMatrix::getFeature: feature %i has type int but requested as float\n");


  return (float *)(_cont_p + n);
}

// get disc. feature number 'n' in current frame


Int32 *
ObservationMatrix::getDiscFea(unsigned short n) {

  assert (n < numFeatures);

  if (n > (numDiscrete-1))
    warning("ObservationMatrix::getFeature: feature %i has type int but requested as float\n");

  return (Int32 *)(_disc_p + n);
}

/* get info for segment 'sentno' from pfile stream 'f' */

size_t
ObservationMatrix::openPFile(StreamInfo *f, size_t sentno) {

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


size_t
ObservationMatrix::openBinaryFile(StreamInfo *f, size_t sentno) {

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

size_t
ObservationMatrix::openAsciiFile(StreamInfo *f,size_t sentno) {

  char *fname = f->dataNames[sentno];
  size_t n_samples = 0;
  char ch;

  if (fname == NULL) {
    warning("ObservationMatrix::openAsciiFile: Filename is NULL for segment %li\n",sentno);
    return 0;
  }

  if ((f->curDataFile = fopen(fname,"r")) == NULL) {
    warning("ObservationMatrix::openAsciiFile: Can't open '%s' for input\n",fname);
    return 0;
  }

  /* for ascii, newline is record delimiter - additional or missing nl's will cause error messages */

  while ((ch = fgetc(f->curDataFile)) != EOF) {
    if (ch == '\n')
      n_samples++;
  }
 rewind(f->curDataFile);
 f->curNumFrames = n_samples;
 
 return n_samples;
}

/* open htk file for segment 'sentno' */

size_t
ObservationMatrix::openHTKFile(StreamInfo *f, size_t sentno) {

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

  if (pk <= WAVEFORM || pk > ANON) {
    warning("Undefined parameter kind for HTK feature file: %i\n",pk);
    return 0;
  }

  int n_fea;

  // parameter kind DISCRETE = all discrete features

  if (parm_kind == DISCRETE) {
    n_fea = samp_size / sizeof(int) ;
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

/* close individual data files */

void
ObservationMatrix::closeDataFiles() {

  for (int i = 0; i < _numStreams; i++)
    if (_inStreams[i]->dataFormat != PFILE)
      fclose(_inStreams[i]->curDataFile);
}



void
ObservationMatrix::printSegmentInfo() {

  printf("Processing segment # %d. Number of frames = %d.\n",
	segmentNumber,numFrames);
}
