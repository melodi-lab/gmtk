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


/* constructor for data buffer (= feature matrix)
 * n_frames: max number of frames in buffer
 * n_floats: number of continuous features
 * n_ints: number of discrete features
 * max_floats: max number of floats per frame in any of the input streams
 * max_ints: max number of ints per frame in any of the input streams
 */


ObservationMatrix::ObservationMatrix(size_t n_frames, 
				     unsigned n_floats, 
				     unsigned n_ints,
				     unsigned max_floats,
				     unsigned max_ints) {
				     

  _numFrames = n_frames;
  _numContinuous = n_floats;
  _numDiscrete = n_ints;
  _numFeatures = _numContinuous + _numDiscrete;
  _maxContinuous = max_floats;
  _maxDiscrete = max_ints;
  
  _stride = _numFeatures; 
  
  features.resize(_numFrames * _numFeatures); // actual observation buffer
  
  cont_fea.resize(_maxContinuous); // temporary buffers for 1 frame of input
  disc_fea.resize(_maxDiscrete);   

  _bufSize = _numFrames;

  // structure of observation matrix: in each frame, all continuous
  // features come first, followed by all discrete features
  
  _cont_p = features.ptr;  // pointer to continuous block 
  _disc_p = features.ptr + _numContinuous; // pointer to discrete block
}


ObservationMatrix::~ObservationMatrix() {
  
  ;
}


/* prepares data buffer for next frame to be read in
 * (advances pointers to current location in buffer)
 */

void
ObservationMatrix::next() {
  
  _cont_p += _numDiscrete;
  _disc_p += _numContinuous;
}

/* read binary floats into data buffer
 * n_floats: number of floats to be read
 * f: file to read from
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

  if (cont_fea.ptr == NULL) {
    warning("ObservationMatrix::readBinFloats: Temporary buffer is NULL\n");
    return 0;
  }

  // read complete frame
  
  n_read = fread((float *)cont_fea.ptr,sizeof(float),n_floats,f);
  
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
      warning("OvservationMatrix::readBinFloats: memory error reading %i'th feature",i);
      return 0;
    }

    if (bswap) 
      swapb_vi32_vi32(1,(const int *)&cont_fea.ptr[i],(int *)fp);
    else
      *fp = cont_fea.ptr[i];
    
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
  
  if (cont_fea.ptr == NULL) {
    warning("ObservationMatrix::readPFLoats: Temporary buffer is NULL\n");
    return 0;
  }
  
  if ((n_read = str->read_ftrslabs((size_t)1,(float *)cont_fea.ptr,NULL)) != 1) {
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

    *fp = cont_fea.ptr[i];
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

  if (cont_fea.ptr == NULL) {
    warning("ObservationMatrix::readAscFloats: Temporary buffer is NULL\n");
    return 0;
  }
	
  float *fea_p = cont_fea.ptr;
  
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

    *fp = cont_fea.ptr[i];
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
  
  if (disc_fea.ptr == NULL) {
    warning("ObservationMatrix::readPInts: Temporary buffer is NULL");
    return 0;
  }
  
  if (str->read_labs((size_t)1,(UInt32 *)disc_fea.ptr) != 1) {
    warning("ObservationMatrix::readPInts: read failed\n");
    return 0;
  }
  for (BP_Range::iterator it = disc_rng->begin(); it <= disc_rng->end(); it++,_disc_p++) {
    int i = (int)*it;

    if ( _disc_p == NULL) {
      warning("ObservationMatrix::readPInts: memory error trying to read %i'the feature",i);
      return 0;
    }

    (Int32 *) _disc_p = disc_fea.ptr[i];
  }
  return 1;
}

/* read binary integers 
 * n_ints: number of features to read
 * f: file to read from
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
  if (disc_fea.ptr == NULL) {
    warning("ObservationMatrix::readBinInts: Temporary buffer is NULL");
    return 0;
  }
  
  n_read = fread((Int32 *)disc_fea.ptr,sizeof(Int32),n_ints,f);

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
      *ip  = swapb_i32_i32(disc_fea.ptr[i]);
    else 
      *ip = disc_fea.ptr[i];
  }
  return 1;
}

/* read ascii integers 
 * n_ints: number of features to read
 * f: file to read from
 * disc_rng: range of features to use
 * returns 1 if successful, else 0
 */

bool
ObservationMatrix::readAscInts(unsigned n_ints, FILE *f, BP_Range *disc_rng) {


  if (_disc_p == NULL) {
    warning("ObservationMatrix::readAscInts: Data buffer is NULL");
    return 0;
  }

   if (disc_fea.ptr == NULL) {
    warning("ObservationMatrix::readAscInts: Temporary buffer is NULL");
    return 0;
  }

  Int32 *dp = disc_fea.ptr;
 
  for (int a = 0; a < n_ints; a++,dp++) {
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
    *ip = disc_fea.ptr[i];
  }
  return 1;
}


/* read single input frame, possibly using input from different files  
 * frameno: number of frame in file to read
 * fds: array of file descriptions
 * numFiles: number of files (length of fds)
 */

void
ObservationMatrix::readFrame(size_t frameno, 
			     FileDescription **fds, 
                             int numFiles) {

  FileDescription *f;
  unsigned num_floats,num_ints;
  int i;

  
  for (i = 0; i < numFiles; i++) {
    
    f = fds[i];	
    
    if (f == NULL)
      error("ObservationMatrix::readFrame: Cannot access file descriptor %i\n",i);

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
  next(); 
}


/* reset pointers to beginning of data buffer */

void
ObservationMatrix::reset() {

 _cont_p = features.ptr;
 _disc_p = features.ptr + _numContinuous;
}


/* resize data buffer
 * n_frames: new number of frames
 */

void
ObservationMatrix::resize(size_t n_frames) {

   features.growIfNeeded(n_frames * _stride);
   _bufSize = n_frames;
}

/* print frame from data buffer
 * stream: output stream
 * frameno: frame number
 */

void
ObservationMatrix::printFrame(FILE *stream, size_t frameno) {
  
  unsigned f;

  Data32 *p = features.ptr + _stride*frameno;

  for (f = 0; f < _numContinuous; f++,p++) {
	float *fp = (float *)p;
	fprintf(stream,"%e ",*fp);
  }
  for (f = 0; f < _numDiscrete; f++,p++)
     fprintf(stream,"%i ",(Int32)*p);

  fprintf(stream,"\n");

}
  
	
float *
ObservationMatrix::getContFea(unsigned short n) {

  if (n > (_numContinuous-1))
    warning("ObservationMatrix::getFeature: feature %i has type int but requested as float\n");

  return (float *)(_cont_p + n);
}

Int32 *
ObservationMatrix::getDiscFea(unsigned short n) {

  if (n > (_numDiscrete-1))
    warning("ObservationMatrix::getFeature: feature %i has type int but requested as float\n");

  return (Int32 *)(_disc_p + n);
}

/* set pointers to beginning of frame f */

void
ObservationMatrix::gotoFrame(size_t f) {

  _cont_p = features.ptr + (f*_stride);
  _disc_p = _cont_p + _numDiscrete;

}

