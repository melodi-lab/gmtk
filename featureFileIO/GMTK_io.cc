//static char *rcsid = "$Id$";
/*
 *
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
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include "GMTK_utils.h"
#include "GMTK_io.h"

/* create data input stream  
 * n_files: number of input files
 * fof: list of files-of-filenames
 * cont_range_str: array of continuous range specification strings
 * disc_range_str: array of discrete range specification strings
 * n_floats: array of numbers of floats for each file
 * n_ints: array of numbers of ints for each file
 * formats: array of format specification for each file
 * swapflag: array of boolean flags indicating whether byteswapping is required
 */

DataInStream::DataInStream(int n_files, const char **fof_names,
                           const char **cont_range_str,
                           const char **disc_range_str,
			   unsigned *n_floats, unsigned *n_ints, 
			   unsigned *formats, bool *swapflag) {

  _numFiles = n_files;

  assert(_numFiles > 0);


  if (fof_names == NULL)
    error("ERROR DataInStream: need list of filenames\n");

  if (formats == NULL) 
    error("ERROR DataInStream: file formats must be specified\n");

  if (n_floats == NULL)
    error("DataInStream: number of floats must be specified\n");

  if (n_ints == NULL)
    error("DataInStream: number of ints must be specified\n");

  inp_buf = NULL;

  _totalFloats = 0;
  _totalInts = 0;

  fd = new FileDescription*[_numFiles];

  for (int i = 0; i < _numFiles; i++) {
      fd[i] = new FileDescription(fof_names[i],cont_range_str[i],
		disc_range_str[i], &n_floats[i],&n_ints[i],
                &formats[i],&swapflag[i],i);


      if (i > 0 ) {

          size_t a = fd[i-1]->fofSize; 
	  size_t b = fd[i]->fofSize;

        if (a != b) {
   	 warning("WARNING: DataInStream: different number of files in '%s' (%li) and '%s' (%li) - will only read minimum number\n",fof_names[i-1],a,
	fof_names[i], b);

	a < b ? fd[i]->fofSize = a : fd[i-1]->fofSize = b; 
       }
     }

    _totalFloats += fd[i]->nFloatsUsed;
    _totalInts += fd[i]->nIntsUsed;
  }
  _numSegs = fd[0]->fofSize;

  inp_buf = new ObservationMatrix(MAXFRAMES,_totalFloats,_totalInts);
}


DataInStream::~DataInStream() {

  for (int i = 0; i < _numFiles; i++) 
    delete fd[i];

  delete [] fd;
  delete inp_buf;

}


/* read input for single segment (utterance)
 * sent_no: segment number
*/


size_t
DataInStream::readFile(size_t sent_no) {

  int i;
  size_t n,n_samps;
  FileDescription *f;
  char *tmp_name;

  _curPos = sent_no;

  size_t max_sents =fd[0]->fofSize; 

  if (sent_no < 0 || sent_no >= max_sents) {
      error("GMTK_io::readFile: Sentence number (%li) outside range of 0 - %li\n",
	sent_no,max_sents);
  }


  for (i = 0; i < _numFiles; i++) {

    f = fd[i];     
         	
    if (f == NULL) 
       error("GMTK_io::readFile: file descriptor for file %i is zero\n",i);

    if (f->dataFormat != PFILE) {
	
    if (f->dataNames == NULL)
      error("List of data files for stream %i is NULL\n",i);

     tmp_name = f->dataNames[sent_no];

      if (tmp_name == NULL)
	error("GMTK_io::readFile: File %i does not exist in list %s\n",
            sent_no,f->fofName);
    }
    else {
      if (f->curDataFile == NULL) 
	error("GMTK_io::readFile: PFile data stream is NULL for file %s\n",
        	f->fofName);

      tmp_name = f->fofName;
    }
    
    switch(f->dataFormat) {
      
    case RAWBIN:
     n_samps = openBinaryFile(f,sent_no);
      break;
    case RAWASC:
      n_samps = openAsciiFile(f,sent_no);
      break;
    case HTK:
      n_samps = openHTKFile(f,sent_no);
      break;
    case PFILE:
      n_samps = openPFile(f,sent_no);
      break;
    default:
      error("GMTK_io::readFile: Invalid file format specified for file %i\n",i);
    }

    if (n_samps == 0) 
      error("GMTK_io::readFile: failure to read file '%s'\n",tmp_name);

    if (i > 0 && n_samps != fd[i-1]->curDataSize) {
	error("GMTK_io::readFile: Number of samples for sentence %i don't match for streams %s and %s (%li vs. %li)\n",sent_no,fd[i-1]->fofName,f->fofName,fd[i-1]->curDataSize,n_samps);
   }
 }

    // resize buffer if necessary
    
    if (n_samps > inp_buf->getBufSize()) 
      inp_buf->resize(n_samps*2);
    
    inp_buf->reset(); 
    
    for (n = 0; n < n_samps; n++) {
      inp_buf->readFrame(n,fd,_numFiles);
   }

   inp_buf->setSegNo(sent_no);
   inp_buf->setNumFrames(n);

   closeDataFiles();
    
   return n; // num frames read
}


/* initialize segment in pfile
 * f: file descriptor
 * sentno: segment numer 
 */

size_t
DataInStream::openPFile(FileDescription *f, size_t sentno) {

  if (f->pfile_istr == NULL) {
    error("GMTK_io::openPFile: stream is NULL");
    return 0;
  }
  	
  if (f->pfile_istr->set_pos(sentno,0) == SEGID_BAD) {
    warning("GMTK_io::openPFile: Can't skip to sent %li frame 0",
		 sentno);
    return 0;
  }

  f->curDataSize = f->pfile_istr->num_frames(sentno);

  return f->curDataSize;
}

/* open binary file
 * f: file descriptor
 * sentno: segment number
 */

size_t
DataInStream::openBinaryFile(FileDescription *f, size_t sentno) {

  char *fname = f->dataNames[sentno];
  int nfloats = f->nFloats;
  int nints = f->nInts;
  size_t fsize;	

  if (fname == NULL) {
    warning("GMTK_io::openBinaryFile: Filename is NULL for segment %li\n",sentno);
    return 0;
  }

  if ((f->curDataFile = fopen(fname,"rb")) == NULL) {
    error("GMTK_io::openBinaryFile: Can't open '%s' for input\n",
		 fname);
  }

  //sanity check on number of bytes

  if (fseek(f->curDataFile,0L,SEEK_END) == -1) {
    error("GMTK_io::openBinaryFile: Can't skip to end of file",
		 fname);
  }

  fsize = ftell(f->curDataFile);

  rewind(f->curDataFile);

  int rec_size = nfloats * sizeof(float) + nints * sizeof(int);
  
  if ((fsize % rec_size) > 0) 
    error("GMTK_io::openBinaryFile: odd number of bytes in file %s\n",fname);


  int n_samples = fsize / rec_size;

  f->curDataSize = n_samples;

  return n_samples;
}

/* open ascii file
 * f: file descriptor
 * sentno: segment number
 */

size_t
DataInStream::openAsciiFile(FileDescription *f,size_t sentno) {

  char *fname = f->dataNames[sentno];
  size_t n_samples = 0;
  char ch;

  if (fname == NULL) {
    warning("GMTK_io::openAsciiFile: Filename is NULL for segment %li\n",sentno);
    return 0;
  }

  if ((f->curDataFile = fopen(fname,"r")) == NULL) {
    warning("GMTK_io::openAsciiFile: Can't open '%s' for input\n",fname);
    return 0;
  }

  // for ascii, newline is record delimiter - additional or missing nl's will cause error messages

  while ((ch = fgetc(f->curDataFile)) != EOF) {
    if (ch == '\n')
      n_samples++;
  }
  rewind(f->curDataFile);
  f->curDataSize = n_samples;	

  return n_samples;
}


/* open HTK file 
 * f: file descriptor
 * sentno: segment number
 */

size_t
DataInStream::openHTKFile(FileDescription *f, size_t sentno) {

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
    warning("GMTK_io::openHTKFile: Filename is NULL for segment %li\n",sentno);
    return 0;
  }

  if ((f->curDataFile = fopen(fname,"rb")) == NULL) {
    warning("GMTK_io::openHTKFile: Can't open '%s' for input\n",fname);
    return 0;
  }

  if (fread(&tmp1,sizeof(Int32),1,f->curDataFile) != 1) {
    warning("GMTK_io::openHTKFile: Can't read number of samples\n");
    return 0;
  }

  if (fread((short *)&tmp2,sizeof(Int32),1,f->curDataFile) != 1) {
    warning("GMTK_io::openHTKFile: Can't read sample period\n");
    return 0;
  }

  if (fread((short *)&stmp1,sizeof(short),1,f->curDataFile) != 1) {
    warning("GMTK_io::openHTKFile: Can't read sample size\n");
    return 0;
  }

  if (fread(&stmp2,sizeof(short),1,f->curDataFile) != 1) {
    warning("GMTK_io::openHTKFile: Can't read parm kind\n");
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

  if (n_samples <= 0)  {
     warning("GMTK_io::openHTKFile: number of samples is %i\n",n_samples);   
     return 0;
  }

  if (samp_period <= 0 || samp_period > 1000000)  {
     warning("GMTK_io::openHTKFile: sample period is %i - must be between 0 and 1000000\n", samp_period); 
     return 0;
  }

  if (samp_size <= 0 || samp_size > 5000) {
     warning("GMTK_io::openHTKFile: sample size is %i - must be between 0 and 5000\n",samp_size);
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
      warning("GMTK_io::openHTKFile:  Number of features in file (%i) does not match number of ints specified (%i)\n", n_fea,nints);
      return 0;
    }
  }

  // otherwise all continuous features

  else {
    n_fea = samp_size / sizeof(float);
    if (n_fea != nfloats) {
      warning("GMTK_io::openHTKFile:  Number of features in file (%i) does not match number of floats specified (%i)\n", n_fea,nfloats);
      return 0;
    }
  }

  f->curDataSize = n_samples;

  return n_samples;
}

/* close individual data files */

void
DataInStream::closeDataFiles() {

  for (int i = 0; i < _numFiles; i++) 
   if (fd[i]->dataFormat != PFILE) 
      fclose(fd[i]->curDataFile);

}


/* get file names associated with individual utterance
 * sentno: segment number
 * buf: buffer for filename pointers 
 */

void
DataInStream::getFileNames(size_t sentno,char **buf) {

  for (int i = 0; i < _numFiles; i++) {
    if (sentno >= fd[i]->fofSize)
      error("DataInStream::getFileNames: sent no %i exceeds number of file name\n",
	    sentno);
    if (buf[i] == NULL)
      error("DataInStream::getFileNames: memory error\n");
    strcpy(buf[i],fd[i]->dataNames[sentno]);
  }
}     	            	


/* get numbers of floats associated with different input streams
 * buf: buffer to store numbers of floats
 */

void
DataInStream::getNumFloats(unsigned *buf) {


  unsigned *p = buf;

  for (int i = 0; i < _numFiles; i++,p++)
        *p = fd[i]->nFloats;
}



/* get numbers of floats associated with different input streams
 * buf: buffer to store numbers of ints
 */

void
DataInStream::getNumInts(unsigned *buf) {


  unsigned *p = buf;

  for (int i = 0; i < _numFiles; i++,p++)
        *p = fd[i]->nInts;

}


/* get pointer to observation matrix */

ObservationMatrix *
DataInStream::getObs() {

  return inp_buf;

}		

/* output stream constructor */


DataOutStream::DataOutStream(int n_files, 
                             unsigned *n_floats, 
                             unsigned *n_ints,
                             unsigned *formats, 
                             bool *swapflags,
                             int *htkPK,
                             unsigned *samp_period,
			     char **pfile_names) {


   assert(n_files > 0);

   _numFiles = n_files;

   _bufSize = MAXFRAMES;

   // obligatory items

   if (n_floats == NULL) 
    error("GMTK_io::DataOutStream: num floats needs to be specified\n");
  
   if (n_ints == NULL)
    error("GMTK_io::DataOutStream: num ints needs to be specified\n");
      
   if (formats == NULL)
    error("GMTK_io::DataOutStream: output data formats need to be specified\n");

   outNames = new char*[_numFiles];
   outFiles = new FILE*[_numFiles];
   dataFormats = new unsigned[_numFiles];
   bswap = new bool[_numFiles];
   nFloats = new unsigned[_numFiles];
   nInts = new unsigned[_numFiles];
   pfile_ostr = new OutFtrLabStream_PFile*[_numFiles];
   sampRate = new unsigned[_numFiles];
   parKind = new short[_numFiles];

   unsigned sp;
   bool sw;
   int pk;
   char *name;


   // null arrays explicitly to avoid portability problems

   for (int i = 0; i < _numFiles; i++) {

     if (samp_period == NULL)
        sp = NULL;
     else
       sp = samp_period[i];

     if (htkPK == NULL)
       pk = NULL;
     else
       pk = htkPK[i];

     if (swapflags == NULL)
       sw = NULL;
     else
       sw = swapflags[i];

     if (pfile_names == NULL)
       name = NULL;
     else
       name = pfile_names[i];

     pfile_ostr[i] = NULL;

     
       
     initOutStream(&n_ints[i],&n_floats[i],&formats[i],
		   &sw,&pk,&sp,&name,i);
   }
}



void
DataOutStream::initOutStream(unsigned *n_ints, unsigned *n_floats, unsigned *format,
			     bool *swapflag,int *pk,unsigned *samp_period,
			     char **pfile_name,unsigned num) {
     


     if (format == NULL) {
       warning("No output file format for file stream %i - writing headerless binary\n",num); 
       dataFormats[num] = RAWBIN;
     }
     else 
       dataFormats[num] = *format;

     _maxInts = 0;
     _maxFloats = 0;

     ftr_buf = NULL;
     lab_buf = NULL;

     if (n_ints == NULL)
       error("DataOutStream::initOutStream: Number of ints for output stream %i needs to be specified\n",num);
     
     nInts[num] = *n_ints;

     if (nInts[num] > _maxInts)
       _maxInts = nInts[num];

     if (n_floats == NULL)
       error("DataOutStream:initOutStream: Number of floats for output stream %i needs to be specified\n",num);
     
     nFloats[num] = *n_floats;

     if (nFloats[num] > _maxFloats)
       _maxFloats = nFloats[num];

     
     sampRate[num] = *samp_period;
     
     parKind[num] = *pk;
     
     if (swapflag == NULL)
       bswap[num] = false;  // default: no byte-swapping on output
     else
       bswap[num] = *swapflag;

     if (dataFormats[num] == PFILE) {

       // allocate temporary feature buffers for pfile stream
     
       ftr_buf = new float[_bufSize*_maxFloats];
       lab_buf = new UInt32[_bufSize*_maxInts];

       if (pfile_name == NULL)
	 error("GMTK_io::DataOutStream: pfile name must be specified\n");

       if ((outFiles[num] = fopen(*pfile_name,"wb")) == NULL)
	 error("GMTK_io::DataOutStream: can't create pfile '%s'\n",pfile_name);

       pfile_ostr[num] = new OutFtrLabStream_PFile(0,*pfile_name,outFiles[num],
						 nFloats[num],nInts[num],1,
						 bswap[num]);

     }
     
     if (nFloats[num] > 0 && nInts[num] > 0 && dataFormats[num] != PFILE)
       error("GMTK_io::DataOutStream: Can't write file with both floats and ints unless it is a pfile\n");

}


void
DataOutStream::initFiles(char **file_names,char **out_ext,size_t n_frames) {


  _numFrames = n_frames;

  if (file_names == NULL)
    error("DataOutStream::initFiles: File names need to be specified\n");

  for (int i = 0; i < _numFiles; i++) {

    if (dataFormats[i] != PFILE) {

      char *ip = file_names[i];
      char *xp = out_ext[i];
      if (ip == NULL)
	error("DataOutStream::initFiles: output file name for stream %i needs to be specified\n",i);
      if (xp != NULL) {
	outNames[i] = new char[strlen(ip)+strlen(xp)+1];
	strcpy(outNames[i],ip);

	for (int c = strlen(outNames[i]); c >= 0; c--)
	  if (outNames[i][c] == '.')
	    outNames[i][c] = '\0';
	strcat(outNames[i],".");
	strcat(outNames[i],xp);
      }
      else {
	outNames[i] = new char[strlen(ip+1)];
	char *cp = outNames[i];
	strcpy(cp,ip);
      }
      
      if ((outFiles[i] = fopen(outNames[i],"w")) == NULL) 
	error("DataOutStream::initFiles: Cannot open '%s' for output\n",outNames[i]);
      
      switch(dataFormats[i]) {
      case HTK:
	writeHTKHeader(outFiles[i],nFloats[i],nInts[i],sampRate[i],
		       _numFrames,parKind[i]);
	break;
      case RAWBIN:
      case RAWASC:
      case PFILE:
	break;
      default:
	error("DataOutStream::initFiles: illegal file format for output stream %s\n",
	      outNames[i]);
      }
    }
  }
}

void
DataOutStream::closeFiles() {

  for (int i = 0;i < _numFiles; i++)
    if (dataFormats[i] != PFILE) {
      fclose(outFiles[i]);
      delete [] outNames[i];
    }

}

       
void   
DataOutStream::writeHTKHeader(FILE *f, 
                             unsigned n_floats,
                             unsigned n_ints,
                             unsigned samp_rate,
                             size_t n_frames,
                             short pk) {

  
  if (f == NULL) 
    error("DataOutStream::writeHTKHeader: output file is NULL\n");

  // checking has been done but just in case

  if (samp_rate == NULL) {
    warning("DataOutStream::writeHTKHEader: sampling rate is NULL");
    warning(" - assuming 10ms frame rate");
    samp_rate = 100000;
  }
  

  if (pk == NULL) {
    warning("DataOutStream::writeHTKHeader: no parameter kind was specified ");

    if (n_floats > 0) {
      warning("- assuming MFCC\n");
      pk = MFCC;
    }
    else if (n_ints > 0) {
      warning("-assuming DISCRETE\n");
      pk = DISCRETE;
    }
  }
  else {

    // we don't write waveforms, only feature files

   if (pk < LPC || pk > ANON)
      error("DataOutStream::writeHTKHeader: Invalid HTK parameter kind\n");

   // discrete data can be DISCRETE, USER or ANON

   if (n_ints > 0 && pk < USER)
      error("DataOutStream::writeHTKHeader: data is discrete but parameter kind is not\n");
   if (n_floats > 0 && pk == DISCRETE)
      error("DataOutStream::writeHTKHeader: parameter kind is discrete but data is not\n");

  }

  if (fwrite((int *)&n_frames,sizeof(int),1,f) != 1) 
    error("DataOutStream::writeHTKHeader: error writing to file\n");

  if (fwrite((int *)&samp_rate,sizeof(int),1,f) != 1)
        error("DataOutStream::writeHTKHeader: error writing to file\n");

  short samp_size;

  if (n_floats > 0)
	samp_size = n_floats * sizeof(float);
  else if (n_ints > 0)
	samp_size = n_ints * sizeof(int);

  if (fwrite((short *)&samp_size,sizeof(short),1,f) != 1)
      error("DataOutStream::writeHTKHeader: error writing to file\n");

  if (fwrite((short *)&pk, sizeof(short),1,f) != 1)
     error("DataOutStream::writeHTKHeader: error writing to file\n");
}
	

/* write feature output
 * sent_no: utterance number
 * crng_str: array of continuous range strings
 * drng_str: array of discrete range strings
 * frng_str: frame range string
 * obs: observation matrix 
 */      	
	

void
DataOutStream::writeData(size_t sent_no,
			 char **crng_str, 
			 char **drng_str,
                         char *frng_str, 
			 ObservationMatrix *obs) {

     BP_Range *cr, *dr, *fr;
     unsigned num_floats, num_ints;
     char *s;

     for (int i = 0; i < _numFiles; i++) {
       if (crng_str == NULL || crng_str[i] == NULL)
	 s = "all"; // default for range strings: all items
       else
	 s = crng_str[i];

       cr = new BP_Range(s,0,nFloats[i]);

       if (drng_str == NULL || drng_str[i] == NULL)
	 s = "all";
       else
	 s  = drng_str[i];

       dr = new BP_Range(s,0,nInts[i]);

       if (frng_str == NULL)
	 s = "all";
       else
	 s = frng_str;

       fr = new BP_Range(s,0,obs->getNumFrames());

       num_floats = cr->length();
       num_ints = dr->length();
 
       if (num_floats > 0 && num_ints > 0 && dataFormats[i] != PFILE)
	 error("DataOutStream::writeData: stream %i: cannot write both discrete and continuous features unless it is a pfile\n");
       
       switch(dataFormats[i]) {
       case PFILE:
	 writePFileFeatures(pfile_ostr[i],obs,cr,dr,fr,sent_no);
	 break;
       case RAWASC:
	 writeAscFeatures(outFiles[i],obs,cr,dr,fr);
	 break;
       case HTK:
       case RAWBIN:
	 writeBinFeatures(outFiles[i],obs,cr,dr,fr);
	 break;
       default: 
	 error("DataOutStream::writeData: undefined file format\n");
       }
     }

     delete cr;
     delete dr;
     delete fr;
}      


/* write binary output
 * fp: output file
 * o: observation matrix
 * cr: continuous feature range
 * dr: discrete feature range
 * fr: frame range
 */

void
DataOutStream::writeBinFeatures(FILE *fp, 
				ObservationMatrix *o, 
				BP_Range *cr, 
				BP_Range *dr,
                                BP_Range *fr) {

  // for all frames in range

  for (BP_Range::iterator fit  = fr->begin(); fit <= fr->max(); fit++) {

    o->gotoFrame(*fit);
    
    // print continuous features
    
    for (BP_Range::iterator cit = cr->begin(); cit <= cr->max(); cit++) {
      
      float *f = o->getContFea(*cit);
      
      if (fwrite((float *)f,sizeof(float),1,fp) != 1)
	error("DataOutStream::writeBinFeatures: couldn't write %i'th feature in frame \n",*cit,*fit);
      
    }
    
    // or discrete features 
    
    for (BP_Range::iterator dit = dr->begin(); dit <= dr->max(); dit++) {
      
      Int32 *i = o->getDiscFea(*dit);
      
      if (fwrite((Int32 *)i,sizeof(Int32),1,fp) != 1)
	error("DataOutStream::writeBinFeatures: couldn't write %i'th feature in frame \n",*dit,*fit);
      
    }
  }
}

/* write features from observation buffer to pfile
 * str: pfile output stream
 * o: observation matrix
 * cr: range of continuous features to be written
 * dr: range of discrete features to be written
 * fr: range of frames to process
 * sent_no: segment number of this utterance 
 */


void
DataOutStream::writePFileFeatures(OutFtrLabStream_PFile *str,
				  ObservationMatrix *o,
				  BP_Range *cr,
				  BP_Range *dr,
				  BP_Range *fr,
				  size_t sent_no) {

  size_t n_frames = fr->length();
  unsigned n_floats = cr->length();
  unsigned n_ints = dr->length();

  if (n_frames > _bufSize) 
    resizeBuffers(n_frames*2);

  float *fp = ftr_buf;      
  UInt32 *lp = lab_buf;

  for (BP_Range::iterator fit = fr->begin(); fit <= fr->max(); fit++) {
    
    o->gotoFrame(*fit);
    
    if (n_floats > 0) {
      
      // fill temporary feature buffer
      
      for (BP_Range::iterator cit = cr->begin(); cit <= cr->max(); cit++,fp++) {
	float *tmp_fp = o->getContFea(*cit);
	if (tmp_fp == NULL)
	  error("DataOutStream::writePFileFeatures: couldn't read %i'th feature for frame %i from buffer",*cit,*fit);
	*fp = *tmp_fp;  
      }
    }
    else if (n_ints > 0) {
      for (BP_Range::iterator cit = cr->begin(); cit <= cr->max(); cit++,lp++){
	UInt32 *tmp_lp = (UInt32 *)o->getDiscFea(*cit);
	if (tmp_lp == NULL)
	  error("DataOutStream::writePFileFeatures: couldn't read %i'th feature for frame %i from buffer",*cit,*fit);
	*lp = (UInt32)*tmp_lp; 
      } 
    }
  }

  // write out feature buffers


  if (n_floats > 0)
    str->write_ftrs(n_frames,ftr_buf);
  else if (n_ints > 0)
    str->write_labs(n_frames,lab_buf);

  str->doneseg(sent_no);

}

/* write ascii output
 * fp: file to write to
 * o: observation matrix 
 * cr: range of continuous features
 * dr: range of discrete features
 * fr: frame range
 */  

void
DataOutStream::writeAscFeatures(FILE *fp,
				  ObservationMatrix *o,
				  BP_Range *cr,
				  BP_Range *dr,
				  BP_Range *fr) {


  if (fp == NULL)
    error("DataOutStream::writeAscFeatures: output file is NULL\n");

  // for all frames in range


  for (BP_Range::iterator fit  = fr->begin(); fit <= fr->max(); fit++) {
    
    o->gotoFrame(*fit);
    
    // print continuous features    

    for (BP_Range::iterator cit = cr->begin(); cit <= cr->max(); cit++) {
      
      float *f = o->getContFea(*cit);
      
      if (f == NULL)
	error("DataOutStream::writeAscFeatures: couldn't read %i'th feature in frame %i from buffer\n",*cit,*fit);
      
      if (!(fprintf(fp,"%e ",*f)))
	error("DataOutStream::writeAscFeatures: couldn't write %i'th feature in frame %i\n",*cit,*fit);
      
    }
    
    // or print discrete features
    
    for (BP_Range::iterator dit = dr->begin(); dit <= dr->max(); dit++) {
      
      Int32 *i = o->getDiscFea(*dit);

      printf("%i\n",*dit);

      if (i == NULL)
	error("DataOutStream::writeAscFeatures: couldn't read %i'th feature in frame %i from buffer\n",*dit,*fit);
      
      if (!(fprintf(fp,"%i ",*i)))
	error("DataOutStream::writeAscFeatures: couldn't write %i'th feature in frame \n",*dit,*fit);
    }
    fprintf(fp,"\n");
  }
}


DataOutStream::~DataOutStream() {
  

  for (int i = 0; i < _numFiles; i++) {

    // this closes and writes out pfile stream

    if (pfile_ostr[i] != NULL)
      delete pfile_ostr[i];
  }

  if (ftr_buf != NULL)
    delete [] ftr_buf;
  if (lab_buf != NULL)
    delete [] lab_buf;
  delete [] outNames;
  delete [] pfile_ostr;
  delete [] outFiles;
  delete [] nInts;
  delete [] nFloats;
  delete [] sampRate;
  delete [] dataFormats;
  delete [] parKind;
  delete [] bswap;
}

/* resize temporary buffers 
 * n: number of frames (new size)
 */

  void
    DataOutStream::resizeBuffers(size_t n) {

  delete ftr_buf;
  delete lab_buf;

  _bufSize = n;

  ftr_buf = new float[_bufSize*_maxFloats];
  lab_buf = new UInt32[_bufSize*_maxInts];

}


