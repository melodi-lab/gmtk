// Very Simple PFILE interface. Does error checking for you.
// 
// Written by: Jeff Bilmes
//             bilmes@icsi.berkeley.edu
//
// $Header$


#ifndef spi_h
#define spi_h

#include <time.h>
#include <stdio.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <math.h>
#include <sys/types.h>
#include <stdlib.h>
#include <math.h>
#include "sArray.h"
#include "QN_PFile.h"


class SPI_base {

public:

  virtual size_t n_ftrs() = 0;
  virtual size_t n_labs() = 0;
  virtual size_t n_segs() = 0;

  // returns total number of frames
  virtual size_t n_frms() = 0;
  // returns number of frames for segment 'pos'
  virtual size_t n_frms(const int pos) = 0;

  // buffer reads.
  virtual size_t read_ftrs(const size_t pos,
		   float*& fb) = 0;
  virtual size_t read_labs(const size_t pos,
		   size_t*& lb) = 0;  
  virtual size_t read_ftrslabs(const size_t pos,
		       float*& fb,
		       QNUInt32*& lb) = 0;
  virtual int read_ftrs(const int pos,
		float*& fb) = 0;
  virtual int read_labs(const int pos,
		int*& lb)  = 0;
  virtual int read_ftrslabs(const int pos,
		    float*& fb,
		    int*& lb) = 0;
};

class SPI : public SPI_base {

  size_t buf_size;

  sArray<float> ftr_buf;
  sArray<QNUInt32> lab_buf;

  char *local_fname;

  FILE *inf;
  QN_InFtrLabStream* in_streamp;

public:

  SPI(const char *const fname);
  virtual ~SPI();

  virtual size_t n_ftrs() { return in_streamp->num_ftrs(); }
  virtual size_t n_labs() { return in_streamp->num_labs(); }
  virtual size_t n_segs() { return in_streamp->num_segs(); }
  virtual size_t n_frms() { return in_streamp->num_frames(); }
  virtual size_t n_frms(const int pos) { return in_streamp->num_frames((size_t)pos); }

  // buffer reads.
  virtual size_t read_ftrs(const size_t pos,
		   float*& fb);
  virtual size_t read_labs(const size_t pos,
		   size_t*& lb);  
  virtual size_t read_ftrslabs(const size_t pos,
		       float*& fb,
		       QNUInt32*& lb);
  virtual int read_ftrs(const int pos,
		float*& fb) 
    { return (int) read_ftrs((size_t)pos,fb); }
  virtual int read_labs(const int pos,
		int*& lb)  
    { return (int) read_labs((size_t)pos,(size_t*&)lb); }
  virtual int read_ftrslabs(const int pos,
		    float*& fb,
		    int*& lb)
    { return (int) read_ftrslabs((size_t)pos,fb,(size_t*&)lb); }
};


// 
// A simple pfile inteface that reads two pfiles
// at the same time, merges them together so that
// the other programs need only know at the beginning
// that there are two pfiles. Merging occurs in the
// order that the file names are presented to the constructor 
// argument.
class SPI2 : public SPI_base {

  sArray<float> ftr_buf1;
  sArray<float> ftr_buf2;

  size_t buf_size;

  sArray<float> ftr_buf;


  sArray<QNUInt32> lab_buf;

  char *local_fname;
  char *local_fname2;

  FILE *inf;
  FILE *inf2;
  QN_InFtrLabStream* in_streamp;
  QN_InFtrLabStream* in_streamp2;

public:

  SPI2(const char *const fname,
      const char *const fname2);
  virtual ~SPI2();

  virtual size_t n_ftrs() { return in_streamp->num_ftrs() + 
			     in_streamp2->num_ftrs(); }
  virtual size_t n_labs() { return in_streamp->num_labs(); }
  virtual size_t n_segs() { return in_streamp->num_segs(); }
  virtual size_t n_frms() { return in_streamp->num_frames(); }
  virtual size_t n_frms(const int pos) { return in_streamp->num_frames((size_t)pos); }

  // buffer reads.
  virtual size_t read_ftrs(const size_t pos,
		   float*& fb);
  // For an SPI2, read_labs returns only the labels for the
  // first pfile.
  virtual size_t read_labs(const size_t pos,
		   size_t*& lb);  
  virtual size_t read_ftrslabs(const size_t pos,
		       float*& fb,
		       QNUInt32*& lb);



  virtual int read_ftrs(const int pos,
		float*& fb) 
    { return (int) read_ftrs((size_t)pos,fb); }

  virtual int read_labs(const int pos,
		int*& lb)  
    { return (int) read_labs((size_t)pos,(size_t*&)lb); }

  virtual int read_ftrslabs(const int pos,
		    float*& fb,
		    int*& lb)
    { return (int) read_ftrslabs((size_t)pos,fb,(size_t*&)lb); }
};



class SPO {

  char *local_fname;

  FILE *ouf;
  QN_OutFtrLabStream* out_streamp;
  
  size_t pos;

public:

  SPO(const char *const fname,
      const size_t n_ftrs,
      const size_t n_labs);
  ~SPO();

  // buffer write
  void write_ftrslabs(const size_t len,
		      float* ftr_buf,
		      QNUInt32* lab_buf);
  void write_ftrslabs(const int len,
		      float* ftr_buf,
		      int* lab_buf)
    {  write_ftrslabs((size_t)len,ftr_buf,(QNUInt32*)lab_buf); }

};


#endif
