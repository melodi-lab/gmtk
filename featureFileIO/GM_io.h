//static char *rcsid = "$Id$";

#include <stdio.h>
#include <stdlib.h>
#include "GM_type.h"
#include "sArray.h"
#include "Range.H"

#define BUFSIZE 300 // initial buffer size (in frames)

//#include "pfile.h"
//#include "htk.h"
//#include ....

// what we need:
// pfile support: read featues and labels from 1 file
// all others: read features and labels from different files
// memory allocation
// error checking (truncated files etc., mismatch feature-label files)
// support for mixes (int and float) files
// support feature ranges for all inputs

// floats and ints should be contiguous
// format for labels?

// supported file formats: raw binary, raw ascii, pfile, htk



enum {
  RAWBIN, 
  RAWASC, 
  PFILE,
  HTK,
  GM,
};

//HTK parameter kinds

enum {
  WAVEFORM,            
  LPC,LPREFC,LPCEPSTRA,LPDELCEP,   
  IREFC,                           
  MFCC,                            
  FBANK,                           
  MELSPEC,                         
  USER,                           
  DISCRETE,                        
  ANON
};

class GM_rec {

  sArray<float> fval; // float part of input record
  sArray<int> ival;   // int part of input record

 public:

  GM_rec(size_t n_frames, size_t rec_size,int n_floats, n_ints);
  ~GM_rec();
  
  print_record(FILE *);
//  copy_record();
//  resize_rec();
  

};
  
class GM_input {

private:

  int num_files;  // number of input files
  int *formats;   // each can have different format

  const char **fnames; // data filename(s)
  FILE **fp;    // data file pointer(s)

  int *n_floats;  // number of floats
  int *n_ints;  // number of ints
  size_t *n_frames; // # frames per file

  float *buf_p;			   

  size_t record_size;  // size of a record in bytes

  int *swap_flags;
  
  Range **ftr_rng; // ftr ranges (can be different for each file)
  Range *sr_rng;   // sentence range (same for all files)
  Range *pr_rng;   // per_sentence range (same for all files)

  GM_rec *inp_buf; // input data buffer
  GM_rec *data;    // data presented to user (range-selected)		       

  int open_binary_file(const char *fname, FILE *fp, int nfloats, int nint);
  int open_ascii_file(const char *fname, FILE *fp, int nfloats, int nint);
  int open_htk_file(const char *fname, FILE *fp, int nfloats, int nint);
  int open_gm_file(const char *fname, FILE *fp, int nfloats, int nint);
  

//int open_pfile(char *fname, FILE *fp, int nfloats, int nint);

public:

  GM_input(int n_files,const char **names, int *nfloats, 
	  int *nints, int *formats);
  ~GM_input();

  int read_data(int num_recs);
};


// output interface

class GM_output {


  int num_files;
  char **fnames;
  FILE **fp;
  int *formats;

  int *n_floats;
  int *n_ints;

public:
  GM_output();
  ~GM_output();

};






