static char rcsid = "$Id$";

#include <iostream.h>
#include <stdlib.h>
#include "GM_type.h"
#include "sArray.h"
#include "range.h"


//#include "pfile.h"
//#include "htk.h"
//#include ....

#define DATA_CHUNK_SIZE 300000 // number of bytes to read as one block

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
};
  
class GM_input {

  int num_files;  // number of input files
  int *formats;   // each can have different format

  char **fname; // data filename(s)
  FILE **fp;    // data file pointer(s)

  Range **ranges;  // range specifications for different files
  
  int n_floats;  // number of floats
  int n_ints;  // number of ints

  size_t record_size;  // size of a record in bytes

public:

  GM_file(int n_files,const char **names,FILE **files, int *nfloats, 
	  int *nints, int *formats);
  ~GM_file();

  int read_data(GM_rec *buf,int num_recs);
}

// output interface

class GM_output {

  int num_files;
  char **fnames;
  FILE **fp;

  int *n_floats;
  int *n_ints;

}
