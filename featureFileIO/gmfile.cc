static char *rcsid = "$Id$";

#include <iostream.h>
#include "GM_io.h"


GM_input(int n_files, const char **filename, FILE **f, int *nfloats, 
	int *nints, int *format, Range **range_specs) {

  num_files = n_files;

  fnames = new char*[num_files];
  fp = new FILE*[num_files];
  n_floats = new int[num_files];
  n_ints = new int[num_files];
  formats = new int[num_files];
  ranges = range_specs;
  
  record_size = 0;		

  for (int i = 0; i < num_files; i++) {
    fnames[i] = filename[i];
    fp[i] = f[i];
    n_floats[i] = nfloats[i];
    record_size += n_floats[i] * sizeof(float); 
    n_ints[i] = nints[i];
    record_size += n_ints[i] * sizeof(int);
    formats[i] = format[i];
  }
};

 
~GM_input() {

  delete [] n_ints;
  delete [] n_floats;
  delete [] formats;
};  



GM_output() {

};

~GM_output() {

};

// main function for reading in data; note that buffer
// needs to be allocated externally 

int
GM_input::read_data(GM_rec *buf, int n_records) {

  // first read floats

  for (int i = 0; i < num_files; i++) {
    if (n_floats[i] > 0) {
      for (Range::iterator f = ranges[i].begin(); !f.at_end(); ++f) {
	fread(bp,sizeof(float),n_elems,file_p[i]);
	bp += n_elems;
      }
      // then read ints
      
}
