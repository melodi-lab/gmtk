//static char *rcsid = "$Id$";

#include <stdio.h>
#include "GM_io.h"


GM_buf::GM_buf(size_t n_records, size_t rec_size) {

  
  
}


GM_rec::GM_rec(int n_floats, int n_ints) {

  
  fval.resize(n_floats);
  ival.resize(n_ints);

}

GM_rec::~GM_rec() {

  delete fval;
  delete ival;

}

GM_rec::print_rec(FILE *f) {

  for (int i=0;i<fval.len();i++) 
    fprintf(f,"%e ",fval[i].val());
  for (int j = 0; j < ival.len(); j++) 
    fprintf(f,"%i ",ival[j].val());
}


GM_input::GM_input(int n_files, const char **filename, 
		   int *nfloats, int *nints, int *format) {

  int i;

  num_files = n_files;
  record_size = 0; 

  fnames = new const char*[num_files];
  fp = new FILE*[num_files];
  n_floats = new int[num_files];
  n_ints = new int[num_files];
  formats = new int[num_files];
  n_frames = new size_t[num_files];

  for ( i = 0; i < num_files; i++) {
    fnames[i] = filename[i];
    n_floats[i] = nfloats[i];
    record_size += n_floats[i] * sizeof(float); 
    n_ints[i] = nints[i];
    record_size += n_ints[i] * sizeof(int);
    formats[i] = format[i];

    // open feature files and check format specifications

    printf("file: %s\n",fnames[i]);
    printf("format %i %i\n",i,format[i]);

    switch(formats[i]) {
    case RAWBIN:
      n_frames[i] = open_binary_file(fnames[i],fp[i],n_floats[i],n_ints[i]);
      break;
    case RAWASC:
      n_frames[i] = open_ascii_file(fnames[i],fp[i],n_floats[i],n_ints[i]);
      break;
    case HTK:
      n_frames[i] = open_htk_file(fnames[i],fp[i],n_floats[i],n_ints[i]);
      break;
    case GM:
      n_frames[i] = open_gm_file(fnames[i],fp[i],n_floats[i],n_ints[i]);
      break;
    default:
      GM_error("ERROR: GM_io::Invalid file format specified for file %i\n",i);
      exit(-1);
    }
    if (n_frames[i] == 0) 
      GM_error("ERROR: GM_io::Failure opening file '%s'\n",fnames[i]);
  }
  for (i = 1; i < num_files; i++) {
    if (n_frames[i] != n_frames[i-1]) {
      GM_error("ERROR: GM_io::Number of samples don't match \n file %s %i\n file %s %i\n",
	       fnames[i-1],n_frames[i-1],fnames[i],n_frames[i]);
    }
  }

  // allocate data buffers

  inp_buf = GM_rec(bufsize,record_size);
  
  //  range_size = parse_ranges();

  //  data = GM_rec(bufsize,range_size);
  
}

// open binary file and check number of bytes

int
GM_input::open_binary_file(const char *fname, FILE *fp, int nfloats, int nints) {

  if ((fp = fopen(fname,"rb")) == NULL) {
    GM_warning("GM_io::open_binary_file: Can't open '%s' for input\n",fname);
    return 0;
  }

  if (fseek(fp,0L,SEEK_END) == -1) {
    GM_warning("GM_io::open_binary_file: Can't skip to end of file",fname);
    return 0;
  }

  size_t file_size = ftell(fp);

  rewind(fp);

  int rec_size = nfloats * sizeof(float) + nints * sizeof(int);
  
  if ((file_size % rec_size) > 0) {
    GM_warning("GM_io::open_binary_file: odd number of bytes in file %s\n",fname);
    return 0;
  }

  int n_samples = file_size / rec_size;
  
  return n_samples;
}


int
GM_input::open_ascii_file(const char *fname,FILE *fp, int nfloats, int nint) {

  size_t n_samples = 0;
  char ch;

  if ((fp = fopen(fname,"r")) == NULL) {
    GM_warning("GM_io::open_ascii_file: Can't open '%s' for input\n",fname);
    return 0;
  }

  // for ascii, newline is record delimiter

  while ((ch = fgetc(fp)) != EOF) {
    if (ch == '\n')
      n_samples++;
  }

  rewind(fp);

  return n_samples;
}


int 
GM_input::open_gm_file(const char *fname,FILE *fp, int nfloats, int nint) {

  // TODO 

  ;

}


// open htk file

int
GM_input::open_htk_file(const char *fname,FILE *fp, int nfloats, int nints) {

  int n_samples;
  int samp_period;
  short samp_size;
  short parm_kind;
  
  if ((fp = fopen(fname,"rb")) == NULL) {
    GM_warning("GM_io::open_htk_file: Can't open '%s' for input\n",fname);
    return 0;
  }
  if (fseek(fp,0L,SEEK_END) == -1) {
    GM_warning("GM_io::open_binary_file: Can't skip to end of file",fname);
    return 0;
  }
  size_t filesize = ftell(fp);

  rewind(fp);

  if (fread(&n_samples,sizeof(int),1,fp) != 1) {
    GM_warning("GM_io::open_htk_file: Can't read number of samples\n");
    return 0;
  }

  if (fread(&samp_period,sizeof(int),1,fp) != 1) {
    GM_warning("GM_io::open_htk_file: Can't read sample period\n");
    return 0;
  }

  if (fread(&samp_size,sizeof(short),1,fp) != 1) {
    GM_warning("GM_io::open_htk_file: Can't read sample size\n");
    return 0;
  }

  if (fread(&parm_kind,sizeof(short),1,fp) != 1) {
    GM_warning("GM_io::open_htk_file: Can't read parm kind\n");
    return 0;
  }

  int n_fea;

  if (parm_kind == DISCRETE) {
    n_fea = samp_size / sizeof(int) ;
    if (n_fea != nints) {
      GM_warning("GM_io::open_htk_file:  Number of features in file (%i) does not match number of ints specified (%i)\n", n_fea,n_ints);
      return 0;
    }
  }
  else {
    n_fea = samp_size / sizeof(float);
    if (n_fea != nfloats) {
      GM_warning("GM_io::open_htk_file:  Number of features in file (%i) does not match number of floats specified (%i)\n", n_fea,n_floats);
      return 0;
    }
  }

  printf("nfea: %i\n",n_fea);
  return n_samples;
}

 
GM_input::~GM_input() {

  delete [] n_ints;
  delete [] n_floats;
  delete [] formats;
}



GM_output::GM_output() {
;
}

GM_output::~GM_output() {
;
}

// read a single frame

int
GM_input::read_frame() {

  int n_read = 0;
  float f;

  for (int i = 0; i < num_files; i++) {
    if (n_floats[i] > 0) {
      if (formats[i] == ASCII) {
	while (n_read < n_floats[i]) {
	  if (fscanf(fp[i],"%e",&f) == 1)
	    inp_buf.read_float(f);
	  n_read++;
	}
      }
      else 
	fread(,sizeof(float),n_floats[i],fp[i]);
    }
    if (n_ints[i] > 0) {
      
  }
}

// read a single file

int
GM_input::read_file() {

  // first read all data in buffer

  // resize buffer if necessary

  if (n_frames[0] > buf_size) {
    if (!(resize_buffer(n_frames[0])))
      GM_error("ERROR: GM_io::Memory failure, file too large\n");
    
  for (size_t i = 0;i < n_frames[0]; i++)
    if (!(read_record()))
      return 0;

  return n_frames[0];
  
}





