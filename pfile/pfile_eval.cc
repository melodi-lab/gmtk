#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <iostream.h>
#include <string.h>
#include <math.h>

#include "pfile.h"
#include "Range.H"
#include "error.h"

#define MINLOGARG 2.45E-308

#define PROB_LOW 0.00001
#define PROB_MIN 0.0001

#define USAGESTRING "USAGE: pfile_eval [options]\n
Option                                     default\n
-i   input pfile                           none\n
-e   compute entropy                       off\n
-bswap do byte-swapping                    off\n
"

FILE* 
iopen(const char *fname)
{
  FILE *file;
  char ch;

  if ((file=fopen(fname,"r")) != NULL) 
    {
      fprintf(stderr,"File %s exists, overwrite? (y/n) ",fname);
      if ((ch=getchar()) == 'y') 
	{
	  remove(fname);
	  file=fopen(fname,"wb");
	}
      else
	exit(03);
      rewind(stdin);
    }
  else 
    file=fopen(fname,"wb");
  return file;
}


int
GetMax(float *vec, int size) {

  float *fp;
  float m;
  int mInd = -1;
  int i;

  for (i=0,fp = vec;i < size;i++,fp++) {
    if (i == 0) {
      m = *fp;
      mInd = 0;
    }
    else if (*fp > m) {
      m = *fp;
      mInd = i;
    }
  }
  return mInd;
}

float
entrp(float *vec, int size) {

  float out = 0.0;
  float *fp;
  int i;
  

  for (i=0,fp=vec;i<size;i++,fp++) {
    if (*fp < MINLOGARG)
      continue;
    out += (*fp * log(*fp));
  }
  return -out;
}



void 
pfile_eval(InFtrLabStream_PFile *data,
		int entropy) {

  const size_t n_sents = data->num_segs();
  const size_t n_labs = data->num_labs();
  const size_t n_ftrs =  data->num_ftrs();
  const size_t n_frames = data->num_frames();
  float en_corr=0.0, en_false=0.0;
  float *ecorr, *efalse;
  float max_corr = 0.0, min_corr = 0.0;
  float max_false = 0.0, min_false = 0.0;
  float var = 0.0, corrmean = 0.0,falsemean = 0.0, diff;
  float corrvar = 0.0, falsevar = 0.0;
  size_t corr = 0, fals=0;
  int a;
  size_t new_id = 0;
  float *fp;
  UInt32 *lp;


  cout << "Input file has " << n_ftrs << " features, " << n_sents << \
    " sentences, " << n_frames << " frames\n";

  ecorr = new float[n_frames];
  efalse = new float[n_frames];

  long *priors = new long[n_ftrs];

  for (int p=0;p<n_ftrs;p++)
    priors[p] = 0;

  size_t buf_size = 300;
  float *ftr_buf = new float[n_ftrs * buf_size];
  UInt32* lab_buf = new UInt32[n_labs * buf_size];

  size_t frames_processed = 0;

  for (int s = 0; s < n_sents;s++) {

    int seg_frames = data->num_frames(s);
    if (seg_frames == SIZET_BAD) 
      {
	fprintf(stderr,
		"Couldn't find number of frames "
		"at sentence %lu in input pfile.\n",
		(unsigned long) s);
	error("Aborting.");
      }

    const SegID seg_id = data->set_pos(s,0);
    
    if (seg_id == SEGID_BAD) {
      fprintf(stderr,
	      "Couldn't seek to start of sentence %lu "
	      "in input pfile.",
	      (unsigned long) s);
      error("Aborting.");
    }
    if (seg_frames > buf_size) {
      delete ftr_buf;
      delete lab_buf;
      buf_size = seg_frames * 2;
      ftr_buf = new float[buf_size * n_ftrs];
      lab_buf = new UInt32[buf_size * n_labs];
    }
    const size_t n_read = data->read_ftrslabs(seg_frames,ftr_buf,lab_buf);
    
    if ( n_read != seg_frames) {
      fprintf(stderr, " At sentence %lu in input pfile, "
	      "only read %lu frames when should have read %lu.\n",
	      (unsigned long) s,
	      (unsigned long) n_read, (unsigned long) n_frames);
      error("Aborting.");
    }
    fp = ftr_buf;
    lp = lab_buf;
    new_id = 0;
    int cframe = 0;
    int fframe = 0;

    for (int f=0; f < n_read;f++) {
      int m = GetMax(fp,n_ftrs);
      int n = lab_buf[f];
      if (n > n_ftrs) {
	fprintf(stderr,"Warning: Label exceeds range: %i\n",n);
	fp += n_ftrs;
	lp++;
	continue;
      }
      if (m == -1)
	error("Can't find max output\n");

      if ( n == m) {
	if (entropy) {
	  ecorr[corr] = entrp(fp,n_ftrs);
	  if (ecorr[corr] > max_corr)
	    max_corr = ecorr[corr];
	  if (min_corr == 0.0) 
	    min_corr = ecorr[corr];
	  else if (ecorr[corr] < min_corr)
	    min_corr = ecorr[corr];
	  en_corr += ecorr[corr];
	}
	corr++;
      }
      else if (n != m && m != -1){
	if (entropy) {
	  efalse[fals] = entrp(fp,n_ftrs);
	  en_false += efalse[fals];
	  if (efalse[fals] > max_false)
	    max_false = efalse[fals];
	  if (min_false == 0.0)
	    min_false = efalse[fals];
	  else if (efalse[fals] < min_false)
	    min_false = efalse[fals];
	}
	fals++;
      }
      else {
	fprintf(stderr,
		"Failure finding max output unit at sent %li, frame %li\n",
		s,f);
	exit(01);
      }
      frames_processed++;
      fp += n_ftrs;
      lp++;
    }
    
  }
  // entropy means
  if (entropy) {
    corrmean = en_corr/(float)corr;
    falsemean = en_false/(float)fals;

  // variances
    for (a = 0; a < corr; a++) {
      diff = ecorr[a]-corrmean;
      corrvar += diff * diff;
    }
    corrvar /= corr;
    for (a = 0; a < fals; a++) {
      diff = efalse[a]-falsemean;
      falsevar += diff * diff;
    }
    falsevar /= fals;
  }
  

  // Output
  printf("correct: %li out of %li, %.2f%% percent\n",corr,frames_processed,
	 corr/(float)frames_processed*100.0);

  if (entropy) {
    printf("entropy corr: mean %f, var %f, min: %f, max: %f\n", 
	   corrmean,corrvar,min_corr, max_corr);
    printf("entropy false: %f, var %f, min %f, max %f\n",
	   falsemean,falsevar,min_false, max_false);
    printf("average entropy: %f\n",(corrmean + falsemean)/2);
    printf("entropy ratio: %f\n",corrmean/falsemean);
  }
  delete ftr_buf;
  delete lab_buf;
  return;
}


void 
main(int argc, const char **argv) {

  const char *inname=NULL,*outname=NULL,*outpname=NULL;      
  const char *program_name = *argv++;
  const char *pname=NULL;
  OutFtrLabStream_PFile *output = NULL;
  int entropy = 0;
  bool bswap = false;
  argc--;
  
  if (argc == 0)
    error(USAGESTRING);
  
  while (argc--)
    {
      const char *argp = *argv++;
      
      if (strcmp(argp, "-help")==0)
	{
	  printf("%s\n",USAGESTRING);
	}
      else if (strcmp(argp, "-i")==0)
	{
	  if (argc>0)
	    {
	      inname = *argv++;
	      argc--;
	    }
	  else
	    error("No input pfile name specified\n");
	}
      else if (strcmp(argp,"-e") == 0)
	{
	  entropy = 1;
	}
       else if (strcmp(argp,"-bswap") == 0)
       {
          bswap = true;
       }
      else {
	error("Unrecognized argument (%s).",argp);
      }
    }
  if (inname==0)
    fprintf(stderr,"No input pfile name supplied.");

  FILE *infile = fopen(inname, "r");
  if (infile == NULL)
    error("Couldn't open input pfile for input\n");

  FILE *priorsfile=stdout;
  if (pname != NULL)
    if ((priorsfile = fopen(pname,"w")) == NULL)
      error("Couldn't open priors file for output\n");

  int debug_level = 0;
  InFtrLabStream_PFile* input
    = new InFtrLabStream_PFile(debug_level, "intput stream", infile, 1,bswap);

  
  FILE *outfile = stdout;
  if (outname != NULL)
    outfile = iopen(outname);
  if (outfile == NULL)
    error("Couldn't open output pfile\n");

  
  pfile_eval(input,entropy);
  
  delete input;
  delete output;
  fclose(infile);
  fclose(outfile);
  return 0;
}
    





