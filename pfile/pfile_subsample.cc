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

#define USAGESTRING "USAGE: pfile_subsample [options]\n
Option                                     default\n
-i   input pfile                           none\n
-o   output pfile                          none\n
-n   # labels                              0\n
-iswap  byteswap input pfile               off\n
-oswap  byteswap output pfile              off\n
-w   # window frames on each side          0\n
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

void
print_vec(FILE *f,float *vec, int size,int frameno, int sentno) {

  int i;
  float *fp;
  fprintf(f,"%i %i ",frameno,sentno);
  for (i = 0, fp = vec; i < size; i++, fp++)
    fprintf(f,"%e ",*fp);
  fprintf(f,"\n");
}


void 
pfile_subsample(InFtrLabStream_PFile *data, OutFtrLabStream_PFile *outdata,
		int total_labs,int window)
{
  
  const size_t n_sents = data->num_segs();
  const size_t num_labs = data->num_labs();
  const size_t n_ftrs =  data->num_ftrs();
  const size_t n_frames = data->num_frames();

  float *fp;
  UInt32 *lp, *right,*left;
  int offset;


  cout << "Input file has " << n_ftrs << " features, " << n_sents << \
    " sentences, " << n_frames << " frames\n";


  size_t buf_size = 300;
  float *ftr_buf = new float[n_ftrs * buf_size];
  float *tmp = new float[n_ftrs * buf_size];
  UInt32* lab_buf = new UInt32[num_labs * buf_size];

  UInt32 *lab_counts = new UInt32[total_labs];

  for (int s = 0; s < n_sents;s++) {
    if (s % 100 == 0)
      printf("Processing sentence %i\n",s);
    
    int seg_frames = data->num_frames(s);
    
    if (seg_frames == SIZET_BAD) 
      {
	fprintf(stderr,
		"Couldn't find number of frames "
		"at sentence %lu in input pfile.\n",
		(unsigned long) s);
	//	error("Aborting.");
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
      lab_buf = new UInt32[buf_size * num_labs];
      tmp = new float[buf_size * n_ftrs];
    }
    const size_t n_read = data->read_ftrslabs(seg_frames,ftr_buf,lab_buf);
    
    if (n_read != seg_frames) {
      fprintf(stderr, " At sentence %lu in input pfile, "
	      "only read %lu frames when should have read %lu.\n",
	      (unsigned long) s,
	      (unsigned long) n_read, (unsigned long) n_frames);
      error("Aborting.");
    }
    fp = ftr_buf;
    lp = lab_buf;

    for (int f=0; f < n_read;f++) {
      int m = GetMax(fp,n_ftrs);
      int n = *lp;
      lab_counts[n]++;
    }
  }
  size_t min_count = lab_counts[0];
  int min_lab = 0;
  for (int i = 1; i < total_labs; i++) {
    if (lab_counts[i]< min_count) {
      min_count = lab_counts[i];
      min_lab = i;
    }
  }
  printf("most frequent label: %i\n",max_lab);
  for (int s = 0; s < n_sents;s++) {
    if (s % 100 == 0)
      printf("Processing sentence %i\n",s);
    
   const size_t n_read = data->read_ftrslabs(seg_frames,ftr_buf,lab_buf);
   fp = ftr_buf;
   lp = lab_buf;
   
   for (int f = 0; f < n_read; f++) {
     if (lp == 
}

/*
       if (m != n) {
	 if (f - offset < 0)
	   left = &lab_buf[0];
	 else
	   left = &lab_buf[f-offset];
	 if (f + offset >= n_read)
	   right = &lab_buf[n_read-1];
	 else
	   right = &lab_buf[f+offset];
	m = check_labwindow(m,n,left,right,window);
       }
      if (m == n) {
	if (corr == 1 && nc_samp == 0) {
	  tmp = 1;
	  outdata->write_ftrslabs(1,fp,&tmp);
	}
	else if (corr == 1 && nc_samp > 0) {
	  if (f % c_samp_per_sent == 0) {
	    tmp = 1;
	    outdata->write_ftrslabs(1,fp,&tmp);
	  }
       }
      }
      if (m != n) {
	if (fls == 1 && nf_samp == 0) {
	  tmp = 0;
	  outdata->write_ftrslabs(1,fp,&tmp);
	}
	else if (fls == 1 && nf_samp > 0) {
	  if (f % f_samp_per_sent == 0) {
	    tmp = 0;
	    outdata->write_ftrslabs(1,fp,&tmp);
	  }
	}
      }
      fp += n_ftrs;
      lp++;
    }
       outdata->doneseg((SegID) s);
  }

  

  delete tmp;
  delete ftr_buf;
  delete lab_buf;
  return;
  } */


void 
main(int argc, const char **argv) {

  const char *inname=NULL,*outname=NULL,*outpname=NULL;      
  const char *program_name = *argv++;
  const char *pname=NULL;
  int n_labs = 0, window = 0;

  bool iswap = true, oswap = true;

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
      else if (strcmp(argp,"-iswap") == 0)
      {
         iswap = true;
      }
      else if (strcmp(argp,"-oswap") == 0)
      {
         oswap = true;
      }
      else if (strcmp(argp,"-o")==0)
	{
	  if (argc>0)
	    {
	      outname = *argv++;
	      argc--;
	    }
	  else
	    error("No output pfile name specified\n");
	}
      else if (strcmp(argp,"-n")==0)
	{
	  n_labs = atoi(*argv++);
	  argc--;
	}
      else if (strcmp(argp,"-w") == 0)
	{
	  window = atoi(*argv++);
	  argc--;
	}
      else
	error("Unrecognized argument\n");
    }
  if (inname==NULL)
    fprintf(stderr,"No input pfile name supplied.");
  if (outname == NULL)
    fprintf(stderr,"No output pfile name supplied.");
  
  FILE *infile = fopen(inname, "r");
  if (infile == NULL)
    error("Couldn't open input pfile for input\n");

  FILE *outfile = fopen(outname, "w");
  if (outfile == NULL)
    error("Couldn't open '%s' for output\n",outname);

  InFtrLabStream_PFile* input
    = new InFtrLabStream_PFile(0, "intput stream", infile, 1,iswap);

  OutFtrLabStream_PFile* output
    = new OutFtrLabStream_PFile(0,
				   "output stream",
				   outfile,
				   input->num_ftrs(),
				   input->num_labs(),
				   1,oswap);
  
  pfile_subsample(input,output,n_labs,window);
  
  delete input;
  delete output;
  fclose(infile);
  fclose(outfile);
  return 0;
}
    





