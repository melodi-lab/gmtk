#include <stdio.h>
#include <stddef.h>
#include <iostream.h>
#include <string.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "pfile.h"
#include "Range.H"
#include "error.h"

#define USAGESTRING "USAGE: pfile_entropy [options]\n
pOption                                     default\n
-i   input pfile                           none\n
-l   list of labels                        none\n
-bswap  byte-swap input                    off\n
-debug  debug level                        0\n
"

#define ZERO  0.0
#define PROB_FLOOR 0.000001


//int atoi(const char *nptr);


FILE * iopen(const char *fname)
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



void pfile_entropy(InFtrLabStream_PFile *data,
		Range &sr_rng,
		Range &fr_rng,
		UInt32 lab)
{

  float sum = 0.0;
  const size_t n_sents = data->num_segs();
  const size_t n_labs = data->num_labs();
  const size_t n_ftrs =  data->num_ftrs();
  const size_t n_frames = data->num_frames();

  size_t buf_size = 300;
 
  float *ftr_buf = new float[n_ftrs * buf_size];
 
  UInt32* lab_buf = new UInt32[n_labs * buf_size];
 
  cout << "Input file has " << n_ftrs << " features, " << n_sents << \
    " sentences, " << n_frames << " frames\n";
  
  for (int srit=0;srit<n_sents;srit++) {
	  const size_t seg_frames = data->num_frames(srit);

	  	  if ((srit % 100) == 0)
	  	    printf("Processing sentence %d\n",srit);

	  if (seg_frames == SIZET_BAD) 
	    {
	      fprintf(stderr,
		      "Couldn't find number of frames "
		      "at sentence %lu in input pfile.\n",
		       (unsigned long) srit);
	      error("Aborting.");
	    }
	  const SegID seg_id = data->set_pos(srit, 0);
	  
	  if (seg_id == SEGID_BAD) {
	    fprintf(stderr,
		    "Couldn't seek to start of sentence %lu "
		    "in input pfile.",
		    (unsigned long) srit);
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
		    (unsigned long) srit,
		    (unsigned long) n_read, (unsigned long) n_frames);
	    error("Aborting.");
	  }
	 
	  int pos = 0;  

	  for (int s = 0; s < n_read; s++) {
	    float  H = 0.0;
	    for (int f = pos; f < pos+n_ftrs; f++) {
	      if (ftr_buf[f] == ZERO)
		ftr_buf[f] = PROB_FLOOR;
	      H += (ftr_buf[f] * log(ftr_buf[f]));
	    }
	    sum += -(H);
	    pos += n_ftrs;
	  }
	  //	    printf("%li %f\n",srit,sum/(float)n_read);
  }
  printf("average entropy: %f\n",sum/n_frames);
  delete ftr_buf;
  delete lab_buf;
  return;
}


void main(int argc, const char **argv) {

  const char *inname;      
  const char *sr_str = NULL;   // sentence range string
  Range *sr_rng;
  const char *fr_str = NULL;   // feature range string
  Range *fr_rng;
  UInt32 lab;
  int debug_level = 0;
  bool bswap = false;
    
  const char *program_name = *argv++;
  argc--;
  
  if (argc == 0)
    error(USAGESTRING);
  
  while (argc--)
    {
      char buf[BUFSIZ];
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
     
     
      else if (strcmp(argp, "-debug")==0)
	{
	  // debug level
	  if (argc>0)
	    {
	      debug_level = atoi(*argv++);
	      argc--;
	    }
	  else
	    error("No debug level given\n");
	}
      else if (strcmp(argp, "-sr")==0)
	{
	  if (argc>0)
	    {
	      sr_str = *argv++;
	      argc--;
	    }
	  else
	    error("No range given\n");
	}
      else if (strcmp(argp, "-fr")==0)
	{
	  if (argc>0)
	    {
	      fr_str = *argv++;
	      argc--;
	    }
	  else
	    error("No range given\n");
	}
      else if (strcmp(argp,"-bswap") == 0) 
      {
         iswap = true;
      }
      else if (strcmp(argp, "-l")==0)
	{
	  if (argc>0)
	    {
	      lab = atoi(*argv++);
	      argc--;
	    }
	  else
	    error("No range given\n");
	}
      
      else {
	error("Unrecognized argument (%s).",argp);
      }
    }
  if (inname==0)
    fprintf(stderr,"No input pfile name supplied.");

  FILE *infile = fopen(inname, "r");
  if (infile == NULL)
    error("Couldn't open input pfile for reading.");
  
  InFtrLabStream_PFile* input
    = new InFtrLabStream_PFile(debug_level, "", infile, 1,bswap);
  
 

  sr_rng = new Range(sr_str,0,input->num_segs());
  fr_rng = new Range(fr_str,0,input->num_ftrs()); 

  
 
  
  pfile_entropy(input,*sr_rng,*fr_rng,lab);
  
  delete input;
  delete sr_rng;
  delete fr_rng;
  fclose(infile);
  return 0;
}
    





