#include <stdio.h>
#include <stddef.h>
#include <iostream.h>
#include <math.h>
#include <stdlib.h>
#include "pfile.h"
#include "error.h"


#define USAGESTRING "USAGE: pfile2rapbin [options]\n
Option                                 default\n
-i1   input pfile                       none\n
-o    output file                       none\n
-bswap do byteswapping                  off\n
-debug  debug level                     0\n
"

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


void pfile2rapbin(InFtrLabStream_PFile *data,
	          FILE *outfile)

{

  const int sep = -1; // rapbin sentence separator
  int num_phones = 0;
  const size_t n_sents = data->num_segs();
  const size_t n_labs = data->num_labs();
  const size_t n_ftrs =  data->num_ftrs();
  const size_t n_frames = data->num_frames();
  const size_t head = n_ftrs;

  cout << "Input file has " << n_ftrs << " features, " << n_sents << \
      " sentences, " << n_frames << " frames\n";

  printf("Writing rapbin file\n");

  size_t buf_size = 300;
  float *ftr_buf = new float[n_ftrs * buf_size];
  UInt32* lab_buf = new UInt32[n_labs * buf_size];

  fwrite(&head,sizeof(head),1,outfile); // write rapbin file header

  for (int r=0;r < data->num_segs(); r++) {

    const size_t seg_frames = data->num_frames(r);
	

    if (seg_frames == SIZET_BAD) 
      {
	error("Couldn't find number of frames "
		"at sentence %lu in input pfile.\n",
		(unsigned long) r);
      }
    const SegID seg_id = data->set_pos(r, 0);
    if (seg_id == SEGID_BAD) 
      error("Couldn't seek to start of sentence %lu in input pfile.\n",r);

    if (seg_frames > buf_size) {
      delete ftr_buf;
      delete lab_buf;
      buf_size = seg_frames * 2;
      ftr_buf = new float[buf_size * n_ftrs];
      lab_buf = new UInt32[buf_size * n_labs];
    }

    size_t n_read = data->read_ftrslabs(seg_frames,ftr_buf,lab_buf);
    
    if ( n_read != seg_frames) {
      error("At sentence %lu in input pfile, "
	      "only read %lu frames when should have read %lu.\n",
	      (unsigned long) r,
	    (unsigned long) n_read, (unsigned long) n_frames);
    }
    float *fp = ftr_buf;
    for (int s = 0; s < n_read; s++) {
      fwrite(&s,sizeof(s),1,outfile);   // output framenumber
      if ((fwrite((float *)fp,sizeof(float),n_ftrs,outfile)) != n_ftrs)
	error("Fewer items written than requested at sent %li frame %li\n",r,s); 
      fp += n_ftrs;
    }
    fwrite((int *)&sep,sizeof(sep),1,outfile); 
  }

  delete ftr_buf;
  delete lab_buf;
  return;
}


main(int argc, const char **argv) {

  const char *inname = NULL,*outname = NULL;
  int debug_level = 0;
  FILE *infile=NULL, *outfile = NULL;
  const char *program_name = *argv++;
  bool swap = false;
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
      else if (strcmp(argp,"-o") == 0)
	{
	  if (argc>0)
	    {
	      outname = *argv++;
	      argc--;
	    }
	  else
	    error("Second input pfile name unspecified\n");
	}
       else if (strcmp(argp,"-bswap") == 0) {
          swap = true;
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
      else {
	error("Unrecognized argument (%s).",argp);
      }
    }

  if (inname == NULL) 
    error("No input pfile name supplied\n");
  if (outname == NULL)
    error("No output pfile name supplied - can't write binary to stdout\n");

  if ((infile = fopen(inname,"r")) == NULL)
      error("Couldn't open input pfile for reading: %s",inname);
  if ((outfile = fopen(outname, "w")) == NULL)
    error("Couldn't open output pfile %s for writing.",outname);

  InFtrLabStream_PFile* input
    = new InFtrLabStream_PFile(debug_level,"",infile,1,swap);
  pfile2rapbin(input,outfile);
  
  delete input;
  fclose(infile);
  fclose(outfile);
  return 0;
}
    





