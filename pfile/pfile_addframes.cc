#include <stdio.h>
#include <stddef.h>
#include <iostream.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "pfile.h"
#include "error.h"

#define USAGESTRING "USAGE: pfile_addframes [options]\n
Option                                     default\n
-i   input pfile                           none\n
-o   output pfile                          none\n
-c   context                               none\n
-l   list                                  none\n
-b   label                                 none\n
-iswap do byteswapping on input            off\n
-oswap do byteswapping on output           off\n 
-debug  debug level                        0\n
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



void pfile_addframes(InFtrLabStream_PFile *data,
		OutFtrLabStream_PFile *outdata,
		size_t sr_rng,
		size_t fr_rng,
		size_t context)
{
  const size_t n_sents = data->num_segs();
  const size_t n_labs = data->num_labs();
  const size_t n_ftrs =  data->num_ftrs();
  const size_t n_frames = data->num_frames();

  cout << "Input file has " << n_ftrs << " features, " << n_sents << \
    " sentences, " << n_frames << " frames\n";

  size_t buf_size = 300;
  float *ftr_buf = new float[n_ftrs * buf_size];
  UInt32* lab_buf = new UInt32[n_labs * buf_size];

  
    size_t outbuf_size = buf_size + (2 * context);
    float *out_ftrbuf = new float[n_ftrs * outbuf_size];
    UInt32* out_labbuf = new UInt32[n_labs * outbuf_size];

    for (long srit=0;srit<sr_rng;srit++) {
      const size_t seg_frames = data->num_frames(srit);
      if (srit % 100 == 0)
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
	delete out_labbuf;
	delete out_ftrbuf;
	buf_size = seg_frames * 2;
	outbuf_size = buf_size + (2 * context);
	ftr_buf = new float[buf_size * n_ftrs];
	lab_buf = new UInt32[buf_size * n_labs];
	out_labbuf = new UInt32[outbuf_size * n_labs];
	out_ftrbuf = new float[outbuf_size * n_ftrs];
      }
      
      const size_t n_read = data->read_ftrslabs(seg_frames,ftr_buf,lab_buf);
      
      if ( n_read != seg_frames) {
	fprintf(stderr, " At sentence %lu in input pfile, "
		"only read %lu frames when should have read %lu.\n",
		(unsigned long)srit,
		(unsigned long) n_read, (unsigned long) n_frames);
	error("Aborting.");
      }
      
      int pos = 0;  int opos=0; int nout=0;
      for (int c = 0; c < context; c++ )
	outdata->write_ftrslabs(1,ftr_buf,lab_buf);
      outdata->write_ftrslabs(n_read,ftr_buf,lab_buf);
      int nf = n_ftrs*n_read;
      int end = nf-n_ftrs;
      pos =  end;
      opos = n_read-1;
      for (int c = 0; c < context; c++ )
	outdata->write_ftrslabs(1,&ftr_buf[pos],&lab_buf[opos]);
      outdata->doneseg((SegID) seg_id);
    }
    delete ftr_buf;
    delete lab_buf;
    delete out_ftrbuf;
    delete out_labbuf;
    return;
}


int main(int argc, const char **argv) {

  const char *inname=NULL,*outname=NULL;      
  int debug_level = 0;
  size_t context = 0;   
  bool iswap = false, oswap = false;
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
      
      else if (strcmp(argp,"-c") == 0)
	{
	  if (argc>0)
	    {
	      context = atoi(*argv++);
	      argc--;
	    }
	  else
	    error("Context size unspecified.");
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
        else if (strcmp(argp,"-iswap") == 0)
         {
           iswap = true;
         }
        else if (strcmp(argp,"-oswap") == 0)
          {
           oswap = true;
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
    = new InFtrLabStream_PFile(debug_level, "", infile, 1,iswap);
  
  FILE *outfile = stdout;
  if (outname != NULL)
    outfile = iopen(outname);
  if (outfile == NULL)
    error("Couldn't open output pfile\n");

  
  if (context == 0)
    error("No context range specified\n");
  

  OutFtrLabStream_PFile *output = 
    new OutFtrLabStream_PFile(debug_level,"",outfile,
				 input->num_ftrs(),input->num_labs(),1,	
	                         oswap);
  
 
  
  pfile_addframes(input,output,input->num_segs(),input->num_ftrs(),context);
  
  delete input;
  delete output;
  fclose(infile);
  fclose(outfile);
  return 0;
}
    





