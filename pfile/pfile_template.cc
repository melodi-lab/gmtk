/*
    $Header$
  
    This program selects an arbitrary subset set of features from each
    frame of a pfile and creates a new pfile with that
    subset in each frame.
    Jeff Bilmes <bilmes@cs.berkeley.edu>
*/

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <values.h>
#include <math.h>
#include <assert.h>

#include "pfile.h"
#include "parse_subset.h"
#include "error.h"
#include "Range.H"
#include "spi.h"

#define MAXHISTBINS 1000


#define MIN(a,b) ((a)<(b)?(a):(b))

static const char* program_name;


extern size_t bin_search(float *array,
			 size_t length, // the length of the array
			 float val);     // the value to search for.


static void
usage(const char* message = 0)
{
    if (message)
        fprintf(stderr, "%s: %s\n", program_name, message);
    fprintf(stderr, "Usage: %s <options>\n", program_name);
    fprintf(stderr,
	    "Where <options> include:\n"
	    "-help           print this message\n"
	    "-i <file-name>  input pfile\n"
	    "[-i2 <file-name>]  optional second input feature-only pfile to be merged online with first\n"
	    "-o <file-name>  output pfile (default stdout) \n"
	    "-q              quiet mode\n"
	    "-sr <range>     sentence range\n"
	    "-fr <range>     feature range\n"
	    "-pr <range>     per-sentence range\n"
	    "-lr <range>     label range\n"
	    "-b              print raw binary data (ints and floats)\n"
            "-iswap1         byte-swap first input pfile\n"
            "-iswap2         byte-swap second input pfile\n"
            "-oswap          byte-swap output pfile\n"
            "-ns             Don't print the frame IDs (i.e., sent and frame #)\n"
	    "-debug <level>  number giving level of debugging output to produce 0=none.\n"
    );
    exit(EXIT_FAILURE);
}


static void
pfile_process(SPI_base *in_stream,
	      SPO *out_stream,
	      Range& srrng,
	      Range& frrng,
	      Range& lrrng,
	      const char * pr_str,
	      const bool print_frameid,
	      const bool quiet,
	      const bool binary)

{
    // Feature and label buffers are dynamically grown as needed.
    size_t buf_size = 300;      // Start with storage for 300 frames.
    const size_t n_labs = in_stream.n_labs();
    const size_t n_ftrs = in_stream.n_ftrs();
    float *ftr_buf_p;

    UInt32* lab_buf;
    UInt32* lab_buf_p;

    //
    // Go through input pfile to get the initial statistics,
    // i.e., max, min, mean, std, etc.
    for (Range::iterator srit=srrng.begin();!srit.at_end();srit++) {

        const size_t n_frames = in_stream.n_frames(*srit);

	Range prrng(pr_str,0,n_frames);

	if (!quiet) {
	  if (*srit % 100 == 0)
	    printf("Processing sentence %d\n",*srit);
	}


	in_stream.read_ftrslabs((*srit), ftr_buf, lab_buf);

	for (Range::iterator prit=prrng.begin();
	     !prit.at_end() ; ++prit) {

	  bool ns = false;
	  ftr_buf_p = ftr_buf + (*prit)*n_ftrs;
	  lab_buf_p = lab_buf + (*prit)*n_labs;

	  if (print_frameid) {
	    if (binary) {
	      const size_t sent_no = *srit;
	      const size_t frame_no = *prit;
	      fwrite(&sent_no,sizeof(sent_no),1,out_fp);
	      fwrite(&frame_no,sizeof(frame_no),1,out_fp);
	    } else {
	      fprintf(out_fp,"%d %u",*srit,*prit);
	      ns = true;
	    }
	  }

	  for (Range::iterator frit=frrng.begin();
	       !frit.at_end(); ++frit) {
	    if (binary) {
	      fwrite(&ftr_buf_p[*frit],
		     sizeof(ftr_buf_p[*frit]),
		     1,out_fp);
	    } else {
	      if (ns) fprintf(out_fp," ");
	      fprintf(out_fp,"%f",ftr_buf_p[*frit]);
	      ns = true;
	    }
	  }
	  for (Range::iterator lrit=lrrng.begin();
	       !lrit.at_end(); ++lrit) {
	    if (binary) {
	      fwrite(&lab_buf_p[*lrit],
		     sizeof(lab_buf_p[*lrit]),
		     1,out_fp);
	    } else {
	      if (ns) fprintf(out_fp," ");	    
	      fprintf(out_fp,"%d",lab_buf_p[*lrit]);
	      ns = true;
	    }
	  }
	  if (!binary)
	    fprintf(out_fp,"\n");
	}

    }
    delete ftr_buf;
    delete lab_buf;
}


static long
parse_long(const char*const s)
{
    size_t len = strlen(s);
    char *ptr;
    long val;

    val = strtol(s, &ptr, 0);

    if (ptr != (s+len))
        error("Not an integer argument.");

    return val;
}

static float
parse_float(const char*const s)
{
    size_t len = strlen(s);
    char *ptr;
    double val;
    val = strtod(s, &ptr);
    if (ptr != (s+len))
        error("Not an floating point argument.");
    return val;
}


main(int argc, const char *argv[])
{
    //////////////////////////////////////////////////////////////////////
    // TODO: Argument parsing should be replaced by something like
    // ProgramInfo one day, with each extract box grabbing the pieces it
    // needs.
    //////////////////////////////////////////////////////////////////////

    const char *input_fname = 0;  // Input pfile name.
    const char *input_fname2 = 0;   // 2nd input pfile name, if provided
    const char *output_fname = 0; // Output pfile name.

    const char *sr_str = 0;   // sentence range string
    Range *sr_rng;
    const char *fr_str = 0;   // feature range string    
    Range *fr_rng;
    const char *lr_str = 0;   // label range string    
    Range *lr_rng;
    const char *pr_str = 0;   // per-sentence range string

    int debug_level = 0;
    bool print_frameid = true;
    bool quiet = false;
    bool binary=false;

    bool iswap1 = true, iswap2 = true, oswap = false;

    program_name = *argv++;
    argc--;

    while (argc--)
    {
        char buf[BUFSIZ];
        const char *argp = *argv++;

        if (strcmp(argp, "-help")==0)
        {
            usage();
        }
        else if (strcmp(argp, "-ns")==0)
        {
            print_frameid = false;
        }
        else if (strcmp(argp, "-q")==0)
        {
            quiet = true;
        }
        else if (strcmp(argp, "-b")==0)
        {
	  binary = true;
        }
        else if (strcmp(argp, "-i")==0)
        {
            // Input file name.
            if (argc>0)
            {
                // Next argument is input file name.
                input_fname = *argv++;
                argc--;
            }
            else
                usage("No input filename given.");
        }
        else if (strcmp(argp, "-i2")==0)
        {
            // Input file name.
            if (argc>0)
            {
                // Next argument is input file name.
                input_fname2 = *argv++;
                argc--;
            }
            else
                usage("No (2nd) input filename given.");
        }
        else if (strcmp(argp, "-o")==0)
        {
            // Output file name.
            if (argc>0)
            {
                // Next argument is output file name.
                output_fname = *argv++;
                argc--;
            }
            else
                usage("No output filename given.");
        }
        else if (strcmp(argp,"-iswap1") == 0)
        {
           iswap1 = true;
        }
        else if (strcmp(argp,"-iswap2") == 0)
        {
           iswap2 = true;
        }
        else if (strcmp(argp,"-oswap") == 0) 
        {
          oswap = true;
        }
        else if (strcmp(argp, "-debug")==0)
        {
            if (argc>0)
            {
                // Next argument is debug level.
                debug_level = (int) parse_long(*argv++);
                argc--;
            }
            else
                usage("No debug level given.");
        } 
        else if (strcmp(argp, "-sr")==0)
        {
            if (argc>0)
            {
	      sr_str = *argv++;
	      argc--;
            }
            else
                usage("No range given.");
        }
        else if (strcmp(argp, "-fr")==0)
        {
            if (argc>0)
            {
	      fr_str = *argv++;
	      argc--;
            }
            else
                usage("No range given.");
        }
        else if (strcmp(argp, "-pr")==0)
        {
            if (argc>0)
            {
	      pr_str = *argv++;
	      argc--;
            }
            else
                usage("No range given.");
        }
        else if (strcmp(argp, "-lr")==0)
        {
            if (argc>0)
            {
	      lr_str = *argv++;
	      argc--;
            }
            else
                usage("No range given.");
        }
        else {
	  sprintf(buf,"Unrecognized argument (%s).",argp);
	  usage(buf);
	}
    }

    //////////////////////////////////////////////////////////////////////
    // Check all necessary arguments provided before creating objects.
    //////////////////////////////////////////////////////////////////////


    SPI_base* in_streamp;
    SPO *out_streamp;
     
    if (input_fname == NULL) {
       usage("No input pfile name supplied");
     }

    if (input_fname2 == NULL) 
      in_streamp = new SPI(input_fname,iswap1);
    else
      in_streamp = new SPI2(input_fname,input_fname2,iswap1,iswap2);

    out_streamp = new SPI(output_fname,n_ftrs,n_labs,oswap);

    //////////////////////////////////////////////////////////////////////
    // Create objects.
    //////////////////////////////////////////////////////////////////////

    InFtrLabStream_PFile* in_streamp
        = new InFtrLabStream_PFile(debug_level, "", in_fp, 1,iswap1);



    sr_rng = new Range(sr_str,0,in_streamp->num_segs());
    fr_rng = new Range(fr_str,0,in_streamp->num_ftrs());
    lr_rng = new Range(lr_str,0,in_streamp->num_labs());

    //////////////////////////////////////////////////////////////////////
    // Do the work.
    //////////////////////////////////////////////////////////////////////

    pfile_process(*in_streamp,
		  *out_streamp,
		  *sr_rng,*fr_rng,*lr_rng,pr_str,
		  print_frameid,quiet,binary);

    //////////////////////////////////////////////////////////////////////
    // Clean up and exit.
    //////////////////////////////////////////////////////////////////////

    delete in_streamp;
    delete sr_rng;
    delete fr_rng;
    delete lr_rng;


    if (fclose(in_fp))
        error("Couldn't close input pfile.");
    if (fclose(out_fp))
        error("Couldn't close output file.");

    return EXIT_SUCCESS;
}
