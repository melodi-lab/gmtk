/*
    $Header$
  
    Finds if two observation files are different.  Written by Jeff
    Bilmes (bilmes@ee.washington.edu) and modified by Karim Filali
    (karim@cs.washington.edu) to handle different types of files.
  
*/

#include <stdlib.h>
#include <cstdio>
#include <cerrno>
#include <cstring>
#ifndef __CYGWIN__
#include <values.h>
#endif
#include <cmath>
#include <cassert>
#include "pfile.h"
#include "error.h"
#include "arguments.h"
#include "GMTK_WordOrganization.h"

#include "GMTK_ObsPrint.h"
#include "GMTK_ObsKLT.h"
#include "GMTK_ObsStats.h"
#include "GMTK_ObsNorm.h"
#include "GMTK_ObsGaussianNorm.h"

#ifdef DEBUG
#define DBGFPRINTF(_x_) fprintf _x_
#else
#define DBGFPRINTF(_x_)
#endif

#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))

ObservationMatrix Obs_Mat1,Obs_Mat2;

extern size_t bin_search(float *array,
			 size_t length, // the length of the array
			 float val);     // the value to search for.


static void obsDiff(ObservationMatrix * obs_mat1,
		    ObservationMatrix * obs_mat2,
		    Range& sr1rng,Range& sr2rng,
		    char *pr1_str, char *pr2_str,
		    const float tolerance,
		    const bool tolerance_percent,
		    bool stop)
{

    // Feature and label buffers are dynamically grown as needed.
    size_t buf_size = 300;      // Start with storage for 300 frames.
    size_t n1_labs = obs_mat1->numDiscrete();
    size_t n2_labs = obs_mat2->numDiscrete();
    size_t n1_ftrs = obs_mat1->numContinuous();
    size_t n2_ftrs = obs_mat2->numContinuous();

    if(n1_ftrs > n2_ftrs) {
      n1_ftrs = n2_ftrs;
    }
    else {
      n2_ftrs = n1_ftrs;
    }
    if(n1_labs > n2_labs) {
      n1_labs = n2_labs;
    }
    else {
      n2_labs = n1_labs;
    }

    float *ftr1_buf = new float[buf_size * n1_ftrs];
    float *ftr2_buf = new float[buf_size * n2_ftrs];
    float *ftr1_buf_p;
    float *ftr2_buf_p;

    UInt32* lab1_buf = new UInt32[buf_size * n1_labs];
    UInt32* lab2_buf = new UInt32[buf_size * n2_labs];
    UInt32* lab1_buf_p;
    UInt32* lab2_buf_p;

    size_t print_count = 0;


    // Go through input pfile to get the initial statistics,
    // i.e., max, min, mean, std, etc.
    Range::iterator srit1=sr1rng.begin();
    Range::iterator srit2=sr2rng.begin();
    for (;!srit1.at_end();srit1++,srit2++) {
       obs_mat1->loadSegment((const unsigned)(*srit1));
       const size_t n1_frames = obs_mat1->numFrames();

       obs_mat2->loadSegment((const unsigned)(*srit2));
       const size_t n2_frames = obs_mat2->numFrames();

	if (print_count ++ % 100 == 0) 
	  printf("Processing sentence %d and %d\n",(*srit1),(*srit2));

	Range pr1rng(pr1_str,0,n1_frames);
	Range pr2rng(pr2_str,0,n2_frames);

	if (pr1rng.length() != pr2rng.length()) {
	  error("Num frames of per-sentence ranges in sentences %d of file1 and %d of file2 different. pf1 = %d, pf2 = %d\n", 
		(*srit1),(*srit2),pr1rng.length(),pr2rng.length());
	}

        // Increase size of buffers if needed.
        if (n1_frames > buf_size || n2_frames > buf_size)
        {
            // Free old buffers.
            delete ftr1_buf;
            delete ftr2_buf;
            delete lab1_buf;
            delete lab2_buf;

            // Make twice as big to cut down on future reallocs.
            buf_size = MAX(n1_frames,n2_frames) * 2;

            // Allocate new larger buffers.
            ftr1_buf = new float[buf_size * n1_ftrs];
            ftr2_buf = new float[buf_size * n2_ftrs];
            lab1_buf = new UInt32[buf_size * n1_labs];
            lab2_buf = new UInt32[buf_size * n2_labs];

        }

	for(unsigned frame_no = 0;  frame_no < n1_frames; ++frame_no) {
	  const float* start_of_frame1 = obs_mat1->floatVecAtFrame(frame_no);
	  const float* start_of_frame2 = obs_mat2->floatVecAtFrame(frame_no);
	  const UInt32* start_of_unsigned_frame1 = obs_mat1->unsignedAtFrame(frame_no);
	  const UInt32* start_of_unsigned_frame2 = obs_mat2->unsignedAtFrame(frame_no);
	  for(unsigned feat_no = 0;  feat_no < n1_ftrs; ++feat_no) {
	    ftr1_buf[frame_no * n1_ftrs + feat_no] = *(start_of_frame1  + feat_no);
	    ftr2_buf[frame_no * n2_ftrs + feat_no] = *(start_of_frame2  + feat_no);
	  }
	  for(unsigned unsigned_feat_no = 0;  unsigned_feat_no < n1_labs; ++unsigned_feat_no) {
	    lab1_buf[frame_no*n1_labs + unsigned_feat_no] = *(start_of_unsigned_frame1+unsigned_feat_no);
	    lab2_buf[frame_no*n2_labs + unsigned_feat_no] = *(start_of_unsigned_frame2+unsigned_feat_no);
	  }
	}

	// construct the output set of features.

	ftr1_buf_p = ftr1_buf;
	ftr2_buf_p = ftr2_buf;

	for (unsigned prit=0; prit<n1_frames; ++prit) {

	  ftr1_buf_p = ftr1_buf + (prit)*n1_ftrs;
	  ftr2_buf_p = ftr2_buf + (prit)*n2_ftrs;

	  for (unsigned frit=0;frit<n1_ftrs; ++frit) {
	    const float diff = ftr1_buf_p[frit] - ftr2_buf_p[frit];
	    if (!tolerance_percent) {
	      if (fabs(diff) > tolerance) {
		printf("sent,pfrm,ftr of file1(%d,%d,%d) and file2(%d,%d,%d) diff > tolerance\n",(*srit1),(prit),(frit),(*srit2),(prit),(frit));
		if (stop)
		  goto done;
	      }
	    } else {
	      float dem = (fabs(ftr1_buf_p[frit]) + fabs(ftr2_buf_p[frit]));
	      if (dem == 0.0) {
		printf("Note: sent,pfrm,ftr of file1(%d,%d,%d) and file2(%d,%d,%d) are both zero\n",(*srit1),(prit),(frit),(*srit2),(prit),(frit));
	      } else if (200.0*fabs(diff)/dem > tolerance) {
		printf("sent,pfrm,ftr of file1(%d,%d,%d) and file2(%d,%d,%d) diff > tolerance\n",
		       (*srit1),(prit),(frit),(*srit2),(prit),(frit));
		if (stop)
		  goto done;
	      }
	    }
	  }
	  
	  lab1_buf_p = lab1_buf + (prit)*n1_labs;
	  lab2_buf_p = lab2_buf + (prit)*n2_labs;

	  for (unsigned lrit=0;lrit<n1_labs; ++lrit) {
	    if (lab1_buf_p[lrit] != lab2_buf_p[lrit]) {
	      printf("sent,pfrm lab of file1(%d,%d) and file2(%d,%d) differ.\n",
		     (*srit1),(prit),(*srit2),(prit));
	      if (stop)
		goto done;
	    }
	  }
	}
    }

done:

    delete ftr1_buf;
    delete ftr2_buf;
    delete lab1_buf;
    delete lab2_buf;
}



char *input_fname[2] = {NULL,NULL};  // Input pfile name.

const char * ifmtStr[2]={"pfile","pfile"};
unsigned ifmt[2];

unsigned int nis[2];
unsigned int nfs[2];

char *sr_str = NULL;
char *sr1_str = 0;   // sentence range string
Range *sr1_rng;
char *fr_str = 0;   // first feature range string
char *fr1_str = 0;   // first feature range string
Range *fr1_rng;
char *lr_str = 0;   // label range string    
char *lr1_str = 0;   // label range string    
Range *lr1_rng;

char *pr_str = 0;   // per-sentence range string
char *pr1_str = 0;   // per-sentence range string

char *sr2_str = 0;   // sentence range string
Range *sr2_rng;
char *fr2_str = 0;   // first feature range string
Range *fr2_rng;
char *lr2_str = 0;   // label range string    
Range *lr2_rng;
char *pr2_str = 0;   // per-sentence range string

float tolerance = 0.0;
bool tolerance_percent = false;
bool stop = false;
bool iswap1 = false, iswap2 = false;

char* perStreamTransformsBoth = NULL;   // 
char* perStreamTransforms[2] = {NULL,NULL};   // 
char* postTransforms                   = NULL;
char* postTransforms1                   = NULL;
char* postTransforms2                   = NULL;


bool quiet = false;
#ifdef INTV_WORDS_BIGENDIAN
bool iswap[2]={true,true};
bool oswap = true;
#else
bool iswap[2]= {false,false};
bool oswap             = false;
#endif 

bool     cppIfAscii        = true;
char*    cppCommandOptions = NULL;

int debug_level = 0;
bool help = false;

Arg Arg::Args[] = {
  Arg("i1",      Arg::Req, input_fname[0],"First input file"),
  Arg("i2",      Arg::Req, input_fname[1],"Second input file"),
  Arg("nf1",   Arg::Opt, nfs[0],"Number of floats in 1st input file"),
  Arg("nf2",   Arg::Opt, nfs[1],"Number of floats in 2nd input file"),
  Arg("ni1",   Arg::Opt, nis[0],"Number of ints (labels) in 1st input file"),
  Arg("ni2",   Arg::Opt, nis[1],"Number of ints (labels) in 2nd input file"),
  Arg("ifmt1", Arg::Opt,ifmtStr[0] ,"Format of 1st input file"),
  Arg("ifmt2", Arg::Opt,ifmtStr[1] ,"Format of 2nd input file"),
  Arg("iswap1",  Arg::Tog, iswap[0],"Byte-swap file 1"),
  Arg("iswap2",  Arg::Tog, iswap[1],"Byte-swap file 2"),
  Arg("sr",     Arg::Opt, sr_str,"Sentence range for both files"),
  Arg("sr1",     Arg::Opt, sr1_str,"Sentence range for 1st file"),
  Arg("sr2",     Arg::Opt, sr2_str,"sentence range for 2nd file"),
  Arg("fr",     Arg::Opt, fr_str,"Range of features selected.  Applies to both files"),
  Arg("fr1",     Arg::Opt, fr1_str,"Range of features selected from 1st file"),
  Arg("fr2",     Arg::Opt, fr2_str,"Range of features selected from 2nd file"),
  Arg("pr",     Arg::Opt, pr_str,"Per-sentence frame range for both files"),
  Arg("pr1",     Arg::Opt, pr1_str,"Per-sentence frame range for 1st file"),
  Arg("pr2",     Arg::Opt, pr2_str,"Per-sentence frame range for 2nd file"),
  Arg("lr",     Arg::Opt, lr_str,"Label range for both files"),
  Arg("lr1",     Arg::Opt, lr1_str,"Label range for 1st file"),
  Arg("lr2",     Arg::Opt, lr2_str,"Label range for 2nd file"),
  Arg("tol",     Arg::Opt, tolerance,"Tolerance absolute value"),
  Arg("tol%",    Arg::Opt, tolerance_percent,"Make tolarance a percentage [0-100]"),
  Arg("stop",     Arg::Opt, stop,"Do *not* continue after first difference found"),
  Arg("trans",  Arg::Opt,perStreamTransformsBoth ,"Transformations string for both files"),
  Arg("trans1",  Arg::Opt,perStreamTransforms[0] ,"Transformations string for 1st file"),
  Arg("trans2",  Arg::Opt,perStreamTransforms[1] ,"Transformations string for 2nd file"),
  Arg("posttrans",  Arg::Opt,postTransforms ,"Final global transformations string for both files"),
  Arg("posttrans1",  Arg::Opt,postTransforms1 ,"Final global transformations string for 1st file"),
  Arg("posttrans2",  Arg::Opt,postTransforms2 ,"Final global transformations string for 2nd file"),
  Arg("cppifascii",Arg::Opt, cppIfAscii,"Pre-process ASCII files using CPP"),
  Arg("cppCommandOptions",Arg::Opt,cppCommandOptions,"Additional CPP command line"),
  Arg("debug",     Arg::Opt, debug_level,"Number giving level of debugging output to produce 0=none"),
  Arg("help",     Arg::Opt, help,"Print this message"),
  Arg()
};

int main(int argc, const char *argv[]) {
    bool doWeSwap;
    
    ByteEndian byteEndian = getWordOrganization();
    switch(byteEndian) {
    case BYTE_BIG_ENDIAN:
      doWeSwap=false;
    break;
    case BYTE_LITTLE_ENDIAN:
    doWeSwap=true;
    break;
    default:
    // We weren't able to figure the Endian out.  Leave the swap defaults as they are.
#ifdef INTV_WORDS_BIGENDIAN
      doWeSwap=true;
#else
      doWeSwap=false;
#endif
    }
    
    for(int i=0; i<2; ++i) {
      iswap[i]=doWeSwap;
    }
  ///////////////////////////////////////////

  bool parse_was_ok = Arg::parse(argc,(char**)argv);

  if(help) {
    Arg::usage();
    exit(0);
  }

  if(!parse_was_ok) {
    Arg::usage(); exit(-1);
  }

  //////////////////////////////////////////////////////////////////////
  // Check all necessary arguments provided before creating objects.
  //////////////////////////////////////////////////////////////////////


    if (input_fname[0]==NULL)
        error("No first input pfile name supplied.");
    FILE *in1_fp = fopen(input_fname[0], "r");
    if (in1_fp==NULL)
        error("Couldn't open input pfile %s for reading.",input_fname[0]);

    if (input_fname[1]==0)
        error("No second input pfile name supplied.");
    FILE *in2_fp = fopen(input_fname[1], "r");
    if (in2_fp==NULL)
        error("Couldn't open input pfile %s for reading.",input_fname[1]);


 for (int i=0;i<2;i++) {
    if (strcmp(ifmtStr[i],"htk") == 0)
      ifmt[i] = HTK;
    else if (strcmp(ifmtStr[i],"binary") == 0)
      ifmt[i] = RAWBIN;
    else if (strcmp(ifmtStr[i],"ascii") == 0)
      ifmt[i] = RAWASC;
    else if (strcmp(ifmtStr[i],"pfile") == 0)
      ifmt[i] = PFILE;
    else
      error("ERROR: Unknown observation file format type: '%s'\n",ifmtStr[i]);
  }

    //////////////////////////////////////////////////////////////////////
    // Create objects.
    //////////////////////////////////////////////////////////////////////

    // If we have a pfile, we can extract the number if features from the file directly
 for(int i=0; i < 2; ++i) {
     if(ifmt[i]==PFILE) {
       FILE *in_fp = fopen(input_fname[i], "r");
       if (in_fp==NULL) error("Couldn't open input pfile for reading.");
       InFtrLabStream_PFile* in_streamp = new InFtrLabStream_PFile(debug_level,"",in_fp,1,iswap[i]);
       nis[i]=in_streamp->num_labs();
       nfs[i]=in_streamp->num_ftrs();
       if (fclose(in_fp)) error("Couldn't close input pfile.");
       delete in_streamp;
     }
     
     if(nis[i]==0 && nfs[i]==0) {
       error("The number of floats and the number of ints cannot be both zero.");
     }
   }
 
    if(fr_str != NULL) {
      fr1_str = fr2_str = fr_str;
    }

    if(lr_str != NULL) {
      lr1_str = lr2_str = lr_str;
    }

    if(pr_str != NULL) {
      pr1_str = pr2_str = fr_str;
    }

    if(sr_str != NULL) {
      sr1_str = sr2_str = sr_str;
    }

    if(perStreamTransformsBoth != NULL) {
      perStreamTransforms[0] = perStreamTransforms[1] = perStreamTransformsBoth;
    }

    if(postTransforms != NULL) {
      postTransforms1=postTransforms2=postTransforms;
    }



    Obs_Mat1.openFiles(1,  // number of files.   For now we use only one
		    (const char**)&input_fname[0],
		    (const char**)&fr1_str,
		    (const char**)&lr1_str,
		    (unsigned*)&nfs[0],
		    (unsigned*)&nis[0],
		    (unsigned*)&ifmt[0],
		    (bool*)&iswap[0],
		    0,  // startSkip
		    0,  // endSkip
		    cppIfAscii,
		    cppCommandOptions,
		    (const char**)&pr1_str,
		    NULL,//actionIfDiffNumFrames,
		    NULL,//actionIfDiffNumSents,
		    &perStreamTransforms[0],
		    postTransforms1);   


     sr1_rng = new Range(sr1_str,0,Obs_Mat1.numSegments());

Obs_Mat2.openFiles(1,  // number of files.   For now we use only one
		   (const char**)&input_fname[1],
		   (const char**)&fr2_str,
		   (const char**)&lr2_str,
		   (unsigned*)&nfs[1],
		   (unsigned*)&nis[1],
		   (unsigned*)&ifmt[1],
		   (bool*)&iswap[1],
		   0,  // startSkip
		   0,  // endSkip
		   cppIfAscii,
		   cppCommandOptions,
		   (const char**)&pr2_str,
		   NULL,//actionIfDiffNumFrames,
		   NULL,//actionIfDiffNumSents,
		   &perStreamTransforms[1],
		   postTransforms2);   


     sr2_rng = new Range(sr2_str,0,Obs_Mat2.numSegments());

    //////////////////////////////////////////////////////////////////////
    // Do the work.
    //////////////////////////////////////////////////////////////////////

    obsDiff(&Obs_Mat1,
	    &Obs_Mat2,
	    *sr1_rng,*sr2_rng,
	    pr1_str,pr2_str,
	    tolerance,
	    tolerance_percent,
	    stop);

    //////////////////////////////////////////////////////////////////////
    // Clean up and exit.
    //////////////////////////////////////////////////////////////////////


    delete sr1_rng;
    delete fr1_rng;
    delete lr1_rng;
    delete sr2_rng;
    delete fr2_rng;
    delete lr2_rng;

    if (fclose(in1_fp))
        error("Couldn't close input1 pfile.");
    if (fclose(in2_fp))
        error("Couldn't close input2 pfile.");

    return EXIT_SUCCESS;
}
