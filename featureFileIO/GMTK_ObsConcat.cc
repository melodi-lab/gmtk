/*
    $Header$
  
    This program reads one or more input files (in pfile, htk, binary,
    or ascii formats) and concatenates them "longitudenaly"
    i.e. writes ou their frames end-to-end, the feature dimension
    staying fixed.  The program is based on Dan Ellis's pfile_concat
    (1998may07 dpwe@icsi.berkeley.edu).

    2003-12-20 karim@cs.washington.edu

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
//#include "parse_subset.h"
#include "error.h"
#include "arguments.h"
#include "GMTK_WordOrganization.h"

#include "GMTK_ObsPrint.h"
#include "GMTK_ObsKLT.h"
#include "GMTK_ObsStats.h"
#include "GMTK_ObsNorm.h"
#include "GMTK_ObsGaussianNorm.h"

#include "range.h"

#ifdef DEBUG
#define DBGFPRINTF(_x_) fprintf _x_
#else
#define DBGFPRINTF(_x_)
#endif


char *output_fname = 0; // Output pfile name.
char * outputList = NULL;
FILE * outputListFp=NULL;
char * outputNameSeparatorStr="_";


void (*copy_swap_func_ptr)(size_t, const intv_int32_t*, intv_int32_t*)=NULL;

ObservationMatrix globalObservationMatrix;

void printHTKHeader(FILE* ofp, bool oswap, int numInts, int numFloats, int numSamples) {

  DBGFPRINTF((stderr,"obsPrint: Printing HTK header.\n"));
  // structure of HTK header
  //  Int32 numSamples;
  Int32 samplePeriod=1;
  short parameterKind;
  short sampleSize;
  if(numFloats > 0) {
    sampleSize=numFloats*sizeof(float);
    parameterKind=USER;  // enum in GMTK_Stream.h
  }
  else {
    sampleSize=numInts*sizeof(short);
    parameterKind=DISCRETE;  // enum in GMTK_Stream.h
  }

  if (oswap) {
    numSamples = swapb_i32_i32(numSamples);
    samplePeriod = swapb_i32_i32(samplePeriod);
    sampleSize = swapb_short_short(sampleSize);
    parameterKind = swapb_short_short(parameterKind);
  }

  if (fwrite(&numSamples,sizeof(Int32),1,ofp) != 1) {
    error("Cannot write HTK number of samples\n");
  }
  if (fwrite((short *)&samplePeriod,sizeof(Int32),1,ofp) != 1) {
    error("Cannot write HTK sample period\n");
  }

  if (fwrite((short *)&sampleSize,sizeof(short),1,ofp) != 1) {
    error("Cannot write HTK sample size\n");
  }

  if (fwrite(&parameterKind,sizeof(short),1,ofp) != 1) {
    error("Cannot write HTK parm kind\n");
  }

  DBGFPRINTF((stderr,"obsPrint: Finished printing HTK header.\n"));

}

void printSegment(unsigned sent_no, FILE* out_fp, float* cont_buf, unsigned num_continuous, UInt32* disc_buf, unsigned num_discrete, unsigned num_frames, const bool dontPrintFrameID,const bool quiet,unsigned ofmt,int debug_level,bool oswap, OutFtrLabStream_PFile* out_stream) {


    if(ofmt==RAWASC || ofmt==RAWBIN || ofmt==HTK) {
      char* current_output_fname = new char[strlen(output_fname)+strlen(outputNameSeparatorStr)+50];
      sprintf(current_output_fname,"%s%s%d",output_fname,outputNameSeparatorStr,sent_no);
      if ((out_fp = fopen(current_output_fname, "w")) == NULL) {
	error("Couldn't open output file (%s) for writing.",current_output_fname);
      }

      if(outputListFp != NULL) {
	fprintf(outputListFp,"%s\n",current_output_fname);
      }

      delete []  current_output_fname;
    }
    
    if(ofmt==HTK) {
      DBGFPRINTF((stderr,"obsPrint: Calling printHTKHeader with numSamples=%d.\n",num_frames));
      printHTKHeader(out_fp,oswap,num_discrete,num_continuous,num_frames);
    }

    float* cont_buf_p = cont_buf;
    UInt32* disc_buf_p = disc_buf;
    //////////  Print the frames ////////////////////////////////////
    for (unsigned frame_no=0; frame_no < num_frames ; ++frame_no) {
      bool ns = false;
      if (!dontPrintFrameID) {
	if (ofmt==FLATBIN || ofmt==RAWBIN) {
	  copy_swap_func_ptr(1,(int*)&sent_no,(int*)&sent_no);
	  copy_swap_func_ptr(1,(int*)&frame_no,(int*)&frame_no);
	  fwrite(&sent_no,sizeof(sent_no),1,out_fp);
	  fwrite(&frame_no,sizeof(frame_no),1,out_fp);
	} else if(ofmt==FLATASC || ofmt==RAWASC ){
	  fprintf(out_fp,"%d %u",sent_no,frame_no);
	  ns = true;
	}
      }
      
      /// Print continuous part of frame /////////////////////////////
      for (unsigned frit=0;frit< num_continuous; ++frit) {
	if (ofmt==FLATBIN || ofmt==RAWBIN || ofmt==HTK) {
	  DBGFPRINTF((stderr,"obsPrint: Printing HTK float %f.\n",cont_buf_p[frit]));
	  copy_swap_func_ptr(1,(int*)&cont_buf_p[frit],(int*)&cont_buf_p[frit]);
	  fwrite(&cont_buf_p[frit], sizeof(cont_buf_p[frit]),  1,out_fp);
	} 
	else if(ofmt==FLATASC || ofmt==RAWASC){
	  if (ns) fprintf(out_fp," ");
	  fprintf(out_fp,"%f",cont_buf_p[frit]);
	  ns = true;
	}
      }
      ///////////////////////////////////////////////////////////////

      /// Print discrete part of the frame ///////////////////////////
      for (unsigned lrit=0;lrit<num_discrete; ++lrit) {
	if (ofmt==FLATBIN || ofmt==RAWBIN || (ofmt==HTK && num_continuous>0) ) {
	  copy_swap_func_ptr(1,(int*)&disc_buf_p[lrit],(int*)&disc_buf_p[lrit]);
	  fwrite(&disc_buf_p[lrit],  sizeof(disc_buf_p[lrit]), 1,out_fp);
	} 
	else if(ofmt==HTK && num_continuous==0) { // in the HTK format we
  // cannot mix floats with discrete data; that's why if there is at
  // least one float component everyting is written out as a float.
	  short short_lab_buf_p=(short)disc_buf_p[lrit];
	  DBGFPRINTF((stderr,"obsPrint: Printing HTK short %d.\n",short_lab_buf_p));
	   if (oswap) {
	     short_lab_buf_p = swapb_short_short(short_lab_buf_p);
	   }
	  fwrite(&short_lab_buf_p,  sizeof(short_lab_buf_p), 1,out_fp);
	}
	else if(ofmt==FLATASC || ofmt==RAWASC) {
	  if (ns) fprintf(out_fp," ");	    
	  fprintf(out_fp,"%d",disc_buf_p[lrit]);
	  ns = true;
	}
      } // end of for (unsigned lrit=0;lrit<num_discrete; ++lrit)
      ////////////////////////////////////////////////////////////////

      if (ofmt==FLATASC || ofmt==RAWASC)
	fprintf(out_fp,"\n");

      cont_buf_p += num_continuous;
      disc_buf_p += num_discrete;
    }  // end of for (unsigned frame_no=0; frame_no < num_frames ; ++frame_no)
    ////////////////////////////////////////////////////////////////////////////

    if(ofmt==PFILE) {
      out_stream->write_ftrslabs(num_frames, cont_buf, disc_buf);
      out_stream->doneseg((SegID) sent_no);
    }
    
    if(ofmt==RAWASC || ofmt==RAWBIN || ofmt==HTK) {
      if (fclose(out_fp)) error("Couldn't close output file.");
    }
}


static void obsConcat(FILE *out_fp, 
		      char *sr_str[], 
		      char *fr_str[], 
		      char *lr_str[], 
		      char *spr_str[], 
		      char *input_fname[],
		      int n_input_fnames, 
		      unsigned nfs[],
		      unsigned nis[],
		      unsigned ifmt[],
		      bool iswap[],
		      unsigned ofmt,
		      unsigned startSkip,
		      unsigned endSkip,
		      bool cppIfAscii,
		      char* cppCommandOptions,
		      char* perStreamTransforms[],
		      char* postTransforms,
		      bool dontPrintFrameID,
		      const int debug_level, 
		      const bool quiet,
		      bool oswap)

{
    // Feature and label buffers are dynamically grown as needed.
    size_t buf_size = 300;      // Start with storage for 300 frames.
    char errmsg[1024];

    // pre-open the first input file so we can check its size
    const char *file_name;
    assert(n_input_fnames > 0);
    file_name = input_fname[0];
    FILE *in_fp = fopen(file_name, "r");
    if (in_fp==NULL) {
	sprintf(errmsg, "Couldn't open input pfile %s for reading.",file_name);
	error(errmsg);
    }

    globalObservationMatrix.openFiles(1,  // number of files
				      (const char**)&file_name,
				      (const char**)&fr_str,
				      (const char**)&lr_str,
				      (unsigned*)&nfs,
				      (unsigned*)&nis,
				      (unsigned*)&ifmt,
				      (bool*)&iswap,
				      startSkip,
				      endSkip,
				      cppIfAscii,
				      cppCommandOptions,
				      (const char**)&spr_str,
				      NULL, // actionIfDiffNumFrames,
				      NULL, //actionIfDiffNumSents,
				      perStreamTransforms,
				      postTransforms);

    unsigned n_ftrs = globalObservationMatrix.numContinuous();
    unsigned n_labs = globalObservationMatrix.numDiscrete();
    
    OutFtrLabStream_PFile* out_stream=NULL;
    if(ofmt==PFILE) {
      out_stream = new OutFtrLabStream_PFile(debug_level,"",out_fp,n_ftrs,n_labs,1,oswap);
    }

    float *ftr_buf = new float[buf_size * n_ftrs];
    UInt32 *lab_buf = new UInt32[buf_size * n_labs];

    // Outer loop around input filenames
    unsigned global_seg_num=0;
    for (int file_ix = 0; file_ix < n_input_fnames; ++file_ix) {
	file_name = input_fname[file_ix];
	// in_fp actually already open first time into the loop (i.e. for file_ix==0)
	if (in_fp == NULL) {
	    in_fp = fopen(file_name, "r");
	    if (in_fp==NULL) {
		error("Couldn't open input file %s for reading.",file_name);
	    }

	    
	    globalObservationMatrix.openFiles(1,  // number of files
					      (const char**)&file_name,
					      (const char**)&fr_str[file_ix],
					      (const char**)&lr_str[file_ix],
					      (unsigned*)&nfs[file_ix],
					      (unsigned*)&nis[file_ix],
					      (unsigned*)&ifmt[file_ix],
					      (bool*)&iswap[file_ix],
					      startSkip,
					      endSkip,
					      cppIfAscii,
					      cppCommandOptions,
					      (const char**)&spr_str[file_ix],
					      NULL, // actionIfDiffNumFrames,
					      NULL, //actionIfDiffNumSents,
					      &perStreamTransforms[file_ix],
					      postTransforms);
	    
	    unsigned check_n_ftrs = globalObservationMatrix.numContinuous();
	    unsigned check_n_labs = globalObservationMatrix.numDiscrete();


	    // Check that this input pfile is the same size as predecessors
	    if(check_n_ftrs != n_ftrs || check_n_labs != n_labs) {
	      sprintf(errmsg, "Features/labels of %s (%d/%d) don't match first input file (%d/%d)", file_name, check_n_labs, check_n_ftrs, n_labs, n_ftrs);
	      error(errmsg);
	    }
	}
	
	// create the sentence-range iterator specifically for this 
	// input file
	Range sr_rng(sr_str[file_ix],0,globalObservationMatrix.numSegments());
	
	for (Range::iterator srit=sr_rng.begin();!srit.at_end();srit++,global_seg_num++) {
	  globalObservationMatrix.loadSegment(*srit);
	  const size_t n_frames = globalObservationMatrix.numFrames();
	  
	  if (!quiet && (*srit) % 100 == 0)
	    printf("Processing sentence %d of file %s\n",
		   (*srit), file_name);
	  
	  // Increase size of buffers if needed.
	  if (n_frames > buf_size){
	    // Free old buffers.
	    delete lab_buf;
	    delete ftr_buf;
	    // Make twice as big to cut down on future reallocs.
	    buf_size = n_frames * 2;
	      // Allocate new larger buffers.
	    ftr_buf = new float[buf_size * n_ftrs];
	    lab_buf = new UInt32[buf_size * n_labs];
	  }

	  for(unsigned frame_no = 0;  frame_no < n_frames; ++frame_no) {
	    const float* start_of_frame = globalObservationMatrix.floatVecAtFrame(frame_no);
	    const UInt32* start_of_unsigned_frame = globalObservationMatrix.unsignedAtFrame(frame_no);
	    for(unsigned feat_no = 0;  feat_no < n_ftrs; ++feat_no) {
	      ftr_buf[frame_no * n_ftrs + feat_no] = *(start_of_frame  + feat_no);
	    }
	    for(unsigned unsigned_feat_no = 0;  unsigned_feat_no < n_labs; ++unsigned_feat_no) {
	      lab_buf[frame_no*n_labs + unsigned_feat_no] = *(start_of_unsigned_frame+unsigned_feat_no);
	    }
	      
	    // Write output.
	    printSegment(global_seg_num, out_fp, ftr_buf,n_ftrs,lab_buf,n_labs,n_frames, dontPrintFrameID,quiet, ofmt, debug_level, oswap, out_stream);
	  }
	  
	  fclose(in_fp);
	  // make sure top of loop knows to open next file
	  in_fp = NULL;
	}
    }
    
    // All done; close output
    // must delete pfile object first, because it needs to rewrite
    if(ofmt==PFILE)
      delete out_stream;
    if (fclose(out_fp))
      error("Couldn't close output file.");
	
    delete lab_buf;
    delete ftr_buf;
}

#define MAX_OBJECTS 5

char *input_fname[MAX_OBJECTS] = {NULL,NULL,NULL,NULL,NULL};  // Input file name.
char * ifmtStr[MAX_OBJECTS]={"pfile","pfile","pfile","pfile","pfile"};
unsigned ifmt[MAX_OBJECTS];



char * ofmtStr="flatasc";
unsigned ofmt;

unsigned int nis[MAX_OBJECTS];
unsigned int nfs[MAX_OBJECTS];

//char  *sr_str               = 0;   // sentence range string
char  *sr_str[MAX_OBJECTS]  = {NULL,NULL,NULL,NULL,NULL};
char  *fr_str[MAX_OBJECTS]  = {NULL,NULL,NULL,NULL,NULL};   // feature range string    
char  *lr_str[MAX_OBJECTS]  = {NULL,NULL,NULL,NULL,NULL};   // label range string  
char  *spr_str[MAX_OBJECTS] = {NULL,NULL,NULL,NULL,NULL};   // per stream per sentence range string 

int  debug_level = 0;
bool dontPrintFrameID = false;
bool quiet = false;
#ifdef INTV_WORDS_BIGENDIAN
bool iswap[MAX_OBJECTS]={true,true,true,true,true};
bool oswap = true;
#else
bool iswap[MAX_OBJECTS]= {false,false,false,false,false};
bool oswap             = false;
#endif 

bool     cppIfAscii        = true;
char*    cppCommandOptions = NULL;

unsigned startSkip = 0;
unsigned endSkip   = 0;

char* perStreamTransforms[MAX_OBJECTS] = {NULL,NULL,NULL,NULL,NULL};   // 
char* postTransforms                   = NULL;

bool help=false;

Arg Arg::Args[] = {
  Arg("i",      Arg::Req, input_fname,"input file",Arg::ARRAY,MAX_OBJECTS),
  Arg("ifmt",      Arg::Opt,ifmtStr ,"format of input file",Arg::ARRAY,MAX_OBJECTS),
  Arg("o",      Arg::Opt, output_fname,"output file"),
  Arg("ofmt",      Arg::Opt, ofmtStr,"format of output file"),
  Arg("nf",   Arg::Opt, nfs,"number of floats in input file",Arg::ARRAY,MAX_OBJECTS),
  Arg("ni",   Arg::Opt, nis,"number of ints (labels) in input file",Arg::ARRAY,MAX_OBJECTS),
  Arg("q",  Arg::Tog, quiet,"quiet mode"),
  Arg("sr",   Arg::Opt, sr_str,"sentence range"),
  Arg("fr",   Arg::Opt, fr_str,"feature range",Arg::ARRAY,MAX_OBJECTS),
  Arg("lr",   Arg::Opt, lr_str,"label range",Arg::ARRAY,MAX_OBJECTS),
  Arg("spr",  Arg::Opt, spr_str,"per stream per-sentence range",Arg::ARRAY,MAX_OBJECTS),
  Arg("startskip",   Arg::Opt, startSkip,"start skip"),
  Arg("endskip",   Arg::Opt, endSkip,"end skip"),
  Arg("iswap",Arg::Opt, iswap,"do byte swapping on the input file",Arg::ARRAY,MAX_OBJECTS),
  Arg("oswap",Arg::Opt, oswap,"do byte swapping on the output file"),
  Arg("cppifascii",Arg::Opt, cppIfAscii,"Pre-process ASCII files using CPP"),
  Arg("ns",    Arg::Opt, dontPrintFrameID,"Don't print the frame IDs (i.e., sent and frame #)"),
  Arg("sep",      Arg::Opt, outputNameSeparatorStr,"String to use as separator when outputting raw ascii or binary files (one sentence per file)."),
  Arg("olist",      Arg::Opt, outputList,"output list-of-files name.  Only meaningful if used with the RAW or HTK formats."),
  Arg("cppCommandOptions",Arg::Opt,cppCommandOptions,"Additional CPP command line"),
  Arg("trans",  Arg::Opt,perStreamTransforms ,"per stream transformations string",Arg::ARRAY,MAX_OBJECTS),
  Arg("posttrans",  Arg::Opt,postTransforms ,"Final global transformations string"),
  Arg("debug",  Arg::Opt, debug_level,"number giving level of debugging output to produce 0=none"),
  Arg("help",   Arg::Tog, help,"print this message"),
  // The argumentless argument marks the end of the above list.
  Arg()
};



int main(int argc, const char *argv[]) {

    //////////////////////////////////////////////////////////////////////
    // Check all necessary arguments provided before creating objects.
    //////////////////////////////////////////////////////////////////////
      
    int numFiles=0;

    // Figure out the Endian of the machine this is running on and set the swap defaults accordingly
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

  oswap=doWeSwap;
  for(int i=0; i<MAX_OBJECTS; ++i) {
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


  // oswap might been assigned a new value on the command line
  if(oswap) {
    copy_swap_func_ptr=&swapb_vi32_vi32;
  }
  else {
    copy_swap_func_ptr=&copy_vi32_vi32;
  }

  //////////////////////////////////////////////////////////////////////
  // Check all necessary arguments provided before creating objects.
  //////////////////////////////////////////////////////////////////////

 for (int i=0;i<MAX_OBJECTS;i++) {
    numFiles += (input_fname[i] != NULL);
    if (strcmp(ifmtStr[i],"htk") == 0)
      ifmt[i] = HTK;
    else if (strcmp(ifmtStr[i],"binary") == 0)
      ifmt[i] = RAWBIN;
    else if (strcmp(ifmtStr[i],"ascii") == 0)
      ifmt[i] = RAWASC;
    else if (strcmp(ifmtStr[i],"pfile") == 0)
      ifmt[i] = PFILE;
    //    else if (strcmp(ifmtStr[i],"flatbin") == 0)
    //  ifmt[i]=FLATBIN;
    //else if (strcmp(ifmtStr[i],"flatasc") == 0)
    //  ifmt[i]=FLATASC;
    else
      error("ERROR: Unknown observation file format type: '%s'\n",ifmtStr[i]);
 }
 
 for(int i=0; i < numFiles; ++i) {
   // first check that we don't have gaps
   if (input_fname[i] == NULL) {
     error("ERROR: Cannot skip an observation file number.");
   }

   if(output_fname!=NULL && strcmp(input_fname[i],output_fname)==0) {
     error("Input and output filenames cannot be the same.");
   }
 }
 if (strcmp(ofmtStr,"htk") == 0)
   ofmt = HTK;
 else if (strcmp(ofmtStr,"binary") == 0)
   ofmt = RAWBIN;
 else if (strcmp(ofmtStr,"ascii") == 0)
   ofmt = RAWASC;
 else if (strcmp(ofmtStr,"pfile") == 0)
   ofmt = PFILE;
 else if (strcmp(ofmtStr,"flatbin") == 0)
   ofmt=FLATBIN;
 else if (strcmp(ofmtStr,"flatasc") == 0)
   ofmt=FLATASC;
 else
   error("ERROR: Unknown observation file format type: '%s'\n",ofmtStr);


 FILE *out_fp=NULL;
 if (output_fname==0 || !strcmp(output_fname,"-")) {
   if(ofmt==RAWASC || ofmt==RAWBIN || ofmt==HTK) {
     error("Need to specify output filename when output type is RAWBIN, RAWASC, or HTK.\n");
   }
   out_fp = stdout;
 }
 else {
   if(ofmt != RAWASC && ofmt != RAWBIN && ofmt != HTK) {
     if ((out_fp = fopen(output_fname, "w")) == NULL) {
       error("Couldn't open output file for writing.\n");
     }
   }
 }
 

 if(outputList != NULL) {
   if ((outputListFp = fopen(outputList, "w")) == NULL) {
     error("Couldn't open output list (%s) for writing.\n",outputList);
   }
 }
 
 //////////////////////////////////////////////////////////////////////
 // Create objects.
 //////////////////////////////////////////////////////////////////////
 // If we have a pfile, we can extract the number if features from the file directly
 for(int i=0; i < numFiles; ++i) {
   if(input_fname[i]!=NULL) {
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
 }
 
 //////////////////////////////////////////////////////////////////////
    // Create objects.
    //////////////////////////////////////////////////////////////////////
    
    // done in subroutine now, after we open the first input file

    //////////////////////////////////////////////////////////////////////
    // Do the work.
    //////////////////////////////////////////////////////////////////////

    obsConcat(out_fp, 
	      sr_str, 
	      fr_str, 
	      lr_str, 
	      spr_str, 
	      input_fname, 
	      numFiles,
	      nfs,
	      nis,
	      ifmt,
	      iswap,
	      ofmt,
	      startSkip,
	      endSkip,
	      cppIfAscii,
	      cppCommandOptions,
	      perStreamTransforms,
	      postTransforms,
	      dontPrintFrameID,
	      debug_level, 
	      quiet,
	      oswap);

    //////////////////////////////////////////////////////////////////////
    // Clean up and exit.
    //////////////////////////////////////////////////////////////////////

    return EXIT_SUCCESS;
}
