/*
    $Header$
  
    This program reads one or more input files (in pfile, htk, binary,
    or ascii formats) and concatenates them "longitudenaly"
    i.e. writes ou their frames end-to-end, the feature dimension
    staying fixed.  The program is based on Dan Ellis's pfile_concat
    (1998may07 dpwe@icsi.berkeley.edu).

    2003-12-20 karim@cs.washington.edu

*/


#ifdef HAVE_CONFIG_H

#include <config.h>
static const char * gmtk_version_id = PACKAGE_STRING;
#ifdef HAVE_HG_H
#include "hgstamp.h"
#endif

#else 
// TODO: automate the process of updating this string.
static const char * gmtk_version_id = "GMTK Version 0.2b Tue Jan 20 22:59:41 2004";
#endif



#include <stdlib.h>
#include <cstdio>
#include <cerrno>
#include <cstring>
#include <limits.h>
#include <float.h>
#include <cmath>
#include <cassert>
#include "pfile.h"
//#include "parse_subset.h"
#include "debug.h"
#include "error.h"
#include "arguments.h"
#include "GMTK_WordOrganization.h"

#include "GMTK_ObsPrint.h"
#include "GMTK_ObsKLT.h"
#include "GMTK_ObsStats.h"
#include "GMTK_ObsNorm.h"
#include "GMTK_ObsGaussianNorm.h"

#include "GMTK_ObservationFile.h"
#include "GMTK_FilterFile.h"
#include "GMTK_FileSource.h"
#include "GMTK_FileSourceNoCache.h"
#include "GMTK_ASCIIFile.h"
#include "GMTK_FlatASCIIFile.h"
#include "GMTK_PFileFile.h"
#include "GMTK_HTKFile.h"
#include "GMTK_HDF5File.h"
#include "GMTK_BinaryFile.h"
#include "GMTK_Filter.h"
#include "GMTK_ObservationArguments.h"
#include "GMTK_FIRFilter.h"
#include "GMTK_ARMAFilter.h"
#include "GMTK_AffineFilter.h"
#include "GMTK_UpsampleFilter.h"
#include "GMTK_UpsampleSmoothFilter.h"

#include "range.h"
#include "vbyteswapping.h"

#ifdef DEBUG
#define DBGFPRINTF(_x_) fprintf _x_
#else
#define DBGFPRINTF(_x_)
#endif



#define GMTK_ARG_OBS_FILES
#define GMTK_ARG_CPP_CMD_OPTS
#define GMTK_ARG_OBS_MATRIX_XFORMATION
#define GMTK_ARG_FILE_RANGE_OPTIONS
#define GMTK_ARG_START_END_SKIP
#define GMTK_ARG_HELP
#define GMTK_ARG_VERSION

#define GMTK_ARGUMENTS_DEFINITION
#include "ObsArguments.h"
#undef GMTK_ARGUMENTS_DEFINITION


char *output_fname = NULL; // Output pfile name.
char * outputList = NULL;
FILE * outputListFp=NULL;
const char * outputNameSeparatorStr="_";


void (*copy_swap_func_ptr)(size_t, const intv_int32_t*, intv_int32_t*)=NULL;

#if 0
ObservationMatrix globalObservationMatrix;
#else
FileSource *gomFS;
#endif

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
  if (fwrite((Int32 *)&samplePeriod,sizeof(Int32),1,ofp) != 1) {
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
    size_t fwrite_result;
    //////////  Print the frames ////////////////////////////////////
    for (unsigned frame_no=0; frame_no < num_frames ; ++frame_no) {
      bool ns = false;
      if (!dontPrintFrameID) {
	if (ofmt==FLATBIN || ofmt==RAWBIN) {
	  copy_swap_func_ptr(1,(int*)&sent_no,(int*)&sent_no);
	  copy_swap_func_ptr(1,(int*)&frame_no,(int*)&frame_no);
	  fwrite_result = fwrite(&sent_no,sizeof(sent_no),1,out_fp);
          if (fwrite_result != 1) {
            error("Error writing to output file");
          }
	  fwrite_result = fwrite(&frame_no,sizeof(frame_no),1,out_fp);
          if (fwrite_result != 1) {
            error("Error writing to output file");
          }
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
	  fwrite_result = fwrite(&cont_buf_p[frit], sizeof(cont_buf_p[frit]),  1,out_fp);
          if (fwrite_result != 1) {
            error("Error writing to output file");
          }
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
	  fwrite_result = fwrite(&disc_buf_p[lrit],  sizeof(disc_buf_p[lrit]), 1,out_fp);
          if (fwrite_result != 1) {
            error("Error writing to output file");
          }
	} 
	else if(ofmt==HTK && num_continuous==0) { // in the HTK format we
  // cannot mix floats with discrete data; that's why if there is at
  // least one float component everyting is written out as a float.
	  short short_lab_buf_p=(short)disc_buf_p[lrit];
	  DBGFPRINTF((stderr,"obsPrint: Printing HTK short %d.\n",short_lab_buf_p));
	   if (oswap) {
	     short_lab_buf_p = swapb_short_short(short_lab_buf_p);
	   }
	  fwrite_result = fwrite(&short_lab_buf_p,  sizeof(short_lab_buf_p), 1,out_fp);
          if (fwrite_result != 1) {
            error("Error writing to output file");
          }
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


#ifndef MEBIBYTE
#define MEBIBYTE (1048576)  
#endif

#define ALLOREMPTY(var) \
  ( ( strncasecmp(var, "all",  4) == 0 ) || \
    ( strncasecmp(var, "nil",  4) == 0 ) || \
    ( strncasecmp(var, "none", 5) == 0 ) || \
    ( strncasecmp(var, "full", 5) == 0 ) || \
    ( strlen(var) == 0 ) )


FileSource *
openOneSourceFile(char *ofs,
		  unsigned ifmt,
		  unsigned nfs,
		  unsigned nis,
		  bool iswp,
		  bool Cpp_If_Ascii,
		  char *cppCommandOptions,
		  const char *prefrs,
		  const char *preirs,
		  const char *prepr,
		  const char *sr,
		  char *Per_Stream_Transforms,
		  const char *frs,
		  const char *irs,
		  char *postpr,
		  char *Post_Transforms,
		  char *gpr_str,
		  bool constantSpace,
		  int startSkip,
		  int endSkip,
		  unsigned fileBufferSize,
		  unsigned fileWindowSize,
		  unsigned fileWindowDelta,
		  unsigned justification,
		  unsigned fileNum)
{
  // range selection is much more efficient if "all" is replaced with NULL
  // since the logical <-> physical mapping step can be skipped
  if (    frs && ALLOREMPTY(frs))         frs = NULL;
  if (    irs && ALLOREMPTY(irs))         irs = NULL;
  if ( prefrs && ALLOREMPTY(prefrs))   prefrs = NULL;
  if ( preirs && ALLOREMPTY(preirs))   preirs = NULL;
  if (     sr && ALLOREMPTY(sr))           sr = NULL;
  if (  prepr && ALLOREMPTY(prepr))     prepr = NULL;
  if ( postpr && ALLOREMPTY(postpr))   postpr = NULL;
  if (gpr_str && ALLOREMPTY(gpr_str)) gpr_str = NULL;

  ObservationFile *obsFile;
  
  obsFile = instantiateFile(ifmt, ofs, nfs, nis, fileNum, iswp,
			    Cpp_If_Ascii, cppCommandOptions, prefrs, preirs,
			    prepr, sr);
  assert(obsFile);
  if (Per_Stream_Transforms || frs || irs || postpr) {
    Filter *fileFilter = instantiateFilters(Per_Stream_Transforms,
					    obsFile->numContinuous(),
					    obsFile->numDiscrete());
    assert(fileFilter);
    obsFile = new FilterFile(fileFilter, obsFile, frs, irs, postpr);
  }
  ObservationFile *mf = obsFile;

  ObservationFile *ff;
  if (Post_Transforms || gpr_str) {
    ff = new FilterFile(instantiateFilters(Post_Transforms, 
                                           obsFile->numLogicalContinuous(),
			                   obsFile->numLogicalDiscrete()),
                        mf, NULL, NULL, gpr_str);
  } else {
    ff = mf;
  }
  unsigned windowBytes = fileWindowSize * MEBIBYTE;
  infoMsg(IM::ObsFile, IM::Low, "windowBytes = %u MiB = %u B\n", fileWindowSize, windowBytes);
  infoMsg(IM::ObsFile, IM::Low, "fileBufferSize = %u\n", fileBufferSize);
  if (constantSpace) {
    return new FileSource(ff, windowBytes, fileWindowDelta, fileBufferSize, 
                          startSkip, endSkip, justification, constantSpace);
  } else {
    return new FileSourceNoCache(ff, windowBytes, fileWindowDelta, fileBufferSize, 
                                 startSkip, endSkip, justification);
  }
}


static void obsConcat(FILE *out_fp, 
		      int n_input_fnames, 
		      unsigned ofmt,
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
  file_name = ofs[0];
  gomFS = openOneSourceFile(ofs[0], ifmts[0], nfs[0], nis[0], iswp[0], 
			    Cpp_If_Ascii, cppCommandOptions,
			    prefrs[0], preirs[0], prepr[0], sr[0],
			    Per_Stream_Transforms[0], frs[0], irs[0],
			    postpr[0], Post_Transforms, gpr_str,
			    constantSpace, startSkip, endSkip,
			    fileBufferSize, fileWindowSize, fileWindowDelta, 
			    justification, 0);
  unsigned n_ftrs = gomFS->numContinuous();
  unsigned n_labs = gomFS->numDiscrete();
    
  OutFtrLabStream_PFile* out_stream=NULL;
  if(ofmt==PFILE) {
    out_stream = new OutFtrLabStream_PFile(debug_level,"",out_fp,n_ftrs,n_labs,1,oswap);
  }

  float *ftr_buf = new float[buf_size * n_ftrs];
  UInt32 *lab_buf = new UInt32[buf_size * n_labs];

  // Outer loop around input filenames
  unsigned global_seg_num=0;
  for (int file_ix = 0; file_ix < n_input_fnames; ++file_ix) {
    file_name = ofs[file_ix];
    delete gomFS;
    gomFS = openOneSourceFile(ofs[file_ix], ifmts[file_ix], nfs[file_ix], nis[file_ix], iswp[file_ix], 
			      Cpp_If_Ascii, cppCommandOptions,
			      prefrs[file_ix], preirs[file_ix], prepr[file_ix], sr[file_ix],
			      Per_Stream_Transforms[file_ix], frs[file_ix], irs[file_ix],
			      postpr[file_ix], Post_Transforms, gpr_str,
			      constantSpace, startSkip, endSkip,
			      fileBufferSize, fileWindowSize, fileWindowDelta, 
			      justification, file_ix);
    unsigned check_n_ftrs = gomFS->numContinuous();
    unsigned check_n_labs = gomFS->numDiscrete();

    // Check that this input pfile is the same size as predecessors
    if(check_n_ftrs != n_ftrs || check_n_labs != n_labs) {
      sprintf(errmsg, "Features/labels of %s (%d/%d) don't match first input file (%d/%d)", 
	      file_name, check_n_labs, check_n_ftrs, n_labs, n_ftrs);
      error(errmsg);
    }
		
    // create the sentence-range iterator specifically for this 
    // input file
    Range sr_rng(sr[file_ix],0,gomFS->numSegments());
		
    for (Range::iterator srit=sr_rng.begin();!srit.at_end();srit++,global_seg_num++) {
      gomFS->openSegment(*srit);
      const size_t n_frames = gomFS->numFrames();
			
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
	const float* start_of_frame = gomFS->floatVecAtFrame(frame_no);
	const UInt32* start_of_unsigned_frame = gomFS->unsignedVecAtFrame(frame_no);
	for(unsigned feat_no = 0;  feat_no < n_ftrs; ++feat_no) {
	  ftr_buf[frame_no * n_ftrs + feat_no] = *(start_of_frame  + feat_no);
	}
	for(unsigned unsigned_feat_no = 0;  unsigned_feat_no < n_labs; ++unsigned_feat_no) {
	  lab_buf[frame_no*n_labs + unsigned_feat_no] = *(start_of_unsigned_frame+unsigned_feat_no);
	}
				
      }
      // Write output.
      printSegment(global_seg_num, out_fp, ftr_buf,n_ftrs,lab_buf,n_labs,n_frames, dontPrintFrameID,
		   quiet, ofmt, debug_level, oswap, out_stream);
    }
  }
    
  // All done; close output
  // must delete pfile object first, because it needs to rewrite
  if(ofmt==PFILE)
    delete out_stream;
  if (out_fp)
    if (fclose(out_fp))
      error("Couldn't close output file.");
	
  delete lab_buf;
  delete ftr_buf;
}

const char * ofmtStr="flatasc";
unsigned ofmt;

int  debug_level = 0;
bool dontPrintFrameID = false;
bool quiet = false;
#ifdef INTV_WORDS_BIGENDIAN
bool oswap = true;
#else
bool oswap             = false;
#endif 

Arg Arg::Args[] = {

#define GMTK_ARGUMENTS_DOCUMENTATION
#include "ObsArguments.h"
#undef GMTK_ARGUMENTS_DOCUMENTATION

  Arg("\n*** Output arguments ***\n"),

  Arg("o",      Arg::Opt, output_fname,"output file"),
  Arg("ofmt",      Arg::Opt, ofmtStr,"format of output file (htk, binary, ascii, pfile, flatbin, flatasc)"),
  Arg("olist",      Arg::Opt, outputList,"output list-of-files name.  Only meaningful if used with the RAW or HTK formats."),
  Arg("sep",      Arg::Opt, outputNameSeparatorStr,"String to use as separator when outputting raw ascii or binary files (one sentence per file)."),
  Arg("oswp",Arg::Opt, oswap,"do byte swapping on the output file"),
  Arg("ns",    Arg::Opt, dontPrintFrameID,"Don't print the frame IDs (i.e., sent and frame #)"),

  Arg("\n*** Misc arguments ***\n"),

  Arg("debug",  Arg::Opt, debug_level,"number giving level of debugging output to produce 0=none"),
  Arg("q",  Arg::Tog, quiet,"quiet mode"),
  // The argumentless argument marks the end of the above list.
  Arg()
};



int main(int argc, const char *argv[]) {

  //////////////////////////////////////////////////////////////////////
  // Check all necessary arguments provided before creating objects.
  //////////////////////////////////////////////////////////////////////
  
  int numFiles=0;
  
  CODE_TO_COMPUTE_ENDIAN
  oswap=doWeSwap;
  
  ///////////////////////////////////////////


  bool parse_was_ok = Arg::parse(argc,(char**)argv);

  if(!parse_was_ok) {
    Arg::usage(); exit(-1);
  }

#define GMTK_ARGUMENTS_CHECK_ARGS
#include "ObsArguments.h"
#undef GMTK_ARGUMENTS_CHECK_ARGS


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
 
  for(int i=0; i < MAX_NUM_OBS_FILES; i+=1) 
    if (ofs[i]) numFiles += 1;

  for(int i=0; i < numFiles; ++i) {
    // first check that we don't have gaps
    if (ofs[i] == NULL) {
      error("ERROR: Cannot skip an observation file number.");
    }

    if(output_fname!=NULL && strcmp(ofs[i],output_fname)==0) {
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
  } else {
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
 
  
  if(nis[0]==0 && nfs[0]==0) {
    error("The number of floats and the number of ints cannot be both zero.");
  }
  
  //////////////////////////////////////////////////////////////////////
  // Create objects.
  //////////////////////////////////////////////////////////////////////
    
  // done in subroutine now, after we open the first input file

  //////////////////////////////////////////////////////////////////////
  // Do the work.
  //////////////////////////////////////////////////////////////////////

  obsConcat(out_fp, 
	    numFiles,
	    ofmt,
	    dontPrintFrameID,
	    debug_level, 
	    quiet,
	    oswap);

  //////////////////////////////////////////////////////////////////////
  // Clean up and exit.
  //////////////////////////////////////////////////////////////////////

  return EXIT_SUCCESS;
}
