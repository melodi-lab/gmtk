/*  Generated header
 *  File Name : GMTK_ObsInfo.cc
 *
 *  Created   : 2003-12-03 11:59:48 karim
 *  Author    : Karim Filali (karim@cs.washington.edu)
 *  Time-stamp: <>
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

#if 0
#  include "GMTK_ObservationMatrix.h"
#else
//#  include "GMTK_ObservationSource.h"
#  include "GMTK_FileSource.h"
#  include "GMTK_CreateFileSource.h"
#  include "GMTK_ASCIIFile.h"
#  include "GMTK_FlatASCIIFile.h"
#  include "GMTK_PFileFile.h"
#  include "GMTK_HTKFile.h"
#  include "GMTK_HDF5File.h"
#  include "GMTK_BinaryFile.h"
#  include "GMTK_Filter.h"
#  include "GMTK_Stream.h"
#endif


#if 0
ObservationMatrix globalObservationMatrix;
#else
FileSource *gomFS;
#endif

void obsInfo(FILE* out_fp, FileSource* obs_mat, bool dont_print_info, bool print_sent_frames, bool print_stream_info) {

  unsigned num_segments      = obs_mat->numSegments();
  unsigned num_streams       = obs_mat->numFiles();
  unsigned total_num_frames  = 0;
  //  StreamInfo* current_stream = NULL;

  for (unsigned seg_no=0; seg_no < num_segments; ++seg_no) {
    obs_mat->openSegment(seg_no);
    total_num_frames += obs_mat->numFrames();
  }
  
  if (!dont_print_info) {
    fprintf(out_fp,"%d stream(s), %d sentences, %d frames, %d discrete feature(s), %d continuous feature(s)\n",
	    num_streams,
	    num_segments,
            total_num_frames,
	    obs_mat->numDiscrete(),
	    obs_mat->numContinuous());
  }
  
  if(print_stream_info) {
#if 0
    for (unsigned stream_no=0; stream_no < num_streams; ++stream_no) {
      current_stream = obs_mat->getStream(stream_no);
      assert(current_stream != NULL);
      fprintf(out_fp,"stream %d: %d discrete feature(s), %d continuous feature(s)\n",stream_no,current_stream->getNumInts(),current_stream->getNumFloats());
    }
#else
    error("printing stream info is not supported with the new observation implementation\n");
#endif
  }

  if (print_sent_frames) {
      for (unsigned seg_no=0; seg_no < num_segments; ++seg_no) {
	obs_mat->openSegment(seg_no);
	fprintf(out_fp,"%d %d\n",seg_no,obs_mat->numFrames());
      }
  }

}


// #ifdef MAIN

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <assert.h>

//#include "pfile.h"
#include "error.h"
#include "arguments.h"
#include "GMTK_WordOrganization.h"



//#define MAX_OBJECTS 60

char *   output_fname      = NULL;
bool     Print_Stream_Info = false;
bool     Print_Sent_Frames = false;
bool     quiet = false;

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

Arg Arg::Args[] = {
  Arg("\n*** ObsInfo options ***\n"),
  Arg("o",        Arg::Opt, output_fname,"output file"),
  Arg("s",        Arg::Opt, Print_Stream_Info,"Also print individual stream info."),
  Arg("p",        Arg::Opt, Print_Sent_Frames,"Also print # frames for each sentence."),
  Arg("q",        Arg::Tog, quiet,"Don't print the normal info (i.e., useful with -p option)."),
#define GMTK_ARGUMENTS_DOCUMENTATION
#include "ObsArguments.h"
#undef GMTK_ARGUMENTS_DOCUMENTATION

  // The argumentless argument marks the end of the above list.
  Arg()
};


int main(int argc, const char *argv[]) {

  CODE_TO_COMPUTE_ENDIAN

  //////////////////////////////////////////////////////////////////////
  // Check all necessary arguments provided before creating objects.
  //////////////////////////////////////////////////////////////////////
  
   
  int numFiles=0;

  bool parse_was_ok = Arg::parse(argc,(char**)argv);

  if(!parse_was_ok) {
    Arg::usage(); exit(-1);
  }

#define GMTK_ARGUMENTS_CHECK_ARGS
#include "ObsArguments.h"
#undef GMTK_ARGUMENTS_CHECK_ARGS


  //////////////////////////////////////////////////////////////////////
  // Check all necessary arguments provided before creating objects.
  //////////////////////////////////////////////////////////////////////


 for(int i=0; i < numFiles; ++i)
   if(output_fname!=NULL && strcmp(ofs[i],output_fname)==0) {
     error("Input and output filenames cannot be the same.");
   }


  FILE *out_fp=NULL;
  if (output_fname==0 || !strcmp(output_fname,"-"))
    out_fp = stdout;
  else {
    if ((out_fp = fopen(output_fname, "w")) == NULL) {
      error("Couldn't open output file (%s) for writing.", output_fname);
    }
  }


 //////////////////////////////////////////////////////////////////////
 // Create objects.
 //////////////////////////////////////////////////////////////////////
 // If we have a pfile, we can extract the number if features from the file directly

#if 0
 globalObservationMatrix.openFiles(numFiles,  // number of files.   For now we use only one
				   (const char**)&input_fname,
				   NULL,      // all feature range
				   NULL,      // all lab range
				   (unsigned*)&nfs,
				   (unsigned*)&nis,
				   (unsigned*)&ifmt,
				   (bool*)&iswap,
				   0,  // no startSkip
				   0,  // no endSkip
				   cppIfAscii,
				   cppCommandOptions,
				   NULL,     // no per stream range
				   actionIfDiffNumFrames,
				   actionIfDiffNumSents
				   );   
#else
 gomFS = instantiateFileSource();
#endif



    //////////////////////////////////////////////////////////////////////
    // Do the work.
    //////////////////////////////////////////////////////////////////////
    
    obsInfo(out_fp, gomFS,quiet,Print_Sent_Frames,Print_Stream_Info);

    if (fclose(out_fp))
      error("Couldn't close output file.");

    return 0;
}
// #endif
