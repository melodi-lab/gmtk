/*  Generated header
 *  File Name : GMTK_ObsInfo.cc
 *
 *  Created   : 2003-12-03 11:59:48 karim
 *  Author    : Karim Filali (karim@cs.washington.edu)
 *  Time-stamp: <>
*/

#include "GMTK_ObservationMatrix.h"


ObservationMatrix globalObservationMatrix;

void obsInfo(FILE* out_fp, ObservationMatrix* obs_mat, bool dont_print_info, bool print_sent_frames, bool print_stream_info) {

  unsigned num_segments      = obs_mat->numSegments();
  unsigned num_streams       = obs_mat->numStreams();
  unsigned total_num_frames  = 0;
  StreamInfo* current_stream = NULL;

  for (unsigned seg_no=0; seg_no < num_segments; ++seg_no) {
    total_num_frames += obs_mat->numFrames(seg_no);
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
    for (unsigned stream_no=0; stream_no < num_streams; ++stream_no) {
      current_stream = obs_mat->getStream(stream_no);
      assert(current_stream != NULL);
      fprintf(out_fp,"stream %d: %d discrete feature(s), %d continuous feature(s)\n",stream_no,current_stream->getNumInts(),current_stream->getNumFloats());
    }
  }

  if (print_sent_frames) {
      for (unsigned seg_no=0; seg_no < num_segments; ++seg_no) {
	fprintf(out_fp,"%d %d\n",seg_no,obs_mat->numFrames(seg_no));
      }
  }

}


#ifdef MAIN

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <values.h>
#include <math.h>
#include <assert.h>
#include "pfile.h"
#include "error.h"
#include "arguments.h"
#include "GMTK_WordOrganization.h"



#define MAX_OBJECTS 5

char *   input_fname[MAX_OBJECTS] = {NULL,NULL,NULL,NULL,NULL};  // Input file name.
char *   ifmtStr[MAX_OBJECTS]     = {"pfile","pfile","pfile","pfile","pfile"};
unsigned ifmt[MAX_OBJECTS];

char *   output_fname      = NULL;
bool     Print_Stream_Info = false;
bool     Print_Sent_Frames = false;

unsigned nis[MAX_OBJECTS];
unsigned nfs[MAX_OBJECTS];

char*    actionIfDiffNumFramesStr[MAX_OBJECTS]={"er","er","er","er","er"};   // 
unsigned actionIfDiffNumFrames[MAX_OBJECTS]={FRAMEMATCH_ERROR,FRAMEMATCH_ERROR,FRAMEMATCH_ERROR,FRAMEMATCH_ERROR,FRAMEMATCH_ERROR};   // 
char*    actionIfDiffNumSentsStr[MAX_OBJECTS]={"te","te","te","te","te"}; 
unsigned actionIfDiffNumSents[MAX_OBJECTS]={SEGMATCH_TRUNCATE_FROM_END,SEGMATCH_TRUNCATE_FROM_END,SEGMATCH_TRUNCATE_FROM_END,SEGMATCH_TRUNCATE_FROM_END,SEGMATCH_TRUNCATE_FROM_END};   // 

bool     quiet = false;
#ifdef INTV_WORDS_BIGENDIAN
bool     iswap[MAX_OBJECTS]={true,true,true,true,true};
#else
bool     iswap[MAX_OBJECTS]={false,false,false,false,false};
#endif 

bool     cppIfAscii        = true;
char*    cppCommandOptions = NULL;

bool     help              = false;

Arg Arg::Args[] = {
  Arg("i",        Arg::Req, input_fname,"input file",Arg::ARRAY,MAX_OBJECTS),
  Arg("ifmt",     Arg::Opt, ifmtStr ,"format of input file",Arg::ARRAY,MAX_OBJECTS),
  Arg("o",        Arg::Opt, output_fname,"output file"),
  Arg("s",        Arg::Opt, Print_Stream_Info,"Also print individual stream info."),
  Arg("p",        Arg::Opt, Print_Sent_Frames,"Also print # frames for each sentence."),
  Arg("nf",       Arg::Opt, nfs,"number of floats in input file",Arg::ARRAY,MAX_OBJECTS),
  Arg("ni",       Arg::Opt, nis,"number of ints (labels) in input file",Arg::ARRAY,MAX_OBJECTS),
  Arg("q",        Arg::Tog, quiet,"Don't print the normal info (i.e., useful with -p option)."),
  Arg("fdiffact", Arg::Opt, actionIfDiffNumFramesStr ,"Action if different number of frames in streams: error (er), repeat last frame (rl), first frame (rf), segmentally expand (se), truncate from start (ts), truncate from end (te)",Arg::ARRAY,MAX_OBJECTS),
  Arg("sdiffact", Arg::Opt, actionIfDiffNumSentsStr ,"Action if different number of sentences in streams: error (er), truncate from end (te), repeat last sent (rl), and wrap around (wa).",Arg::ARRAY,MAX_OBJECTS),
  Arg("iswap",    Arg::Opt, iswap,"do byte swapping on the input file",Arg::ARRAY,MAX_OBJECTS),
  Arg("cppifascii",        Arg::Opt, cppIfAscii,"Pre-process ASCII files using CPP"),
  Arg("cppCommandOptions", Arg::Opt, cppCommandOptions,"Additional CPP command line"),
  Arg("help",     Arg::Tog, help,"print this message"),
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

 for(int i=0; i < numFiles; ++i)
   if(output_fname!=NULL && strcmp(input_fname[i],output_fname)==0) {
     error("Input and output filenames cannot be the same.");
   }

 for(int i=0; i < MAX_OBJECTS; ++i) {
   if(input_fname[i]!=NULL) {
     if (strcmp(actionIfDiffNumFramesStr[i],"er") == 0)
       actionIfDiffNumFrames[i] = FRAMEMATCH_ERROR;
     else if (strcmp(actionIfDiffNumFramesStr[i],"rl") == 0)
       actionIfDiffNumFrames[i] = FRAMEMATCH_REPEAT_LAST;
     else if (strcmp(actionIfDiffNumFramesStr[i],"rf") == 0)
       actionIfDiffNumFrames[i] = FRAMEMATCH_REPEAT_FIRST;
     else if (strcmp(actionIfDiffNumFramesStr[i],"se") == 0)
       actionIfDiffNumFrames[i] = FRAMEMATCH_EXPAND_SEGMENTALLY;
     else if (strcmp(actionIfDiffNumFramesStr[i],"ts") == 0)
       actionIfDiffNumFrames[i] = FRAMEMATCH_TRUNCATE_FROM_START;
     else if (strcmp(actionIfDiffNumFramesStr[i],"te") == 0)
       actionIfDiffNumFrames[i] = FRAMEMATCH_TRUNCATE_FROM_END;
     else
       error("ERROR: Unknown action when diff num of frames: '%s'\n",actionIfDiffNumFramesStr[i]);
   }
 }

for(int i=0; i < MAX_OBJECTS; ++i) {
   if(input_fname[i]!=NULL) {
     if (strcmp(actionIfDiffNumSentsStr[i],"er") == 0)
       actionIfDiffNumSents[i] = SEGMATCH_ERROR;
     else if (strcmp(actionIfDiffNumSentsStr[i],"rl") == 0)
       actionIfDiffNumSents[i] = SEGMATCH_REPEAT_LAST;
     else if (strcmp(actionIfDiffNumSentsStr[i],"wa") == 0)
       actionIfDiffNumSents[i] = SEGMATCH_WRAP_AROUND;
     else if (strcmp(actionIfDiffNumSentsStr[i],"te") == 0)
       actionIfDiffNumSents[i] = SEGMATCH_TRUNCATE_FROM_END;
     else
       error("ERROR: Unknown action when diff num of sentences: '%s'\n",actionIfDiffNumSentsStr[i]);
   }
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
 for(int i=0; i < MAX_OBJECTS; ++i) {
   if(input_fname[i]!=NULL) {
     if(ifmt[i]==PFILE) {
       FILE *in_fp = fopen(input_fname[i], "r");
       if (in_fp==NULL) error("Couldn't open input pfile for reading.");
       InFtrLabStream_PFile* in_streamp = new InFtrLabStream_PFile(0,"",in_fp,1,iswap[i]);
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




    //////////////////////////////////////////////////////////////////////
    // Do the work.
    //////////////////////////////////////////////////////////////////////
    
    obsInfo(out_fp, &globalObservationMatrix,quiet,Print_Sent_Frames,Print_Stream_Info);

    if (fclose(out_fp))
      error("Couldn't close output file.");

    return 0;
}
#endif
