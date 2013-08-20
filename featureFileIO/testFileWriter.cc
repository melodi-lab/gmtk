
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <assert.h>


/*
 * Written by Richard Rogers <rprogers@uw.edu>
 *
 * Copyright (C) 2012 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */


#ifdef HAVE_CONFIG_H
#  include <config.h>
static const char * gmtk_version_id = PACKAGE_STRING;
#  ifdef HAVE_HG_H
#    include "hgstamp.h"
#  endif

#else 
// TODO: automate the process of updating this string.
static const char * gmtk_version_id = "GMTK Version 0.2b Tue Jan 20 22:59:41 2004";
#endif

#include "general.h"
#include "error.h"
#include "debug.h"
#include "arguments.h"

#include "GMTK_ObservationSource.h"
#include "GMTK_FileSource.h"
#include "GMTK_ASCIIFile.h"
#include "GMTK_BinaryFile.h"
#include "GMTK_FlatASCIIFile.h"
#include "GMTK_HDF5File.h"
#include "GMTK_HTKFile.h"
#include "GMTK_PFileFile.h"
#include "GMTK_CreateFileSource.h"
#include "GMTK_Stream.h"
unsigned chunkSize = 1;


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
const char * ofmtStr="flatascii";
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
  // final one to signal the end of the list
  Arg("chunkSize", Arg::Opt,chunkSize, "# of frames to load/print at once"),

  Arg("\n*** Output arguments ***\n"),

  Arg("o",      Arg::Req, output_fname,"output file"),
  Arg("ofmt",      Arg::Opt, ofmtStr,"format of output file (htk, binary, ascii, pfile, flatbin, flatasc)"),
  Arg("olist",      Arg::Opt, outputList,"output list-of-files name.  Only meaningful if used with the RAW or HTK formats."),
  Arg("sep",      Arg::Opt, outputNameSeparatorStr,"String to use as separator when outputting raw ascii or binary files (one sentence per file)."),
  Arg("oswp",Arg::Opt, oswap,"do byte swapping on the output file"),
  Arg("ns",    Arg::Opt, dontPrintFrameID,"Don't print the frame IDs (i.e., sent and frame #)"),

  Arg()
};


int 
main(int argc, char *argv[]) {

  bool parse_was_ok = Arg::parse(argc,(char**)argv);

  if(help) {
    Arg::usage();
    exit(0);
  }
  if(!parse_was_ok) {
    Arg::usage(); exit(-1);
  }
  if (print_version_and_exit) {
#ifdef HAVE_CONFIG_H
    printf("%s (Mercurial id: %s)\n",gmtk_version_id,HGID);
#else
    printf("%s\n", gmtk_version_id);
#endif
    exit(0);
  }
  infoMsg(IM::Max,"Finished parsing arguments\n");
#define GMTK_ARGUMENTS_CHECK_ARGS
#include "ObsArguments.h"
#undef GMTK_ARGUMENTS_CHECK_ARGS

  infoMsg(IM::Max,"Opening Files ...\n");

  FileSource *f = instantiateFileSource();

  ObservationFile *obsFile = NULL;
  ofmt = formatStrToNumber(ofmtStr);
  switch (ofmt) {
  case FLATASC:
    obsFile = new FlatASCIIFile(output_fname, f->numContinuous(), f->numDiscrete());
    break;
  case RAWASC:
    obsFile = new ASCIIFile(outputList, output_fname, outputNameSeparatorStr, f->numContinuous(), f->numDiscrete());
    break;
  case PFILE:
    obsFile = new PFileFile(output_fname, f->numContinuous(), f->numDiscrete(), oswap);
    break;
  case HTK:
    obsFile = new HTKFile(outputList, output_fname, outputNameSeparatorStr, oswap, f->numContinuous(), f->numDiscrete());
    break;
  case HDF5:
    error("ERROR: Unsupported observation file format type: '%s'\n", ofmtStr);
  case RAWBIN:
    obsFile = new BinaryFile(outputList, output_fname, outputNameSeparatorStr, oswap, f->numContinuous(), f->numDiscrete());
    break;
  default:
    error("ERROR: Unknown observation file format type: '%s'\n", ofmtStr);
  }

  for (unsigned j=0; j < f->numSegments(); j+=1) {
    assert(f->openSegment(j));
    printf("Processing sentence %u\n", j);
#if 0
    for (unsigned i=0; i < f->numFrames(); i+=1) {
#if 0
      obsFile->writeFrame(f->loadFrames(i, 1));
#else
      Data32 const *frame = f->loadFrames(i, 1);
      for (unsigned j=0; j < f->numContinuous(); j+=1) {
	obsFile->writeFeature((Data32)frame[j]);
      }
      for (unsigned j=f->numContinuous(); j < f->numFeatures(); j+=1) {
	obsFile->writeFeature((Data32)frame[j]);
      }
#endif
    }
    obsFile->endOfSegment();
#else
    obsFile->writeSegment(f->loadFrames(0, f->numFrames()), f->numFrames());
#endif
  }
  delete obsFile;
  exit(0);
}
