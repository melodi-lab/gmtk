
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <assert.h>

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
#include "vbyteswapping.h"

#include "GMTK_WordOrganization.h"
#include "GMTK_ObservationSource.h"
#include "GMTK_ASCIIStream.h"
#include "GMTK_BinStream.h"
#include "GMTK_FileSource.h"
#include "GMTK_ASCIIFile.h"
#include "GMTK_FlatASCIIFile.h"
#include "GMTK_PFileFile.h"
#include "GMTK_HTKFile.h"
#include "GMTK_HDF5File.h"
#include "GMTK_BinaryFile.h"
#include "GMTK_Filter.h"
#include "GMTK_FilterFile.h"
#include "GMTK_FIRFilter.h"
#include "GMTK_AffineFilter.h"
#include "GMTK_Stream.h"


#define OPTIONAL_OBSERVATION_FILES
#define GMTK_ARG_OBS_FILES
#define GMTK_ARG_CPP_CMD_OPTS
#define GMTK_ARG_OBS_MATRIX_XFORMATION
#define GMTK_ARG_FILE_RANGE_OPTIONS
#define GMTK_ARG_START_END_SKIP
#define GMTK_ARG_HELP
#define GMTK_ARG_VERSION
#define GMTK_ARG_STREAMING_INPUT
#define GMTK_ARG_STREAMING_OUTPUT
#define GMTK_ARGUMENTS_DEFINITION
#include "ObsArguments.h"
#undef GMTK_ARGUMENTS_DEFINITION

Arg Arg::Args[] = {
#define GMTK_ARGUMENTS_DOCUMENTATION
#include "ObsArguments.h"
#undef GMTK_ARGUMENTS_DOCUMENTATION
  // final one to signal the end of the list
  Arg()
};


void
sendFrame(Data32 const *frame, unsigned segNum, unsigned frameNum, 
	  unsigned numCont, unsigned numDisc, bool binary, bool swap)
{
  if (binary) {
    printf("F");
    Data32 const *buf = frame;
    for (unsigned ff=0; ff < numCont; ff+=1) {
      float obs = *((float *)(buf++));
      if (swap) obs = swapb_f32_f32(obs);
      if (fwrite(&obs, sizeof(obs), 1, stdout) != 1) {
	error("ERROR: failed to write the %uth float in segment %u frame %u",ff, segNum, frameNum);
      }
    }
    for (unsigned ff=0; ff < numDisc; ff+=1) {
      unsigned obs = *((unsigned *)(buf++));
      if (swap) obs = swapb_i32_i32(obs);
      if (fwrite(&obs, sizeof(obs), 1, stdout) != 1) {
	error("ERROR: failed to write the %uth int in segment %u frame %u",ff, segNum, frameNum);
      }
    }
  } else {

    // FIXME - matlab format support
    if (prettyPrintStream) {
      printf("%u %u: ", segNum, frameNum);
    }
    printf("F");
    Data32 const *buf = frame;
    for (unsigned ff=0; ff < numCont; ff+=1)
      printf(" %f", *((float *)(buf++)));
    for (unsigned ff=0; ff < numDisc; ff+=1)
      printf(" %d", *((int *)(buf++)));
    printf("\n");
  }
}



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

  
  bool machineBigEndian = getWordOrganization() == BYTE_BIG_ENDIAN;
  bool needOutputSwap   = machineBigEndian != outputNetByteOrder;

  if (binaryOutputStream) {
    printf("%s%s", GMTK_BIN_PROTOCOL_COOKIE, GMTK_BIN_PROTOCOL_VERSION);
  } else {
    printf("%s%s", GMTK_ASC_PROTOCOL_COOKIE, GMTK_ASC_PROTOCOL_VERSION);
  }
  if (streamSource) {

    // get input from another stream

    FILE *inFile;

    if (strcmp("-", os)) {
      inFile = fopen(os, binaryInputStream ? "rb" : "r");
    } else {
      inFile = stdin;
    }

    if (!inFile) {
      error("ERROR: '%s' %s", os, strerror(errno));
    }
    
    // FIXME - use StreamSource & support posttrans ?

    // FIXME - add -sfr and -sir for feature range selection ?
    ObservationStream *inStream;
    if (binaryInputStream) {
      inStream = new BinaryStream(inFile, snf, sni, NULL /* sfr */ , NULL /* sir */, inputNetByteOrder);
    } else {
      inStream = new  ASCIIStream(inFile, snf, sni, NULL /* sfr */ , NULL /* sir */);
    }
    assert(inStream);

    unsigned segNum = 0;
    unsigned frmNum = 0;

    Data32 const *frame;
    for (; !inStream->EOS(); ) {
      frame = inStream->getNextLogicalFrame();
      if (!frame) {
	segNum += 1;
	frmNum = 0;
	printf("%s", binaryOutputStream ? "E" : "E\n");
	continue;
      }
      sendFrame(frame, segNum, frmNum, inStream->numContinuous(), inStream->numDiscrete(),
		binaryOutputStream, needOutputSwap);
      frmNum += 1;
    }

  } else {

    // stream data out of random-access file(s)
    
    infoMsg(IM::Max,"Opening Files ...\n");

    FileSource globalObservationMatrix;
    ObservationFile *obsFile[MAX_NUM_OBS_FILES];
    unsigned nFiles=0;
    unsigned nCont = 0;
    for (unsigned i=0; i < MAX_NUM_OBS_FILES && ofs[i] != NULL; i+=1, nFiles+=1) {

      obsFile[i] = instantiateFile(ifmts[i], ofs[i], nfs[i], nis[i], i, iswp[i],
				   Cpp_If_Ascii, cppCommandOptions, prefrs[i], preirs[i],
				   prepr[i], sr[i]);
      assert(obsFile[i]);
      Filter *fileFilter = instantiateFilters(Per_Stream_Transforms[i],
					      obsFile[i]->numContinuous());
      if (fileFilter) {
	obsFile[i] = new FilterFile(fileFilter, obsFile[i], frs[i], irs[i], postpr[i]);
	nCont += obsFile[i]->numContinuous();
      } else
	error("current implementation requires filter");
    }
    globalObservationMatrix.initialize(nFiles, obsFile, 
				       1024*1024 /* FIXME - argument */, 
				       Action_If_Diff_Num_Sents,
				       Action_If_Diff_Num_Frames,
				       gpr_str, startSkip, endSkip,
				       instantiateFilters(Post_Transforms, nCont));

    FileSource *f = &globalObservationMatrix;
    for (unsigned j=0; j < f->numSegments(); j+=1) {
      assert(f->openSegment(j));
      for (unsigned k=0; k < f->numFrames(); k+=1) {
	Data32 const *buf = f->loadFrames(k,1);
	sendFrame(buf, j, k, f->numContinuous(), f->numDiscrete(),
		  binaryOutputStream, needOutputSwap);
	
      }
      printf("%s", binaryOutputStream ? "E" : "E\n");
    }
  }
  exit(0);
}




#if 0

    FileSource *f = &globalObservationMatrix;
    for (unsigned j=0; j < f->numSegments(); j+=1) {
      assert(f->openSegment(j));
      if (binaryOutputStream) {
	for (unsigned k=0; k < f->numFrames(); k+=1) {
	  // binary output
	  printf("F");
	  Data32 const *buf = f->loadFrames(k,1);
	  for (unsigned ff=0; ff < f->numContinuous(); ff+=1) {
	    float obs = *((float *)(buf++));
	    if (needOutputSwap) obs = swapb_f32_f32(obs);
	    if (fwrite(&obs, sizeof(obs), 1, stdout) != 1) {
	      error("ERROR: failed to write the %uth float in segment %u frame %u",ff, j, k);
	    }
	  }
	  for (unsigned ff=0; ff < f->numDiscrete(); ff+=1) {
	    unsigned obs = *((unsigned *)(buf++));
	    if (needOutputSwap) obs = swapb_i32_i32(obs);
	    if (fwrite(&obs, sizeof(obs), 1, stdout) != 1) {
	      error("ERROR: failed to write the %uth int in segment %u frame %u",ff, j, k);
	    }
	  }
	}
	printf("E");
      } else {
	for (unsigned k=0; k < f->numFrames(); k+=1) {
	  // ASCII output
	  if (prettyPrintStream) {
	    printf("%u %u: ", j, k);
	  }
	  printf("F");
	  Data32 const *buf = f->loadFrames(k,1);
	  for (unsigned ff=0; ff < f->numContinuous(); ff+=1)
	    printf(" %f", *((float *)(buf++)));
	  for (unsigned ff=0; ff < f->numDiscrete(); ff+=1)
	    printf(" %d", *((int *)(buf++)));
	  printf("\n");
	}
	printf("E\n");
      }
    }
#endif
