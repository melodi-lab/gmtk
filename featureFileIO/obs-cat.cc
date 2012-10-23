
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
#include "GMTK_FileStream.h"
#include "GMTK_FileSrcStream.h"
#include "GMTK_StreamSource.h"
#include "GMTK_FileSource.h"
#include "GMTK_CreateFileSource.h"
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


#define GMTK_ARG_CPP_CMD_OPTS
#define GMTK_ARG_OBS_MATRIX_XFORMATION
#define GMTK_ARG_FILE_RANGE_OPTIONS
#define GMTK_ARG_START_END_SKIP
#define GMTK_ARG_HELP
#define GMTK_ARG_VERSION
#define GMTK_ARG_STREAM_AND_FILE_INPUT
#define GMTK_ARG_STREAMING_OUTPUT
#define GMTK_ARGUMENTS_DEFINITION
#include "ObsArguments.h"
#undef GMTK_ARGUMENTS_DEFINITION

unsigned verb = IM::Default;

Arg Arg::Args[] = {
#define GMTK_ARGUMENTS_DOCUMENTATION
#include "ObsArguments.h"
#undef GMTK_ARGUMENTS_DOCUMENTATION
  Arg("verbosity",Arg::Opt, verb, "Level of debugging output, [0-100]"),
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


void
sendFeatureCounts(ObservationSource *src, bool binaryOutputStream, bool needOutputSwap) {
  unsigned nCont = src->numContinuous();
  unsigned nDisc = src->numDiscrete();
  if (binaryOutputStream) {
    if (needOutputSwap) {
      nCont = swapb_i32_i32(nCont);
      nDisc = swapb_i32_i32(nDisc);
    }
    if (fwrite(&nCont, sizeof(nCont), 1, stdout) != 1) {
      error("ERROR: failed to write the number of continous features");
    }
    if (fwrite(&nDisc, sizeof(nDisc), 1, stdout) != 1) {
      error("ERROR: failed to write the number of discrete features");
    }
  } else {
    printf("%u %u\n", src->numContinuous(), src->numDiscrete());
  }
}


ObservationStream * 
makeFileStream(unsigned ifmts[], unsigned i) {
  ObservationFile *obsFile = 
    instantiateFile(ifmts[i], ofs[i], nfs[i], nis[i], i, iswp[i],
		    Cpp_If_Ascii, cppCommandOptions, prefrs[i], preirs[i],
		    prepr[i], sr[i]);
  assert(obsFile);
  if (Per_Stream_Transforms[i]) {
    Filter *fileFilter = instantiateFilters(Per_Stream_Transforms[i],
					    obsFile->numLogicalContinuous(),
					    obsFile->numLogicalDiscrete());
    if (fileFilter) {
      obsFile = new FilterFile(fileFilter, obsFile, frs[i], irs[i], postpr[i]);
    } else {
      error("ERROR: failed to create filter for '%s'",Per_Stream_Transforms[i]);
    }
  }
  return new FileStream(obsFile);
}


ObservationStream *
makeStream(unsigned ifmts[], unsigned i) {
  FILE *inFile;
  
  if (strcmp("-", oss[i])) {
    inFile = fopen(oss[i], ifmts[i] == RAWBIN ? "rb" : "r");
  } else {
    inFile = stdin;
  }
  
  if (!inFile) {
    error("ERROR: '%s' %s", oss[i], strerror(errno));
  }
  
  if (ifmts[i] == RAWBIN) {
    return new BinaryStream(inFile, nfs[i], nis[i], frs[i], irs[i], inputNetByteOrder[i]);
  } else if (ifmts[i] == RAWASC) {
    return new  ASCIIStream(inFile, nfs[i], nis[i], frs[i], irs[i]);
  } else {
    error("ERROR: -fmt%u must be 'binary' or 'ascii', got '%s'", i, fmts[i]);
  }
  return NULL; // should never reach here
}


ObservationStream *
makeFileSource() {
#if 0
  ObservationFile   *obsFile[MAX_NUM_OBS_FILES] = {NULL,NULL,NULL,NULL,NULL};
  unsigned nFiles = 0;
  unsigned nCont  = 0;
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
      error("current implementation requires filter\n");
  }
  if (nFiles == 0) {
    error("ERROR: no input files specified");
  }
  FileSource *fileSrc = new FileSource();
  fileSrc->initialize(nFiles, obsFile, fileBufferSize,
		      Action_If_Diff_Num_Sents,
		      Action_If_Diff_Num_Frames,
		      gpr_str, startSkip, endSkip,
		      instantiateFilters(Post_Transforms, nCont));
  return new FileSrcStream(fileSrc);
#else
  return new FileSrcStream(instantiateFileSource());
#endif
}


int 
main(int argc, char *argv[]) {

  CODE_TO_COMPUTE_ENDIAN

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

  IM::setGlbMsgLevel(IM::ObsStream, verb);
  
  bool machineBigEndian = getWordOrganization() == BYTE_BIG_ENDIAN;
  bool needOutputSwap   = machineBigEndian != outputNetByteOrder;

  if (binaryOutputStream) {
    printf("%s%s", GMTK_BIN_PROTOCOL_COOKIE, GMTK_BIN_PROTOCOL_VERSION);
  } else {
    printf("%s%s", GMTK_ASC_PROTOCOL_COOKIE, GMTK_ASC_PROTOCOL_VERSION);
  }
  
  ObservationStream *obsStream[MAX_NUM_OBS_STREAMS] = {NULL,NULL,NULL,NULL,NULL};
  unsigned           nStreams = 0;

  // If there are only files, we'll use a FileSrcStream wrapper around
  // a FileSource so that we can do -fdiffact etc. Otherwise, we wrap
  // FileStreams around the individual ObservationFiles.
 
  bool allFiles = true;
  for (unsigned i=0; i < MAX_NUM_OBS_STREAMS; i+=1) {
    allFiles = allFiles && (oss[i] == NULL);
  }
  if (allFiles) {
    obsStream[0] = makeFileSource();
    if (!obsStream[0]) {
      error("ERROR: no input sources specified\n");
    }
    nStreams = 1;
  } else {
    for (unsigned i=0; i < MAX_NUM_OBS_STREAMS && (ofs[i] || oss[i]); i+=1) {
      if (ofs[i]) {
	obsStream[i] = makeFileStream(ifmts, i);
      } else if (oss[i]) {
	obsStream[i] = makeStream(ifmts, i);
      } else {
	error("ERROR: no input stream or file specified for -os%u or -of%u", i,i);
      }
      assert(obsStream[i]);
      nStreams += 1;
    }
  }
  if (nStreams == 0) {
    error("ERROR: no input sources specified (use -ofX or -osX)");
  }

  StreamSource *source = new StreamSource(nStreams, obsStream, streamBufferSize, Post_Transforms, startSkip);
  sendFeatureCounts(source, binaryOutputStream, needOutputSwap);

  unsigned segNum = 0;
  unsigned frmNum = 0;
  
  Data32 const *frame;
  for (; !source->EOS(); segNum += 1) {
    source->preloadFrames(1);
    for (frmNum=0; source->numFrames() == 0 || frmNum < source->numFrames(); frmNum += 1) {
      frame = source->loadFrames(frmNum, 1);
      if (!frame) {
	error("ERROR: incomplete segment");
      }
      if (source->numFrames() != 0 && frmNum >= source->numFrames()) {
	assert(frame == NULL);
	continue;
      }
      sendFrame(frame, segNum, frmNum, source->numContinuous(), source->numDiscrete(),
		binaryOutputStream, needOutputSwap);
      source->enqueueFrames(1);
    }
    printf("%s", binaryOutputStream ? "E" : "E\n");
  }

  exit(0);
}

