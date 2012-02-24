
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

#include "GMTK_ObservationSource.h"
#include "GMTK_FileSource.h"
#include "GMTK_ASCIIFile.h"
#include "GMTK_FlatASCIIFile.h"
#include "GMTK_PFileFile.h"
#include "GMTK_HTKFile.h"
#include "GMTK_HDF5File.h"
#include "GMTK_BinaryFile.h"
#include "GMTK_Stream.h"

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
#define GMTK_ARGUMENTS_DOCUMENTATION
#include "ObsArguments.h"
#undef GMTK_ARGUMENTS_DOCUMENTATION
  // final one to signal the end of the list
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

  FileSource globalObservationMatrix;
  ObservationFile *obsFile[MAX_NUM_OBS_FILES];
  unsigned nFiles=0;
  for (unsigned i=0; i < MAX_NUM_OBS_FILES && ofs[i] != NULL; i+=1, nFiles+=1) {
    switch (ifmts[i]) {
    case RAWASC:
      obsFile[i] = new ASCIIFile(ofs[i], nfs[i], nis[i], i,
				 Cpp_If_Ascii, cppCommandOptions,
				 frs[i], irs[i], prepr[i], sr[i]);
      break;
    case PFILE:
      obsFile[i] = new PFileFile(ofs[i], nfs[i], nis[i], i, iswp[i], frs[i], irs[i], prepr[i], sr[i]);
      break;
    case HTK:
      obsFile[i] = new HTKFile(ofs[i], nfs[i], nis[i], i, iswp[i], Cpp_If_Ascii, cppCommandOptions,
			       frs[i], irs[i], prepr[i], sr[i]);
      break;
    case HDF5:
      obsFile[i] = new HDF5File(ofs[i], i, Cpp_If_Ascii, cppCommandOptions,
				frs[i], irs[i], prepr[i], sr[i]);
      break;
    case FLATASC:
      obsFile[i] = new FlatASCIIFile(ofs[i], nfs[i], nis[i], i, Cpp_If_Ascii, cppCommandOptions,
				     frs[i], irs[i], prepr[i], sr[i]);
      break;
    case RAWBIN:
      obsFile[i] = new BinaryFile(ofs[i], nfs[i], nis[i], i, iswp[i], Cpp_If_Ascii, cppCommandOptions,
				  frs[i], irs[i], prepr[i], sr[i]);
      break;
    default:
      error("ERROR: Unknown observation file format type: '%s'\n", fmts[i]);
    }
  }
  globalObservationMatrix.initialize(nFiles, obsFile, gpr_str, startSkip, endSkip);
  FileSource *f = &globalObservationMatrix;
  for (unsigned j=0; j < f->numSegments(); j+=1) {
    assert(f->openSegment(j));
    for (unsigned k=0; k < f->numFrames(); k+=1) {
      printf("%03u %03u", j, k);
      Data32 const *buf = f->loadFrames(k,1);
      for (unsigned ff=0; ff < f->numContinuous(); ff+=1)
	printf(" %f", *((float *)(buf++)));
      for (unsigned ff=0; ff < f->numDiscrete(); ff+=1)
	printf(" %d", *((int *)(buf++)));
      printf("\n");
    }
  }
  
  exit(0);
}
