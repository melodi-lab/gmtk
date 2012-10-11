
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
#include "GMTK_Filter.h"
#include "GMTK_FilterFile.h"
#include "GMTK_MergeFile.h"
#include "GMTK_FIRFilter.h"
#include "GMTK_AffineFilter.h"
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
//  FilterFile *filteredFile[MAX_NUM_OBS_FILES];

  unsigned nFiles=0;
  unsigned nCont = 0;
  for (unsigned i=0; i < MAX_NUM_OBS_FILES && ofs[i] != NULL; i+=1, nFiles+=1) {

    obsFile[i] = instantiateFile(ifmts[i], ofs[i], nfs[i], nis[i], i, iswp[i],
				 Cpp_If_Ascii, cppCommandOptions, prefrs[i], preirs[i],
				 prepr[i], sr[i]);
    assert(obsFile[i]);
    Filter *fileFilter = instantiateFilters(Per_Stream_Transforms[i],
					    obsFile[i]->numLogicalContinuous(),
                                            obsFile[i]->numLogicalDiscrete());
    if (fileFilter) {
      obsFile[i] = new FilterFile(fileFilter, obsFile[i], frs[i], irs[i], postpr[i]);
      nCont += obsFile[i]->numContinuous();
    } else
      error("current implementation requires filter\n");
  }

  MergeFile *mf = new MergeFile(nFiles, obsFile, 
				Action_If_Diff_Num_Sents,
				Action_If_Diff_Num_Frames,
				Ftr_Combo);
#if 0
  globalObservationMatrix.initialize(nFiles, obsFile, 1024*1024 /* FIXME */,
				     Action_If_Diff_Num_Sents,
				     Action_If_Diff_Num_Frames,
				     gpr_str, startSkip, endSkip,
				     instantiateFilters(Post_Transforms, nCont),
				     justification, Ftr_Combo);
  FileSource *f = &globalObservationMatrix;
#endif

  for (unsigned j=0; j < mf->numSegments(); j+=1) {
    assert(mf->openSegment(j));
printf("Processing sentence %u\n", j);
    for (unsigned k=0; k < mf->numFrames(); k+=1) {
      //      printf("%03u %03u", j, k);
      printf("%u %u", j, k);
      Data32 const *buf = mf->getFrames(k,1);
      for (unsigned ff=0; ff < mf->numContinuous(); ff+=1)
	printf(" %f", *((float *)(buf++)));
      for (unsigned ff=0; ff < mf->numDiscrete(); ff+=1)
	printf(" %d", *((int *)(buf++)));
      printf("\n");
    }
  }
  delete mf;
  exit(0);
}
