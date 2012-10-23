
/*
 * GMTK_CreateFileSource.cc
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2012, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <string.h>

#include "error.h"
#include "general.h"
#include "machine-dependent.h"
#include "debug.h"

#include "GMTK_ObservationFile.h"
#include "GMTK_FilterFile.h"
#include "GMTK_FileSource.h"
#include "GMTK_FileSourceNoCache.h"
#include "GMTK_Filter.h"
#include "GMTK_MergeFile.h"
#include "GMTK_ObservationArguments.h"


extern char    *ofs[];
extern unsigned nfs[];
extern unsigned nis[];
extern const char   *fmts[];
extern const char    *frs[];
extern const char    *irs[];
extern const char    *prefrs[];
extern const char    *preirs[];
extern const char    *sr[];
extern char  *prepr[];
extern char *postpr[];
extern char *gpr_str;
extern unsigned justification;
extern bool iswp[];
extern unsigned ifmts[];
extern bool Cpp_If_Ascii;
extern char *cppCommandOptions;
extern char    *Per_Stream_Transforms[];
extern char    *Post_Transforms;
extern int startSkip;
extern int endSkip;
extern unsigned    Action_If_Diff_Num_Frames[];
extern unsigned    Action_If_Diff_Num_Sents[];
extern unsigned Ftr_Combo;
extern unsigned fileBufferSize;
extern unsigned fileWindowSize;
extern unsigned fileWindowDelta;
extern bool constantSpace;

#ifndef MAX_NUM_OBS_FILES
#define MAX_NUM_OBS_FILES 10
#endif

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
instantiateFileSource() {

  if (ofs == NULL)
    error("CreateFileSource: list of file names is NULL\n");
  if (fmts == NULL)
    error("CreateFileSource: list of file formats is NULL\n");
  if (nfs == NULL)
    error("CreateFileSource: list of number of floats is NULL\n");
  if (nis == NULL)
    error("CreateFileSource: list of number of ints is NULL\n");

  // range selection is much more efficient if "all" is replaced with NULL
  // since the logical <-> physical mapping step can be skipped
  for (unsigned i=0; i < MAX_NUM_OBS_FILES; i+=1) {
    if (   frs[i] && ALLOREMPTY(frs[i]))       frs[i] = NULL;
    if (   irs[i] && ALLOREMPTY(irs[i]))       irs[i] = NULL;
    if (prefrs[i] && ALLOREMPTY(prefrs[i])) prefrs[i] = NULL;
    if (preirs[i] && ALLOREMPTY(preirs[i])) preirs[i] = NULL;
    if (    sr[i] && ALLOREMPTY(sr[i]))         sr[i] = NULL;
    if ( prepr[i] && ALLOREMPTY(prepr[i]))   prepr[i] = NULL;
    if (postpr[i] && ALLOREMPTY(postpr[i])) postpr[i] = NULL;
  }
  if (gpr_str && ALLOREMPTY(gpr_str)) gpr_str = NULL;

  ObservationFile *obsFile[MAX_NUM_OBS_FILES];
  unsigned nFiles=0;
  for (unsigned i=0; i < MAX_NUM_OBS_FILES && ofs[i] != NULL; i+=1, nFiles+=1) {
    obsFile[i] = instantiateFile(ifmts[i], ofs[i], nfs[i], nis[i], i, iswp[i],
                                 Cpp_If_Ascii, cppCommandOptions, prefrs[i], preirs[i],
                                 prepr[i], sr[i]);
    assert(obsFile[i]);
    if (Per_Stream_Transforms[i] || frs[i] || irs[i] || postpr[i]) {
      Filter *fileFilter = instantiateFilters(Per_Stream_Transforms[i],
					      obsFile[i]->numContinuous(),
					      obsFile[i]->numDiscrete());
      assert(fileFilter);
      obsFile[i] = new FilterFile(fileFilter, obsFile[i], frs[i], irs[i], postpr[i]);
    }
  }

  if (nFiles == 0) return NULL;

  for (unsigned i=nFiles+1; i < MAX_NUM_OBS_FILES; i+=1) {
    if (ofs[i]) error("instantiateFileSource: Observation files [%u,%u] are missing\n", nFiles,i-1);
  }

  ObservationFile *mf;
  if (nFiles > 1) {
    mf = new  MergeFile(nFiles, obsFile,
			Action_If_Diff_Num_Sents,
			Action_If_Diff_Num_Frames,
			Ftr_Combo);
  } else {
    mf = obsFile[0];
  }
  ObservationFile *ff;
  if (Post_Transforms || gpr_str) {
    ff = new FilterFile(instantiateFilters(Post_Transforms,
			mf->numLogicalContinuous(), mf->numLogicalDiscrete()), 
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

