
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
#include "GMTK_Filter.h"
#include "GMTK_MergeFile.h"


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

FileSource *
instantiateFileSource() {
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
    assert(fileFilter);
    obsFile[i] = new FilterFile(fileFilter, obsFile[i], frs[i], irs[i], postpr[i]);
    nCont += obsFile[i]->numContinuous();
  }
  MergeFile *mf  = new  MergeFile(nFiles, obsFile,
				  Action_If_Diff_Num_Sents,
				  Action_If_Diff_Num_Frames,
				  Ftr_Combo);
  FilterFile *ff = new FilterFile(instantiateFilters(Post_Transforms, nCont), 
				  mf, NULL, NULL, gpr_str);
  unsigned windowBytes = fileWindowSize * MEBIBYTE;
  return new FileSource(ff, windowBytes, fileWindowDelta, fileBufferSize, 
			startSkip, endSkip, justification, constantSpace);
}


void 
instantiateFileSource(FileSource *source) {
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
    assert(fileFilter);
    obsFile[i] = new FilterFile(fileFilter, obsFile[i], frs[i], irs[i], postpr[i]);
    nCont += obsFile[i]->numContinuous();
  }
  MergeFile *mf  = new  MergeFile(nFiles, obsFile,
				  Action_If_Diff_Num_Sents,
				  Action_If_Diff_Num_Frames,
				  Ftr_Combo);
  FilterFile *ff = new FilterFile(instantiateFilters(Post_Transforms, nCont), 
				  mf, NULL, NULL, gpr_str);
  unsigned windowBytes = fileWindowSize * MEBIBYTE;
  source->initialize(ff, windowBytes, fileWindowDelta, fileBufferSize, 
		     startSkip, endSkip, justification, constantSpace);
}

