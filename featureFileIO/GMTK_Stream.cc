/*
 *
 * Written by Katrin Kirchhoff <katrin@ee.washington.edu>
 *
 * Modified by Karim Filali <karim@cs.washington.edu> to add the
 * option to pipe the list of file names through CPP.  Made a few
 * other minor bug fixes.
 *
 * $Header$
 *
 * Copyright (c) 2001
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about the suitability of this soft
 * for any purpose. It is provided "as is" without express or implied warranty.
 * */

#include <stdio.h>
#include <string.h>

#include "GMTK_Stream.h"

#include "error.h"
#include "general.h"

// karim - 29aug2003
#include "fileParser.h"
#ifdef PIPE_ASCII_FILES_THROUGH_CPP
#ifndef DECLARE_POPEN_FUNCTIONS_EXTERN_C
extern "C" {
  FILE     *popen(const char *, const char *) __THROW;
  int pclose(FILE *stream) __THROW;
};
#endif
#endif
#ifdef PIPE_ASCII_FILES_THROUGH_CPP
#define CPP_DIRECTIVE_CHAR '#'
#endif

StreamInfo::StreamInfo(const char *name, const char *crng_str,
		       const char *drng_str,
		       unsigned *nfloats, unsigned *nints, 
		       unsigned *format, bool swap, unsigned num,bool cppIfAscii,char* cppCommandOptions) 
  : cppIfAscii(cppIfAscii),cppCommandOptions(cppCommandOptions),numFileNames(0),pfile_istr(NULL),dataNames(NULL),cont_rng(NULL),disc_rng(NULL)
{

  if (name == NULL) 	
    error("StreamInfo: File name is NULL for stream %i\n",num);	
  
  // local copy of file name
  fofName = new char[strlen(name)+1];
  strcpy(fofName,name);
  
  if (format == NULL)
    error("StreamInfo: Data format unspecified for stream %i\n",num);
   
   dataFormat = *format;
   
   if (nfloats == NULL)
     error("StreamInfo: Number of floats unspecified for stream %i\n",num);

   nFloats = *nfloats;

   if (nints == NULL) 
      error("StreamInfo: Number of ints unspecified for stream %i\n",num);

   nInts = *nints;

   // Update: this is already taken car of when {c|d}rng_str eq "nil" or "none"
   ////// Special cases {c|d}rng_str eq "-1"
   //if(crng_str != NULL && strcmp(crng_str,"-1")==0) {
     // empty range
   //cont_rng = new BP_Range(NULL,0,0);
   //}
   //else {
   cont_rng = new BP_Range(crng_str,0,nFloats);
   //}
   //if(drng_str != NULL && strcmp(drng_str,"-1")==0) {
   //disc_rng = new BP_Range(NULL,0,0);  // empty range
   //}
   //else {
   disc_rng = new BP_Range(drng_str,0,nInts);
   //}

   // If in the future we want to have the option to append deltas and
   // double deltas, it should be taken into account below.

   nFloatsUsed = cont_rng->length();
   nIntsUsed = disc_rng->length();

   if( (nFloatsUsed + nIntsUsed) == 0) {
       warning("WARNING: No features were selected for stream %i\n",num);
   }

   bswap = swap;

   if (dataFormat == PFILE) {

     if ((curDataFile = fopen(fofName,"rb")) == NULL)
       error("StreamInfo: Can't open '%s' for input\n", fofName);

     pfile_istr = new InFtrLabStream_PFile(0,fofName,curDataFile,1,bswap);
     
     fofSize = pfile_istr->num_segs();
     
     if (pfile_istr->num_ftrs() != nFloats) 
       error("StreamInfo: File %s has %i floats, expected %i",
	     fofName,
	     pfile_istr->num_ftrs(),
	     nFloats);
     
     // in a pfile, only the labs can be ints
     
     if (pfile_istr->num_labs() != nInts)
       error("StreamInfo: File %s has %i ints, expected %i",
	     fofName,
	     pfile_istr->num_labs(),
	     nInts);
   }
   else {
     pfile_istr = NULL;
#ifdef PIPE_ASCII_FILES_THROUGH_CPP     
     if(cppIfAscii) { 
       string cppCommand = string("cpp");
       if (cppCommandOptions != NULL) {
	 cppCommand = cppCommand + string(" ") + string(cppCommandOptions);
       }
       // make sure the file  exists first.
       if ((fofFile = ::fopen(fofName,"r")) == NULL) {
	 error("ERROR: unable to open file (%s) for reading",fofName);
       }
       fclose(fofFile);
       cppCommand = cppCommand + string(" ") + string(fofName);
       fofFile = ::popen(cppCommand.c_str(),"r");    
       if (fofFile == NULL)
	 error("ERROR, can't open file stream from (%s)",fofName);
     }
     else {
       if ((fofFile = fopen(fofName,"r")) == NULL)
	 error("StreamInfo: Can't open '%s' for input\n",fofName);
     }
#else
      if ((fofFile = fopen(fofName,"r")) == NULL)
	 error("StreamInfo: Can't open '%s' for input\n",fofName);
#endif

     fofSize = readFof(fofFile);
#ifdef PIPE_ASCII_FILES_THROUGH_CPP     
     if(cppIfAscii) {
       if (pclose(fofFile) != 0)
	 warning("WARNING: Can't close pipe 'cpp %s'.",fofName);
     }
     else
       fclose(fofFile);
#else
     fclose(fofFile);
#endif
   }
}


StreamInfo::~StreamInfo() {

  if (fofName != NULL)
    delete [] fofName;

  if (dataFormat != PFILE) {
    if (dataNames != NULL) {
      for (unsigned i = 0; i < fofSize; i++)
	delete [] dataNames[i];
      delete [] dataNames;
    }
  }
  else 
    delete pfile_istr;

  if (cont_rng != NULL)
    delete cont_rng;   
  if (disc_rng != NULL)
    delete disc_rng;
}



/**
 *  calcNumFileNames -- calculate the number of file names in the file pointed to by the file handle f
 *
 *  pre-conditions: the file name, fofName, must be initialized 
 *
 *  side effects: if f is a CPP pipe, the file is closed and re-opened
 *  in order to achieve the same effect as a rewind
 * */
size_t StreamInfo::calcNumFileNames(FILE* &f) {
  char line[MAXSTRLEN];

  numFileNames = 0;
#ifdef PIPE_ASCII_FILES_THROUGH_CPP     
  if(cppIfAscii) {
    while (fgets(line,sizeof(line),f) != NULL) {
      int l = strlen(line);
      if(l==0 || l==1 && line[0]=='\n') continue;
      if(line[0]==CPP_DIRECTIVE_CHAR) continue;  // lines that start with # are CPP directives 
      numFileNames++;
    }
    // since it's a pipe we need to close it and reopen it
    fclose(f);
    string cppCommand = string("cpp");
    if (cppCommandOptions != NULL) {
      cppCommand = cppCommand + string(" ") + string(cppCommandOptions);
    }
    cppCommand = cppCommand + string(" ") + string(fofName);
    f = ::popen(cppCommand.c_str(),"r");    
    if (f == NULL)
      error("ERROR, can't open file stream from (%s)",fofName);
  } 
  else {
    while (fgets(line,sizeof(line),f) != NULL) {
      int l = strlen(line);
      if(l==0 || (l==1 && line[0]=='\n') || line[0]==CPP_DIRECTIVE_CHAR) continue;
      numFileNames++;
    }
    rewind(f);
  }
#else
   while (fgets(line,sizeof(line),f) != NULL) {
      int l = strlen(line);
      if(l==0 || l==1 && line[0]=='\n') continue;
      numFileNames++;
    }
    rewind(f);
#endif

  return numFileNames;

}	

size_t StreamInfo::readFof(FILE *f) {

  size_t n_lines = 0;

  numFileNames = 0;
  char line[MAXSTRLEN];

  numFileNames=calcNumFileNames(f);
  dataNames = new char*[numFileNames];

  n_lines = 0;
  while (fgets(line,sizeof(line),f) != NULL) {
    int l = strlen(line);
    if(l==0|| l==1 && line[0]=='\n') continue;
    if(line[0]==CPP_DIRECTIVE_CHAR) continue;
    if (line[l-1] != '\n') {
      if (n_lines < numFileNames-1) 
	error("StreamInfo::readFof: line %i too long in file '%s' - increase MAXSTRLEN\n",
	      n_lines+1,fofName);
    }
    else
      line[l-1] = '\0';
    
    dataNames[n_lines] = new char[l+1];
    strcpy(dataNames[n_lines],line);
    n_lines++;
  }
  assert(numFileNames==n_lines);
  return n_lines;
}




	








