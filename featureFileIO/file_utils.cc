
/*
 * file_utils.cc
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2011, < fill in later >
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
#include "file_utils.h"

using namespace std;

/**
 * openCPPableFile -- open an ASCII file that may need to be preprocessed by CPP
 *
 * returns NULL on failure
 */
FILE *
openCPPableFile(char const *filename, bool cppIfAscii, 
		char const *cppCommandOptions) 
{
  FILE *f = NULL;
#ifdef PIPE_ASCII_FILES_THROUGH_CPP     
  if(cppIfAscii) { 
    string cppCommand = CPP_Command();
    if (cppCommandOptions != NULL) {
      cppCommand = cppCommand + string(" ") + string(cppCommandOptions);
    }
    // make sure the file  exists first.
    if ((f = ::fopen(filename,"r")) == NULL) {
      warning("openCPPableFile: unable to open file (%s) for reading",filename);
    }
    fclose(f);
    cppCommand = cppCommand + string(" ") + string(filename);
    f = ::popen(cppCommand.c_str(),"r");    
    if (f == NULL)
      warning("openCPPableFile: can't open file stream from (%s)",filename);
  } else {
    if ((f = fopen(filename,"r")) == NULL)
      warning("openCPPableFile: Can't open '%s' for input\n",filename);
  }
#else
  if ((f = fopen(filename,"r")) == NULL)
    warning("oepnCPPableFile: Can't open '%s' for input\n",filename);
#endif
  return f;
}


/**
 * closeCPPableFile -- close a file opened by openCPPableFile
 *
 */
void
closeCPPableFile(FILE * &f, bool cppIfAscii) {
#ifdef PIPE_ASCII_FILES_THROUGH_CPP     
  if(cppIfAscii) {
    // first, scan until end of file since sometimes it appears
    // that cosing a pipe when not at the end causes an error (e.g., mac osx).
    freadUntilEOF(f);
    if (pclose(f) != 0) {
      // we don' give a warning here since sometimes 'cpp' might return
      // with an error that we really don't care about. TODO: the proper
      // thing to do here is no to use 'cpp' as a pre-processor and use
      // some other macro preprocessor (such as m4).
      // warning("WARNING: Can't close pipe '%s %s'.",CPP_Command());
    }
  }
  else
    fclose(f);
#else
  fclose(f);
#endif
  f = NULL;
}


/**
 *  calcNumFileNames -- calculate the number of file names in the 
 *                      file pointed to by the file handle f
 *
 *  pre-conditions: the file f and its file name, fofName, must be initialized 
 *
 *  side effects: if f is a CPP pipe, the file is closed and re-opened
 *  in order to achieve the same effect as a rewind
 *
 */
unsigned 
calcNumFileNames(FILE* &f, char const *fofName, 
		 bool cppIfAscii, char const *cppCommandOptions) 
{
  char line[MAXSTRLEN];

  unsigned numFileNames = 0;
#ifdef PIPE_ASCII_FILES_THROUGH_CPP     
  if(cppIfAscii) {
    while (fgets(line,sizeof(line),f) != NULL) {
      int l = strlen(line);
      if(l==0 || (l==1 && line[0]=='\n')) continue;
      if(line[0]==CPP_DIRECTIVE_CHAR) continue;  // lines that start with # are CPP directives 
      numFileNames++;
    }
    // since it's a pipe we need to close it and reopen it
    pclose(f);
    string cppCommand = CPP_Command();
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
      if(l==0 || (l==1 && line[0]=='\n') ) continue;
      numFileNames++;
    }
    rewind(f);
#endif

  return numFileNames;
}	


/**
 * readFof -- reads the list of filenames in f into dataNames,
 *            stores the filename count in numFileNames
 * 
 *  pre-conditions: the file f and its file name, fofName, must be initialized 
 *
 *  side effects: if f is a CPP pipe, the file is closed and re-opened
 *  in order to achieve the same effect as a rewind
 *
 */
unsigned
readFof(FILE * &f, char const *fofName, unsigned &numFileNames, char **&dataNames, 
	bool cppIfAscii, char const *cppCommandOptions) 
{

  unsigned n_lines = 0;

  numFileNames = 0;
  char line[MAXSTRLEN];

  numFileNames=calcNumFileNames(f, fofName, cppIfAscii, cppCommandOptions);
  dataNames = new char*[numFileNames];

  n_lines = 0;
  while (fgets(line,sizeof(line),f) != NULL) {
    int l = strlen(line);
    if(l==0|| (l==1 && line[0]=='\n')) continue;
    if(line[0]==CPP_DIRECTIVE_CHAR) continue;
    if (line[l-1] != '\n') {
      if (n_lines < numFileNames-1) 
	error("readFof: line %i too long in file '%s' - increase MAXSTRLEN (currently set to %d) or decrease your line lengths.\n",
	      n_lines+1,fofName,MAXSTRLEN);
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

