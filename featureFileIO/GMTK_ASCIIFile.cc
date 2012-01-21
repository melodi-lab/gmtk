
/*
 * GMTK_ASCIIFile.cc
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

#include "GMTK_ASCIIFile.h"

using namespace std;

/**
 *  calcNumFileNames -- calculate the number of file names in the file pointed to by the file handle f
 *
 *  pre-conditions: the file name, fofName, must be initialized 
 *
 *  side effects: if f is a CPP pipe, the file is closed and re-opened
 *  in order to achieve the same effect as a rewind
 * */
unsigned
ASCIIFile::calcNumFileNames(FILE* &f) {
  char line[MAXSTRLEN];

  numFileNames = 0;
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

unsigned
ASCIIFile::readFof(FILE * &f) {

  unsigned n_lines = 0;

  numFileNames = 0;
  char line[MAXSTRLEN];

  numFileNames=calcNumFileNames(f);
  dataNames = new char*[numFileNames];

  n_lines = 0;
  while (fgets(line,sizeof(line),f) != NULL) {
    int l = strlen(line);
    if(l==0|| (l==1 && line[0]=='\n')) continue;
    if(line[0]==CPP_DIRECTIVE_CHAR) continue;
    if (line[l-1] != '\n') {
      if (n_lines < numFileNames-1) 
	error("ASCIFile::readFof: line %i too long in file '%s' - increase MAXSTRLEN (currently set to %d) or decrease your line lengths.\n",
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

ASCIIFile::ASCIIFile(const char *name, unsigned nfloats, unsigned nints, 
		     unsigned num, bool cppIfAscii, char* cppCommandOptions)
  : cppIfAscii(cppIfAscii),cppCommandOptions(cppCommandOptions)
{
  buffer = NULL;
  if (name == NULL) 	
    error("ASCIIFile: File name is NULL for stream %i\n",num);	
  if (nfloats == 0 && nints == 0)
    error("ASCIIFile: number of float and int features cannot both be zero");

  // local copy of file name
  fofName = new char[strlen(name)+1];
  strcpy(fofName,name);
  nFloats = nfloats;
  nInts = nints;
  
#ifdef PIPE_ASCII_FILES_THROUGH_CPP     
  if(cppIfAscii) { 
    string cppCommand = CPP_Command();
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
  } else {
    if ((fofFile = fopen(fofName,"r")) == NULL)
      error("ASCIIFile: Can't open '%s' for input\n",fofName);
  }
#else
  if ((fofFile = fopen(fofName,"r")) == NULL)
    error("ASCIIFile: Can't open '%s' for input\n",fofName);
#endif

  // for some reason this asserts while the following doesn't
  // assert(numFileNames == readFof(fofFile));
  
  unsigned rfof = readFof(fofFile);
  assert(numFileNames == rfof);

#ifdef PIPE_ASCII_FILES_THROUGH_CPP     
  if(cppIfAscii) {
    // first, scan until end of file since sometimes it appears
    // that cosing a pipe when not at the end causes an error (e.g., mac osx).
    freadUntilEOF(fofFile);
    if (pclose(fofFile) != 0) {
      // we don' give a warning here since sometimes 'cpp' might return
      // with an error that we really don't care about. TODO: the proper
      // thing to do here is no to use 'cpp' as a pre-processor and use
      // some other macro preprocessor (such as m4).
      // warning("WARNING: Can't close pipe '%s %s'.",CPP_Command());
    }
  }
  else
    fclose(fofFile);
#else
  fclose(fofFile);
#endif
  fofFile = NULL;
}


// Begin sourcing data from the requested segment.
// Must be called before any other operations are performed on a segment.
bool
ASCIIFile::openSegment(unsigned seg) {
  assert(seg < numFileNames);
  bool found_new_line=true;
  char *fname = dataNames[seg];
  unsigned n_samples = 0;
  char ch;
  FILE *curDataFile = NULL;

  if (fname == NULL) {
    warning("ASCIIFile::openSegment: Filename is NULL for segment %u\n",seg);
    return false;
  }

#ifdef PIPE_ASCII_FILES_THROUGH_CPP
  if(cppIfAscii) {
    string cppCommand = string("cpp");
    if (cppCommandOptions != NULL) {
      cppCommand = cppCommand + string(" ") + string(cppCommandOptions);
    }
    // make sure the file  exists first.
    if ((curDataFile = ::fopen(fname,"r")) == NULL) {
      error("ERROR: unable to open segment no %u (file '%s') for reading",seg,fname);
    }
    fclose(curDataFile);
    cppCommand = cppCommand + string(" ") + string(fname);
    curDataFile = ::popen(cppCommand.c_str(),"r");
    if (curDataFile == NULL)
      error("ERROR, can't open file (%s)",curDataFile);
    int tmp;
    while ((tmp = fgetc(curDataFile)) != EOF) {
      ch=tmp;
      if(ch=='\n') { found_new_line=true; continue;}
      if(ch==CPP_DIRECTIVE_CHAR) { found_new_line=false; continue; }
      if(found_new_line) {
	found_new_line=false;
	n_samples++;
      }
    }
    pclose(curDataFile);
    curDataFile = ::popen(cppCommand.c_str(),"r");
  } else {
    if ((curDataFile = fopen(fname,"r")) == NULL) {
      warning("ASCIIFile::openSegment: Can't open '%s' for input\n",fname);
      return false;
    }
    /* for ascii, newline is record delimiter - additional or missing nl's will cause error messages.  Fixed this -- karim 07sep2003*/
    int tmp;
    while ((tmp = fgetc(curDataFile)) != EOF) {
      ch = tmp;
      if(ch=='\n') { found_new_line=true; continue;}
      if(ch==CPP_DIRECTIVE_CHAR) { found_new_line=false; continue; }
      if(found_new_line) {
	found_new_line=false;
	n_samples++;
      }
    }
    rewind(curDataFile);
  }
#else
  if ((curDataFile = fopen(fname,"r")) == NULL) {
    warning("ASCIIFile::openSegment: Can't open '%s' for input\n",fname);
    return false;
  }
  /* for ascii, newline is record delimiter - additional or missing nl's will cause error messages */
  int tmp;
  while ((tmp = fgetc(curDataFile)) != EOF) {
    ch = tmp;
    if(ch=='\n') { found_new_line=true; continue;}
    if(ch==CPP_DIRECTIVE_CHAR) { found_new_line=false; continue; }
    if(found_new_line) {
      found_new_line=false;
      n_samples++;
    }
  }
  rewind(curDataFile);
#endif
  nFrames = n_samples;
  
  if (buffer) delete [] buffer;
  buffer = new Data32[n_samples * (nFloats + nInts)];
  Data32 *dest = buffer;

  int lineNum=0;


  // consume CPP special directives if any
#ifdef PIPE_ASCII_FILES_THROUGH_CPP
  int tmp;
  if(cppIfAscii) {
    while((tmp=fgetc(curDataFile))==CPP_DIRECTIVE_CHAR) {
      while((tmp=fgetc(curDataFile))!='\n');
      lineNum++;
    }
    ungetc(tmp,curDataFile);
  }
#endif



  // could be made a bit more efficient since we check whether
  // num_floats and num_ints > 0 for each frame.
  for(unsigned s=0; s < n_samples; ++s) {
    lineNum++;
    for (unsigned n = 0; n < nFloats; n+=1) {
      if (fscanf(curDataFile,"%e", (float *)(dest++)) != 1) {
	error("ERROR: ASCIIFile::openSegment: couldn't read %u'th item in frame %u\n",n,s);
      }
    }
    for (unsigned n = 0; n < nInts; n+=1) {
      if (fscanf(curDataFile,"%d", (Int32 *)(dest++)) != 1) {
	error("ERROR: ASCIIFile::openSegment: couldn't read %u'th item in frame %u\n",n,s);
      }
    }
  }
  fclose(curDataFile);
  curDataFile = NULL;
  return true;
}

