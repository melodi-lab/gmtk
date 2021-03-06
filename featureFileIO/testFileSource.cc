
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

#include "GMTK_FileSource.h"
#include "GMTK_ASCIIFile.h"
#include "arguments.h"

#define MAX_OBJECTS 10

char *input_fname[MAX_OBJECTS] = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};  // Input file name.
unsigned int numInt;
unsigned int numFloat;
char *cppCommandOptions = NULL;
bool printVersion = false;
unsigned help=0;  // help=0...5 depending on the amount of info we want printed
Arg Arg::Args[] = {

  Arg("\n*** Input arguments ***\n"),

  Arg("i",    Arg::Req, input_fname,"input file. Replace X with the file number",Arg::ARRAY,MAX_OBJECTS),
  Arg("nf",   Arg::Req, numFloat,"number of floats in input file(s)"),
  Arg("ni",   Arg::Req, numInt,"number of ints (labels) in input file(s)"),
  Arg("cppCommandOptions",Arg::Opt,cppCommandOptions,"Command line options to give to 'cpp'"),

  Arg("\n*** Misc arguments ***\n"),

  Arg("help",  Arg::Help, help,  "Print this message. Add an argument from 1 to 5 for increasing help info."),
  Arg("version", Arg::Tog, printVersion, "Print GMTK version and exit."),
  // The argumentless argument marks the end of the above list.
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
  if (printVersion) {
#ifdef HAVE_CONFIG_H
    printf("%s (Mercurial id: %s)\n",gmtk_version_id,HGID);
#else
    printf("%s\n", gmtk_version_id);
#endif
    exit(0);
  }
  ObservationFile *file[MAX_OBJECTS];

  unsigned nfiles =0;
  unsigned i;
  for (i=0; i < MAX_OBJECTS && input_fname[i] != NULL; i+=1) {
    file[i] = new ASCIIFile(input_fname[i], numFloat, numInt, i, cppCommandOptions);
    fprintf(stderr, "file %u: %s\n", i, input_fname[i]);
  }
  nfiles = i;
  fprintf(stderr, "%u files\n", nfiles);
  FileSource fs(nfiles, file);

#define chunksize 5

  for (unsigned j=0; j < fs.numSegments(); j+=1) {
    assert(fs.openSegment(j));
    unsigned numFrames = fs.numFrames();

    unsigned k;
    for (k=0; k + chunksize < numFrames; ) {
      Data32 const *buf = fs.loadFrames(k,chunksize);
      if (buf) {
	for (i=0; i < chunksize; i+=1, k+=1) {
	  printf("%03u %03u", j, k);
	  for (unsigned f=0; f < fs.numContinuous(); f+=1)
	    printf(" %f", *((float *)(buf++)));
	  for (unsigned f=0; f < fs.numDiscrete(); f+=1)
	    printf(" %d", *((int *)(buf++)));
	  printf("\n");
	}
	//	printf("--------------------------------\n");
      } else {
	printf(" NULL");
      }
    }
    for ( ; k < numFrames; k+=1) {
      Data32 const *buf = fs.loadFrames(k,1);
      printf("%03u %03u", j, k);
      if (buf) {
	for (unsigned f=0; f < fs.numContinuous(); f+=1)
	  printf(" %f", *((float *)(buf++)));
	for (unsigned f=0; f < fs.numDiscrete(); f+=1)
	  printf(" %d", *((int *)(buf++)));
      } else {
	printf(" NULL");
      }
      printf("\n");
    }

  }
  exit(0);
}
