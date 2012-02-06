
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

#include "GMTK_PFileFile.h"
#include "arguments.h"

#define MAX_OBJECTS 10

char *input_fname[MAX_OBJECTS] = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};  // Input file name.
bool iswp[MAX_OBJECTS] = {false,false,false,false,false};
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
  Arg("iswp",Arg::Opt,iswp,"Endian swap condition for observation file X",Arg::ARRAY,MAX_OBJECTS),

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
  for (unsigned i=0; i < MAX_OBJECTS && input_fname[i] != NULL; i+=1) {
    PFileFile f(input_fname[i], numFloat, numInt, i, iswp[i]);
    printf("reading %s - %u\n", input_fname[i], f.numSegments());
    for (unsigned j=0; j < f.numSegments(); j+=1) {
      assert(f.openSegment(j));
      for (unsigned k=0; k < f.numFrames(); k+=1) {
	printf("%03u %03u", j, k);
	Data32 const *buf = f.getFrames(k,1);
	for (unsigned f=0; f < numFloat; f+=1)
	  printf(" %f", *((float *)(buf++)));
	for (unsigned f=0; f < numInt; f+=1)
	  printf(" %d", *((int *)(buf++)));
	printf("\n");
      }
    }
  }
  exit(0);
}
