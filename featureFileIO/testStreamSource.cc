
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

#include "GMTK_StreamSource.h"
#include "GMTK_ASCIIStream.h"
#include "GMTK_Filter.h"
#include "GMTK_FilterFile.h"
#include "GMTK_FIRFilter.h"
#include "GMTK_AffineFilter.h"
#include "arguments.h"

unsigned int numInt;
unsigned int numFloat;
char *strans  = NULL;
char *frs_str = NULL;
char *irs_str = NULL;

bool printVersion = false;
unsigned help=0;  // help=0...5 depending on the amount of info we want printed
Arg Arg::Args[] = {

  Arg("\n*** Input arguments ***\n"),

  Arg("nf",   Arg::Req, numFloat,"number of floats in input file(s)"),
  Arg("ni",   Arg::Req, numInt,"number of ints (labels) in input file(s)"),
  Arg("strans", Arg::Opt, strans, "Stream transform"),
  Arg("frs", Arg::Opt, frs_str, "Float feature range"),
  Arg("irs", Arg::Opt, irs_str, "Int feature range"),
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

  //  Filter *filt = NULL;
  ASCIIStream as(stdin, numFloat, numInt, frs_str, irs_str);
  ObservationStream *os[1];
  os[0] = &as;
  //  if (strans) filt = instantiateFilters(strans, as.numLogicalContinuous());
  StreamSource ss(1U, os, 100U, strans);

  unsigned segNum, frameNum;
  for (segNum=0; !ss.EOS(); segNum+=1) {
    ss.preloadFrames(3);
    frameNum = 0;
    do {
      Data32 const *frame = ss.loadFrames(frameNum, 1);
      if (ss.numFrames() > 0 && frameNum >= ss.numFrames()) {
	break;
      }
      
      printf("%03u %03u", segNum, frameNum);
      for (unsigned f=0; f < ss.numContinuous(); f+=1)
	printf(" %f", *((float *)(frame++)));
      for (unsigned f=0; f < ss.numDiscrete(); f+=1)
	printf(" %d", *((int *)(frame++)));
      printf("\n");
      frameNum += 1;
    } while (ss.numFrames() == 0 || frameNum < ss.numFrames());
  }
  exit(0);
}
