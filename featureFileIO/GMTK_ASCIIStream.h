
/*
 * GMTK_ASCIIStream.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2011, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */

#ifndef GMTK_ASCIISTREAM_H
#define GMTK_ASCIISTREAM_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
using namespace std;

#include "machine-dependent.h"
#include "error.h"

#include "GMTK_ObservationStream.h"

class ASCIIStream: public ObservationStream {
  FILE *f;   // file to read data from
  
 public:

  ASCIIStream() {f=NULL;}
  
  ASCIIStream(FILE *file, unsigned nFloat, unsigned nInt, 
	     char const *contFeatureRangeStr=NULL, char const *discFeatureRangeStr=NULL) 
    : ObservationStream(nFloat, nInt, contFeatureRangeStr, discFeatureRangeStr), f(file)
  {}


  ~ASCIIStream() {
    if (f) fclose(f);
  }


  // FIXME - see if it needs to fail on a read to set eof
  bool EOS() {return feof(f);}

  Data32 const *getNextFrame();

};

#endif

