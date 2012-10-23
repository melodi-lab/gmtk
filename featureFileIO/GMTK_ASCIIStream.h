
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

#include <string.h>
#include <stdio.h>
using namespace std;

#include "machine-dependent.h"
#include "error.h"

#include "GMTK_ObservationStream.h"

#define GMTK_ASC_COOKIE_LENGTH    6
#define GMTK_ASC_VERSION_LENGTH   6
#define GMTK_ASC_PROTOCOL_COOKIE  "GMTa\n"
#define GMTK_ASC_PROTOCOL_VERSION "0000\n"


class ASCIIStream: public ObservationStream {
  FILE *f;   // file to read data from

  char version[GMTK_ASC_VERSION_LENGTH]; // protocol version #
  
 public:

  ASCIIStream() {f=NULL;}
  
  ASCIIStream(FILE *file, unsigned nFloat, unsigned nInt, 
	     char const *contFeatureRangeStr=NULL, char const *discFeatureRangeStr=NULL) 
    : ObservationStream(nFloat, nInt, contFeatureRangeStr, discFeatureRangeStr), f(file)
  {
    char cookie[GMTK_ASC_COOKIE_LENGTH];
    if (fgets(cookie, GMTK_ASC_COOKIE_LENGTH, f) != cookie) {
      error("ERROR: ASCIIStream did not begin with 'GMTa\\n'\n");
    }
    if (strcmp(cookie, GMTK_ASC_PROTOCOL_COOKIE)) {
      error("ERROR: ASCIIStream did not begin with 'GMTa\\n'\n");
    }
    if (fgets(version, GMTK_ASC_VERSION_LENGTH, f) != version) {
      error("ERROR: ASCIIStream couldn't read protocol version\n");
    }
    if (strcmp(version, GMTK_ASC_PROTOCOL_VERSION) > 0) {
      version[GMTK_ASC_VERSION_LENGTH-2] = 0;
      error("ERROR: input ASCIIStream version %s is newer than this implementation's version %s",
	      version, GMTK_ASC_PROTOCOL_VERSION);
    }
    unsigned nCont, nDisc;
    if (fscanf(f, "%u %u", &nCont, &nDisc) != 2) {
      error("ERROR: ASCIIStream failed to read the number of continuous and discreate features");
    }
    if (nCont != nFloat) {
      error("ERROR: ASCIIStream contains %u continuous features, but expected %u",
	    nCont, nFloat);
    }
    if (nDisc != nInt) {
      error("ERROR: ASCIIStream contains %u discrete features, but expected %u",
	    nDisc, nInt);
    }
  }


  ~ASCIIStream() {
    if (f) fclose(f);
  }


  // Note that the end-of-file indicator is only set on a read
  // past the last byte in the file. The fscanf(" %c ",...) in
  // ASCIIStream::getNextFrame() seems to try to read past the
  // last E\n in the stream, thus setting the eof indicator. The
  // binary protocol has no white-space between the flag characters
  // and the data, so BinaryStream::EOS() is a little more
  // complicated

  bool EOS() {return feof(f);}

  Data32 const *getNextFrame();

};

#endif

