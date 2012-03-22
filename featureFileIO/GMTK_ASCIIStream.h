
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
#include "GMTK_StreamCookie.h"

class ASCIIStream: public ObservationStream {
  FILE *f;   // file to read data from

  char version[GMTK_VERSION_LENGTH]; // protocol version #
  
 public:

  ASCIIStream() {f=NULL;}
  
  ASCIIStream(FILE *file, unsigned nFloat, unsigned nInt, 
	     char const *contFeatureRangeStr=NULL, char const *discFeatureRangeStr=NULL) 
    : ObservationStream(nFloat, nInt, contFeatureRangeStr, discFeatureRangeStr), f(file)
  {
    char cookie[GMTK_COOKIE_LENGTH];
    if (fgets(cookie, GMTK_COOKIE_LENGTH, f) != cookie) {
      error("ERROR: ASCIIStream did not begin with 'GMTK\\n'\n");
    }
    if (strcmp(cookie, GMTK_PROTOCOL_COOKIE)) {
      error("ERROR: ASCIIStream did not begin with 'GMTK\\n'\n");
    }
    if (fgets(version, GMTK_VERSION_LENGTH, f) != version) {
      error("ERROR: ASCIIStream couldn't read protocol version\n");
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

