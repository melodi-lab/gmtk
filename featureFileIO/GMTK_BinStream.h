
/*
 * GMTK_BinStream.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2012, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */

#ifndef GMTK_BINSTREAM_H
#define GMTK_BINSTREAM_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include <stdio.h>
using namespace std;

#include "machine-dependent.h"
#include "vbyteswapping.h"
#include "error.h"

#include "GMTK_WordOrganization.h"
#include "GMTK_ObservationStream.h"

#define GMTK_BIN_COOKIE_LENGTH    6
#define GMTK_BIN_VERSION_LENGTH   6
#define GMTK_BIN_PROTOCOL_COOKIE  "GMTb\n"
#define GMTK_BIN_PROTOCOL_VERSION "0000\n"

class BinaryStream: public ObservationStream {
  FILE *f;     // file to read data from
  bool swap;   // true if we need to swap to match the requested byte order

  char version[GMTK_BIN_VERSION_LENGTH]; // protocol version #
  
 public:

  BinaryStream() {f=NULL;}
  
  BinaryStream(FILE *file, unsigned nFloat, unsigned nInt,
	       char const *contFeatureRangeStr=NULL, char const *discFeatureRangeStr=NULL,
	       bool netByteOrder=true) 
    : ObservationStream(nFloat, nInt, contFeatureRangeStr, discFeatureRangeStr), f(file)
  {
    char cookie[GMTK_BIN_COOKIE_LENGTH];
    if (fgets(cookie, GMTK_BIN_COOKIE_LENGTH, f) != cookie) {
      error("ERROR: BinaryStream did not begin with 'GMTb\\n'");
    }
    if (strcmp(cookie, GMTK_BIN_PROTOCOL_COOKIE)) {
      error("ERROR: BinaryStream did not begin with 'GMTb\\n'");
    }
    if (fgets(version, GMTK_BIN_VERSION_LENGTH, f) != version) {
      error("ERROR: BinaryStream couldn't read protocol version");
    }
    if (strcmp(version, GMTK_BIN_PROTOCOL_VERSION) > 0) {
      version[GMTK_BIN_VERSION_LENGTH-2] = 0;
      error("ERROR: input BinaryStream version %s is newer than this implementation's version %s",
              version, GMTK_BIN_PROTOCOL_VERSION);
    }
    // TODO - send BOM instead of risking being wrong?
    swap = ( netByteOrder && getWordOrganization() != BYTE_BIG_ENDIAN)    ||
           (!netByteOrder && getWordOrganization() != BYTE_LITTLE_ENDIAN);
    unsigned nCont, nDisc;
    if (fread(&nCont, sizeof(nCont), 1, f) != 1) {
      error("ERROR: BinaryStream::couldn't read the number of continuous features");
    }
    if (swap) nCont = swapb_i32_i32(nCont);
    if (nCont != nFloat) {
      error("ERROR: BinaryStream contains %u continuous features, but expected %u",
	    nCont, nFloat);
    }
    if (fread(&nDisc, sizeof(nDisc), 1, f) != 1) {
      error("ERROR: BinaryStream::couldn't read the number of continuous features");
    }
    if (swap) nDisc = swapb_i32_i32(nDisc);
    if (nDisc != nInt) {
      error("ERROR: BinaryStream contains %u discrete features, but expected %u",
	    nDisc, nInt);
    }
  }


  ~BinaryStream() {
    if (f) fclose(f);
  }


  bool EOS();

  Data32 const *getNextFrame();

};

#endif

