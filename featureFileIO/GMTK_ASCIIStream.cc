
/*
 * GMTK_ASCIIStream.cc
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
#include <string>
using namespace std;


#include "GMTK_ASCIIStream.h"

Data32 const*
ASCIIStream::getNextFrame() {
  char tag;
  if (fscanf(f," %c ", &tag) != 1) {
    error("ERROR: ASCIIStream::getNextFrame: couldn't read stream tag\n");
  }
  if (tag == 'E' || tag == 'e') return NULL;
  if (tag != 'F' && tag != 'f') {
    error("ERROR: ASCIIStream::getNextFrame: expected tag E or F, got '%c'\n", tag);
  }
  float *fdest = (float *)frameData;
  for (unsigned n=0; n < nFloat; n+=1) {
    if (fscanf(f, " %e ", fdest++) != 1) {
      error("ERROR: ASCIIStream::getNextFrame: couldn't read the %u'th continuous feature\n", n);
    }
  }
  Int32 *idest = (Int32 *)fdest;
  for (unsigned n=0; n < nInt; n+=1) {
    if (fscanf(f, " %d ", idest++) != 1) {
      error("ERROR: ASCIIStream::getNextFrame: couldn't read the %u'th discrete feature\n", n);
    }
  }
  return frameData;
}

