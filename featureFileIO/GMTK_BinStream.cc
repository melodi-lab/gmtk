
/*
 * GMTK_BinStream.cc
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2012, < fill in later >
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

#include "vbyteswapping.h"
#include "GMTK_BinStream.h"

Data32 const*
BinaryStream::getNextFrame() {
  char tag;
  if (fscanf(f,"%c", &tag) != 1) {
#if 0
    if (feof(f)) return NULL;
#endif
    error("ERROR: BinaryStream::getNextFrame: couldn't read stream tag");
  }
  if (tag == 'E' || tag == 'e') return NULL;
  if (tag != 'F' && tag != 'f') {
    error("ERROR: BinaryStream::getNextFrame: expected tag E or F, got '%c'", tag);
  }
  float *fdest = (float *)frameData;
  if (nFloat > 0) {
    if (fread(fdest, sizeof(Data32), nFloat, f) != nFloat) {
      error("ERROR: BinaryStream::getNextFrame: couldn't read the continuous features");
    }
    if (swap) swapb_vf32_vf32(nFloat, fdest, fdest);
  }
  Int32 *idest = (Int32 *)(fdest + nFloat);
  if (nInt > 0) {
    if (fread(idest, sizeof(Data32), nInt, f) != nInt) {
      error("ERROR: BinaryStream::getNextFrame: couldn't read the discrete features");
    }
    if (swap) swapb_vi32_vi32(nInt, idest, idest);
  }
  return frameData;
}


bool
BinaryStream::EOS() {
  assert(f);
  int c = fgetc(f);         // read a byte to possilby set eof indicator
  bool eof = feof(f) != 0;  // eof?
  ungetc(c, f);             // return the read byte to the stream. may fail if @ eof, but so what?
  return eof;
}
