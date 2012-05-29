
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
    error("ERROR: ASCIIStream::getNextFrame: couldn't read stream tag in segment %u frame %u\n", currentSegNum, currentFrameNum);
  }
  if (tag == 'E' || tag == 'e') {
    currentFrameNum = 0;
    currentSegNum += 1;
    return NULL;
  }
  if (tag != 'F' && tag != 'f') {
    error("ERROR: ASCIIStream::getNextFrame: expected tag E or F, got '%c' in segment %u frame %u\n", tag, currentSegNum, currentFrameNum);
  }
  float *fdest = (float *)frameData;
  for (unsigned n=0; n < nFloat; n+=1) {
    if (fscanf(f, " %e ", fdest++) != 1) {
      error("ERROR: ASCIIStream::getNextFrame: couldn't read the %u'th continuous feature in segment %u frame %u\n", n, currentSegNum, currentFrameNum);
    }
  }
  Int32 *idest = (Int32 *)fdest;
  for (unsigned n=0; n < nInt; n+=1) {
    if (fscanf(f, " %d ", idest++) != 1) {
      error("ERROR: ASCIIStream::getNextFrame: couldn't read the %u'th discrete feature in segment %u frame %u\n", n, currentSegNum, currentFrameNum);
    }
  }
  currentFrameNum += 1;
  return frameData;
}

