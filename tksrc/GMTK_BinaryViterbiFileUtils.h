

/*
 * GMTK_BinaryViterbiFileUtils.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (C) 2013 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 * 
 * A few utility routines to simplify reading the binary Viterbi files
 *
 */

#ifndef GMTK_BINARYVITERBIFILEUTILS_H
#define GMTK_BINARYVITERBIFILEUTILS_H

#include "vbyteswapping.h"
#include "GMTK_SectionScheduler.h"

// only call this with integer types: [unsigned] {char,short,int,long}
template<typename Z>
int readVitZ(Z &x) {
  if (fread(&x, sizeof(Z), 1, SectionScheduler::binaryViterbiFile) != 1) {
    char *err = strerror(errno);
    error("ERROR: Unable to read binary Viterbi file '%s': %s\n",
	  SectionScheduler::binaryViterbiFilename, err);
  }
  if (SectionScheduler::binaryViterbiSwap) {
    switch (sizeof(Z)) {
    case 1: return 1;
    case 2: x = swapb_short_short(x); break;
    case 4: x = swapb_i32_i32(x); break;
    case 8: x = swapb_i64_i64(x); break;
    default: assert(0); // should be impossible
    }
  }
  return 1;
}


// only call this with real types: {float,double}
template<typename R>
int readVitR(R &x) {
  if (fread(&x, sizeof(R), 1, SectionScheduler::binaryViterbiFile) != 1) {
    char *err = strerror(errno);
    error("ERROR: Unable to read binary Viterbi file '%s': %s\n",
	  SectionScheduler::binaryViterbiFilename, err);
  }
  if (SectionScheduler::binaryViterbiSwap) {
    switch (sizeof(R)) {
    case 4: x = swapb_f32_f32(x); break;
    case 8: x = swapb_f64_f64(x); break;
    default: assert(0); // should be impossible
    }
  }
  return 1;
}

void readVitIntVector(size_t len, unsigned *ptr);

#endif
