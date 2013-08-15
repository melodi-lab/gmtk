

/*
 * GMTK_BinaryViterbiFileUtils.cc
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (C) 2013 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 * 
 *
 */


#include <stdio.h>
#include "GMTK_BinaryViterbiFileUtils.h"

void readVitIntVector(size_t len, unsigned *ptr) {
  if (fread(ptr, sizeof(unsigned), len, JunctionTree::binaryViterbiFile) != len) {
    char *err = strerror(errno);
    error("ERROR: faild to read Viterbi values from '%s': %s\n", 
	  JunctionTree::binaryViterbiFilename, err);
  }
  if (JunctionTree::binaryViterbiSwap) {
    swapb_vi32_vi32(len, (int *)ptr, (int *)ptr);
  }
}
