
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "GMTK_ASCIIStream.h"


/*
 * Written by Richard Rogers <rprogers@uw.edu>
 *
 * Copyright (C) 2012 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

RAND rnd(false);

// testASCIIStream file nf ni
int
main(int argc, char *argv[]) {
  if (argc != 4) {
    fprintf(stderr, "testASCIIStream file nf ni\n");
    exit (1);
  }
  FILE *f;

  if (strcmp("-", argv[1]))
    f = fopen(argv[1], "r");
  else
    f = stdin;

  if (!f) {
    perror(argv[1]);
    exit(1);
  }
  unsigned nf = atoi(argv[2]);
  unsigned ni = atoi(argv[3]);

  ASCIIStream as(f, nf, ni);

  Data32 const *frame;
  for (; !as.EOS(); ) {
    frame = as.getNextLogicalFrame();
    if (!frame) {
      printf("eos\n");
      continue;
    }
    for (unsigned i=0; i < nf; i+=1) {
      printf("%f ", ((float *)frame)[i]);
    }
    for (unsigned i=0; i < ni; i+=1) {
      printf("%d ", ((int *)frame)[as.numContinuous() + i]);
    }
    printf("\n");
  }
  exit(0);
}

