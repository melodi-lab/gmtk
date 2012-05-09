
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "GMTK_BinStream.h"

// testBinStream file nf ni netorder
int
main(int argc, char *argv[]) {
  if (argc != 5) {
    fprintf(stderr, "testBinStream file nf ni no\n");
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
  bool     no = atoi(argv[4]) != 0;

  BinaryStream bs(f, nf, ni, NULL, NULL, no);

  Data32 const *frame;
  for (; !bs.EOS(); ) {
    frame = bs.getNextLogicalFrame();
    if (!frame) {
      printf("eos\n");
      continue;
    }
    for (unsigned i=0; i < nf; i+=1) {
      printf("%f ", ((float *)frame)[i]);
    }
    for (unsigned i=0; i < ni; i+=1) {
      printf("%d ", ((int *)frame)[bs.numContinuous() + i]);
    }
    printf("\n");
  }
  exit(0);
}

