
#include <stdlib.h>
#include <stdio.h>

#include "GMTK_ASCIIStream.h"

// testASCIIStream file nf ni
int
main(int argc, char *argv[]) {
  if (argc != 4) {
    fprintf(stderr, "testASCIIStream file nf ni\n");
    exit (1);
  }
  FILE *f = fopen(argv[1], "r");
  if (!f) {
    perror(argv[1]);
    exit(1);
  }
  unsigned nf = atoi(argv[2]);
  unsigned ni = atoi(argv[3]);

  ASCIIStream as(f, nf, ni);

  Data32 const *frame;
  for (; !as.EOS(); ) {
    frame = as.getNextFrame();
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

