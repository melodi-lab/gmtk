
#include <stdio.h>

#include "sArray.h"
#include "rand.h"

#define testSize 20

int
main(int argc, char *argv[]) {
  RAND r;
  sArray<int> a(testSize), b(testSize), c(testSize);
  for (unsigned i=0; i < testSize; i+=1) {
    int x = r.uniform(30);
    a[i] = x; b[i] = x; c[i] = x;
  }

  b.qsort();
  c.sort();

  for (unsigned i=0; i < testSize; i+=1) {
    printf("%2d   %2d   %2d\n", a[i], b[i], c[i]);
  }
  exit(0);
}
