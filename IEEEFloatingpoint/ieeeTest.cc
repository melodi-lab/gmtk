
#include <stdlib.h>
#include <cstdio>

#include "ieeeFPsetup.h"

int
main (int argc, char *argv[]) {
  ieeeFPsetup();
  float x = ((float)rand()/(float)RAND_MAX) / 0.0;
  printf("x = %f\n", x);
  exit(0);
}
