
#include <stdlib.h>
#include <cstdio>
#include <fenv.h>

#include "ieeeFPsetup.h"

volatile float zero = 0.0;
volatile float five = 5.0;

int
main (int argc, char *argv[]) {
  ieeeFPsetup();
#if 0
  feraiseexcept(FE_OVERFLOW);
#else
  float x = five / zero;
  printf("x = %f\n", x);
#endif
  exit(0);
}
