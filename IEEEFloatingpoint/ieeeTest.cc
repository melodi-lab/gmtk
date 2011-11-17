
#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <cstdio>

#if HAVE_FENV_H
#include <fenv.h>
#endif

#include "ieeeFPsetup.h"

volatile float zero  = 0.0f;
volatile float one   = 1.0f;
volatile float three = 3.0f;

int
main (int argc, char *argv[]) {
  ieeeFPsetup();
  int which = 0;
  if (argc > 1)
    which = atoi(argv[1]);
  float x;
  switch (which) {
  case 1: x = one / three; break;        // FE_INEXACT
  case 2: x = one / zero;  break;        // FE_DIVZERO
  case 3: for (x=0.5f; 1; x*=x) ; break; // FE_UNDERFLOW
  case 4: for (x=2.0f; 1; x*=x) ; break; // FE_OVERFLOW
  case 5: x = zero / zero; break;        // FE_INVALID
  default: x=0.0f;
  }
  printf("x = %f\n", x);
  exit(0);
}
