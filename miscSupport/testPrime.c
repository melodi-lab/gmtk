
/*  
 * Written by Richard Rogers <rprogers@ee.washington.edu> 
 *
 * Copyright (C) 2013 Richard Rogers 
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "prime.h"

#define SIEVE_MAX (UINT32_C(1)<<27)
#if 0
#define SIEVE_MAX (UINT32_C(1)<<22)
#define SIEVE_MAX (UINT32_C(1)<<24)
#define SIEVE_MAX (UINT32_C(3000000000))
#endif

int
main(int argc, char *argv[]) {
  uint8_t  *sieve   = (uint8_t *) calloc(SIEVE_MAX / 8 + 1, 1);
  uint8_t   mask[8] = {1,2,4,8,16,32,64,128};
  uint32_t  n;
  uint64_t  i;
  uint32_t  root_max = 1 + (uint32_t)sqrt(SIEVE_MAX);
  int       pass = true;

  if (!sieve) {
    perror ("failed to allocate sieve array");
    exit (1);
  }
  for (n=2; n < SIEVE_MAX; n+=1) {
    if (n % (SIEVE_MAX / 78) == 0) {
      printf ("*"); fflush(stdout);
    }
    if (sieve[n / 8] & mask[n % 8]) {
      if (prime32(n)) {
        printf ("\n%"PRIu32" falsely reported prime by prime32()!\n", n);
        pass = false;
      }
    } else {
      if (!prime32(n)) {
        printf ("\n%"PRIu32" not reported prime by prime32()!\n", n);
        pass = false;
      }
      if (n < root_max)
        for (i = n * 2; i < SIEVE_MAX; i+=n)
          sieve[i / 8] |= mask[i % 8];
    }
  }
  printf ("\n");
  free (sieve);
  return pass ? 0 : 1;
}

