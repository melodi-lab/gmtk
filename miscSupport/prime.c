
/*  
 * Written by Richard Rogers <rprogers@ee.washington.edu> 
 *
 * Copyright (C) 2013 Richard Rogers 
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 */

#include "prime.h"

uint32_t 
mod_exp (uint32_t b, uint64_t n, uint32_t m) {
  uint64_t n_bits   = n, 
           b_2j_r   = b % m;
  uint32_t residue  = 1;
  
  for ( ; n_bits; n_bits >>= 1, b_2j_r = (b_2j_r * b_2j_r) % m) {
    if (n_bits & 1) {
      residue = (residue * b_2j_r) % m;
    }
  }
  return residue;
}

int
miller (uint32_t n, uint32_t b) {
  uint32_t s, t, j, n_bits;
  uint64_t pow;   /* pow is 2^j */

  for (s=0, n_bits = n-1; !(n_bits & 1); s += 1, n_bits >>= 1)
    ;
  t = n_bits;
  if (mod_exp (b, t, n) == 1) 
    return true;
  for (j = 0, pow = 1; j < s; j += 1, pow *= 2) {
    if ( mod_exp (b, pow * t, n) == n - 1 )
      return true;
  }
  return false;
}

int 
prime32 (uint32_t n) {
  if (n == 2 || n == 7 || n == 61) return true;  /* bases */
  if ( n < 2 || !(n & 1) )         return false; /* 0, 1, or even */
  if (!miller (n,  2))             return false; /* composite */
  if (!miller (n,  7))             return false;
  if (!miller (n, 61))             return false;

  return true; /* strong pseudoprime to bases 2, 7, & 61, thus prime
                  since all such 32-bit numbers are prime */
}
