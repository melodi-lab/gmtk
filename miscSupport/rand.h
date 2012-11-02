// 
// Random number interface.
// 
// Written by: Jeff Bilmes
//             bilmes@icsi.berkeley.edu
//
// $Header$


#ifndef rand_h
#define rand_h

#include <time.h>
#include <stdio.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <math.h>
#include <sys/types.h>
#include <stdlib.h>
#include <math.h>

#include "general.h"

extern "C" {
  extern double drand48(void) __THROW;
  unsigned short *seed48(unsigned short seed16v[3]) __THROW;
}

class RAND {
 public:
  static unsigned short seedv[3];

  double inverse_error_func(double p);
  double inverse_normal_func(double p);

  RAND(bool seed = false) { 
    if (seed) {
      seedv[0] = time(0); seedv[1] = time(0); seedv[2] = time(0);
      seed48(seedv); 
    }
  }

  void seed()
    {
      seedv[0] = time(0); seedv[1] = time(0); seedv[2] = time(0);
      seed48(seedv); 
    }

  // seed with a double 
  void seed(double* d) {
    unsigned short* tmp = (unsigned short*)d;
    seed48(tmp);
  }
  // seed with a float
  void seed(float* f) {
    double d = *f;
    seed(&d);
  }


  /* return an integer random number uniformly in [l,u] inclusive, l < u */
  int uniform(int l, int u) 
    { return (l + (int)((1+u-l)*::drand48())); }

  /* return an integer random number uniformly in [0,u] inclusive */
  int uniform(int u) 
    { return ((int)((u+1)* ::drand48())); }

  /* return an integer random number uniformly in [0,u-1] inclusive */
  int uniformOpen(int u) 
    { return ((int)((u)* ::drand48())); }

  /* return 1 with probability p and 0 with probability 1-p */
  int coin(float p) 
    { return (int)
	(drand48() < p); }

  // drand48 plus epsilon, i.e., (0,1)
  double drand48pe() {
    double rc;
    do {
      rc = ::drand48();
    } while (rc == 0.0);
    return rc;
  }

  // simple interface to drand
  double drand48() { return ::drand48(); }

  // choose integers from [0:(u-1)] according to a length u 'dist', 
  // assume sum(dist) = 1.0
  int sample(const int u,const float *const dist);

  // return a sample from a N(0,1) zero mean unity variance
  // Gaussian
  double normal() {
    return inverse_normal_func(drand48pe());
  }

  // randomly permute the integers in vec of lenght len
  void rpermute(int * vec, const unsigned len);

  // randomly permute the unsigned integers in vec of length len
  void rpermute(unsigned * vec, const unsigned len);

};


extern RAND rnd;
#endif
