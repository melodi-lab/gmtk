//
//  Random Generation support.
//
//
// Written by: Jeff Bilmes
//             bilmes@icsi.berkeley.edu


#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#ifndef __CYGWIN__
#include <values.h>
#endif
#include <math.h>
#include <float.h>
#include <assert.h>

#include "general.h"
VCID("$Header$")
#include "error.h"

#include "rand.h"

#ifndef M_SQRT2
#define M_SQRT2        1.41421356237309504880  /* sqrt(2) */
#endif

// declare a default global rnd object.
// RAND rnd(false); // change to true at some point.
//////////////////////////////
// data needed for rand objects.
unsigned short RAND::seedv[3];


/*-
 *-----------------------------------------------------------------------
 * RAND::sample
 *      Sample an integer in the range [0:u] according to the
 *      probability mass function given by 'dist'
 *
 * Results:
 *      Returns the resulting integer.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
int
RAND::sample(const int u, const float *const dist)
{

  int i;
  double sum = 0.0;
  for (i=0;i<u;i++) {
    sum += dist[i];
  }
  double tmp = ::drand48()*sum;
  const float *arr_p = dist;
  i = 0;
  // do a dumb linear search for now.
  do {
    tmp -= *arr_p++;
    i++;
  } while (tmp > 0.0 && i < u);
  return i-1;
}


/*-
 *-----------------------------------------------------------------------
 * inverse_error_func
 *      This function implements an approximation to the inverse
 *      error function. The code was taken from the source code
 *      of gnuplot. See Abramovitz & Stegan (or some stat book)
 *      for the definition of the inverse error function.
 *
 * Results:
 *      returns the inverse error function for value p
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
double 
RAND::inverse_error_func(double p) 
{
	/* 
           Source: This routine was derived (using f2c) from the 
                   FORTRAN subroutine MERFI found in 
                   ACM Algorithm 602 obtained from netlib.

                   MDNRIS code contains the 1978 Copyright 
                   by IMSL, INC. .  Since MERFI has been 
                   submitted to netlib, it may be used with 
                   the restriction that it may only be 
                   used for noncommercial purposes and that
                   IMSL be acknowledged as the copyright-holder
                   of the code.
        */



	/* Initialized data */
	static double a1 = -.5751703;
	static double a2 = -1.896513;
	static double a3 = -.05496261;
	static double b0 = -.113773;
	static double b1 = -3.293474;
	static double b2 = -2.374996;
	static double b3 = -1.187515;
	static double c0 = -.1146666;
	static double c1 = -.1314774;
	static double c2 = -.2368201;
	static double c3 = .05073975;
	static double d0 = -44.27977;
	static double d1 = 21.98546;
	static double d2 = -7.586103;
	static double e0 = -.05668422;
	static double e1 = .3937021;
	static double e2 = -.3166501;
	static double e3 = .06208963;
	static double f0 = -6.266786;
	static double f1 = 4.666263;
	static double f2 = -2.962883;
	static double g0 = 1.851159e-4;
	static double g1 = -.002028152;
	static double g2 = -.1498384;
	static double g3 = .01078639;
	static double h0 = .09952975;
	static double h1 = .5211733;
	static double h2 = -.06888301;

	/* Local variables */
	static double a, b, f, w, x, y, z, sigma, z2, sd, wi, sn;

	x = p;

	/* determine sign of x */
	if (x > 0)
		sigma = 1.0;
	else
		sigma = -1.0;

	/* Note: -1.0 < x < 1.0 */

	z = fabs(x);

	/* z between 0.0 and 0.85, approx. f by a 
	   rational function in z  */

	if (z <= 0.85) {
		z2 = z * z;
		f = z + z * (b0 + a1 * z2 / (b1 + z2 + a2 
		    / (b2 + z2 + a3 / (b3 + z2))));

	/* z greater than 0.85 */
	} else {
		a = 1.0 - z;
		b = z;

		/* reduced argument is in (0.85,1.0), 
		   obtain the transformed variable */

		w = sqrt(-(double)log(a + a * b));

		/* w greater than 4.0, approx. f by a 
		   rational function in 1.0 / w */

		if (w >= 4.0) {
			wi = 1.0 / w;
			sn = ((g3 * wi + g2) * wi + g1) * wi;
			sd = ((wi + h2) * wi + h1) * wi + h0;
			f = w + w * (g0 + sn / sd);

		/* w between 2.5 and 4.0, approx. 
		   f by a rational function in w */

		} else if (w < 4.0 && w > 2.5) {
			sn = ((e3 * w + e2) * w + e1) * w;
			sd = ((w + f2) * w + f1) * w + f0;
			f = w + w * (e0 + sn / sd);

		/* w between 1.13222 and 2.5, approx. f by 
		   a rational function in w */
		} else if (w <= 2.5 && w > 1.13222) {
			sn = ((c3 * w + c2) * w + c1) * w;
			sd = ((w + d2) * w + d1) * w + d0;
			f = w + w * (c0 + sn / sd);
		}
	}
	y = sigma * f;
	return(y);
}




/*-
 *-----------------------------------------------------------------------
 * inverse_normal_func
 *      This function implements an approximation to the inverse
 *      normal(0,1) function. The code was taken from the source code
 *      of gnuplot. See Abramovitz & Stegan (or some stat book)
 *      for the definition of the inverse error function.
 *      So, to get a sample from a Gaussian N(0,1), call
 *      inverse_normal_func(rnd.drand48pe()).
 *
 * Results:
 *      returns the inverse Normal function for value p, where
 *      'p' is presummably uniformly distributed on (0,1).
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
double 
RAND::inverse_normal_func(double p)
{
	/* 
           Source: This routine was derived (using f2c) from the 
                   FORTRAN subroutine MDNRIS found in 
                   ACM Algorithm 602 obtained from netlib.

                   MDNRIS code contains the 1978 Copyright 
                   by IMSL, INC. .  Since MDNRIS has been 
                   submitted to netlib it may be used with 
                   the restriction that it may only be 
                   used for noncommercial purposes and that
                   IMSL be acknowledged as the copyright-holder
                   of the code.
        */

	/* Initialized data */
	static double eps = 1e-10;
	static double g0 = 1.851159e-4;
	static double g1 = -.002028152;
	static double g2 = -.1498384;
	static double g3 = .01078639;
	static double h0 = .09952975;
	static double h1 = .5211733;
	static double h2 = -.06888301;
	static double sqrt2 = M_SQRT2; /* 1.414213562373095; */

	/* Local variables */
	static double a, w, x;
	static double sd, wi, sn, y;

	// double inverse_error_func(double p);

	/* Note: 0.0 < p < 1.0 */
	/* assert ( 0.0 < p && p < 1.0 ); */

	/* p too small, compute y directly */
	if (p <= eps) {
		a = p + p;
		w = sqrt(-(double)log(a + (a - a * a)));

		/* use a rational function in 1.0 / w */
		wi = 1.0 / w;
		sn = ((g3 * wi + g2) * wi + g1) * wi;
		sd = ((wi + h2) * wi + h1) * wi + h0;
		y = w + w * (g0 + sn / sd);
		y = -y * sqrt2;
	} else {
		x = 1.0 - (p + p);
		y = inverse_error_func(x);
		y = -sqrt2 * y;
	}
	return(y);
} 





/*-
 *-----------------------------------------------------------------------
 * permute
 *      Randomly and uniformly permute the integers in vec of lenght len.
 * 
 *  Preconditions:
 *     vec must be of length len.
 *
 * Results:
 *      returns the permuted vector in place
 *
 * Side Effects:
 *      changes the argument
 *
 *-----------------------------------------------------------------------
 */
void
RAND::rpermute(int * vec, const unsigned len)
{
  const unsigned lenm1 = len-1;
  for (unsigned i=0;i<lenm1;i++) {
    unsigned val = uniform(i,lenm1);
    swap(vec[i],vec[val]);
  }
}
void
RAND::rpermute(unsigned * vec, const unsigned len)
{
  if (len == 0) return;
  const unsigned lenm1 = len-1;
  for (unsigned i=0;i<lenm1;i++) {
    unsigned val = uniform(i,lenm1);
    swap(vec[i],vec[val]);
  }
}



#ifdef MAIN

#include <stdio.h>

int
main()
{
  RAND rnd(true);

  printf("Uniform int between %d and %d is %d\n",
	 0,10,rnd.uniform(0,10));

  printf("Coin flip with p = %f is %d\n",
	 0.6,rnd.coin(0.6));

  printf("normal drand48pe = %f, uniform on (0,1)\n",
	 rnd.drand48pe());

  float dist[6] = { 0.1, 0.1, 0.1, 0.2, 0.3, 0.2 }; 
  
  printf("sample from distribution is %d\n",
	 rnd.sample(5,&dist[0]));
}


#endif
