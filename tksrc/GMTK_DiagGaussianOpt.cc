/*-
 * GMTK_DiagGaussianOpt.cc
 *      -  Code for plain vanilla diagonal Gaussians.
 *      -  These routines might benefit from separate optimization (such as separate loop unrolling).
 *     
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2001, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle, and Jeff Bilmes make no representations about
 * the suitability of this software for any purpose.  It is provided
 * "as is" without express or implied warranty.
 *
 */


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <float.h>
#include <assert.h>
#include <ctype.h>

#include <string>

#include "general.h"
VCID("$Header$")
#include "error.h"
#include "rand.h"

#include "GMTK_DiagGaussian.h"
#include "GMTK_GMParms.h"
#include "GMTK_MixtureCommon.h"



/*-
 *-----------------------------------------------------------------------
 * log_p()
 *      Computes the probability of this Gaussian.
 * 
 * Preconditions:
 *      preCompute() must have been called on covariance matrix before this.
 *
 * Postconditions:
 *      nil
 *
 * Side Effects:
 *      nil, other than possible FPEs if the values are garbage
 *
 * Results:
 *      Returns the probability.
 *
 *-----------------------------------------------------------------------
 */

logpr
DiagGaussian::log_p(const float *const x,
		    const Data32* const base,
		    const int stride)
{
  assert ( basicAllocatedBitIsSet() );

  //////////////////////////////////////////////////////////////////
  // The local accumulator type in this routine.
  // This can be changed from 'float' to 'double' to
  // provide extra range for temporary accumulators. Alternatively,
  // decreasing the program's mixCoeffVanishRatio at the beginning
  // of training should eliminate any component that produces
  // such low scores.
#define DIAG_GAUSSIAN_TMP_ACCUMULATOR_TYPE double

  ////////////////////
  // note: 
  // covariances must have been precomputed for this
  // to work.
  const float *xp = x;
  const float *const x_endp = x + _dim;
  const float *mean_p = mean->basePtr();
  const float *var_inv_p = covar->baseVarInvPtr();
  DIAG_GAUSSIAN_TMP_ACCUMULATOR_TYPE d=0.0;
  do {
    const DIAG_GAUSSIAN_TMP_ACCUMULATOR_TYPE tmp
      = (*xp - *mean_p);
    d += (tmp*(*var_inv_p))*tmp;

    xp++;
    mean_p++;
    var_inv_p++;
  } while (xp != x_endp);
  d *= -0.5;
  return logpr(0,(covar->log_inv_normConst() + d));
}



/////////////////
// EM routines //
/////////////////



/*-
 *-----------------------------------------------------------------------
 * emIncrementMeanDiagCovar
 *      Simultaneously increments a mean and a diagonal covariance vector
 *      with one loop rather than doing each separately with two loops.
 * 
 * Preconditions:
 *      Vectors must be allocated and pointing to appropriately sized
 *      arrays. No other assumptions are made (e.g., such as like prob
 *      is large enough).
 *
 * Postconditions:
 *      Vectors have been accumulated by f.
 *
 * Side Effects:
 *      Changes meanAccumulator and diagCovarAccumulator arrays.
 *
 * Results:
 *      nothing.
 *
 *-----------------------------------------------------------------------
 */
void 
DiagGaussian::emIncrementMeanDiagCovar(const float fprob,
				       const float * const f,
				       const unsigned len,
				       float *meanAccumulator,
				       float *diagCovarAccumulator)
{
  register const float * f_p = f;
  register const float *const f_p_endp = f + len;
  register float *meanAccumulator_p = meanAccumulator;
  register float *diagCovarAccumulator_p = diagCovarAccumulator;
  do {

#if 0
    // this code has aliasing so is commented out in favor of the
    // code below.
    register float tmp = (*f_p)*fprob;
    *meanAccumulator_p += tmp;
    tmp *= (*f_p);
    *diagCovarAccumulator_p += tmp;
    meanAccumulator_p++;
    diagCovarAccumulator_p++;
    f_p ++;
#endif

    // a version of the above code that avoids aliasing of f_p and
    // meanAccumulator_p so the compiler can probably optimize better.

    register float tmp = (*f_p)*fprob;
    register float tmp2 = tmp*(*f_p);

    *meanAccumulator_p += tmp;
    *diagCovarAccumulator_p += tmp2;

    meanAccumulator_p++;
    diagCovarAccumulator_p++;
    f_p ++;

  } while (f_p != f_p_endp);
}




/*-
 *-----------------------------------------------------------------------
 * fkIncrementMeanDiagCovar
 *      Simultaneously increments a mean and a diagonal covariance vector
 *      with one loop rather than doing each separately with two loops.
 *      Here we incement the mean and covariance according to what is needed
 *      to produce the Fisher kernel vector (i.e., what is needed to produce the
 *      Fisher kernel of the DBN). The resulting feature space parameters are
 *      stored in the very same EM accumulators.
 * 
 * Preconditions:
 *      Vectors must be allocated and pointing to appropriately sized
 *      arrays. No other assumptions are made (e.g., such as like prob
 *      is large enough).
 *
 * Postconditions:
 *      Vectors have been accumulated by f.
 *
 * Side Effects:
 *      Changes meanAccumulator and diagCovarAccumulator arrays.
 *
 * Results:
 *      nothing.
 *
 *-----------------------------------------------------------------------
 */
void 
DiagGaussian::fkIncrementMeanDiagCovar(const float fprob,
				       const float * const f,
				       const unsigned len,
				       float *curMeans,
				       float *curDiagCovars,
				       float *meanAccumulator,
				       float *diagCovarAccumulator)
{
  register const float * f_p = f;
  register const float *const f_p_endp = f + len;

  register float *mean_p = curMeans;
  register float *diagCovar_p = curDiagCovars;

  register float *meanAccumulator_p = meanAccumulator;
  register float *diagCovarAccumulator_p = diagCovarAccumulator;

  do {

    register float tmp = (*f_p - *mean_p)/(*diagCovar_p);

    // store the values in temporaries so that the
    // compiler knows that there is no aliasing.
    register float mean_val = fprob*tmp;
    tmp = tmp*tmp;
    register float covar_val = -0.5*fprob*(1.0/(*diagCovar_p) + tmp);

    // increment the actual accumulators
    *meanAccumulator_p += mean_val;
    *diagCovarAccumulator_p += covar_val;

    // update all pointers

    mean_p++;
    diagCovar_p++;
    meanAccumulator_p++;
    diagCovarAccumulator_p++;
    f_p ++;
  } while (f_p != f_p_endp);

}


