/*-
 * 
 * GMTK_LinMeanCondDiagGaussianOpt.cc
 * 
 *        Code from LinMeanCondDiagGaussian that will benefit from other optimizations.
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
#include "lineqsolve.h"

#include "GMTK_LinMeanCondDiagGaussian.h"
#include "GMTK_GMParms.h"
#include "GMTK_MixtureCommon.h"
#include "GMTK_DiagGaussian.h"


/*-
 *-----------------------------------------------------------------------
 * log_p()
 *      Computes the probability of this Gaussian.
 * 
 * Preconditions:
 *      preCompute() must have been called before this.
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
LinMeanCondDiagGaussian::log_p(const float *const x,
			       const Data32* const base,
			       const int stride)
{
  assert ( basicAllocatedBitIsSet() );

  logpr rc;
  rc.set_to_zero();
  Dlinks* const dLinks = dLinkMat->dLinks;


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
  assert ( dLinks->preComputedOffsets.len() == dLinkMat->arr.len() );
  const int* lagStrideOffsetsp = dLinks->preComputedOffsets.ptr;
  const float* buryValsp = dLinkMat->arr.ptr;
  DIAG_GAUSSIAN_TMP_ACCUMULATOR_TYPE d=0.0;

  int i=0; do {
    // TODO: vectorize this version for better unrolling 
    float u=0.0;
    const int nLinks = dLinks->numLinks(i);
    if (nLinks > 0) {
      // TODO: this is just a dot-product, should call inlined vectorized version of this.
      const int *lagStrideOffsets_endp = lagStrideOffsetsp+nLinks;
      do {
	u += (*buryValsp) *
	  *((float*)base + *lagStrideOffsetsp);
	lagStrideOffsetsp++;
	buryValsp++;
      } while (lagStrideOffsetsp != lagStrideOffsets_endp);
    }
    u += *mean_p;

    const DIAG_GAUSSIAN_TMP_ACCUMULATOR_TYPE tmp
      = (*xp - u);

    d += tmp*tmp*(*var_inv_p);

    xp++;
    mean_p++;
    var_inv_p++;
    i++;
  } while (xp != x_endp);
  d *= -0.5;
  return logpr(0,(covar->log_inv_normConst() + d));

}



/////////////////
// EM routines //
/////////////////


void
LinMeanCondDiagGaussian::fkIncrementMeanDiagCovarDlinks()
{
  // TODO: go through and implmenet the fk for this class. Will need
  // to go into the dlink structure, as basically we need to compute:
  //      \sum_t p_t(i,l|x_{1:T}) D_{il} (x_t - B_{il}z_t - u_{il})z_t^T
  // See the TR for more information.
  error("LinMeanCondDiagGaussian::fkIncrementMeanDiagCovarDlinks(): not implemented");
}



