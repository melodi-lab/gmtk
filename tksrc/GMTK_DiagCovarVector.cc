/*-
 * GMTK_DiagCovarVector.cc
 *     Trainable shared diagonal covariance matrix (vector)
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
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <ieeefp.h>
#include <float.h>
#include <assert.h>

#include "general.h"
#include "error.h"
#include "rand.h"

#include "GMTK_DiagCovarVector.h"
#include "GMTK_GaussianComponent.h"
#include "GMTK_GMParms.h"

#ifndef M_PI
#define M_PI               3.14159265358979323846  /* pi */
#endif



VCID("$Header$");


//////////////////////////////////
// static member initialization
unsigned DiagCovarVector::numFlooredVariances = 0;

////////////////////////////////////////////////////////////////////
//        General create, read, destroy routines 
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * DiagCovarVector::DiagCovarVector()
 *      Constructor
 *
 * Results:
 *      Constructs the object.
 *
 * Side Effects:
 *      None so far.
 *
 *-----------------------------------------------------------------------
 */
DiagCovarVector::DiagCovarVector() 
{
}


/*-
 *-----------------------------------------------------------------------
 * DiagCovarVector::read()
 *      read in, and make sure it is a valid cov. matrix.
 *
 * Results:
 *      Constructs the object.
 *
 * Side Effects:
 *      None so far.
 *
 *-----------------------------------------------------------------------
 */
void 
DiagCovarVector::read(iDataStreamFile& is)
{
  NamedObject::read(is);
  covariances.read(is); 
  for (int i=0;i<covariances.len();i++) {
    if (covariances[i] < GaussianComponent::varianceFloor()) {
      error("Error: reading diagonal covariance matrix '%s' (from file '%s'), but covariance[%d] = (%e) < current Floor = (%e)",
	    name().c_str(),is.fileName(),
	    i,covariances[i],GaussianComponent::varianceFloor());
    }
  }
  preCompute();
}



////////////////////////////////////////////////////////////////////
//        Misc Support
////////////////////////////////////////////////////////////////////

void
DiagCovarVector::makeRandom()
{
  for (int i=0;i<covariances.len();i++) {
    covariances[i] = 1.0+rnd.drand48pe();
  }
  preCompute();
}

void
DiagCovarVector::makeUniform()
{
  for (int i=0;i<covariances.len();i++) {
    covariances[i] = 1.0;
  }
  preCompute();
}


void
DiagCovarVector::preCompute()
{
  variances_inv.growIfNeeded(covariances.len());
  double det = 1.0;
  for (int i=0;i<covariances.len();i++) {
    if (covariances[i] <= GaussianComponent::varianceFloor()) {
      // Theoretically, this shouldn't happen unless you are reading
      // in a file that was computed from a previous run with a different threshold.
      coredump("ERROR: element %i of diagonal covariance matrix '%s' is at or below floor value",i,name().c_str());
    }
    variances_inv[i] = 1.0/covariances[i];
    det *= covariances[i];
  }
  if (det <= DBL_MIN) {
    coredump("ERROR: determinant of diagonal covariance matrix '%s' has hit minimum",name().c_str());
  }
  const double tmp = (::pow(2*M_PI,covariances.len()/2.0)*::sqrt(det));
  if (tmp <= DBL_MIN)
    coredump("ERROR:  norm const has hit maximum of diagonal covariance matrix '%s'",name().c_str());
  _log_inv_normConst = -0.5*(covariances.len()*::log(2*M_PI) + ::log(det));
}






/////////////////
// EM routines //
/////////////////



void
DiagCovarVector::emStartIteration()
{
  if (!GM_Parms.amTrainingCovars())
    return;

  if(emOnGoingBitIsSet())
    return; // already done

  if (!emEmAllocatedBitIsSet()) {
    // this is presumably the first time
    nextCovariances.resize(covariances.len());
    nextMeans.resize(covariances.len());
    emSetEmAllocatedBit();
  }
  // EM iteration is now going.
  emSetOnGoingBit();
  emSetSwappableBit();
  
  accumulatedProbability = 0.0;
  numFlooredVariances = 0;
  for (int i=0;i<covariances.len();i++) {
    nextCovariances[i] = nextMeans[i] = 0.0;
  }
}


void
DiagCovarVector::emIncrement(logpr prob,
			     const float*f,
			     const Data32* const base,
			     const int stride)
{
  if (!GM_Parms.amTrainingCovars())
    return;

  emStartIteration();

  if (prob < minIncrementProbabilty) {
    missedIncrementCount++;
    return;
    // don't accumulate anything since this one is so small and
    // if we did an unlog() and converted to a single precision
    // floating point number, it would be a denomral or underflow.
  } 

  accumulatedProbability += prob;
  const float f_prob = prob.unlog();

  float * means_p = nextMeans.ptr;
  float * covars_p = nextCovariances.ptr;
  float * means_end_p = nextMeans.ptr + nextMeans.len();
  const float * f_p = f;
  do {
    const float tmp = (*f_p)*f_prob;
    *covars_p += (*f_p)*tmp;
    *means_p += tmp;

    covars_p++;
    means_p++;
    f_p++;
  } while (means_p != means_end_p);

}


void
DiagCovarVector::emEndIteration()
{
  if (!GM_Parms.amTrainingCovars())
    return;

  if (!emOnGoingBitIsSet())
    return;

  if (accumulatedProbability.zero()) {
    // TODO: need to check if this will overflow here
    // when dividing by it. This is more than just checking
    // for zero. Also need to do this in every such EM object.
    warning("WARNING: Diagonal covariance vector named '%s' did not receive any accumulated probability in EM iteration",name().c_str());
  }

  const double invRealAccumulatedProbability = 
    accumulatedProbability.inverse().unlog();
  // finish computing the next means.

  unsigned prevNumFlooredVariances = numFlooredVariances;
  for (int i=0;i<nextMeans.len();i++) {
    nextMeans[i] *= invRealAccumulatedProbability;
    nextCovariances[i] *= invRealAccumulatedProbability;
    
    nextCovariances[i] = 
      nextCovariances[i]  - nextMeans[i]*nextMeans[i];

    
    /////////////////////////////////////////////////////
    // When variances hit zero or their floor:
    // There could be several reasons for the variances hitting the floor:
    //   1) The prediction of means is very good which leads to low
    //      variances.  In this case, we shouldn't drop the component,
    //      instead we should just hard-limit the variance (i.e., here
    //      mixCoeffs[this] is not too small).  
    //   2) Very small quantity of training data (i.e., mixCoeffs[this] is
    //      very small). In this case we should remove the component
    //      completely.
    //   3) If there is only one mixture, and this happens, then it could be
    //      that the prob of this Gaussian is small. In this case, there's 
    //      probably a problem with the graph topology. 
    //      I.e., really, we should remove the RV state leading to
    //      this this Gaussian. For now,
    //      however, if this happens, the variance will be floored like
    //      in case 1.

    if (nextCovariances[i] < GaussianComponent::varianceFloor()) {

      numFlooredVariances++;

      // Don't let variances go less than variance floor. 

      // At this point, we either could keep the old variance 
      // values, or hard limit them to the varianceFloor. 

      // either A or B but not both should be uncommented below.

      // A: keep old variance
      nextCovariances[i] = covariances[i];

      // B: hard limit the variances
      // nextCovariances[i] = GaussianComponent::varianceFloor();
    }
  }
  if (prevNumFlooredVariances < numFlooredVariances) {
    warning("WARNING: covariance vector named '%s' had %d variances floored\n",
	    name().c_str(),numFlooredVariances-prevNumFlooredVariances);
  }
  
  // stop EM
  emClearOnGoingBit();
}


void
DiagCovarVector::emSwapCurAndNew()
{
  if (!GM_Parms.amTrainingCovars())
    return;

  if (!emSwappableBitIsSet())
    return;

  for (int i=0;i<covariances.len();i++) {
    genSwap(covariances[i],nextCovariances[i]);
  }
  // set up new parameters for their potential next use.
  preCompute();

  emClearSwappableBit();
}


void
DiagCovarVector::emStoreAccumulators(oDataStreamFile& ofile)
{
  error("not implemented");
}

void
DiagCovarVector::emLoadAccumulators(iDataStreamFile& ifile)
{
  error("not implemented");
}


void
DiagCovarVector::emAccumulateAccumulators(iDataStreamFile& ifile)
{
  error("not implemented");
}






