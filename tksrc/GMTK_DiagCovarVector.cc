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
#include "GMTK_MixGaussiansCommon.h"

#ifndef M_PI
#define M_PI               3.14159265358979323846  /* pi */
#endif

VCID("$Header$");

////////////////////////////////////////////////////////////////////
//        Static members
////////////////////////////////////////////////////////////////////

unsigned DiagCovarVector::numFlooredVariances = 0;

double DiagCovarVector::cloneSTDfrac = 0.1;

void DiagCovarVector::checkForValidValues()
{
  if (DiagCovarVector::cloneSTDfrac < 0)
    error("ERROR: DiagCovarVector's cloneSTDfrac (%e) must be >= 0",
	  DiagCovarVector::cloneSTDfrac);
}


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
  setBasicAllocatedBit();
  preCompute();
}



////////////////////////////////////////////////////////////////////
//        Misc Support
////////////////////////////////////////////////////////////////////

void
DiagCovarVector::makeRandom()
{
  for (int i=0;i<covariances.len();i++) {
    covariances[i] = 10*(1.0+rnd.drand48pe());
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


/*-
 *-----------------------------------------------------------------------
 * noisyClone()
 *      Create a copy of self, but perturb the mean vector
 *      a bit with some noise.
 * 
 * Preconditions:
 *      The mean must be read in, and basicAllocatedBitIsSet() must be true.
 *
 * Postconditions:
 *      none.
 *
 * Side Effects:
 *      'this' is not changed at all. Allocates new memory though.
 *
 * Results:
 *      returns the new noisy mean.
 *
 *-----------------------------------------------------------------------
 */
DiagCovarVector*
DiagCovarVector::noisyClone()
{
  assert ( basicAllocatedBitIsSet() );

  // first check if self is already cloned, and if so, return that.
  DiagCovarVector* clone;

  map<DiagCovarVector*,DiagCovarVector*>::iterator it = MixGaussiansCommon::diagCovarCloneMap.find(this);
  if (it == MixGaussiansCommon::diagCovarCloneMap.end()) {
    clone = new DiagCovarVector();
    unsigned cloneNo=0; do {
      char buff[256];
      sprintf(buff,"%d",cloneNo);
      clone->_name = _name + string("_cl") + buff;
      cloneNo++;
    } while (GM_Parms.covarsMap.find(clone->_name) != GM_Parms.covarsMap.end());
    clone->refCount = 0;
    clone->covariances.resize(covariances.len());
    for (int i=0;i<covariances.len();i++) {
      clone->covariances[i] = covariances[i] + 
	cloneSTDfrac*covariances[i]*rnd.normal();
    }
    clone->setBasicAllocatedBit();
    MixGaussiansCommon::diagCovarCloneMap[this] = clone;

    // also add self to GMParms object.
    GM_Parms.add(clone);
    clone->preCompute();
  } else {
    clone = (*it).second;
  }
  return clone;
}




/*-
 *-----------------------------------------------------------------------
 * precompute()
 *      Precompute a number of internal variables for speed.
 *      This routine **** MUST BE CALLED ANYTIME *** the paramters
 *      of this object change. If this routine is not called,
 *      any Gaussian using this covariance will produce
 *      invalid results.
 * 
 * Preconditions:
 *      basicAllocatedBitIsSet() must be set
 *
 * Postconditions:
 *      All internal variables are allocated, and the 
 *      object is ready for computing probabilities.
 *
 * Side Effects:
 *      nil
 *
 * Results:
 *      nil
 *
 *-----------------------------------------------------------------------
 */
void
DiagCovarVector::preCompute()
{
  assert ( basicAllocatedBitIsSet() );

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
DiagCovarVector::emStartIteration(sArray<float>& componentsNextCovars)
{
  assert ( basicAllocatedBitIsSet() );

  if (!GM_Parms.amTrainingCovars())
    return;

  /////////////////////////////////////////////
  // make sure our caller has its covar accumulator resized
  // and initialized.
  componentsNextCovars.growIfNeeded(covariances.len());
  for (int i=0;i<covariances.len();i++) {
    componentsNextCovars[i] = 0.0;
  }

  if(emOnGoingBitIsSet()) {
    // EM already on going.
    // Increment the count of number of Gaussian Components using this mean.
    refCount++; 
    return; 
  }

  if (!emEmAllocatedBitIsSet()) {
    // this is presumably the first time
    emSetEmAllocatedBit();
    // allocate the final covars means if needed
    nextCovariances.growIfNeeded(covariances.len());
    // allocate the final next means if needed
    // nextMeans.growIfNeeded(covariances.len());
  }
  // EM iteration is now going.
  emSetOnGoingBit();
  
  accumulatedProbability = 0.0;
  numFlooredVariances = 0;
  refCount=1;
  for (int i=0;i<covariances.len();i++) {
    // nextMeans[i] = 0.0;
    nextCovariances[i] = 0.0;
  }


  // make it swapable, although at this point
  // it would swap in the unaccumulated values.
  emSetSwappableBit();
}




/*-
 *-----------------------------------------------------------------------
 * emIncrement
 *      Add the data item for the current f into the
 *      partial accumulator (2nd moment accumulators) object given by the argument
 *      pointer. 
 *      NOTE: This routine lives in the inner most loop of
 *      EM training, so it is important that this is as
 *      fast as possible.
 * 
 * Preconditions:
 *      basic structures must be allocated.
 *
 * Postconditions:
 *      data has been accumulated
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      nil
 *
 *-----------------------------------------------------------------------
 */

void
DiagCovarVector::emIncrement(const logpr prob,
			     const float fprob,
			     const float*f,
			     const Data32* const base,
			     const int stride,
			     float *const partialAccumulatedNextCovar)
{
  assert ( basicAllocatedBitIsSet() );
  if (!GM_Parms.amTrainingCovars())
    return;


  
  /////////////////////////////////////////////
  // Note: unlike the normal EM mode described
  // in GMTK_EMable.h, we do not call
  // emStartIteration() here and assume that it
  // was called by the Gaussian component that
  // is using this mean. This is because
  // this object keeps a reference count (needed for
  // sharing), and calling that routine repeatedly 
  // would result in an incorrect count. We do
  // make sure that em has been allocated with the
  // following assertion.
  assert ( emEmAllocatedBitIsSet() ); 

  // we assume here that (prob > minIncrementProbabilty),
  // i.e., that this condition has been checked by the caller
  // of this routine (meaning that fprob is valid)
  assert (prob >= minIncrementProbabilty);

  accumulatedProbability += prob;

  float * covars_p = partialAccumulatedNextCovar;
  float * covars_end_p = partialAccumulatedNextCovar + covariances.len();
  const float * f_p = f;
  do {
    const float tmp = (*f_p)*(*f_p)*fprob;
    *covars_p += tmp;

    covars_p++;
    f_p++;
  } while (covars_p != covars_end_p);

}



/*-
 *-----------------------------------------------------------------------
 * emEndIteration
 *      End the current em iteration. The way this works is as follows.
 *      Some number of mixture Gaussian (MG) objects might share this covariance
 *      object. Each time it ends, it (the MG object) calls this object with its 
 *      accumulated covariance (but without having been normalized by the sum
 *      of posteriors, as that is done here). That covariance is then
 *      accumulated, and a reference count (keeping track of how many
 *      MGs are sharing this object) is decremented. If the ref count
 *      hits zero, then we do the final update of the covariance.
 *
 *      NOTE: this version of the routine uses the parents' accumulated first and 
 *      2nd order statistics.
 *
 *      NOTE: any changes here should also be made to other emEndIteration routines 
 *      in this object.
 * 
 * Preconditions:
 *      basic structures must be allocated, EM must be ongoing.
 *
 * Postconditions:
 *      em iteration is ended.
 *
 * Side Effects:
 *      possibly updates all next parameters
 *
 * Results:
 *      nil
 *
 *-----------------------------------------------------------------------
 */
void
DiagCovarVector::emEndIteration(const logpr parentsAccumulatedProbability,
				const float*const partialAccumulatedNextMeans,
				const float *const partialAccumulatedNextCovar)
{
  assert ( basicAllocatedBitIsSet() );
  if (!GM_Parms.amTrainingCovars())
    return;

  if (refCount > 0) {
    // if this isn't the case, something is wrong.
    assert ( emOnGoingBitIsSet() );


    // we just accumulate in the covariance
    // shared by the parent.

    // NOTE: What are we doing here? We need to compute
    // the weighted average of the shared means and covariances.
    // Lets say that C is the shared covariance formula,
    // for all classes in the set J. Cj is the covariance
    // just for class j. The normal update formula for class j's
    // covariance is the following:
    // 
    //        \sum p(j|x_i) (x_i-m_j)^2
    //          i
    // Cj =  ----------------------------------
    //        \sum p(j|x_i)
    //          i
    //
    // Now, however, we are equating Cj with a number of classes,
    // i.e., we are saying that 
    // 
    //     C = \sum p(j) Cj
    //         j in J
    //
    // where p(j) = \sum_i p(j|x_i) / \sum_k \sum_i p(k|x_i)
    // 
    // and we are using for all j in J C in place of Cj
    // This then gives the formula
    // 
    //                              \sum p(j|x_i) (x_i-m_j)^2
    //        1.0                     i
    // C =  ------- \sum E[I(j)] --------------------------------------
    //         N    j in J            \sum p(j|x_i)
    //                                  i
    // 
    // where N = \sum E[I(j)]
    //             j
    // 
    // and where I(j) is the indicator of class j so 
    // E[I(j)] = \sum_i p(j|x_i) are the expected counts, and
    // where p(j) = E[I(j)]/N
    // 
    // The complication is that the means and covariances might be
    // shared.  What we are doing here is, using the EM version of the
    // formula cov(X) = EX^2 - (EX)^2 accumulating each of the class
    // conditional partially weighted EX^2 and (EX)^2 values.  This is
    // valid since in the above formula E[I(j)] cancels out the
    // denominator when we are computing the final C matrix. The
    // entire thing needs to be divided again by N, though, which is
    // done below after refCount hits zero below.

    if ( parentsAccumulatedProbability >
	 GaussianComponent::minAccumulatedProbability()) {
      // Only accumulate here if there is something significant 
      // to accumlate.

      // accumulate the 1st and 2nd order statistics given
      // by the mean object.
      const double invRealAccumulatedProbability = 
	parentsAccumulatedProbability.inverse().unlog();
    
      for (int i=0;i<covariances.len();i++) {
	nextCovariances[i] += 
	  (partialAccumulatedNextCovar[i] - 
	   partialAccumulatedNextMeans[i]*partialAccumulatedNextMeans[i]*invRealAccumulatedProbability);
      }
    }

    refCount--;
  }


  /////////////////////////////////////////////
  // if there is still someone who
  // has not given us his/her 1st order stats,
  // then we return w/o finishing.
  if (refCount > 0)
    return;

  // otherwise, we're ready to finish and
  // compute the next covariances.

  if (accumulatedProbability < GaussianComponent::minAccumulatedProbability()) {
    warning("WARNING: Diag covariance vec '%s' received only %e accumulated log probability in EM iteration, using previous values.",
	    name().c_str(),
	    accumulatedProbability.val());
    for (int i=0;i<covariances.len();i++) 
      nextCovariances[i] = covariances[i];
  } else {

    // TODO: should check for possible overflow here of 
    // accumulatedProbability when we do the inverse and unlog.
    // Ideally this won't happen for a given minAccumulatedProbability().
    const double invRealAccumulatedProbability = 
      accumulatedProbability.inverse().unlog();

    // Finally, divide by N (see the equation above)
    // here computing the final variances.
    unsigned prevNumFlooredVariances = numFlooredVariances;
    double minVar = 0.0;
    for (int i=0;i<covariances.len();i++) {
      nextCovariances[i] *= invRealAccumulatedProbability;


      if (i == 0 || nextCovariances[i] < minVar)
	minVar = nextCovariances[i];
    
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
      warning("WARNING: covariance vector named '%s' had %d variances floored, minimum variance found was %e.\n",
	      name().c_str(),numFlooredVariances-prevNumFlooredVariances,minVar);
    }

  }

  // stop EM
  emClearOnGoingBit();
}





/*-
 *-----------------------------------------------------------------------
 * emEndIteration
 *      End the current em iteration. The way this works is as follows.
 *      Some number of mixture Gaussian (MG) objects might share this covariance
 *      object. Each time it ends, it (the MG object) calls this object with its 
 *      accumulated covariance (but without having been normalized by the sum
 *      of posteriors, as that is done here). That covariance is then
 *      accumulated, and a reference count (keeping track of how many
 *      MGs are sharing this object) is decremented. If the ref count
 *      hits zero, then we do the final update of the covariance.
 *
 *      NOTE: this version of the routine uses the already mostly computed 
 *      parents' next covariance.
 *
 *      NOTE: any changes here should also be made to other emEndIteration routines 
 *      in this object.
 * 
 * Preconditions:
 *      basic structures must be allocated, EM must be ongoing.
 *
 * Postconditions:
 *      em iteration is ended.
 *
 * Side Effects:
 *      possibly updates all next parameters
 *
 * Results:
 *      nil
 *
 *-----------------------------------------------------------------------
 */
void
DiagCovarVector::emEndIteration(const float *const parentsAccumulatedNextCovar)
{
  assert ( basicAllocatedBitIsSet() );
  if (!GM_Parms.amTrainingCovars())
    return;

  if (refCount > 0) {
    // if this isn't the case, something is wrong.
    assert ( emOnGoingBitIsSet() );

    // we just accumulate in the covariance
    // shared by the parent.
    for (int i=0;i<covariances.len();i++)
      nextCovariances[i] += parentsAccumulatedNextCovar[i];

    refCount--;
  }


  /////////////////////////////////////////////
  // if there is still someone who
  // has not given us his/her 1st order stats,
  // then we return w/o finishing.
  if (refCount > 0)
    return;

  // otherwise, we're ready to finish and
  // compute the next covariances.

  if (accumulatedProbability < GaussianComponent::minAccumulatedProbability()) {
    warning("WARNING: Diag covariance vec '%s' received only %e accumulated log probability in EM iteration, using previous values.",accumulatedProbability.val(),name().c_str());
    for (int i=0;i<covariances.len();i++) 
      nextCovariances[i] = covariances[i];
  } else {

    // TODO: should check for possible overflow here of 
    // accumulatedProbability when we do the inverse and unlog.
    // Ideally this won't happen for a given minAccumulatedProbability().
    const double invRealAccumulatedProbability = 
      accumulatedProbability.inverse().unlog();

    // Finally, divide by N (see the equation above)
    // here computing the final variances.
    unsigned prevNumFlooredVariances = numFlooredVariances;
    double minVar = 0.0;
    for (int i=0;i<covariances.len();i++) {
      nextCovariances[i] *= invRealAccumulatedProbability;

      if (i == 0 || nextCovariances[i] < minVar)
	minVar = nextCovariances[i];

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
      warning("WARNING: covariance vector named '%s' had %d variances floored, minimum variance found was %e.\n",
	      name().c_str(),numFlooredVariances-prevNumFlooredVariances,minVar);
    }
  }

  // stop EM
  emClearOnGoingBit();
}




void
DiagCovarVector::emSwapCurAndNew()
{
  assert ( basicAllocatedBitIsSet() );
  if (!GM_Parms.amTrainingCovars())
    return;

  // we should have that the number of calls
  // to emStartIteration and emEndIteration are
  // the same.
  assert ( refCount == 0 );

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
  assert ( basicAllocatedBitIsSet() );
  assert ( emEmAllocatedBitIsSet() );
  EMable::emStoreAccumulators(ofile);
}

void
DiagCovarVector::emLoadAccumulators(iDataStreamFile& ifile)
{
  assert ( basicAllocatedBitIsSet() );
  assert ( emEmAllocatedBitIsSet() );
  EMable::emLoadAccumulators(ifile);
}


void
DiagCovarVector::emAccumulateAccumulators(iDataStreamFile& ifile)
{
  assert ( basicAllocatedBitIsSet() );
  assert ( emEmAllocatedBitIsSet() );
  EMable::emAccumulateAccumulators(ifile);
}






