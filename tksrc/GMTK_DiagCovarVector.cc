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
#include <float.h>
#include <assert.h>

#include "general.h"
#include "error.h"
#include "rand.h"

#include "GMTK_DiagCovarVector.h"
#include "GMTK_GaussianComponent.h"
#include "GMTK_GMParms.h"
#include "GMTK_MixGaussiansCommon.h"
#include "GMTK_MeanVector.h"
#include "GMTK_DlinkMatrix.h"

#ifndef M_PI
#define M_PI               3.14159265358979323846  /* pi */
#endif

VCID("$Header$");

////////////////////////////////////////////////////////////////////
//        Static members
////////////////////////////////////////////////////////////////////

unsigned DiagCovarVector::numFlooredVariances = 0;

bool DiagCovarVector::floorVariancesWhenReadIn = false;

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
  int length;
  is.read(length,"DiagCovarVector read length");
  if (length <= 0)
    error("ERROR: diag covariance matrix %s specifies length (%d) < 0 in input. Must be positive.",
	  name().c_str(),length);
  covariances.resize(length);
  unsigned numFloored=0;
  for (int i=0;i<length;i++) {
    is.read(covariances[i],":reading covar value");
    if (covariances[i] < (float)GaussianComponent::varianceFloor()) {
      if (!floorVariancesWhenReadIn) {
	error("Error: reading diagonal covariance matrix '%s' (from file '%s'), but covariance[%d] = (%e) < current Floor = (%e)",
	      name().c_str(),
	      is.fileName(),
	      i,covariances[i],
	      GaussianComponent::varianceFloor());
      } else {
	numFloored++;
	covariances[i] = GaussianComponent::varianceFloor();
      }
    }
  }
  if (numFloored > 0)
    warning("WARNING: reading diagonal covariance matrix '%s' (from file '%s'), and %d variance values (out of %d) were  < current Floor = (%e), forcing them to floor.",
	    name().c_str(),
	    is.fileName(),
	    numFloored,
	    covariances.len(),
	    GaussianComponent::varianceFloor());
  setBasicAllocatedBit();
  preCompute();
  numTimesShared = 0;
  refCount = 0;
}


/*-
 *-----------------------------------------------------------------------
 * DiagCovarVector::write()
 *      write out
 *
 * Results:
 *
 * Side Effects:
 *      None so far.
 *
 *-----------------------------------------------------------------------
 */
void 
DiagCovarVector::write(oDataStreamFile& os)
{
  NamedObject::write(os);
  os.write(covariances.len(),"diag cov vector write length");
  for (int i=0;i<covariances.len();i++) {
    os.write(covariances[i],"diag cov vector write, values");
  }
  os.nl();
}


////////////////////////////////////////////////////////////////////
//        Misc Support
////////////////////////////////////////////////////////////////////

void
DiagCovarVector::makeRandom()
{
  if (!emAmTrainingBitIsSet())
    return;

  for (int i=0;i<covariances.len();i++) {
    covariances[i] = 10*(1.0+rnd.drand48pe());
  }
  preCompute();
}

void
DiagCovarVector::makeUniform()
{
  if (!emAmTrainingBitIsSet())
    return;

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
    clone->numTimesShared = 0;
    clone->covariances.resize(covariances.len());
    for (int i=0;i<covariances.len();i++) {

      float tmp;
      // try 10 times to get a covariance above the floor.
      for (int j=0;j<10;j++) {
	tmp = covariances[i] + cloneSTDfrac*covariances[i]*rnd.normal();
	if (tmp >= (float)GaussianComponent::varianceFloor())
	  break;
      }
      if (tmp < (float)GaussianComponent::varianceFloor()) {
	// then we get something guaranteed above the floor
	tmp = GaussianComponent::varianceFloor()+
	  cloneSTDfrac*covariances[i]*fabs(rnd.normal());
      }

      clone->covariances[i] = tmp;
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
    if (covariances[i] < (float)GaussianComponent::varianceFloor()) {
      // Theoretically, this shouldn't happen unless you are reading
      // in a file that was computed from a previous run with a different threshold.
      error("ERROR: element %i of diagonal covariance matrix '%s' has value %.16e which is below floor value of %.16e.",
	    i,
	    name().c_str(),
	    covariances[i],
	    GaussianComponent::varianceFloor());
    }
    variances_inv[i] = 1.0/covariances[i];
    det *= covariances[i];
    if (det <= DBL_MIN) {
      warning("WARNING: determinant of diagonal covariance matrix '%s' (=%e) is hiting minimum (=%e) after %d stages. Possible causes include: 1) not enough training segments, or 2) data that is inappropriately scaled, or 3) too much pruning, or 4) impossible or infrequent state configurations, or 5) not large enough varFloor & floor on read command line args.",
	      name().c_str(),det,DBL_MIN,i);
    }
  }
  if (det <= DBL_MIN) {
    error("ERROR: determinant of diagonal covariance matrix '%s' has hit minimum. Possible causes include: 1) not enough training segments, or 2) data that is inappropriately scaled, or 3) too much pruning, or 4) impossible or infrequent state configurations, or 5) not large enough varFloor & floor on read command line args.",
	  name().c_str());
  }
  const double tmp = (::pow(2*M_PI,covariances.len()/2.0)*::sqrt(det));
  if (tmp <= DBL_MIN)
    coredump("ERROR: norm const has hit maximum of diagonal covariance matrix '%s'",name().c_str());
  _log_inv_normConst = -0.5*(covariances.len()*::log(2*M_PI) + ::log(det));
}






/////////////////
// EM routines //
/////////////////



void
DiagCovarVector::emStartIteration(sArray<float>& componentsNextCovars)
{
  assert ( basicAllocatedBitIsSet() );

  // we return if both 1) the not training bit is set
  // and 2) there is no chance that this object will be shared.
  // If 1) is not true, we are training this object so we continue,
  // and if 2) is not true (we are sharing), then the accumulaters
  // created for this object will be needed by the other objects 
  // in this object's Gaussian component, so we'll need to compute them,
  // even though the parameters of this object will not be updated
  // when we swap them in.
  // if (numTimesShared == 1 && !emAmTrainingBitIsSet())
  // return;

  if(emOnGoingBitIsSet()) {
    // EM already on going.
    // Increment the count of number of Gaussian Components using this mean.
    refCount++; 
    // this object therefore is shared, set the bit saying so.
    emSetSharedBit();

    // Make sure our callers accumulators are allocated.  The reason
    // for this is that the caller of this routine is one who is
    // sharing this object with at least one other caller, and this
    // caller is being set up after the first caller.  This caller has
    // therefore not had its own accumulators allocated yet unless
    // this is the second iteration in an internal EM iteration run
    // (e.g., we are not running in parallel), but in any event it
    // should not be calling its emStartIteration() multiple times.
    componentsNextCovars.growIfNeeded(covariances.len());
    for (int i=0;i<covariances.len();i++) {
      componentsNextCovars[i] = 0.0;
    }

    // We return now since we might have already
    // accumulated some probability for this object
    // (which would be stored in accumulatedProbability)
    // but accumulated it for an object that is
    // sharing self but has a different set of its
    // own accumulators.
    return; 
  }

  if (!emEmAllocatedBitIsSet()) {
    // this is presumably the first time
    emSetEmAllocatedBit();
  }
  // EM iteration is now going.
  emSetOnGoingBit();

  // accumulators are not initialized at this point.
  emClearAccInitializedBit();

  accumulatedProbability = 0.0;
  numFlooredVariances = 0;
  refCount=1;
  emClearSharedBit();

  /////////////////////////////////////////////
  // make sure our caller has its covar accumulator resized
  // and initialized.
  componentsNextCovars.growIfNeeded(covariances.len());
  for (int i=0;i<covariances.len();i++) {
    componentsNextCovars[i] = 0.0;
  }

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

  // we return if both 1) the not training bit is set
  // and 2) there is no chance that this object will be shared.
  // If 1) is not true, we are training this object so we continue,
  // and if 2) is not true (we are sharing), then the accumulaters
  // created for this object will be needed by the other objects 
  // in this object's Gaussian component, so we'll need to compute them,
  // even though the parameters of this object will not be updated
  // when we swap them in.
  // if (numTimesShared == 1 && !emAmTrainingBitIsSet())
  // return;



  
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
 * emEndIterationNoSharingAlreadyNormalized()
 *      end the EM iteration for this var in the case that there is no
 *      sharing occuring, either mean sharing or covariance sharing.
 *      This version of the routine assumes that the covariance has
 *      already been normalized for us, but it has not yet been
 *      checked if the variances have fallen below floor (which we
 *      therefore check here).
 * 
 * Preconditions:
 *      see the assertions
 *
 * Postconditions:
 *      EM iteration has been ended.
 *
 * Side Effects:
 *      internal parameters are changed.
 *
 * Results:
 *      nil
 *
 *-----------------------------------------------------------------------
 */
void
DiagCovarVector::emEndIterationNoSharingAlreadyNormalized(const float *const parentsAccumulatedNextCovar)
{
  assert ( basicAllocatedBitIsSet() );

  // we return if both 1) the not training bit is set
  // and 2) there is no chance that this object will be shared.
  // If 1) is not true, we are training this object so we continue,
  // and if 2) is not true (we are sharing), then the accumulaters
  // created for this object will be needed by the other objects 
  // in this object's Gaussian component, so we'll need to compute them,
  // even though the parameters of this object will not be updated
  // when we swap them in.
  //if (numTimesShared == 1 && !emAmTrainingBitIsSet())
  // return;


  // if this isn't the case, something is wrong.
  assert ( emOnGoingBitIsSet() );

  // shouldn't be called when sharing occurs.
  assert ( refCount == 1 );
  assert (!emSharedBitIsSet());

  if (!emAccInitializedBitIsSet()) {
    nextCovariances.growIfNeeded(covariances.len());
    for (int i=0;i<covariances.len();i++) {
      nextCovariances[i] = 0.0;
    }
    emSetAccInitializedBit();
  }

  refCount = 0;

  accumulatedProbability.floor();
  if (accumulatedProbability < minContAccumulatedProbability()) {
    warning("WARNING: Diag covariance vec '%s' received only %e accumulated log probability in EM iteration, using previous values.",
	    name().c_str(),
	    accumulatedProbability.val());
    
    for (int i=0;i<covariances.len();i++)
      nextCovariances[i] = covariances[i];
  } else {

    // Finally, divide by N (see the equation above)
    // here computing the final variances.
    unsigned prevNumFlooredVariances = numFlooredVariances;
    double minVar = 0.0;
    for (int i=0;i<covariances.len();i++) {
      // "compute" the next variance
      nextCovariances[i] = parentsAccumulatedNextCovar[i];

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

      if (nextCovariances[i] < (float)GaussianComponent::varianceFloor()) {

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
	      name().c_str(),
	      numFlooredVariances-prevNumFlooredVariances,
	      minVar);
    }
  }

  // make it swapable
  emSetSwappableBit();

  // stop EM
  emClearOnGoingBit();
}




/*-
 *-----------------------------------------------------------------------
 * emEndIterationSharedMeansCovars()
 *      end the EM iteration for this var in the case that the
 *      covariances *and* means are shared amongst multiple arbitrary 
 *      Gaussians. 
 *      
 *      Note: this routine allows for an arbitrary set of
 *      Gaussians to share another arbitrary set of Means and a
 *      third arbitrary set of Covariances. It is more general
 *      than tying multiple means & variances together identicaly.
 *      In that case, state tieing should be used.
 *
 *      Note: this routine is a GEM rather than an EM.
 * 
 * Preconditions:
 *      see the assertions
 *
 * Postconditions:
 *      EM iteration has been ended.
 *
 * Side Effects:
 *      internal parameters are changed.
 *
 * Results:
 *      nil
 *
 *-----------------------------------------------------------------------
 */
void
DiagCovarVector::emEndIterationSharedMeansCovars(const logpr parentsAccumulatedProbability,
						 const float*const partialAccumulatedNextMeans,
						 const float *const partialAccumulatedNextCovar,
						 const MeanVector* mean)
{
  assert ( basicAllocatedBitIsSet() );


  // we return if both 1) the not training bit is set
  // and 2) there is no chance that this object will be shared.
  // If 1) is not true, we are training this object so we continue,
  // and if 2) is not true (we are sharing), then the accumulaters
  // created for this object will be needed by the other objects 
  // in this object's Gaussian component, so we'll need to compute them,
  // even though the parameters of this object will not be updated
  // when we swap them in.
  // if (numTimesShared == 1 && !emAmTrainingBitIsSet())
  // return;



  if (!emAccInitializedBitIsSet()) {
    nextCovariances.growIfNeeded(covariances.len());
    for (int i=0;i<covariances.len();i++) {
      nextCovariances[i] = 0.0;
    }
    emSetAccInitializedBit();
  }

  if (refCount > 0) {
    // if this isn't the case, something is wrong.
    assert ( emOnGoingBitIsSet() );

    if ( parentsAccumulatedProbability >
	 minContAccumulatedProbability()) {
      // Only accumulate here if there is something significant 
      // to accumlate.

      // accumulate the 1st and 2nd order statistics given
      // by the mean object.
      const double realAccumulatedProbability = 
	parentsAccumulatedProbability.unlog();

      const float* prev_mean_ptr = mean->means.ptr;
    
      for (int i=0;i<covariances.len();i++) {
	// do it in double precision
	const double tmp = ((double)partialAccumulatedNextCovar[i] 
			    - 2.0*(double)partialAccumulatedNextMeans[i]*(double)prev_mean_ptr[i]
			    + (double)prev_mean_ptr[i]*(double)prev_mean_ptr[i]*(double)realAccumulatedProbability);
	// convert and accumulate
	nextCovariances[i] += tmp;
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

  if (accumulatedProbability < minContAccumulatedProbability()) {
    warning("WARNING: shared diag covariance vec '%s' received only %e accumulated log probability (min is %e) in EM iteration, using previous values.",
	    name().c_str(),
	    accumulatedProbability.val(),
	    minContAccumulatedProbability().val());
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

      if (nextCovariances[i] < (float)GaussianComponent::varianceFloor()) {

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
      warning("WARNING: shared covariance vector named '%s' had %d variances floored, minimum variance found was %e.\n",
	      name().c_str(),
	      numFlooredVariances-prevNumFlooredVariances,
	      minVar);
    }

  }

  // make it swapable
  emSetSwappableBit();

  // stop EM
  emClearOnGoingBit();
}




/*-
 *-----------------------------------------------------------------------
 * emEndIterationSharedCovars
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
DiagCovarVector::emEndIterationSharedCovars(const float *const parentsAccumulatedNextCovar)
{
  assert ( basicAllocatedBitIsSet() );

  // we return if both 1) the not training bit is set
  // and 2) there is no chance that this object will be shared.
  // If 1) is not true, we are training this object so we continue,
  // and if 2) is not true (we are sharing), then the accumulaters
  // created for this object will be needed by the other objects 
  // in this object's Gaussian component, so we'll need to compute them,
  // even though the parameters of this object will not be updated
  // when we swap them in.
  // if (numTimesShared == 1 && !emAmTrainingBitIsSet())
  // return;



  if (!emAccInitializedBitIsSet()) {
    nextCovariances.growIfNeeded(covariances.len());
    for (int i=0;i<covariances.len();i++) {
      nextCovariances[i] = 0.0;
    }
    emSetAccInitializedBit();
  }


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

  accumulatedProbability.floor();
  if (accumulatedProbability < minContAccumulatedProbability()) {
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

      if (nextCovariances[i] < (float)GaussianComponent::varianceFloor()) {

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
	      name().c_str(),
	      numFlooredVariances-prevNumFlooredVariances,
	      minVar);
    }
  }

  // make it swapable
  emSetSwappableBit();

  // stop EM
  emClearOnGoingBit();
}






/*-
 *-----------------------------------------------------------------------
 * emEndIterationSharedMeansCovarsDlinks()
 *      end the EM iteration for this var in the case that the
 *      covariances, means, and dlinks are shared amongst multiple arbitrary 
 *      Gaussians. 
 *      
 *      Note: this routine allows for an arbitrary set of
 *      Gaussians to share another arbitrary set of Means and a
 *      third arbitrary set of Covariances. It is more general
 *      than tying multiple means & variances together identicaly.
 *      In that case, state tieing should be used.
 *
 *      Note: this routine is a GEM rather than an EM.
 * 
 * Preconditions:
 *      see the assertions
 *
 * Postconditions:
 *      EM iteration has been ended.
 *
 * Side Effects:
 *      internal parameters are changed.
 *
 * Results:
 *      nil
 *
 *-----------------------------------------------------------------------
 */
void
DiagCovarVector::emEndIterationSharedMeansCovarsDlinks(const logpr parentsAccumulatedProbability,
						       const float* const xAccumulators,
						       const float* const xxAccumulators,
						       const float* const xzAccumulators,
						       const float* const zAccumulators,
						       const float* const zzAccumulators,
						       const MeanVector* mean,
						       const DlinkMatrix* dLinkMat)
{
  assert ( basicAllocatedBitIsSet() );

  // we return if both 1) the not training bit is set
  // and 2) there is no chance that this object will be shared.
  // If 1) is not true, we are training this object so we continue,
  // and if 2) is not true (we are sharing), then the accumulaters
  // created for this object will be needed by the other objects 
  // in this object's Gaussian component, so we'll need to compute them,
  // even though the parameters of this object will not be updated
  // when we swap them in.
  // if (numTimesShared == 1 && !emAmTrainingBitIsSet())
  // return;


  if (!emAccInitializedBitIsSet()) {
    nextCovariances.growIfNeeded(covariances.len());
    for (int i=0;i<covariances.len();i++) {
      nextCovariances[i] = 0.0;
    }
    emSetAccInitializedBit();
  }

  if (refCount > 0) {
    // if this isn't the case, something is wrong.
    assert ( emOnGoingBitIsSet() );

    if ( parentsAccumulatedProbability >
	 minContAccumulatedProbability()) {
      // Only accumulate here if there is something significant 
      // to accumlate.


      const double realAccumulatedProbability = 
	parentsAccumulatedProbability.unlog();

      // previous mean
      const float* prev_mean_ptr = mean->means.ptr;
      // previous B matrix
      const float *prev_dlinkMatrix_ptr = dLinkMat->arr.ptr;

      // accumulator pointers
      const float*zAccumulators_ptr = zAccumulators;
      const float*xzAccumulators_ptr = xzAccumulators;
      const float*zzAccumulators_ptr = zzAccumulators;

      // expanded full zzmatrix
      sArray <double> zzExp;

      /////////////////////////////////////////////////////////////////
      // the comments in the following loop (i.e., term 1, 2, etc.)
      // refer to the sheet of equations containing this derivation.
      for (int i=0;i<covariances.len();i++) {

	const int nLinks = dLinkMat->numLinks(i);
	const double x = xAccumulators[i];

	// compute the various terms
	// term 0
	const double xx = xxAccumulators[i];

	// term 1
	const double x_u = x*prev_mean_ptr[i];

	// compute terms involving nLinks computations.
	double xz_B = 0.0;  // term 2
	double u_z_B = 0.0;   // term 3
	double B_zz_B = 0.0;   // term 4

	for (int j=0;j<nLinks;j++) {
	  // term 2 
	  xz_B += (double)(xzAccumulators_ptr[j])*
	          (double)(prev_dlinkMatrix_ptr[j]);
	  // term 3
	  u_z_B += (double)(zAccumulators_ptr[j])*
	          (double)(prev_dlinkMatrix_ptr[j]);

	  // term 4
	  double dotproduct = 0.0;
	  for (int k=0;k<nLinks;k++) {
	    double zz_tmp;
	    if (j > k) 
	      zz_tmp = 
		zzAccumulators_ptr[k*nLinks - k*(k+1)/2 + j];
	    else // k >= j
	      zz_tmp = 
		zzAccumulators_ptr[j*nLinks - j*(j+1)/2 + k];
	    dotproduct += 
	      zz_tmp*(double)(prev_dlinkMatrix_ptr[k]);
	  }
	  B_zz_B += dotproduct*prev_dlinkMatrix_ptr[j];
	}
	// finish term 3
	u_z_B *= prev_mean_ptr[i];
	// update pointers
	prev_dlinkMatrix_ptr += nLinks;
	xzAccumulators_ptr += nLinks;
	zAccumulators_ptr += nLinks;
	zzAccumulators_ptr += nLinks*(nLinks+1)/2;

	// term 5
	const double u_u = 
	  (double)prev_mean_ptr[i]*
	  (double)prev_mean_ptr[i]*
	  (double)realAccumulatedProbability;

	// convert and accumulate
	const double final_accumulator = 
	  (xx - 2.0*(x_u + xz_B - u_z_B) + B_zz_B + u_u);
	nextCovariances[i] += final_accumulator;
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

  if (accumulatedProbability < minContAccumulatedProbability()) {
    warning("WARNING: Shared diag covariance vec '%s' received only %e accumulated log probability (min is %e) in EM iteration, using previous values.",
	    name().c_str(),
	    accumulatedProbability.val(),
	    minContAccumulatedProbability().val());
    for (int i=0;i<covariances.len();i++) 
      nextCovariances[i] = covariances[i];
  } else {

    // TODO: should check for possible overflow here of 
    // accumulatedProbability when we do the inverse and unlog.
    // Ideally this won't happen for a given minContAccumulatedProbability().
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

      if (nextCovariances[i] < (float)GaussianComponent::varianceFloor()) {

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
      warning("WARNING: shared covariance vector named '%s' had %d variances floored, minimum variance found was %e.\n",
	      name().c_str(),
	      numFlooredVariances-prevNumFlooredVariances,
	      minVar);
    }

  }

  // make it swapable
  emSetSwappableBit();

  // stop EM
  emClearOnGoingBit();
}





/*-
 *-----------------------------------------------------------------------
 * emEndIterationSharedCovars()
 *      end the EM iteration for this var in the case that the
 *      covariances are shared amongst multiple Gaussians but
 *      where each such Gaussian has its own mean (i.e., no mean sharing).
 *      This routine takes the partially accumulated 1st and 2nd moments,
 *      and the corresponding probability.
 * 
 * Preconditions:
 *      see the assertions
 *
 * Postconditions:
 *      EM iteration has been ended.
 *
 * Side Effects:
 *      internal parameters are changed.
 *
 * Results:
 *      nil
 *
 *-----------------------------------------------------------------------
 */
void
DiagCovarVector::emEndIterationSharedCovars(const logpr parentsAccumulatedProbability,
					    const float*const partialAccumulatedNextMeans,
					    const float *const partialAccumulatedNextCovar)
{
  assert ( basicAllocatedBitIsSet() );


  // we return if both 1) the not training bit is set
  // and 2) there is no chance that this object will be shared.
  // If 1) is not true, we are training this object so we continue,
  // and if 2) is not true (we are sharing), then the accumulaters
  // created for this object will be needed by the other objects 
  // in this object's Gaussian component, so we'll need to compute them,
  // even though the parameters of this object will not be updated
  // when we swap them in.
  // if (numTimesShared == 1 && !emAmTrainingBitIsSet())
  // return;


  if (!emAccInitializedBitIsSet()) {
    nextCovariances.growIfNeeded(covariances.len());
    for (int i=0;i<covariances.len();i++) {
      nextCovariances[i] = 0.0;
    }
    emSetAccInitializedBit();
  }

  if (refCount > 0) {
    // if this isn't the case, something is wrong.
    assert ( emOnGoingBitIsSet() );

    if ( parentsAccumulatedProbability >
	 minContAccumulatedProbability()) {
      // Only accumulate here if there is something significant 
      // to accumlate.

      // accumulate the 1st and 2nd order statistics given
      // by the mean object.
      const double invRealAccumulatedProbability = 
	parentsAccumulatedProbability.inverse().unlog();
    
      for (int i=0;i<covariances.len();i++) {
	const double tmp = 
	  ((double)partialAccumulatedNextCovar[i] - 
	   (double)partialAccumulatedNextMeans[i]*
	   (double)partialAccumulatedNextMeans[i]*
	   (double)invRealAccumulatedProbability);
	nextCovariances[i] += tmp;
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

  if (accumulatedProbability < minContAccumulatedProbability()) {
    warning("WARNING: Diag covariance vec '%s' received only %e accumulated log probability (min is %e) in EM iteration, using previous values.",
	    name().c_str(),
	    accumulatedProbability.val(),
	    minContAccumulatedProbability().val());
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

      if (nextCovariances[i] < (float)GaussianComponent::varianceFloor()) {

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
	      name().c_str(),
	      numFlooredVariances-prevNumFlooredVariances,
	      minVar);
    }

  }

  // make it swapable
  emSetSwappableBit();

  // stop EM
  emClearOnGoingBit();
}



/*-
 *-----------------------------------------------------------------------
 * emEndIterationNoSharing
 *      end the EM iteration for this var in the case that there is no
 *      sharing occuring, either mean sharing or covariance sharing.
 * 
 * Preconditions:
 *      see the assertions
 *
 * Postconditions:
 *      EM iteration has been ended.
 *
 * Side Effects:
 *      internal parameters are changed.
 *
 * Results:
 *      nil
 *
 *-----------------------------------------------------------------------
 */
void
DiagCovarVector::emEndIterationNoSharing(const float*const partialAccumulatedNextMeans,
					 const float *const partialAccumulatedNextCovar)
{
  assert ( basicAllocatedBitIsSet() );

  // we return if both 1) the not training bit is set
  // and 2) there is no chance that this object will be shared.
  // If 1) is not true, we are training this object so we continue,
  // and if 2) is not true (we are sharing), then the accumulaters
  // created for this object will be needed by the other objects 
  // in this object's Gaussian component, so we'll need to compute them,
  // even though the parameters of this object will not be updated
  // when we swap them in.
  // if (numTimesShared == 1 && !emAmTrainingBitIsSet())
  //return;


  // if this isn't the case, something is wrong.
  assert ( emOnGoingBitIsSet() );

  // shouldn't be called when sharing occurs.
  assert ( refCount == 1 );
  assert (!emSharedBitIsSet());

  if (!emAccInitializedBitIsSet()) {
    nextCovariances.growIfNeeded(covariances.len());
    for (int i=0;i<covariances.len();i++) {
      nextCovariances[i] = 0.0;
    }
    emSetAccInitializedBit();
  }

  refCount = 0;
  // we're ready to finish and compute the next covariances.

  if (accumulatedProbability < minContAccumulatedProbability()) {
    warning("WARNING: Diag covariance vec '%s' received only %e accumulated log probability (min is %e) in EM iteration, using previous values.",
	    name().c_str(),
	    accumulatedProbability.val(),
	    minContAccumulatedProbability().val());
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
      const double tmp = 
	((double)partialAccumulatedNextCovar[i] - 
          	 (double)partialAccumulatedNextMeans[i]*
	         (double)partialAccumulatedNextMeans[i]*(double)invRealAccumulatedProbability)
	*(double)invRealAccumulatedProbability;
      nextCovariances[i] = tmp;

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

      if (nextCovariances[i] < (float)GaussianComponent::varianceFloor()) {

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
	      name().c_str(),
	      numFlooredVariances-prevNumFlooredVariances,
	      minVar);
    }

  }

  // make it swapable
  emSetSwappableBit();

  // stop EM
  emClearOnGoingBit();
}







void
DiagCovarVector::emSwapCurAndNew()
{
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
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



/*-
 *-----------------------------------------------------------------------
 *
 * Accumulator loading/storing routines for parallel training support.
 *
 *-----------------------------------------------------------------------
 */


void
DiagCovarVector::emStoreAccumulators(oDataStreamFile& ofile)
{
  assert ( basicAllocatedBitIsSet() );
  if (numTimesShared == 1 && !emAmTrainingBitIsSet()) {
    // then we are not training, because
    // we have turned off training of this object.
    // We write out '0' to state that 
    // there are no values stored for this object.
    unsigned flag = 0;
    ofile.write(flag,"writing acc flag");
    return;
  } else {
    // either the training bit is set, or the training bit is not set
    // but this object is being shared more than once. In the former
    // case, we of course write out the accumulators. In the latter
    // case, while this object won't change, it's accumulators might
    // be needed by another object for which the training bit is set.
    if (accumulatedProbability.zero()) {
      // then we indeed have no probability values, so lets emit a warning
      warning("WARNING: zero accumulator values for %s '%s'\n",
	      typeName().c_str(),
	      name().c_str());
      // We write out '0' to state that 
      // there are no values stored for this object.
      unsigned flag = 0;
      ofile.write(flag,"writing acc flag");
    } else {
      // we write a 1 to indicate that there are accumulators
      // stored for this object.
      unsigned flag = 1;
      ofile.write(flag,"writing acc flag");
      // store the accumulators as normal.
      ofile.write(accumulatedProbability.val(),"EM store accums");
      // call virtual function to do actual work for object.
      emStoreObjectsAccumulators(ofile);
    }
  }
}







