/*-
 * GMTK_MeanVector.cc
 *     A vector used for means of Gaussians.
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

#include "general.h"
#include "error.h"
#include "rand.h"

#include "GMTK_MeanVector.h"
#include "GMTK_GMParms.h"
#include "GMTK_GaussianComponent.h"
#include "GMTK_MixtureCommon.h"
#include "GMTK_DiagCovarVector.h"
#include "GMTK_DlinkMatrix.h"

VCID("$Header$");

////////////////////////////////////////////////////////////////////
//        Static members
////////////////////////////////////////////////////////////////////

double MeanVector::cloneSTDfrac = 0.1;

void MeanVector::checkForValidValues()
{
  if (MeanVector::cloneSTDfrac < 0)
    error("ERROR: MeanVector's cloneSTDfrac (%e) must be >= 0",
	  MeanVector::cloneSTDfrac);
}

////////////////////////////////////////////////////////////////////
//        General create, read, destroy routines 
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * MeanVector::MeanVector()
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
MeanVector::MeanVector()
{}



/*-
 *-----------------------------------------------------------------------
 * MeanVector::read(is)
 *      read in the mean array from file 'is'. 
 * 
 * Results:
 *      No results.
 *
 * Side Effects:
 *
 *-----------------------------------------------------------------------
 */
void MeanVector::read(iDataStreamFile& is) { 
  NamedObject::read(is);
  int length;
  is.read(length,"MeanVector::read, distribution length");
  if (length <= 0)
    error("ERROR: mean vector %s specifies length (%d) < 0 in input. Must be positive.",
	    name().c_str(),length);
  means.resize(length);
  for (int i=0;i<length;i++) {
    is.read(means[i],":reading mean values");
  }
  setBasicAllocatedBit();
  numTimesShared = 0;
  refCount = 0;
}


/*-
 *-----------------------------------------------------------------------
 * MeanVector::write(is)
 *      read in the mean array from file 'is'. 
 * 
 * Results:
 *      No results.
 *
 * Side Effects:
 *
 *-----------------------------------------------------------------------
 */
void MeanVector::write(oDataStreamFile& os) { 
  NamedObject::write(os);
  os.write(means.len(),"mean vector write length");
  for (int i=0;i<means.len();i++) {
    os.write(means[i],"mean vector write, values");
  }
  os.nl();
}


////////////////////////////////////////////////////////////////////
//        Misc Support
////////////////////////////////////////////////////////////////////

void
MeanVector::makeRandom()
{
  if (!emAmTrainingBitIsSet())
    return;

  for (int i=0;i<means.len();i++) {
    means[i] = rnd.drand48pe();
  }
}

void
MeanVector::makeUniform()
{
  if (!emAmTrainingBitIsSet())
    return;

  for (int i=0;i<means.len();i++) {
    means[i] = 0.0;
  }
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
MeanVector*
MeanVector::noisyClone()
{
  assert ( basicAllocatedBitIsSet() );

  // first check if self is already cloned, and if so, return that.
  MeanVector* clone;

  map<MeanVector*,MeanVector*>::iterator it = MixtureCommon::meanCloneMap.find(this);
  if (it == MixtureCommon::meanCloneMap.end()) {
    clone = new MeanVector();
    // make sure we get a unique name
    unsigned cloneNo=0; do {
      char buff[256];
      sprintf(buff,"%d",cloneNo);
      clone->_name = _name + string("_cl") + buff;
      cloneNo++;
    } while (GM_Parms.meansMap.find(clone->_name) != GM_Parms.meansMap.end());
    clone->refCount = 0;
    clone->numTimesShared = 0;
    clone->means.resize(means.len());
    for (int i=0;i<means.len();i++) {
      clone->means[i] = means[i] + 
	cloneSTDfrac*means[i]*rnd.normal();
    }
    clone->setBasicAllocatedBit();
    MixtureCommon::meanCloneMap[this] = clone;

    // also add self to GMParms object.
    GM_Parms.add(clone);

  } else {
    clone = (*it).second;
  }
  return clone;
}


/////////////////
// EM routines //
/////////////////


void
MeanVector::emStartIteration(sArray<float>& componentsNextMeans)
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
    componentsNextMeans.growIfNeeded(means.len());
    for (int i=0;i<means.len();i++) {
      componentsNextMeans[i] = 0.0;
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
  refCount = 1;
  emClearSharedBit();

  /////////////////////////////////////////////
  // make sure our caller has its mean accumulator resized
  // and initialized.
  componentsNextMeans.growIfNeeded(means.len());
  for (int i=0;i<means.len();i++) {
    componentsNextMeans[i] = 0.0;
  }

}


void
MeanVector::emIncrement(const logpr prob,
			const float fprob,
			const float* const f,
			const Data32* const base,
			const int stride,
			float *const partialAccumulatedNextMeans)
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
  //    return;


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

#if 0
  // this accumulation is now done by caller in conjunction
  // with the covar vector.
  // This routine is called often so we use pointer arith.
  float * nextMeans_p = partialAccumulatedNextMeans;
  float * nextMeans_end_p = partialAccumulatedNextMeans + means.len();
  const float * f_p = f;
  do {
    *nextMeans_p += (*f_p)*fprob;
    nextMeans_p++;
    f_p++;
  } while (nextMeans_p != nextMeans_end_p);
#endif

}

/*-
 *-----------------------------------------------------------------------
 * emEndIterationSharedMeansCovarsDlinks()
 *      end the EM iteration for this mean object, where we have 
 *      shared means, shared covariances, and shared dlink matrices.
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
MeanVector::emEndIterationSharedMeansCovarsDlinks(const logpr parentsAccumulatedProbability,
						  const float*const xAccumulators,
						  const float*const zAccumulators,
						  const DlinkMatrix* dLinkMat,
						  const DiagCovarVector* covar)
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
    // make sure next-means are set up
    nextMeans.growIfNeeded(means.len());
    for (int i=0;i<means.len();i++) {
      nextMeans[i] = 0.0;
    }
    sharedMeansDenominator.growIfNeeded(means.len());
    for (int i=0;i<means.len();i++) {
      sharedMeansDenominator[i] = 0.0;
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
      
      // grab a pointer to the previous inverse variances
      // needed for normalization.
      const float* previous_variances_inv_ptr = covar->variances_inv.ptr;

      const double realAccumulatedProbability = 
	parentsAccumulatedProbability.unlog();

      // grab a pointer to the previous B matrix
      const float *prev_dlinkMatrix_ptr = dLinkMat->arr.ptr;

      // also get a pointer to the z accumulators
      const float*zAccumulators_ptr = zAccumulators;
      

      // accumulate in the 1st order statistics
      for (int i=0;i<means.len();i++) {
	
	const int nLinks = dLinkMat->numLinks(i);

	// first compute the B*zAccumulator value
	double tmp=0.0;
	for (int j=0;j<nLinks;j++) {
	  tmp += (double)(*prev_dlinkMatrix_ptr++)*(double)(*zAccumulators_ptr++);
	}

	// subtract off the Bz from x
	tmp = xAccumulators[i] - tmp;

	// and accumulate the result multiplying by the previous iteration's inverse variance
	nextMeans[i] += tmp*previous_variances_inv_ptr[i];

	// and get the denominator
	sharedMeansDenominator[i] += previous_variances_inv_ptr[i]*realAccumulatedProbability;
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

  accumulatedProbability.floor();
  if (accumulatedProbability < minContAccumulatedProbability()) {
    infoMsg(IM::Warning,"WARNING: Shared mean vec '%s' received only a total of %e accumulated log probability (min is %e) in EM iteration, using previous means",name().c_str(),
	    accumulatedProbability.val(),
	    minContAccumulatedProbability().val());
    for (int i=0;i<nextMeans.len();i++)
      nextMeans[i] = means[i];
  } else {
    unsigned previousMeansUsed = 0;
    // finish computing the next means.
    for (int i=0;i<means.len();i++) {
      // first make sure denominator is well behaved
      if (sharedMeansDenominator[i] <= DBL_MIN) {
	// use previous mean for this iteration.
	nextMeans[i] = means[i];
	previousMeansUsed++;
      } else {
	// make sure the division is done in double precision.
	const double tmp = (double) nextMeans[i] / (double) sharedMeansDenominator[i];
	if (tmp >= FLT_MAX) {
	  // use previous mean for this iteration.
	  nextMeans[i] = means[i];
	  previousMeansUsed++;
	} else {
	  // then covert it to nextMeans type.
	  nextMeans[i] = tmp;
	}
      }
    }
    if (previousMeansUsed > 0) 
      infoMsg(IM::Warning,"WARNING: Shared mean vec '%s' used %d previous means values because of low counts.",
	      name().c_str(),
	      previousMeansUsed);
  }

  // make it swapable
  emSetSwappableBit();

  // stop EM
  emClearOnGoingBit();
}




/*-
 *-----------------------------------------------------------------------
 * emEndIterationSharedMeansCovars()
 *      end the EM iteration for this mean object, where we have both
 *      shared means and shared covariances.
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
MeanVector::emEndIterationSharedMeansCovars(const logpr parentsAccumulatedProbability,
						  const float*const partialAccumulatedNextMeans,
						  const DiagCovarVector* covar)
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
    // make sure next-means are set up
    nextMeans.growIfNeeded(means.len());
    for (int i=0;i<means.len();i++) {
      nextMeans[i] = 0.0;
    }
    sharedMeansDenominator.growIfNeeded(means.len());
    for (int i=0;i<means.len();i++) {
      sharedMeansDenominator[i] = 0.0;
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

      // grab a pointer to the previous inverse variances
      // needed for normalization.
      const float* previous_variances_inv_ptr = covar->variances_inv.ptr;

      const double realAccumulatedProbability = 
	parentsAccumulatedProbability.unlog();

      // accumulate in the 1st order statistics given
      // by the mean object.
      for (int i=0;i<means.len();i++) {
	nextMeans[i] += 
	  (partialAccumulatedNextMeans[i]*previous_variances_inv_ptr[i]);
	sharedMeansDenominator[i] += previous_variances_inv_ptr[i]*realAccumulatedProbability;
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

  accumulatedProbability.floor();
  if (accumulatedProbability < minContAccumulatedProbability()) {
    infoMsg(IM::Warning,"WARNING: Shared mean vec '%s' received only a total of %e accumulated log probability (min is %e) in EM iteration, using previous means",name().c_str(),
	    accumulatedProbability.val(),
	    minContAccumulatedProbability().val());
    for (int i=0;i<nextMeans.len();i++)
      nextMeans[i] = means[i];
  } else {
    unsigned previousMeansUsed = 0;
    // finish computing the next means.
    for (int i=0;i<means.len();i++) {
      // first make sure denominator is well behaved
      if (sharedMeansDenominator[i] <= DBL_MIN) {
	// use previous mean for this iteration.
	nextMeans[i] = means[i];
	previousMeansUsed++;
      } else {
	// make sure the division is done in double precision.
	const double tmp = (double)nextMeans[i] / (double)sharedMeansDenominator[i];
	// then covert it to nextMeans type.
	nextMeans[i] = tmp;
      }
    }
    if (previousMeansUsed > 0) 
      infoMsg(IM::Warning,"WARNING: Shared mean vec '%s' used %d previous means values because of low counts.",
	      name().c_str(),
	      previousMeansUsed);
  }

  // make it swapable
  emSetSwappableBit();

  // stop EM
  emClearOnGoingBit();
}



/*-
 *-----------------------------------------------------------------------
 * emEndIterationNoSharing()
 *      end the EM iteration for this mean object, where we have no
 *      sharing.
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
MeanVector::emEndIterationNoSharing(const float*const partialAccumulatedNextMeans)
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

  // if this isn't the case, something is wrong.
  assert ( emOnGoingBitIsSet() );

  // shouldn't be called when sharing occurs.
  assert ( refCount == 1 );
  assert (!emSharedBitIsSet());


  if (!emAccInitializedBitIsSet()) {
    // make sure next-means are set up
    nextMeans.growIfNeeded(means.len());
    for (int i=0;i<means.len();i++) {
      nextMeans[i] = 0.0;
    }
    emSetAccInitializedBit();
  }
  
  refCount = 0;

  accumulatedProbability.floor();
  if (accumulatedProbability < minContAccumulatedProbability()) {
    infoMsg(IM::Warning,"WARNING: Mean vec '%s' received only %e accumulated log probability (min is %e) in EM iteration, using previous means",name().c_str(),
	    accumulatedProbability.val(),
	    minContAccumulatedProbability().val());
    for (int i=0;i<nextMeans.len();i++)
      nextMeans[i] = means[i];
  } else {
    const double invRealAccumulatedProbability =
      accumulatedProbability.inverse().unlog();
    // finish computing the next means.
    for (int i=0;i<means.len();i++) {
      nextMeans[i] = partialAccumulatedNextMeans[i]*invRealAccumulatedProbability;
    }
  }

  // make it swapable
  emSetSwappableBit();

  // stop EM
  emClearOnGoingBit();
}



/*-
 *-----------------------------------------------------------------------
 * emEndIterationNoSharingAlreadyNormalized()
 *      end the EM iteration for this mean object, where we have no
 *      sharing, but the mean has already been normalized. We still
 *      need to check for a small accumulator probability in which case we 
 *      just use the previous mean (rather than the new one).
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
MeanVector::emEndIterationNoSharingAlreadyNormalized(const float*const accumulatedNextMeans)
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


  // if this isn't the case, something is wrong.
  assert ( emOnGoingBitIsSet() );

  // shouldn't be called when sharing occurs.
  assert ( refCount == 1 );
  assert (!emSharedBitIsSet());


  if (!emAccInitializedBitIsSet()) {
    // make sure next-means are set up
    nextMeans.growIfNeeded(means.len());
    for (int i=0;i<means.len();i++) {
      nextMeans[i] = 0.0;
    }
    emSetAccInitializedBit();
  }

  refCount = 0;

  accumulatedProbability.floor();
  if (accumulatedProbability < minContAccumulatedProbability()) {
    infoMsg(IM::Warning,"WARNING: Mean vec '%s' received only %e accumulated log probability (min is %e) in EM iteration, using previous means",name().c_str(),
	    accumulatedProbability.val(),
	    minContAccumulatedProbability().val());
    for (int i=0;i<nextMeans.len();i++)
      nextMeans[i] = means[i];
  } else {
    for (int i=0;i<means.len();i++) {
      nextMeans[i] = accumulatedNextMeans[i];
    }
  }


  // make it swapable
  emSetSwappableBit();

  // stop EM
  emClearOnGoingBit();
}




void
MeanVector::emSwapCurAndNew()
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
  for (int i=0;i<means.len();i++) {
    genSwap(means[i],nextMeans[i]);
  }
  // make no longer swappable
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
MeanVector::emStoreAccumulators(oDataStreamFile& ofile)
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
      infoMsg(IM::SoftWarning,"WARNING: zero accumulator values for %s '%s'\n",
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

