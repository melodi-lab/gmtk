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
#include "GMTK_MixGaussiansCommon.h"
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

  map<MeanVector*,MeanVector*>::iterator it = MixGaussiansCommon::meanCloneMap.find(this);
  if (it == MixGaussiansCommon::meanCloneMap.end()) {
    clone = new MeanVector();
    // make sure we get a unique name
    unsigned cloneNo=0; do {
      char buff[256];
      sprintf(buff,"%d",cloneNo);
      clone->_name = _name + string("_cl") + buff;
      cloneNo++;
    } while (GM_Parms.meansMap.find(clone->_name) != GM_Parms.meansMap.end());
    clone->refCount = 0;
    clone->means.resize(means.len());
    for (int i=0;i<means.len();i++) {
      clone->means[i] = means[i] + 
	cloneSTDfrac*means[i]*rnd.normal();
    }
    clone->setBasicAllocatedBit();
    MixGaussiansCommon::meanCloneMap[this] = clone;

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
  // if (!emAmTrainingBitIsSet())
  // return;

  /////////////////////////////////////////////
  // make sure our caller has its mean accumulator resized
  // and initialized.
  componentsNextMeans.growIfNeeded(means.len());
  for (int i=0;i<means.len();i++) {
    componentsNextMeans[i] = 0.0;
  }

  if(emOnGoingBitIsSet()) {
    // EM already on going.
    // Increment the count of number of Gaussian Components using this mean.
    refCount++; 
    // this object therefore is shared, set the bit saying so.
    emSetSharedBit();
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
  // if (!emAmTrainingBitIsSet())
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

  // This routine is called often so we use pointer arith.
  float * nextMeans_p = partialAccumulatedNextMeans;
  float * nextMeans_end_p = partialAccumulatedNextMeans + means.len();
  const float * f_p = f;
  do {
    *nextMeans_p += (*f_p)*fprob;
    nextMeans_p++;
    f_p++;
  } while (nextMeans_p != nextMeans_end_p);

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
  // if (!emAmTrainingBitIsSet())
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
    warning("WARNING: Shared mean vec '%s' received only a total of %e accumulated log probability (min is %e) in EM iteration, using previous means",name().c_str(),
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
      warning("WARNING: Shared mean vec '%s' used %d previous means values because of low counts.",
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
  // if (!emAmTrainingBitIsSet())
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
    warning("WARNING: Shared mean vec '%s' received only a total of %e accumulated log probability (min is %e) in EM iteration, using previous means",name().c_str(),
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
      warning("WARNING: Shared mean vec '%s' used %d previous means values because of low counts.",
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
  // if (!emAmTrainingBitIsSet())
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
    warning("WARNING: Mean vec '%s' received only %e accumulated log probability (min is %e) in EM iteration, using previous means",name().c_str(),
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
  // if (!emAmTrainingBitIsSet())
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
    warning("WARNING: Mean vec '%s' received only %e accumulated log probability (min is %e) in EM iteration, using previous means",name().c_str(),
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


void
MeanVector::emStoreAccumulators(oDataStreamFile& ofile)
{
  assert ( basicAllocatedBitIsSet() );
  // if (!emAmTrainingBitIsSet())
  // return;

  if ( !emEmAllocatedBitIsSet() ) {
    warning("WARNING: storing zero accumulators for mean '%s'\n",
	    name().c_str());
    emStoreZeroAccumulators(ofile);
    return;
  }
  EMable::emStoreAccumulators(ofile);
}


void
MeanVector::emStoreZeroAccumulators(oDataStreamFile& ofile)
{
  // if (!emAmTrainingBitIsSet())
  // return;

  assert ( basicAllocatedBitIsSet() );
  EMable::emStoreZeroAccumulators(ofile);
}


void
MeanVector::emLoadAccumulators(iDataStreamFile& ifile)
{
  // if (!emAmTrainingBitIsSet())
  // return;

  assert ( basicAllocatedBitIsSet() );
  assert ( emEmAllocatedBitIsSet() );
  EMable::emLoadAccumulators(ifile);
}


void
MeanVector::emAccumulateAccumulators(iDataStreamFile& ifile)
{
  // if (!emAmTrainingBitIsSet())
  // return;

  assert ( basicAllocatedBitIsSet() );
  assert ( emEmAllocatedBitIsSet() );
  EMable::emAccumulateAccumulators(ifile);
}


