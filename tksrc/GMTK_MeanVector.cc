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
#include <ieeefp.h>
#include <float.h>
#include <assert.h>

#include "general.h"
#include "error.h"
#include "rand.h"

#include "GMTK_MeanVector.h"
#include "GMTK_GMParms.h"
#include "GMTK_GaussianComponent.h"
#include "GMTK_MixGaussiansCommon.h"

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
  for (int i=0;i<means.len();i++) {
    means[i] = rnd.drand48pe();
  }
}

void
MeanVector::makeUniform()
{
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

  if (!GM_Parms.amTrainingMeans())
    return;

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
    return;
  }

  if (!emEmAllocatedBitIsSet()) {
    // this is presumably the first time
    emSetEmAllocatedBit();
    // allocate the final next means if needed
    nextMeans.growIfNeeded(means.len());
  }

  // EM iteration is now going.
  emSetOnGoingBit();

  accumulatedProbability = 0.0;
  refCount = 1;
  for (int i=0;i<means.len();i++) {
    nextMeans[i] = 0.0;
  }

  // make it swapable, although at this point
  // it would swap in the unaccumulated values.
  emSetSwappableBit();
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

  if (!GM_Parms.amTrainingMeans())
    return;
  
  /////////////////////////////////////////////
  // Note: unlike the normal EM mode described
  // in GMTK_EMable.h, we do not call
  // emStartIteration() here and assume that it
  // was called by the Gaussian component that
  // is using this mean.

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


void
MeanVector::emEndIteration(const float*const partialAccumulatedNextMeans)
{
  assert ( basicAllocatedBitIsSet() );

  if (!GM_Parms.amTrainingMeans())
    return;

  if (refCount > 0) {
    // if this isn't the case, something is wrong.
    assert ( emOnGoingBitIsSet() );

    // accumulate in the 1st order statistics given
    // by the mean object.
    for (int i=0;i<means.len();i++) {
      nextMeans[i] += partialAccumulatedNextMeans[i];
    }

    refCount--;
  }

  /////////////////////////////////////////////
  // if there is still someone who
  // has not given us his/her 1st order stats,
  // then we return w/o finishing.
  if (refCount > 0)
    return;

  if (accumulatedProbability < GaussianComponent::minAccumulatedProbability()) {
    warning("WARNING: Mean vec '%s' received only %e accumulated log probability in EM iteration, using previous means",
	    accumulatedProbability.val(),name().c_str());
    for (int i=0;i<nextMeans.len();i++)
      nextMeans[i] = means[i];
  } else {
    const double invRealAccumulatedProbability =
      accumulatedProbability.inverse().unlog();
    // finish computing the next means.
    float * nextMeans_p = nextMeans.ptr;
    float * nextMeans_end_p = nextMeans.ptr + nextMeans.len();
    do {
      *nextMeans_p *= invRealAccumulatedProbability;
      nextMeans_p++;
    } while (nextMeans_p != nextMeans_end_p);
  }

  // stop EM
  emClearOnGoingBit();
}


void
MeanVector::emSwapCurAndNew()
{
  assert ( basicAllocatedBitIsSet() );

  if (!GM_Parms.amTrainingMeans())
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
  error("not implemented");
}

void
MeanVector::emLoadAccumulators(iDataStreamFile& ifile)
{
  assert ( basicAllocatedBitIsSet() );
  error("not implemented");
}


void
MeanVector::emAccumulateAccumulators(iDataStreamFile& ifile)
{
  assert ( basicAllocatedBitIsSet() );
  error("not implemented");
}


