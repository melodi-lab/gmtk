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

VCID("$Header$");


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




/////////////////
// EM routines //
/////////////////



void
MeanVector::emStartIteration(sArray<float>& componentsNextMeans)
{
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

  if (accumulatedProbability.zero()) {
    // TODO: need to check if this will overflow here
    // when dividing by it. This is more than just checking
    // for zero. Also need to do this in every such EM object.
    warning("WARNING: Mean vector named '%s' did not receive any accumulated probability in EM iteration, using previous means",name().c_str());
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
  error("not implemented");
}

void
MeanVector::emLoadAccumulators(iDataStreamFile& ifile)
{
  error("not implemented");
}


void
MeanVector::emAccumulateAccumulators(iDataStreamFile& ifile)
{
  error("not implemented");
}


