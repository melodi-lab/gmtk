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
    means[i] = rnd.drand48();
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
MeanVector::emStartIteration()
{
  if (!GM_Parms.amTrainingMeans())
    return;

  if(emOnGoingBitIsSet())
    return; // already done

  if (!emEmAllocatedBitIsSet()) {
    // this is presumably the first time
    emSetEmAllocatedBit();
    nextMeans.resize(means.len());
  }
  // EM iteration is now going.
  emSetOnGoingBit();
  emSetSwappableBit();

  accumulatedProbability = 0.0;
  for (int i=0;i<means.len();i++) {
    nextMeans[i] = 0.0;
  }
}


void
MeanVector::emIncrement(logpr prob,
			const float* const f,
			const Data32* const base,
			const int stride)
{
  if (!GM_Parms.amTrainingMeans())
    return;

  emStartIteration();

  if (prob.val() < log_FLT_MIN) {
    return;
    // don't accumulate anything since this one is so small and
    // if we did an unlog() and converted to a single precision
    // floating point number, it would be a denomral.
  } 

  accumulatedProbability += prob;
  const float f_prob = prob.unlog();

  float * means_p = nextMeans.ptr;
  float * means_end_p = nextMeans.ptr + nextMeans.len();
  const float * f_p = f;
  do {
    *means_p += (*f_p)*f_prob;
    means_p++;
    f_p++;
  } while (means_p != means_end_p);

}


void
MeanVector::emEndIteration()
{
  if (!GM_Parms.amTrainingMeans())
    return;

  if ( !emOnGoingBitIsSet() )
    return; 

  if (accumulatedProbability.zero()) {
    warning("WARNING: Mean vector named '%s' did not receive any accumulated probability in EM iteration",name().c_str());
  }

  const double invRealAccumulatedProbability = 
    accumulatedProbability.inverse().unlog();
  // finish computing the next means.
  float * means_p = nextMeans.ptr;
  float * means_end_p = nextMeans.ptr + nextMeans.len();
  do {
    *means_p *= invRealAccumulatedProbability;
    means_p++;
  } while (means_p != means_end_p);

  // stop EM
  emClearOnGoingBit();
}


void
MeanVector::emSwapCurAndNew()
{
  if (!GM_Parms.amTrainingMeans())
    return;

  if (!emSwappableBitIsSet())
    return;
  for (int i=0;i<means.len();i++) {
    genSwap(means[i],nextMeans[i]);
  }
  emSetSwappableBit();
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

