/*-
 * GMTK_LinMeanCondDiagGaussian.cc
 *        Code for plain vanilla diagonal Gaussians.
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
#include <ieeefp.h>
#include <float.h>
#include <assert.h>
#include <ctype.h>

#include <string>

#include "general.h"
VCID("$Header$");
#include "error.h"
#include "rand.h"

#include "GMTK_LinMeanCondDiagGaussian.h"
#include "GMTK_GMParms.h"




void
LinMeanCondDiagGaussian::read(iDataStreamFile& is)
{
  // The index in the global mean array of this mean.
  int meanIndex; 
  // The index in the global variance array of this variance vector
  int covarIndex;

  // read name
  NamedObject::read(is);

  // read mean vector
  string str;
  is.read(str);

  if (GM_Parms.meansMap.find(str) ==  GM_Parms.meansMap.end()) 
      error("Error: LinMeanCondDiagGaussian '%s' in file '%s' specifies mean name '%s' that does not exist",
	    _name.c_str(),is.fileName(),str.c_str());
  meanIndex = GM_Parms.meansMap[str];
  mean = GM_Parms.means[meanIndex];


  // read covariance vector
  is.read(str);
  if (GM_Parms.covarsMap.find(str) == GM_Parms.covarsMap.end())
    error("Error: LinMeanCondDiagGaussian '%s' in file '%s' specifies covar name '%s' that does not exist",
	  _name.c_str(),is.fileName(),str.c_str());
  
  covarIndex = GM_Parms.covarsMap[str];
  covar = GM_Parms.covars[covarIndex];


  // read the dlink structure 
  is.read(str);
  if (GM_Parms.dLinksMap.find(str) == GM_Parms.dLinksMap.end())
    error("Error: LinMeanCondDiagGaussian '%s' in file '%s' specifies dlink name '%s' that does not exist",
	  _name.c_str(),is.fileName(),str.c_str());
  
  dLinks = GM_Parms.dLinks[GM_Parms.dLinksMap[str]];


  // read the dlink matrix parameter values
  is.read(str);
  if (GM_Parms.dLinkMatsMap.find(str) == GM_Parms.dLinkMatsMap.end())
    error("Error: LinMeanCondDiagGaussian '%s' in file '%s' specifies dlink matrix name '%s' that does not exist",
	  _name.c_str(),is.fileName(),str.c_str());
  
  dLinkMat = GM_Parms.dLinkMats[GM_Parms.dLinkMatsMap[str]];

  if (!dLinks->compatibleWith(*dLinkMat))
    error("Error: LinMeanCondDiagGaussian '%s' in file '%s' specifices a dlink structure '%s' and matrix '%s' that are not compatible with each other\n",
	  _name.c_str(),is.fileName(),
	  dLinks->name().c_str(),
	  dLinkMat->name().c_str());

  // check that lengths match, etc.
  if (covar->dim() != mean->dim()) {
    error("Error: LinMeanCondDiagGaussian '%s' in file '%s' specifices a mean '%s' with dim %d and covariance '%s' with dim '%d'\n",
	  _name.c_str(),is.fileName(),
	  mean->name().c_str(),
	  mean->dim(),
	  covar->name().c_str(),
	  covar->dim());
  }
  if ((unsigned)covar->dim() != dim()) {
    error("Error: LinMeanCondDiagGaussian '%s' in file '%s' of dim %d does not match its mean '%s' with dim %d or covariance '%s' with dim '%d'\n",
	  _name.c_str(),is.fileName(),
	  _dim,
	  mean->name().c_str(),
	  mean->dim(),
	  covar->name().c_str(),
	  covar->dim());
  }


  if ((unsigned)dLinks->dim() != _dim)
    error("Error: LinMeanCondDiagGaussian '%s' in file '%s' specifies a dlink structure '%s'(matrix '%s') that does not match its mean and covariance with dim '%d'\n",
	  _name.c_str(),is.fileName(),
	  dLinks->name().c_str(),
	  dLinkMat->name().c_str(),
	  _dim);

  setBasicAllocatedBit();
}


void
LinMeanCondDiagGaussian::write(oDataStreamFile& os)
{
  assert ( basicAllocatedBitIsSet() );

  // write the type of self and the name
  os.write((int)GaussianComponent::Diag);
  NamedObject::write(os);
  os.nl();

  // write mean vector
  os.write(mean->name()); 
  os.write(covar->name());
  os.write(dLinks->name());
  os.write(dLinkMat->name());
  os.nl();
}


void
LinMeanCondDiagGaussian::makeRandom()
{
  assert ( basicAllocatedBitIsSet() );

  mean->makeRandom();
  covar->makeRandom();
  dLinkMat->makeRandom();
}


void
LinMeanCondDiagGaussian::makeUniform()
{
  assert ( basicAllocatedBitIsSet() );

  mean->makeUniform();
  covar->makeUniform();
  dLinkMat->makeUniform();  
}


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
    float u=0.0;
    const int nComs = dLinks->numLinks(i);
    if (nComs > 0) {
      const int *lagStrideOffsets_endp = lagStrideOffsetsp+nComs;
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
LinMeanCondDiagGaussian::emStartIteration()
{
  assert ( basicAllocatedBitIsSet() );

  if (!GM_Parms.amTrainingDiagGaussians())
    return;

  if (emOnGoingBitIsSet())
    return;

  if (!emEmAllocatedBitIsSet()) {
    // this is presumably the first time
    emSetEmAllocatedBit();
  }

  // EM iteration is now going.
  emSetOnGoingBit();
  emSetSwappableBit();

  accumulatedProbability = 0.0;
  mean->emStartIteration(nextMeans);
  // dLinkMat->emStartIteration(nextDiagCovars);
  covar->emStartIteration(nextDiagCovars);
}


void
LinMeanCondDiagGaussian::emIncrement(logpr prob,
			  const float*f,
			  const Data32* const base,
			  const int stride)
{

  assert ( basicAllocatedBitIsSet() );

  if (!GM_Parms.amTrainingDiagGaussians())
    return;

  emStartIteration();

  if (prob < minIncrementProbabilty) {
    missedIncrementCount++;
    return;
    // don't accumulate anything since this one is so small and
    // if we did an unlog() and converted to a single precision
    // floating point number, it would be a denomral.
  }

  accumulatedProbability += prob;
  // prob.unlog() here so it doesn't need to be done
  // twice by the callees.
  const float fprob = prob.unlog();
  mean->emIncrement(prob,fprob,f,base,stride,nextMeans.ptr);
  // dLinkMat->emIncrement(prob,fprob,f,base,stride,
  // dLinks,nextDlinkMat);
  covar->emIncrement(prob,fprob,f,base,stride,nextDiagCovars.ptr);

}


void
LinMeanCondDiagGaussian::emEndIteration()
{
  assert ( basicAllocatedBitIsSet() );

  if (!GM_Parms.amTrainingDiagGaussians())
    return;

  if (!emOnGoingBitIsSet())
    return;

  if (accumulatedProbability == 0.0) {
    // TODO: need to check if this will overflow here
    // when dividing by it. This is more than just checking
    // for zero. Also need to do this in every such EM object.
    warning("WARNING: Gaussian Component named '%s' did not receive any accumulated probability in EM iteration",name().c_str());
  }

  mean->emEndIteration(nextMeans.ptr);
  covar->emEndIteration(accumulatedProbability,nextMeans.ptr,nextDiagCovars.ptr);

  emClearOnGoingBit();
}


void
LinMeanCondDiagGaussian::emSwapCurAndNew()
{
  assert ( basicAllocatedBitIsSet() );

  if (!GM_Parms.amTrainingDiagGaussians())
    return;

  if (!emSwappableBitIsSet())
    return;

  mean->emSwapCurAndNew();
  covar->emSwapCurAndNew();

  emClearSwappableBit();
}


void
LinMeanCondDiagGaussian::emStoreAccumulators(oDataStreamFile& ofile)
{
  assert ( basicAllocatedBitIsSet() );
  error("not implemented");
}

void
LinMeanCondDiagGaussian::emLoadAccumulators(iDataStreamFile& ifile)
{
  assert ( basicAllocatedBitIsSet() );
  error("not implemented");
}


void
LinMeanCondDiagGaussian::emAccumulateAccumulators(iDataStreamFile& ifile)
{
  assert ( basicAllocatedBitIsSet() );
  error("not implemented");
}


////////////////////////////////////////////////////////////
// Sample generation
////////////////////////////////////////////////////////////



void LinMeanCondDiagGaussian::sampleGenerate(float *const sample,
				  const Data32* const base,
				  const int stride)
{
  error("not implemented");
}









