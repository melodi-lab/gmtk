/*-
 * GMTK_DiagGaussian.cc
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
#include <float.h>
#include <assert.h>
#include <ctype.h>

#include <string>

#include "general.h"
VCID("$Header$");
#include "error.h"
#include "rand.h"

#include "GMTK_DiagGaussian.h"
#include "GMTK_GMParms.h"
#include "GMTK_MixGaussiansCommon.h"



void
DiagGaussian::read(iDataStreamFile& is)
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
      error("Error: DiagGaussian '%s' in file '%s' specifies mean name '%s' that does not exist",
	    _name.c_str(),is.fileName(),str.c_str());
  meanIndex = GM_Parms.meansMap[str];
  mean = GM_Parms.means[meanIndex];


  // read covariance vector
  is.read(str);
  if (GM_Parms.covarsMap.find(str) == GM_Parms.covarsMap.end())
    error("Error: DiagGaussian '%s' in file '%s' specifies covar name '%s' that does not exist",
	  _name.c_str(),is.fileName(),str.c_str());
  
  covarIndex = GM_Parms.covarsMap[str];
  covar = GM_Parms.covars[covarIndex];

  // check that lengths match, etc.
  if (covar->dim() != mean->dim()) {
    error("Error: LinMeanCondDiagGaussian '%s' in file '%s' specifices a mean '%s' with dim %d and covariance '%s' with dim '%d'\n",
	  _name.c_str(),is.fileName(),
	  mean->name().c_str(),
	  mean->dim(),
	  covar->name().c_str(),
	  covar->dim());
  }
  if ((unsigned)covar->dim() != _dim) {
    error("Error: LinMeanCondDiagGaussian '%s' in file '%s' of dim %d does not match its mean '%s' with dim %d or covariance '%s' with dim '%d'\n",
	  _name.c_str(),is.fileName(),
	  _dim,
	  mean->name().c_str(),
	  mean->dim(),
	  covar->name().c_str(),
	  covar->dim());
  }

  setBasicAllocatedBit();
}


void
DiagGaussian::write(oDataStreamFile& os)
{
  assert ( basicAllocatedBitIsSet() );

  // write the type of self and the name
  os.write((int)GaussianComponent::Diag);
  NamedObject::write(os);
  os.nl();

  // write mean vector
  os.write(mean->name());
  os.write(covar->name());
  os.nl();
}


void
DiagGaussian::makeRandom()
{
  assert ( basicAllocatedBitIsSet() );

  mean->makeRandom();
  covar->makeRandom();
}


void
DiagGaussian::makeUniform()
{
  assert ( basicAllocatedBitIsSet() );

  mean->makeUniform();
  covar->makeUniform();
}


/*-
 *-----------------------------------------------------------------------
 * noisyClone()
 *      Create a copy of self, but with cloned perturbed the mean/variance vectors
 * 
 * Preconditions:
 *      Obj must be read in, and basicAllocatedBitIsSet() must be true.
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
GaussianComponent*
DiagGaussian::noisyClone()
{
  assert ( basicAllocatedBitIsSet() );

  map<GaussianComponent*,GaussianComponent*>::iterator it = MixGaussiansCommon::gcCloneMap.find(this);
  // first check if self is already cloned, and if so, return that.
  if (it == MixGaussiansCommon::gcCloneMap.end()) {
    DiagGaussian* clone;
    clone = new DiagGaussian(dim());

    unsigned cloneNo=0; do {
      char buff[256];
      sprintf(buff,"%d",cloneNo);
      clone->_name = _name + string("_cl") + buff;
      cloneNo++;
    } while (GM_Parms.gaussianComponentsMap.find(clone->_name) 
	     != GM_Parms.gaussianComponentsMap.end());

    clone->mean = mean->noisyClone();
    clone->covar = covar->noisyClone();

    clone->setBasicAllocatedBit();
    MixGaussiansCommon::gcCloneMap[this] = clone;

    // also add self to GMParms object.
    GM_Parms.add(clone);

    return clone;
  } else {
    return (*it).second;
  }
}




/*-
 *-----------------------------------------------------------------------
 * log_p()
 *      Computes the probability of this Gaussian.
 * 
 * Preconditions:
 *      preCompute() must have been called on covariance matrix before this.
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
DiagGaussian::log_p(const float *const x,
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
  DIAG_GAUSSIAN_TMP_ACCUMULATOR_TYPE d=0.0;
  do {
    const DIAG_GAUSSIAN_TMP_ACCUMULATOR_TYPE tmp
      = (*xp - *mean_p);
    d += tmp*tmp*(*var_inv_p);

    xp++;
    mean_p++;
    var_inv_p++;
  } while (xp != x_endp);
  d *= -0.5;
  return logpr(0,(covar->log_inv_normConst() + d));

}



/////////////////
// EM routines //
/////////////////



void
DiagGaussian::emStartIteration()
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
  covar->emStartIteration(nextDiagCovars);
}


void
DiagGaussian::emIncrement(logpr prob,
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
  covar->emIncrement(prob,fprob,f,base,stride,nextDiagCovars.ptr);
}


void
DiagGaussian::emEndIteration()
{
  assert ( basicAllocatedBitIsSet() );
  if (!GM_Parms.amTrainingDiagGaussians())
    return;

  // if EM not ongoing, we just return
  // here since, if this object is shared,
  // we don't want to end the epoch > 1 time
  if (!emOnGoingBitIsSet())
    return; 

  if (accumulatedProbability == 0.0) {
    // Note: we assume here that the mean and covar object will check for us that 
    // its accumulated probability is above threshold. Here, we just
    // check for zero, and then issue a warning if needed. This might not
    // indicate a problem as the mean and covar of this object might be shared
    // and might have received plenty of count.
    warning("WARNING: Gaussian Component named '%s' did not receive any accumulated probability in EM iteration, check child mean '%s' and covar '%s'",name().c_str(),mean->name().c_str(),covar->name().c_str());
  }

  mean->emEndIteration(nextMeans.ptr);
  covar->emEndIteration(accumulatedProbability,nextMeans.ptr,nextDiagCovars.ptr);

  emClearOnGoingBit();
}


void
DiagGaussian::emSwapCurAndNew()
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


////////////////////////////////////////////////////////////
// Parallel EM support
////////////////////////////////////////////////////////////

void
DiagGaussian::emStoreAccumulators(oDataStreamFile& ofile)
{
  assert ( basicAllocatedBitIsSet() );
  assert ( emEmAllocatedBitIsSet() );
  EMable::emStoreAccumulators(ofile);
  for (int i=0;i<nextMeans.len();i++) {
    ofile.write(nextMeans[i],"nxm");
  }
  for (int i=0;i<nextDiagCovars.len();i++) {
    ofile.write(nextDiagCovars[i],"ndc");
  }
}

void
DiagGaussian::emLoadAccumulators(iDataStreamFile& ifile)
{
  assert ( basicAllocatedBitIsSet() );
  assert ( emEmAllocatedBitIsSet() );
  EMable::emLoadAccumulators(ifile);
  for (int i=0;i<nextMeans.len();i++) {
    ifile.read(nextMeans[i],"nxm");
  }
  for (int i=0;i<nextDiagCovars.len();i++) {
    ifile.read(nextDiagCovars[i],"ndc");
  }
}


void
DiagGaussian::emAccumulateAccumulators(iDataStreamFile& ifile)
{
  assert ( basicAllocatedBitIsSet() );
  assert ( emEmAllocatedBitIsSet() );
  EMable::emAccumulateAccumulators(ifile);
  for (int i=0;i<nextMeans.len();i++) {
    float tmp;
    ifile.read(tmp,"nxm");
    nextMeans[i] += tmp;
  }
  for (int i=0;i<nextDiagCovars.len();i++) {
    float tmp;
    ifile.read(tmp,"ndc");
    nextDiagCovars[i] += tmp;
  }
}


////////////////////////////////////////////////////////////
// Sample generation
////////////////////////////////////////////////////////////



void DiagGaussian::sampleGenerate(float *const sample,
				  const Data32* const base,
				  const int stride)
{
  error("not implemented");
}









