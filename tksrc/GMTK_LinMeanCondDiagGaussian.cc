/*-
 * 
 * GMTK_LinMeanCondDiagGaussian.cc
 * 
 *        Code for linear mean conditioanl diagonal Gaussians.
 *        This code also supports
 *             - linear BMMs (linear buried Markov models)
 *                 (with dependencies coming from the past and/or the future)
 *             - full covariance matrices
 *             - semi-tied shared factored sparse inverse covariance matrices
 *             - covariances viewed as a directed graphical model
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
#include "lineqsolve.h"

#include "GMTK_LinMeanCondDiagGaussian.h"
#include "GMTK_GMParms.h"
#include "GMTK_MixGaussiansCommon.h"


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

  // read the dlink matrix parameter values
  is.read(str);
  if (GM_Parms.dLinkMatsMap.find(str) == GM_Parms.dLinkMatsMap.end())
    error("Error: LinMeanCondDiagGaussian '%s' in file '%s' specifies dlink matrix name '%s' that does not exist",
	  _name.c_str(),is.fileName(),str.c_str());
  
  dLinkMat = GM_Parms.dLinkMats[GM_Parms.dLinkMatsMap[str]];

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


  if ((unsigned)dLinkMat->dim() != _dim)
    error("Error: LinMeanCondDiagGaussian '%s' in file '%s' specifies a dlink matrix '%s' that does not match its mean and covariance having dim '%d'\n",
	  _name.c_str(),is.fileName(),
	  dLinkMat->name().c_str(),
	  _dim);

  setBasicAllocatedBit();
}


void
LinMeanCondDiagGaussian::write(oDataStreamFile& os)
{
  assert ( basicAllocatedBitIsSet() );

  // write the type of self and the name
  os.write((int)GaussianComponent::LinMeanCondDiag);
  NamedObject::write(os);
  os.nl();

  // write mean vector
  os.write(mean->name()); 
  os.write(covar->name());
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
  Dlinks* const dLinks = dLinkMat->dLinks;


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
    const int nLinks = dLinks->numLinks(i);
    if (nLinks > 0) {
      const int *lagStrideOffsets_endp = lagStrideOffsetsp+nLinks;
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
LinMeanCondDiagGaussian::noisyClone()
{
  assert ( basicAllocatedBitIsSet() );

  map<GaussianComponent*,GaussianComponent*>::iterator it = MixGaussiansCommon::gcCloneMap.find(this);
  // first check if self is already cloned, and if so, return that.
  if (it == MixGaussiansCommon::gcCloneMap.end()) {
    LinMeanCondDiagGaussian* clone;
    clone = new LinMeanCondDiagGaussian(dim());

    unsigned cloneNo=0; do {
      char buff[256];
      sprintf(buff,"%d",cloneNo);
      clone->_name = _name + string("_cl") + buff;
      cloneNo++;
    } while (GM_Parms.gaussianComponentsMap.find(clone->_name) 
	     != GM_Parms.gaussianComponentsMap.end());

    clone->mean = mean->noisyClone();
    clone->covar = covar->noisyClone();
    clone->dLinkMat = dLinkMat->noisyClone();

    clone->setBasicAllocatedBit();
    MixGaussiansCommon::gcCloneMap[this] = clone;

    // also add self to GMParms object.
    GM_Parms.add(clone);

    return clone;
  } else {
    return (*it).second;
  }
}


/////////////////
// EM routines //
/////////////////



void
LinMeanCondDiagGaussian::emStartIteration()
{
  assert ( basicAllocatedBitIsSet() );
  if (!GM_Parms.amTrainingLinMeanCondDiagGaussians())
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
  mean->emStartIteration(xAccumulators);
  dLinkMat->emStartIteration(xzAccumulators,zzAccumulators,zAccumulators);
  covar->emStartIteration(xxAccumulators);
}


void
LinMeanCondDiagGaussian::emIncrement(logpr prob,
				     const float*f,
				     const Data32* const base,
				     const int stride)
{
  
  assert ( basicAllocatedBitIsSet() );
  if (!GM_Parms.amTrainingLinMeanCondDiagGaussians())
    return;

  emStartIteration();

  if (prob < minIncrementProbabilty) {
    missedIncrementCount++;
    return;
    // don't accumulate anything since this one is so small and
    // if we did an unlog() and converted to a single precision
    // floating point number, it might be a denomral.
  }

  accumulatedProbability += prob;
  // prob.unlog() here so it doesn't need to be done
  // multiple times by the callees.
  const float fprob = prob.unlog();
  mean->emIncrement(prob,fprob,f,base,stride,xAccumulators.ptr);
  dLinkMat->emIncrement(prob,fprob,f,base,stride,
			xzAccumulators.ptr,
			zzAccumulators.ptr,
			zAccumulators.ptr);
  covar->emIncrement(prob,fprob,f,base,stride,xxAccumulators.ptr);

}




/*-
 *-----------------------------------------------------------------------
 * emEndIteration	
 *      This routine ends the EM iteration. This routine "merges"
 *      together the means and dlink coefficients into a single matrix because 
 *      it is in that domain where we must compute a matrix inverse.
 * 
 * Preconditions:
 *      basic structures must be allocated.
 *
 * Postconditions:
 *      next B matrix (consisting of burying coefficients and means) 
 *      and covariances are computed.
 *
 * Side Effects:
 *      allocates and then frees memory from the heap.
 *
 * Results:
 *      nil
 *
 *-----------------------------------------------------------------------
 */
void
LinMeanCondDiagGaussian::emEndIteration()
{
  assert ( basicAllocatedBitIsSet() );
  if (!GM_Parms.amTrainingLinMeanCondDiagGaussians())
    return;


  if (!emOnGoingBitIsSet())
    return;

  accumulatedProbability.floor();
  if (accumulatedProbability < minContAccumulatedProbability()) {
    error("ERROR: Gaussian Component named '%s' received only %e accumulated log probability (min is %e) in EM iteration, also check child mean '%s', covar '%s', and dlink matrix '%s'",name().c_str(),
	  accumulatedProbability.val(),
	  minContAccumulatedProbability().val(),
	  mean->name().c_str(),
	  covar->name().c_str(),
	  dLinkMat->name().c_str());
  }

  const double realAccumulatedProbability =
    accumulatedProbability.unlog();

  // Now we need the fully expanded forms of 
  // xz, zz, and B. These are used because we need to
  // compute the inverse matrix in the formula B = (XZ)(ZZ)^(-1).
  // It is assumed that since this routine is not called often (relative
  // to emIncrement), we can do all the memory allocation/reclaimation here.
  
  /////////////////////////////////////////////////////////
  // xzExpAccumulators (expanded accumulators) contains the same
  // information as xzAccumulators but includes the accumulated mean
  // xAccumulators in the right-most position of each vector.
  sArray<double> xzExpAccumulators;

  xzExpAccumulators.resize( dLinkMat->dLinks->totalNumberLinks() + mean->dim() );

  // now copy xz and x over to expanded xz
  double *xzExpAccumulators_p = xzExpAccumulators.ptr;
  float *xzAccumulators_p = xzAccumulators.ptr;
  for (int feat=0;feat<mean->dim();feat++) {
    const int nLinks = dLinkMat->numLinks(feat);

    for (int dlink=0;dlink<nLinks;dlink++) {
      *xzExpAccumulators_p++ =
	*xzAccumulators_p++;
    }

    *xzExpAccumulators_p++ =
      xAccumulators[feat];
  }

  ///////////////////////////////////////////////////////////////////////////
  // zzExpAccumulators contains the same information as
  // zzAccumulators, but includes the extra right most column and
  // bottom most row containing the value accumulatedProbability.
  sArray<double> zzExpAccumulators;  
  // resize matrix, going from a size of just upper triangular
  // representation to one where we include the extra z=1 variable at
  // the end and represent a full matrix (needed for computing matrix
  // inverse below). 7/2/01 notes.
  zzExpAccumulators.resize (
			    2*dLinkMat->zzAccumulatorLength()
			    + dLinkMat->totalNumberLinks()
			    + mean->dim()
			    );
  // now expand out. This code does two things simultaneously:
  // 1) it turns the triangular matrix into a full matrix, and
  // 2) it adds an extra row&col to the full matrix for the final z value
  // which wasn't represented during the EM increment stage.
  float *zzAccumulators_p = zzAccumulators.ptr;
  float *zAccumulators_p = zAccumulators.ptr;
  double *zzExpAccumulators_p = zzExpAccumulators.ptr;
  for (int feat=0;feat<mean->dim();feat++) {
    const int nLinks = dLinkMat->numLinks(feat);
    
    float *zzp = zzAccumulators_p; // ptr to current zz
    double *zzep = zzExpAccumulators_p; // ptr to current exp zz

    for (int dlink=0;dlink<=nLinks;dlink++) {
      double *zze_rp = zzep; // row ptr to expanded zz
      double *zze_cp = zzep; // col ptr to expanded zz
      for (int j=0;j<(nLinks-dlink);j++) {
	*zze_rp = *zze_cp = *zzp++;
	zze_rp ++;            // increment by one value
	zze_cp += (nLinks+1); // increment by stride
      }
      if (dlink < nLinks) 
	*zze_rp = *zze_cp = *zAccumulators_p++;
      else 
	*zze_rp = *zze_cp = 
	  realAccumulatedProbability;

      zzep += (nLinks+2);
    }

    zzAccumulators_p += nLinks*(nLinks+1)/2;
    zzExpAccumulators_p += (nLinks+1)*(nLinks+1);

    // sanity assertions
    assert ( zzp == zzAccumulators_p );
    assert ( (zzep - nLinks - 1) == zzExpAccumulators_p );
  }

  ////////////////////////////////////////////////////////////////////////////
  // Now go through and compute the next dlink coefficients (which
  // includes the coefficients AND the means which will be contained
  // in the right most position of the sparse array.)
  sArray<double> nextDlinkMat;
  nextDlinkMat.resize( dLinkMat->dLinks->totalNumberLinks() + mean->dim() );
  double *nextDlinkMat_p =  nextDlinkMat.ptr;
  zzExpAccumulators_p = zzExpAccumulators.ptr;
  xzExpAccumulators_p = xzExpAccumulators.ptr;
  for (int feat=0;feat<mean->dim();feat++) {
    const int nLinks = dLinkMat->numLinks(feat);

    // Solve for the burying coefficients using
    // a routine which solves Ax = b for x where
    // A is nXn, x is nX1, and b is nX1
    // here,
    //     A = zzExpAccumulators_p,
    //     x = the output
    //     b = xzExpAccumulators_p (which gets destroyed)

    // First, copy xzAccumulators over to nextDlinkMat since
    // we will need the values of xzAccumulators later.
    ::memcpy(nextDlinkMat_p,xzExpAccumulators_p,sizeof(double)*(nLinks+1));
    // Finally, solve for the link values putting the results
    // in nextDlinkMat_p (the current values get destroyed).
    ::lineqsolve(nLinks+1,1,
		 zzExpAccumulators_p,nextDlinkMat_p);

    // now solve for the variances
    double tmp = 0.0;
    for (int i=0;i<(nLinks+1);i++) {
      tmp += ( nextDlinkMat_p[i] * xzExpAccumulators_p[i]);
    }

    // finally solve for the variance, putting the result back in
    // the xx accumulator. Do *NOT* normalize by the posterior
    // accumulator here as that will be done by the covariance
    // object when we give this to it, below. Note also that
    // we do not check for variances being too small here, that
    // is done below as well, since when sharing occurs, the
    // individual covars might be small, but the shared covariances
    // might be fine.
    xxAccumulators[feat] = 
      (xxAccumulators[feat] - tmp);

    xzExpAccumulators_p += (nLinks+1);
    nextDlinkMat_p += (nLinks+1);
    zzExpAccumulators_p += (nLinks+1)*(nLinks+1);

  }

  // finally, copy out the means and dlinks in 
  // *NON* normalized form. Renormalization will occur
  // (with the (shared if any) sum posteriors) by the objects
  // that handle those guys, below.

  // Store the next means in xAccumulators and
  // store the dlinks in xzAccumulators.

  // now copy it and means over.
  nextDlinkMat_p = nextDlinkMat.ptr;
  xzAccumulators_p = xzAccumulators.ptr;
  for (int feat=0;feat<mean->dim();feat++) {
    const int nLinks = dLinkMat->numLinks(feat);
    for (int dlink=0;dlink<nLinks;dlink++) {
      *xzAccumulators_p++ = // unnormalize
	(*nextDlinkMat_p++) * realAccumulatedProbability; 
    }
    xAccumulators[feat] = // unnormalize
      (*nextDlinkMat_p++) * realAccumulatedProbability; 
  }

  // Finally, incorporate our hard work into the 
  // (respectively) shared objects who will do any
  // normalization
  mean->emEndIteration(xAccumulators.ptr);
  dLinkMat->emEndIteration(xzAccumulators.ptr);
  covar->emEndIteration(xxAccumulators.ptr);

  // Finally, end the EM epoch.
  emClearOnGoingBit();
}


void
LinMeanCondDiagGaussian::emSwapCurAndNew()
{
  assert ( basicAllocatedBitIsSet() );
  if (!GM_Parms.amTrainingLinMeanCondDiagGaussians())
    return;

  if (!emSwappableBitIsSet())
    return;

  mean->emSwapCurAndNew();
  dLinkMat->emSwapCurAndNew();
  covar->emSwapCurAndNew();

  emClearSwappableBit();
}


////////////////////////////////////////////////////////////
// Parallel EM support
////////////////////////////////////////////////////////////


void
LinMeanCondDiagGaussian::emStoreAccumulators(oDataStreamFile& ofile)
{
  assert ( basicAllocatedBitIsSet() );
  if ( !emEmAllocatedBitIsSet() ) {
    warning("WARNING: storing zero accumulators for lin mean cond Gaussian '%s'\n",
	    name().c_str());
    emStoreZeroAccumulators(ofile);
    return;
  }
  EMable::emStoreAccumulators(ofile);
  for (int i=0;i<xAccumulators.len();i++) {
    ofile.write(xAccumulators[i],"LMDG store accums x.");
  }
  for (int i=0;i<xxAccumulators.len();i++) {
    ofile.write(xxAccumulators[i],"LMDG store accums xx.");
  }
  for (int i=0;i<xzAccumulators.len();i++) {
    ofile.write(xzAccumulators[i],"LMDG store accums xz.");
  }
  for (int i=0;i<zzAccumulators.len();i++) {
    ofile.write(zzAccumulators[i],"LMDG store accums zz.");
  }
  for (int i=0;i<zAccumulators.len();i++) {
    ofile.write(zAccumulators[i],"LMDG store accums z.");
  }
}



void
LinMeanCondDiagGaussian::emStoreZeroAccumulators(oDataStreamFile& ofile)
{
  assert ( basicAllocatedBitIsSet() );
  EMable::emStoreZeroAccumulators(ofile);
  for (int i=0;i<mean->dim();i++) {
    ofile.write((float)0.0,"LMDG store zero accums x.");
  }
  for (int i=0;i<covar->dim();i++) {
    ofile.write((float)0.0,"LMDG store zero accums xx.");
  }
  for (int i=0;i<(int)dLinkMat->totalNumberLinks();i++) {
    ofile.write((float)0.0,"LMDG store zero accums xz.");
  }
  for (int i=0;i<(int)dLinkMat->zzAccumulatorLength();i++) {
    ofile.write((float)0.0,"LMDG store zero accums zz.");
  }
  for (int i=0;i<covar->dim();i++) {
    ofile.write((float)0.0,"LMDG store zero accums z.");
  }
}


void
LinMeanCondDiagGaussian::emLoadAccumulators(iDataStreamFile& ifile)
{
  assert ( basicAllocatedBitIsSet() );
  assert ( emEmAllocatedBitIsSet() );
  EMable::emLoadAccumulators(ifile);
  for (int i=0;i<xAccumulators.len();i++) {
    ifile.read(xAccumulators[i],"LMDG load accums x.");
  }
  for (int i=0;i<xxAccumulators.len();i++) {
    ifile.read(xxAccumulators[i],"LMDG load accums xx.");
  }
  for (int i=0;i<xzAccumulators.len();i++) {
    ifile.read(xzAccumulators[i],"LMDG load accums xz.");
  }
  for (int i=0;i<zzAccumulators.len();i++) {
    ifile.read(zzAccumulators[i],"LMDG load accums zz.");
  }
  for (int i=0;i<zAccumulators.len();i++) {
    ifile.read(zAccumulators[i],"LMDG load accums z.");
  }

}


void
LinMeanCondDiagGaussian::emAccumulateAccumulators(iDataStreamFile& ifile)
{
  assert ( basicAllocatedBitIsSet() );
  assert ( emEmAllocatedBitIsSet() );
  EMable::emAccumulateAccumulators(ifile);
  for (int i=0;i<xAccumulators.len();i++) {
    float tmp;
    ifile.read(tmp,"LMDG accumulate accums x.");
    xAccumulators[i] += tmp;
  }
  for (int i=0;i<xxAccumulators.len();i++) {
    float tmp;
    ifile.read(tmp,"LMDG accumulate accums xx.");
    xxAccumulators[i] += tmp;
  }
  for (int i=0;i<xzAccumulators.len();i++) {
    float tmp;
    ifile.read(tmp,"LMDG accumulate accums xz.");
    xzAccumulators[i] += tmp;
  }
  for (int i=0;i<zzAccumulators.len();i++) {
    float tmp;
    ifile.read(tmp,"LMDG accumulate accums zz.");
    zzAccumulators[i] += tmp;
  }
  for (int i=0;i<zAccumulators.len();i++) {
    float tmp;
    ifile.read(tmp,"LMDG accumulate accums z.");
    zAccumulators[i] += tmp;
  }
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









