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
VCID("$Header$")
#include "error.h"
#include "rand.h"
#include "lineqsolve.h"

#include "GMTK_LinMeanCondDiagGaussian.h"
#include "GMTK_GMParms.h"
#include "GMTK_MixtureCommon.h"
#include "GMTK_DiagGaussian.h"

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
      error("Error: LinMeanCondDiagGaussian '%s' in file '%s' line %d, specifies mean name '%s' that does not exist",
	    _name.c_str(),is.fileName(),is.lineNo(),str.c_str());
  meanIndex = GM_Parms.meansMap[str];
  mean = GM_Parms.means[meanIndex];
  mean->numTimesShared++;


  // read covariance vector
  is.read(str);
  if (GM_Parms.covarsMap.find(str) == GM_Parms.covarsMap.end())
    error("Error: LinMeanCondDiagGaussian '%s' in file '%s' line %d, specifies covar name '%s' that does not exist",
	  _name.c_str(),is.fileName(),is.lineNo(),str.c_str());
  
  covarIndex = GM_Parms.covarsMap[str];
  covar = GM_Parms.covars[covarIndex];
  covar->numTimesShared++;

  // read the dlink matrix parameter values
  is.read(str);
  if (GM_Parms.dLinkMatsMap.find(str) == GM_Parms.dLinkMatsMap.end())
    error("Error: LinMeanCondDiagGaussian '%s' in file '%s' line %d, specifies dlink matrix name '%s' that does not exist",
	  _name.c_str(),is.fileName(),is.lineNo(),str.c_str());
  
  dLinkMat = GM_Parms.dLinkMats[GM_Parms.dLinkMatsMap[str]];
  dLinkMat->numTimesShared++;


  // check that lengths match, etc.
  if (covar->dim() != mean->dim()) {
    error("Error: LinMeanCondDiagGaussian '%s' in file '%s' line %d, specifices a mean '%s' with dim %d and covariance '%s' with dim '%d'\n",
	  _name.c_str(),is.fileName(),is.lineNo(),
	  mean->name().c_str(),
	  mean->dim(),
	  covar->name().c_str(),
	  covar->dim());
  }
  if ((unsigned)covar->dim() != dim()) {
    error("Error: LinMeanCondDiagGaussian '%s' in file '%s' line %d of dim %d does not match its mean '%s' with dim %d or covariance '%s' with dim '%d'\n",
	  _name.c_str(),
	  is.fileName(),is.lineNo(),
	  _dim,
	  mean->name().c_str(),
	  mean->dim(),
	  covar->name().c_str(),
	  covar->dim());
  }


  if ((unsigned)dLinkMat->dim() != _dim)
    error("Error: LinMeanCondDiagGaussian '%s' in file '%s' line %d specifies a dlink matrix '%s' that does not match its mean and covariance having dim '%d'\n",
	  _name.c_str(),
	  is.fileName(),is.lineNo(),
	  dLinkMat->name().c_str(),
	  _dim);

  setBasicAllocatedBit();
}


void
LinMeanCondDiagGaussian::write(oDataStreamFile& os)
{
  assert ( basicAllocatedBitIsSet() );

  // write the type of self and the name
  os.write((int)Component::LinMeanCondDiagGaussian);
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
  if (!emAmTrainingBitIsSet())
    return;

  mean->makeRandom();
  covar->makeRandom();
  dLinkMat->makeRandom();
}


void
LinMeanCondDiagGaussian::makeUniform()
{
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;

  mean->makeUniform();
  covar->makeUniform();
  dLinkMat->makeUniform();  
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
Component*
LinMeanCondDiagGaussian::noisyClone()
{
  assert ( basicAllocatedBitIsSet() );

  map<Component*,Component*>::iterator it = MixtureCommon::mcCloneMap.find(this);
  // first check if self is already cloned, and if so, return that.
  if (it == MixtureCommon::mcCloneMap.end()) {
    LinMeanCondDiagGaussian* clone;
    clone = new LinMeanCondDiagGaussian(dim());

    unsigned cloneNo=0; do {
      char buff[256];
      sprintf(buff,"%d",cloneNo);
      clone->_name = _name + string("_cl") + buff;
      cloneNo++;
    } while (GM_Parms.componentsMap.find(clone->_name) 
	     != GM_Parms.componentsMap.end());

    if (cloneShareMeans && cloneShareCovars && cloneShareDlinks) {
      warning("WARNING: Dlink Gaussian component '%s' is cloning, and was asked to share both means, covariances, and dlinks. No sharing is occuring instead.",name().c_str());
      clone->mean = mean->noisyClone();
      clone->covar = covar->noisyClone();
      clone->dLinkMat = dLinkMat->noisyClone();
    } else {
      if (cloneShareMeans)
	clone->mean = mean;
      else 
	clone->mean = mean->noisyClone();
      if (cloneShareCovars)
	clone->covar = covar;
      else
	clone->covar = covar->noisyClone();
      if (cloneShareDlinks)
	clone->dLinkMat = dLinkMat;
      else
	clone->dLinkMat = dLinkMat->noisyClone();
    }
    
    // need to tell mean, covar, and dlink that either
    //    1) if this is a new mean,covar,dlink, that a
    //       a parent object is using them, or
    //    2) if this is a shared mean,covar,dlink, that
    //       an additional parent object is using them.
    clone->mean->numTimesShared++;
    clone->covar->numTimesShared++;
    clone->dLinkMat->numTimesShared++;

    clone->setBasicAllocatedBit();
    MixtureCommon::mcCloneMap[this] = clone;

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
  if (!emAmTrainingBitIsSet())
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
  if (!emAmTrainingBitIsSet())
    return;

  if (!emOnGoingBitIsSet())
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
  DiagGaussian::emIncrementMeanDiagCovar(fprob,f,xAccumulators.size(),xAccumulators.ptr,xxAccumulators.ptr);
}


/*-
 *-----------------------------------------------------------------------
 * emEndIteration
 *      This routine ends the EM iteration. It figures out the
 *      current sharing structure, and calls the appropriate
 *      routine (EM or GEM) accordingly.
 *
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
  assert (basicAllocatedBitIsSet());
  if (!emAmTrainingBitIsSet())
    return;

  if (!emOnGoingBitIsSet())
    return;

  const unsigned shareBits =
    (mean->emSharedBitIsSet()?0x1:0x0)
    |
    (covar->emSharedBitIsSet()?0x2:0x0)
    |
    (dLinkMat->emSharedBitIsSet()?0x4:0x0);


  /////////////////////////////////////////////
  // use a switch as we might in the future
  // want to implement other forms.
  switch (shareBits) {
  case 0x0:
    // no sharing at all
    {
      emEndIterationNoSharing();
    }
    break;
  case 0x2:
    // tied diagonal covariance matrices
    {
      emEndIterationSharedCovars();
    }
    break;
  default:
    // the most general case, a GEM
    {
      emEndIterationSharedAll();
    }
    break;
  }

  // Finally, end the EM epoch.
  emClearOnGoingBit();
}


/*-
 *-----------------------------------------------------------------------
 * emEndIterationNoSharing
 *      This routine ends the EM iteration. This routine "merges"
 *      together the means and dlink coefficients into a single matrix because 
 *      it is in that domain where we must compute a matrix inverse in the no sharing
 *      case.
 *
 *   NOTE: this routine is almost identical to emEndIterationSharedCovars();
 * 
 * Preconditions:
 *      - basic structures must be allocated.
 *      - Must be called only by emEndIteration()
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
LinMeanCondDiagGaussian::emEndIterationNoSharing()
{
  // accumulatedProbability.floor();
  if (accumulatedProbability < minContAccumulatedProbability()) {
    warning("WARNING: Lin Mean-Cond Diag Gaussian Component named '%s' received only %e accumulated log probability (min is %e) in EM iteration, Global missed increment count is %d. Also check child mean '%s', covar '%s', and dlink matrix '%s'",
	  name().c_str(),
	  accumulatedProbability.val(),
	  minContAccumulatedProbability().val(),
	  missedIncrementCount,
	  mean->name().c_str(),
	  covar->name().c_str(),
	  dLinkMat->name().c_str());
    //////////////////////////////////////////////////////////
    // Since the probability is so small, it is likely
    // that the accumulators are tiny or zero anyway. We
    // pass them on in this form to the child object accumulators
    // (mean, covar, dlinkmat), which will increment them as is (but it
    // shouldn't do much since they are zero). We expect that the child objects
    // should notice that their accumulated probability is small
    // and hence use the previous iterations parameter values.
    goto finishup;
  }

  /////////////////////////////////////////////////////////
  // make this a {} block so that we can jump over it
  // above in the goto.
  { 
    const double realAccumulatedProbability =
      accumulatedProbability.unlog();
    const double invRealAccumulatedProbability =
      accumulatedProbability.inverse().unlog();


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
    // inverse below). See 7/2/01 bilmes hand-written notes.
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
    // includes the coefficients AND the means, the later
    // of which which will be contained in the right most position of the 
    // sparse array.)
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
      // the xx accumulator. Normalize by the posterior
      // accumulator here as that will not be done by the covariance
      // object when we give this to it, below. Note also that
      // we do not check for variances being too small here, that
      // is done in the variance object itself.
      xxAccumulators[feat] = 
	(xxAccumulators[feat] - tmp)*invRealAccumulatedProbability;

      xzExpAccumulators_p += (nLinks+1);
      nextDlinkMat_p += (nLinks+1);
      zzExpAccumulators_p += (nLinks+1)*(nLinks+1);

    }

    // finally, copy out the means and dlinks in 
    // *ALREADY NORMALIZED* form. 

    // Store the next means in xAccumulators and
    // store the dlinks in xzAccumulators.

    // now copy it and means over.
    nextDlinkMat_p = nextDlinkMat.ptr;
    xzAccumulators_p = xzAccumulators.ptr;
    for (int feat=0;feat<mean->dim();feat++) {
      const int nLinks = dLinkMat->numLinks(feat);
      for (int dlink=0;dlink<nLinks;dlink++) {
	*xzAccumulators_p++ = *nextDlinkMat_p++;
      }
      xAccumulators[feat] = *nextDlinkMat_p++;
    }
  }

 finishup:
  // Finally, incorporate our hard work into the 
  // (respectively) shared objects who will do any
  // normalization
  mean->emEndIterationNoSharingAlreadyNormalized(xAccumulators.ptr);
  covar->emEndIterationNoSharingAlreadyNormalized(xxAccumulators.ptr);
  dLinkMat->emEndIterationNoSharingAlreadyNormalized(xzAccumulators.ptr);

  // Finally, end the EM epoch.
  emClearOnGoingBit();
}



/*-
 *-----------------------------------------------------------------------
 * emEndIterationSharedCovars
 *      This routine ends the EM iteration. This routine "merges"
 *      together the means and dlink coefficients into a single matrix because 
 *      it is in that domain where we must compute a matrix inverse. It
 *      is assumed that the covariance is shared.
 *
 *   NOTE: this routine is almost identical to emEndIterationNoSharing();
 *

 * 
 * Preconditions:
 *      - basic structures must be allocated.
 *      - Must be called only by emEndIteration()
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
LinMeanCondDiagGaussian::emEndIterationSharedCovars()
{
  // accumulatedProbability.floor();
  if (accumulatedProbability < minContAccumulatedProbability()) {
    warning("WARNING: Lin Mean-Cond Diag Gaussian Component named '%s' received only %e accumulated log probability (min is %e) in EM iteration, Global missed increment count is %d. Also check child mean '%s', covar '%s', and dlink matrix '%s'",
	  name().c_str(),
	  accumulatedProbability.val(),
	  minContAccumulatedProbability().val(),
	  missedIncrementCount,
	  mean->name().c_str(),
	  covar->name().c_str(),
	  dLinkMat->name().c_str());
    //////////////////////////////////////////////////////////
    // Since the probability is so small, it is likely
    // that the accumulators are tiny or zero anyway. We
    // pass them on in this form to the child object accumulators
    // (mean, covar, dlinkmat), which will increment them as is (but it
    // shouldn't do much since they are zero). We expect that the child objects
    // should notice that their accumulated probability is small
    // and hence use the previous iterations parameter values.
    goto finishup;
  }

  /////////////////////////////////////////////////////////
  // make this a {} block so that we can jump over it
  // above in the goto.
  { 
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
    // *ALREADY NORMALIZED* form. 

    // Store the next means in xAccumulators and
    // store the dlinks in xzAccumulators.

    // now copy it and means over.
    nextDlinkMat_p = nextDlinkMat.ptr;
    xzAccumulators_p = xzAccumulators.ptr;
    for (int feat=0;feat<mean->dim();feat++) {
      const int nLinks = dLinkMat->numLinks(feat);
      for (int dlink=0;dlink<nLinks;dlink++) {
	*xzAccumulators_p++ = *nextDlinkMat_p++;
      }
      xAccumulators[feat] = *nextDlinkMat_p++;
    }
  }

 finishup:
  // Finally, incorporate our hard work into the 
  // (respectively) shared objects who will do any
  // normalization
  mean->emEndIterationNoSharingAlreadyNormalized(xAccumulators.ptr);
  covar->emEndIterationSharedCovars(xxAccumulators.ptr);
  dLinkMat->emEndIterationNoSharingAlreadyNormalized(xzAccumulators.ptr);

  // Finally, end the EM epoch.
  emClearOnGoingBit();
}



/*-
 *-----------------------------------------------------------------------
 * emEndIterationSharedAll()
 *      This routine ends the EM iteration. This routine assumes
 *      that arbitrary sharing might be taking place and so uses the
 *      a GEM.
 *
 * 
 * Preconditions:
 *      - basic structures must be allocated.
 *      - Must be called only by emEndIteration()
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
LinMeanCondDiagGaussian::emEndIterationSharedAll()
{

  // accumulatedProbability.floor();
  if (accumulatedProbability < minContAccumulatedProbability()) {
    warning("WARNING: Lin Mean-Cond Diag Gaussian Component named '%s' received only %e accumulated log probability (min is %e) in EM iteration, Global missed increment count is %d. Also check child mean '%s', covar '%s', and dlink matrix '%s'",
	  name().c_str(),
	  accumulatedProbability.val(),
	  minContAccumulatedProbability().val(),
	  missedIncrementCount,
	  mean->name().c_str(),
	  covar->name().c_str(),
	  dLinkMat->name().c_str());
    //////////////////////////////////////////////////////////
    // Since the probability is so small, it is likely
    // that the accumulators are tiny or zero anyway. We
    // pass them on in this form to the child object accumulators
    // (mean, covar, dlinkmat), which will increment them as is (but it
    // shouldn't do much since they are zero). We expect that the child objects
    // should notice that their accumulated probability is small
    // and hence use the previous iterations parameter values.
  }

  mean->emEndIterationSharedMeansCovarsDlinks(accumulatedProbability,
					      xAccumulators.ptr,
					      zAccumulators.ptr,
					      dLinkMat,
					      covar);

  covar->emEndIterationSharedMeansCovarsDlinks(accumulatedProbability,
					       xAccumulators.ptr,
					       xxAccumulators.ptr,
					       xzAccumulators.ptr,
					       zAccumulators.ptr,
					       zzAccumulators.ptr,
					       mean,
					       dLinkMat);

  dLinkMat->emEndIterationSharedMeansCovarsDlinks(xzAccumulators.ptr,
						  zAccumulators.ptr,
						  zzAccumulators.ptr,
						  mean,
						  covar);
						  

  // Finally, end the EM epoch.
  emClearOnGoingBit();
}









void
LinMeanCondDiagGaussian::emSwapCurAndNew()
{
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
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
LinMeanCondDiagGaussian::emStoreObjectsAccumulators(oDataStreamFile& ofile)
{
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
LinMeanCondDiagGaussian::emLoadObjectsDummyAccumulators(iDataStreamFile& ifile)
{
  // ASSUME ACCUMULATOR TYPE IS 'float'
  float tmp;
  for (int i=0;i<mean->dim();i++) {
    ifile.read(tmp,"LMDG load accums x.");
  }
  for (int i=0;i<covar->dim();i++) {
    ifile.read(tmp,"LMDG load accums xx.");
  }
  for (int i=0;i<(int)dLinkMat->totalNumberLinks();i++) {
    ifile.read(tmp,"LMDG load accums xz.");
  }
  for (int i=0;i<(int)dLinkMat->zzAccumulatorLength();i++) {
    ifile.read(tmp,"LMDG load accums zz.");
  }
  for (int i=0;i<(int)dLinkMat->totalNumberLinks();i++) {
    ifile.read(tmp,"LMDG load accums z.");
  }

}


void
LinMeanCondDiagGaussian::emZeroOutObjectsAccumulators()
{
  for (int i=0;i<xAccumulators.len();i++) {
    xAccumulators[i] = 0.0;
  }
  for (int i=0;i<xxAccumulators.len();i++) {
    xxAccumulators[i] = 0.0;
  }
  for (int i=0;i<xzAccumulators.len();i++) {
    xzAccumulators[i] = 0.0;
  }
  for (int i=0;i<zzAccumulators.len();i++) {
    zzAccumulators[i] = 0.0;
  }
  for (int i=0;i<zAccumulators.len();i++) {
    zAccumulators[i] = 0.0;
  }
}


void
LinMeanCondDiagGaussian::emLoadObjectsAccumulators(iDataStreamFile& ifile)
{
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
LinMeanCondDiagGaussian::emAccumulateObjectsAccumulators(iDataStreamFile& ifile)
{
  // ASSUME ACCUMULATOR TYPE IS 'float'
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









