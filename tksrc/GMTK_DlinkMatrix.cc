/*-
 * GMTK_DlinkMatrix.cc
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2001, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle, make no representations about
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
#include "lineqsolve.h"

#include "GMTK_DlinkMatrix.h"
#include "GMTK_Dlinks.h"
#include "GMTK_GMParms.h"
#include "GMTK_GaussianComponent.h"
#include "GMTK_MixGaussiansCommon.h"
#include "GMTK_MeanVector.h"
#include "GMTK_DiagCovarVector.h"

VCID("$Header$");


////////////////////////////////////////////////////////////////////
//        Static members
////////////////////////////////////////////////////////////////////

double DlinkMatrix::cloneSTDfrac = 0.0;

void DlinkMatrix::checkForValidValues()
{
  if (DlinkMatrix::cloneSTDfrac < 0)
    error("ERROR: DlinkMatrix's cloneSTDfrac (%e) must be >= 0",
	  DlinkMatrix::cloneSTDfrac);
}



////////////////////////////////////////////////////////////////////
//        General create, read, destroy routines 
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * Dlinkmatrix::Dlinkmatrix()
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
DlinkMatrix::DlinkMatrix() 
{
}




////////////////////////////////////////////////////////////////////
//        Misc Support
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * read(is)
 *      read in the array from file 'is'. 
 *      The data probs are stored as doubles, but when they are read in
 *      they are converted to the log domain. Also, they are read
 *      in as a single array for speed reasons.
 *
 * Preconditions:
 *      none
 *
 * Postconditions:
 *      object is read in.
 *
 * Side Effects:
 *      Changes the pmf member function in the object.
 *
 * Results:
 *      nil
 *
 *-----------------------------------------------------------------------
 */
void
DlinkMatrix::read(iDataStreamFile& is)
{
  NamedObject::read(is);
  string str;

  // read the dlink structure 
  is.read(str);
  if (GM_Parms.dLinksMap.find(str) == GM_Parms.dLinksMap.end())
    error("Error: Dlink matrix '%s' in file '%s' specifies dlink structure name '%s' that does not exist",
	  name().c_str(),
	  is.fileName(),
	  str.c_str());
  
  dLinks = GM_Parms.dLinks[GM_Parms.dLinksMap[str]];

  int _dim;

  is.read(_dim,"DlinkMatrix::read, _dim");
  if (_dim <= 0)
    error("ERROR: reading DlinkMatrix '%s' from file '%s', dim (%d) must be positive",
	  name().c_str(),
	  is.fileName(),
	  _dim);

  if (_dim != dLinks->dim())
    error("Error: Dlink matrix '%s' in file '%s' specifices a dlink structure '%s' with incompatible dimensions\n",
	  name().c_str(),
	  is.fileName(),
	  dLinks->name().c_str());


  for (int i=0;i<_dim;i++) {
    int nlinks;
    is.read(nlinks,"DlinkMatrix::read, nlinks");

    if (nlinks < 0) 
      error("ERROR: reading DlinkMatrix '%s' from file '%s', # dlinks (%d) must be >= 0",
	    name().c_str(),
	    is.fileName(),
	    nlinks);

    if (nlinks != dLinks->numLinks(i))
      error("Error: Dlink matrix '%s' in file '%s' specifices a dlink structure '%s' with incompatible number of links\n",
	    name().c_str(),
	    is.fileName(),
	    dLinks->name().c_str());

    int oldLen = arr.len();
    arr.resizeAndCopy(oldLen+nlinks);

    for (int j=0;j<nlinks;j++) {
      is.read(arr[oldLen+j],"DlinkMatrix::read, v");
    }
  }

  setBasicAllocatedBit();
}




/*-
 *-----------------------------------------------------------------------
 * DlinkMatrix::write(os)
 *      write out data to file 'os'. 
 * 
 * Results:
 *      No results.
 *
 * Side Effects:
 *      No effects other than  moving the file pointer of os.
 *
 *-----------------------------------------------------------------------
 */
void
DlinkMatrix::write(oDataStreamFile& os)
{
  NamedObject::write(os);
  os.nl();

  os.write(dLinks->name(),"DlinkMatrix::write, dlm");
  os.nl();

  os.write(dim(),"DlinkMatrix::write, dim()");
  os.nl();
  int ptr = 0;
  for (int i=0;i<dim();i++) {
    os.write(dLinks->numLinks(i),"DlinkMatrix::write, nlinks");
    for (int j=0;j<dLinks->numLinks(i);j++) {
      os.write(arr[ptr++],"DlinkMatrix::write val ");
    }
    os.nl();
  }
}




/*-
 *-----------------------------------------------------------------------
 * makeRandom()
 *      assign random values to all elements
 * 
 * Preconditions:
 *      Object must be allocated.
 *
 * Postconditions:
 *      Object has random values.
 *
 * Side Effects:
 *      destroys previous values.
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
void
DlinkMatrix::makeRandom()
{
  if (!emAmTrainingBitIsSet())
    return;

  for (int i=0;i<arr.len();i++)
    arr[i] = rnd.drand48pe();
}



/*-
 *-----------------------------------------------------------------------
 * makeUniform()
 *      assign uniform (i.e., in this case 0) values to all elements
 * 
 * Preconditions:
 *      Object must be allocated.
 *
 * Postconditions:
 *      Object has random values.
 *
 * Side Effects:
 *      destroys previous values.
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
void
DlinkMatrix::makeUniform()
{
  if (!emAmTrainingBitIsSet())
    return;

  for (int i=0;i<arr.len();i++)
    arr[i] = 0;
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
DlinkMatrix*
DlinkMatrix::noisyClone()
{
  assert ( basicAllocatedBitIsSet() );

  // first check if self is already cloned, and if so, return that.
  DlinkMatrix* clone;

  map<DlinkMatrix*,DlinkMatrix*>::iterator it = 
    MixGaussiansCommon::dLinkMatCloneMap.find(this);
  if (it == MixGaussiansCommon::dLinkMatCloneMap.end()) {
    clone = new DlinkMatrix();
    // make sure we get a unique name
    unsigned cloneNo=0; do {
      char buff[256];
      sprintf(buff,"%d",cloneNo);
      clone->_name = _name + string("_cl") + buff;
      cloneNo++;
    } while (GM_Parms.dLinkMatsMap.find(clone->_name) != GM_Parms.dLinkMatsMap.end());
    clone->refCount = 0;
    clone->dLinks = dLinks;

    clone->arr.resize(arr.len());
    for (int i=0;i<arr.len();i++) {
      clone->arr[i] = arr[i] + 
	cloneSTDfrac*arr[i]*rnd.normal();
    }

    clone->setBasicAllocatedBit();
    MixGaussiansCommon::dLinkMatCloneMap[this] = clone;

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
DlinkMatrix::emStartIteration(sArray<float>& xzAccumulators,
			      sArray<float>& zzAccumulators,
			      sArray<float>& zAccumulators)
{
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;

  /////////////////////////////////////////////
  // make sure our caller has its accumulator resized
  // and initialized.
  xzAccumulators.growIfNeeded(dLinks->totalNumberLinks());
  for (int i=0;i<xzAccumulators.len();i++) {
    xzAccumulators[i] = 0.0;
  }
  zzAccumulators.growIfNeeded(dLinks->zzAccumulatorLength());
  for (int i=0;i<zzAccumulators.len();i++) {
    zzAccumulators[i] = 0.0;
  }
  zAccumulators.growIfNeeded(dLinks->totalNumberLinks());
  for (int i=0;i<zAccumulators.len();i++) {
    zAccumulators[i] = 0.0;
  }

  sharedZZDenominator.growIfNeeded(zzAccumulators.len());
  for (int i=0;i<sharedZZDenominator.len();i++) {
    sharedZZDenominator[i] = 0.0;
  }

  if(emOnGoingBitIsSet()) {
    // EM already on going.
    // Increment the count of number of Gaussian Components using this mean.
    refCount++;
    // this object therefore is shared, set the bit saying so.
    emSetSharedBit();

    if (refCount == 2) {
      // initialize the denominator if this is the
      // first time that we find out that this mean is shared
      // (known by the fact that refCount == 2).
      sharedZZDenominator.growIfNeeded(zzAccumulators.len());
      for (int i=0;i<sharedZZDenominator.len();i++) {
	sharedZZDenominator[i] = 0.0;
      }
    }
    return;
  }

  if (!emEmAllocatedBitIsSet()) {
    // this is presumably the first time
    emSetEmAllocatedBit();
    // allocate the final next means if needed
    nextArr.growIfNeeded(dLinks->totalNumberLinks());
  }

  // EM iteration is now going.
  emSetOnGoingBit();

  accumulatedProbability = 0.0;
  refCount = 1;
  emClearSharedBit();
  for (int i=0;i<nextArr.len();i++) {
    nextArr[i] = 0.0;
  }

  // make it swapable, although at this point
  // it would swap in the unaccumulated values.
  emSetSwappableBit();
}


/*-
 *-----------------------------------------------------------------------
 * emIncrement
 *      Add the data item for the current f into the
 *      partial accumulators for this dlink object.
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
DlinkMatrix::emIncrement(const logpr prob,
			 const float fprob,
			 const float* const f,
			 const Data32* const base,
			 const int stride,
			 float* xzAccumulators,
			 float* zzAccumulators,
			 float* zAccumulators)
{
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;

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
  assert ( prob >= minIncrementProbabilty );

  accumulatedProbability += prob;

  // compute the (possibly) shared cache arrays
  if (dLinks->totalNumberLinks() == 0)
    return;
  dLinks->cacheArrays(base,f);

  // This routine is called often, so we use pointer arithmetic

  // accumulate the zx and z counters together since they have
  // the same length.
  float *xzAccumulators_p = xzAccumulators;
  float *xzArrayCache_p = dLinks->xzArrayCache.ptr;
  float *zAccumulators_p = zAccumulators;
  float *zArrayCache_p = dLinks->zArrayCache.ptr;
  float *xzArrayCache_endp = xzArrayCache_p + 
    dLinks->totalNumberLinks();
  do {
    (*xzAccumulators_p++) += (*xzArrayCache_p++) * fprob;
    (*zAccumulators_p++) += (*zArrayCache_p++) * fprob;
  } while (xzArrayCache_p != xzArrayCache_endp);

  // accumulate the zz counters  
  float *zzAccumulators_p = zzAccumulators;
  float *zzArrayCache_p = dLinks->zzArrayCache.ptr;
  float *zzArrayCache_endp = zzArrayCache_p + 
    dLinks->zzAccumulatorLength();
  do {
    (*zzAccumulators_p++) += (*zzArrayCache_p++) * fprob;
  } while (zzArrayCache_p != zzArrayCache_endp);


}

#if 0

void
DlinkMatrix::emEndIteration(const float*const xzAccumulators)
{
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;

  if (refCount > 0) {
    // if this isn't the case, something is wrong.
    assert ( emOnGoingBitIsSet() );

    for (int i=0;i<nextArr.len();i++) {
      nextArr[i] += xzAccumulators[i];
    }

    refCount--;
  }

  /////////////////////////////////////////////
  // if there is still someone who
  // has not given us his/her accumulators
  // then we return w/o finishing.
  if (refCount > 0)
    return;

  // accumulatedProbability.floor();
  if (accumulatedProbability < minContAccumulatedProbability()) {
    warning("WARNING: dLink matrx '%s' received only %e accumulated log probability in EM iteration, using previous matrix",
	    name().c_str(),
	    accumulatedProbability.val());
    for (int i=0;i<nextArr.len();i++)
      nextArr[i] = arr[i];
  } else {
    const double invRealAccumulatedProbability =
      accumulatedProbability.inverse().unlog();
    // finish computing the next means.
    for (int i=0;i<nextArr.len();i++) {
      nextArr[i] *= invRealAccumulatedProbability;
    }
  }

  // stop EM
  emClearOnGoingBit();
}

#endif



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
 *      Note: this routine is probably a GEM rather than an EM.
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
DlinkMatrix::emEndIterationSharedMeansCovarsDlinks(const float*const xzAccumulators,
						   const float*const zAccumulators,
						   const float*const zzAccumulators,
						   const MeanVector* mean,
						   const DiagCovarVector* covar)
{
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;

  if (refCount > 0) {
    // if this isn't the case, something is wrong.
    assert ( emOnGoingBitIsSet() );

    // grab a pointer to the previous inverse variances
    // needed for normalization.
    const float* previous_variances_inv_ptr = covar->variances_inv.ptr;
    // previous mean
    const float* prev_mean_ptr = mean->means.ptr;

    const float* xzAccumulators_ptr = xzAccumulators;
    const float* zAccumulators_ptr = zAccumulators;
    const float* zzAccumulators_ptr = zzAccumulators;    

    float *nextArr_ptr = nextArr.ptr;

    double *sharedZZDenominator_ptr = sharedZZDenominator.ptr;

    for (int i=0;i<dim();i++) {
      
      const int nLinks = numLinks(i);

      for (int j=0;j<nLinks;j++) {
	double tmp = 
	  (double)previous_variances_inv_ptr[i]*
	  (
	   (double)(*xzAccumulators_ptr++)
	   -
	   ((double)(prev_mean_ptr[i])*(*zAccumulators_ptr++))
	   );
	*nextArr_ptr++ += tmp;

	for (int k=0;k<nLinks;k++) {
	  *sharedZZDenominator_ptr++ 
	    +=
	    previous_variances_inv_ptr[i]*
	    (*zzAccumulators_ptr++);
	}
      }
    }

    refCount--;
  }

  /////////////////////////////////////////////
  // if there is still someone who
  // has not given us his/her accumulators
  // then we return w/o finishing.
  if (refCount > 0)
    return;

  // accumulatedProbability.floor();
  if (accumulatedProbability < minContAccumulatedProbability()) {
    warning("WARNING: shared dLink matrix '%s' received only %e accumulated log probability in EM iteration, using previous matrix",
	    name().c_str(),
	    accumulatedProbability.val());
    for (int i=0;i<nextArr.len();i++)
      nextArr[i] = arr[i];
  } else {

    sArray<double> nextDlinkMat;

    double *sharedZZDenominator_ptr = sharedZZDenominator.ptr;
    float *nextArr_ptr = nextArr.ptr;

    unsigned numFlooredDlinks = 0;
    double minDlinkValue = DBL_MAX;
    for (int i=0;i<dim();i++) {
      const int nLinks = numLinks(i);

      nextDlinkMat.growIfNeeded(nLinks);
      // copy and convert to double
      for (int j=0;j<nLinks;j++) {
	nextDlinkMat[j] = nextArr_ptr[j];
      }
      
      // solve for the link values putting the results in
      // nextDlinkMat destroying the old values.
      ::lineqsolve(nLinks,1,sharedZZDenominator_ptr,nextDlinkMat.ptr);

      // copy out and convert back to single precision.
      for (int j=0;j<nLinks;j++) {
	if (nextDlinkMat[j] > FLT_MIN)
	  nextArr_ptr[j] = nextDlinkMat[j];
	else {
	  nextArr_ptr[j] = FLT_MIN;
	  numFlooredDlinks ++;
	}
	if (nextDlinkMat[j] < minDlinkValue)
	  minDlinkValue = nextDlinkMat[j];
      }

      nextArr_ptr += nLinks;
      sharedZZDenominator_ptr += (nLinks*nLinks);
    }

    if (numFlooredDlinks > 0) {
      warning("WARNING: shared dlink matrix '%s' had %d dlinks floored to %e, minimum dlink value found was %e.\n",
	      name().c_str(),
	      numFlooredDlinks,
	      FLT_MIN,
	      minDlinkValue);
    }
  }

  // stop EM
  emClearOnGoingBit();
}




/*-
 *-----------------------------------------------------------------------
 * emEndIterationNoSharingAlreadyNormalized()
 *      end the EM iteration for this mean object, where we have no
 *      sharing, but the mean has already been normalized. We still
 *      need to check for a small accumulator probability in which case we 
 *      just use the previous dlinks (rather than the new ones).
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
DlinkMatrix::emEndIterationNoSharingAlreadyNormalized(const float*const xzAccumulators)
{
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;

  // if this isn't the case, something is wrong.
  assert ( emOnGoingBitIsSet() );

  // shouldn't be called when sharing occurs, ensure this.
  assert ( refCount == 1 );
  assert (!emSharedBitIsSet());
  refCount = 0;

  // accumulatedProbability.floor();
  if (accumulatedProbability < minContAccumulatedProbability()) {
    warning("WARNING: dLink matrx '%s' received only %e accumulated log probability in EM iteration, using previous matrix",
	    name().c_str(),
	    accumulatedProbability.val());
    for (int i=0;i<nextArr.len();i++)
      nextArr[i] = arr[i];
  } else {
    // finish computing the next means.
    for (int i=0;i<nextArr.len();i++) {
      nextArr[i] = xzAccumulators[i];
    }
  }

  // stop EM
  emClearOnGoingBit();
}


void
DlinkMatrix::emSwapCurAndNew()
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
  for (int i=0;i<arr.len();i++) {
    genSwap(arr[i],nextArr[i]);
  }
  // make no longer swappable
  emClearSwappableBit();
}

void
DlinkMatrix::emStoreAccumulators(oDataStreamFile& ofile)
{
  assert (basicAllocatedBitIsSet());
  if (!emAmTrainingBitIsSet())
    return;
  if ( !emEmAllocatedBitIsSet() ) {
    warning("WARNING: storing zero accumulators for dlink matrix '%s'\n",
	    name().c_str());
    emStoreZeroAccumulators(ofile);
    return;
  }
  EMable::emStoreAccumulators(ofile);
}


void
DlinkMatrix::emStoreZeroAccumulators(oDataStreamFile& ofile)
{
  assert (basicAllocatedBitIsSet());
  if (!emAmTrainingBitIsSet())
    return;
  EMable::emStoreZeroAccumulators(ofile);
}


void
DlinkMatrix::emLoadAccumulators(iDataStreamFile& ifile)
{
  assert (basicAllocatedBitIsSet());
  if (!emAmTrainingBitIsSet())
    return;
  assert (emEmAllocatedBitIsSet());
  EMable::emLoadAccumulators(ifile);
}


void
DlinkMatrix::emAccumulateAccumulators(iDataStreamFile& ifile)
{
  assert (basicAllocatedBitIsSet());
  if (!emAmTrainingBitIsSet())
    return;
  assert (emEmAllocatedBitIsSet());
  EMable::emAccumulateAccumulators(ifile);
}



