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
#include <ieeefp.h>
#include <float.h>
#include <assert.h>

#include "general.h"
#include "error.h"
#include "rand.h"

#include "GMTK_DlinkMatrix.h"
#include "GMTK_Dlinks.h"

VCID("$Header$");


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
	  name().c_str(),is.fileName(),str.c_str());
  
  dLinks = GM_Parms.dLinks[GM_Parms.dLinksMap[str]];

  int _dim;

  is.read(_dim,"DlinkMatrix::read, _dim");
  if (_dim <= 0)
    error("ERROR: reading DlinkMatrix '%s' from file '%s', dim (%d) must be positive",name().c_str(),is.fileName(),_dim);

  _numLinks.resize(_dim);

  for (int i=0;i<_dim;i++) {
    int nlinks;
    is.read(nlinks,"DlinkMatrix::read, nlinks");

    if (nlinks < 0) 
      error("ERROR: reading DlinkMatrix '%s' from file '%s', # dlinks (%d) must be >= 0",name().c_str(),is.fileName(),nlinks);

    _numLinks[i] = nlinks;

    int oldLen = arr.len();
    arr.resizeAndCopy(oldLen+nlinks);

    for (int j=0;j<nlinks;j++) {
      is.read(arr[oldLen+j],"DlinkMatrix::read, v");
    }
  }

  if (!compatibleWith(*dLinks))
    error("Error: Dlink matrix '%s' in file '%s' specifices a dlink structure '%s' is incompatible\n",
	  name().c_str(),is.fileName(),
	  dLinks->name().c_str());

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
  os.write(dim(),"DlinkMatrix::write, dim()");
  os.nl();
  int ptr = 0;
  for (int i=0;i<dim();i++) {
    os.write(_numLinks[i],"DlinkMatrix::write, nlinks");
    for (int j=0;j<_numLinks[i];j++) {
      os.write(arr[ptr++],"DlinkMatrix::write val ");
    }
    os.nl();
  }
}


/*-
 *-----------------------------------------------------------------------
 * preCompute()
 *      precompute the offset array.
 * 
 * Preconditions:
 *      basic object should be read in.
 *
 * Postconditions:
 *      offset array is computed.
 *
 * Side Effects:
 *      destroys old offset array
 *
 * Results:
 *      retunrs nil
 *
 *-----------------------------------------------------------------------
 */
void
DlinkMatrix::preCompute(const unsigned stride)
{
  dLinks->preCompute(stride);

  unsigned len = 0;
  zzAccumulatorLength = 0;
  unsigned maxDlinks = 0;
  for (int i=0;i<dim();i++) {
    if (numLinks(i) > maxDlinks)
      maxDlinks = numLinks(i);
    len += numLinks(i);
    // only storing upper triangular portion of symetric matrix
    zzAccumulatorLength += numLinks(i)*(numLinks(i)+1)/2;
  }

  // now allocate the cache
  zArray.resize(maxDlinks);
  zzArrayCache.resize(zzAccumulatorLength);
  xzArrayCache.resize(len);
  clearArrayCache();
}



/*-
 *-----------------------------------------------------------------------
 * compatibleWith()
 *      returns true of this object is compatible with the argument function.
 * 
 * Preconditions:
 *      Object must be read in.
 *
 * Postconditions:
 *      --
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      true only if compatibility holds.
 *
 *-----------------------------------------------------------------------
 */
bool 
DlinkMatrix::compatibleWith(Dlinks& d)
{
  if (dim() != d.dim())
    return false;
  for (int i=0;i<dim();i++) {
    if (numLinks(i) != d.numLinks(i))
      return false;
  }
  return true;
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
  for (int i=0;i<arr.len();i++)
    arr[i] = 0;
}



/*-
 *-----------------------------------------------------------------------
 * Function
 *      What the function does.
 * 
 * Preconditions:
 *      What must be true before the function is called.
 *
 * Postconditions:
 *      What is true after the function is called.
 *
 * Side Effects:
 *      Any exernal side effects
 *
 * Results:
 *      What does the function return, if anything. 
 *
 *-----------------------------------------------------------------------
 */
void
DlinkMatrix::clearArrayCache()
{
  // assume that a NULL pointer won't match anything.
  arrayCacheTag = NULL;
}


/*-
 *-----------------------------------------------------------------------
 * Function
 *      What the function does.
 * 
 * Preconditions:
 *      What must be true before the function is called.
 *
 * Postconditions:
 *      What is true after the function is called.
 *
 * Side Effects:
 *      Any exernal side effects
 *
 * Results:
 *      What does the function return, if anything. 
 *
 *-----------------------------------------------------------------------
 */
void
DlinkMatrix::cacheArrays(const Data32* const base, const float *const f)
{
  if (base == arrayCacheTag)
    return;

  // refill the cache.

  const int* lagStrideOffsetsp = dLinks->preComputedOffsets.ptr;
  float *zzArrayCachep = zzArrayCache.ptr;
  float *xzArrayCachep = xzArrayCache.ptr;
  const float *fp = f;
  int i=0; do {
    const int nLinks = numLinks(i);
    if (nLinks > 0) {
      const float feat = *fp;
      float *zArray_p = zArray.ptr;
      const int *const lagStrideOffsets_endp = 
	lagStrideOffsetsp + nLinks;
      do { 
	const float zval = *((float*)base + *lagStrideOffsetsp++);
	*zArray_p++ = zval;
	*xzArrayCachep++ = feat*zval;
      } while (lagStrideOffsetsp != lagStrideOffsets_endp);
      matrixSelfOuterProduct(zArray,nLinks,zzArrayCachep);
      zzArrayCachep += nLinks*(nLinks+1)/2;
    }
    fp++;
  } while (++i != dim());

  arrayCacheTag = base;
}


/////////////////
// EM routines //
/////////////////



void
DlinkMatrix::emStartIteration(sArray<float>& xzAccumulators,
			      sArray<float>& zzAccumulators)
{
  assert ( basicAllocatedBitIsSet() );

  if (!GM_Parms.amTrainingMeans())
    return;

  /////////////////////////////////////////////
  // make sure our caller has its accumulator resized
  // and initialized.
  xzAccumulators.growIfNeeded(dLinks->totalNumberLinks());
  for (int i=0;i<xzAccumulators.len();i++) {
    xzAccumulators = 0.0;
  }
  zzAccumulators.growIfNeeded(dLinks->zzAccumulatorLength);
  for (int i=0;i<zzAccumulators.len();i++) {
    zzAccumulators = 0.0;
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
    nextArr.growIfNeeded(dLinks->totalNumberLinks());
  }

  // EM iteration is now going.
  emSetOnGoingBit();

  accumulatedProbability = 0.0;
  refCount = 1;
  for (int i=0;i<nextArr.len();i++) {
    nextArr[i] = 0.0;
  }

  // make it swapable, although at this point
  // it would swap in the unaccumulated values.
  emSetSwappableBit();
}


void
DlinkMatrix::emIncrement(const logpr prob,
			 const float fprob,
			 const float* const f,
			 const Data32* const base,
			 const int stride,
			 float* xzAccumulators,
			 float* zzAccumulators)
{
  assert ( basicAllocatedBitIsSet() );

  if (!GM_Parms.amTrainingMeans())
    return;
  
  /////////////////////////////////////////////
  // Note: unlike the normal EM mode described
  // in GMTK_EMable.h, we do not call
  // emStartIteration() here and assume that it
  // was called by the Gaussian component that
  // is using this mean. This is because
  // this object keeps a reference count (needed for
  // sharing), and calling that routine repeatedly 
  // would result in an incorrect count.

  // we assume here that (prob > minIncrementProbabilty),
  // i.e., that this condition has been checked by the caller
  // of this routine (meaning that fprob is valid)
  assert ( prob >= minIncrementProbabilty );

  accumulatedProbability += prob;

  // compute the possibly shared cache arrays
  cacheArrays(base,f);

  // This routine is called often so we use pointer arithmetic where possible.
  const int* lagStrideOffsetsp = dLinks->preComputedOffsets.ptr;
  float *xzAccumulators_p = xzAccumulators;
  float *zzAccumulators_p = zzAccumulators;




}


void
DlinkMatrix::emEndIteration(const float*const partialAccumulatedNextMeans)
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
DlinkMatrix::emSwapCurAndNew()
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
DlinkMatrix::emStoreAccumulators(oDataStreamFile& ofile)
{
  assert ( basicAllocatedBitIsSet() );
  error("not implemented");
}

void
DlinkMatrix::emLoadAccumulators(iDataStreamFile& ifile)
{
  assert ( basicAllocatedBitIsSet() );
  error("not implemented");
}


void
DlinkMatrix::emAccumulateAccumulators(iDataStreamFile& ifile)
{
  assert ( basicAllocatedBitIsSet() );
  error("not implemented");
}



