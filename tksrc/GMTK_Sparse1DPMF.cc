/*-
 * GMTK_Sparse1DPMF.cc
 *     Trainable (with say EM) sparse 1D discrete probability
 *     distributions.
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

#include "GMTK_Sparse1DPMF.h"
#include "GMTK_GMParms.h"

VCID("$Header$");


////////////////////////////////////////////////////////////////////
//        General create, read, destroy routines 
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * Sparse1DPMF::Sparse1DPMF()
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
Sparse1DPMF::Sparse1DPMF() 
{

}



/*-
 *-----------------------------------------------------------------------
 * Sparse1DPMF::read(is)
 *      read in a distribution from file 'is'. 
 *      The data probs are stored as doubles, but when they are read in
 *      they are converted into the log domain (so discrete data on disk
 *      is NOT represented as log probabilities. This is because 
 *          1) discrete data typically doesn't need such a huge dynamic range
 *              (like Gaussian probabilties do).
 *          2) it is easier to examine data on disk when it is not in log domain.
 * 
 * Results:
 *      No results.
 *
 * Side Effects:
 *      Changes the pmf member function in the object.
 *
 *-----------------------------------------------------------------------
 */
void
Sparse1DPMF::read(iDataStreamFile& is)
{
  int len;
  NamedObject::read(is);
  is.read(_card,"Sparse1DPMF::read, card");
  if (_card <= 0)
    error("Sparse1DPMF: read length (%d) < 0 in input",_card);

  is.read(len,"Sparse1DPMF::read, len");
  if (len <= 0)
    error("Sparse1DPMF: read length (%d) < 0 in input",len);
  if (len > _card)
    error("Sparse1DPMF: read length (%d) > card (%d) in input",len,_card);

  pmf.resize(len);

  int prev_val = -1;
  for (int i=0;i<len;i++) {
    int val;
    double prob;

    is.read(val,"Sparse1DPMF::read, reading value");
    if (val < 0 || val > _card-1)
      error("Sparse1DPMF::read, bad value = %d, must be in range [0:%d]",val,
	    _card-1);
    if (val <= prev_val)
      error("Sparse1DPMF::read, values must be unique and in sorted order");
    prev_val = val;

    is.readDouble(prob,"Sparse1DPMF::read, reading prob");
    if (prob < 0.0 || prob > 1.0)
      error("Sparse1DPMF: read, invalid pmf value (%g)",val);
    pmf[i].val = val;
    pmf[i].prob = prob;
  }
}




/*-
 *-----------------------------------------------------------------------
 * Sparse1DPMF::write(is)
 *      write out distribution to file 'os'. 
 *      the data probs are stored on disk as doubles,  NOT in log domain.
 * 
 * Results:
 *      No results.
 *
 * Side Effects:
 *      No effectcs other than  moving the file pointer of os.
 *
 *-----------------------------------------------------------------------
 */
void
Sparse1DPMF::write(oDataStreamFile& os)
{
  NamedObject::write(os);
  os.write(_card,"Sparse1DPMF::write, card");
  os.write(pmf.len(),"Sparse1DPMF::write, len");
  for (int i=0;i<pmf.len();i++) {
    os.write(pmf[i].val,"Sparse1DPMF::write, writing value");
    os.writeDouble(pmf[i].prob.unlog(),"Sparse1DPMF::write, writing prob");
  }
  os.nl();
}





////////////////////////////////////////////////////////////////////
//        Misc Support
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * prob
 *      Return the probability for value val
 * 
 * Preconditions:
 *      Must be filled in.
 *
 * Postconditions:
 *      same as preconditions.
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      returns the probability.
 *
 *-----------------------------------------------------------------------
 */

logpr
Sparse1DPMF::prob(const int val)
{
  assert ( pmf.len() > 0 );
  assert ( val >= 0 && val < _card );
  const int last = pmf.len()-1;
  if (val > pmf[last].val)
    return LZERO;

  int l,u;
  l = 0;
  u = last;
  while (l<=u) {
    // pmf[l] <= val <= pmf[u]
    int m = (l+u)/2;
    // if l = u, 
    //      then  m = l = u
    // if even(u), even(l)
    //      then l < m < u
    // if one even, one odd, and l+1 = u
    //      then  m = l = u-1
    // if one even, one odd, and l+1 < u
    //      then l < m < u

    if (val < pmf[m].val)
      u = m-1;
    else if (val > pmf[m].val)
      l = m+1;
    else
      return pmf[m].prob;
  }


  return LZERO;
}




////////////////////////////////////////////////////////////////////
//        Misc Support
////////////////////////////////////////////////////////////////////



/*-
 *-----------------------------------------------------------------------
 * normalize, makeRandom, and makeUniform
 *      normalizes, makes random, and makes uniform (resp.) the values. 
 * 
 * Preconditions:
 *      internals must be allocated
 *
 * Postconditions:
 *      values are normalized, randomized, or uniformized
 *
 * Side Effects:
 *      Changes all values of internal probabilities.
 *
 * Results:
 *      nothing
 *
 *-----------------------------------------------------------------------
 */
void
Sparse1DPMF::normalize()
{
  logpr sum = 0.0;
  for (int i=0;i<pmf.len();i++) {
    sum += pmf[i].prob;
  }
  for (int i=0;i<pmf.len();i++) {
    pmf[i].prob = pmf[i].prob / sum;
  }
}

void
Sparse1DPMF::makeRandom()
{
  logpr sum = 0.0;
  for (int i=0;i<pmf.len();i++) {
    logpr tmp = rnd.drand48();
    sum += tmp;
    pmf[i].prob = tmp;
  }
  for (int i=0;i<pmf.len();i++) {
    pmf[i].prob = pmf[i].prob / sum;
  }
}

void
Sparse1DPMF::makeUniform()
{
  // NOTE: this doesn't assign non-zero values to "holes"
  // in the table (which are forced to be zero).
  logpr p = 1.0/pmf.len();
  for (int i=0;i<pmf.len();i++) {
    pmf[i].prob = p;
  }
}



////////////////////////////////////////////////////////////////////
//        EM Routines
////////////////////////////////////////////////////////////////////

 
void
Sparse1DPMF::emStartIteration()
{
  if (!GM_Parms.amTrainingSparse1DPMFs())
    return;

  if(emOnGoingBitIsSet())
    return; // already done

  if (!emEmAllocatedBitIsSet()) {
    nextPmf.resize(pmf.len());
    emSetEmAllocatedBit();
  }
  // EM iteration is now going.
  emSetOnGoingBit();
  emSetSwappableBit();

  accumulatedProbability = 0.0;  

  for (int i=0;i<nextPmf.len();i++) {
    nextPmf[i].set_to_zero();
  }

}


void
Sparse1DPMF::emIncrement(logpr prob,const int val)
{
  if (!GM_Parms.amTrainingSparse1DPMFs())
    return;

  emStartIteration();

  assert ( pmf.len() > 0 );
  assert ( val >= 0 && val < _card );
  const int last = pmf.len()-1;

  /////////////////////////////////////////////////////////
  // we shouln't be trying to increment an impossible value.
  // something must be wrong upstairs.
  assert (val <= pmf[last].val);

  // find the entry using bin search
  int l,u;
  l = 0;
  u = last;
  while (l<=u) {
    // pmf[l] <= val <= pmf[u]
    int m = (l+u)/2;
    // if l = u, 
    //      then  m = l = u
    // if even(u), even(l)
    //      then l < m < u
    // if one even, one odd, and l+1 = u
    //      then  m = l = u-1
    // if one even, one odd, and l+1 < u
    //      then l < m < u

    if (val < pmf[m].val)
      u = m-1;
    else if (val > pmf[m].val)
      l = m+1;
    else {
      // found the value
      accumulatedProbability += prob;
      nextPmf[m] += prob;
    }
  }

  /////////////////////////////////////////////////////////
  // again, we shouln't be trying to increment an impossible value.
  // something must be wrong.
  assert (0);

}

void
Sparse1DPMF::emEndIteration()
{
  if (!GM_Parms.amTrainingSparse1DPMFs())
    return;

  if ( !emOnGoingBitIsSet() )
    return; 

  if (accumulatedProbability.zero()) {
    warning("WARNING: Sparse 1D PMF named '%s' did not receive any accumulated probability in EM iteration",name().c_str());
  }

  for (int i=0;i<nextPmf.len();i++) {
    nextPmf[i] /= accumulatedProbability;
  }

  // stop EM
  emClearOnGoingBit();

}

void
Sparse1DPMF::emSwapCurAndNew()
{
  if (!GM_Parms.amTrainingSparse1DPMFs())
    return;

  if (!emSwappableBitIsSet())
    return;

  for (int i=0;i<nextPmf.len();i++) {
    genSwap(nextPmf[i],pmf[i].prob);
  }

  emClearSwappableBit();
}


void
Sparse1DPMF::emStoreAccumulators(oDataStreamFile& ofile)
{
  error("not implemented");
}

void
Sparse1DPMF::emLoadAccumulators(iDataStreamFile& ifile)
{
  error("not implemented");
}


void
Sparse1DPMF::emAccumulateAccumulators(iDataStreamFile& ifile)
{
  error("not implemented");
}


