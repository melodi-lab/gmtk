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
    error("ERROR: invalid read cardinality (%d) < 0 in SPMF '%s' in file '%s'",
	  _card,name().c_str(),is.fileName());

  is.read(len,"Sparse1DPMF::read, len");
  if (len <= 0)
    error("ERROR: invalid read length (%d) < 0 in SPMF '%s' of cardinality %d in file '%s'",
	  len,name().c_str(),_card,is.fileName());

  pmf.resize(len);

  int prev_val = -1;
  for (int i=0;i<len;i++) {
    int val;

    is.read(val,"Sparse1DPMF::read, reading value");
    if (val < 0 || val > _card-1)
      error("ERROR: invalid value (%d) in SPMF '%s' of cardinality %d in file '%s'",
	  val,name().c_str(),_card,is.fileName());
    if (val <= prev_val)
      error("ERROR: values must be sorted (%d) in SPMF '%s' of cardinality %d in file '%s'",
	  val,name().c_str(),_card,is.fileName());

    prev_val = val;

    pmf[i] = val;
    
  }

  // read name of dense 1d PMF to use
  string str;
  is.read(str);
  if (GM_Parms.dPmfsMap.find(str) == GM_Parms.dPmfsMap.end()) {
    error("ERROR: in SPMF '%s', can't find DPMF named '%s' when reading file '%s'\n",
	  name().c_str(),str.c_str(),is.fileName());
  }

  dense1DPMF = GM_Parms.dPmfs[
			      GM_Parms.dPmfsMap[str]
  ];

  // now make sure length is right.
  if ((unsigned)len != dense1DPMF->length()) {
    error("ERROR: in SPMF '%s', DPMF named '%s' has %d elements but SPMF needs %d\n",
	  name().c_str(),
	  str.c_str(),dense1DPMF->length(),len);
  }
  setBasicAllocatedBit();
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
  assert ( basicAllocatedBitIsSet() );
  NamedObject::write(os);
  os.write(_card,"Sparse1DPMF::write, card");
  os.write(pmf.len(),"Sparse1DPMF::write, len");
  for (int i=0;i<pmf.len();i++) {
    os.write(pmf[i],"Sparse1DPMF::write, writing value");
  }
  os.nl();
  os.write(dense1DPMF->name());
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
  assert ( basicAllocatedBitIsSet() );
  assert ( pmf.len() > 0 );
  assert ( val >= 0 && val < _card );
  const int last = pmf.len()-1;
  if (val > pmf[last]) {
    logpr p; // make a log(zero) value.
    return p; // and return log(0).
  }

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

    if (val < pmf[m])
      u = m-1;
    else if (val > pmf[m])
      l = m+1;
    else
      return dense1DPMF->p(m);
  }

  logpr p;
  return p;

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
  dense1DPMF->normalize();
}

void
Sparse1DPMF::makeRandom()
{
  if (!emAmTrainingBitIsSet())
    return;

  dense1DPMF->makeRandom();
}

void
Sparse1DPMF::makeUniform()
{
  // NOTE: this does NOT assign non-zero values to "holes"
  // in the table (which are forced to be zero).
  if (!emAmTrainingBitIsSet())
    return;

  dense1DPMF->makeUniform();
}



////////////////////////////////////////////////////////////////////
//        EM Routines
////////////////////////////////////////////////////////////////////

 
void
Sparse1DPMF::emStartIteration()
{
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;

  if(emOnGoingBitIsSet())
    return; // already done

  if (!emEmAllocatedBitIsSet()) {
    emSetEmAllocatedBit();
  }
  // EM iteration is now going.
  emSetOnGoingBit();
  emSetSwappableBit();

  accumulatedProbability = 0.0;  

  if (dense1DPMF->length() != (unsigned)length()) {
    error("ERROR: Sparse PMF '%s' with '%d' possible values is trying to start an EM iteration using a dense PMF '%s' of length '%d'\n",
	  name().c_str(),length(),dense1DPMF->name().c_str(),
	  dense1DPMF->length());
  }

  dense1DPMF->emStartIteration();

}


void
Sparse1DPMF::emIncrement(logpr prob,const int val)
{
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;

  if (!emOnGoingBitIsSet())
    emStartIteration();

  assert ( pmf.len() > 0 );
  assert ( val >= 0 && val < _card );
  const int last = pmf.len()-1;

  /////////////////////////////////////////////////////////
  // we shouln't be trying to increment an impossible value.
  // something must be wrong upstairs.
  assert (val <= pmf[last]);

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

    if (val < pmf[m])
      u = m-1;
    else if (val > pmf[m])
      l = m+1;
    else {
      // found the value
      accumulatedProbability += prob;
      dense1DPMF->emIncrement(prob,m);
      return;
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
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;

  if (!emOnGoingBitIsSet())
    return; 

  accumulatedProbability.floor();
  if (accumulatedProbability.zero()) {
    warning("WARNING: Sparse 1D PMF named '%s' did not receive any accumulated probability in EM iteration",name().c_str());
  }

  dense1DPMF->emEndIteration();

  // stop EM
  emClearOnGoingBit();

}

void
Sparse1DPMF::emSwapCurAndNew()
{

  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;

  if (!emSwappableBitIsSet())
    return;

  dense1DPMF->emSwapCurAndNew();

  emClearSwappableBit();
}


void
Sparse1DPMF::emStoreAccumulators(oDataStreamFile& ofile)
{
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;
  if ( !emEmAllocatedBitIsSet() ) {
    warning("WARNING: storing zero accumulators for SPMF '%s'\n",
	    name().c_str());
    emStoreZeroAccumulators(ofile);
    return;
  }
  EMable::emStoreAccumulators(ofile);
}


void
Sparse1DPMF::emStoreZeroAccumulators(oDataStreamFile& ofile)
{
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;
  EMable::emStoreZeroAccumulators(ofile);
}


void
Sparse1DPMF::emLoadAccumulators(iDataStreamFile& ifile)
{
  assert (basicAllocatedBitIsSet());
  if (!emAmTrainingBitIsSet())
    return;
  assert (emEmAllocatedBitIsSet());
  EMable::emLoadAccumulators(ifile);
}


void
Sparse1DPMF::emAccumulateAccumulators(iDataStreamFile& ifile)
{
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;
  assert ( emEmAllocatedBitIsSet() );
  EMable::emAccumulateAccumulators(ifile);
}


