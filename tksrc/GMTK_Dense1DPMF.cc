/*-
 * GMTK_Dense1DPMF.cc
 *     Trainable (with say EM) 1D discrete probability
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

#include "GMTK_Dense1DPMF.h"
#include "GMTK_GMParms.h"

VCID("$Header$");


////////////////////////////////////////////////////////////////////
//        General create, read, destroy routines 
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * Dense1DPMF::Dense1DPMF()
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
Dense1DPMF::Dense1DPMF() 
{

}


/*-
 *-----------------------------------------------------------------------
 * Dense1DPMF::read(is)
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
Dense1DPMF::read(iDataStreamFile& is)
{

  NamedObject::read(is);
  int length;
  is.read(length,"Dense1DPMF::read, distribution length");
  if (length <= 0)
    error("ERROR: DPMF '%s', file '%s', bad length (%d) < 0 in input",
	  name().c_str(),
	  is.fileName(),length);
  pmf.resize(length);
  for (int i=0;i<length;i++) {
    double prob;
    is.readDouble(prob,"Dense1DPMF::read, reading prob");
    if (prob < 0 || prob > 1)
      error("ERROR: DPMF '%s', file '%s', bad pmf value (%g)",
	    name().c_str(),
	    is.fileName(),prob);
    pmf[i] = prob;
  }
}




/*-
 *-----------------------------------------------------------------------
 * Dense1DPMF::write(is)
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
Dense1DPMF::write(oDataStreamFile& os)
{
  NamedObject::write(os);
  os.write(pmf.len(),"Dense1DPMF::write, distribution length");
  for (int i=0;i<pmf.len();i++) {
    // convert out of log domain and write out.
    os.writeDouble(pmf[i].unlog(),"Dense1DPMF::write, writing prob");
  }
  os.nl();
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
Dense1DPMF::normalize()
{
  logpr sum = 0.0;
  for (int i=0;i<pmf.len();i++) {
    sum += pmf[i];
  }
  for (int i=0;i<pmf.len();i++) {
    pmf[i] = pmf[i] / sum;
  }
}

void
Dense1DPMF::makeRandom()
{
  logpr sum = 0.0;
  for (int i=0;i<pmf.len();i++) {
    logpr tmp = rnd.drand48();
    sum += tmp;
    pmf[i] = tmp;
  }
  for (int i=0;i<pmf.len();i++) {
    pmf[i] = pmf[i] / sum;
  }
}

void
Dense1DPMF::makeUniform()
{
  logpr val = 1.0/pmf.len();
  for (int i=0;i<pmf.len();i++) {
    pmf[i] = val;
  }
}



/////////////////
// EM routines //
/////////////////


void
Dense1DPMF::emStartIteration()
{
  if (!GM_Parms.amTrainingDense1DPMFs())
    return;
  if(emOnGoingBitIsSet())
    return; 

  if (!emEmAllocatedBitIsSet()) {
    // this is presumably the first time
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


/*-
 *-----------------------------------------------------------------------
 * emIncrement()
 *      Increment postDist into the next parameters for this PMF.
 * 
 * Preconditions:
 *      prob must be valid probability, and postDist
 *      a posterior distribution that has *ALREADY* had 
 *      'prob' multiplied in by the caller.
 *
 * Postconditions:
 *      postDist has been accumulated in.
 *
 * Side Effects:
 *      Changes internal variables.
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */

void
Dense1DPMF::emIncrement(logpr prob,
			sArray<logpr>& postDist)
{
  if (!GM_Parms.amTrainingDense1DPMFs())
    return;
  emStartIteration();

  assert ( postDist.len() == nextPmf.len() );

  if (prob < minIncrementProbabilty) {
    missedIncrementCount++;
    return;
  } 
  accumulatedProbability+= prob;
  for (int i=0;i<nextPmf.len();i++) {
    // we assume here that 'prob' has already
    // been multipled into postDist by our friendly
    // caller. If this is not the case, then
    // this code would be wrong.
    nextPmf[i] += postDist[i];
  }
}



/*-
 *-----------------------------------------------------------------------
 * emIncrement()
 *      Increment a single value count into the next parameters for this PMF.
 * 
 * Preconditions:
 *      prob must be valid probability, and val
 *      must indicate a valid entry into the array.
 *
 * Postconditions:
 *      val of post. val has been accumulated in.
 *
 * Side Effects:
 *      Changes internal variables.
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */

void
Dense1DPMF::emIncrement(logpr prob,
			const int val)
{
  assert ( val >= 0 && val < nextPmf.len() );

  if (!GM_Parms.amTrainingDense1DPMFs())
    return;
  emStartIteration();

  if (prob < minIncrementProbabilty) {
    missedIncrementCount++;
    return;
  } 
  accumulatedProbability+= prob;
  nextPmf[val] += prob;
}


void
Dense1DPMF::emEndIteration()
{
  if (!GM_Parms.amTrainingDense1DPMFs())
    return;

  if ( !emOnGoingBitIsSet() )
    return; 

  if (accumulatedProbability.zero()) {
    warning("WARNING: Dense1DPMF named '%s' did not receive any accumulated probability in EM iteration",name().c_str());
  }

  for (int i=0;i<nextPmf.len();i++) {
    nextPmf[i] /= accumulatedProbability;
  }

  // stop EM
  emClearOnGoingBit();
}



void
Dense1DPMF::emSwapCurAndNew()
{
  if (!GM_Parms.amTrainingDense1DPMFs())
    return;

  if (!emSwappableBitIsSet())
    return;

  for (int i=0;i<nextPmf.len();i++) {
    genSwap(nextPmf[i],pmf[i]);
  }

  emClearSwappableBit();
}




void
Dense1DPMF::emStoreAccumulators(oDataStreamFile& ofile)
{
  error("not implemented");
}

void
Dense1DPMF::emLoadAccumulators(iDataStreamFile& ifile)
{
  error("not implemented");
}


void
Dense1DPMF::emAccumulateAccumulators(iDataStreamFile& ifile)
{
  error("not implemented");
}


////////////////////////////////////////////////////////////
// Sample generation
////////////////////////////////////////////////////////////


