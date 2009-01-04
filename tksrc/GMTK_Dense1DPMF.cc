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
#include <float.h>
#include <assert.h>

#include "general.h"
#include "error.h"
#include "rand.h"

#include "GMTK_Dense1DPMF.h"
#include "GMTK_GMParms.h"
#include "GMTK_MixtureCommon.h"
#include "GMTK_CPT.h"
#include "tieSupport.h"


VCID("$Header$")

/*
 * This routine copies the magnitude of x, the sign of y, and returns the result, i.e.,
 *  copysign(x,y) = fabs(x)*sign(y)
 */      
extern "C" double copysign(double x, double y) __THROW;

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
  is.read(length,"Can't read Dense1DPMF's distribution length");
  if (length <= 0)
    error("ERROR: reading file '%s' line %d, DPMF '%s' has a bad length of (%d) <= 0 in input",
	  is.fileName(),is.lineNo(),name().c_str(),length);
  pmf.resize(length);


  // read optional smoothing parameters. We support either
  //   1) constant dirichlet priors, where a constant alpha
  //      is given right here, where we have a Dirichlet 
  //      with hyperparameter alpha >= 0 which is constant for
  //      all rv values. The syntax for this is:
  //         
  //         DirichletConst <alpha>
  // 
  //      and where the accumulators now have E[counts] + alpha rather
  //      than just alpha. Note that \alpha is given as a real fractional
  //      count value, so it is not a log prob.
  //   2) specify a counts object, where we have a full set of 'numValues'
  //      counts for all values of the random variable for all possible
  //      parent values. This is much more general than the above, as
  //      the counts object can specify a different Dirichlet hyperparameter
  //      for each value of the RV. Syntax for this is:
  //   
  //        DirichletTable table-object
  //
  //      where table-object is a previously defined Dirichlet Table object in the
  //      master file which is compatible with this table. Note that it
  //      has to be compatible at read time. Since a Dense1DPMF can increase
  //      in cardinality when used for Gaussian mixture responsibilities, once
  //      the cardinality increases, the Dirichlet Table is no longer valid and
  //      it is not any longer used. To get it to work again in this case, you
  //      must specify another table and re-read in the object.

  
  if (is.readIfMatch(DirichletConstStr,"DenseCPT double value or Dirichlet const spec")) {
    // so we should have a single constant alpha value next.
    is.read(dirichletAlpha,"Can't read DPMF Dirichlet hyperparameter");
    smoothingType = DirichletConstVal;
  } else if (is.readIfMatch(DirichletTableStr,"DenseCPT double value or Dirichlet table spec")) {
    // so we should have a pointer to a previously existing count table.
    string dirichletTableName;
    is.read(dirichletTableName);
    if (GM_Parms.dirichletTabsMap.find(dirichletTableName) == GM_Parms.dirichletTabsMap.end()) {
	error("ERROR: reading file '%s' line %d, DPMF '%s' specified Dirichlet Table (%s) that does not exist.",
	      is.fileName(),is.lineNo(),
	      name().c_str(),
	      dirichletTableName.c_str());

    }
    dirichletTable = GM_Parms.dirichletTabs[GM_Parms.dirichletTabsMap[dirichletTableName]];
    smoothingType = DirichletTableVal;
    // next check that the table matches the DPMF.
    if (dirichletTable->numDimensions() != 1) {
	error("ERROR: reading file '%s' line %d, DPMF '%s' (a 1-D table), but Dirichlet Table '%s' has dimensionality %d",
	      is.fileName(),is.lineNo(),
	      name().c_str(),
	      dirichletTable->name().c_str(),
	      dirichletTable->numDimensions());
    }

    // check self cardinality
    if (pmf.size() != dirichletTable->lastDimension()) {
      error("ERROR: reading file '%s' line %d, in DPMF '%s', has length %d, but Dirichlet Table '%s' has its last dimension of size %d",
	    is.fileName(),is.lineNo(),
	    name().c_str(),
	    pmf.size(),
	    dirichletTable->name().c_str(),
	    dirichletTable->lastDimension());
    }

    // everything matches up, but include last sanity check
    assert ( dirichletTable->tableSize() == (unsigned)length );

  }


  cachedMaxValue.set_to_zero();
  logpr sum;
  for (int i=0;i<length;i++) {

    double prob;
    is.readDouble(prob,"Can't read Dense1DPMF's prob");

    // we support reading in both regular probability values
    // (in the range [+0,1] inclusive) and log probability 
    // values (in the range (-infty,-0] inclusive. These
    // ranges give distinct values for probabilties, except for
    // the value 0 which can either be real probability zero (impossible
    // event) or it could be log(1) = 0 (the certain event). Since
    // the IEEE FP standard supports both +0 and -0, and since the
    // ASCII read routines preserve ASCII string '-0.0' to be negative zero,
    // we consider -0.0 as log(1) , and +0.0 as real zero.
    if (prob > 1)
      error("ERROR: reading file '%s' line %d, DPMF '%s' has invalid probability value (%e), entry %d",
	    is.fileName(),is.lineNo(),
	    name().c_str(),
	    prob,
	    i);
    if (prob > 0) {
      // regular probability
      pmf[i] = prob;
    } else if (prob < 0) {
      // log base e probability
      pmf[i].setFromLogP(prob);
    } else {
      // is zero, so need to check sign bit for
      // either -0 (log(1)) or +0 (true zero prob)
      if (copysign(1.0,prob)==1.0) {      
	// regular zero probability
	pmf[i].set_to_zero();
      } else {
	// prob == -0, so set to log(1)
	pmf[i].set_to_one();
      }
    }
    sum += pmf[i];
    if (pmf[i] > cachedMaxValue)
      cachedMaxValue = pmf[i];

  }


  if (CPT::normalizationThreshold != 0) {
    double abs_diff = fabs(sum.unlog() - 1.0);
    // be more forgiving as cardinality increases
    if (abs_diff > length*CPT::normalizationThreshold) 
      error("ERROR: reading file '%s' line %d, DPMF '%s' has probabilities that sum to %e but should sum to unity, absolute difference = %e, current normalization threshold = %f.",
	    is.fileName(),is.lineNo(),
	    name().c_str(),
	    sum.unlog(),
	    abs_diff,
	    CPT::normalizationThreshold);
  }

  setBasicAllocatedBit();
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
  assert ( basicAllocatedBitIsSet() );
  NamedObject::write(os);
  os.write(pmf.len(),"Dense1DPMF::write, distribution length");

  if (smoothingType == DirichletConstVal) {
    os.write(DirichletConstStr);
    os.write(dirichletAlpha);
    os.nl();
  } else if (smoothingType == DirichletTableVal) {
    os.write(DirichletTableStr);
    os.write(dirichletTable->name());
    os.nl();
  }


  normalize();
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
  assert ( basicAllocatedBitIsSet() );


  logpr sum = 0.0;
  for (int i=0;i<pmf.len();i++) {
    sum += pmf[i];
  }
  cachedMaxValue.set_to_zero();
  for (int i=0;i<pmf.len();i++) {
    pmf[i] = pmf[i] / sum;
    if (pmf[i] > cachedMaxValue)
      cachedMaxValue = pmf[i];
  }
}

void
Dense1DPMF::makeRandom()
{
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;

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
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;

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
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
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

  if (smoothingType == NoneVal || !useDirichletPriors) {
    for (int i=0;i<nextPmf.len();i++) {
      nextPmf[i].set_to_zero();
    }
  } else if (smoothingType == DirichletConstVal) {
    // const dirichlet priors 
    logpr alpha(dirichletAlpha);
    for (int i=0;i<nextPmf.len();i++) {
      nextPmf[i] = alpha;
    }
  } else {
    // table dirichlet priors 
    for (int i=0;i<nextPmf.len();i++) {
      nextPmf[i] = dirichletTable->tableValue(i);
    }
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
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;

  if (!emOnGoingBitIsSet()) 
    emStartIteration();

  assert ( postDist.len() == nextPmf.len() );

  if (prob < minIncrementProbabilty) {
    missedIncrementCount++;
    return;
  }
  accumulatedProbability+= prob;
  if (!fisherKernelMode) {
    for (int i=0;i<nextPmf.len();i++) {
      // we assume here that 'prob' has already
      // been multipled into postDist by our friendly
      // caller. If this is not the case, then
      // this code would be wrong.
      nextPmf[i] += postDist[i];
    }
  } else {
    // do the fisher kernel score case, the same as above but divide by
    // the previous parameters.
    for (int i=0;i<nextPmf.len();i++) {
      if (!pmf[i].zero())
	nextPmf[i] += postDist[i]/pmf[i];
    }
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
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;
  if (!emOnGoingBitIsSet()) 
    emStartIteration();

  assert ( val >= 0 && val < nextPmf.len() );

  if (prob < minIncrementProbabilty) {
    missedIncrementCount++;
    return;
  } 
  accumulatedProbability+= prob;
  if (!fisherKernelMode) {
    nextPmf[val] += prob;
  } else {
    // do Fisher score accumulation case.
    if (!pmf[val].zero())
      nextPmf[val] += prob/pmf[val];
  }

}


void
Dense1DPMF::emEndIteration()
{
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;

  if ( !emOnGoingBitIsSet() )
    return; 

  accumulatedProbability.floor();
  if (accumulatedProbability < minDiscAccumulatedProbability()) {
    infoMsg(IM::SoftWarning,"WARNING: Dense1DPMF named '%s' received only %e accumulated probability in EM iteration. Using previous values.",name().c_str(),accumulatedProbability.val());
    for (int i=0;i<nextPmf.len();i++) {
      nextPmf[i] = pmf[i];
    }
  } else {
    for (int i=0;i<nextPmf.len();i++) {
      nextPmf[i] /= accumulatedProbability;
    }
  }

  // stop EM
  emClearOnGoingBit();
}



void
Dense1DPMF::emSwapCurAndNew()
{
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;

  if (!emSwappableBitIsSet())
    return;

  const unsigned oldLen = nextPmf.len();
  unsigned newLen = nextPmf.len();
  unsigned numVanished = 0;
  unsigned numSplit = 0;
  for (unsigned i=0;i<(unsigned)nextPmf.len();i++) {
    if (MixtureCommon::vanishingComponentSet.
	find(pair<Dense1DPMF*,unsigned>(this,i))
	!= MixtureCommon::vanishingComponentSet.end()) {
      numVanished++;
      newLen--;
    } else if (MixtureCommon::splittingComponentSet.
	       find(pair<Dense1DPMF*,unsigned>(this,i))
	       != MixtureCommon::splittingComponentSet.end()) {
      numSplit++;
      newLen++;
    }
  }
  if (numVanished > 0 || numSplit > 0)
    warning("NOTE: DPMF '%s' has %d/%d elements vanishing/splitting",
	      name().c_str(),numVanished,numSplit);
  
  // command line check should ensure this
  assert ( newLen > 0 );

  pmf.resizeIfDifferent(newLen);

  unsigned newIndex = 0;
  for (unsigned i=0;i<(unsigned)nextPmf.len();i++) {
    if (MixtureCommon::vanishingComponentSet.
	find(pair<Dense1DPMF*,unsigned>(this,i))
	!= MixtureCommon::vanishingComponentSet.end()) {
      // do nothing, don't copy it over
      ;
    } else if (MixtureCommon::splittingComponentSet.
	       find(pair<Dense1DPMF*,unsigned>(this,i))
	       != MixtureCommon::splittingComponentSet.end()) {

      // copy it and clone over
      pmf[newIndex++] = nextPmf[i]/2.0;
      pmf[newIndex++] = nextPmf[i]/2.0;
    } else {
      pmf[newIndex++] = nextPmf[i];
    }
  }
  assert ( newIndex == newLen );

  nextPmf.resizeIfDifferent(newLen);
  
  if (useDirichletPriors) {
    if (smoothingType == DirichletTableVal) {
      if (newLen != dirichletTable->lastDimension()) {
	// if length changed, we can no longer use the table.
	// rather than die with an error, we instead turn off Dirichlet smoothing using this table.
	smoothingType = NoneVal;
	infoMsg(IM::Med,"WARNING: Dense1DPMF named '%s' changed from length %d to length %d, so is no longer using DirichletTable '%s' of length %d.",
		name().c_str(),oldLen,newLen,
		dirichletTable->name().c_str(),
		dirichletTable->lastDimension());
      }
    }
  }

  // renormalize & re-compute max value.
  normalize();
  emClearSwappableBit();
}


/*-
 *-----------------------------------------------------------------------
 *
 * Accumulator loading/storing routines for parallel training support.
 * These routines are virtual, and are called from the EMable object
 * which has the code containing the logic and checks about which one
 * to call when.
 *
 *----------------------------------------------------------------------- 
*/


void Dense1DPMF::emStoreObjectsAccumulators(oDataStreamFile& ofile,
					    bool writeLogVals,
					    bool writeZeros)
{
  if (writeZeros) {
    if (writeLogVals) {
      logpr z;
      z.set_to_zero();
      for (int i=0;i<nextPmf.len();i++) {
	ofile.write(z.val(),"DPMF store accums");
      }   
    } else {
      for (int i=0;i<nextPmf.len();i++) {
	ofile.write(0.0,"DPMF store accums");
      }   
    } 
  } else {
    if (writeLogVals) {
      for (int i=0;i<nextPmf.len();i++) {
	ofile.write(nextPmf[i].val(),"DPMF store accums");
      }
    } else {
      for (int i=0;i<nextPmf.len();i++) {
	ofile.write(nextPmf[i].unlog(),"DPMF store accums");
      }
    }
  }
}

void Dense1DPMF::emLoadObjectsDummyAccumulators(iDataStreamFile& ifile)
{
  logpr tmp;
  for (int i=0;i<pmf.len();i++) {
    ifile.read(tmp.valref(),"DPMF load accums");
  }
}

void Dense1DPMF::emZeroOutObjectsAccumulators()
{
  for (int i=0;i<nextPmf.len();i++) {
    nextPmf[i].set_to_zero();
  }
}

void Dense1DPMF::emLoadObjectsAccumulators(iDataStreamFile& ifile)
{
  for (int i=0;i<nextPmf.len();i++) {
    ifile.read(nextPmf[i].valref(),"DPMF load accums");
  }
}

void Dense1DPMF::emAccumulateObjectsAccumulators(iDataStreamFile& ifile)
{
  for (int i=0;i<nextPmf.len();i++) {
    logpr tmp;
    ifile.read(tmp.valref(),"DPMF accumulate accums");
    nextPmf[i] += tmp;
  }
}



/*-
 *-----------------------------------------------------------------------
 * Dense1DPMF::identicalIndependentClone
 *      creates an exact copy of this object that shares nothing with
 *      the original
 *
 * Preconditions:
 *      1) object being copied should be allocated
 *      2) GM_Parms should contain all parameters, so that a unique name
 *         for the new object can be generated
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      the new object is added to GM_Parms
 *
 * Results:
 *      a pointer the new object
 *
 *-----------------------------------------------------------------------
 */
Dense1DPMF* 
Dense1DPMF::identicalIndependentClone()
{
  Dense1DPMF* newDPMF = new Dense1DPMF();

  newDPMF->pmf.resize(pmf.len());
  for (int i=0;i<pmf.len();i++) 
    newDPMF->pmf[i] = pmf[i];

  newDPMF->setName(new_name(name(),&GM_Parms.dPmfsMap));
  newDPMF->setBasicAllocatedBit();
  GM_Parms.add(newDPMF);

  return newDPMF;
} 

#if 0

void
Dense1DPMF::emStoreAccumulators(oDataStreamFile& ofile)
{
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet()) {
    // then we are not training, either because
    // 1) this object obtained no probability during training or
    // 2) we have turned off training of this object.
    // In either case, we write out '0' to state that 
    // there are no values stored for this object.
    unsigned flag = 0;
    ofile.write(flag,"DPMF writing flag");
    if ( !emEmAllocatedBitIsSet() ) {
      // then we indeed have no probability values, so lets emit a warning
      warning("WARNING: zero accumulator values for DPMF '%s'\n",
	      name().c_str());
    }
    return;
  }
  // in this case, we are guaranteed that emEmAllocatedBitIsSet() is true
  // since if we are training, the em structures have at least been allocated.
  // store the accumulators as normal.
  EMable::emStoreAccumulators(ofile);
  for (int i=0;i<pmf.len();i++) {
    ofile.write(nextPmf[i].val(),"DPMF store accums");
  }
}


void
Dense1DPMF::emLoadAccumulators(iDataStreamFile& ifile)
{
  assert (basicAllocatedBitIsSet()); 
  unsigned flag;
  ifile.read(flag,"DPMF load flag");
  if (!emAmTrainingBitIsSet()) {
    // then we are not adjusting this object here.
    if (flag == 0) {
      // then not only are we not training, but
      // we also have nothing here to load, so we move on.
      return;
    } else {
      // this means that the accumulator file has data for this object
      // (which must have occured during the production of the
      // accumulator files) , but the user must have specified on the
      // final command line during the accumulate-the-accumulator
      // stage that this object should not be trained. Therefore,
      // we emit a warning, read in the accumulators into dummy locations.
      warning("WARNING: loading into dummy accumulators for fixed DPMF '%s'\n",
	      name().c_str());
      EMable::emLoadDummyAccumulators(ifile);
      logpr tmp;
      for (int i=0;i<pmf.len();i++) {
	ifile.read(tmp.valref(),"DPMF load accums");
      }
      return;
    }
  }
  // so we are training. This next assertion
  // is needed, since it shouldn't be possible that
  // we're training, but the EM data structures have not
  // been allocated.
  assert (emEmAllocatedBitIsSet());

  
  if (flag == 0) {
    // then we don't have any accumulator values for this object, but
    // still need to initizlie the accumulators.
    emZeroOutAccumulators();
    for (int i=0;i<nextPmf.len();i++) {
      nextPmf[i].set_to_zero();
    }
  } else {
    // load up the real accumulators.
    EMable::emLoadAccumulators(ifile);
    for (int i=0;i<nextPmf.len();i++) {
      ifile.read(nextPmf[i].valref(),"DPMF load accums");
    }
  }
}


void
Dense1DPMF::emAccumulateAccumulators(iDataStreamFile& ifile)
{
  assert ( basicAllocatedBitIsSet() );
  unsigned flag;
  ifile.read(flag,"DPMF load flag");
  if (!emAmTrainingBitIsSet()) {
    // then we are not adjusting this object here.
    if (flag == 0) {
      // then not only are we not training, but
      // we also have nothing here to load, so we move on.
      return;
    } else {
      // this means that the accumulator file has data for this object
      // (which must have occured during the production of the
      // accumulator files) , but the user must have specified on the
      // final command line during the accumulate-the-accumulator
      // stage that this object should not be trained. Therefore, we
      // emit a warning, read in the accumulators into dummy
      // locations.
      warning("WARNING: accumulating into dummy accumulators for fixed DPMF '%s'\n",
	      name().c_str());
      EMable::emLoadDummyAccumulators(ifile);
      logpr tmp;
      for (int i=0;i<pmf.len();i++) {
	ifile.read(tmp.valref(),"DPMF load accums");
      }
      return;
    }
  }
  // so we are training. This next assertion
  // is needed, since it shouldn't be possible that
  // we're training, but the EM data structures have not
  // been allocated.
  assert (emEmAllocatedBitIsSet());

  if (flag == 0) {
    // then we don't have any accumulator values for this object, and
    // there is nothing to do.
  } else {
    // accumulate up the real accumulators.
    EMable::emAccumulateAccumulators(ifile);
    for (int i=0;i<nextPmf.len();i++) {
      logpr tmp;
      ifile.read(tmp.valref(),"DPMF accumulate accums");
      nextPmf[i] += tmp;
    }
  }
}

#endif

////////////////////////////////////////////////////////////
// Sample generation
////////////////////////////////////////////////////////////


