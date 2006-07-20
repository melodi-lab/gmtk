/*-
 * GMTK_EMable.cc
 *     A vector used for means of Gaussians.
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
#include "logp.h"

#include "GMTK_EMable.h"
#include "GMTK_GMParms.h"

VCID("$Header$")

logpr
// EMable::minIncrementProbabilty = logpr((void*)NULL,log_FLT_MIN);
EMable::minIncrementProbabilty = logpr((void*)NULL,log_FLT_MIN);

unsigned long
EMable::missedIncrementCount = 0;

bool EMable::useDirichletPriors = false;


////////////////////////////////////////////////////
// The minimum accumulated probability of mean and covariance -like
// objects. If the accumulated probability falls below this
// value, then the mean or variance like object will not
// update its values.
logpr EMable::_minContAccumulatedProbability = 
EMable::setMinContAccumulatedProbability(logpr((void*)NULL, (double)-600.0));
logpr EMable::setMinContAccumulatedProbability(const logpr floor) 
{ 
  // We hard limit to be no less than exp(-700).
  // we use -700 because if it is much smaller than this
  // we won't be able to take the exp of the inverse (i.e.,
  // exp(700) is to large to represent in non-log space.
  // We need to do this, however, because continuous parameters
  // need to be updated not in log probability.
  if (floor.val() < -700)
    _minContAccumulatedProbability = logpr((void*)NULL, (double)-700.0);
  else 
    _minContAccumulatedProbability = floor; 
  return _minContAccumulatedProbability;
}

logpr EMable::_minDiscAccumulatedProbability =
EMable::setMinDiscAccumulatedProbability(logpr((void*)NULL, (double)LSMALL/2.0));

logpr EMable::setMinDiscAccumulatedProbability(const logpr floor) 
{ 
  // Hard limit to be no less than LSMALL, the smallest 
  // possible non-zero log probability.
  if (floor.val() < LSMALL)
    _minDiscAccumulatedProbability = logpr((void*)NULL, (double)LSMALL);
  else 
    _minDiscAccumulatedProbability = floor; 
  return _minDiscAccumulatedProbability;
}




////////////////////////////////////////////////////////////////////
//        Support for parallel training/accumulator load/store.
////////////////////////////////////////////////////////////////////



/*-
 *-----------------------------------------------------------------------
 *
 * Accumulator loading/storing routines.
 *
 *-----------------------------------------------------------------------
 */


/*
 * emStoreAccumulators:
 *   store the accumulators to file ofile.
 *   we perform very simple compression (i.e., don't write long vectors of zeros) on the output.
 *   When possible, write the accumulators in log probability space.
 *
 */
void
EMable::emStoreAccumulators(oDataStreamFile& ofile) 
{
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet()) {
    // then we are not training, because
    // we have turned off training of this object.
    // We write out '0' to state that 
    // there are no values stored for this object.
    unsigned flag = 0;
    ofile.write(flag,"writing acc flag");
    return;
  } else {
    // the training bit is set.
    if (accumulatedProbability.zero()) {
      // then we are training, but indeed have no probability values, so lets emit a warning
      infoMsg(IM::SoftWarning,"WARNING: zero accumulator values for %s '%s'\n",
	      typeName().c_str(),
	      name().c_str());
      // We write out a special flag '0' to state that there are no
      // values stored for this object. This saves space compared to
      // writing out all the accumulators with zero values. This also
      // makes parameters and accumulator files compatible between
      // different structures (so we can use different structures to
      // train different parts of the parameters, or in different
      // ways).
      unsigned flag = 0;
      ofile.write(flag,"writing acc flag");
    } else {
      // if we have accumulated probability, em structures
      // should have been allocated. This ensures
      // that when we write a 1 before the accumulators,
      // we will always have the full set of accumulators written.
      assert (emEmAllocatedBitIsSet());

      // We first write a 1 to indicate that there are real
      // accumulators stored for this object.
      unsigned flag = 1;
      ofile.write(flag,"writing acc flag");

      // store the accumulators as normal.
      ofile.write(accumulatedProbability.val(),"EM store accums");
      // call virtual function to do actual work for object.
      emStoreObjectsAccumulators(ofile);
    }
  }
}


void
EMable::emLoadAccumulators(iDataStreamFile& ifile)
{
  assert (basicAllocatedBitIsSet()); 
  // first read the stored flag.
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
      infoMsg(IM::Warning,"WARNING: loading into dummy accumulators for fixed %s '%s'\n",
	      typeName().c_str(),
	      name().c_str());
      // EMable::emLoadDummyAccumulators(ifile);
      // load a dummy accummulator
      logpr tmp;
      ifile.read(tmp.valref(),"EM load accums");
      // call virtual function to do actual work for object.
      emLoadObjectsDummyAccumulators(ifile);
      return;
    }
  } else {
    // so we are training. 
  
    if (flag == 0) {
      // then we don't have any accumulator values for this object, but
      // still need to initialize the accumulators.
      // EMable::emZeroOutAccumulators();
      accumulatedProbability.set_to_zero();
      // call virtual function to do actual work for object.
      emZeroOutObjectsAccumulators();
    } else {

      // This next assertion is needed, since it shouldn't be possible
      // that we're training and have non-zero accumulation, but the
      // EM data structures have not been allocated.
      assert (emEmAllocatedBitIsSet());

      // load up the real accumulators.
      // EMable::emLoadAccumulators(ifile);
      ifile.read(accumulatedProbability.valref(),"EM load accums");
      // call virtual function to do actual work for object.
      emLoadObjectsAccumulators(ifile);
    }
  }
}


void
EMable::emAccumulateAccumulators(iDataStreamFile& ifile)
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
      infoMsg(IM::Warning,"WARNING: accumulating into dummy accumulators for fixed %s '%s'\n",
	      typeName().c_str(),
	      name().c_str());
      // EMable::emLoadDummyAccumulators(ifile);
      // load a dummy accummulator
      logpr tmp;
      ifile.read(tmp.valref(),"EM load accums");
      // call virtual function to do actual work for object.
      emLoadObjectsDummyAccumulators(ifile);
      return;
    }
  } else {
    // so we are training. 

    if (flag == 0) {
      // then we don't have any accumulator values for this object, and
      // there is nothing to do.
    } else {
      // This next assertion is needed, since it shouldn't be possible
      // that we're training and have non-zero accumulation, but the EM
      // data structures have not been allocated.
      assert (emEmAllocatedBitIsSet());

      // accumulate up the real accumulators.
      // EMable::emAccumulateAccumulators(ifile);
      logpr tmp;
      ifile.read(tmp.valref(),"EM accumulate accums");
      accumulatedProbability += tmp;
      // call virtual function to do actual work for object.
      emAccumulateObjectsAccumulators(ifile);
    }
  }
}


void
EMable::emInitAccumulators()
{
  assert (basicAllocatedBitIsSet()); 

  if (!emAmTrainingBitIsSet()) {
    // we are not "training" here, so we just move on.
    return;
  } else {
    // so we are training or using these accumulators. 
    accumulatedProbability.set_to_zero();
    // call virtual function to do actual work for object.
    emZeroOutObjectsAccumulators();
  }
}


/*
 * emWriteBasicAccumulators:
 *  write accumulators in a simple, easy, unencoded, and fixed
 *  length way to read by humans, or to be used as a input to feature
 *  transform for a kernel machine. Note, we do not write out
 *  log probabilities if logp is currently using them.
 *
 * 
 */

void
EMable::emWriteUnencodedAccumulators(oDataStreamFile& ofile,
				     bool writeLogVals)
{
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet()) {
    // then we are not training, because
    // we have turned off training of this object.
    // We write nothing out.
    return;
  } else {

    // store the accumulators as normal values      
    if (writeLogVals) {
      ofile.write(accumulatedProbability.val(),"EM store accums");
    } else {
      ofile.write(accumulatedProbability.unlog(),"EM store accums");
    }
    if (accumulatedProbability.zero()) {
      // then we have no probability values, so we need to write out 'zero' accumulators.
      // call virtual function to do actual work for object.
      emStoreObjectsAccumulators(ofile,
				 writeLogVals, // log value writing
				 true // set to true to just write zero vector of correct len
				 );
    } else {
      // the training bit is set.  EM structures should have been
      // allocated. This ensures
      // that when we write a 1 before the accumulators,
      // we will always have the full set of accumulators written.
      assert (emEmAllocatedBitIsSet());

      // call virtual function to do actual work for object.
      emStoreObjectsAccumulators(ofile,
				 writeLogVals // log value writing
				 );
    }
  }
}





