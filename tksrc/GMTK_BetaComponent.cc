/*-
 * GMTK_BetaComponent.cc
 *        Code for Beta observation distribution.
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2008, < fill in later >
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

#include "GMTK_BetaComponent.h"
#include "GMTK_GMParms.h"
#include "GMTK_MixtureCommon.h"
#include "tieSupport.h"


////////////////////////////////////////////////////
// The default value of the minimum possible variance of any
// Beta Distribution. This must be >= FLT_MIN for numeric stability,
// and it should be made availalbe as a command line parameter at some point.
double BetaComponent::_varianceFloor = BetaComponent::setVarianceFloor(1e-10);

double BetaComponent::setVarianceFloor(const double floor) 
{ 
  if (floor < FLT_MIN) 
    _varianceFloor = FLT_MIN; 
  else 
    _varianceFloor = floor; 
  return _varianceFloor;
}


// TODO: these values should not be exported to the command line
// until (and if) sharing is done.
// true if when a clone occurs, we use the same alpha
// (i.e., share the alpha and clone other things)
bool BetaComponent::cloneShareAlpha = false;
// true if when a clone occurs, we use the same beta (i.e.,
// share the beta and only clone other things).
bool BetaComponent::cloneShareBeta = false;


/*-
 *-----------------------------------------------------------------------
 * recomputeNormalizer()
 *      Recompute pre-computed values for probabilit evaluation.
 * 
 * Preconditions:
 *      Basic alpha and beta parameter matrices should be allocated.
 *
 * Postconditions:
 *      normalizer is recomputed. 
 *
 * Side Effects:
 *      modifies variable 'normalizer' and the alpha/beta parameters.
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */

void
BetaComponent::recomputeNormalizer()
{
  assert ( basicAllocatedBitIsSet() );
  normalizer.resizeIfDifferent(_dim);
  for (unsigned i=0; i< _dim; i++) {
    normalizer.ptr[i] =  lgamma( alpha->values.ptr[i] + beta->values.ptr[i] )
      - lgamma ( alpha->values.ptr[i] ) - lgamma ( beta->values.ptr[i] );
  }
}




/*-
 *-----------------------------------------------------------------------
 * log_p()
 *      Computes the probability of this Beta distribution vector (multiplying them together)
 * 
 * Preconditions:
 *      none
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
BetaComponent::log_p(const float *const x,
		    const Data32* const base,
		    const int stride)
{
  assert ( basicAllocatedBitIsSet() );

  double val = 0;
  for (unsigned i = 0; i < _dim ; i++ ) {

    // enforce strict limits
    if (x[i] <= lower || x[i] >= upper) {
      // unfortnately there is no way here to give a better error
      // message.  I suppose if we threw an exception, it could be
      // caught where more information is availale regarding the
      // observation location and file name, etc.
      error("ERROR: obsevation value is %f (element %d), but Beta distribution %s with range (lower,upper)=(%f,%f) cannot have such a value.",x[i],i,_name.c_str(),lower,upper);
    }

    double x_norm = rangeScale*(x[i] - lower);
    
    // Compute the probability in the log domain. 
    val += normalizer.ptr[i] + 
      (alpha->values.ptr[i]-1.0)*log(x_norm) + (beta->values.ptr[i]-1.0)*log(1 - x_norm);

  }
  return logpr(0,val);
}


// Slightly perturb either the alpha & beta parameters.
void
BetaComponent::perturbAlpha()
{
  assert ( basicAllocatedBitIsSet() );
  for (unsigned i=0; i< _dim ; i++) {
    alpha->values.ptr[i] +=  (alpha->values.ptr[i]/8.0)*(rnd.drand48()-0.5); 
  }  
}
void
BetaComponent::perturbBeta()
{
  assert ( basicAllocatedBitIsSet() );
  for (unsigned i=0;i<_dim;i++) {
    beta->values.ptr[i] +=  (beta->values.ptr[i]/8.0)*(rnd.drand48()-0.5); 
  }
}




void
BetaComponent::read(iDataStreamFile& is)
{
  // The index in the global alpha 
  int alphaIndex; 
  // The index in the global beta
  int betaIndex;

  // read name
  NamedObject::read(is);

  // read range (lower,upper)
  is.read(lower,"beta lower");
  is.read(upper,"beta upper");  
  if (upper <= lower) {
    error("Error: BetaComponent '%s' in file '%s' line %d specifies invalid (lower,upper)=(%f,%d) range",
	  _name.c_str(),is.fileName(),is.lineNo(),lower,upper);
  }
  rangeScale = 1.0/(upper - lower);

  // read alpha vector
  string str;
  is.read(str);

  if (GM_Parms.realMatsMap.find(str) ==  GM_Parms.realMatsMap.end()) 
      error("Error: BetaComponent '%s' in file '%s' line %d specifies alpha name '%s' that does not exist",
	    _name.c_str(),is.fileName(),is.lineNo(),str.c_str());
  alphaIndex = GM_Parms.realMatsMap[str];
  alpha = GM_Parms.realMats[alphaIndex];

  if (alpha->rows() != 1 || (unsigned)alpha->cols() != _dim)
    error("Error: BetaComponent '%s' in file '%s' line %d specifices a alpha real matrix '%s' with rows,cols = %d, %d but rows must be 1 and columns must be %d.\n",
	  _name.c_str(),is.fileName(),is.lineNo(),
	  alpha->name().c_str(),
	  alpha->rows(),
	  alpha->cols(),
	  _dim);


  // read beta vector
  is.read(str);
  if (GM_Parms.realMatsMap.find(str) == GM_Parms.realMatsMap.end())
    error("Error: BetaComponent '%s' in file '%s' line %d specifies beta name '%s' that does not exist",
	  _name.c_str(),is.fileName(),is.lineNo(),str.c_str());
  
  betaIndex = GM_Parms.realMatsMap[str];
  beta = GM_Parms.realMats[betaIndex];

  if (beta->rows() != 1 || (unsigned)beta->cols() != _dim)
    error("Error: BetaComponent '%s' in file '%s' line %d specifices a beta real matrix '%s' with rows,cols = %d, %d but rows must be 1 and columns must be %d.\n",
	  _name.c_str(),is.fileName(),is.lineNo(),
	  beta->name().c_str(),
	  beta->rows(),
	  beta->cols(),
	  _dim);

  setBasicAllocatedBit();
  recomputeNormalizer();


}


void
BetaComponent::write(oDataStreamFile& os)
{
  assert ( basicAllocatedBitIsSet() );

  // write the type of self and the name
  os.write((int)Component::BetaComponent);
  NamedObject::write(os);
  os.nl();
  os.write(lower);  os.write(upper); os.nl();

  // write object names
  os.write(alpha->name());
  os.write(beta->name());
  os.nl();
}


/*-
 *-----------------------------------------------------------------------
 * BetaComponent::makeRandom
 *      calls makeRandom for the parameters
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
BetaComponent::makeRandom()
{
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;
  for (unsigned i = 0; i < _dim ; i++ ) {
    // start with values > 1.0
    alpha->values.ptr[i] = 1.0+rnd.drand48pe();
    beta->values.ptr[i] = 1.0+rnd.drand48pe();
  }
}


/*-
 *-----------------------------------------------------------------------
 * makeUniform()
 *      calls makeUniform for the parameters.
 * 
 * Preconditions:
 *      Object must be allocated.
 *
 * Postconditions:
 *      Object has uniform values.
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
BetaComponent::makeUniform()
{
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;

  for (unsigned i = 0; i < _dim; i++ ) {
    alpha->values.ptr[i] = 1.0;
    beta->values.ptr[i] = 1.0;
  }

}


/*-
 *-----------------------------------------------------------------------
 * noisyClone()
 *      Create a copy of self, but with cloned perturbed the alpha/beta vectors
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
 *      returns the new noisy distribution.
 *
 *-----------------------------------------------------------------------
 */
Component*
BetaComponent::noisyClone()
{
  assert ( basicAllocatedBitIsSet() );

  map<Component*,Component*>::iterator it = MixtureCommon::mcCloneMap.find(this);

  // first check if self is already cloned, and if so, return that.
  if (it == MixtureCommon::mcCloneMap.end()) {
    BetaComponent* clone;
    clone = new BetaComponent(_dim);

    clone->lower = lower;
    clone->upper = upper;
    clone->rangeScale = rangeScale;

    unsigned cloneNo=0; do {
      // TODO: change not to use finite length buffer like this.
      char buff[256];
      sprintf(buff,"%d",cloneNo);
      clone->_name = _name + string("_cl") + buff;
      cloneNo++;
    } while (GM_Parms.componentsMap.find(clone->_name) 
	     != GM_Parms.componentsMap.end());

    // sharing during training is not implemented (yet). Make sure the
    // user doesn't change code above and set these values.
    assert (!cloneShareBeta && !cloneShareAlpha);

    if (cloneShareBeta && cloneShareAlpha) {
      warning("WARNING: Beta Component '%s' is cloning, and was asked to share both alpha and beta. No sharing is occuring instead.",name().c_str());
      clone->alpha = alpha->cleanClone();
      clone->beta = beta->cleanClone();
      clone->perturbAlpha();
      clone->perturbBeta();
    } else {
      if (cloneShareBeta)
	clone->beta = beta;
      else {
	clone->beta = beta->cleanClone();
	clone->perturbBeta();
      }

      if (cloneShareAlpha)
	clone->alpha = alpha;
      else {
	clone->alpha = alpha->cleanClone();
	clone->perturbAlpha();
      }
    }
    
    // need to tell alpha and beta that either
    //    1) if this is a new alpha,beta that a
    //       a parent object is using them, or
    //    2) if this is a shared alpha,beta that
    //       an additional parent object is using them.
    clone->alpha->numTimesShared++;
    clone->beta->numTimesShared++;

    clone->setBasicAllocatedBit();
    MixtureCommon::mcCloneMap[this] = clone;

    // also add self to GMParms object.
    GM_Parms.add(clone);

    return clone;
  } else {
    return (*it).second;
  }
}



/*-
 *-----------------------------------------------------------------------
 * BetaComponent::identicalIndependentClone
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
Component* 
BetaComponent::identicalIndependentClone()
{
  BetaComponent* newDG = new BetaComponent(_dim);

  newDG->beta = beta->cleanClone();
  newDG->alpha = alpha->cleanClone();
  newDG->lower = lower;
  newDG->upper = upper;
  newDG->rangeScale = rangeScale;

  // don't change usage counts here - do it in calling function,
  // because only that knows how the sharing is arranged

  newDG->_name = new_name(name(),&GM_Parms.componentsMap);
  newDG->setBasicAllocatedBit();
  GM_Parms.add(newDG);

  return newDG;
}

/////////////////
// EM routines //
/////////////////



void
BetaComponent::emStartIteration()
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

  sumx.resizeIfDifferent(_dim);
  sumxx.resizeIfDifferent(_dim);
  for (unsigned i = 0 ; i < _dim; i ++ ) {
    sumx.ptr[i] = sumxx.ptr[i] = 0;

  }

  alpha->emStartIteration();
  beta->emStartIteration();

}


void
BetaComponent::emIncrement(logpr prob,
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
    // floating point number, it would be a denomral.
  }

  accumulatedProbability += prob;
  // prob.unlog() here so it doesn't need to be done
  // twice by the callees.
  const double fprob = prob.unlog();

  // we next do an estimate of the stats for a Beta distribuion.
  // NOTE: this is currently just a hack, this is *NOT* the
  // MLE. To do an MLE of a Beta requires an additional iterative
  // procedure.

  if (!fisherKernelMode) {
    // do normal EM increment
    // this call is optimized to do the 1st and 2nd moment stats simultaneously
    for (unsigned i = 0; i < _dim; i++ ) {
      // enforce strict positivity.
      if (f[i] <= 0.0) {
	// unfortnately there is no way here to give a better error
	// message.  I suppose if we threw an exception, it could be
	// caught where more information is availale regarding the
	// observation location and file name, etc.
	error("ERROR: em training, obsevation value is %f, but a Beta distribution cannot have a non-positive observation.",f[i]);
      }
      sumx.ptr[i] += f[i]*fprob;
      sumxx.ptr[i] += f[i]*f[i]*fprob;
    }
  } else {
    error("ERROR: fisher kernel not yet implemented with beta components");
  }

  alpha->emIncrement(prob);
  beta->emIncrement(prob);
}

void
BetaComponent::emEndIteration()
{

  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;

  // if EM not ongoing, we just return
  // here since, if this object is shared,
  // we don't want to end the epoch > 1 time
  if (!emOnGoingBitIsSet())
    return; 

  accumulatedProbability.floor();
  if (accumulatedProbability.zero()) {
    // Note: we assume here that the alpha/beta object will check for us that 
    // its accumulated probability is above threshold. Here, we just
    // check for zero, and then issue a warning if needed. This might not
    // indicate a problem as the alpha/beta of this object might be shared
    // and might have received plenty of count (TODO: implement sharing).
    warning("WARNING: Beta Component named '%s' did not receive any accumulated probability in EM iteration. Global missed increment count is %d. Also check child alpha '%s' and beta '%s'",
	    name().c_str(),
	    missedIncrementCount,
	    alpha->name().c_str(),
	    beta->name().c_str());
  }

  // TODO: implement the sharing case properly.
  const double inv_denom = accumulatedProbability.inverse().unlog();

  for (unsigned i = 0; i < _dim; i++ ) {
    double mean = sumx.ptr[i]*inv_denom;
    double variance = sumxx.ptr[i]*inv_denom - mean*mean;

    // normalize for range.
    mean = (mean - lower)*rangeScale;
    variance = variance*rangeScale*rangeScale;

    if (variance < _varianceFloor)
      variance = _varianceFloor;

    double tmp = mean * (1.0 - mean )/variance - 1.0;
    beta->nextValues.ptr[i] = (1.0 - mean ) * tmp;
    alpha->nextValues.ptr[i] = mean * tmp;

  }

  alpha->emEndIteration();
  beta->emEndIteration();
  emClearOnGoingBit();

}


void
BetaComponent::emSwapCurAndNew()
{
  assert ( basicAllocatedBitIsSet() );

  if (!emAmTrainingBitIsSet())
    return;

  if (!emSwappableBitIsSet())
    return;

  alpha->emSwapCurAndNew();
  beta->emSwapCurAndNew();
  recomputeNormalizer();

  emClearSwappableBit();
}


////////////////////////////////////////////////////////////
// Parallel EM support
////////////////////////////////////////////////////////////


void
BetaComponent::emStoreObjectsAccumulators(oDataStreamFile& ofile,
					 bool writeLogVals,
					 bool writeZeros)
{
  assert (emEmAllocatedBitIsSet());
    

  // since this is a Beta, we ignore the writeLogVals
  // argument since it doesn't make sense to take log of
  // these values since they are continuous. etc.
  if (writeZeros) {
    for (unsigned i=0;i<2*_dim;i++) {
      ofile.write(0.0,"Beta Component store accum.");
    }
  } else {
    for (unsigned i = 0; i < _dim; i++) {
      ofile.write(sumx.ptr[i],"Beta Component store accum.");
      ofile.write(sumxx.ptr[i],"Beta Component store accum.");
    }
  }
}


void
BetaComponent::emLoadObjectsDummyAccumulators(iDataStreamFile& ifile)
{

  // assume alpha and beta are of type float
  for (unsigned i = 0; i < _dim; i++) {
    float tmp;
    ifile.read(tmp,"Beta load accums al.");
    ifile.read(tmp,"Beta load accums bt.");
  }
}


void
BetaComponent::emZeroOutObjectsAccumulators()
{
  assert (emEmAllocatedBitIsSet());
  for (unsigned i = 0; i < _dim; i++) {
    sumx.ptr[i] = sumxx.ptr[i] = 0;
  }
}

void
BetaComponent::emLoadObjectsAccumulators(iDataStreamFile& ifile)
{
  assert (emEmAllocatedBitIsSet());
  for (unsigned i = 0; i < _dim; i++) {
    ifile.read(sumx.ptr[i],"Beta Component load accums x");
    ifile.read(sumxx.ptr[i],"Beta Component load accums xx.");
  }
}


void
BetaComponent::emAccumulateObjectsAccumulators(iDataStreamFile& ifile)
{
  //
  // assume alphas and betas are of type float
  assert (emEmAllocatedBitIsSet());
  for (unsigned i = 0; i < _dim; i++) {
    float tmp; // here is the float assumption.
    ifile.read(tmp,"Beta component accumulate accums sumx.");
    sumx.ptr[i] += tmp;
    ifile.read(tmp,"Diag Gaussian accumulate accums sumxx.");
    sumxx.ptr[i] += tmp;
  }
}


////////////////////////////////////////////////////////////
// Sample generation
////////////////////////////////////////////////////////////



void BetaComponent::sampleGenerate(float *const sample,
				  const Data32* const base,
				  const int stride)
{
  error("BetaComponent sample generate not implemented");
}

