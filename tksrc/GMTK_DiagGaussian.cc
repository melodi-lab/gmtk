/*-
 * GMTK_DiagGaussian.cc
 *        Code for plain vanilla diagonal Gaussians.
 *
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (C) 2001 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
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
#if HAVE_CONFIG_H
#include <config.h>
#endif
#if HAVE_HG_H
#include "hgstamp.h"
#endif
VCID(HGID)

#include "error.h"
#include "rand.h"

#include "GMTK_DiagGaussian.h"
#include "GMTK_GMParms.h"
#include "GMTK_MixtureCommon.h"
#include "tieSupport.h"



void
DiagGaussian::read(iDataStreamFile& is)
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
      error("Error: DiagGaussian '%s' in file '%s' line %d specifies mean name '%s' that does not exist",
	    _name.c_str(),is.fileName(),is.lineNo(),str.c_str());
  meanIndex = GM_Parms.meansMap[str];
  mean = GM_Parms.means[meanIndex];
  mean->numTimesShared++;


  // read covariance vector
  is.read(str);
  if (GM_Parms.covarsMap.find(str) == GM_Parms.covarsMap.end())
    error("Error: DiagGaussian '%s' in file '%s' line %d specifies covar name '%s' that does not exist",
	  _name.c_str(),is.fileName(),is.lineNo(),str.c_str());
  
  covarIndex = GM_Parms.covarsMap[str];
  covar = GM_Parms.covars[covarIndex];
  covar->numTimesShared++;

  // check that lengths match, etc.
  if (covar->dim() != mean->dim()) {
    error("Error: DiagGaussian '%s' in file '%s' line %d specifices a mean '%s' with dim %d and covariance '%s' with dim '%d'\n",
	  _name.c_str(),is.fileName(),is.lineNo(),
	  mean->name().c_str(),
	  mean->dim(),
	  covar->name().c_str(),
	  covar->dim());
  }
  if ((unsigned)covar->dim() != _dim) {
    error("Error: CondDiagGaussian '%s' in file '%s' line %d of dim %d does not match its mean '%s' with dim %d or covariance '%s' with dim '%d'\n",
	  _name.c_str(),is.fileName(),is.lineNo(),
	  _dim,
	  mean->name().c_str(),
	  mean->dim(),
	  covar->name().c_str(),
	  covar->dim());
  }

  setBasicAllocatedBit();
}


void
DiagGaussian::write(oDataStreamFile& os)
{
  assert ( basicAllocatedBitIsSet() );

  // write the type of self and the name
  os.write((int)Component::DiagGaussian);
  NamedObject::write(os);
  os.nl();

  // write mean vector
  os.write(mean->name());
  os.write(covar->name());
  os.nl();
}


/*-
 *-----------------------------------------------------------------------
 * DiagGaussian::makeRandom
 *      calls makeRandom for the mean and covar
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
DiagGaussian::makeRandom()
{
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;
  mean->makeRandom();
  covar->makeRandom();
}


/*-
 *-----------------------------------------------------------------------
 * makeUniform()
 *      calls makeUniform for the mean and covar
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
DiagGaussian::makeUniform()
{
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;

  mean->makeUniform();
  covar->makeUniform();
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
DiagGaussian::noisyClone()
{
  assert ( basicAllocatedBitIsSet() );

  map<Component*,Component*>::iterator it = MixtureCommon::mcCloneMap.find(this);
  // first check if self is already cloned, and if so, return that.
  if (it == MixtureCommon::mcCloneMap.end()) {
    DiagGaussian* clone;
    clone = new DiagGaussian(dim());

    unsigned cloneNo=0; do {
      char buff[256];
      sprintf(buff,"%d",cloneNo);
      clone->_name = _name + string("_cl") + buff;
      cloneNo++;
    } while (GM_Parms.componentsMap.find(clone->_name) 
	     != GM_Parms.componentsMap.end());

    if (cloneShareMeans && cloneShareCovars) {
      // the issue here is that if the user wants to clone a gaussian that shares
      // both mean and covar, it would be identical and no reason to clone in the first
      // place.
      warning("WARNING: Diagonal Gaussian component '%s' is cloning, and was asked to share both means and covariances. No sharing is occuring instead.",name().c_str());
      clone->mean = mean->noisyClone();
      clone->covar = covar->noisyClone();
    } else {
      if (cloneShareMeans)
	clone->mean = mean;
      else
	clone->mean = mean->noisyClone();
      if (cloneShareCovars)
	clone->covar = covar;
      else
	clone->covar = covar->noisyClone();
    }
    
    // need to tell mean, and covar that either
    //    1) if this is a new mean,covar that a
    //       a parent object is using them, or
    //    2) if this is a shared mean,covar that
    //       an additional parent object is using them.
    clone->mean->numTimesShared++;
    clone->covar->numTimesShared++;

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
 * DiagGaussian::identicalIndependentClone
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
DiagGaussian::identicalIndependentClone()
{

  DiagGaussian* newDG = new DiagGaussian(dim());

  newDG->mean = mean->identicalIndependentClone();
  newDG->covar = covar->identicalIndependentClone();

  // don't change usage counts here - do it in calling function,
  // because only that knows how the sharing is arranged

  //newDG->mean->numTimesShared++;
  //newDG->covar->numTimesShared++;

  //mean->numTimesShared--;
  //covar->numTimesShared--;

  newDG->_name = new_name(name(),&GM_Parms.componentsMap);
  newDG->setBasicAllocatedBit();
  GM_Parms.add(newDG);

  return newDG;
}

/////////////////
// EM routines //
/////////////////



void
DiagGaussian::emStartIteration()
{
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;

  if (emOnGoingBitIsSet())
    return;

  if (!emEmAllocatedBitIsSet()) {
    // this is presumably the first time
    trMembers.allocateIfNeeded();
    emSetEmAllocatedBit();
  }

  // EM iteration is now going.
  emSetOnGoingBit();
  emSetSwappableBit();

  accumulatedProbability = 0.0;
  mean->emStartIteration(trMembers->nextMeans);
  covar->emStartIteration(trMembers->nextDiagCovars);
}


void
DiagGaussian::emIncrement(logpr prob,
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
  const float fprob = prob.unlog();


  // these next calls are presumed to be valid for both EM and FK (fisher kernel accumulation).
  mean->emIncrement(prob,fprob,f,base,stride,trMembers->nextMeans.ptr);
  covar->emIncrement(prob,fprob,f,base,stride,trMembers->nextDiagCovars.ptr);

  if (!fisherKernelMode) {
    // do normal EM increment
    // this call is optimized to do the 1st and 2nd moment stats simultaneously
    emIncrementMeanDiagCovar(fprob,f,trMembers->nextMeans.size(),trMembers->nextMeans.ptr,trMembers->nextDiagCovars.ptr);
  } else {
    // do a fisher score increment
    fkIncrementMeanDiagCovar(fprob,f,trMembers->nextMeans.size(),
			     mean->means.ptr,
			     covar->covariances.ptr,
			     trMembers->nextMeans.ptr,
			     trMembers->nextDiagCovars.ptr);

  }
}

void
DiagGaussian::emEndIteration()
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
    // Note: we assume here that the mean and covar object will check for us that 
    // its accumulated probability is above threshold. Here, we just
    // check for zero, and then issue a warning if needed. This might not
    // indicate a problem as the mean and covar of this object might be shared
    // and might have received plenty of count.
    warning("WARNING: Gaussian Component named '%s' did not receive any accumulated probability in EM iteration. Global missed increment count is %d. Also check child mean '%s' and covar '%s'",
	    name().c_str(),
	    missedIncrementCount,
	    mean->name().c_str(),
	    covar->name().c_str());
  }


  if (covar->emSharedBitIsSet()) {
    // covariance is shared
    if (mean->emSharedBitIsSet()) {
      // mean and covariance are both shared
      mean->emEndIterationSharedMeansCovars
	(accumulatedProbability,trMembers->nextMeans.ptr,covar);
      covar->emEndIterationSharedMeansCovars
	(accumulatedProbability,
	 trMembers->nextMeans.ptr,
	 trMembers->nextDiagCovars.ptr,
	 mean);
    } else {
      // Mean is not shared, but covariance is shared
      // can still use normal EM
      mean->emEndIterationNoSharing(trMembers->nextMeans.ptr);
      covar->emEndIterationSharedCovars(accumulatedProbability,
					trMembers->nextMeans.ptr,
					trMembers->nextDiagCovars.ptr);
    }
  } else {
    // covariance is not shared
    if (mean->emSharedBitIsSet()) {
      // covariance is not shared, mean is shared,
      // but use the "all shared" versions since we don't
      // have the new mean at this point.
      mean->emEndIterationSharedMeansCovars
	(accumulatedProbability,trMembers->nextMeans.ptr,covar);
      covar->emEndIterationSharedMeansCovars
	(accumulatedProbability,
	 trMembers->nextMeans.ptr,
	 trMembers->nextDiagCovars.ptr,
	 mean);
    } else {
      // nothing is shared, use EM
      mean->emEndIterationNoSharing(trMembers->nextMeans.ptr);
      covar->emEndIterationNoSharing(trMembers->nextMeans.ptr,trMembers->nextDiagCovars.ptr);
    }
  }

  // mean->emEndIteration(nextMeans.ptr);
  // covar->emEndIteration(accumulatedProbability,nextMeans.ptr,nextDiagCovars.ptr);

  emClearOnGoingBit();
}


void
DiagGaussian::emSwapCurAndNew()
{
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;

  if (!emSwappableBitIsSet())
    return;

  mean->emSwapCurAndNew();
  covar->emSwapCurAndNew();

  emClearSwappableBit();
}


////////////////////////////////////////////////////////////
// Parallel EM support
////////////////////////////////////////////////////////////


void
DiagGaussian::emStoreObjectsAccumulators(oDataStreamFile& ofile,
					 bool writeLogVals,
					 bool writeZeros)
{
  // since this is a Gaussian, we ignore the writeLogVals
  // argument since it doesn't make sense to take log of
  // these values since they are continuous, could be negative, etc.
  if (writeZeros) {
    const unsigned totalLen = trMembers->nextMeans.len() + trMembers->nextDiagCovars.len();
    for (unsigned i=0;i<totalLen;i++) {
      ofile.write(0.0,"Diag Gaussian store accums nm + nc.");
    }
  } else {
    for (int i=0;i<trMembers->nextMeans.len();i++) {
      ofile.write(trMembers->nextMeans[i],"Diag Gaussian store accums nm.");
    }
    for (int i=0;i<trMembers->nextDiagCovars.len();i++) {
      ofile.write(trMembers->nextDiagCovars[i],"Diag Gaussian store accums nc.");
    }
  }
}


void
DiagGaussian::emLoadObjectsDummyAccumulators(iDataStreamFile& ifile)
{
  //
  // ASSUME MEANS AND COVARIANCES ARE OF TYPE FLOAT
  // See the MeanVector.h and DiagcovarVector.h for specifics.
  float tmp;
  for (int i=0;i<mean->dim();i++) {
    ifile.read(tmp,"Diag Gaussian load accums nm.");
  }
  for (int i=0;i<covar->dim();i++) {
    ifile.read(tmp,"Diag Gaussian load accums nc.");
  }
}


void
DiagGaussian::emZeroOutObjectsAccumulators()
{
  for (int i=0;i<trMembers->nextMeans.len();i++) {
    trMembers->nextMeans[i] = 0.0;
  }
  for (int i=0;i<trMembers->nextDiagCovars.len();i++) {
    trMembers->nextDiagCovars[i] = 0.0;
  }
}

void
DiagGaussian::emLoadObjectsAccumulators(iDataStreamFile& ifile)
{
  for (int i=0;i<trMembers->nextMeans.len();i++) {
    ifile.read(trMembers->nextMeans[i],"Diag Gaussian load accums nm.");
  }
  for (int i=0;i<trMembers->nextDiagCovars.len();i++) {
    ifile.read(trMembers->nextDiagCovars[i],"Diag Gaussian load accums nc.");
  }
}


void
DiagGaussian::emAccumulateObjectsAccumulators(iDataStreamFile& ifile)
{
  //
  // ASSUME MEANS AND COVARIANCES ARE OF TYPE FLOAT
  // See the MeanVector.h and DiagcovarVector.h for specifics.
  for (int i=0;i<trMembers->nextMeans.len();i++) {
    float tmp;
    ifile.read(tmp,"Diag Gaussian accumulate accums nm.");
    trMembers->nextMeans[i] += tmp;
  }
  for (int i=0;i<trMembers->nextDiagCovars.len();i++) {
    float tmp;
    ifile.read(tmp,"Diag Gaussian accumulate accums nc.");
    trMembers->nextDiagCovars[i] += tmp;
  }
}


////////////////////////////////////////////////////////////
// Sample generation
////////////////////////////////////////////////////////////



void DiagGaussian::sampleGenerate(float *const sample,
				  const Data32* const base,
				  const int stride)
{
  error("not implemented");
}


