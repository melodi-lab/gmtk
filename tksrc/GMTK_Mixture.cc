/*-
 * GMTK_Mixture.cc
 *
 *        Code for mixtures of components, where each component can be
 *        one of a variety of different types (e.g., diagonal covariance Gaussian, 
 *        linear conditional Gaussian, sparse inverse covariance Gaussian, etc.).
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
VCID("$Header$");
#include "error.h"
#include "rand.h"
#include "debug.h"

#include "GMTK_Mixture.h"
#include "GMTK_GMParms.h"
#include "GMTK_ObservationMatrix.h"



void
Mixture::read(iDataStreamFile& is)
{
  // read name of self
  NamedObject::read(is);
  
  // read number of mixtures
  is.read(numComponents,"Mixture::read numComponents");
  if (numComponents <= 0) 
    error("Error: Mixture '%s' has negative or zero number of components = %d\n",
	  _name.c_str(),numComponents);

  // read name of dense 1d PMF to use for mixtures
  string str;
  is.read(str);
  if (GM_Parms.dPmfsMap.find(str) == GM_Parms.dPmfsMap.end()) {
    error("ERROR: Mixture '%s' in file '%s', can't find PMF named '%s'\n",
	  _name.c_str(),
	  is.fileName(),str.c_str());
  }

  dense1DPMF = GM_Parms.dPmfs[
			      GM_Parms.dPmfsMap[str]
  ];

  // now make sure that this one matches the number of components.
  if (numComponents != dense1DPMF->length()) {
    error("ERROR: Mixture '%s' in file '%s', PMF named '%s' has %d elements but we need %d\n",
	  _name.c_str(),is.fileName(),
	  str.c_str(),dense1DPMF->length(),numComponents);
  }

  // now read 'name' pointers to all the Components
  components.resize(numComponents);
  for (unsigned i=0;i<numComponents;i++) {
    is.read(str);
    if (GM_Parms.componentsMap.find(str) == GM_Parms.componentsMap.end()) {
      error("ERROR: Mixture '%s' in file '%s', can't find Component named '%s'\n",_name.c_str(),is.fileName(),str.c_str());
    }
    Component*gc = GM_Parms.components [
	GM_Parms.componentsMap[str]
    ];
    components[i] = gc;
    if (gc->dim() != dim()) {
      error("ERROR: Mixture '%s' in file '%s' of dim %d trying to use component '%s' of dim %d\n",
	    name().c_str(),is.fileName(),dim(),gc->name().c_str(),gc->dim());
    }
  }

  // make ready for probability evaluation.
  componentCache.resize(10);
  setBasicAllocatedBit();
}

void
Mixture::write(oDataStreamFile& os)
{
  assert ( basicAllocatedBitIsSet() );

  NamedObject::write(os);
  os.nl();
  // read number of mixture components
  os.write(numComponents,"Mixture::write numComponents");

  // write name of dense 1d PMF to use for mixtures
  os.write(dense1DPMF->name()); os.nl();

  // and the component names
  for (unsigned i=0;i<numComponents;i++) {
    os.write(components[i]->name());
  }
  os.nl();
}



/*-
 *-----------------------------------------------------------------------
 * totalNumberParameters()
 *      return total number of parameters used by this mixture.
 * 
 * Preconditions:
 *      none
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
unsigned Mixture::totalNumberParameters()
{
  unsigned sum=0;
  for (unsigned i=0;i<components.size();i++) 
    sum += components[i]->totalNumberParameters();
  return (sum + dense1DPMF->totalNumberParameters());
}


void
Mixture::makeUniform()
{
  if (!emAmTrainingBitIsSet())
    return;

  dense1DPMF->makeUniform();
  for (unsigned i=0;i<numComponents;i++) {
    components[i]->makeUniform();
  }
}


void
Mixture::makeRandom()
{
  if (!emAmTrainingBitIsSet())
    return;

  dense1DPMF->makeRandom();
  for (unsigned i=0;i<numComponents;i++) {
    components[i]->makeRandom();
  }
}


/*-
 *-----------------------------------------------------------------------
 * emptyComponentCache()
 *      For all allocated entries in the component cache, set the
 *      cache so that the element is empty.
 * 
 * Preconditions:
 *      none (i.e., component cache can be allocated or deallocated).
 *
 * Postconditions:
 *      component cache is essentially empty (i.e., mixture probabilities will be re-computed)
 *
 * Side Effects:
 *      changes component cache status, but does no memory allocation.
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
void
Mixture::emptyComponentCache()
{
  for (unsigned i=0;i<componentCache.size();i++) {
    componentCache[i].prob.valref() = (-LZERO);
  }
}

//
// compute the log probability of x with stride 'stride'
// 
logpr
Mixture::log_p(const float *const x,
	       const Data32* const base,
	       const int stride)
{
  // printf("Computing prob for name '%s'\n",
  // name().c_str());
  logpr rc;
  for (unsigned i=0;i<numComponents;i++) {
    rc += dense1DPMF->p(i)* components[i]->log_p(x,base,stride);
  }
  return rc;
}


//
// compute the log probability of x with stride 'stride'
// This version uses the current global observatio matrix directly,
// and caches the result for EM.
logpr
Mixture::log_p(const unsigned frameIndex, 
	       const unsigned firstFeatureElement)
{
  assert ( basicAllocatedBitIsSet() );

  const float *const x = globalObservationMatrix.floatVecAtFrame(frameIndex,firstFeatureElement);
  const Data32* const base = globalObservationMatrix.baseAtFrame(frameIndex);
  const int stride =  globalObservationMatrix.stride();

  if (cacheComponentsInEmTraining) {

    if (componentCache.size() < (frameIndex+1)) {
      // never have more than 25% more frames than needed while still
      // having log number of allocations in the ultimate size of the observation vectors.
      componentCache.resize( ((frameIndex+1)*5)>>2 );
    }
    // TODO: this stuff is needed only for EM, don't cache components
    // when just doing decoding.
    if (componentCache[frameIndex].cmpProbArray.size() < numComponents)
      componentCache[frameIndex].cmpProbArray.resize(numComponents);

    if (componentCache[frameIndex].prob.valref() != (-LZERO)) {
      // already done for this mixture
      // infoMsg(IM::Mega,"Using cached value for Gaussian component\n");
      // TODO: fix bug here where if firstFeatureElement is not used to check the cache entry
      //       and if this object is shared over multiple firstFeatureElements.
      return componentCache[frameIndex].prob;
    }

    logpr rc;
    // TODO: this stuff is needed only for EM, don't cache components
    // when just doing decoding.
    for (unsigned i=0;i<numComponents;i++) {
      logpr tmp = dense1DPMF->p(i)* components[i]->log_p(x,base,stride);
      // store each component prob value
      componentCache[frameIndex].cmpProbArray[i].prob = tmp;
      rc += tmp;
    }

    // and store the sum as well.
    componentCache[frameIndex].prob = rc;
    return rc;
  } else {
    // don't cache our probabilities.
    logpr rc;
    for (unsigned i=0;i<numComponents;i++) {
      logpr tmp = dense1DPMF->p(i)* components[i]->log_p(x,base,stride);
      rc += tmp;
    }
    return rc;
  }
}



/////////////////
// EM routines //
/////////////////

///////////////////////////////////////////////////////////////////
// NOTES: on algorithm for splitting and vanishing Components with sharing.
// 
// J. Bilmes, $Header$
// 
// If a component component is split, it gets cloned meaning that its
// parameters (e.g., mean and covariances) are cloned.  If another component is
// split, and it shares either the parameters (mean/variance) of another one that
// was split, it will continue to share the parameters (mean/variance) of the
// cloned component.  The cloning relationships are forgotten once all
// is done at the end of the epoch.
//
// First, two sets are used and cleared out at the end of each epoch.
// vanish set, saying a DPMF's particular i'th entry is to be removed.
// split set, saying a DPMF's particular i'th entry is to be split
//
// Next, maps are used for cloning, so that if a GC, mean, variance,
// etc. is to be cloned twice, the clone will be the same the 2nd time
// it is around. This is so that if an object (GC, mean, variance,
// etc.) is shared, its cloned object should be shared as well.  These
// maps are cleared out after the end of each epoch (so that cloning
// will start up again next time), and live in the global GM_Parms
// object.
//       
// The reason for using maps here rather than additional members in
// each object is that splitting/vanishing is probably the exception
// rather than the norm, therefore we don't want to use up more memory
// by adding a field in each object pointing to its current clone when
// that pointer won't be used very often. The map implementation
// uses almost no extra memory unless things start splitting/vanishing.
//
// Here are sketches of the algorithms.
//
// EM end of a Mixture
//   end the DPMF
//   for each entry i in the pmf
//         end comp epoch i
//                // need to end comp epoch in any case because
//                // otherwise cov ref count won't hit zero.
//         if pmf i is below threshold
//                 add (DPMF ptr,i) to vanishing set
//         if pmf i is above threshold
//                 add (DPMF ptr,i) to splitting set
//
// most of the work gets done in the swapping routines.
//
//  mixgauss's swap
//     check on DPMF's new length
//     for each entry i of PMF
//         if pmf&i is in vanishing set
//             skip it, don't swap, move to last position
//         if pmf&i ptr in split set
//             clone cmp i and add it in.
//
//
//  densePMF's swap
//
//     newLen = len;
//     go through pmf, entry i
//        if i is in delete set
//           newLen --;
//        if i is in clone set
//           newLen ++
//     make new len array
//     go through pmf, entry i
//        if i is in delete set
//           skip j
//        if i is in clone set
//           add i an j with split prob
//        else
//           just skip
//     normalize
//
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


void
Mixture::emStartIteration()
{
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;

  if(emOnGoingBitIsSet())
    return; // already done

  if (!emEmAllocatedBitIsSet()) {
    // this is presumably the first time
    weightedPostDistribution.resize(components.size());
    emSetEmAllocatedBit();
  }
  
  // check the length here because
  //  1) it is cheap
  //  2) we want to make sure that this Mixture DPMF length
  //     hasn't changed while this mixture wasn't active (i.e.,
  //     it is possible that during pruning, this mixture
  //     never went through swap below, which is the routine
  //     which adjusts the components.
  if (dense1DPMF->length() != numComponents) {
    error("ERROR: mixture '%s' with '%d' components is trying to start an EM iteration with a dense PMF '%s' of length '%d'\n",
	  name().c_str(),numComponents,dense1DPMF->name().c_str(),
	  dense1DPMF->length());
  }

  // EM iteration is now going.
  emSetOnGoingBit();
  emSetSwappableBit();

  accumulatedProbability = 0.0;
  dense1DPMF->emStartIteration();
  for (unsigned i=0;i<numComponents;i++) {
    components[i]->emStartIteration();
  }
}



/*-
 *-----------------------------------------------------------------------
 * emIncrement()
 *      increment the EM counts
 * 
 * Preconditions:
 *      IMPORTANT: log_p() must have been called on this frameIndex and firstFeatureElement
 *      for the CURRENT observation matrix for this function to be valid. Note that
 *      there is no tagging of the cache to determine if it is valid. This rountine simply assumes
 *      that log_p() has been called, and the cache for the appropriate frame has been filled in.
 *
 * Postconditions:
 *      The EM counts have been added in.
 *
 * Side Effects:
 *      Possible changes to internal structures.
 *
 * Results:
 *      none.
 *
 *-----------------------------------------------------------------------
 */
void
Mixture::emIncrement(logpr prob,
			  const unsigned frameIndex, 
			  const unsigned firstFeatureElement)
{
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;

  if(!emOnGoingBitIsSet())
    emStartIteration();

  if (prob < minIncrementProbabilty) {
    missedIncrementCount++;
    return;
  } 
  accumulatedProbability+= prob;

  const float *const x = globalObservationMatrix.floatVecAtFrame(frameIndex,firstFeatureElement);
  const Data32* const base = globalObservationMatrix.baseAtFrame(frameIndex);
  const int stride = globalObservationMatrix.stride();

  if (cacheComponentsInEmTraining) {
    logpr tmp = prob/componentCache[frameIndex].prob;
    for (unsigned i=0;i<numComponents;i++) {
      weightedPostDistribution[i] =
	componentCache[frameIndex].cmpProbArray[i].prob*tmp;
    }
  } else { // no caching
    // first compute the local mixture posterior distribution.
    logpr sum;
    sum.set_to_zero();
    for (unsigned i=0;i<numComponents;i++) {
      weightedPostDistribution[i] = 
	dense1DPMF->p(i)* components[i]->log_p(x,base,stride);
      sum += weightedPostDistribution[i];
    }
    logpr tmp = prob/sum;
    for (unsigned i=0;i<numComponents;i++) {
      weightedPostDistribution[i] =
	weightedPostDistribution[i]*tmp;
    }
  }

  // increment the mixture weights
  dense1DPMF->emIncrement(prob,weightedPostDistribution);

  // and the components themselves.
  for (unsigned i=0;i<numComponents;i++) {
    components[i]->emIncrement(weightedPostDistribution[i],
			       x,base,stride);
  }
}


/*-
 *-----------------------------------------------------------------------
 * emEndIteration
 *      end the current iteration of the mixture.
 *      First this will end the component mixture.
 *      Next, we compute the mixture specific ratios.
 *      If we find any components that satisfy the tresholds,
 *      we enter them into the map so that other routines
 *      can identify which components are to be added/deleted.
 * 
 * Preconditions:
 *      Basic stuff must be allocated.
 *
 * Postconditions:
 *      iteration has been ended, no longer valid to increment.
 *
 * Side Effects:
 *      changes object values.
 *
 * Results:
 *      nil
 *
 *-----------------------------------------------------------------------
 */
void
Mixture::emEndIteration()
{
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;

  if ( !emOnGoingBitIsSet() )
    return; // done already

  // floor to zero of small.

  accumulatedProbability.floor();
  if (accumulatedProbability.zero()) {
    warning("WARNING: mixture named '%s' did not receive any accumulated probability in EM iteration",name().c_str());
  }

  dense1DPMF->emEndIteration();

  if (dense1DPMF->emAmTrainingBitIsSet()) {
    /////////////////////////////////////////////////////////////////
    // The next bunch of code does splitting/vanishing
    // using both he MCVR/MCSR stuff but also works
    // with the force top split and force bottom vanish
    // framework.
    // The intension is that MCVR will go ahead and vanish
    // everyone according to that value. If there are any additional
    // ones in the bottom N that need to be vanish, that will
    // occur on top of the MCVR.
    // Similarly, MCSR will split anyone that satisfies that ratio.
    // The force top split will then split any additional onces
    // that need to be split according to the top N value.
  

    /////////////////////////////////////////////////////////
    // first split/vanish based on ratio.

    logpr mixCoeffVanishThreshold =
      logpr((double)1.0/numComponents)/logpr(mixCoeffVanishRatio);
    logpr mixCoeffSplitThreshold =
      logpr(mixCoeffSplitRatio)/logpr((double)numComponents);

    assert ( mixCoeffVanishThreshold < mixCoeffSplitThreshold );

    unsigned numVanishedSoFar = 0;
    unsigned numSplitSoFar = 0;
    for (unsigned i=0;i<numComponents;i++) {
      if (numVanishedSoFar < (numComponents-1) && 
	  dense1DPMF->np(i) < mixCoeffVanishThreshold) {
	// make sure not to vanish everyone.
	numVanishedSoFar++;
	MixtureCommon::vanishingComponentSet.insert(pair<Dense1DPMF*,unsigned>(dense1DPMF,i));
      } else if (dense1DPMF->np(i) >= mixCoeffSplitThreshold) {
	numSplitSoFar++;
	MixtureCommon::splittingComponentSet.insert(pair<Dense1DPMF*,unsigned>(dense1DPMF,i));
      }
    }

    ///////////////////////////////////////////////////////////
    // next split/vanish based on force split/vanish

    unsigned localNumTopToForceSplit = numTopToForceSplit;
    unsigned localNumBottomToForceVanish = numBottomToForceVanish;
    if ((localNumTopToForceSplit + localNumBottomToForceVanish) >= 
	dense1DPMF->length()) {
      // need to adjust the top/bottom ones to match
      // the local number of mixtures. Do it "roughly" proportionally.
      double splitFrac = localNumTopToForceSplit/
	(localNumTopToForceSplit + localNumBottomToForceVanish);
      double vanishFrac = 1.0 - splitFrac;
      localNumTopToForceSplit = 
	(unsigned)floor(splitFrac*(dense1DPMF->length()));
      localNumBottomToForceVanish = 
	(unsigned)floor(vanishFrac*(dense1DPMF->length()-1));
    }

    assert ( localNumTopToForceSplit <= dense1DPMF->length() );
    assert ( localNumBottomToForceVanish < dense1DPMF->length() );


    if (localNumTopToForceSplit > 0 
	||
	localNumBottomToForceVanish > 0) {
      // more work

    // if request is such that everything should be split,
    // we make sure that nothing is vanished.
      if (localNumTopToForceSplit == dense1DPMF->length())
	localNumBottomToForceVanish = 0;      

      vector< pair<logpr,unsigned> > coefs;
      coefs.resize(dense1DPMF->length());
      for (unsigned i=0;i<coefs.size();i++) {
	coefs[i].first = dense1DPMF->np(i);
	coefs[i].second = i;
      }

      // need to sort in *accending* order.
      sort(coefs.begin(),coefs.end(),LogpUnsignedPairCompare());

      ///////////////////////////////////////////////////////////////////////
      // Add new split/vanish occurances to appropriate set.
      // Make sure we don't vanish too many, as the vanishing ratio might have
      // already vanished some of the ones that we want to vanish here.
      if (localNumBottomToForceVanish > 0) {
	for (unsigned i=0;
	     (i<localNumBottomToForceVanish)
	       &&
	       (numVanishedSoFar < (numComponents-1))
	       ;i++) {
	  const bool alreadyVanished =
	    (MixtureCommon::vanishingComponentSet.find(pair<Dense1DPMF*,unsigned>(dense1DPMF,coefs[i].second))
	     != MixtureCommon::vanishingComponentSet.end());
	  const bool alreadySplit =
	    (MixtureCommon::splittingComponentSet.find(pair<Dense1DPMF*,unsigned>(dense1DPMF,coefs[i].second))
	     != MixtureCommon::splittingComponentSet.end());

	  //////////////////////////////////////////////////////////
	  // Make sure not to vanish if already requested to split
	  // and don't bother vanishing again if vanishing again.
	  if (!alreadyVanished && !alreadySplit) {

	    numVanishedSoFar++;
	    MixtureCommon::vanishingComponentSet.insert(pair<Dense1DPMF*,unsigned>(dense1DPMF,coefs[i].second));
	  }
	}
      }

      if (localNumTopToForceSplit > 0) {

	for (unsigned j=0,i=(coefs.size()-1);
	     j<localNumTopToForceSplit;
	     j++,i--) {
	  const bool alreadyVanished =
	    (MixtureCommon::vanishingComponentSet.find(pair<Dense1DPMF*,unsigned>(dense1DPMF,coefs[i].second))
	     != MixtureCommon::vanishingComponentSet.end());
	  const bool alreadySplit =
	    (MixtureCommon::splittingComponentSet.find(pair<Dense1DPMF*,unsigned>(dense1DPMF,coefs[i].second))
	     != MixtureCommon::splittingComponentSet.end());

	  //////////////////////////////////////////////////////////
	  // Make sure not to split if already requested to split
	  // and don't bother vanishing again if vanishing again.
	  if (!alreadyVanished && !alreadySplit) {
	    numSplitSoFar++;
	    MixtureCommon::splittingComponentSet.insert(pair<Dense1DPMF*,unsigned>(dense1DPMF,coefs[i].second));
	  }
	}
      }
    }
  }

  // finally end the components iteration.
  for (unsigned i=0;i<numComponents;i++) {
    components[i]->emEndIteration();
  }

  // stop EM
  emClearOnGoingBit();
}



/*-
 *-----------------------------------------------------------------------
 * emSwapCurAndNew
 *      swaps in the new parameters and makes them current.
 *      This routine has extra logic than the standard swap
 *      because of sharing, and the desire to split and remove
 *      mixture components.
 * 
 * Preconditions:
 *      basic allocated, and em iteration finished.
 *
 *
 * Postconditions:
 *      New parameters are in. Components might be
 *      split or removed.
 *
 * Side Effects:
 *      Might drop the link to a component
 *
 * Results:
 *      nil
 *
 *-----------------------------------------------------------------------
 */
void
Mixture::emSwapCurAndNew()
{
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;

  if (!emSwappableBitIsSet())
      return;

  dense1DPMF->emSwapCurAndNew();

  const unsigned newNumComponents = dense1DPMF->length();
  vector < Component* > newComponents;
  newComponents.resize(newNumComponents);

  unsigned newIndex = 0;  
  for (unsigned i=0;i<numComponents;i++) {
    if (MixtureCommon::vanishingComponentSet.
	find(pair<Dense1DPMF*,unsigned>(dense1DPMF,i))
	!= MixtureCommon::vanishingComponentSet.end()) {
      // Do we swap in new values?? No, we keep
      // old values. If someone else needs them, however,
      // the new ones will get swapped in by whoever
      // has not eliminated this component.
      // components[i]->emSwapCurAndNew();
      // Next, do nothing, i.e., don't copy it over to new components.
      ;
    } else if (MixtureCommon::splittingComponentSet.
	       find(pair<Dense1DPMF*,unsigned>(dense1DPMF,i))
	       != MixtureCommon::splittingComponentSet.end()) {
      // first swap in new values
      components[i]->emSwapCurAndNew();
      // next copy it over,
      newComponents[newIndex] = components[i];
      // next copy a cloned copy over,
      newComponents[newIndex+1] = components[i]->noisyClone();
      // two components added
      newIndex += 2;
    } else {
      // first swap in new values
      components[i]->emSwapCurAndNew();
      // next copy it over,
      newComponents[newIndex] = components[i];
      newIndex++;
    }
  }
  // make sure everything is as expected
  assert ( newIndex == newNumComponents );

  // finally, get ready for next iteration.
  components = newComponents;
  numComponents = newNumComponents;
  weightedPostDistribution.resizeIfDifferent(components.size());

  emClearSwappableBit();
}



////////////////////////////////////////////////////////////
// Sample generation
////////////////////////////////////////////////////////////



void Mixture::sampleGenerate(float *const sample,
				  const Data32* const base,
				  const int stride)
{
  assert ( basicAllocatedBitIsSet() );
  error("not implemented");
}









