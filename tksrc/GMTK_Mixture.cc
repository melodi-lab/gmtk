/*-
 * GMTK_MixGaussians.cc
 *        Code for mixtures of Gaussians of a variety of types.
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

#include "GMTK_MixGaussians.h"
#include "GMTK_GMParms.h"
#include "GMTK_ObservationMatrix.h"



void
MixGaussians::read(iDataStreamFile& is)
{
  // read name of self
  NamedObject::read(is);
  
  // read number of mixtures
  is.read(numComponents,"MixGaussians::read numComponents");
  if (numComponents <= 0) 
    error("Error: MixGaussians '%s' has number of components = %d\n",
	  _name.c_str(),numComponents);

  // read name of dense 1d PMF to use for mixtures
  string str;
  is.read(str);
  if (GM_Parms.dPmfsMap.find(str) == GM_Parms.dPmfsMap.end()) {
    error("ERROR: Mix Gaussian '%s' in file '%s', can't find PMF named '%s'\n",
	  is.fileName(),
	  _name.c_str(),str.c_str());
  }

  dense1DPMF = GM_Parms.dPmfs[
			      GM_Parms.dPmfsMap[str]
  ];

  // now make sure that this one matches the number of components.
  if (numComponents != dense1DPMF->length()) {
    error("ERROR: MixGaussians '%s' in file '%s', PMF named '%s' has %d elements but we need %d\n",
	  _name.c_str(),is.fileName(),
	  str.c_str(),dense1DPMF->length(),numComponents);
  }

  // now read 'name' pointers to all the Gaussians.
  components.resize(numComponents);
  for (unsigned i=0;i<numComponents;i++) {
    is.read(str);
    if (GM_Parms.gaussianComponentsMap.find(str) == GM_Parms.gaussianComponentsMap.end()) {
      error("ERROR: MixGaussians '%s' in file '%s', can't find Gaussian Component named '%s'\n",_name.c_str(),is.fileName(),str.c_str());
    }
    GaussianComponent*gc = GM_Parms.gaussianComponents [
	GM_Parms.gaussianComponentsMap[str]
    ];
    components[i] = gc;
    if (gc->dim() != dim()) {
      error("ERROR: MixGaussians '%s' in file '%s' of dim %d trying to use component '%s' of dim %d\n",
	    name().c_str(),is.fileName(),dim(),gc->name().c_str(),gc->dim());
    }
  }

  // make ready for probability evaluation.
  componentCache.resize(10);
  setBasicAllocatedBit();
}

void
MixGaussians::write(oDataStreamFile& os)
{
  assert ( basicAllocatedBitIsSet() );

  NamedObject::write(os);
  os.nl();
  // read number of mixture components
  os.write(numComponents,"MixGaussians::write numComponents");

  // write name of dense 1d PMF to use for mixtures
  os.write(dense1DPMF->name()); os.nl();

  // and the component names
  for (unsigned i=0;i<numComponents;i++) {
    os.write(components[i]->name());
  }
  os.nl();
}



void
MixGaussians::makeUniform()
{
  dense1DPMF->makeUniform();
  for (unsigned i=0;i<numComponents;i++) {
    components[i]->makeUniform();
  }
}


void
MixGaussians::makeRandom()
{
  dense1DPMF->makeRandom();
  for (unsigned i=0;i<numComponents;i++) {
    components[i]->makeRandom();
  }
}


//
// compute the log probability of x with stride 'stride'
// 
logpr
MixGaussians::log_p(const float *const x,
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
MixGaussians::log_p(const unsigned frameIndex, 
		    const unsigned firstFeatureElement)

{
  assert ( basicAllocatedBitIsSet() );

  const float *const x = globalObservationMatrix.floatVecAtFrame(frameIndex,firstFeatureElement);
  const Data32* const base = globalObservationMatrix.baseAtFrame(frameIndex);
  const int stride =  globalObservationMatrix.stride;

  if (componentCache.size() < (frameIndex+1)) {
    // never have more than 25% more frames than needed while still
    // having log number of allocations in the ultimate size of the observation vectors.
    componentCache.resize( ((frameIndex+1)*5)>>2 );
  }
  if (componentCache[frameIndex].cmpProbArray.size() < numComponents)
    componentCache[frameIndex].cmpProbArray.resize(numComponents);

  logpr rc;
  for (unsigned i=0;i<numComponents;i++) {
    logpr tmp = dense1DPMF->p(i)* components[i]->log_p(x,base,stride);
    // store each component prob value
    componentCache[frameIndex].cmpProbArray[i].prob = tmp;
    rc += tmp;
  }
  // and store the sum as well.
  componentCache[frameIndex].prob = rc;
  return rc;
}



/////////////////
// EM routines //
/////////////////

///////////////////////////////////////////////////////////////////
// NOTES: on algorithm for splitting and vanishing Gaussians with sharing.
// 
// J. Bilmes, $Header$
// 
// If a Gaussian component is split, it gets cloned meaning that its
// mean and covariances are cloned.  If another Gaussian component is
// split, and it shares either the mean/variance of another one that
// was split, it will continue to share the mean/variance of the
// cloned Gaussian.  The cloning relationships are forgotten once all
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
// EM end of a MixGaussian
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
MixGaussians::emStartIteration()
{
  assert ( basicAllocatedBitIsSet() );

  if (!GM_Parms.amTrainingMixGaussians())
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
  //  2) we want to make sure that this Gaussians DPMF length
  //     hasn't changed while this mixture wasn't active (i.e.,
  //     it is possible that during pruning, this mixture
  //     never went through swap below, which is the routine
  //     which adjusts the components.
  if (dense1DPMF->length() != numComponents) {
    error("ERROR: Gaussian mixture '%s' with '%d' components is trying to start an EM iteration with a dense PMF '%s' of length '%d'\n",
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
MixGaussians::emIncrement(logpr prob,
			  const unsigned frameIndex, 
			  const unsigned firstFeatureElement)
{
  assert ( basicAllocatedBitIsSet() );

  if (!GM_Parms.amTrainingMixGaussians())
    return;  

  emStartIteration();

  if (prob < minIncrementProbabilty) {
    missedIncrementCount++;
    return;
  } 
  accumulatedProbability+= prob;

  const float *const x = globalObservationMatrix.floatVecAtFrame(frameIndex,firstFeatureElement);
  const Data32* const base = globalObservationMatrix.baseAtFrame(frameIndex);
  const int stride = globalObservationMatrix.stride;

  logpr tmp = prob/componentCache[frameIndex].prob;
  for (unsigned i=0;i<numComponents;i++) {
    weightedPostDistribution[i] =
      componentCache[frameIndex].cmpProbArray[i].prob*tmp;
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
 *      end the current iteration of the gaussian mixtures.
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
MixGaussians::emEndIteration()
{
  assert ( basicAllocatedBitIsSet() );

  if (!GM_Parms.amTrainingMixGaussians())
    return;

  if ( !emOnGoingBitIsSet() )
    return; // done already

  // floor to zero of small.

  accumulatedProbability.floor();
  if (accumulatedProbability.zero()) {
    warning("WARNING: Gaussian mixture named '%s' did not receive any accumulated probability in EM iteration",name().c_str());
  }

  dense1DPMF->emEndIteration();

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
      MixGaussiansCommon::vanishingComponentSet.insert(pair<Dense1DPMF*,unsigned>(dense1DPMF,i));
    } else if (dense1DPMF->np(i) >= mixCoeffSplitThreshold) {
      numSplitSoFar++;
      MixGaussiansCommon::splittingComponentSet.insert(pair<Dense1DPMF*,unsigned>(dense1DPMF,i));
    }
    // need to end component iteration in both cases.
    components[i]->emEndIteration();
  }

  ///////////////////////////////////////////////////////////
  // next split/vanish based on force split/vanish

  unsigned localNumTopToForceSplit = numTopToForceSplit;
  unsigned localNumBottomToForceVanish = numBottomToForceVanish;
  if ((localNumTopToForceSplit + localNumBottomToForceVanish) >= 
      dense1DPMF->length()) {
    // need to adjust the top/bottom ones to match
    // the local number of mixtures. Do it proportionally.
    double splitFrac = localNumTopToForceSplit/
      (localNumTopToForceSplit + localNumBottomToForceVanish);
    double vanishFrac = 1.0 - splitFrac;
    localNumTopToForceSplit = 
      (unsigned)floor(splitFrac*(dense1DPMF->length()-1));
    localNumBottomToForceVanish = 
      (unsigned)floor(vanishFrac*(dense1DPMF->length()-1));
  }

  assert ( localNumTopToForceSplit < dense1DPMF->length() );
  assert ( localNumBottomToForceVanish < dense1DPMF->length() );

  if (localNumTopToForceSplit > 0 
      ||
      localNumBottomToForceVanish > 0) {
    // more work

    vector< pair<logpr,unsigned> > coefs;
    coefs.resize(dense1DPMF->length());
    for (unsigned i=0;i<coefs.size();i++) {
      coefs[i].first = dense1DPMF->np(i);
      coefs[i].second = i;
    }

    // need to sort in *accending* order.
    LogpUnsignedPairCompare lupc;
    sort(coefs.begin(),coefs.end(),lupc);

    /////////////////////////////////////
    // Add new split/vanish occurances to appropriate set.
    // Make sure we don't vanish too many, as the vanishing ratio might have
    // already vanished some of the ones that we want to vanish here.
    if (localNumBottomToForceVanish > 0) {
      for (unsigned i=0;
	   (i<localNumBottomToForceVanish)
	     &&
	     (numVanishedSoFar < (numComponents-1))
	     ;i++) {
	bool alreadyVanished =
	  (MixGaussiansCommon::vanishingComponentSet.find(pair<Dense1DPMF*,unsigned>(dense1DPMF,coefs[i].second))
	   != MixGaussiansCommon::vanishingComponentSet.end());
	bool alreadySplit =
	  (MixGaussiansCommon::splittingComponentSet.find(pair<Dense1DPMF*,unsigned>(dense1DPMF,coefs[i].second))
	   != MixGaussiansCommon::splittingComponentSet.end());
	if (!alreadyVanished && !alreadySplit) {
	  // then ok to vanish
	  numVanishedSoFar++;
	  MixGaussiansCommon::vanishingComponentSet.insert(pair<Dense1DPMF*,unsigned>(dense1DPMF,coefs[i].second));
	}
      }
    }

    if (localNumTopToForceSplit > 0) {

      for (unsigned j=0,i = (coefs.size()-1);
	   ;j<localNumTopToForceSplit;j++,i--) {
	bool alreadyVanished =
	  (MixGaussiansCommon::vanishingComponentSet.find(pair<Dense1DPMF*,unsigned>(dense1DPMF,coefs[i].second))
	   != MixGaussiansCommon::vanishingComponentSet.end());
	bool alreadySplit =
	  (MixGaussiansCommon::splittingComponentSet.find(pair<Dense1DPMF*,unsigned>(dense1DPMF,coefs[i].second))
	   != MixGaussiansCommon::splittingComponentSet.end());

	if (!alreadyVanished && !alreadySplit) {
	  numSplitSoFar++;
	  MixGaussiansCommon::splittingComponentSet.insert(pair<Dense1DPMF*,unsigned>(dense1DPMF,coefs[i].second));
	}
      }
    }
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
 *      Might drop the link to a Gaussian component
 *
 * Results:
 *      nil
 *
 *-----------------------------------------------------------------------
 */
void
MixGaussians::emSwapCurAndNew()
{
  assert ( basicAllocatedBitIsSet() );

  if (!GM_Parms.amTrainingMixGaussians())
    return;

  if (!emSwappableBitIsSet())
      return;

  dense1DPMF->emSwapCurAndNew();

  const unsigned newNumComponents = dense1DPMF->length();
  vector < GaussianComponent* > newComponents;
  newComponents.resize(newNumComponents);

  unsigned newIndex = 0;  
  for (unsigned i=0;i<numComponents;i++) {
    if (MixGaussiansCommon::vanishingComponentSet.
	find(pair<Dense1DPMF*,unsigned>(dense1DPMF,i))
	!= MixGaussiansCommon::vanishingComponentSet.end()) {
      // Do we swap in new values?? No, we keep
      // old values. If someone else needs them, however,
      // the new ones will get swapped in by whoever
      // has not eliminated this component.
      // components[i]->emSwapCurAndNew();
      // Next, do nothing, i.e., don't copy it over to new components.
      ;
    } else if (MixGaussiansCommon::splittingComponentSet.
	       find(pair<Dense1DPMF*,unsigned>(dense1DPMF,i))
	       != MixGaussiansCommon::splittingComponentSet.end()) {
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



void
MixGaussians::emStoreAccumulators(oDataStreamFile& ofile)
{
  assert ( basicAllocatedBitIsSet() );
  if ( !emEmAllocatedBitIsSet() ) {
    warning("WARNING: storing zero accumulators for mix gaussian '%s'\n",
	    name().c_str());
    emStoreZeroAccumulators(ofile);
    return;
  }
  EMable::emStoreAccumulators(ofile);
}

void
MixGaussians::emStoreZeroAccumulators(oDataStreamFile& ofile)
{
  assert ( basicAllocatedBitIsSet() );
  EMable::emStoreZeroAccumulators(ofile);
}

void
MixGaussians::emLoadAccumulators(iDataStreamFile& ifile)
{
  assert ( basicAllocatedBitIsSet() );
  assert ( emEmAllocatedBitIsSet() );
  EMable::emLoadAccumulators(ifile);
}


void
MixGaussians::emAccumulateAccumulators(iDataStreamFile& ifile)
{
  assert ( basicAllocatedBitIsSet() );
  assert ( emEmAllocatedBitIsSet() );
  EMable::emAccumulateAccumulators(ifile);
}


////////////////////////////////////////////////////////////
// Sample generation
////////////////////////////////////////////////////////////



void MixGaussians::sampleGenerate(float *const sample,
				  const Data32* const base,
				  const int stride)
{
  assert ( basicAllocatedBitIsSet() );
  error("not implemented");
}









