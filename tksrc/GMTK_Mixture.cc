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
#include <ieeefp.h>
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
}

void
MixGaussians::write(oDataStreamFile& os)
{
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


void
MixGaussians::emStartIteration()
{
  if (!GM_Parms.amTrainingMixGaussians())
    return;
  if(emOnGoingBitIsSet())
    return; // already done

  if (!emEmAllocatedBitIsSet()) {
    // this is presumably the first time
    weightedPostDistribution.resize(components.size());
    emSetEmAllocatedBit();
  }

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

#if 0
void
MixGaussians::emIncrement(logpr prob,
			  const float *f,
			  const Data32* const base,
			  const int stride)
{

  // printf("In emIncrement, calling with log prob = %g, float = 0x%X\n",
  // prob.val(),(void*)f);

  if (!GM_Parms.amTrainingMixGaussians())
    return;  

  emStartIteration();

  if (prob < minIncrementProbabilty) {
    missedIncrementCount++;
    return;
  } 
  accumulatedProbability+= prob;

  // first compute the local Gaussian mixture posterior distribution.
  logpr sum;
  for (unsigned i=0;i<numComponents;i++) {
    postDistribution[i] = 
      dense1DPMF->p(i)* components[i]->log_p(f,base,stride);
    sum += postDistribution[i];
  }

  logpr tmp = prob/sum;
  for (unsigned i=0;i<numComponents;i++) {
    postDistribution[i] =
      postDistribution[i]*tmp;
  }

  // increment the mixture weights
  dense1DPMF->emIncrement(prob,postDistribution);

  // and the components themselves.
  for (unsigned i=0;i<numComponents;i++) {
    components[i]->emIncrement(postDistribution[i],
			       f,base,stride);
  }

}
#endif


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

  const logpr mixCoeffVanishThreshold =   
    logpr((double)1.0/numComponents)/logpr(mixCoeffVanishRatio);

  // and the components themselves.
  for (unsigned i=0;i<numComponents;i++) {
    components[i]->emIncrement(weightedPostDistribution[i],
			       x,base,stride);
  }
}


void
MixGaussians::emEndIteration()
{

  if (!GM_Parms.amTrainingMixGaussians())
    return;

  if ( !emOnGoingBitIsSet() )
    return; // done already

  if (accumulatedProbability.zero()) {
    warning("WARNING: Gaussian mixture named '%s' did not receive any accumulated probability in EM iteration",name().c_str());
  }

  dense1DPMF->emEndIteration();

  for (unsigned i=0;i<numComponents;i++)
    components[i]->emEndIteration();

  // stop EM
  emClearOnGoingBit();
}



void
MixGaussians::emSwapCurAndNew()
{

  if (!GM_Parms.amTrainingMixGaussians())
    return;

  if (!emSwappableBitIsSet())
      return;

  dense1DPMF->emSwapCurAndNew();
  for (unsigned i=0;i<numComponents;i++) {
    components[i]->emSwapCurAndNew();
  }
  emClearSwappableBit();
}




void
MixGaussians::emStoreAccumulators(oDataStreamFile& ofile)
{
  error("not implemented");
}

void
MixGaussians::emLoadAccumulators(iDataStreamFile& ifile)
{
  error("not implemented");
}


void
MixGaussians::emAccumulateAccumulators(iDataStreamFile& ifile)
{
  error("not implemented");
}


////////////////////////////////////////////////////////////
// Sample generation
////////////////////////////////////////////////////////////



void MixGaussians::sampleGenerate(float *const sample,
				  const Data32* const base,
				  const int stride)
{
  error("not implemented");
}









