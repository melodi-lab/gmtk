/*-
 * GMTK_ZeroScoreMixGaussian
 *        Mixture of gaussian that always returns probability 0 regardless
 *        of its parents or the corresponding feature values. This object
 *        can be useful in certain situations. The name of this mixture
 *        is internal, and is "internal::ZeroScoreGaussianMixture"
 *        Note that this object has dim = 0. It is up to the other checking
 *        routines when checking the dimension to make sure that zero matches
 *        anything else, as zero is not normally a valid specifyable dimension.
 *
 *  Written by Jeff Bilmes <bilmes@ee.washington.edu>
 * 
 *  $Header$
 * 
 * Copyright (c) 2001, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about
 * the suitability of this software for any purpose.  It is provided
 * "as is" without express or implied warranty.
 *
 */


#ifndef GMTK_ZEROSCOREMIXGAUSSIANS_H
#define GMTK_ZEROSCOREMIXGAUSSIANS_H

#include "fileParser.h"
#include "logp.h"
#include "machine-dependent.h"


#include "GMTK_MixGaussians.h"

// name to use in collections when refering to one of these objects.
#define ZEROSCOREMIXGAUSSIAN_NAME "internal:ZeroScore"

class ZeroScoreMixGaussian : public MixGaussians {
 
public:

  ZeroScoreMixGaussian() 
    : MixGaussians(0,ci_zeroScoreMixGaussian)
  { _name = ZEROSCOREMIXGAUSSIAN_NAME; setBasicAllocatedBit(); }
  ~ZeroScoreMixGaussian() {}

  //////////////////////////////////////////////
  // read/write basic parameters
  void read(iDataStreamFile& is) { }
  void write(oDataStreamFile& os) { }

  unsigned totalNumberParameters() { return 0; }

  // these routines are used to not save gaussian
  // components (and their means, variances, etc.) 
  // that are not actively used in a parameter file (such
  // as those that have vanished away).
  void recursivelyClearUsedBit() { }
  void recursivelySetUsedBit() { }

  //////////////////////////////////
  // set all current parameters to valid but random values
  void makeRandom() {}
  // set all current parameters to valid but "uniform" values 
  // (for Gaussians this means N(0,1))
  void makeUniform() {}
  //////////////////////////////////

  //////////////////////////////////
  // probability evaluation
  logpr log_p(const float *const x,    // real-valued scoring obs at time t
	      const Data32* const base, // ptr to base obs at time t
	      const int stride)       // stride
    { // a no-argument logpr returns a zero value.
      return logpr(); }
  // a version that uses the current global obervation matrix directly.
  logpr log_p(const unsigned frameIndex,
	      const unsigned firstFeatureElement)
    { // a no-argument logpr returns a zero value.
      return logpr(); }
  //////////////////////////////////


  //////////////////////////////////
  // Full Baum-Welch EM training  //
  //////////////////////////////////
  void emStartIteration() {}
  void emIncrement(logpr prob, 
		   const unsigned frameIndex, 
		     const unsigned firstFeatureElement) {}
  void emEndIteration() {}
  void emSwapCurAndNew() {}


  // parallel training
  void emStoreAccumulators(oDataStreamFile& ofile) {}
  void emLoadAccumulators(iDataStreamFile& ifile) {}
  void emAccumulateAccumulators(iDataStreamFile& ifile) {}

  void emStoreObjectsAccumulators(oDataStreamFile& ofile) {}
  void emLoadObjectsDummyAccumulators(iDataStreamFile& ifile) {}
  void emZeroOutObjectsAccumulators() {}
  void emLoadObjectsAccumulators(iDataStreamFile& ifile) {}
  void emAccumulateObjectsAccumulators(iDataStreamFile& ifile) {}
  const string typeName() { return "Zero Score Gaussian mixture"; }
  //////////////////////////////////

  //////////////////////////////////
  // Sample Generation            //
  //////////////////////////////////
  void sampleGenerate(float *sample,
		      const Data32* const base,
		      const int stride) {}
  //////////////////////////////////

};


#endif


