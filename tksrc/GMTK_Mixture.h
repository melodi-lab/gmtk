/*-
 * GMTK_Mixture
 *        Mixture of components of various types.
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


#ifndef GMTK_MIXTURES_H
#define GMTK_MIXTURES_H

#include "fileParser.h"
#include "logp.h"
#include "machine-dependent.h"
#include "cArray.h"

#include "GMTK_MixtureCommon.h"
#include "GMTK_Component.h"

#include "GMTK_DiagGaussian.h"

#include "GMTK_Dense1DPMF.h"


class Mixture : public MixtureCommon {

  friend class GMTK_Tie;

  // functions in tieSupport.h
  friend MeanVector* find_MeanVector_of_Mixture(Mixture *mixture);
  friend bool all_DiagGaussian(Mixture* mixture);
  friend std::vector<MeanVector*> find_MeanVectors_of_Mixture(Mixture *mixture);

  ///////////////////////////////////////////
  // the (possibly) shared components
  vector < Component* > components;

  ///////////////////////////////////////////
  // For EM, the posteriors
  sArray < logpr > weightedPostDistribution;
  // For EM training,
  // create a 2D component array cache.
  // the first index is by frame number for the current
  // utterance, and the second is by mixture
  // component value (which is stored as we might
  // be doing pruning to save memory).
  struct CompProb {
    logpr prob;
    CompProb() : prob((void*)0) {
      // include constructor to avoid unnecessary initialization.
    }
    // unsigned componentNum; to be used soon.
  };
  struct CompCacheArray {
    cArray < CompProb > cmpProbArray;
    logpr prob;
    unsigned firstFeatureElement;
    CompCacheArray() : prob((void*)0) {
      // initial value is special value indicating that this
      // entry is empty.
      prob.valref() = (-LZERO);
      firstFeatureElement = ~0x0;
    }
    // might want to add more fields later.
  };

  cArray< CompCacheArray > componentCache;

  ///////////////////////////////////////////

  ///////////////////////////////////////////
  // the (possibly) shared 1DPMFs used for the mixture weights.
  Dense1DPMF* dense1DPMF;


  ////////////////////////////////////////////////////////////
  // this dummy variable apparently needs to be here so that gdb 5.0 on
  // Solaris can print out *this. If this is removed, that version of
  // gdb can't do that.
  //   int _dummy1;
  //   int _dummy2;
  //   int _dummy3;
  //   int _dummy4;
  //   Dense1DPMF* dcopy;


public:

  Mixture(const int dim,ContinuousImplementation mtype=ci_mixture)
    : MixtureCommon(dim,mtype)
  { }
  //  Mixture() {}
  ~Mixture() {}

  //////////////////////////////////////////////
  // read/write basic parameters
  void read(iDataStreamFile& is);
  void write(oDataStreamFile& os);

  unsigned totalNumberParameters();

  // component cache support
  // make all component cache entries effectively empty.
  void emptyComponentCache();
  // clear the memory associated with the component cache.
  void freeComponentCache();

  // these routines are used to not save 
  // components (and their means, variances, parms, etc.) 
  // that are not actively used in a parameter file (such
  // as those that have vanished away).
  void recursivelyClearUsedBit() { 
    emClearUsedBit();
    for (unsigned i=0;i<components.size();i++)
      components[i]->recursivelyClearUsedBit();
  }
  void recursivelySetUsedBit() {
    emSetUsedBit();    
    for (unsigned i=0;i<components.size();i++)
      components[i]->recursivelySetUsedBit();
  }

  //////////////////////////////////
  // set all current parameters to valid but random values
  void makeRandom();
  // set all current parameters to valid but "uniform" values 
  // (for Gaussians this means N(0,1))
  void makeUniform();
  //////////////////////////////////

  //////////////////////////////////
  // probability evaluation
  logpr log_p(const float *const x,    // real-valued scoring obs at time t
	      const Data32* const base, // ptr to base obs at time t
	      const int stride);       // stride
  // a version that uses the current global obervation matrix directly.
  logpr log_p(const unsigned frameIndex,
	      const unsigned firstFeatureElement);
  logpr maxValue();
  //////////////////////////////////


  //////////////////////////////////
  // Full Baum-Welch EM training  //
  //////////////////////////////////
  void emStartIteration();
  void emIncrement(logpr prob, 
		   const unsigned frameIndex, 
		   const unsigned firstFeatureElement);
  void emEndIteration();
  void emSwapCurAndNew();


  // parallel training
  void emStoreObjectsAccumulators(oDataStreamFile& ofile,
				  bool writeLogVals = true,
				  bool writeZeros = false) {}
  void emLoadObjectsDummyAccumulators(iDataStreamFile& ifile) {}
  void emZeroOutObjectsAccumulators() {}
  void emLoadObjectsAccumulators(iDataStreamFile& ifile) {}
  void emAccumulateObjectsAccumulators(iDataStreamFile& ifile) {}
  const string typeName() { return "Mixture"; }
  //////////////////////////////////

  //////////////////////////////////
  // Sample Generation            //
  //////////////////////////////////
  void sampleGenerate(float *sample,
		      const Data32* const base,
		      const int stride);
  //////////////////////////////////


};


#endif


