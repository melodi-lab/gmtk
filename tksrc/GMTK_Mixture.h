/*-
 * GMTK_MixGaussians
 *        Mixture of gaussians of various types.
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


#ifndef GMTK_MIXGAUSSIANS_H
#define GMTK_MIXGAUSSIANS_H

#include "fileParser.h"
#include "logp.h"
#include "machine-dependent.h"


#include "GMTK_MixGaussiansCommon.h"
#include "GMTK_GaussianComponent.h"

#include "GMTK_DiagGaussian.h"

#include "GMTK_Dense1DPMF.h"


class MixGaussians : public MixGaussiansCommon {

  ///////////////////////////////////////////
  // the (possibly) shared components
  vector < GaussianComponent* > components;

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
    // unsigned componentNum; to be used soon.
  };
  struct CompCacheArray {
    vector < CompProb > cmpProbArray;
    logpr prob;
    // might want to add more fields later.
  };

  vector< CompCacheArray > componentCache;

  ///////////////////////////////////////////

  ///////////////////////////////////////////
  // the (possibly) shared 1DPMFs used for the mixture weights.
  Dense1DPMF* dense1DPMF;


  ////////////////////////////////////////////////////////////
  // this dummy variable apparently needs to be here so that gdb 5.0 on
  // Solaris can print out *this. If this is removed, that version of
  // gdb can't do that.
  int _dummy;
 
public:

  MixGaussians(const int dim) 
    : MixGaussiansCommon(dim,ci_mixGaussian)
  { }
  ~MixGaussians() {}

  //////////////////////////////////////////////
  // read/write basic parameters
  void read(iDataStreamFile& is);
  void write(oDataStreamFile& os);

  unsigned totalNumberParameters();

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
  void emStoreAccumulators(oDataStreamFile& ofile);
  void emStoreZeroAccumulators(oDataStreamFile& ofile);
  void emLoadAccumulators(iDataStreamFile& ifile);
  void emAccumulateAccumulators(iDataStreamFile& ifile);
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


