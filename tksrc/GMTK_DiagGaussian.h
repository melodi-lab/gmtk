/*-
 * GMTK_DiagGaussian
 *        Diagonal Covariance Gaussian Distribution, no
 *        additional dependencies to other observation.
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


#ifndef GMTK_DIAGGAUSSIAN_H
#define GMTK_DIAGGAUSSIAN_H

#include "fileParser.h"
#include "logp.h"
#include "machine-dependent.h"
#include "sArray.h"

#include "GMTK_GaussianComponent.h"
#include "GMTK_MeanVector.h"
#include "GMTK_DiagCovarVector.h"

class DiagGaussian : public GaussianComponent {

  friend class GMTK_Tie;
  friend MeanVector* find_MeanVector_of_DiagGaussian(DiagGaussian *diag_gaussian);

  ///////////////////////////////////////////////////////
  // The means. 
  // This might be tied with multiple other distributions.
  MeanVector* mean;

  // For EM Training: Local copy of mean & diagCov accumulators for this DiagGaussian,
  // which is needed for sharing.
  sArray<float> nextMeans;
  sArray<float> nextDiagCovars;

  ///////////////////////////////////////////////////////
  // The diagonal of the covariance matrix
  // This might be tied with multiple other distributions.
  DiagCovarVector* covar;

 
public:

  
  DiagGaussian(const int dim) : GaussianComponent(dim) { }
  ~DiagGaussian() {}

  //////////////////////////////////////////////
  // read/write basic parameters
  void read(iDataStreamFile& is);
  void write(oDataStreamFile& os);

  // create a copy of self, but with slightly perturbed
  // means/variance values.
  Component* noisyClone();

  /////////////////////////////////////////////////
  // create a copy of self, with entirely new parameters with
  // identical values; NOTHING is shared
  Component* identicalIndependentClone();

  //////////////////////////////////
  // set all current parameters to valid but random values
  void makeRandom();
  // set all current parameters to valid but "uniform" values 
  // (for Gaussians this means N(0,1))
  void makeUniform();
  ///////////////////////////////////

  ///////////////////////////////////
  unsigned totalNumberParameters() { return 
				       mean->totalNumberParameters() +
				       covar->totalNumberParameters(); }


  // these routines are used to not save gaussian
  // components (and their means, variances, etc.) 
  // that are not actively used in a parameter file (such
  // as those that have vanished away).
  void recursivelyClearUsedBit() { 
    emClearUsedBit();
    mean->recursivelyClearUsedBit();
    covar->recursivelyClearUsedBit();
  }
  void recursivelySetUsedBit() {
    emSetUsedBit();    
    mean->recursivelySetUsedBit();
    covar->recursivelySetUsedBit();
  }




  //////////////////////////////////
  // probability evaluation
  logpr log_p(const float *const x,    // real-valued scoring obs at time t
	      const Data32* const base, // ptr to base obs at time t
	      const int stride);       // stride
  //////////////////////////////////


  //////////////////////////////////
  // Full Baum-Welch EM training  //
  //////////////////////////////////
  void emStartIteration();
  void emIncrement(logpr prob,
		   const float*f,
		   const Data32* const base,
		   const int stride);
  static void emIncrementMeanDiagCovar(const float prob,
				       const float *const f,
				       const unsigned len,
				       float *meanAccumulator,
				       float *diagCovarAccumulator);

  void emEndIteration();
  void emSwapCurAndNew();



  // parallel training
  void emStoreObjectsAccumulators(oDataStreamFile& ofile);
  void emLoadObjectsDummyAccumulators(iDataStreamFile& ifile);
  void emZeroOutObjectsAccumulators();
  void emLoadObjectsAccumulators(iDataStreamFile& ifile);
  void emAccumulateObjectsAccumulators(iDataStreamFile& ifile);
  const string typeName() { return "Diag Gaussian"; }
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
