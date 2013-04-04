/*-
 * GMTK_MissingFeatureScaledDiagGaussian
 *        Missing Feaure Scaled Diagonal Covariance Gaussian Distribution.
 *
 *  Written by Jeff Bilmes <bilmes@ee.washington.edu>
 * 
 *  $Header$
 * 
 * Copyright (C) 2010 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 *
 *
 */


#ifndef GMTK_MISSINGFEATURESCALEDDIAGGAUSSIAN_H
#define GMTK_MISSINGFEATURESCALEDDIAGGAUSSIAN_H

#include "fileParser.h"
#include "logp.h"
#include "machine-dependent.h"
#include "sArray.h"




#include "GMTK_GaussianComponent.h"
#include "GMTK_MeanVector.h"
#include "GMTK_DiagCovarVector.h"
#include "GMTK_Component.h"
#include "GMTK_RealMatrix.h"

class MissingFeatureScaledDiagGaussian : public GaussianComponent {

  friend class GMTK_Tie;
  friend MeanVector* find_MeanVector_of_MissingFeatureScaledDiagGaussian(MissingFeatureScaledDiagGaussian *diag_gaussian);

  ///////////////////////////////////////////////////////
  // The means. 
  // This might be tied with multiple other distributions.
  MeanVector* mean;

  // For EM Training: Local copy of mean & diagCov accumulators for this MissingFeatureScaledDiagGaussian,
  // which is needed for sharing.
  sArray<float> nextMeans;
  sArray<float> nextDiagCovars;

  ///////////////////////////////////////////////////////
  // The diagonal of the covariance matrix
  // This might be tied with multiple other distributions.
  DiagCovarVector* covar;

  // The exponential scale vector
  RealMatrix* scale;

  ////////////////////////////////////////////////
  // the local accumulated probability needed for
  // each element of the vector.
  sArray<logpr> elementAccumulatedProbability;

  /////////////////////////////////////////////////
  // modify the usage counts of any members that use them; typically
  // called with amount=1 or -1
  void adjustNumTimesShared(int amount){
    mean->numTimesShared += amount;
    covar->numTimesShared += amount;
  };

public:

  
  MissingFeatureScaledDiagGaussian(const int dim) : GaussianComponent(dim) { }
  ~MissingFeatureScaledDiagGaussian() {}

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

  // return the maximum possible value of this component. For a Gaussian,
  // just return the value if the argument to the exponent is zero.
  virtual logpr maxValue() { return logpr(0,covar->log_inv_normConst()); }

  //////////////////////////////////


  //////////////////////////////////
  // Full Baum-Welch EM training  //
  //////////////////////////////////
  void emStartIteration();
  void emIncrement(logpr prob,
		   const float*f,
		   const Data32* const base,
		   const int stride);
  void emIncrementMeanDiagCovar(logpr prob,
				const float prob,
				const float *const f,
				const unsigned len,
				float *meanAccumulator,
				float *diagCovarAccumulator);

  void emEndIteration();
  void emSwapCurAndNew();

  void fkIncrementMeanDiagCovar(logpr prob,
				const float fprob,
				const float *const f,
				const unsigned len,
				float *curMeans,
				float *curDiagCovars,
				float *meanAccumulator,
				float *diagCovarAccumulator);
  

  // parallel training
  void emStoreObjectsAccumulators(oDataStreamFile& ofile,
				  bool writeLogVals = true,
				  bool writeZeros = false);
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
