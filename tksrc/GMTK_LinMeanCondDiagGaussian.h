/*-
 * GMTK_LinMeanCondDiagGaussian
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


#ifndef GMTK_LINMEANCONDDIAGGAUSSIAN_H
#define GMTK_LINMEANCONDDIAGGAUSSIAN_H

#include "fileParser.h"
#include "logp.h"
#include "machine-dependent.h"
#include "sArray.h"

#include "GMTK_RandomVariable.h"
#include "GMTK_GaussianComponent.h"
#include "GMTK_MeanVector.h"
#include "GMTK_DiagCovarVector.h"
#include "GMTK_DlinkMatrix.h"
#include "GMTK_Dlinks.h"

class LinMeanCondDiagGaussian : public GaussianComponent {

  ///////////////////////////////////////////////////////
  // The means. 
  // This might be tied with multiple other distributions.
  MeanVector* mean;

  ///////////////////////////////////////////////////////
  // The diagonal of the covariance matrix
  // This might be tied with multiple other distributions.
  DiagCovarVector* covar;
  
  ///////////////////////////////////////////////////////
  // parameters for the dlink structure
  DlinkMatrix* dLinkMat;

  // Accumulators for EM training: First, a local mean & diagCov
  // accumulator, needed for sharing.
  // WARNING: If changed from float -> double, accumulator
  // writing routines will also need to change.
  // EX 
  sArray<float> xAccumulators;
  // E[x^2] accumulators
  sArray<float> xxAccumulators;
  // Next, an array containing the accumulators for
  // the feature vector 'x' times the conditioning variables 'z'
  // i.e., this is E[XZ'].
  sArray<float> xzAccumulators;
  // next, E[Z]
  sArray<float> zAccumulators;
  //  E[ZZ']
  sArray<float> zzAccumulators;

public:

  LinMeanCondDiagGaussian(const int dim) : GaussianComponent(dim) { }
  ~LinMeanCondDiagGaussian() {}

  //////////////////////////////////////////////
  // read/write basic parameters
  void read(iDataStreamFile& is);
  void write(oDataStreamFile& os);

  // create a copy of self, but with slightly perturbed
  // means/variance values.
  GaussianComponent* noisyClone();


  //////////////////////////////////
  // set all current parameters to valid but random values
  void makeRandom();
  // set all current parameters to valid but "uniform" values 
  // (for Gaussians this means N(0,1))
  void makeUniform();
  //////////////////////////////////


  ///////////////////////////////////
  unsigned totalNumberParameters() { return 
				       mean->totalNumberParameters()+
				       covar->totalNumberParameters()+
				       dLinkMat->totalNumberParameters(); }




  // these routines are used to not save gaussian
  // components (and their means, variances, etc.) 
  // that are not actively used in a parameter file (such
  // as those that have vanished away).
  void recursivelyClearUsedBit() { 
    emClearUsedBit();
    mean->recursivelyClearUsedBit();
    covar->recursivelyClearUsedBit();
    dLinkMat->recursivelyClearUsedBit();
  }
  void recursivelySetUsedBit() {
    emSetUsedBit();    
    mean->recursivelySetUsedBit();
    covar->recursivelySetUsedBit();
    dLinkMat->recursivelySetUsedBit();
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
  void emEndIteration();
  void emEndIterationNoSharing();
  void emEndIterationSharedCovars();
  void emEndIterationSharedAll();
  void emSwapCurAndNew();


  // parallel training
  void emStoreObjectsAccumulators(oDataStreamFile& ofile);
  void emLoadObjectsDummyAccumulators(iDataStreamFile& ifile);
  void emZeroOutObjectsAccumulators();
  void emLoadObjectsAccumulators(iDataStreamFile& ifile);
  void emAccumulateObjectsAccumulators(iDataStreamFile& ifile);
  const string typeName() { return "Linear mean-conditional Gaussian"; }
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
