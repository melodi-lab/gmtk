/*-
 * GMTK_DiagCovarVector.h
 *
 *  Written by Jeff Bilmes <bilmes@ee.washington.edu>
 * 
 *  $Header$
 * 
 * Copyright (C) 2001 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 *
 */


#ifndef GMTK_DIAGCOVARVECTOR_H
#define GMTK_DIAGCOVARVECTOR_H

#include <list>

#include "fileParser.h"
#include "logp.h"
#include "sArray.h"

#include "GMTK_EMable.h"
#include "GMTK_NamedObject.h"

class MeanVector;
class Clusterable;
class DlinkMatrix;


class DiagCovarVector : public EMable {

  friend class DiagGaussian;
  friend class LinMeanCondDiagGaussian;
  friend class MissingFeatureScaledDiagGaussian;
  friend class MeanVector;
  friend class DlinkMatrix;
  friend class GMTK_Tie;
  friend double cluster_scaled_log_likelihood(std::list<Clusterable*> &items, double* tot_occupancy);

  //////////////////////////////////
  // The actual covariance "matrix"
  sArray<float> covariances;

  /////////////////////////////////////////////////
  // counts the number of gaussian components
  // that are sharing this mean at all. This is a static
  // count, and is computed as any object that
  // might use a mean (such as a Gaussian component)
  // is read in.
  unsigned numTimesShared;

  // Precomputed values for efficiency (computed at the
  // top of every training EM epoch, or during a read-in).
  // Precomputed inverse variances.
  sArray<float> variances_inv;
  // precomputed logged normalization constant.
  float _log_inv_normConst;
  ///////////////////////////////////////////////////////


  ////////////////////////////////////////////////////////////////////
  // Data structures support for EM
  ////////////////////////////////////////////////////////////////////
  struct Training_Members {
    // NOTE: if we change from type float, to type double,
    //   we will need to check the load/store accumulator code
    //   for subclasses of GaussianComponent.
    sArray<float> nextCovariances;

    ///////////////////////////////////////////////////////

    /////////////////////////////////////////////////
    // counts the number of gaussian components
    // that are sharing this covariance at EM training time. This is a dynamic
    // count, and is computed as EM training is run. This
    // value does not necessarily equal the number of
    // objects that have specified this object in
    // the object files.
    unsigned refCount;

    /////////////////////////////////////////////////////////////
    // used by EM to count the number of times variances 
    // became very small.
    unsigned numFlooredVariances;


    ////////////////////////////////////////////////
    // the local accumulated probability needed for each element of
    // the mean vector in a missing feature scaled gaussian. This 
    // is currently (as of 6/21/2013) used only by MissingFeatureScaledDiagGaussian,
    // for other uses this array remains empty.
    sArray<logpr> elementAccumulatedProbability;

    // default constructor
    Training_Members() : refCount(0),numFlooredVariances(0) {}  
  };
  auto_deleting_ptr<struct Training_Members> trMembers;
  ////////////////////////////////////////////////////////////////////
  // End of data structures support for EM
  ////////////////////////////////////////////////////////////////////



public:

  /////////////////////////////////////////////////////////////
  // Floor the variances to the variance floor when they are read in.
  static bool floorVariancesWhenReadIn;

  ///////////////////////////////////////////////////////////  
  // General constructor
  DiagCovarVector();
  ~DiagCovarVector() {}

  // When noisy cloning an object, this gives
  // the fraction to multiply to get the STD of the noise.
  static double cloneSTDfrac;
  static void checkForValidValues();

  ///////////////////////////////////////////////////////////  
  void makeRandom();
  void makeUniform();


  ///////////////////////////////////
  unsigned totalNumberParameters() { return covariances.len(); }
  void recursivelyClearUsedBit() {  emClearUsedBit();  }
  void recursivelySetUsedBit() { emSetUsedBit();  }

  void initElementAccumulatedProbability() {
    assert ( emEmAllocatedBitIsSet() );
    trMembers->elementAccumulatedProbability.growIfNeeded(covariances.len());
    for (int i=0 ; i<covariances.len() ; i++) {
      trMembers->elementAccumulatedProbability[i].set_to_zero();
    }
  }



  ///////////////////////////////////////
  int dim() { return covariances.len(); }

  // easy access to internal pointers for efficient reading.
  const float* basePtr() { return &covariances[0]; }
  const float* baseVarInvPtr() { return &variances_inv[0]; }
  float log_inv_normConst() { return _log_inv_normConst; }

  ///////////////////////////////////////
  // preCompute: should be called anytime the parameters
  // change.
  void preCompute();

  //////////////////////////////////////////////
  // read/write basic parameters
  void read(iDataStreamFile& is);
  void write(oDataStreamFile& os);

  // create a copy of self, but with slightly perturbed values
  DiagCovarVector* noisyClone();

  /////////////////////////////////////////////////
  // create a copy of self, with entirely new parameters with
  // identical values; NOTHING is shared
  DiagCovarVector* identicalIndependentClone();
  
  //////////////////////////////////
  // Public interface support for EM
  //////////////////////////////////
  void emStartIteration(sArray<float>&);
  using EMable::emIncrement;
  void emIncrement(const logpr prob,
		   const float fprob,
		   const float* f,
		   const Data32* const base,
		   const int stride,
		   float *const partialAccumulatedNextCovars);
  using EMable::emEndIteration;
  void emEndIteration(const float*const c);
  void emEndIteration(const logpr prob,const float*const m,const float*const v);
  void emEndIterationNoSharing(const float*const m,const float*const v);
  void emEndIterationNoSharingElementProbabilities(const float*const m,
						   const float*const v,
						   const logpr *const elementAccumulatedProbabilities);

  void emEndIterationNoSharingAlreadyNormalized(const float*const c);

  void emEndIterationSharedMeansCovars(const logpr parentsAccumulatedProbability,
				       const float*const partialAccumulatedNextMeans,
				       const float *const partialAccumulatedNextCovar,
				       const MeanVector* mean);
  void emEndIterationSharedMeansCovarsElementProbabilities(const logpr parentsAccumulatedProbability,
				       const float*const partialAccumulatedNextMeans,
				       const float *const partialAccumulatedNextCovar,
				       const MeanVector* mean,
				       const logpr *const elementAccumulatedProbabilities);


  void emEndIterationSharedMeansCovarsDlinks(const logpr accumulatedProbability,
					     const float* const xAccumulators,
					     const float* const xxAccumulators,
					     const float* const xzAccumulators,
					     const float* const zAccumulators,
					     const float* const zzAccumulators,
					     const MeanVector* mean,
					     const DlinkMatrix* dLinkMat);

  void emEndIterationSharedCovars(const logpr parentsAccumulatedProbability,
				  const float*const partialAccumulatedNextMeans,
				  const float *const partialAccumulatedNextCovar);
  void emEndIterationSharedCovarsElementProbabilities(
		  const logpr parentsAccumulatedProbability,
		  const float*const partialAccumulatedNextMeans,
		  const float *const partialAccumulatedNextCovar,
		  const logpr *const elementAccumulatedProbabilities);


  void emEndIterationSharedCovars(const float*const c);
  void emSwapCurAndNew();


  // parallel training
  void emStoreObjectsAccumulators(oDataStreamFile& ofile,
				  bool writeLogVals = true,
				  bool writeZeros = false) {};
  void emLoadObjectsDummyAccumulators(iDataStreamFile& ifile) {};
  void emZeroOutObjectsAccumulators() {};
  void emLoadObjectsAccumulators(iDataStreamFile& ifile) {};
  void emAccumulateObjectsAccumulators(iDataStreamFile& ifile) {};
  const string typeName() { return "covariance vector"; }
  // need to override parent class's routine in this case.
  void emStoreAccumulators(oDataStreamFile& ofile);

  //////////////////////////////////

};



#endif 
