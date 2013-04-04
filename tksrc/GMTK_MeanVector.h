/*-
 * GMTK_MeanVector.h
 *
 *  Written by Jeff Bilmes <bilmes@ee.washington.edu>
 * 
 *  $Header$
 * 
 * Copyright (C) 2001 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 *
 *
 */


#ifndef GMTK_MEANVECTOR_H
#define GMTK_MEANVECTOR_H

#include <list>

#include "fileParser.h"
#include "logp.h"
#include "sArray.h"

#include "GMTK_EMable.h"
#include "GMTK_NamedObject.h"


class DiagCovarVector;
class DlinkMatrix;
class Clusterable;

class MeanVector : public EMable {

private:

  friend class DiagGaussian;
  friend class LinMeanCondDiagGaussian;
  friend class MissingFeatureScaledDiagGaussian;
  friend class DiagCovarVector;
  friend class DlinkMatrix;
  friend class GMTK_Tie;
  friend class ClusterableMean;
  friend class ClusterableMixture;
  friend double cluster_scaled_log_likelihood(std::list<Clusterable*> &items, double* tot_occupancy);

  //////////////////////////////////
  // The actual mean vector
  sArray<float> means;

  //////////////////////////////////
  // Data structures support for EM
  // NOTE: if we change from type float, to type double,
  //   we will need to check the load/store accumulator code
  //   for subclasses of GaussianComponent.
  sArray<float> nextMeans;

  /////////////////////////////////////////////////////////
  // the denominator used in the case of shared means.
  sArray<double> sharedMeansDenominator;

  /////////////////////////////////////////////////
  // counts the number of gaussian components
  // that are sharing this mean at all. This is a static
  // count, and is computed as any object that
  // might use a mean (such as a Gaussian component)
  // is read in.
  unsigned numTimesShared;

  /////////////////////////////////////////////////
  // counts the number of gaussian components
  // that are sharing this mean at EM training time. This is a dynamic
  // count, and is computed as EM training is run. This
  // value does not necessarily equal the number of
  // objects that have specified this object in
  // the object files.
  unsigned refCount;

  /////////////////////////////////////////////////
  // allow access to internal accumulator pointer for
  // shared covariance object to use to add to its
  // accumulation.
  const float* const accumulatorPtr() { return nextMeans.ptr; }

public:

  ///////////////////////////////////////////////////////////  
  // General constructor
  MeanVector();
  ~MeanVector() { } 

  // When noisy cloning an object, this gives
  // the fraction to multiply to get the STD of the noise.
  static double cloneSTDfrac;
  static void checkForValidValues();

  //////////////////////////////////
  // set all current parameters to random values
  void makeRandom();
  void makeUniform();

  ///////////////////////////////////////
  int dim() { return means.len(); }

  const float *basePtr() { return &means[0]; }

  //////////////////////////////////////////////
  // read/write basic parameters
  void read(iDataStreamFile& is);
  void write(oDataStreamFile& os);

  // create a copy of self, but with slightly perturbed
  // means values.
  MeanVector* noisyClone();

  /////////////////////////////////////////////////
  // create a copy of self, with entirely new parameters with
  // identical values; NOTHING is shared
  MeanVector* identicalIndependentClone();

  unsigned totalNumberParameters() { return means.len(); }
  void recursivelyClearUsedBit() {  emClearUsedBit();  }
  void recursivelySetUsedBit() { emSetUsedBit();  }



  //////////////////////////////////
  // Public interface support for EM
  //////////////////////////////////
  void emStartIteration(sArray<float>& componentsNextMeans);
  void emIncrement(const logpr prob,
		   const float fprob,
		   const float *f,
		   const Data32* const base,
		   const int stride,
		   float *const partialAccumulatedNextMeans);
  void emEndIteration(const float *const partialAccumulatedNextMeans);
  void emEndIterationSharedMeansCovars(const logpr parentsAccumulatedProbability,
					     const float*const partialAccumulatedNextMeans,
					     const DiagCovarVector* covar);
  void emEndIterationNoSharing(const float *const partialAccumulatedNextMeans);

  void emEndIterationNoSharingElementProbabilities(const float *const partialAccumulatedNextMeans,
						   const logpr *const elementAccumulatedProbabilities);

  void emEndIterationNoSharingAlreadyNormalized(const float *const accumulatedNextMeans);
  void emEndIterationSharedMeansCovarsDlinks(const logpr accumulatedProbability,
					     const float*const xAccumulators,
					     const float*const zAccumulators,
					     const DlinkMatrix* dLinkMat,
					     const DiagCovarVector* covar);

  void emSwapCurAndNew();


  // parallel training
  void emStoreObjectsAccumulators(oDataStreamFile& ofile,
				  bool writeLogVals = true,
				  bool writeZeros = false) {};
  void emLoadObjectsDummyAccumulators(iDataStreamFile& ifile) {};
  void emZeroOutObjectsAccumulators() {};
  void emLoadObjectsAccumulators(iDataStreamFile& ifile) {};
  void emAccumulateObjectsAccumulators(iDataStreamFile& ifile) {};
  const string typeName() { return "mean vector"; }
  // need to override parent class's routine in this case.
  void emStoreAccumulators(oDataStreamFile& ofile);
  //////////////////////////////////


};



#endif // defined MEANVECTOR_H
