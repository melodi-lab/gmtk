/*-
 * GMTK_RealArray.h
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


#ifndef GMTK_MEANVECTOR_H
#define GMTK_MEANVECTOR_H

#include "fileParser.h"
#include "logp.h"
#include "sArray.h"

#include "GMTK_RealArray.h"
#include "GMTK_EMable.h"
#include "GMTK_NamedObject.h"


class DiagCovarVector;
class DlinkMatrix;

class MeanVector : public EMable {

private:

  friend class DiagGaussian;
  friend class LinMeanCondDiagGaussian;
  friend class DiagCovarVector;
  friend class DlinkMatrix;

  //////////////////////////////////
  // The acutal mean vector
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

  void emEndIterationNoSharingAlreadyNormalized(const float *const accumulatedNextMeans);
  void emEndIterationSharedMeansCovarsDlinks(const logpr accumulatedProbability,
					     const float*const xAccumulators,
					     const float*const zAccumulators,
					     const DlinkMatrix* dLinkMat,
					     const DiagCovarVector* covar);

  void emSwapCurAndNew();


  // parallel training
  void emStoreObjectsAccumulators(oDataStreamFile& ofile) {};
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
