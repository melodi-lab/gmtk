/*-
 * GMTK_DiagCovarVector.h
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


#ifndef GMTK_DIAGCOVARVECTOR_H
#define GMTK_DIAGCOVARVECTOR_H

#include "fileParser.h"
#include "logp.h"
#include "sArray.h"

#include "GMTK_EMable.h"
#include "GMTK_RealArray.h"
#include "GMTK_NamedObject.h"
#include "GMTK_RandomVariable.h"

class DiagCovarVector : public EMable, public NamedObject {


  //////////////////////////////////
  // The actual covariance "matrix"
  RealArray covariances;

  //////////////////////////////////
  // Data structures support for EM
  sArray<float> nextCovariances;
  // Need the following means because of potential sharing.
  // that is, if multiple Gaussians share different means
  // and different covariances, then the covariance matrices
  // need to keep track of the accumulated means according to
  // all the Gaussians that share this Covariance matrix (which
  // might not be the same sharing as the means.
  sArray<float> nextMeans;

  ///////////////////////////////////////////////////////
  // Precomputed values for efficiency (computed at the
  // top of every EM epoch).
  //    Precomputed inverse variances.
  sArray<float> variances_inv;
  //    precomputed logged normalization constant.
  float _log_inv_normConst;
  ///////////////////////////////////////////////////////

  /////////////////////////////////////////////////
  // counts the number of gaussian components
  // that are sharing this covariance
  unsigned refCount;

  /////////////////////////////////////////////////////////////
  // used by EM to count the number of times variances 
  // became very small.
  static unsigned numFlooredVariances;

public:

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
  void write(oDataStreamFile& os) { 
    NamedObject::write(os);
    covariances.write(os); 
  }

  // create a copy of self, but with slightly perturbed values
  DiagCovarVector* noisyClone();

  //////////////////////////////////
  // Public interface support for EM
  //////////////////////////////////
  void emStartIteration(sArray<float>&);
  void emIncrement(const logpr prob,
		   const float fprob,
		   const float* f,
		   const Data32* const base,
		   const int stride,
		   float *const partialAccumulatedNextCovars);
  void emEndIteration(const float*const c);
  void emEndIteration(const logpr prob,const float*const m,const float*const v);
  void emSwapCurAndNew();
  void emStoreAccumulators(oDataStreamFile& ofile);
  void emLoadAccumulators(iDataStreamFile& ifile);
  void emAccumulateAccumulators(iDataStreamFile& ifile);
  //////////////////////////////////

};



#endif 
