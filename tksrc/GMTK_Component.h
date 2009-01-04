/*-
 * GMTK_Component
 *        Component elements that are shared by all mixture components
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


#ifndef GMTK_COMPONENT_H
#define GMTK_COMPONENT_H

#include "fileParser.h"
#include "logp.h"

#include "machine-dependent.h"

#include "GMTK_NamedObject.h"
#include "GMTK_EMable.h"

class Component :  public EMable {

  friend class GMTK_Tie;

protected:
  ///////////////////////////////////////////////////////
  // dim = dimensionality of this component.
  // Typically, this will be the number of features over which
  // this component applies.
  const unsigned _dim;

  /////////////////////////////////////////////////
  // modify the usage counts of any members that use them (e.g. means,
  // covars); typically called with amount=1 or -1
  virtual void adjustNumTimesShared(int amount) = 0;

public:

  ///////////////////////////////////////////////////////////
  // This object must be aware of types of its derived
  // classes for reading from file purposes. 
  enum ComponentType { 
    // standard diagonal Gaussian
    DiagGaussian = 0, 
    // Linear mean-conditional diagonal Gaussian. 
    // Implements full-covariance Gaussians (via chol. factorization),
    // sparse covariance Gaussians with Bayesian network covariance matrices,
    // and linear BMMs.
    LinMeanCondDiagGaussian = 1,
    // Polynomial non-linear mean-conditional diagonal Gaussians 
    // Implements sparse covariance Gaussians with Bayesian network covariance matrices
    // but non-linear dependencies, and polynomial non-linear BMMs.
    PolyNLinMeanCondDiagGaussian = 2,
    // A N-dimensional Gamma distribution component, for real values that are
    // strictly > 0.
    GammaComponent = 3,
    // A N-dimensional Beta distribution component, for real values that are
    // in a given range
    BetaComponent = 4
    // add others here as they are written...
  };

  Component(const int dim) : _dim(dim) {}
  virtual ~Component() { }

  unsigned dim() const { return _dim; }


  //////////////////////////////////////////////
  // read/write basic parameters
  virtual void read(iDataStreamFile& is) = 0;
  virtual void write(oDataStreamFile& os) = 0;
  //////////////////////////////////////////////

  /////////////////////////////////////////////////
  // create a copy of self, with entirel or partially new parameters
  // (so clone might or might not share, depends on the object), and
  // with slightly (and randomly) perturbed values.
  virtual Component* noisyClone() = 0;

  /////////////////////////////////////////////////
  // create a copy of self, with entirely new parameters with
  // identical values; NOTHING is shared
  virtual Component* identicalIndependentClone() = 0;


  //////////////////////////////////
  // set all current parameters to valid but random values
  virtual void makeRandom() = 0;
  // set all current parameters to valid but "uniform" values 
  // (for Gaussians this means N(0,1))
  virtual void makeUniform() = 0;
  //////////////////////////////////

  virtual void recursivelyClearUsedBit() = 0;
  virtual void recursivelySetUsedBit() = 0;

  //////////////////////////////////
  // probability evaluation
  virtual logpr log_p(const float *const x,     // real-valued scoring obs at time t
		      const Data32* const base, // ptr to base obs at time t
		      const int stride) = 0;    // stride
  double p(const float *const x,
		   const Data32* const base,
		   const int stride)
  { return ::exp(log_p(x,base,stride).val()); }


  // return the maximum possible value of this component. Default
  // return value is 1.0 (which is not even a bound for continuous RVs
  // since the scores for, say, Gaussians can be > 1.0, so this function
  // should be redefined in child classes.
  virtual logpr maxValue() { return 1.0; }
  //////////////////////////////////


  //////////////////////////////////
  // Sample Generation            //
  //////////////////////////////////
  virtual void sampleGenerate(float *sample,
			      const Data32* const base,
			      const int stride) = 0;
  //////////////////////////////////

};


#endif
