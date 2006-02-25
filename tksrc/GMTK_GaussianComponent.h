/*-
 * GMTK_GaussianComponent
 *        Gaussian Component elements that are shared by all Gaussian type component clases.
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


#ifndef GMTK_GAUSSIANCOMPONENT_H
#define GMTK_GAUSSIANCOMPONENT_H

#include "fileParser.h"
#include "logp.h"

#include "machine-dependent.h"

#include "GMTK_NamedObject.h"
#include "GMTK_Component.h"

class GaussianComponent :  public Component {

  ///////////////////////////////////////////////////////
  // The value that, if any variances go below, cause
  // either warnings (and adjustments) or errors to occur
  // depending on how probable this component is.
  static double _varianceFloor;

public:


  GaussianComponent(const int dim);
  virtual ~GaussianComponent() { }

  static double varianceFloor() { return _varianceFloor; }
  static double setVarianceFloor(const double floor);

  static bool cloneShareMeans;
  static bool cloneShareCovars;
  static bool cloneShareDlinks;


  //////////////////////////////////////////////
  // read/write basic parameters
  virtual void read(iDataStreamFile& is) = 0;
  virtual void write(oDataStreamFile& os) = 0;
  //////////////////////////////////////////////

  /////////////////////////////////////////////////
  // create a copy of self, with entirely new parameters
  // (so clone shares nothing), and with slightly (and 
  // randomly) perturbed values.
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
