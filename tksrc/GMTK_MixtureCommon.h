/*-
 * GMTK_GaussianCommon
 *        Common elements that are shared by all Gaussiaan type clases.
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


#ifndef GMTK_GAUSSIANCOMMON_H
#define GMTK_GAUSSIANCOMMON_H

#include "fileParser.h"
#include "logp.h"

#include "machine-dependent.h"

class GaussianCommon {

  ///////////////////////////////////////////////////////
  // dim = dimensionality of this Gaussian.
  // Typically, this will be the number of features over which
  // this Gaussian applies.
  int dim;

  ///////////////////////////////////////////////////////
  // The value that, if any variances go below, cause
  // either warnings (and adjustments) or errors to occur
  // depending on how probable this component is.
  static double _varianceFloor;

  /////////////////////////////////////////////////////////
  // precompute values so that evaluation is efficient.
  virtual void preCompute();

public:

  
  GaussianCommon() {}
  virtual ~GaussianCommon() {}


  static double varianceFloor() { return _varianceFloor; }
  static void setVarianceFloor(double floor) 
    { if (floor < FLT_MIN) floor = FLT_MIN; _varianceFloor = floor; }


  //////////////////////////////////////////////
  // read/write basic parameters
  virtual void read(iDataStreamFile& is);
  virtual void write(oDataStreamFile& os);
  //////////////////////////////////////////////


  //////////////////////////////////
  // set all current parameters to valid but random values
  virtual void makeRandom();
  // set all current parameters to valid but "uniform" values 
  // (for Gaussians this means N(0,1))
  virtual void makeUniform();
  //////////////////////////////////

  //////////////////////////////////
  // probability evaluation
  virtual logpr log_p(const float *const x,    // real-valued scoring obs at time t
		      const Data32* const base, // ptr to base obs at time t
		      const int stride);       // stride
  virtual double p(const float *const x,
		   const Data32* const base,
		   const int stride)
  { return exp(log_p(x,base,stride).val()); }
  //////////////////////////////////


  //////////////////////////////////
  // Full Baum-Welch EM training  //
  //////////////////////////////////
  virtual void emInit();
  virtual void startEmEpoch();
  // WARNING: Interface to this next one will soon change. JB: 4/19/01
  virtual void emAccumulate(const float prob,
			    const float *const obs);
  virtual void endEmEpoch(logpr cmpSop_acc);
  // For parallelism, loading/storing partially completed accumulators.
  virtual void emLoadAccumulators(iDataStreamFile& ifile);
  virtual void emStoreAccumulators(oDataStreamFile& ofile);
  virtual void emAccumulateAccumulators(iDataStreamFile& ifile);
  //////////////////////////////////


  //////////////////////////////////
  // Sample Generation            //
  //////////////////////////////////
  virtual void sampleGenerate(float *const sample,
		      const Data32* const base);
  //////////////////////////////////


};


#endif
