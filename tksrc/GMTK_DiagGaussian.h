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


#include "GMTK_RandomVariable.h"
#include "GMTK_GaussianCommon.h"
#include "GMTK_EMable.h"
#include "GMTK_MeanVector.h"
#include "GMTK_DiagCovarVector.h"

class DiagGaussian : public GaussianCommon, public EMable {

  ///////////////////////////////////////////////////////
  // The means. 
  // This might be tied with multiple other distributions.
  MeanVector* means;
  // The index in the global mean array of this mean.
  int meanIndex; 

  ///////////////////////////////////////////////////////
  // The diagonal of the covariance matrix
  // This might be tied with multiple other distributions.
  DiagCovarVector* variance;
  // The index in the global variance array of this variance vector
  int varianceIndex;

 
public:
  
  DiagGaussian() { bitmask = 0x0; }
  ~DiagGaussian();

  //////////////////////////////////////////////
  // read/write basic parameters
  void read(iDataStreamFile& is);
  void write(oDataStreamFile& os);


  //////////////////////////////////
  // set all current parameters to valid but random values
  void makeRandom();
  // set all current parameters to valid but "uniform" values 
  // (for Gaussians this means N(0,1))
  void makeUniform();
  //////////////////////////////////

  //////////////////////////////////
  // probability evaluation
  logpr log_p(const float *const x,    // real-valued scoring obs at time t
	      const Data32* const base, // ptr to base obs at time t
	      const int stride);       // stride
  //////////////////////////////////


  //////////////////////////////////
  // Full Baum-Welch EM training  //
  //////////////////////////////////
  void emStartIteration() {}
  void emIncrement(RandomVariable*,logpr prob) {}
  void emEndIteration() {}
  void emSwapCurAndNew() {}
  void emStoreAccumulators(oDataStreamFile& ofile) {}
  void emLoadAccumulators(iDataStreamFile& ifile) {}
  void emAccumulateAccumulators(iDataStreamFile& ifile) {}
  //////////////////////////////////

  //////////////////////////////////
  // Sample Generation            //
  //////////////////////////////////
  void sampleGenerate(float *const sample,
		      const Data32* const base);
  //////////////////////////////////


};


#endif
