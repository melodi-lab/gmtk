/*-
 * GMTK_MixGaussians
 *        Mixture of gaussians of various types.
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


#ifndef GMTK_MIXGAUSSIANS_H
#define GMTK_MIXGAUSSIANS_H

#include "fileParser.h"
#include "logp.h"
#include "machine-dependent.h"


#include "GMTK_MixGaussiansCommon.h"
#include "GMTK_GaussianComponent.h"

#include "GMTK_DiagGaussian.h"

#include "GMTK_Dense1DPMF.h"


class MixGaussians : public MixGaussiansCommon {

  vector < GaussianComponent* > components;
 
public:

  MixGaussians(const int dim) : MixGaussiansCommon(dim) { }
  ~MixGaussians() {}

  void preCompute();

  //////////////////////////////////////////////
  // read/write basic parameters
  void read(iDataStreamFile& is);
  void write(oDataStreamFile& os);


  //////////////////////////////////
  // set all current parameters to valid but random values
  void makeRandom() {}
  // set all current parameters to valid but "uniform" values 
  // (for Gaussians this means N(0,1))
  void makeUniform() {}
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
  void emStartIteration();
  void emIncrement(RandomVariable*,logpr prob);
  void emEndIteration();
  void emSwapCurAndNew();
  void emStoreAccumulators(oDataStreamFile& ofile);
  void emLoadAccumulators(iDataStreamFile& ifile);
  void emAccumulateAccumulators(iDataStreamFile& ifile);
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


