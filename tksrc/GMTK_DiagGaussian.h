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


#ifndef DIAGGAUSSIAN_H
#define DIAGGAUSSIAN_H

#include "fileParser.h"
#include "logp.h"

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

  // A bitmask giving the state of the object.
  enum {
    bm_basicAllocated = 0x1,
    bm_emAllocated = 0x2,
    bm_swapped = 0x4,

    bm_initState = 0x0
  };
  unsigned int bitmask;

 
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
  // swap the old and the new parameters.
  void swapCurAndNew();
  //////////////////////////////////

  //////////////////////////////////
  // probability evaluation
  logpr log_p(const float *const x,    // real-valued scoring obs at time t
	      const ptr32* const base, // ptr to base obs at time t
	      const int stride);       // stride
  //////////////////////////////////


  //////////////////////////////////
  // Full Baum-Welch EM training  //
  //////////////////////////////////
  virtual void emInit();
  virtual void startEmEpoch();
  // WARNING: Interface to this will soon change. JB: 4/19/01
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
  void sampleGenerate(float *const sample,
		      const ptr32* const base);
  //////////////////////////////////


};


#endif
