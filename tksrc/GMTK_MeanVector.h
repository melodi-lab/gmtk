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
#include "GMTK_RandomVariable.h"


class MeanVector : public EMable, public NamedObject {

private:

  friend class DiagGaussian;

  //////////////////////////////////
  // The acutal mean vector
  RealArray means;

  //////////////////////////////////
  // Data structures support for EM
  sArray<float> nextMeans;

  /////////////////////////////////////////////////
  // counts the number of gaussian components
  // that are sharing this mean.
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

  //////////////////////////////////
  // set all current parameters to random values
  void makeRandom();
  void makeUniform();

  ///////////////////////////////////////
  int dim() { return means.len(); }

  const float *basePtr() { return &means[0]; }

  //////////////////////////////////////////////
  // read/write basic parameters
  void read(iDataStreamFile& is) { 
    NamedObject::read(is);
    means.read(is); 
    setBasicAllocatedBit();
  }
  void write(oDataStreamFile& os) { 
    NamedObject::write(os);
    means.write(os); 
  }


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
  void emSwapCurAndNew();
  void emStoreAccumulators(oDataStreamFile& ofile);
  void emLoadAccumulators(iDataStreamFile& ifile);
  void emAccumulateAccumulators(iDataStreamFile& ifile);
  //////////////////////////////////


};



#endif // defined MEANVECTOR_H
