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

class DiagCovarVector : public EMable {


  //////////////////////////////////
  // The actual covariance "matrix"
  RealArray covariances;

  //////////////////////////////////
  // Data structures support for EM
  RealArray nextCovariances;

  ///////////////////////////////////////////////////////
  // Precomputed values for efficiency (computed at the
  // top of every EM epoch).
  //    Precomputed inverse variances.
  sArray<float> variances_inv;
  //    precomputed logged normalization constant.
  float log_inv_normConst;
  ///////////////////////////////////////////////////////


public:

  ///////////////////////////////////////////////////////////  
  // General constructor
  DiagCovarVector();


  //////////////////////////////////////////////
  // read/write basic parameters
  void read(iDataStreamFile& is) { 
    covariances.read(is); 
    for (i=0;i<covariances.len();i++) {
      if (covariances[i] < GaussianCommon::varianceFloor()) {
	error("DiagCovarVector:: read, covariance[%d] = (%e) < current Floor = (%e)",
	      i,covariances[i],GaussianCommon::varianceFloor());
      }
    }
  }
  void write(oDataStreamFile& os) { covariances.write(os); }

};



#endif 
