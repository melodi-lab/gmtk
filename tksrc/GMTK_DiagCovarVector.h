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

class DiagCovarVector : public EMable {


  //////////////////////////////////////////////////////  
  // the name
  char *_name;

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
  ~DiagCovarVector() { delete [] _name; }

  ///////////////////////////////////////////////////////////  
  void makeRandom();

  //////////////////////////////////////////////
  // read/write basic parameters
  void read(iDataStreamFile& is);
  void write(oDataStreamFile& os) { 
    os.write(_name,"DiagCovarVector::write name");
    covariances.write(os); 
  }



  //////////////////////////////////
  // Public interface support for EM
  //////////////////////////////////
  void emInit() {}
  void startEmEpoch() {}
  void emAccumulate(const float prob,
		    const float *const oo_array) {}
  void endEmEpoch(logpr cmpSop_acc) {}
  void emLoadAccumulators(iDataStreamFile& ifile) {}
  void emStoreAccumulators(oDataStreamFile& ofile) {}
  void emAccumulateAccumulators(iDataStreamFile& ifile) {}
  void swapCurAndNew() {}
  //////////////////////////////////

};



#endif 
