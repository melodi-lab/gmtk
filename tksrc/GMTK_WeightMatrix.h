/*-
 * GMTK_Weightmatrix.h
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


#ifndef GMTK_WEIGHTMATRIX_H
#define GMTK_WEIGHTMATRIX_H

#include "fileParser.h"
#include "logp.h"
#include "sArray.h"

#include "GMTK_EMable.h"

class WeightMatrix : public EMable {


  ///////////////////////////////////////////////////////////  
  // The data values
  sArray<float> weights;
  ///////////////////////////////////////////////////////////  

  //////////////////////////////////
  // Data structures support for EM
  sArray<float> nextWeights;


  ///////////////////////////////////////////////////////////  
  // dimensons of the distribution (separate from the memory
  // allocated for arr above).
  int _rows; // == number of outputs
  int _cols; // == number of inputs
  ///////////////////////////////////////////////////////////  

public:

  ///////////////////////////////////////////////////////////  
  // General constructor
  WeightMatrix();

  int rows() { return _rows; }
  int cols() { return _cols; }
  int nInputs() { return _cols; }
  int nOutputs() { return _rows; }

  ///////////////////////////////////////////////////////////    
  void read(iDataStreamFile& is);

  ///////////////////////////////////////////////////////////    
  void write(oDataStreamFile& os);


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
