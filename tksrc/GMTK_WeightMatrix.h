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
#include "GMTK_NamedObject.h"
#include "GMTK_RandomVariable.h"


class WeightMatrix : public EMable, public NamedObject  {


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
  ~WeightMatrix() {}

  int rows() { return _rows; }
  int cols() { return _cols; }
  int nInputs() { return _cols; }
  int nOutputs() { return _rows; }

  ///////////////////////////////////////////////////////////    
  void read(iDataStreamFile& is);

  ///////////////////////////////////////////////////////////    
  void write(oDataStreamFile& os);

  ///////////////////////////////////////////////////////////    
  void makeRandom() {}
  void makeUniform() {}

  //////////////////////////////////
  // Public interface support for EM
  //////////////////////////////////
  void emStartIteration() {}
  void emIncrement(RandomVariable*,logpr prob) {}
  void emEndIteration() {}
  void emSwapCurAndNew() {}
  void emStoreAccumulators(oDataStreamFile& ofile) {}
  void emLoadAccumulators(iDataStreamFile& ifile) {}
  void emAccumulateAccumulators(iDataStreamFile& ifile) {}
  //////////////////////////////////



};



#endif 
