/*-
 * GMTK_Weightmatrix.h
 *
 *  Written by Jeff Bilmes <bilmes@ee.washington.edu>
 * 
 *  $Header$
 * 
 * Copyright (C) 2001 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 *
 */


#ifndef GMTK_WEIGHTMATRIX_H
#define GMTK_WEIGHTMATRIX_H

#include "fileParser.h"
#include "logp.h"
#include "sArray.h"

#include "GMTK_EMable.h"
#include "GMTK_NamedObject.h"

class RV;

class WeightMatrix : public EMable  {

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

  unsigned totalNumberParameters() { return 0; }
  void recursivelyClearUsedBit() {  emClearUsedBit();  }
  void recursivelySetUsedBit() { emSetUsedBit();  }


  //////////////////////////////////
  // Public interface support for EM
  //////////////////////////////////
  void emStartIteration() {}
  void emIncrement(RV*,logpr prob) {}
  void emEndIteration() {}
  void emSwapCurAndNew() {}

  // parallel training
  void emStoreObjectsAccumulators(oDataStreamFile& ofile,
				  bool writeLogVals = true,
				  bool writeZeros = false) {}
  void emLoadObjectsDummyAccumulators(iDataStreamFile& ifile) {}
  void emZeroOutObjectsAccumulators() {}
  void emLoadObjectsAccumulators(iDataStreamFile& ifile) {}
  void emAccumulateObjectsAccumulators(iDataStreamFile& ifile) {}
  const string typeName() { return "Weight matrix"; }
  //////////////////////////////////

};



#endif 
