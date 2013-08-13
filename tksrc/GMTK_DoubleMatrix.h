/*-
 * GMTK_DoubleMatrix.h
 *
 *  Written by Jeff Bilmes <bilmes@ee.washington.edu>
 * 
 * Copyright (C) 2013 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 *
 */


#ifndef GMTK_DOUBLEMATRIX_H
#define GMTK_DOUBLEMATRIX_H

#include "fileParser.h"
#include "logp.h"
#include "sArray.h"

#include "GMTK_EMable.h"
#include "GMTK_NamedObject.h"

class RV;

class DoubleMatrix : public EMable  {

  friend class GammaComponent;
  friend class BetaComponent;
  friend class DeepVECPT;

protected:

  ///////////////////////////////////////////////////////////  
  // The data values
  sArray<double> values;
  ///////////////////////////////////////////////////////////  

  //////////////////////////////////
  // Data structures support for EM
  sArray<double> nextValues;

  ///////////////////////////////////////////////////////////  
  // dimensons of the distribution (separate from the memory
  // allocated for arr above).
  int _rows;
  int _cols;
  ///////////////////////////////////////////////////////////  

  /////////////////////////////////////////////////
  // counts the number of objects that use this real matrix for
  // any purpose. This is a static count, computed as
  // the objects are read in.
  unsigned numTimesShared;


  /////////////////////////////////////////////////
  // counts the number of objects that are using this real matrix at
  // EM training time.  This is a dynamic count, and is computed as EM
  // training is run. This value does not necessarily equal the number
  // of objects that have specified this object in the data files.
  unsigned refCount;


public:

  ///////////////////////////////////////////////////////////  
  // General constructor
  DoubleMatrix();
  ~DoubleMatrix() {}

  int rows() { return _rows; }
  int cols() { return _cols; }

  // for use when using this marix 
  // int nInputs() { return _cols; }
  // int nOutputs() { return _rows; }

  ///////////////////////////////////////////////////////////    
  void read(iDataStreamFile& is);

  ///////////////////////////////////////////////////////////    
  void write(oDataStreamFile& os);

  // create an exact copy of self
  virtual DoubleMatrix* cleanClone();

  unsigned totalNumberParameters() { return _rows*_cols; }
  void recursivelyClearUsedBit() {  emClearUsedBit();  }
  void recursivelySetUsedBit() { emSetUsedBit();  }

  // The following are dummy routine values for now. These
  // em routines will be implemented by the objects that 
  // use this marix.

  //////////////////////////////////
  // Public interface support for EM
  //////////////////////////////////

  ///////////////////////////////////////////////////////////    
  // since multiple objects might use a matrix, we let the user of
  // these object types define how to make random and uniform.
  virtual void makeRandom() {}
  virtual void makeUniform() {}

  // Give only the most basic defintions, sub-classes can fill in.
  virtual void emStartIteration();
  virtual void emIncrement(logpr prob);
  virtual void emEndIteration();
  virtual void emSwapCurAndNew(); 

  // parallel training. Since accumulators might be very different
  // for different users of this object, we leave saving accumulators to the
  // user.
  virtual void emStoreObjectsAccumulators(oDataStreamFile& ofile,
				  bool writeLogVals = true,
				  bool writeZeros = false) {}
  virtual void emLoadObjectsDummyAccumulators(iDataStreamFile& ifile) {}
  virtual void emZeroOutObjectsAccumulators() {}
  virtual void emLoadObjectsAccumulators(iDataStreamFile& ifile) {}
  virtual void emAccumulateObjectsAccumulators(iDataStreamFile& ifile) {}
  virtual const string typeName() { return "Double Matrix"; }
  //////////////////////////////////

};



#endif 
