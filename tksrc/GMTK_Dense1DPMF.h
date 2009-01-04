/*-
 * GMTK_Dense1DPMF.h
 *      .h file for GMTK_Dense1DPMF.h, trainable 1D discrete probability
 *      distributions.
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


#ifndef GMTK_DENSE1DPMF_H
#define GMTK_DENSE1DPMF_H


#include "fileParser.h"
#include "logp.h"
#include "sArray.h"

#include "GMTK_EMable.h"
#include "GMTK_NamedObject.h"
#include "GMTK_DirichletTable.h"
#include "GMTK_DirichletPrior.h"
#include "tieSupport.h"


class Dense1DPMF : public EMable, public DirichletPrior {

  friend class GMTK_Tie;

  ///////////////////////////////////////////////////////////  
  // The probability mass function
  sArray <logpr> pmf;
  ///////////////////////////////////////////////////////////  

  //////////////////////////////////
  // Data structures support for EM.
  //   The next probability mass function 
  sArray <logpr> nextPmf;
  ///////////////////////////////////////////////////////////  

  logpr cachedMaxValue;

public:

  ///////////////////////////////////////////////////////////  
  // General constructor
  Dense1DPMF();
  ~Dense1DPMF() { }

  /////////////////////////////////////////////////
  // create a copy of self, with entirely new parameters with
  // identical values; NOTHING is shared
  Dense1DPMF* identicalIndependentClone();

  unsigned length() { return (unsigned)pmf.len(); }
  unsigned card() { return (unsigned)pmf.len(); }

  logpr p(unsigned i) { 
    assert ( i < (unsigned)pmf.len() );
    return pmf.ptr[i]; 
  }
  
  // also give access to next pmf
  logpr np(unsigned i) { 
    assert ( i < (unsigned)pmf.len() );
    return nextPmf[i]; 
  }

  ///////////////////////////////////////////////////////////  
  // Re-normalize the output distributions
  void normalize();
  // set all values to random values.
  void makeRandom();
  // set all values to uniform values.
  void makeUniform();


  ///////////////////////////////////
  unsigned totalNumberParameters() { return pmf.len(); }

  logpr maxValue() { return cachedMaxValue; }


  ///////////////////////////////////////////////////////////    
  // read in the basic parameters, assuming file pointer 
  // is located at the correct position.
  void read(iDataStreamFile& is);

  ///////////////////////////////////////////////////////////    
  // write out the basic parameters, starting at the current
  // file position.
  void write(oDataStreamFile& os);


  //////////////////////////////////
  // Public interface support for EM
  //////////////////////////////////
  void emStartIteration();
  void emIncrement(logpr prob,sArray<logpr>&);
  void emIncrement(logpr prob,const int val);
  void emEndIteration();
  void emSwapCurAndNew();

  // parallel training
  void emStoreObjectsAccumulators(oDataStreamFile& ofile,
				  bool writeLogVals = true,
				  bool writeZeros = false);
  void emLoadObjectsDummyAccumulators(iDataStreamFile& ifile);
  void emZeroOutObjectsAccumulators();
  void emLoadObjectsAccumulators(iDataStreamFile& ifile);
  void emAccumulateObjectsAccumulators(iDataStreamFile& ifile);
  const string typeName() { return "DPMF"; }
  //////////////////////////////////


};



#endif
