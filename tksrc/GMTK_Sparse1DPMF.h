/*-
 * GMTK_Dense1DPMF.h
 *      .h file for GMTK_Sparse1DPMF.h, trainable 1D discrete probability
 *      distributions. This is just a dense1dpmf but with values which
 *      say which rv values are not zero. Zero values are not represented.
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


#ifndef GMTK_SPARSE1DPMF_H
#define GMTK_SPARSE1DPMF_H


#include "fileParser.h"
#include "logp.h"
#include "sArray.h"

#include "GMTK_EMable.h"
#include "GMTK_NamedObject.h"
#include "GMTK_RandomVariable.h"
#include "GMTK_Dense1DPMF.h"


class Sparse1DPMF : public EMable, public NamedObject {


  //////////////////////////////////////
  // The cardinality of this RV, i.e., 
  // values may take value between [0:card-1]
  int _card;

  ///////////////////////////////////////////
  // the (possibly) shared 1DPMFs used for the prob values
  Dense1DPMF* dense1DPMF;

  ///////////////////////////////////////////////////////////  
  // The probability mass function
  sArray <RandomVariable::DiscreteVariableType> pmf;
  ///////////////////////////////////////////////////////////  

  //////////////////////////////////
  // Data structures support for EM.
  //   The previous probability mass function 
  //   (we don't use 'Entrys' here as the vals are the same.
  // sArray <logpr> nextPmf;
  ///////////////////////////////////////////////////////////  

public:

  ///////////////////////////////////////////////////////////  
  // General constructor
  Sparse1DPMF();
  ~Sparse1DPMF() {}

  //////////////////////////////////////
  // Return the number of valid values
  int length() { return pmf.len(); }

  ///////////////////////////////////////
  // return the cardinality of the RV
  int card() { return _card; }

  ///////////////////////////////////////
  // return the probability that for value 'val'
  // This will do a search for RV value val
  logpr prob(const int val);

  ///////////////////////////////////////
  // return the probability at entry i
  logpr probAtEntry(const int i) {
    assert ( i >= 0 && i < pmf.len() );
    return dense1DPMF->p(i);
  }
  ///////////////////////////////////////
  // return the value at entry i
  RandomVariable::DiscreteVariableType valueAtEntry(const int i) {
    assert ( i >= 0 && i < pmf.len() );
    return pmf[i];
  }

  ///////////////////////////////////////////////////////////  
  // Re-normalize the output distributions
  void normalize();
  // set all values to random values.
  void makeRandom();
  // set all values to uniform values.
  void makeUniform();

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
  void emIncrement(logpr prob,const int val);
  void emEndIteration();
  void emSwapCurAndNew();
  void emStoreAccumulators(oDataStreamFile& ofile);
  void emLoadAccumulators(iDataStreamFile& ifile);
  void emAccumulateAccumulators(iDataStreamFile& ifile);
  //////////////////////////////////


};



#endif 
