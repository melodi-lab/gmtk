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
#include "GMTK_DiscRV.h"
#include "GMTK_Dense1DPMF.h"


class Sparse1DPMF : public EMable {


  //////////////////////////////////////
  // The cardinality of this RV, i.e., 
  // values may take value between [0:card-1]
  unsigned _card;

  ///////////////////////////////////////////
  // the (possibly) shared 1DPMFs used for the prob values
  Dense1DPMF* dense1DPMF;

  ///////////////////////////////////////////////////////////  
  // The probability mass function
  sArray < DiscRVType > pmf;
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
  unsigned length() { return pmf.size(); }

  ///////////////////////////////////////
  // return the cardinality of the RV
  unsigned card() { return _card; }

  ///////////////////////////////////////
  // return the probability that for value 'val'
  // This will do a search for RV value val
  logpr prob(const unsigned val);

  ///////////////////////////////////////
  // return the probability at entry i
  logpr probAtEntry(const unsigned i) {
    assert ( i < pmf.size() );
    return dense1DPMF->p(i);
  }
  ///////////////////////////////////////
  // return the value at entry i
  DiscRVType valueAtEntry(const unsigned i) {
    assert ( i < pmf.size() );
    return pmf.ptr[i];
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
  // presumably, get parameters from dense pmfs.
  unsigned totalNumberParameters() { return dense1DPMF->totalNumberParameters();  }

  ///////////////////////////////////////////////////////////    
  // write out the basic parameters, starting at the current
  // file position.
  void write(oDataStreamFile& os);


  //////////////////////////////////
  // Public interface support for EM
  //////////////////////////////////
  void emStartIteration();
  void emIncrement(logpr prob,const unsigned val);
  void emEndIteration();
  void emSwapCurAndNew();

  // parallel training
  void emStoreObjectsAccumulators(oDataStreamFile& ofile,
				  bool writeLogVals = true,
				  bool writeZeros = false) {};
  void emLoadObjectsDummyAccumulators(iDataStreamFile& ifile) {};
  void emZeroOutObjectsAccumulators() {};
  void emLoadObjectsAccumulators(iDataStreamFile& ifile) {};
  void emAccumulateObjectsAccumulators(iDataStreamFile& ifile) {};
  const string typeName() { return "SPMF"; }
  //////////////////////////////////


};



#endif 
