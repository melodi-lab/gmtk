/*-
 * GMTK_Dense1DPMF.h
 *      .h file for GMTK_Sparse1DPMF.h, trainable 1D discrete probability
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


#ifndef GMTK_DISCRETE1DPDF
#define GMTK_DISCRETE1DPDF


#include "fileParser.h"
#include "logp.h"
#include "sArray.h"

class Sparse1DPMF : public EMable {

  //////////////////////////////////////
  // The cardinality of this RV, i.e., 
  // values may take value between [0:card-1]
  int _card;

  struct Entry {
    int val;
    logpr prob;
  };

  ///////////////////////////////////////////////////////////  
  // The probability mass function
  sArray <Entry> pmf
  ///////////////////////////////////////////////////////////  

  //////////////////////////////////
  // Data structures support for EM.
  //   The previous probability mass function 
  //   (we don't use 'Entrys' here as the vals are the same.
  sArray <logpr> nextPmf;
  ///////////////////////////////////////////////////////////  

public:

  ///////////////////////////////////////////////////////////  
  // General constructor
  Sparse1DPMF();

  //////////////////////////////////////
  // Return the number of valid values
  int length() { return pmf.len(); }

  ///////////////////////////////////////
  // return the cardinality of the RV
  int card() { return _card; }

  ///////////////////////////////////////
  // return the probability that for value 'val'
  logpr prob(const int val);

  ///////////////////////////////////////////////////////////  
  // Re-normalize the distribution
  normalize();

  ///////////////////////////////////////////////////////////    
  // read in the basic parameters, assuming file pointer 
  // is located at the correct position.
  void read(iDataStreamFile& is);

  ///////////////////////////////////////////////////////////    
  // write out the basic parameters, starting at the current
  // file position.
  void write(oDataStreamFile& os);


};



#endif // defined DISCRETE1DPDF
