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

class Dense1DPMF : public EMable {

  ///////////////////////////////////////////////////////////  
  // the name
  char *_name;

  ///////////////////////////////////////////////////////////  
  // The probability mass function
  sArray <logpr> pmf;
  ///////////////////////////////////////////////////////////  

  //////////////////////////////////
  // Data structures support for EM.
  //   The previous probability mass function 
  sArray <logpr> nextPmf;
  ///////////////////////////////////////////////////////////  

public:

  ///////////////////////////////////////////////////////////  
  // General constructor
  Dense1DPMF();
  ~Dense1DPMF() { delete [] _name; }

  int length() { return pmf.len(); }

  int card() { return pmf.len(); }

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
