/*-
 * GMTK_Discrete1DPDF.h
 *      .h file for GMTK_Discrete1DPDF.h, trainable 1D discrete probability
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

class Discrete1DPDF : public EMable {


  ///////////////////////////////////////////////////////////  
  // The probability mass function
  sArray <logpr> pmf
  ///////////////////////////////////////////////////////////  

  //////////////////////////////////
  // Data structures support for EM.
  //   The previous probability mass function 
  sArray <logpr> nextPmf;
  ///////////////////////////////////////////////////////////  

public:

  ///////////////////////////////////////////////////////////  
  // General constructor
  Discrete1DPDF();

  int length() { return pmf.len(); }

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
