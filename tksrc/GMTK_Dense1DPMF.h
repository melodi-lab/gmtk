/*-
 * GMTK_1D_Dist.h
 *      .h file for GMTK_1D_Dist.cc, trainable 1D discrete probability
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


#ifndef GMTK_1D_DIST
#define GMTK_1D_DIST


#include "fileParser.h"
#include "logp.h"
#include "sArray.h"

class Discrete1DPDF {


  ///////////////////////////////////////////////////////////  
  // The probability mass function
  sArray <logpr> pmf
  ///////////////////////////////////////////////////////////  

  ///////////////////////////////////////////////////////////  
  // length of the distribution (separate from the memory
  // allocated for pmf above).
  int length;
  ///////////////////////////////////////////////////////////  

public:

  ///////////////////////////////////////////////////////////  
  // General constructor
  Discrete1DPDF() 


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



#endif
