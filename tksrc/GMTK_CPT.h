/*-
 * GMTK_CPT
 *      .h file the .cc file.
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


#ifndef GMTK_CPT
#define GMTK_CPT


#include "fileParser.h"
#include "logp.h"
#include "sArray.h"

class CPT {


  ///////////////////////////////////////////////////////////  
  // The 'rows' probability mass functions
  sArray <logpr> pmf
  ///////////////////////////////////////////////////////////  

  ///////////////////////////////////////////////////////////  
  // cardinality of parents
  int rows
  ///////////////////////////////////////////////////////////  

  ///////////////////////////////////////////////////////////  
  // cardinality of self
  int cols;
  ///////////////////////////////////////////////////////////  

public:

  ///////////////////////////////////////////////////////////  
  // General constructor
  CPT() 


  ///////////////////////////////////////////////////////////  
  // Re-normalize the distributions
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


#endif // defined GMTK_CPT
