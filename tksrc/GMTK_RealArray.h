/*-
 * GMTK_RealArray.h
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


#ifndef GMTK_REALARRAY
#define GMTK_REALARRAY

#include "fileParser.h"
#include "logp.h"
#include "sArray.h"

class RealArray {

  ///////////////////////////////////////////////////////////  
  // The data values
  sArray<float> arr;
  ///////////////////////////////////////////////////////////  

  ///////////////////////////////////////////////////////////  
  // length of the distribution (separate from the memory
  // allocated for arr above).
  int length;
  ///////////////////////////////////////////////////////////  

public:

  ///////////////////////////////////////////////////////////  
  // General constructor
  RealArray() 


  ///////////////////////////////////////////////////////////    
  void read(iDataStreamFile& is);

  ///////////////////////////////////////////////////////////    
  void write(oDataStreamFile& os);

};



#endif // defined REALARRAY
