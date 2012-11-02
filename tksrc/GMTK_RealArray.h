/*-
 * GMTK_RealArray.h
 *    An array of real values with some functionality for 
 *         reading, writing, and so on. This is just an sArray<float> packaged
 *         with  reading and writing routines.
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


#ifndef GMTK_REALARRAY_H
#define GMTK_REALARRAY_H

#include "fileParser.h"
#include "logp.h"
#include "sArray.h"

class RealArray : public sArray<float> {



public:

  ///////////////////////////////////////////////////////////  
  // General constructor
  RealArray();

  ///////////////////////////////////////////////////////////    
  // reading and writing from/to disk.
  void read(iDataStreamFile& is);
  void write(oDataStreamFile& os);

};



#endif // defined REALARRAY_H

