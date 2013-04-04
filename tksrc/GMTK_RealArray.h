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
 * Copyright (C) 2001 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 *
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

