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


#ifndef GMTK_MEANVECTOR_H
#define GMTK_MEANVECTOR_H

#include "fileParser.h"
#include "logp.h"
#include "sArray.h"

#include "GMTK_RealArray.h"
#include "GMTK_EMable.h"
#include "GMTK_NamedObject.h"

class MeanVector : public EMable, public NamedObject {


  //////////////////////////////////
  // The acutal mean vector
  RealArray means;

  //////////////////////////////////
  // Data structures support for EM
  RealArray nextMeans;


public:

  ///////////////////////////////////////////////////////////  
  // General constructor
  MeanVector();
  ~MeanVector() { } 

  //////////////////////////////////
  // set all current parameters to random values
  void makeRandom();

  //////////////////////////////////////////////
  // read/write basic parameters
  void read(iDataStreamFile& is) { 
    NamedObject::read(is);
    means.read(is); 
  }
  void write(oDataStreamFile& os) { 
    NamedObject::write(os);
    means.write(os); 
  }

  //////////////////////////////////
  // Public interface support for EM
  //////////////////////////////////
  void emInit() {}
  void emStartIteration() {}
  void emIncrement(logpr prob) {}
  void emEndIteration() {}
  void emSwapCurAndNew() {}
  void emStoreAccumulators(oDataStreamFile& ofile) {}
  void emLoadAccumulators(iDataStreamFile& ifile) {}
  void emAccumulateAccumulators(iDataStreamFile& ifile) {}
  //////////////////////////////////



};



#endif // defined MEANVECTOR_H
