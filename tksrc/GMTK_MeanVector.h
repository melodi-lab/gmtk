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

#include "GMTK_RealArray.h"

class MeanVector : public EMable {

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

  //////////////////////////////////
  // set all current parameters to random values
  void makeRandom();

  //////////////////////////////////////////////
  // read/write basic parameters
  void read(iDataStreamFile& is) { means.read(is); }
  void write(oDataStreamFile& os) { means.write(os); }

  //////////////////////////////////
  // Public interface support for EM
  //////////////////////////////////
  void emInit();
  void startEmEpoch();
  void emAccumulate(const float prob,
		    const float *const oo_array);
  void endEmEpoch(logpr cmpSop_acc);
  void emLoadAccumulators(iDataStreamFile& ifile);
  void emStoreAccumulators(oDataStreamFile& ofile);
  void emAccumulateAccumulators(iDataStreamFile& ifile);
  //////////////////////////////////



};



#endif // defined REALARRAY
