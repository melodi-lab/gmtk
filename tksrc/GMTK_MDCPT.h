/*-
 * GMTK_MDCPT.h
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


#ifndef GMTK_MDCPT_H
#define GMTK_MDCPT_H

#include "fileParser.h"
#include "logp.h"
#include "sArray.h"
#include "CPT.h"

class MDCPT : public CPT {

  //////////////////////////////////
  // The acutal cpt. This is the table for
  // Pr( variable at mdcpt.len()-1 |  variables from 0 to mdcpt.len()-1 )
  sArray < logpr > mdcpt;

  //////////////////////////////////
  // Support for computing probabilities as below.
  // This indices the section of the table mdcpt above
  // for a particular set of parent values. 
  logpr* mdcpt_ptr; 
  // The accumulative cardinalities to help index into
  // the table above.
  sArray <int> cumulativeCardinalities;

  //////////////////////////////////
  // Data structures support for EM
  sArray < logpr > nextMdcpt;


public:

  ///////////////////////////////////////////////////////////  
  // General constructor
  MDCPT();

  //////////////////////////////////
  // set all current parameters to random values
  void randomize();


  //////////////////////////////////
  // various forms of probability calculation
  void setParentValues( sArray <int>& parentValues );
  logpr probGivenParents(const int val);
  logpr probGivenParents(sArray <int>& parentValues, 
				 const int val);
  int numValsGivenParents();
  // returns an iterator for the first one.
  iterator first();
  // Given a current iterator, return the next one in the sequence.
  iterator next(iterator &);



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



#endif // defined MDCPT
