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

#include "GMTK_CPT.h"
#include "GMTK_EMable.h"
#include "GMTK_RandomVariable.h"
#include "GMTK_NamedObject.h"

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
  vector <int> cumulativeCardinalities;

  //////////////////////////////////
  // Data structures support for EM
  sArray < logpr > nextMdcpt;

  // the overall expected occurence of this CPT
  logpr accumulator;


public:

  ///////////////////////////////////////////////////////////  
  // General constructor
  MDCPT();
  ~MDCPT() { }

  ///////////////////////////////////////////////////////////    
  // Semi-constructors: useful for debugging.
  // See parent class for further documention.
  void setNumParents(const int _nParents);
  void setNumCardinality(const int var, const int card);
  void allocateBasicInternalStructures();

  //////////////////////////////////
  // various forms of probability calculation
  void becomeAwareOfParentValues( vector <int>& parentValues,
				  vector <int>& cards);
  void becomeAwareOfParentValues( vector <RandomVariable *>& parents );

  logpr probGivenParents(const int val) {
    assert ( bitmask & bm_basicAllocated );
    assert ( val >= 0 && val < cardinalities[_numParents] );
    return *(mdcpt_ptr + val);
  }
  logpr probGivenParents(vector <int>& parentValues, 
			 vector <int>& cards, 
			 const int val) {
    assert ( bitmask & bm_basicAllocated );
    becomeAwareOfParentValues(parentValues,cards);
    return probGivenParents(val);
  }
  logpr probGivenParents(vector <RandomVariable *>& parents,
			 const int val) {
    assert ( bitmask & bm_basicAllocated );
    becomeAwareOfParentValues(parents);
    return probGivenParents(val);
  }
  int numValsGivenParents() { 
    assert ( bitmask & bm_basicAllocated );
    return cardinalities[_numParents]; 
  }


  // returns an iterator for the first one.
  iterator begin() {
    assert ( bitmask & bm_basicAllocated );
    iterator it(this);
    it.internalState = 0;
    it.probVal = *mdcpt_ptr;
    return it;
  }

  iterator end() {
    assert ( bitmask & bm_basicAllocated );
    iterator it(this);
    it.internalState = cardinalities[_numParents];
    return it;
  }

  // Given a current iterator, return the next one in the sequence.
  bool next(iterator &it) {
    assert ( bitmask & bm_basicAllocated );
    // don't increment past the last value.
    it.internalState++;
    if (it.internalState == cardinalities[_numParents])
      return false;
    it.probVal = mdcpt_ptr[it.internalState];
    return true;
  }


  //////////////////
  int randomSample();

  
  ///////////////////////////////////////////////////////////  
  // Re-normalize the output distributions
  void normalize();
  // set all values to random values.
  void makeRandom();
  // set all values to uniform values.
  void makeUniform();

  //////////////////////////////////////////////
  // read/write basic parameters
  void read(iDataStreamFile& is);
  void write(oDataStreamFile& os);

  //////////////////////////////////
  // Public interface support for EM
  //////////////////////////////////
  void emStartIteration();
  void emIncrement(logpr p,RandomVariable*rv);
  void emEndIteration();
  void emSwapCurAndNew();
  void emStoreAccumulators(oDataStreamFile& ofile);
  void emLoadAccumulators(iDataStreamFile& ifile);
  void emAccumulateAccumulators(iDataStreamFile& ifile);
  //////////////////////////////////



};



#endif // defined MDCPT
