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

protected:

  // special constructor for subclasses who want to 
  // use a different type.
  MDCPT(DiscreteImplementaton _cptType) : CPT(_cptType) {};

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
    assert ( val >= 0 && val < card() );
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

  // returns an iterator for the first one that is not zero prob.
  iterator begin() {
    assert ( bitmask & bm_basicAllocated );
    iterator it(this);
    it.internalState = 0;
    it.probVal = *mdcpt_ptr;
    if (it.probVal.essentially_zero()) {
      // go to first entry which is not zero.
      do {
	it.internalState++;
	// We keep the following assertion as we
	// must have that at least one entry is non-zero.
	// The read code of the MDCPT should ensure this
	// as sure all parameter update procedures.
	assert (it.internalState < card());
	it.probVal = mdcpt_ptr[it.internalState];
      } while (it.probVal.essentially_zero());
    }
    it.value = it.internalState;
    return it;
  }

  // returns an iterator for the first one that is not zero prob.
  void begin(iterator& it) {
    assert ( bitmask & bm_basicAllocated );
    it.setCPT(this);
    it.internalState = 0;
    it.probVal = *mdcpt_ptr;
    if (it.probVal.essentially_zero()) {
      // go to first entry which is not zero.
      do {
	it.internalState++;
	// We keep the following assertion as we
	// must have that at least one entry is non-zero.
	// The read code of the MDCPT should ensure this
	// as sure all parameter update procedures.
	assert (it.internalState < card());
      } while (mdcpt_ptr[it.internalState].essentially_zero());
      it.probVal = mdcpt_ptr[it.internalState];
    }
    it.value = it.internalState;
  }

  // Given a current iterator, return the next one in the sequence.
  // Skip the zero probability ones.
  bool next(iterator &it) {
    // don't increment past the last value.
    do {
      it.internalState++;
      if (it.internalState == card())
	return false;
      it.probVal = mdcpt_ptr[it.internalState];
    } while (it.probVal.essentially_zero());
    it.value = it.internalState;
    return true;
  }


  bool end(iterator& it) {
    return (it.internalState == card());
  }

  ///////////////////////////////////
  int randomSample();

  ///////////////////////////////////
  unsigned totalNumberParameters() { return mdcpt.len(); }

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


  // parallel training
  void emStoreObjectsAccumulators(oDataStreamFile& ofile);
  void emLoadObjectsDummyAccumulators(iDataStreamFile& ifile);
  void emZeroOutObjectsAccumulators();
  void emLoadObjectsAccumulators(iDataStreamFile& ifile);
  void emAccumulateObjectsAccumulators(iDataStreamFile& ifile);
  const string typeName() { return "DenseCPT"; }
  //////////////////////////////////


};



#endif // defined MDCPT
