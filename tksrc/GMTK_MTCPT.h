/*-
 * GMTK_MTCPT.h
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


#ifndef GMTK_MTCPT_H
#define GMTK_MTCPT_H

#include <vector>

#include "fileParser.h"
#include "logp.h"

#include "GMTK_RngDecisionTree.h"
#include "GMTK_RandomVariable.h"
#include "GMTK_CPT.h"
#include "GMTK_Sparse1DPMF.h"

#include "GMTK_EMable.h"
#include "GMTK_GMParms.h"
#include "GMTK_NamedObject.h"


class MTCPT : public CPT  {

  //////////////////////////////////
  // Index into the world structure
  // of the decision tree
  unsigned dtIndex; 

  ///////////////////////////////////////
  // Direct pointer to the decision tree.
  RngDecisionTree<unsigned>* dt;

  ////////////////
  // the value
  int _val;

public:

  ///////////////////////////////////////////////////////////  
  // General constructor
  MTCPT();
  ~MTCPT() { }

  ///////////////////////////////////////////////////////////    
  void setNumParents(const int _nParents);
  void setNumCardinality(const int var, const int card);
  void allocateBasicInternalStructures();
  ///////////////////////////////////////////////////////////    

  //////////////////////////////////
  // various forms of probability calculation
  void becomeAwareOfParentValues( vector <int>& parentValues,
				  vector <int>& cards ) {
    _val = dt->query(parentValues,cards);
  }
  void becomeAwareOfParentValues( vector <RandomVariable *>& parents ) {
    _val = dt->query(parents);
  }
  logpr probGivenParents(const int val) {
    assert ( bitmask & bm_basicAllocated );
    assert ( val >= 0 && val <= cardinalities[_numParents] );
    if (val == _val)
      return 1.0;
    else
      return 0.0;
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
    iterator it;
    it.internalState = 0;
    it.probVal = 1.0;
    return it;
  }

  iterator end() {
    assert ( bitmask & bm_basicAllocated );
    iterator it(this);
    it.internalState = 1;
    return it;
  }
  bool next(iterator &it) {
    assert ( bitmask & bm_basicAllocated );
    // don't increment past the last value.
    if (it.internalState == 1)
      return false;
    it.internalState = 1;
    it.probVal = 0.0;
    return true;
  }
  virtual int valueAtIt(const int internalState) { 
    assert ( internalState == 0);
    return _val;
  }

  // random sample given current parents value
  int randomSample() { return _val; }
  
  //////////
  // these routines do nothing for an MTCPT since
  // there is only one possible value for a DT given the parents.
  ///////////////////////////////////////////////////////////  
  // Re-normalize the output distributions
  void normalize() {}
  // set all values to random values.
  void makeRandom() {}
  // set all values to uniform values.
  void makeUniform() {}

  //////////////////////////////////////////////
  // read/write basic parameters
  void read(iDataStreamFile& is);
  void write(oDataStreamFile& os);

  //////////////////////////////////
  // Public interface support for EM
  //////////////////////////////////
  void emStartIteration();
  void emIncrement(logpr p,RandomVariable*);
  void emEndIteration();
  void emSwapCurAndNew();
  void emStoreAccumulators(oDataStreamFile& ofile);
  void emLoadAccumulators(iDataStreamFile& ifile);
  void emAccumulateAccumulators(iDataStreamFile& ifile);
  //////////////////////////////////

};



#endif // defined MTCPT
