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
  RngDecisionTree* dt;

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
    if (_val >= card()) {
      warning("ERROR: Deterministic CPT '%s' of card %d querying DT '%s' received value %d",
	      name().c_str(),
	      card(),
	      dt->name().c_str(),
	      _val);
      fprintf(stderr,"Parent values:");
      for (unsigned i=0;i<parentValues.size();i++) {
	fprintf(stderr," %d", parentValues[i]);
      }
      error("");
    }
  }
  void becomeAwareOfParentValues( vector <RandomVariable *>& parents ) {
    _val = dt->query(parents);
    if (_val >= card()) {
      warning("ERROR: Deterministic CPT '%s' of card %d querying DT '%s' received value %d",
	      name().c_str(),
	      card(),
	      dt->name().c_str(),
	      _val);
      fprintf(stderr,"Parents configuration :");
      for (unsigned i=0;i<parents.size();i++) {
	fprintf(stderr,"%s(%d)=%d,",
		parents[i]->name().c_str(),
		parents[i]->timeIndex,
		parents[i]->val);
      }
      error("");
    }
  }
  logpr probGivenParents(const int val) {
    assert ( bitmask & bm_basicAllocated );
    assert ( val >= 0 && val <= card() );
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

  // returns an iterator for the first one.  It *must* be that
  // becomeAwareOfParentValues has already been called
  iterator begin() {
    assert ( bitmask & bm_basicAllocated );
    iterator it(this);
    it.internalState = 0;
    it.probVal = 1.0;
    it.value = _val;
    return it;
  }


  // returns an iterator for the first one.  It *must* be that
  // becomeAwareOfParentValues has already been called
  void begin(iterator& it) {
    assert ( bitmask & bm_basicAllocated );
    it.setCPT(this);
    // indicates internal state
    it.internalState = 0;
    it.probVal = 1.0;
    it.value = _val;
  }

  inline bool next(iterator &it) {
    // this is an MTCPT so we end here immediately.
    it.internalState = 1;
    return false;
  }

  bool end(iterator& it) {
    return (it.internalState == 1);
  }


  virtual int valueAtIt(const int internalState) { 
    assert ( internalState == 0);
    return _val;
  }

  // used for elimination/triangulation
  virtual unsigned averageCardinality() { return 1; }
  virtual unsigned maxCardinality() { return 1; }


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


  ///////////////////////////////////
  // get parameters from dense pmfs
  unsigned totalNumberParameters() { return 0; }

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


  // parallel training
  void emStoreObjectsAccumulators(oDataStreamFile& ofile) {}
  void emLoadObjectsDummyAccumulators(iDataStreamFile& ifile) {}
  void emZeroOutObjectsAccumulators() {}
  void emLoadObjectsAccumulators(iDataStreamFile& ifile) {}
  void emAccumulateObjectsAccumulators(iDataStreamFile& ifile) {}
  const string typeName() { return "DeterministicCPT"; }
  //////////////////////////////////

};



#endif // defined MTCPT
