/*-
 * GMTK_MSCPT.h
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


#ifndef GMTK_MSCPT_H
#define GMTK_MSCPT_H

#include <vector>

#include "fileParser.h"
#include "logp.h"

#include "GMTK_RngDecisionTree.h"
#include "GMTK_RandomVariable.h"
#include "GMTK_CPT.h"
#include "GMTK_Sparse1DPMF.h"
#include "GMTK_NameCollection.h"

#include "GMTK_EMable.h"
#include "GMTK_GMParms.h"
#include "GMTK_NamedObject.h"
#include "GMTK_NameCollection.h"

class MSCPT : public CPT {

  //////////////////////////////////
  // Index into the world structure
  // of the decision tree
  unsigned dtIndex; 

  ///////////////////////////////////////
  // Direct pointer to the decision tree.
  RngDecisionTree* dt;

  ///////////////////////////////////////
  // Direct pointer to collection (indirect
  // mapping to spmfs).
  NameCollection* ncl;

  ///////////////////////////////////////
  // mapping from DT leaves to SPMFs
  NameCollection* spmfCollection;

  ///////////////////////////////////////
  // Index of world's sparse mass function,
  // cached for current value of parents.
  unsigned spmfIndex;

  ///////////////////////////////////////
  // Direct index to the Sparse PMF
  Sparse1DPMF* spmf;



public:

  ///////////////////////////////////////////////////////////  
  // General constructor
  MSCPT();
  ~MSCPT() { }

  ///////////////////////////////////////////////////////////    
  void setNumParents(const int _nParents);
  void setNumCardinality(const int var, const int card);
  void allocateBasicInternalStructures();
  ///////////////////////////////////////////////////////////    

  //////////////////////////////////
  // various forms of probability calculation
  void becomeAwareOfParentValues( vector <int>& parentValues,
				  vector <int>& cards ) {
    spmfIndex = dt->query(parentValues,cards);
    if (spmfIndex < 0 || spmfIndex >=  GM_Parms.sPmfs.size()) {
      error("ERROR: MSCPT '%s' uses DT '%s' with invalid SPMF index '%d'\n",
	    name().c_str(),dt->name().c_str(),spmfIndex);
    }
    spmf = GM_Parms.sPmfs[spmfIndex];
    if (spmf->card() != card()) {
      warning("ERROR: MSCPT '%s' of card %d querying DT '%s' received index %d of SPMF '%s' having card %d",
	      name().c_str(),
	      card(),
	      dt->name().c_str(),
	      spmfIndex,
	      spmf->name().c_str(),
	      spmf->card());
      fprintf(stderr,"Parent values:");
      for (unsigned i=0;i<parentValues.size();i++) {
	fprintf(stderr," %d", parentValues[i]);
      }
      error("");
    }
  }
  void becomeAwareOfParentValues( vector <RandomVariable *>& parents ) {
    spmfIndex = dt->query(parents);
    if (spmfIndex < 0 || spmfIndex >=  GM_Parms.sPmfs.size()) {
      error("ERROR: MSCPT '%s' uses DT '%s' with invalid SPMF index '%d'\n",
	    name().c_str(),dt->name().c_str(),spmfIndex);
    }
    spmf = GM_Parms.sPmfs[spmfIndex];
    if (spmf->card() != card()) {
      warning("ERROR: MSCPT '%s' of card %d querying DT '%s' received index %d of SPMF '%s' having card %d",
	      name().c_str(),
	      card(),
	      dt->name().c_str(),
	      spmfIndex,
	      spmf->name().c_str(),
	      spmf->card());
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
    return spmf->prob(val);
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

  // returns an iterator for the first one.
  iterator begin() {
    assert ( bitmask & bm_basicAllocated );
    iterator it(this);
    it.internalState = 0;
    it.probVal = spmf->probAtEntry(0);
    return it;
  }

  iterator end() {
    assert ( bitmask & bm_basicAllocated );
    iterator it(this);
    it.internalState = spmf->length();
    return it;
  }
  bool next(iterator &it) {
    assert ( bitmask & bm_basicAllocated );
    it.internalState++;
    // don't increment past the last value.
    if (it.internalState >= spmf->length()) {
      it.internalState = spmf->length();
      return false;
    }
    it.probVal = spmf->probAtEntry(it.internalState);
    return true;
  }
  virtual int valueAtIt(const int internalState) { 
    assert ( internalState >= 0 && internalState < spmf->length() );
    return spmf->valueAtEntry(internalState); 
  }

  int randomSample();

  ///////////////////////////////////////////////////////////  
  // Re-normalize the output distributions
  void normalize();
  // set all values to random values.
  void makeRandom();
  // set all values to uniform values.
  void makeUniform();

  ///////////////////////////////////////////////////////////
  // get parameters from dense pmfs, since this heavily
  // depends on the DT mapping.
  unsigned totalNumberParameters() { return 0;  }

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
  void emStoreZeroAccumulators(oDataStreamFile& ofile);
  void emLoadAccumulators(iDataStreamFile& ifile);
  void emAccumulateAccumulators(iDataStreamFile& ifile);
  //////////////////////////////////

};



#endif // defined MSCPT
