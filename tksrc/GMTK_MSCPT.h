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

#include "fileParser.h"
#include "logp.h"
#include "sArray.h"

#include "GMTK_CPT.h"
#include "GMTK_EMable.h"
#include "GMTK_RandomVariable.h"
#include "GMTK_RngDecisionTree.h"

class MSCPT : public EMable, public CPT {

  //////////////////////////////////

  //////////////////////////////////
  // Index into the world structure
  // of the decision tree
  int dtIndex; 

  ///////////////////////////////////////
  // Direct pointer to the decision tree.
  RngDecisionTree* dt;

  ///////////////////////////////////////
  // Index of world's sparse mass function,
  // cached for current value of parents.
  int spmfIndex;

  ///////////////////////////////////////
  // Direct index to the Sparse PMF
  Sparse1DPMF* spmf;

public:

  ///////////////////////////////////////////////////////////  
  // General constructor
  MSCPT();

  ///////////////////////////////////////////////////////////    
  void setNumParents(const int _nParents) 
     { error("Not Implemented"); }
  void setNumCardinality(const int var, const int card)
     { error("Not Implemented"); }
  void allocateBasicInternalStructures()
     { error("Not Implemented"); }
  ///////////////////////////////////////////////////////////    

  //////////////////////////////////
  // various forms of probability calculation
  void becomeAwareOfParentValues( sArray <int>& parentValues ) {
    dtIndex = dt->query(parentValues);
    spmf = GM_Parms.dts[dtIndex];
  }
  void becomeAwareOfParentValues( sArray <randomVariable *>& parents ) {
    dtIndex = dt->query(parents);
    spmf = GM_Parms.dts[dtIndex];
  }
  logpr probGivenParents(const int val) {
    assert ( bitmask & bm_basicAllocated );
    assert ( val >= 0 && val <= cardinalities[numParents] );
    return spmf->prob(val);
  }
  logpr probGivenParents(sArray <int>& parentValues, 
			 const int val) {
    assert ( bitmask & bm_basicAllocated );
    becomeAwareOfParentValues(parentValues);
    return probGivenParents(val);
  }
  logpr probGivenParents(sArray <randomVariable *>& parents,
			 const int val) {
    assert ( bitmask & bm_basicAllocated );
    becomeAwareOfParentValues(parents);
    return probGivenParents(val);
  }
  int numValsGivenParents() { 
    assert ( bitmask & bm_basicAllocated );
    return cardinalities[numParents]; 
  }

  // returns an iterator for the first one.
  iterator first() {
    assert ( bitmask & bm_basicAllocated );
    iterator it;
    it.val = 0;
    it.probVal = *mscpt_ptr;
    return it;
  }

  // Given a current iterator, return the next one in the sequence.
  bool next(iterator &it) {
    assert ( bitmask & bm_basicAllocated );

    if (it.val == cardinalities[numParents]-1)
      return false;
    it.val++;
    it.probVal = mscpt_ptr[it.val];
    return true;
  }
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
  void emInit() { error("not implemented"); }
  void startEmEpoch()  { error("not implemented"); }
  void emAccumulate(const float prob,
		    const float *const oo_array) { error("not implemented"); }
  void endEmEpoch(logpr cmpSop_acc)  { error("not implemented"); }
  void emLoadAccumulators(iDataStreamFile& ifile)  { error("not implemented"); }
  void emStoreAccumulators(oDataStreamFile& ofile)  { error("not implemented"); }
  void emAccumulateAccumulators(iDataStreamFile& ifile)  { error("not implemented"); }
  void swapCurAndNew() { error("not implemented"); }
  //////////////////////////////////



};



#endif // defined MSCPT
