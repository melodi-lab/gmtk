/* 
 * GMTK_DiscreteRandomVariable.h
 * A specific case of the RandomVariable class
 *
 * Written by Geoffrey Zweig <gzweig@us.ibm.com> & Jeff Bilmes <bilmes@ee.washington.edu>
 * 
 * Copyright (c) 2001, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about the suitability of this software 
 * for any purpose. It is provided "as is" without express or implied warranty.
 *
 */ 

#ifndef GMTK_DISCRETERANDOMVARIABLE_H
#define GMTK_DISCRETERANDOMVARIABLE_H

#include <vector>

#include "GMTK_RandomVariable.h"
#include "GMTK_CPT.h"

class DiscreteRandomVariable : public RandomVariable
{
private:

  //////////////////////////////////////////////////////////////////////
  // CPT array, one for each set of possible parents we might
  // have (so size of this array is the number of different
  // possible conditional parents).
  vector < CPT* > conditionalCPTs;

  // the current CPT after findConditionalParents() is called.
  // This variable is NULL if there are no parents of this node.
  CPT* curCPT;

  // iterator used between clamp functions.
  CPT::iterator it;

public:

  DiscreteRandomVariable(string _label, int card);

  ////////////////////////////////////////////////
  // Assuming the parents have been allocated, this forces
  // the internal CPT structures to 1) be allocated and 
  // 2) match the cardinalities of the parents.
  void allocateProbabiltyTables();

  void setCpts(vector<CPT*> &cpts);

  ////////////////////////////////////////////////////////////////
  // Set up conditional parents pointers and other tables.
  void findConditionalParents();

  ////////////////////////////////////////////////////////////////
  // These next several routines use the switching state and
  // parent set determined when clampFirstValue() is called.
  // clampFirstValue() invokes findConditionalParents(), and  
  // calls to clampNextValue() and probGivenParents() rely on being
  // called in the appropriate context. i.e. the values of the switching
  // parents must not have changed since the last call to clampFirstValue().
  // This imposes some constraints -- which do happen to be satisfied by
  // the inference loops.
  // 
  // compute the probability
  logpr probGivenParents() {
    return curCPT->probGivenParents(*curConditionalParents,val);
  }
  // clamp this RV to its "first" value
  void clampFirstValue() 
  { findConditionalParents(); it = curCPT->begin(); val = it.val(); }
  // continue on
  bool clampNextValue() { it++; return (it != curCPT->end()); }
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  // Value caching support.
  DISCRETE_VARIABLE_TYPE cached_val;
  void cacheValue() {cached_val=val;}
  void restoreCachedValue() {val=cached_val;}

    /////////////////////////////////////////////////////////////////////////
    // stores a variable's value elsewhere
    virtual void storeValue(VariableValue &vv) {vv.ival = val;}

    /////////////////////////////////////////////////////////////////////////
    // sets a variables value as specified
    virtual void setValue(VariableValue &vv) {val = vv.ival;}

  void makeRandom() { 
    for (unsigned i=0;i<conditionalCPTs.size();i++) 
      conditionalCPTs[i]->makeRandom();
  }

  void makeUniform()  { 
    for (unsigned i=0;i<conditionalCPTs.size();i++) 
      conditionalCPTs[i]->makeUniform();
  }

  ////////////////////////////////////////////////////////////////
  // Sample, set value.
  void instantiate() { val = curCPT->randomSample(); }

  ///////////////////////////////////////////////////
  // EM Support
  void emStartIteration() {error("zeroAccumulators undefined");}
  void emIncrement(logpr posterior) {error("increment undefined");}
  void emEndIteration() {error("update undefined");}
  ///////////////////////////////////////////////////


};

#endif
