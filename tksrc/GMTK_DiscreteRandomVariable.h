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
  CPT* curCPT;

  // iterator used between clamp functions.
  CPT::iterator it;

  // Cached value of findConditionalParents(). Can reuse
  // this value w/o needing to do the integer map lookup
  // again.
  unsigned cachedIntFromSwitchingState;

public:

  DiscreteRandomVariable(char * _label, vartype vt, int card);
/*
  DiscreteRandomVariable(char * _label, vartype vt, int card)
    : RandomVariable(_label, vt, card) {;}
*/

  ////////////////////////////////////////////////
  // Assuming the parents have been allocated, this forces
  // the internal CPT structures to 1) be allocated and 
  // 2) match the cardinalities of the parents.
  void allocateProbabiltyTables();


  ////////////////////////////////////////////////////////////////
  // Set up conditional parents pointers and other tables.
  void findConditionalParents();

  ////////////////////////////////////////////////////////////////
  // These next several routines use the switching state and
  // parent set determined the last time findConditionalParents
  // was called.
  // 
  // compute the probability
  logpr probGivenParents() {
    return curCPT->probGivenParents(*curConditionalParents,val);
  }
  // clamp this RV to its "first" value
  void clampFirstValue() { it = curCPT->begin(); val = it.val(); }
  // continue on
  bool clampNextValue() { it++; return (it != curCPT->end()); }
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  // Value caching support.
  DISCRETE_VARIABLE_TYPE cached_val;
  void cacheValue() {cached_val=val;}
  void restoreCachedValue() {val=cached_val;}

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

  ////////////////////////////////////////////////////////////////
  // tie self set of parameters with those of rv
  void tieWith(RandomVariable *rv);


  ///////////////////////////////////////////////////
  // EM Support
  void zeroAccumulators();
  void increment(logpr posterior);
  void update();
  ///////////////////////////////////////////////////


};

#endif
