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

#ifndef GMTK_DISCRETERANDOMVARIABLE
#define GMTK_DISCRETERANDOMVARIABLE

#include "GMTK_RandomVariable.h"
#include "GMTK_CPT.h"

class DiscreteRandomVariable : public RandomVariable
{
private:

  //////////////////////////////////////////////////////////////////////
  // CPT array, one for each set of possible parents we might
  // have (so size of this array is the number of different
  // possible conditional parents).
  sArray < CPT* > conditionalCPTs;

  // the current CPT after findConditionalParents() is called.
  CPT* curCPT;

  // iterator used between clamp functions.
  CPT::iterator it;

  // Cached value of findConditionalParents(). Can reuse
  // this value w/o needing to do the integer map lookup
  // again.
  int cachedIntFromSwitchingState;

public:

  DiscreteRandomVariable(string _label, vartype vt, int card=0)
  {RandomVariable::RandomVariable(string _label, vartype vt, int card=0);}

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
    return curCPT->probGivenParents(curConditionalParents,val);
  }
  // clamp this RV to its "first" value
  void clampFirstValue() { it = curCPT->first(); val = it.val; }
  // continue on
  bool clampNextValue() { curCPT->next(it); val = it.val; }
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  // Value caching support.
  DISCRETE_VARIABLE_TYPE cached_val;
  void cacheValue() {cached_val=val;}
  void restoreCachedValue() {val=cached_val;}

  void makeRandom() { 
    for (int i=0;i<conditionalCPTs.len();i++) 
      conditionalCPTs[i].makeRandom();
  }

  void makeUniform()  { 
    for (int i=0;i<conditionalCPTs.len();i++) 
      conditionalCPTs[i].makeUniform();
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
