/* 
 * GMTK_ContinuousRandomVariable.h
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

#ifndef GMTK_CONTINUOUSRANDOMVARIABLE_H
#define GMTK_CONTINUOUSRANDOMVARIABLE_H

#include <vector>

#include "GMTK_RandomVariable.h"
#include "GMTK_CPT.h"


class ContinuousRandomVariable : public RandomVariable
{
private:


public:

  ContinuousRandomVariable(string _label) 
    : RandomVariable(_label,Continuous) {}

  ////////////////////////////////////////////////////////////////
  // Set up conditional parents pointers and other tables.
  void findConditionalParents();

  // compute the probability
  logpr probGivenParents();

  ////////////////////////////////////////////////
  // clamp this RV to its "first" value,
  // presumably this is an observation and contains only one value.
  void clampFirstValue();
  { findConditionalParents(); it = curCPT->begin(); val = it.val(); }
  // always the last value.
  bool clampNextValue() { return true; }
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  // Value caching support.
  CONTINUOUS_VARIABLE_TYPE cached_val;
  void cacheValue() {cached_val=val;}
  void restoreCachedValue() {val=cached_val;}

  /////////////////////////////////////////////////////////////////////////
  // stores a variable's value elsewhere
  virtual void storeValue(VariableValue &vv) {vv.ival = val;}

  /////////////////////////////////////////////////////////////////////////
  // sets a variables value as specified
  virtual void setValue(VariableValue &vv) {val = vv.ival;}

  void makeRandom();
  void makeUniform();

  ////////////////////////////////////////////////////////////////
  // Sample, set value.
  void instantiate();

  ///////////////////////////////////////////////////
  // EM Support
  void emStartIteration();
  void emIncrement(logpr posterior);
  void emEndIteration();
  ///////////////////////////////////////////////////


};

#endif
