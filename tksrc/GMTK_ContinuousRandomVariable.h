/* 
 * GMTK_ContinuousRandomVariable.h
 * A specific case of the RandomVariable class
 *
 * Jeff Bilmes <bilmes@ee.washington.edu>
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

#include "logp.h"

#include "GMTK_RandomVariable.h"
#include "GMTK_CPT.h"
#include "GMTK_FileParser.h"
#include "GMTK_MixGaussians.h"
#include "GMTK_GMParms.h"


class ContinuousRandomVariable : public RandomVariable
{
private:

  // include an enum with the crv type.
  class MappingOrDirect {
  public:
    MappingOrDirect() {
      direct = false;
      gaussian = NULL;
      dtMapper = NULL;
    };
    bool direct;
    union { 
      MixGaussians* gaussian;
      RngDecisionTree<unsigned>* dtMapper;
    };
  };

  //////////////////////////////////////////////////////////////////////
  // array, one for each set of possible conditional parents we might
  // have (so size of this array is the number of different
  // possible conditional parents), and that is determined
  // by the number of regions carved out in the state space
  // of the switching parents.
  vector <MappingOrDirect> conditionalGaussians;

  // the current Gaussian after findConditionalParents() is called.
  MappingOrDirect* curMappingOrDirect;

  //////////////////////////////////////////////////////////////
  // the feature range to which this random variable corresponds.
  // it corresponds to [firstFeature:lastFeature] inclusive
  unsigned firstFeature;
  unsigned lastFeature;


public:

  ContinuousRandomVariable(string _label);

  ////////////////////////////////////////////////////////////////
  // Set up conditional parents pointers and other tables.
  void findConditionalParents();

  // compute the probability
  logpr probGivenParents() { return 0.0; }

  ////////////////////////////////////////////////
  // clamp this RV to its "first" value,
  // presumably this is an observation and contains only one value.
  void clampFirstValue() {}
  // always the last value.
  bool clampNextValue() { return true; }
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  // Value caching support.
  void cacheValue() { error("not implemented"); }
  void restoreCachedValue() { error("not implemented"); }
  /////////////////////////////////////////////////////////////////////////
  // stores a variable's value elsewhere
  void storeValue(VariableValue &vv) { error("not implemented"); }
  /////////////////////////////////////////////////////////////////////////
  // sets a variables value as specified
  void setValue(VariableValue &vv) { error("not implemented"); }


  void makeRandom() {}
  void makeUniform() {}

  void tieParametersWith(RandomVariable*const other) { error("not implemented"); }

  ////////////////////////////////////////////////////////////////
  // Sample, set value.
  void instantiate() {}

  ///////////////////////////////////////////////////
  // EM Support
  void emStartIteration() {}
  void emIncrement(logpr posterior) {}
  void emEndIteration() {}
  void emClearAllocatedBit() { } 
  void emClearSwappedBit() { }
  void emSwapCurAndNew() { }
  ///////////////////////////////////////////////////


  ///////////////////////////////////////////////////
  // reproduction routines.

  RandomVariable *create() { 
    return new ContinuousRandomVariable(label);
  }
  RandomVariable *clone() { error("not implemented"); return this; }

  ///////////////////////////////////////////////////


};

#endif
