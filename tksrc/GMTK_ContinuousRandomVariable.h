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
#include "GMTK_MixGaussiansCommon.h"
#include "GMTK_GMParms.h"
#include "GMTK_NameCollection.h"

#include "GMTK_ObservationMatrix.h"

class ContinuousRandomVariable : public RandomVariable
{
private:
  friend class FileParser;

  /////////////////////////////////////////////
  // A MappingOrDirect object is used to store either
  // a direct pointer to a Gaussian Mixture object
  // (which occurrs when there are no ('nil') conditional parents)
  // or a pointer to a decision tree which is used to map from
  // the current set of conditional parent values an integer
  // which then indexes into the GM_Parms Gaussian Mixture array
  // to locate a particular gaussian.
  class MappingOrDirect {
  public:
    // include an enum with the crv type.
    MappingOrDirect() {
      direct = false;
      gaussian = NULL;
      mapping.dtMapper = NULL;
      mapping.collection = NULL;
    };
    bool direct;
    union { 
      // if direct, a direct pointer to a mixture Gaussian
      MixGaussiansCommon* gaussian;
      // if not direct, a DT and a MG collection object.
      struct MappingStruct {
	// DT to map from parent values to an integer
	RngDecisionTree* dtMapper;
	// the resulting integer is an ofset in this table
	// which points directly to one of the GaussianMixture objects.
	NameCollection* collection;
      } mapping;
    };
  };

  //////////////////////////////////////////////////////////////////////
  // This array has one element for each set of possible conditional parents 
  // we might have (so size of this array is the number of different
  // possible conditional parents), and that equivalent to 
  // the number of regions carved out in the state space
  // of the switching parents.
  vector <MappingOrDirect> conditionalGaussians;


  ////////////////////////////////////////////////////////////////////
  // the current Gaussian after findConditionalParents() is called.
  // It is "current" in the sence that it valid for the set
  // of parent values that are clamped. If the parent values
  // change, this CPT will no longer be valid until another
  // findConditionalParents() is called.
  MappingOrDirect* curMappingOrDirect;

  // cached probability, not used at the moment
  // bool probIsCached;
  // logpr _cachedProb;

  /////////////////////////////////////////////////////////////////
  // the feature range to which this random variable corresponds
  // in a data file. it corresponds to 
  // [firstFeatureElement:lastFeatureElement] inclusive.
  unsigned firstFeatureElement;
  unsigned lastFeatureElement;


public:

  ContinuousRandomVariable(string _label);

  ////////////////////////////////////////////////////////////////
  // Set up conditional parents pointers and other tables.
  void findConditionalParents();

  // compute the probability
  logpr probGivenParents();

  ////////////////////////////////////////////////
  // clamp this RV to its "first" value,
  // presumably this is an observation and contains only one value.
  void clampFirstValue() {
    // we do not support hidden continuous variables just yet.
    assert ( !hidden );
    if (!globalObservationMatrix.active())
      warning("WARNING: clamping value of observation variable w/o observation matrix");
    findConditionalParents();
  }
  ////////////////////////////////////////////////////////////
  // always the last value, since we do not support hidden
  // continuous variables at this time.
  bool clampNextValue() { return false; }
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


  void makeRandom();
  void makeUniform();

  void tieParametersWith(RandomVariable*const other,
			 bool checkStructure=true);

  ////////////////////////////////////////////////////////////////
  // Sample, set value.
  void instantiate() { error("not implemented"); }

  // what is the dimensionality of this variable?
  // intended for indexing fval() after sampling 
  int dimensionality() 
      {error("dimensionality function not implemented"); return 0;}

  // return a pointer to the actual array of values
  // intended for use after sampling
  // will probably return a pointer into the global observation matrix
  float *fval() {error("fval function not implemented"); return NULL;}

  ///////////////////////////////////////////////////
  // EM Support
  void emIncrement(logpr prob);

  ///////////////////////////////////////////////////


  ///////////////////////////////////////////////////
  // reproduction routines.

  RandomVariable *create() { 
    return new ContinuousRandomVariable(label);
  }
  RandomVariable *clone();
  RandomVariable *cloneWithoutParents();

  ///////////////////////////////////////////////////


};

#endif
