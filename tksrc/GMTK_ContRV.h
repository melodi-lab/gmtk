/*
 * GMTK_ContRV.h
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *  $Header$
 *
 * Copyright (C) 2001 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 *
 *
 *
 *
 *
 *
 */

#ifndef GMTK_CONT_RV_H
#define GMTK_CONT_RV_H

#include <vector>
#include <string>
#include <set>

#include "GMTK_RV.h"
#include "GMTK_CPT.h"
#if 0
#  include "GMTK_ObservationMatrix.h"
#else
#  include "GMTK_ObservationSource.h"
#endif
#include "GMTK_MixtureCommon.h"
#include "GMTK_GMParms.h"
#include "GMTK_NameCollection.h"


class FileParser;

class ContRV : public RV {
  friend class FileParser;

protected:

  /////////////////////////////////////////////
  // Note that a non-switching ContRV must have parents
  struct MappingStruct {
    // DT to map from parent values to an integer
    RngDecisionTree* dtMapper;
    // the resulting integer is an ofset in this table
    // which points directly to one of the Mixture objects.
    NameCollection* collection;
  };

  //////////////////////////////////////////////////////////////////////////
  // A MappingOrDirect object is used to store either a direct pointer
  // to a Mixture object (which occurs when there are 'nil'
  // conditional parents) or a pointer to a decision tree which is
  // used to map from the current set of conditional parent values an
  // integer which then indexes into the current collection array that
  // locates a particular mixture in the global GM_Parms Mixture
  // array. Note that for non-switching ContRVs, direct is almost always
  // false (but you still need to check). In the switching case
  // (subclasses of this class) we need to check the direct field of
  // this object.
  class MappingOrDirect {
  public:
    // include an enum with the crv type.
    MappingOrDirect() {
      direct = false;
      mixture = NULL;
      mapping.dtMapper = NULL;
      mapping.collection = NULL;
    };
    bool direct;
    union { 
      // if direct, a direct pointer to a mixture
      MixtureCommon* mixture;
      // if not direct, a DT and a MG collection object.
      struct MappingStruct mapping;
    };
  };

  //////////////////////////////////////////////////////////////////////
  // This array has one element for each set of possible conditional parents 
  // we might have (so size of this array is the number of different
  // possible conditional parents), and that equivalent to 
  // the number of regions carved out in the state space
  // of the switching parents. For a non-switching object (such
  // as the present object), this array always has length 1. For
  // subclassses with switching functionality , the array might be longer.
  vector <MappingOrDirect> conditionalMixtures;

  ////////////////////////////////////////////////////////////////////
  // the current mixture after setConditionalParents() is called.  It
  // is "current" in the sence that it valid for the current set of
  // parent values. If the parent values change, however, this pointer
  // will no longer be valid until another setConditionalParents() is
  // called.
  MappingOrDirect* curMappingOrDirect;

public:
  
  /////////////////////////////////////////////////////////////////////////
  // constructor: Initialize with the variable type.  The default
  // timeIndex value of ~0x0 indicates a static network.  Discrete
  // nodes must be specified with their cardinalities.
  ContRV(RVInfo& _rv_info,
	 unsigned _timeFrame = ~0x0)
    : RV(_rv_info,_timeFrame)
  {
  }

  virtual ~ContRV() {;}


  // printing routines.
  virtual void printNameFrameValue(FILE *f,bool nl=true) {
    RV::printNameFrame(f,false);
    fprintf(f,"=C%s",nls(nl));
  }

  virtual void printSelf(FILE *f,bool nl=true);
  virtual void printSelfVerbose(FILE *f) = 0;

  // what is the dimensionality of this variable?
  // intended for indexing fval() after sampling 
  int dimensionality() 
  {  return (rv_info.rvFeatureRange.lastFeatureElement - 
	     rv_info.rvFeatureRange.firstFeatureElement + 1); }


  //////////////////////////////////////////////////////////
  // computing with probabilities
  //////////////////////////////////////////////////////////

  // no need to refine the parent class functions.
  // virtual void probGivenParents(logpr& p) = 0;
  // virtual logpr maxValue() = 0;

  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  // Parameter and Structure adjustment/checking       /////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////

  void makeRandom() = 0;
  void makeUniform() = 0;


  //////////////////////////////////////////////////////////
  // Iterate through values of this RV                   ///
  //////////////////////////////////////////////////////////

  void begin(logpr& p) = 0;
  bool next(logpr& p) = 0;
  void randomSample() = 0;


  /////////////////////////////////////////
  // EM Learning                         //
  /////////////////////////////////////////

  void emIncrement(logpr posterior) = 0;

  /////////////////////////////////////////
  // Cloning Support                     //
  /////////////////////////////////////////

  virtual ContRV* cloneRVShell();

};


// Define a few down-cast routine here for when we have a RV that is
// known to be a ContRV. Note that this is not entirely type-safe, but
// speed is more important than type safety here, so we need to be
// able to do static down casts.
inline ContRV* RV2CRV(RV* rv) { 
  return ((ContRV*)rv);
}
inline ContRV& RV2CRV(RV& rv) { 
  return *((ContRV*)&rv);
}


#endif


