/*
 * GMTK_Sw_HidDiscRV.h
 *
 *  Observed Discrete Random Variables with switching parent functionality.
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *  $Header$
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
 *
 * The discrete random variable type.
 *
 *
 *
 */

#ifndef GMTK_SW_HID_DISC_RV_H
#define GMTK_SW_HID_DISC_RV_H

#include <vector>

#include "GMTK_HidDiscRV.h"
#include "GMTK_SwDiscRV.h"

class Sw_HidDiscRV : public HidDiscRV, public SwDiscRV {
  friend class FileParser;
  friend class CPT;
  friend class MDCPT;
  friend class MSCPT;
  friend class MTCPT;


  virtual void setParents(vector<RV *> &sparents,vector<vector<RV *> > &cpl) {
    setSwitchingConditionalParents(sparents,cpl,this,allParents);
  }
  virtual void setDTMapper(RngDecisionTree *dt) {
    dtMapper = dt;
  }
  virtual vector< RV* >& condParentsVec(unsigned j) {
    assert ( j < conditionalParentsList.size() );
    return conditionalParentsList[j];
  }
  virtual vector< RV* >& switchingParentsVec() {
    return switchingParents;
  }
  virtual void setCpts(vector<CPT*> &cpts) {
    conditionalCPTs = cpts;
  }

protected:

public:

  /////////////////////////////////////////////////////////////////////////
  // constructor: Initialize with the variable type.  The default
  // timeIndex value of ~0x0 indicates a static network.  Discrete
  // nodes must be specified with their cardinalities.
  Sw_HidDiscRV(RVInfo& _rv_info,
	       unsigned _timeFrame = ~0x0,
	       unsigned _cardinality = 0)
    : HidDiscRV(_rv_info,_timeFrame,_cardinality)
  {
    t.hidden = 1;
    t.switching = 1;
    t.scale = 0;
    t.penalty = 0;
    t.shift = 0;
  }

  void printSelf(FILE*f,bool nl=true);
  void printSelfVerbose(FILE*f);

  inline virtual void probGivenParents(logpr& p) {
    setCurrentConditionalParents(this);
    curCPT = conditionalCPTs[cachedSwitchingState];
    p = curCPT->probGivenParents(*curConditionalParents,this);
  }

  inline virtual void begin(logpr& p) {
    setCurrentConditionalParents(this);
    curCPT = conditionalCPTs[cachedSwitchingState];
    curCPT->becomeAwareOfParentValuesAndIterBegin(*curConditionalParents,it,this,p);
    return;
  }

  virtual bool next(logpr& p) { 
    return curCPT->next(it,p);
  }


  void emIncrement(logpr posterior) {
    setCurrentConditionalParents(this);
    curCPT->emIncrement(posterior,*curConditionalParents,this);
  }


  unsigned averageCardinality() { return SwDiscRV::averageCardinality(rv_info); }
  unsigned maxCardinality() { return SwDiscRV::maxCardinality(rv_info); }


  virtual Sw_HidDiscRV* cloneRVShell();
  virtual Sw_HidDiscRV* create() {
    Sw_HidDiscRV*rv = new Sw_HidDiscRV(rv_info,frame(),cardinality);
    rv->t = t;
    return rv;
  }

  virtual void randomSample() {
    setCurrentConditionalParents(this);
    curCPT->becomeAwareOfParentValues( *curConditionalParents, this );
    curCPT->randomSample(this); 
  }


  ////////////////////////////////////////////////////////////////////////
  // Increment the statistics with probabilty 'posterior' for the
  // current random variable's parameters, for the case where
  // the random variable and its parents are set to their current
  // values (i.e., the increment corresponds to the currently set
  // parent/child values). 


};


#endif
