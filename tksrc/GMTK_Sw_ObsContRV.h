/*
 * GMTK_Sw_ObsContRV.h
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

#ifndef GMTK_SW_OBS_CONT_RV_H
#define GMTK_SW_OBS_CONT_RV_H

#include <vector>

#include "GMTK_ObsContRV.h"
#include "GMTK_SwContRV.h"

class Sw_ObsContRV : public ObsContRV, public SwContRV {
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

protected:

public:

  /////////////////////////////////////////////////////////////////////////
  // constructor: Initialize with the variable type.  The default
  // timeIndex value of ~0x0 indicates a static network.  Discrete
  // nodes must be specified with their cardinalities.
  Sw_ObsContRV(RVInfo& _rv_info,
	       unsigned _timeFrame = ~0x0)
    : ObsContRV(_rv_info,_timeFrame)
  {
  }

  void printSelf(FILE*f,bool nl=true);
  void printSelfVerbose(FILE*f);

  void probGivenParents(logpr& p);

  void begin(logpr& p) {
    probGivenParents(p);
    return;
  }
  bool next(logpr& p) { return false; }
  void emIncrement(logpr posterior);

  virtual Sw_ObsContRV* cloneRVShell();
  virtual Sw_ObsContRV* create() {
    Sw_ObsContRV*rv = new Sw_ObsContRV(rv_info,frame());
    return rv;
  }

  bool iterable() const {
    return iterableSw();
  }


};


#endif
