/*
 * GMTK_ObsContRV.h
 *
 *  Scale/Penalty/Shift Observed Discrete Random Variables.
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
 *
 */

#ifndef GMTK_SC_PN_SH_OBS_CONT_RV_H
#define GMTK_SC_PN_SH_OBS_CONT_RV_H

#include <vector>

#include "GMTK_ObsContRV.h"
#include "GMTK_ScPnShRV.h"

class ScPnSh_ObsContRV : public ObsContRV, public ScPnShRV {
  friend class FileParser;

protected:

public:

  /////////////////////////////////////////////////////////////////////////
  // constructor: Initialize with the variable type.  The default
  // timeIndex value of ~0x0 indicates a static network.  Discrete
  // nodes must be specified with their cardinalities.
  ScPnSh_ObsContRV(RVInfo& _rv_info,
		   unsigned _timeFrame = ~0x0)
    : ObsContRV(_rv_info,_timeFrame)
  {
  }

  virtual void printSelf(FILE *f,bool nl=true) {
    ObsContRV::printSelf(f,false);
    ScPnShRV::printSelf(rv_info.rvWeightInfo[0],f,nl);
  }

  virtual void printSelfVerbose(FILE *f) 
  {
    ObsContRV::printSelfVerbose(f);
    ScPnShRV::printSelfVerbose(rv_info.rvWeightInfo[0],f);
  }

  virtual void begin(logpr& p) {
    ObsContRV::begin(p);
    modifyProbability(p,rv_info.rvWeightInfo[0],this);
    return;
  }

  inline virtual void probGivenParents(logpr& p) {
    ObsContRV::probGivenParents(p);
    modifyProbability(p,rv_info.rvWeightInfo[0],this);
  }
  inline virtual logpr maxValue() {
    logpr p = ObsContRV::maxValue();
    // TODO: a bug here, if we call modifyProbability at the beginning
    // before we have any observations, and if the modification parameters
    // come from an observation file, then this will produce random results.
    // I.e., we might be calling maxValue of this rv to do some pre-pruning
    // stuff here.
    // This is an inherent problem as the max probability of a random variable
    // is very time-dependent (it depends on the obs file). 
    modifyProbability(p,rv_info.rvWeightInfo[0],this);
    return p;
  }


  virtual ScPnSh_ObsContRV* cloneRVShell() {
    return (ScPnSh_ObsContRV*)ObsContRV::cloneRVShell();
  }
  virtual ScPnSh_ObsContRV* create() {
    ScPnSh_ObsContRV*rv = new ScPnSh_ObsContRV(rv_info,frame());
    return rv;
  }

};


#endif
