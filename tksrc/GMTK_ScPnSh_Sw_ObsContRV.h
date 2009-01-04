/*
 * GMTK_ScPnSh_Sw_ObsContRV.h
 *
 *  Scale/Penalty/Shift Switching Observed Continuous Random Variables.
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

#ifndef GMTK_SC_PN_SH_SW_OBS_CONT_RV_H
#define GMTK_SC_PN_SH_SW_OBS_CONT_RV_H

#include <vector>

#include "GMTK_Sw_ObsContRV.h"
#include "GMTK_ScPnShRV.h"

class ScPnSh_Sw_ObsContRV : public Sw_ObsContRV, public ScPnShRV {
  friend class FileParser;

protected:

public:

  /////////////////////////////////////////////////////////////////////////
  // constructor: Initialize with the variable type.  The default
  // timeIndex value of ~0x0 indicates a static network.  Discrete
  // nodes must be specified with their cardinalities.
  ScPnSh_Sw_ObsContRV(RVInfo& _rv_info,
		      unsigned _timeFrame = ~0x0)
    : Sw_ObsContRV(_rv_info,_timeFrame)
  {
  }

  virtual void printSelf(FILE *f,bool nl=true) {
    Sw_ObsContRV::printSelf(f,false);
    ScPnShRV::printSelf(rv_info.rvWeightInfo[0],f,nl);
  }

  virtual void printSelfVerbose(FILE *f) 
  {
    Sw_ObsContRV::printSelfVerbose(f);
    ScPnShRV::printSelfVerbose(rv_info.rvWeightInfo[0],f);
  }

  virtual void begin(logpr& p) {
    Sw_ObsContRV::begin(p);
    if (rv_info.rvWeightInfo.size() > 1) 
      modifyProbability(p,rv_info.rvWeightInfo[cachedSwitchingState],this);
    else 
      modifyProbability(p,rv_info.rvWeightInfo[0],this);
    return;
  }

  inline virtual void probGivenParents(logpr& p) {
    Sw_ObsContRV::probGivenParents(p);
    if (rv_info.rvWeightInfo.size() > 1) 
      modifyProbability(p,rv_info.rvWeightInfo[cachedSwitchingState],this);
    else 
      modifyProbability(p,rv_info.rvWeightInfo[0],this);
  }

  virtual logpr maxValue() {
    logpr p = Sw_ObsContRV::maxValue();
    if (rv_info.rvWeightInfo.size() > 1) 
      modifyProbability(p,rv_info.rvWeightInfo[cachedSwitchingState],this);
    else 
      modifyProbability(p,rv_info.rvWeightInfo[0],this);
    return p;
  }

  virtual ScPnSh_Sw_ObsContRV* cloneRVShell() {
    return (ScPnSh_Sw_ObsContRV*)Sw_ObsContRV::cloneRVShell();
  }
  virtual ScPnSh_Sw_ObsContRV* create() {
    ScPnSh_Sw_ObsContRV*rv = new ScPnSh_Sw_ObsContRV(rv_info,frame());
    return rv;
  }

};


#endif
