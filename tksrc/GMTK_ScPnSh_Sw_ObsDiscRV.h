/*
 * GMTK_ScPnSh_Sw_ObsDiscRV.h
 *
 *  Scale/Penalty/Shift Switching Hidden Discrete Random Variables.
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

#ifndef GMTK_SC_PN_SH_SW_OBS_DISC_RV_H
#define GMTK_SC_PN_SH_SW_OBS_DISC_RV_H

#include <vector>

#include "GMTK_Sw_ObsDiscRV.h"
#include "GMTK_ScPnShRV.h"

class ScPnSh_Sw_ObsDiscRV : public Sw_ObsDiscRV, public ScPnShRV {
  friend class FileParser;

protected:

public:

  /////////////////////////////////////////////////////////////////////////
  // constructor: Initialize with the variable type.  The default
  // timeIndex value of ~0x0 indicates a static network.  Discrete
  // nodes must be specified with their cardinalities.
  ScPnSh_Sw_ObsDiscRV(RVInfo& _rv_info,
		   unsigned _timeFrame = ~0x0,
		   unsigned _cardinality = 0)
    : Sw_ObsDiscRV(_rv_info,_timeFrame,_cardinality)
  {
  }

  void printSelf(FILE *f,bool nl=true) {
    printNameFrameValue(f,false);
    ScPnShRV::printSelf(rv_info.rvWeightInfo[0],f,false);
    fprintf(f,"switching observed discrete cardinality = %d%s",cardinality,nls(nl));
  }

  void printSelfVerbose(FILE *f) {
    fprintf(f,"Switching Observed Discrete Random variable:\n");
    printNameFrameValue(f,true);
    fprintf(f,"From line %d in file %s\n",rv_info.fileLineNumber,rv_info.rvFileName.c_str());
    fprintf(f,"RV has cardinality = %d\n",cardinality);
    ScPnShRV::printSelfVerbose(rv_info.rvWeightInfo[0],f);
  }
  
  void begin(logpr& p) {
    Sw_ObsDiscRV::begin(p);
    if (rv_info.rvWeightInfo.size() > 1) 
      modifyProbability(p,rv_info.rvWeightInfo[cachedSwitchingState],this);
    else 
      modifyProbability(p,rv_info.rvWeightInfo[0],this);
    return;
  }

  void probGivenParents(logpr& p) {
    setCurrentConditionalParents(this);
    curCPT = conditionalCPTs[cachedSwitchingState];
    p = curCPT->probGivenParents(*curConditionalParents,this);
    if (rv_info.rvWeightInfo.size() > 1) 
      modifyProbability(p,rv_info.rvWeightInfo[cachedSwitchingState],this);
    else 
      modifyProbability(p,rv_info.rvWeightInfo[0],this);
  }

  ScPnSh_Sw_ObsDiscRV* cloneRVShell() {
    return (ScPnSh_Sw_ObsDiscRV*)Sw_ObsDiscRV::cloneRVShell();
  }
  ScPnSh_Sw_ObsDiscRV* create() {
    ScPnSh_Sw_ObsDiscRV*rv = new ScPnSh_Sw_ObsDiscRV(rv_info,frame(),cardinality);
    return rv;
  }

};


#endif
