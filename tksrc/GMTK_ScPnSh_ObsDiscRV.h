/*
 * GMTK_ObsDiscRV.h
 *
 *  Scale/Penalty/Shift Observed Discrete Random Variables.
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *  $Header$
 *
 * Copyright (C) 2001 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 *
 *
 *
 */

#ifndef GMTK_SC_PN_SH_OBS_DISC_RV_H
#define GMTK_SC_PN_SH_OBS_DISC_RV_H

#include <vector>

#include "GMTK_ObsDiscRV.h"
#include "GMTK_ScPnShRV.h"

class ScPnSh_ObsDiscRV : public ObsDiscRV, public ScPnShRV {
  friend class FileParser;

protected:

public:

  /////////////////////////////////////////////////////////////////////////
  // constructor: Initialize with the variable type.  The default
  // timeIndex value of ~0x0 indicates a static network.  Discrete
  // nodes must be specified with their cardinalities.
  ScPnSh_ObsDiscRV(RVInfo& _rv_info,
		   unsigned _timeFrame = ~0x0,
		   unsigned _cardinality = 0)
    : ObsDiscRV(_rv_info,_timeFrame,_cardinality)
  {
  }

  virtual void printSelf(FILE *f,bool nl=true) {
    ObsDiscRV::printSelf(f,false);
    ScPnShRV::printSelf(rv_info.rvWeightInfo[0],f,nl);
  }

  virtual void printSelfVerbose(FILE *f) {
    fprintf(f,"Observed Discrete Random variable:\n");
    printNameFrameValue(f,true);
    fprintf(f,"From line %d in file %s\n",rv_info.fileLineNumber,rv_info.rvFileName.c_str());
    fprintf(f,"RV has cardinality = %d\n",cardinality);
    ScPnShRV::printSelf(rv_info.rvWeightInfo[0],f);
  }
  
  virtual void begin(logpr& p) {
    ObsDiscRV::begin(p);
    modifyProbability(p,rv_info.rvWeightInfo[0],this);
    return;
  }

  inline virtual void probGivenParents(logpr& p) {
    p = curCPT->probGivenParents(allParents,this);
    modifyProbability(p,rv_info.rvWeightInfo[0],this);
  }

  // See https://j.ee.washington.edu/trac/gmtk/ticket/6#comment:25
  // maxValue() is inherited from a base class, so the scale, shift, and penalty
  // are not applied. I'm not sure if that matters since the RV is observed. - RR

  virtual ScPnSh_ObsDiscRV* cloneRVShell() {
    return (ScPnSh_ObsDiscRV*)ObsDiscRV::cloneRVShell();
  }
  virtual ScPnSh_ObsDiscRV* create() {
    ScPnSh_ObsDiscRV*rv = new ScPnSh_ObsDiscRV(rv_info,frame(),cardinality);
    return rv;
  }

};


#endif
