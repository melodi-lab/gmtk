/*
 * GMTK_ScPnSh_Sw_ObsContRV.h
 *
 *  Scale/Penalty/Shift Switching Observed Continuous Random Variables.
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

#ifndef GMTK_SC_PN_SH_SW_OBS_CONT_RV_H
#define GMTK_SC_PN_SH_SW_OBS_CONT_RV_H

#include <vector>

#include "GMTK_Sw_ObsContRV.h"
#include "GMTK_ScPnShRV.h"
#include "GMTK_MixtureCommon.h"
#include "GMTK_Mixture.h"

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
    // ticket 6: This code is copied from ObsContRV::maxValue(), but 
    // modified to include the scale/penalty/shift if they are available.
    logpr mval, tmp;
    for (unsigned i=0; i< conditionalMixtures.size(); i++) { // i is the switching state
      RVInfo::WeightInfo &wi = (rv_info.rvWeightInfo.size() > 1) ?
	  rv_info.rvWeightInfo[i]
	: rv_info.rvWeightInfo[0];
      if (conditionalMixtures[i].direct) {
	tmp = conditionalMixtures[i].mixture->maxValue();
	if (safeToModifyProbability(wi)) 
	  modifyProbability(tmp, wi, this);
	if (tmp > mval)
	  mval = tmp;
      } else {
	for (unsigned j=0; j < conditionalMixtures[i].mapping.collection->mxSize(); j++) {
	  tmp = conditionalMixtures[i].mapping.collection->mx(j)->maxValue();
	  if (safeToModifyProbability(wi)) 
	    modifyProbability(tmp, wi, this);
	  if (tmp > mval)
	    mval = tmp;
	}
      }
    }
    return mval;
#if 0
      // ticket 6: cachedSwitchingState may not be initialized before maxValue() is called
      // modifyProbability() may require observation data that's not available yet
      // modifyProbability() should be done before taking the max
      logpr p = Sw_ObsContRV::maxValue();
      if (rv_info.rvWeightInfo.size() > 1) 
	modifyProbability(p,rv_info.rvWeightInfo[cachedSwitchingState],this);
      else 
	modifyProbability(p,rv_info.rvWeightInfo[0],this);
      return p;
#endif
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
