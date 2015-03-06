/*
 * GMTK_ObsContRV.h
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

#ifndef GMTK_SC_PN_SH_OBS_CONT_RV_H
#define GMTK_SC_PN_SH_OBS_CONT_RV_H

#include <vector>

#include "GMTK_ObsContRV.h"
#include "GMTK_ScPnShRV.h"
#include "GMTK_MixtureCommon.h"
#include "GMTK_Mixture.h"

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
    // ticket 6: This code is copied from ObsContRV::maxValue(), but 
    // modified to include the scale/penalty/shift if they are available.
    logpr mval, tmp;
    assert(rv_info.rvWeightInfo.size() == 1); // this isn't a switched RV
    for (unsigned i=0; i< conditionalMixtures.size(); i++) { // i is the switching state
      // this loop should just iterate 1 time, since this isn't a switched RV
      RVInfo::WeightInfo &wi = rv_info.rvWeightInfo[0];
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
    // ticket 6:
    // modifyProbability() may require observation data that's not available yet
    // modifyProbability() should be done before taking the max
    logpr p = ObsContRV::maxValue();
    modifyProbability(p,rv_info.rvWeightInfo[0],this);
    return p;
#endif
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
