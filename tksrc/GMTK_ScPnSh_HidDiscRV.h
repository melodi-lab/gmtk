/*
 * GMTK_HidDiscRV.h
 *
 *  Scale/Penalty/Shift Hidden Discrete Random Variables.
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

#ifndef GMTK_SC_PN_SH_HID_DISC_RV_H
#define GMTK_SC_PN_SH_HID_DISC_RV_H

#include <vector>

#include "GMTK_HidDiscRV.h"
#include "GMTK_ScPnShRV.h"

class ScPnSh_HidDiscRV : public HidDiscRV, public ScPnShRV {
  friend class FileParser;

protected:

public:

  /////////////////////////////////////////////////////////////////////////
  // constructor: Initialize with the variable type.  The default
  // timeIndex value of ~0x0 indicates a static network.  Discrete
  // nodes must be specified with their cardinalities.
  ScPnSh_HidDiscRV(RVInfo& _rv_info,
		   unsigned _timeFrame = ~0x0,
		   unsigned _cardinality = 0)
    : HidDiscRV(_rv_info,_timeFrame,_cardinality)
  {
  }

  virtual void printSelf(FILE *f,bool nl=true) {
    HidDiscRV::printSelf(f,false);
    ScPnShRV::printSelf(rv_info.rvWeightInfo[0],f,nl);
  }

  virtual void printSelfVerbose(FILE *f) {
    fprintf(f,"Hidden Discrete Random variable:\n");
    printNameFrameValue(f,true);
    fprintf(f,"From line %d in file %s\n",rv_info.fileLineNumber,rv_info.rvFileName.c_str());
    fprintf(f,"RV has cardinality = %d\n",cardinality);
    ScPnShRV::printSelfVerbose(rv_info.rvWeightInfo[0],f);
  }
  
  virtual void begin(logpr& p) {
    HidDiscRV::begin(p);
    modifyProbability(p,rv_info.rvWeightInfo[0],this);
    return;
  }

  virtual bool next(logpr& p) { 
      if (HidDiscRV::next(p)) {
	modifyProbability(p,rv_info.rvWeightInfo[0],this);
	return true;
      } else
	return false;
  }

  inline virtual void probGivenParents(logpr& p) {
    p = curCPT->probGivenParents(allParents,this);
    modifyProbability(p,rv_info.rvWeightInfo[0],this);
  }


  virtual ScPnSh_HidDiscRV* cloneRVShell() {
    return (ScPnSh_HidDiscRV*)HidDiscRV::cloneRVShell();
  }
  virtual ScPnSh_HidDiscRV* create() {
    ScPnSh_HidDiscRV*rv = new ScPnSh_HidDiscRV(rv_info,frame(),cardinality);
    return rv;
  }

};


#endif
