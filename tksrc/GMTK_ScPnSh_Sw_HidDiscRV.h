/*
 * GMTK_ScPnSh_Sw_HidDiscRV.h
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

#ifndef GMTK_SC_PN_SH_SW_HID_DISC_RV_H
#define GMTK_SC_PN_SH_SW_HID_DISC_RV_H

#include <vector>

#include "GMTK_Sw_HidDiscRV.h"
#include "GMTK_ScPnShRV.h"

class ScPnSh_Sw_HidDiscRV : public Sw_HidDiscRV, public ScPnShRV {
  friend class FileParser;

protected:

public:

  /////////////////////////////////////////////////////////////////////////
  // constructor: Initialize with the variable type.  The default
  // timeIndex value of ~0x0 indicates a static network.  Discrete
  // nodes must be specified with their cardinalities.
  ScPnSh_Sw_HidDiscRV(RVInfo& _rv_info,
		   unsigned _timeFrame = ~0x0,
		   unsigned _cardinality = 0)
    : Sw_HidDiscRV(_rv_info,_timeFrame,_cardinality)
  {
  }

  virtual void printSelf(FILE *f,bool nl=true) {
    Sw_HidDiscRV::printSelf(f,false);
    ScPnShRV::printSelf(rv_info.rvWeightInfo[0],f,nl);
  }

  virtual void printSelfVerbose(FILE *f) {
    fprintf(f,"Switching Hidden Discrete Random variable:\n");
    printNameFrameValue(f,true);
    fprintf(f,"From line %d in file %s\n",rv_info.fileLineNumber,rv_info.rvFileName.c_str());
    fprintf(f,"RV has cardinality = %d\n",cardinality);
    ScPnShRV::printSelfVerbose(rv_info.rvWeightInfo[0],f);
  }
  
  virtual void begin(logpr& p) {
    Sw_HidDiscRV::begin(p);
    if (rv_info.rvWeightInfo.size() > 1) 
      modifyProbability(p,rv_info.rvWeightInfo[cachedSwitchingState],this);
    else 
      modifyProbability(p,rv_info.rvWeightInfo[0],this);
    return;
  }

  virtual bool next(logpr& p) { 
    if (Sw_HidDiscRV::next(p)) {
      if (rv_info.rvWeightInfo.size() > 1) 
	modifyProbability(p,rv_info.rvWeightInfo[cachedSwitchingState],this);
      else 
	modifyProbability(p,rv_info.rvWeightInfo[0],this);
      return true;
    } else
      return false;
  }

  inline virtual void probGivenParents(logpr& p) {
    setCurrentConditionalParents(this);
    curCPT = conditionalCPTs[cachedSwitchingState];
    p = curCPT->probGivenParents(*curConditionalParents,this);
    if (rv_info.rvWeightInfo.size() > 1) 
      modifyProbability(p,rv_info.rvWeightInfo[cachedSwitchingState],this);
    else 
      modifyProbability(p,rv_info.rvWeightInfo[0],this);
  }

  virtual ScPnSh_Sw_HidDiscRV* cloneRVShell() {
    return (ScPnSh_Sw_HidDiscRV*)Sw_HidDiscRV::cloneRVShell();
  }
  virtual ScPnSh_Sw_HidDiscRV* create() {
    ScPnSh_Sw_HidDiscRV*rv = new ScPnSh_Sw_HidDiscRV(rv_info,frame(),cardinality);
    return rv;
  }

};

/*

* To reduce branchs further, we can do something like the following for
* each possible condition:

class ScPnSh_Sw_HidDiscRV : public Sw_HidDiscRV, public ScPnShRV {
  friend class FileParser;

protected:

public:

  /////////////////////////////////////////////////////////////////////////
  // constructor: Initialize with the variable type.  The default
  // timeIndex value of ~0x0 indicates a static network.  Discrete
  // nodes must be specified with their cardinalities.
  ScPnSh_Sw_HidDiscRV(RVInfo& _rv_info,
		   unsigned _timeFrame = ~0x0,
		   unsigned _cardinality = 0)
    : Sw_HidDiscRV(_rv_info,_timeFrame,_cardinality)
  {
  }

  virtual void printSelf(FILE *f,bool nl=true) {
    Sw_HidDiscRV::printSelf(f,false);
    ScPnShRV::printSelf(rv_info.rvWeightInfo[0],f,nl);
  }

  virtual void printSelfVerbose(FILE *f) {
    fprintf(f,"Switching Hidden Discrete Random variable:\n");
    printNameFrameValue(f,true);
    fprintf(f,"From line %d in file %s\n",rv_info.fileLineNumber,rv_info.rvFileName.c_str());
    fprintf(f,"RV has cardinality = %d\n",cardinality);
    ScPnShRV::printSelfVerbose(rv_info.rvWeightInfo[0],f,true);
  }
  
  virtual void begin(logpr& p) = 0; 
  virtual bool next(logpr& p) = 0;
  virtual void probGivenParents(logpr& p) = 0;

  virtual ScPnSh_Sw_HidDiscRV* cloneRVShell() {
    return (ScPnSh_Sw_HidDiscRV*)Sw_HidDiscRV::cloneRVShell();
  }

};


#define DEFINE_SCPNSH_SW_HIDDISCRV_CLASS(_COND_) \
\
class ScPnSh_## _COND_ ## _Sw_HidDiscRV : public ScPnSh_Sw_HidDiscRV, public ScPnShRV { \
  friend class FileParser;  \
protected:  \
public:  \
  ScPnSh_ ## _COND_ ## _Sw_HidDiscRV(RVInfo& _rv_info,  \
		   unsigned _timeFrame = ~0x0,  \
		   unsigned _cardinality = 0)  \
    : Sw_HidDiscRV(_rv_info,_timeFrame,_cardinality)  \
  {  \
  }  \
  \
  \
  virtual void begin(logpr& p) {  \
    Sw_HidDiscRV::begin(p);  \
    if (rv_info.rvWeightInfo.size() > 1)   \
      modifyProbability ## _COND_ (p,rv_info.rvWeightInfo[cachedSwitchingState],this);  \
    else   \
      modifyProbability ## _COND_ (p,rv_info.rvWeightInfo[0],this);  \
    return;  \
  }  \
  \
  virtual bool next(logpr& p) {   \
    if (Sw_HidDiscRV::next(p)) {  \
      if (rv_info.rvWeightInfo.size() > 1)   \
	modifyProbability ## _COND_ (p,rv_info.rvWeightInfo[cachedSwitchingState],this);  \
      else   \
	modifyProbability ## _COND_(p,rv_info.rvWeightInfo[0],this);  \
      return true;  \
    } else  \
      return false;  \
  }  \
  \
  inline virtual void probGivenParents(logpr& p) {  \
    setCurrentConditionalParents(this);  \
    curCPT = conditionalCPTs[cachedSwitchingState];  \
    p = curCPT->probGivenParents(*curConditionalParents,this);  \
    if (rv_info.rvWeightInfo.size() > 1)   \
      modifyProbability ## _COND_ (p,rv_info.rvWeightInfo[cachedSwitchingState],this);  \
    else   \
      modifyProbability ## _COND_ (p,rv_info.rvWeightInfo[0],this);  \
  }  \
  \
  virtual ScPnSh_Sw_HidDiscRV* create() {  \
    ScPnSh_Sw_HidDiscRV*rv = new ScPnSh_Sw_HidDiscRV(rv_info,frame(),cardinality);  \
    rv->t = t;  \
    return rv;  \
  }  \
  \
};

DEFINE_SCPNSH_SW_HIDDISCRV_CLASS(CO)

*/

#endif
