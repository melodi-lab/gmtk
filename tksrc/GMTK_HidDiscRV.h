/*
 * GMTK_HidDiscRV.h
 *
 *  Hidden Discrete Random Variables.
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

#ifndef GMTK_HID_DISC_RV_H
#define GMTK_HID_DISC_RV_H

#include <vector>

#include "GMTK_DiscRV.h"


class HidDiscRV : public DiscRV {
  friend class FileParser;
  friend class CPT;
  friend class MDCPT;
  friend class MSCPT;
  friend class MTCPT;

protected:

public:

  /////////////////////////////////////////////////////////////////////////
  // constructor: Initialize with the variable type.  The default
  // timeIndex value of ~0x0 indicates a static network.  Discrete
  // nodes must be specified with their cardinalities.
  HidDiscRV(RVInfo& _rv_info,
	    unsigned _timeFrame = ~0x0,
	    unsigned _cardinality = 0)
    : DiscRV(_rv_info,_timeFrame,_cardinality)
  {
  }

  virtual void printSelf(FILE *f,bool nl=true);
  virtual void printSelfVerbose(FILE *f);
  
  virtual void begin(logpr& p) {
    curCPT->becomeAwareOfParentValuesAndIterBegin(allParents,it,this,p);
    return;
  }
  virtual bool next(logpr& p) { 
    return curCPT->next(it,p);
  }

  // This function assumes that::
  //   1) deterministic() is true
  //   2) the curCPT is a deterministic CPT
  // If these conditions are not true, calling this function will yield a run-time error.
  void assignDeterministicChild() { curCPT->assignDeterministicChild(allParents,this); }

  virtual HidDiscRV* cloneRVShell();
  virtual HidDiscRV* create() {
    HidDiscRV*rv = new HidDiscRV(rv_info,frame(),cardinality);
    return rv;
  }


};


#endif
