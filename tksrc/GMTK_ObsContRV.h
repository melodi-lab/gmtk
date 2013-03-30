/*
 * GMTK_ObsContRV.h
 *
 * Observed Continuous Random Variable
 *  
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *  $Header$
 *
 * Copyright (C) 2001 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 *
 *
 *
 */

#ifndef GMTK_OBS_CONT_RV_H
#define GMTK_OBS_CONT_RV_H

#include <vector>
#include <string>
#include <set>

#include "GMTK_ContRV.h"
#include "GMTK_CPT.h"
#if 0
#  include "GMTK_ObservationMatrix.h"
#else
#  include "GMTK_FileSource.h"
#endif
#include "GMTK_MixtureCommon.h"
#include "GMTK_GMParms.h"
#include "GMTK_NameCollection.h"

class FileParser;

class ObsContRV : public ContRV {
  friend class FileParser;

protected:


public:
  
  /////////////////////////////////////////////////////////////////////////
  // Constructor: Initialize with the variable type.  The default
  // timeIndex value of ~0x0 indicates a static network.  Discrete
  // nodes must be specified with their cardinalities.
  ObsContRV(RVInfo& _rv_info,
	    unsigned _timeFrame = ~0x0)
    : ContRV(_rv_info,_timeFrame)
  {
  }
  virtual ~ObsContRV() {;}

  /////////////////////////////////////////////////////////////////
  // The feature range to which this random variable corresponds in a
  // data file. it corresponds to
  // [firstFeatureElement:lastFeatureElement] inclusive.
  unsigned firstFeatureElement() { return rv_info.rvFeatureRange.firstFeatureElement; }
  unsigned lastFeatureElement() { return rv_info.rvFeatureRange.lastFeatureElement; }

  virtual void printNameFrameValue(FILE *f,bool nl=true);
  virtual void printSelf(FILE *f,bool nl=true);
  virtual void printSelfVerbose(FILE *f);
  
  void probGivenParents(logpr& p);
  logpr maxValue();

  void begin(logpr& p) {
    ObsContRV::probGivenParents(p);
    return;
  }
  bool next(logpr& p) { return false; }

  void makeRandom() {  error("not implemented, this should be called to somewhere else"); }
  void makeUniform(){  error("not implemented, this should be called to somewhere else"); }

  ///////////////////////////////////////////////////
  // EM Support
  void emIncrement(logpr prob);

  void randomSample() { error("ERROR: random sampling from observed continous random variable not yet implemented"); }

  virtual ObsContRV* cloneRVShell();
  virtual ObsContRV* create() {
    ObsContRV*rv = new ObsContRV(rv_info,frame());
    return rv;
  }



};



#endif


