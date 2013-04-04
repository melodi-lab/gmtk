/*
 * GMTK_ObsDiscRV.h
 *
 *  Observed Discrete Random Variables.
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *  $Header$
 *
 * Copyright (C) 2001 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 *
 *
 *
 * The discrete random variable type.
 *
 *
 *
 */

#ifndef GMTK_OBS_DISC_RV_H
#define GMTK_OBS_DISC_RV_H

#include <vector>

#include "GMTK_DiscRV.h"
#if 1
#  include "GMTK_FileSource.h"
#endif

class ObsDiscRV : public DiscRV {
  friend class FileParser;
  friend class CPT;
  friend class MDCPT;
  friend class MSCPT;
  friend class MTCPT;

protected:

  // the feature file element corresponding to this RV in the
  // integer portion of the global observation matrix.
  unsigned featureElement() { return rv_info.rvFeatureRange.firstFeatureElement; }


public:

  /////////////////////////////////////////////////////////////////////////
  // constructor: Initialize with the variable type.  The default
  // timeIndex value of ~0x0 indicates a static network.  Discrete
  // nodes must be specified with their cardinalities.
  ObsDiscRV(RVInfo& _rv_info,
	    unsigned _timeFrame = ~0x0,
	    unsigned _cardinality = 0)
    : DiscRV(_rv_info,_timeFrame,_cardinality)
  {
  }
  virtual ~ObsDiscRV() {;}

  virtual void printSelf(FILE *f,bool nl=true);
  virtual void printSelfVerbose(FILE *f);

  
  virtual void begin(logpr& p) {
    DiscRV::probGivenParents(p);
    return;
  }
  virtual bool next(logpr& p) { return false; }

  virtual ObsDiscRV* cloneRVShell();
  virtual ObsDiscRV* create() {
    ObsDiscRV*rv = new ObsDiscRV(rv_info,frame(),cardinality);
    return rv;
  }

  // set the RV (which must be observed) to the observed value in a
  // file (this should be done once outside of the inference inner
  // loops).
  void setToObservedValue();

  // only valid when this var is non-switching, observed, with a deterministic implementation.
  void computeParentsSatisfyingChild(
	    // input arguments
	    unsigned par, // parent number
	    vector <RV*> & parents, 
	    vector <RV*> & hiddenParents,
	    PackCliqueValue& hiddenParentPacker,
	    sArray < DiscRVType*>& hiddenNodeValPtrs,
	    RV* child,
	    // output arguments
	    sArray < unsigned >& packedParentVals,
	    unsigned& num);

};

void setObservedRVs(vector <RV*>& rvs);
void setObservedRVs(set <RV*>& rvs);

#endif
