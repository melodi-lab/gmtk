/*
 * GMTK_SwDiscRV.h
 *
 *  Switching functionality for a RV, but adds parameters that exist
 *  both in hidden and obseved versions of this variable.
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
 * The top level GMTK random variable object for the RV class hierarchy.
 *
 *
 *
 */

#ifndef GMTK_SW_DISC_RV_H
#define GMTK_SW_DISC_RV_H

#include <vector>
#include <string>
#include <set>

#include "logp.h"
#include "GMTK_RVInfo.h"
#include "GMTK_NamedObject.h"
#include "GMTK_RngDecisionTree.h"

#include "GMTK_SwRV.h"

class SwDiscRV : public SwRV {

  friend class FileParser;

protected:

  //////////////////////////////////////////////////////////////////////
  // CPT array, one for each set of possible parents we might
  // have (so size of this array is the number of different
  // possible conditional parents).
  vector < CPT* > conditionalCPTs;

public:

  SwDiscRV() {}
  virtual ~SwDiscRV() {}

  ////////////////////////////////////////////////////////////////////////
  // Ties the parameters of 'this' with whatever those of 'other' are. 
  void tieParametersWith(SwDiscRV* other) {
    SwRV::tieParametersWith(other);
    conditionalCPTs = other->conditionalCPTs;
  }

  unsigned averageCardinality(RVInfo& rv_info);
  unsigned maxCardinality(RVInfo& rv_info);

  bool iterableSw() const { 
    for (unsigned i=0;i<conditionalCPTs.size();i++)
      if (conditionalCPTs[i]->iterable())
	return false;
    return dtMapper->iterable();
  }

  logpr maxValue() {
    logpr res;
    for (unsigned i=0;i<conditionalCPTs.size();i++) {
      logpr tmp = conditionalCPTs[i]->maxValue();
      if (tmp > res)
	res = tmp;
    }
    return res;
  }

};

#endif
