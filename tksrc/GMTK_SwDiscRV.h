/*
 * GMTK_SwDiscRV.h
 *
 *  Switching functionality for a RV, but adds parameters that exist
 *  both in hidden and obseved versions of this variable.
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


};

#endif
