/*
 * GMTK_SwRV.h
 *
 *  Switching functionality for a RV.
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

#ifndef GMTK_SW_RV_H
#define GMTK_SW_RV_H

#include <vector>
#include <string>
#include <set>

#include "logp.h"
#include "GMTK_RVInfo.h"
#include "GMTK_NamedObject.h"
#include "GMTK_RngDecisionTree.h"

#include "GMTK_DiscRV.h"

class SwRV {

  friend class FileParser;

protected:

  ////////////////////////////////////////////////////////////////////////
  // The set of (assumed to be necessarily discrete) switching parents
  // of a RV that has switching parents.
  vector< RV *> switchingParents;

  ////////////////////////////////////////////////////////////////////////
  // For each possible different list of conditional parents that
  // might exist for all possible values of the switching parents,
  // this array gives that list of appropriate conditional
  // parents. For example, suppose that S is the set of conditional
  // parents, and that 0 <= S <= 5 corresponds to one set of
  // conditional parents, and 6 <= S < = 10 corresponds to another set
  // of conditional parents, and those are the only two set of
  // conditional parents that exist for all values of the switching
  // parents, this list is of size two.
  vector< vector < RV* > > conditionalParentsList;

  ////////////////////////////////////////////////////////////////////////
  // This variable is assigned a pointer to the current set of
  // conditional parents, which is dependent on the current value of
  // the switching parents.  Note that this points to one of the
  // entries in conditionalParentsList
  vector<RV *> *curConditionalParents;

  ////////////////////////////////////////////////////////////////////////
  // This decision tree is used to map from the set of conditional
  // parents to the integer used to select the curret set of
  // conditional parents via the DT query of the switching parents.
  // The routine maps from the current set of switching parent values
  // to an integer, which indicates which set of conditional parents
  // should be active for those switching parent values.
  RngDecisionTree *dtMapper;

  // set the switching and conditional parents of this object,
  // and also set the union into allParents.
  void setSwitchingConditionalParents(vector<RV *> &sparents,
				      vector<vector<RV *> > &cpl,
				      RV* self,
				      vector<RV*>& allParents);

public:

  SwRV() {}
  virtual ~SwRV() {}

  // Cached value of the most recent swithcing parents query via the
  // DT. We can reuse this value w/o needing to do the DT integer map
  // lookup again.
  unsigned cachedSwitchingState;

  ////////////////////////////////////////////////////////////////////////
  // Based on the values of the switching parents, this routine determines 
  // the appropriate value of the conditionalParents, the current set
  // of conditional parents.
  void setCurrentConditionalParents(RV* rv);

  void tieParametersWith(SwRV* other) {
    dtMapper = other->dtMapper;
  }

};



#endif
