/*
 * GMTK_SwRV.cc
 *
 * Switching support functionality for random variables with switching parents.
 * 
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

#include "GMTK_SwRV.h"
#include "GMTK_DiscRV.h"

#include "GMTK_RngDecisionTree.h"


/*-
 *-----------------------------------------------------------------------
 * Function
 *
 *     setSwitchingConditionalParents(): set the switching and
 *     conditional parents of this object, and also set the union into
 *     allParents.
 *
 * Preconditions:
 *      'this' Variable must have switching.
 *
 * Postconditions:
 *      parent are re-set.
 *
 * Side Effects:
 *      parent are re-set.
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
void SwRV::setSwitchingConditionalParents(vector<RV *> &sparents,
					  vector<vector<RV *> > &cpl,
					  RV* self,
					  vector<RV*>& allParents)
{
  // make sure we give this swithcing information.
  assert ( sparents.size() > 0 );
  assert ( cpl.size() > 1 );


  switchingParents = sparents;
  conditionalParentsList = cpl;

  // now set this as a child of all parents, making sure to avoid
  // duplicates by creating a temporary set.
  set<RV *> parentSet;
  for (unsigned i=0;i<sparents.size();i++) {
    parentSet.insert(sparents[i]);
  }
  for (unsigned i=0;i<cpl.size();i++) {
    for (unsigned j=0;j<cpl[i].size();j++) {
      parentSet.insert(cpl[i][j]);
    }
  }
  set<RV *>::iterator si;
  allParents.clear();
  for (si = parentSet.begin(); si != parentSet.end(); si++) {
    RV* rv = (*si);
    allParents.push_back(rv);
    rv->allChildren.push_back(self);
  }
}


/*-
 *-----------------------------------------------------------------------
 * Function
 *
 *     setConditionalParents(): sets the current set of conditional
 *     parents based on the switching parents. Note that this routine
 *     takes as an argument the child random variable since even
 *     though 'this' is the child, we can not cast from 'this' to the
 *     true base class RV (due to multiple inheritance and a desire
 *     not to use RV as a virtual base class).
 * 
 * Preconditions:
 *      all variable must be appropriately assigned.
 *
 * Postconditions:
 *      The current set of condtional parents are set appropriately
 *
 * Side Effects:
 *      Changes member function.
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
void
SwRV::setCurrentConditionalParents(RV* rv)
{
  assert (dtMapper != NULL);
  cachedSwitchingState = dtMapper->query(switchingParents,rv);
  if ( cachedSwitchingState >= conditionalParentsList.size()) {
    warning("ERROR: Random Variable %s:%d using Decision Tree '%s' yielded an invalid switching position %d. Must be between 0 and %d.\n",
	    rv->name().c_str(),rv->frame(), (dtMapper == NULL?"NULL":dtMapper->name().c_str()),cachedSwitchingState,conditionalParentsList.size());
    fprintf(stderr,"Current values of switching parents: ");
    // TODO: use vector RV print routine.
    for (unsigned i=0;i<switchingParents.size();i++) {
      switchingParents[i]->printNameFrameValue(stderr,false);
    }
    error("\n");
  }
  curConditionalParents = & conditionalParentsList[cachedSwitchingState];
}


