/*
 * GMTK_SwDiscRV.cc
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


#include "general.h"
VCID("$Header$")

#include <iostream>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <string.h>

#include "GMTK_SwDiscRV.h"


// This routine returns the average cardinality (average number of
// possible child random variable values with non-zero probability)
// for the current random variable.  This routines assume that the
// RV's parents are in the same clique (or an earlier clique
// relative to the JT root) so that we only consider the child
// values with non-zero probability.
unsigned SwDiscRV::averageCardinality(RVInfo& rv_info) {
  unsigned res = 0;
  for (unsigned i=0;i<rv_info.discImplementations.size();i++) {
    if (rv_info.discImplementations[i] == CPT::di_MDCPT)
	res += rv_info.rvCard;
    else if (rv_info.discImplementations[i] == CPT::di_MTCPT)
	res += 1;
    else { // sparce CPT
	// we only get the sparse CPT cardinality if it has
	// been read in, otherwise, we just take card.
	unsigned tmp;
	if (i < conditionalCPTs.size())
	  tmp = conditionalCPTs[i]->averageCardinality();
	else 
	  tmp = rv_info.rvCard;
	res += tmp;
    }
  }
  if (res == 0) {
    // if there are no CPTs (i.e., if we didn't load in the
    // parameters during say a triangulation procedure, then we just
    // return the cardinality.
    return rv_info.rvCard;
  } else {
    // otherwise, we return what we've got.
    return res/rv_info.discImplementations.size();
  }
}

// This routine returns the maximum possible cardinality (max number
// of possible child random variable values with non-zero
// probability) for the current random variable.  This routines
// assume that the RV's parents are in the same clique (or an
// earlier clique relative to the JT root) so that we only consider
// the child values with non-zero probability.
unsigned SwDiscRV::maxCardinality(RVInfo& rv_info) {
  unsigned res = 0;
  for (unsigned i=0;i<rv_info.discImplementations.size();i++) {
    if (rv_info.discImplementations[i] == CPT::di_MDCPT) {
        res = rv_info.rvCard;
	// it'll never get bigger than this so we 
	// might as well leave the loop
	break;
    } else if (rv_info.discImplementations[i] == CPT::di_MTCPT) {
	if (res == 0) res = 1;
    } else { // sparce CPT
	// we only get the sparse CPT cardinality if it has
	// been read in, otherwise, we just take card.
	unsigned tmp;
	if (i < conditionalCPTs.size())
	  tmp = conditionalCPTs[i]->maxCardinality();
	else 
	  tmp = rv_info.rvCard;
	if (tmp > res) res = tmp;
    }
  }
  if (res == 0) {
    // if there are no CPTs (i.e., if we didn't load in the
    // parameters during say a triangulation procedure, then we just
    // return the cardinality.
    return rv_info.rvCard;
  } else {
    // otherwise, we return the maximum.
    return res;
  }
}


