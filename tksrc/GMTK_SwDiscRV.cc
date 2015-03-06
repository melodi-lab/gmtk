/*
 * GMTK_SwDiscRV.cc
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


#include "general.h"
#if HAVE_CONFIG_H
#include <config.h>
#endif
#if HAVE_HG_H
#include "hgstamp.h"
#endif
VCID(HGID)


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


