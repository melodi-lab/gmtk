/*-
 * GMTK_Mixgaussiancommon.cc
 *        Any of the common code for the family of Gaussian-like
 *        classes.
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2001, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle, and Jeff Bilmes make no representations about
 * the suitability of this software for any purpose.  It is provided
 * "as is" without express or implied warranty.
 *
 */


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <float.h>
#include <assert.h>

#include "general.h"
VCID("$Header$");
#include "error.h"
#include "rand.h"


#include "GMTK_MixGaussiansCommon.h"
#include "GMTK_Dense1DPMF.h"
#include "GMTK_MeanVector.h"
#include "GMTK_DiagCovarVector.h"
#include "GMTK_GaussianComponent.h"

//////////////////////////////////////////////////////////////////
// set the mcvr. By default it is set to a large
// value meaning that mixtures are not removed.
// A reasonable value, to eagerly start removing
// mixtures, is about 50.0 or so.
double
MixGaussiansCommon::mixCoeffVanishRatio = 1e20;

//////////////////////////////////////////////////////////////////
// set the mcsr. By default it is set to a large
// value meaning that mixtures are not split.
double
MixGaussiansCommon::mixCoeffSplitRatio = 1e10;


void
MixGaussiansCommon::checkForValidRatioValues() {
  // this next check guarantees that we will never eliminate
  // all components
  if (mixCoeffVanishRatio < 1.0)
    error("ERROR: must have mixCoeffVanishRatio >= 1.0");
  if (mixCoeffSplitRatio <= 0)
    error("ERROR: must have mixCoeffSplitRatio > 0.0");
  if (1.0/mixCoeffVanishRatio >= mixCoeffSplitRatio) 
    error("ERROR: must have 1.0/mixCoeffVanishRatio < mixCoeffSplitRatio");
}


//////////////////////////////////////////////////////////////////
// support for component vanishing & splitting with multiple
// potentially shared objects.
set<pair<Dense1DPMF*,unsigned> > MixGaussiansCommon::vanishingComponentSet;
set<pair<Dense1DPMF*,unsigned> > MixGaussiansCommon::splittingComponentSet;

map<MeanVector*,MeanVector*> MixGaussiansCommon::meanCloneMap;
map<DiagCovarVector*,DiagCovarVector*> MixGaussiansCommon::diagCovarCloneMap;
map<DlinkMatrix*,DlinkMatrix*> MixGaussiansCommon::dLinkMatCloneMap;
map<GaussianComponent*,GaussianComponent*> MixGaussiansCommon::gcCloneMap;

