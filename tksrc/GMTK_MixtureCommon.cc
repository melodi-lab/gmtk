/*-
 * GMTK_MixtureCommon.cc
 *        Any of the common code for the family of mixture-like
 *        classes.
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (C) 2001 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
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
#if HAVE_CONFIG_H
#include <config.h>
#endif
#if HAVE_HG_H
#include "hgstamp.h"
#endif
VCID(HGID)

#include "error.h"
#include "rand.h"


#include "GMTK_MixtureCommon.h"
#include "GMTK_Dense1DPMF.h"
#include "GMTK_MeanVector.h"
#include "GMTK_DiagCovarVector.h"
#include "GMTK_Component.h"

//////////////////////////////////////////////////////////////////
// set the mcvr. By default it is set to a large
// value meaning that mixtures are not removed.
// A reasonable value, to eagerly start removing
// mixtures, is about 50.0 or so.
double
MixtureCommon::mixCoeffVanishRatio = 1e20;

//////////////////////////////////////////////////////////////////
// set the mcsr. By default it is set to a large
// value meaning that mixtures are not split.
double
MixtureCommon::mixCoeffSplitRatio = 1e10;

unsigned
MixtureCommon::numTopToForceSplit = 0;

unsigned 
MixtureCommon::numBottomToForceVanish = 0;

bool
MixtureCommon::cacheComponentsInEmTraining = true;

bool
MixtureCommon::cacheMixtureProbabilities = true;

void
MixtureCommon::checkForValidRatioValues() {
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
set<pair<Dense1DPMF*,unsigned> > MixtureCommon::vanishingComponentSet;
set<pair<Dense1DPMF*,unsigned> > MixtureCommon::splittingComponentSet;

map<MeanVector*,MeanVector*> MixtureCommon::meanCloneMap;
map<DiagCovarVector*,DiagCovarVector*> MixtureCommon::diagCovarCloneMap;
map<DlinkMatrix*,DlinkMatrix*> MixtureCommon::dLinkMatCloneMap;
map<Component*,Component*> MixtureCommon::mcCloneMap;

