/*-
 * GMTK_EMable.cc
 *     A vector used for means of Gaussians.
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
#include "error.h"
#include "rand.h"
#include "logp.h"

#include "GMTK_EMable.h"
#include "GMTK_GMParms.h"

VCID("$Header$");

logpr
// EMable::minIncrementProbabilty = logpr((void*)NULL,log_FLT_MIN);
EMable::minIncrementProbabilty = logpr((void*)NULL,log_FLT_MIN);

unsigned long
EMable::missedIncrementCount = 0;


////////////////////////////////////////////////////
// The minimum accumulated probability of mean and covariance -like
// objects. If the accumulated probability falls below this
// value, then the mean or variance like object will not
// update its values.
logpr EMable::_minContAccumulatedProbability = 
EMable::setMinContAccumulatedProbability(logpr((void*)NULL, (double)-600.0));
logpr EMable::setMinContAccumulatedProbability(const logpr floor) 
{ 
  // We hard limit to be no less than exp(-700).
  // we use -700 because if it is much smaller than this
  // we won't be able to take the exp of the inverse (i.e.,
  // exp(700) is to large to represent in non-log space.
  // We need to do this, however, because continuous parameters
  // need to be updated not in log probability.
  if (floor.val() < -700)
    _minContAccumulatedProbability = logpr((void*)NULL, (double)-700.0);
  else 
    _minContAccumulatedProbability = floor; 
  return _minContAccumulatedProbability;
}

logpr EMable::_minDiscAccumulatedProbability =
EMable::setMinDiscAccumulatedProbability(logpr((void*)NULL, (double)LSMALL/2.0));

logpr EMable::setMinDiscAccumulatedProbability(const logpr floor) 
{ 
  // Hard limit to be no less than LSMALL, the smallest 
  // possible non-zero log probability.
  if (floor.val() < LSMALL)
    _minDiscAccumulatedProbability = logpr((void*)NULL, (double)LSMALL);
  else 
    _minDiscAccumulatedProbability = floor; 
  return _minDiscAccumulatedProbability;
}





////////////////////////////////////////////////////////////////////
//        General create, read, destroy routines 
////////////////////////////////////////////////////////////////////

