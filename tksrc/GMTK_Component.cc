/*-
 * GMTK_GaussianComponent.cc
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
#include <ieeefp.h>
#include <float.h>
#include <assert.h>

#include "general.h"
VCID("$Header$");
#include "error.h"

#include "GMTK_GaussianComponent.h"
#include "rand.h"

////////////////////////////////////////////////////
// The default value of the minimum possible variance of any
// Gaussian. This must be >= FLT_MIN for numeric stability,
// and it should be made availalbe as a command line parameter at some point.
double GaussianComponent::_varianceFloor = GaussianComponent::setVarianceFloor(1e-15);

double GaussianComponent::setVarianceFloor(const double floor) 
{ 
  if (floor < FLT_MIN) 
    _varianceFloor = FLT_MIN; 
  else 
    _varianceFloor = floor; 
  return _varianceFloor;
}


////////////////////////////////////////////////////
// The minimum accumulated probability of mean and covariance -like
// objects. If the accumulated probability falls below this
// value, then the mean or variance like object will not
// update its values.
logpr GaussianComponent::_minAccumulatedProbability = 
GaussianComponent::setMinAccumulatedProbability(logpr((void*)NULL, (double)-600.0));


logpr GaussianComponent::setMinAccumulatedProbability(const logpr floor) 
{ 
  // hard limit to be no less than exp(-700).
  if (floor.val() < -700)
    _minAccumulatedProbability = logpr((void*)NULL, (double)-700.0);
  else 
    _minAccumulatedProbability = floor; 
  return _minAccumulatedProbability;
}


GaussianComponent::GaussianComponent(const int dim) : _dim(dim) 
{
}
