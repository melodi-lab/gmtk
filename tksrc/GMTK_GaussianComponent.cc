/*-
 * GMTK_GaussianComponent.cc
 *        Any of the common code for the family of Gaussian-like
 *        components classes.
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
VCID("$Header$")
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

// true if when a clone occurs, we use the same mean (i.e.,
// share the mean and only clone other things)
bool GaussianComponent::cloneShareMeans = false;
// true if when a clone occurs, we use the same covariance
// (i.e., share the covariance and clone other things)
bool GaussianComponent::cloneShareCovars = false;
// true if hwen a clone occurs, we use the same dlink matrix
bool GaussianComponent::cloneShareDlinks = false;

// l2-norm based accuracy regularization ceofficients for training
// regularized Gaussians. 
// The first coefficient applies only to the means of Gaussiasn. 
double GaussianComponent::gmarCoeffL2 = 0.0;
// The second coefficient applies only to the dlink matrices of Gaussians,
// so can be used to produce regularized BMMs or regularized covariance matrices
double GaussianComponent::gdarCoeffL2 = 0.0;



GaussianComponent::GaussianComponent(const int dim) 
  : Component(dim) 
{

}
