/*-
 * GMTK_GaussianCommon.cc
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

#include "GMTK_GaussianCommon.h"
#include "rand.h"

////////////////////////////////////////////////////
// The default value of the minimum possible variance of any
// Gaussian. This must be >= FLT_MIN for numeric stability.
double GaussianCommon::_varianceFloor = 1e-30; 

double GaussianCommon::varianceFloor() 
{ return _varianceFloor; }

void GaussianCommon::setVarianceFloor(const double floor) 
{ 
  if (floor < FLT_MIN) 
    _varianceFloor = FLT_MIN; 
  else 
    _varianceFloor = floor; 
}




