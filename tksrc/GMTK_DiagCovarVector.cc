/*-
 * GMTK_DiagCovarVector.cc
 *     Trainable shared diagonal covariance matrix (vector)
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
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <ieeefp.h>
#include <float.h>
#include <assert.h>

#include "general.h"
#include "error.h"
#include "rand.h"

#include "GMTK_DiagCovarVector.h"
#include "GMTK_GaussianComponent.h"

#ifndef M_PI
#define M_PI               3.14159265358979323846  /* pi */
#endif



VCID("$Header$");


////////////////////////////////////////////////////////////////////
//        General create, read, destroy routines 
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * DiagCovarVector::DiagCovarVector()
 *      Constructor
 *
 * Results:
 *      Constructs the object.
 *
 * Side Effects:
 *      None so far.
 *
 *-----------------------------------------------------------------------
 */
DiagCovarVector::DiagCovarVector() 
{
}


/*-
 *-----------------------------------------------------------------------
 * DiagCovarVector::read()
 *      read in, and make sure it is a valid cov. matrix.
 *
 * Results:
 *      Constructs the object.
 *
 * Side Effects:
 *      None so far.
 *
 *-----------------------------------------------------------------------
 */
void 
DiagCovarVector::read(iDataStreamFile& is)
{
  NamedObject::read(is);
  covariances.read(is); 
  for (int i=0;i<covariances.len();i++) {
    if (covariances[i] < GaussianComponent::varianceFloor()) {
      error("Error: reading diagonal covariance matrix '%s' (file '%s'), but covariance[%d] = (%e) < current Floor = (%e)",
	    name().c_str(),is.fileName(),
	    i,covariances[i],GaussianComponent::varianceFloor());
    }
  }
}





////////////////////////////////////////////////////////////////////
//        Misc Support
////////////////////////////////////////////////////////////////////

void
DiagCovarVector::makeRandom()
{
  for (int i=0;i<covariances.len();i++) {
    covariances[i] = 1.0+rnd.drand48();
  }
}

void
DiagCovarVector::makeUniform()
{
  for (int i=0;i<covariances.len();i++) {
    covariances[i] = 1.0;
  }
}


void
DiagCovarVector::preCompute()
{
  variances_inv.growIfNeeded(covariances.len());
  double det = 1.0;
  for (int i=0;i<covariances.len();i++) {
    if (covariances[i] <= GaussianComponent::varianceFloor()) {
      // Theoretically, this shouldn't happen unless you are reading
      // in a file that was computed from a previous run with a different threshold.
      coredump("ERROR: element %i of diagonal covariance matrix '%s' is at or below floor value",i,name().c_str());
    }
    variances_inv[i] = 1.0/covariances[i];
    det *= covariances[i];
  }
  if (det <= DBL_MIN) {
    coredump("ERROR: determinant of diagonal covariance matrix '%s' has hit minimum",name().c_str());
  }
  const double tmp = (::pow(2*M_PI,covariances.len()/2.0)*::sqrt(det));
  if (tmp <= DBL_MIN)
    coredump("ERROR:  norm const has hit maximum of diagonal covariance matrix '%s'",name().c_str());
  _log_inv_normConst = -0.5*(covariances.len()*::log(2*M_PI) + ::log(det));
}









