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
      error("DiagCovarVector:: read, covariance[%d] = (%e) < current Floor = (%e)",
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







