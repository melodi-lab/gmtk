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


////////////////////////////////////////////////////////////////////
//        General create, read, destroy routines 
////////////////////////////////////////////////////////////////////

