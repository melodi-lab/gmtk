/*-
 * GMTK_CPT.cc
 *     Trainable (with say EM) CPT
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

#include "GMTK_CPT.h"


VCID("$Header$");

////////////////////////////////////////////////////////////////////
//        Static Data
////////////////////////////////////////////////////////////////////


int CPT::warningNumParents = 50;

////////////////////////////////////////////////////////////////////
//        General create, read, destroy routines 
////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////
//        Misc Support
////////////////////////////////////////////////////////////////////




/*-
 *-----------------------------------------------------------------------
 * Function
 *      compareCardinalities: compare the cardinalities of this CPT with that of an other. REturn
 *      true if they are equil false otherwise.
 *
 * Results:
 *      returns true if cards are equal.
 *
 * Side Effects:
 *      none
 *
 *-----------------------------------------------------------------------
 */
bool 
CPT::compareCardinalities(CPT& cpt)
{
  if (cardinalities.len() != cpt.cardinalities.len())
    return false;

  for (int i=0;i<cardinalities.len();i++) {
    if (cardinalities[i] != cpt.cardinalities[i])
      return false;
  }
}


