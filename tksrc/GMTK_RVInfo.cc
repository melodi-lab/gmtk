/*-
 * GMTK_RVInfo.cc
 *     RV generic information
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2007, < fill in later >
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

#include "GMTK_RVInfo.h"


#include "GMTK_RV.h"

VCID("$Header$")


/*-
 *-----------------------------------------------------------------------
 * RVInfo::clear()
 *   clear out the current RV structure to a generic known state, and 
 *   also clear up any memory used by this object. Useful for construction/destruction,
 *   or when we are parsing and encounter a new RV.
 * 
 * Preconditions:
 *      None
 *
 * Postconditions:
 *      None
 *
 * Side Effects:
 *      changes internal member variables.
 *
 * Results:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
void 
RVInfo::clear() {
    name.erase();

    rvType = t_unknown;
    rvDisp = d_unknown;
    rvFeatureRange.clear();
    rvWeightInfo.clear();
    eliminationOrderHint = 0.0;
    variablePositionInStrFile = -1;

    switchingParents.clear();
    switchMapping.clear();

    conditionalParents.clear();
    discImplementations.clear();
    contImplementations.clear();
    listIndices.clear();

    isDeterministic = isSparse = false;

    delete rv;
    rv = NULL;
}
