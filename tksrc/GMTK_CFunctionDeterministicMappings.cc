/*-
 * GMTK_CFunctionDecisionTrees.cc
 *     Code that defines and declares internal pre-defined GMTK deterministic
 *     mapping functions, that map from a set of random variables (nominally a
 *     set of parent variables and one child variable) that are non-negative integer
 *     valued down to a single integer value.
 *     This file serves two purposes:
 *         1) to define and register internal oft-used determinstic
 *            mapping functions that are useful in a variety of
 *            contexts (such as copy parent, etc.).  
 *         2) to allow the user to define their own mapping function. This is useful
 *            when a decision tree is large, complicated, and most importantly slow.
 *            By being able to define it here, at compile time, we can take advantage
 *            of the optimizing C++ comiler producing an efficient implementation of
 *            the desired formula (not to mention that now, not only formulas but
 *            one can also use loops, subroutine calls, local variables, floating
 *            point, and so on).  
 *            In this case, the user must make sure that their function doesn't have
 *            infinite loops, and so on, as otherwise inference will just stall.
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

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <map>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#include "error.h"
#include "general.h"
#include "rand.h"
#include "sArray.h"

#include "GMTK_DiscRV.h"
#include "GMTK_HidDiscRV.h"
#include "GMTK_ObsDiscRV.h"
#include "GMTK_PackCliqueValue.h"
#include "GMTK_RngDecisionTree.h"



// Internal functions, do not modify anything here.


void
cFunctionDeterministicMapping_copyparent(
        const vector< RV* >& parent_variables,
	const RV* const child_rv)
{
  DiscRVType foo;
  return foo + p0;
}

DEFINE_FUNCTION_MAPPER_C_CODE(copysingleparent)
{
  return p0;
}



// Additional user defined DTs. Uncomment, change to the name you want
// and then register this below.




// DT registraiton code.


void
registerAllMappers()
{

  // DO Not change anything here.
  registerDeterministicCMapper("internal:copyparent",1,cFunctionDeterministicMapping_copyparent);
  registerDeterministicCMapper("internal:incrementIfTrue",5,cFunctionDeterministicMapping_increment);

  // Uncomment to register user defined DTs.
  // 

}


























