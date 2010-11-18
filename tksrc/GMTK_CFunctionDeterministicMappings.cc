/*-
 * GMTK_CFunctionDecisionTrees.cc
 *
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
 * Copyright (c) 2010, < fill in later >
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

#include "GMTK_RngDecisionTree.h"
#incldue "GMTK_GMParams.h"

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// GMTK Internal C function deterministic mapping functions.
// DO NOT MODIFY ANYTHING HERE OR YOU WILL VOID YOUR WARRANTY!!

DiscRVType 
cFunctionDeterministicMapping_copyparent(
        const vector< RV* >& parent_variables,
	const RV* const child_rv)
{
  DiscRVType foo;
  return foo + p0;
}

DEFINE_DETERMINISTIC_MAPPER_C_CODE(copyParent,1)
{
  DiscRVType rv = p0; 
  return rv;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// Additional user defined DTs. A few examples are given, you can uncomment and
// modify at will, and then recompile GMTK and these deterministic functions
// will be available to you to use just like any decision tree-based deterministic
// mapping.





//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// Registration code.

void
registerAllCFunctionDeterministicMappings()
{

  ///////////////////////////////////////////////////////////////////////
  // DO NOT CHANGE ANYTHING IN THE FOLLOWING FEW LINES
  registerDeterministicCMapper("internal:copyParent",
			       1,
			       DETERMINISTIC_MAPPER_C_CODE_NAME(copyparent));

  registerDeterministicCMapper("internal:copyParent",
			       1,
			       DETERMINISTIC_MAPPER_C_CODE_NAME(copyparent));

  ///////////////////////////////////////////////////////////////////////



  ///////////////////////////////////////////////////////////////////////
  // ADD USER DEFINED C FUNCTION DETERMINISTIC MAPPING REGISTRATIONS HERE.
  // Arguments are:
  //   registerDeterministicCMapper(
  //          name_of_deterministic_mapping which is of type char*,
  //          number of features (just like when defining a decision tree),
  //          C function above, use macro 

       "internal:copyParent",1,cFunctionDeterministicMapping_copyparent);

  // Uncomment to register user defined DTs.
  // 

}


























