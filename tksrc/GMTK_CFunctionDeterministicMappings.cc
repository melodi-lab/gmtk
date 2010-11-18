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
 *    How to use this file?
 *         1) Search for the tag REGISTRATION_OF_USER_DEFINED_DETERMINISTIC_MAPPERS:
 *            Here you define a bunch of C-like functions that become deterministic mappers.
 *            You can define any C subroutines in this region of the code as well.
 *         2) Next, search for the tag REGISTRATION_OF_USER_DEFINED_DETERMINISTIC_MAPPERS:
 *            below and copy/uncomment necessary code to register your newly defined 
 *            deterministic mappers.
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
#include "GMTK_DiscRV.h"
#incldue "GMTK_GMParams.h"

#define DEFINE_DETERMINISTIC_MAPPER_MACROS
#include "GMTK_CFunctionDecisionTrees.h"


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// GMTK Internal C function deterministic mapping functions.
// DO NOT MODIFY ANYTHING IN THE NEXT BIT OF CODE STARTING HERE.

#define COPYPARENT_NUM_FEATURES 1
DEFINE_DETERMINISTIC_MAPPER_C_CODE(copyParent,COPYPARENT_NUM_FEATURES)
{
  DiscRVType rv = p0; 
  return rv;
}

// DO NOT MODIFY ANYTHING IN THE NEXT BIT OF CODE ENDING HERE.
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// DEFINITION_OF_USER_DEFINED_DETERMINISTIC_MAPPERS: 
// Additional user defined DTs. A few examples are given, you can uncomment and
// modify at will, and then recompile GMTK and these deterministic functions
// will be available to you to use just like any decision tree-based deterministic
// mapping.


#define USERDTBINMAPPING_NUM_FEATURES 10
DEFINE_DETERMINISTIC_MAPPER_C_CODE(userDTbinMapping,USERDTBINMAPPING_NUM_FEATURES)
{ 
  const DiscRVType SHIFT_ZERO=0;
  return (p0 > 0 ? (((p1+1 - 17*p7 - 18*p8 + p9)~/p0) + 1*(p6-SHIFT_ZERO) >= 3 ? ( ((p1 + 1 - 17*p7 - 18*p8 + p9)~/p0) + 1*(p6-SHIFT_ZERO) < 4 ? 1 : 0) : 0) : 0) 
    ||
    (p3 > 0 ? (((p4+19 - 17*p7 - 18*p8 + p9)~/p3) + 1*(p6-SHIFT_ZERO) >= 3 ? ( ((p4+19 - 17*p7 - 18*p8 + p9)~/p3 + 1*(p6-SHIFT_ZERO)) < 4 ? 1 : 0) : 0) : 0);
}



//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// Registration code.

void
registerAllCFunctionDeterministicMappings()
{

  ///////////////////////////////////////////////////////////////////////
  // DO NOT CHANGE ANYTHING IN THE FOLLOWING FEW LINES STARTING HERE.
  registerDeterministicCMapper("internal:copyParent",
			       COPYPARENT_NUM_FEATURES,
			       DETERMINISTIC_MAPPER_C_CODE_NAME(copyparent));
  // DO NOT CHANGE ANYTHING IN THE ABOVE FEW LINES ENDING HERE.
  ///////////////////////////////////////////////////////////////////////




  ///////////////////////////////////////////////////////////////////////
  // REGISTRATION_OF_USER_DEFINED_DETERMINISTIC_MAPPERS: 
  // ADD USER DEFINED C FUNCTION DETERMINISTIC MAPPING REGISTRATIONS HERE.
  // Arguments are:
  //   registerDeterministicCMapper(
  //          name_of_deterministic_mapping which is of type char*,
  //          number of features (just like when defining a decision tree),
  //          C function above, use macro
  //        );

  registerDeterministicCMapper("cmapper:userDTbinMapping",
			       USERDTBINMAPPING_NUM_FEATURES,
			       DETERMINISTIC_MAPPER_C_CODE_NAME(userDTbinMapping));

       "internal:copyParent",1,cFunctionDeterministicMapping_copyparent);

  // Uncomment to register user defined DTs.
  // 

}


























