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

#define DEFINE_DETERMINISTIC_MAPPER_MACROS 1
#include "GMTK_CFunctionDeterministicMappings.h"
#include "GMTK_RngDecisionTree.h"
#include "GMTK_DiscRV.h"
#include "GMTK_GMParms.h"





//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// GMTK Internal C function deterministic mapping functions.
// DO NOT MODIFY ANYTHING IN THE NEXT BIT OF CODE STARTING HERE.
//
#define COPYPARENT_NUM_FEATURES 1
DEFINE_DETERMINISTIC_MAPPER_C_CODE(copyParent,COPYPARENT_NUM_FEATURES)
{
  DiscRVType rv = p0; 
  return rv;
}
//
#define ALWAYSZERO_NUM_FEATURES 0
DEFINE_DETERMINISTIC_MAPPER_C_CODE(alwaysZero,ALWAYSZERO_NUM_FEATURES)
{
  return (DiscRVType) 0;
}
//
#define ALWAYSONE_NUM_FEATURES 0
DEFINE_DETERMINISTIC_MAPPER_C_CODE(alwaysOne,ALWAYSONE_NUM_FEATURES)
{
  return (DiscRVType) 1;
}
//
#define INCREMENT_NUM_FEATURES 1
DEFINE_DETERMINISTIC_MAPPER_C_CODE(increment,INCREMENT_NUM_FEATURES)
{
  return (p0+1);
}
//
#define DECREMENT_NUM_FEATURES 1
DEFINE_DETERMINISTIC_MAPPER_C_CODE(decrement,DECREMENT_NUM_FEATURES)
{
  return (p0-1);
}
//
#define CONDITIONAL_INCREMENT_NUM_FEATURES 2
DEFINE_DETERMINISTIC_MAPPER_C_CODE(conditionalIncrement,CONDITIONAL_INCREMENT_NUM_FEATURES)
{
  // increments p0 if p1 is non-zero, otherwise returns p0.
  if (p1)
    return (p0+1);
  else
    return (p0);
}
//
#define CONDITIONAL_DECREMENT_NUM_FEATURES 2
DEFINE_DETERMINISTIC_MAPPER_C_CODE(conditionalDecrement,CONDITIONAL_DECREMENT_NUM_FEATURES)
{
  // decrements p0 if p1 is non-zero, otherwise returns p0.
  if (p1)
    return (p0-1);
  else
    return (p0);
}
//
#define CONDITIONAL_LIMITED_INCREMENT_NUM_FEATURES 3
DEFINE_DETERMINISTIC_MAPPER_C_CODE(conditionalLimitedIncrement,CONDITIONAL_LIMITED_INCREMENT_NUM_FEATURES)
{
  // increments p0 if p1 is non-zero. Increment up to and including value given by p2
  // but not beyond. Otherwise returns p0.
  if (p1 && p0 < p2)
    return (p0+1);
  else
    return (p0);
}
//
#define CONDITIONAL_LIMITED_DECREMENT_NUM_FEATURES 3
DEFINE_DETERMINISTIC_MAPPER_C_CODE(conditionalLimitedDecrement,CONDITIONAL_LIMITED_DECREMENT_NUM_FEATURES)
{
  // decrements p0 if p1 is non-zero. Decrement down to and including value given by p2
  // but not below. Otherwise returns p0.
  if (p1 && p0 > p2)
    return (p0-1);
  else
    return (p0);
}


//
// The below includes routines with a variable number of parents.
// TODO: get the below working, with variable numbers of parents, we
// need to change DT code to allow variable num features.
//
DEFINE_DETERMINISTIC_MAPPER_C_CODE(allParentsEqual,CDT_VARIABLE_NUMBER_FEATURES)
{
  // returns 1 if all parents are equal, and otherwise returns zero.
  // Note that this works for any number of parents.
  if ( numParents == 0 ) {
    error("CDT allParentsEqual called with zero features. Need to have at least 1.");
  }
  DiscRVType rv = p0; 
  for (unsigned i = 1 ; i < numParents ; i++ ) {
    if (par(i) != rv)
      return (DiscRVType) 0;
  }
  return (DiscRVType) 1;
}
//
DEFINE_DETERMINISTIC_MAPPER_C_CODE(allParentsUnEqual,CDT_VARIABLE_NUMBER_FEATURES)
{
  // returns 1 if all parents are *un*equal, and otherwise returns zero.
  // Note that this works for any number of parents.
  if ( numParents == 0 ) {
    error("CDT allParentsUnEqual called with zero features. Need to have at least 1.");
  }
  // return 0 if we find any two parents that are equal. 
  // TODO: Is there a faster way than O(N^2) to do this?
  for (unsigned i = 0 ; i < numParents ; i++ ) {
    for (unsigned j = i+1 ; j < numParents ; j++ ) {
      if (par(i) == par(j))
	return (DiscRVType) 0;
    }
  }
  return (DiscRVType) 1;
}
// DO NOT MODIFY ANYTHING IN THE ABOVE BIT OF CODE ENDING HERE.
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

// note: when a, b are integers, a/b = floor(a/b).
// To get ceil and round, using only integer ops, we have that:
//    ceil((float)a/(float)b) = (a + a - 1)/b
//    round((float)a/(float)b) = (a + a/2)/b = (a + (a>>1))/b.
// One is free here to convert to/from floating point, but often one need not do that.

/*
 * Here is an example mapping that one might want to define.
 * Uncomment to activate.
#define AJIT_MAPPING_NUM_FEATURES 10
DEFINE_DETERMINISTIC_MAPPER_C_CODE(ajitMapping,AJIT_MAPPING_NUM_FEATURES)
{ 
  const DiscRVType SHIFT_ZERO=0;

  return (p0 > 0 ? ((  ( (p1+1 - 17*p7 - 18*p8 + p9) + (p0>>1))/p0) + 1*(p6-SHIFT_ZERO) >= 3 ? ( ( ((p1 + 1 - 17*p7 - 18*p8 + p9) + (p0>>1))/p0) + 1*(p6-SHIFT_ZERO) < 4 ? 1 : 0) : 0) : 0) 
    ||
    (p3 > 0 ? (( ((p4+19 - 17*p7 - 18*p8 + p9) + (p3>>1))/p3) + 1*(p6-SHIFT_ZERO) >= 3 ? ( ( ((p4+19 - 17*p7 - 18*p8 + p9) + (p3>>1))/p3 + 1*(p6-SHIFT_ZERO)) < 4 ? 1 : 0) : 0) : 0);
}
*/

//
// add more deterministic mapping functions here as desired ...
// 
#ifdef USER_INTERNAL_CFUNC_DTS
#include "user_internal_cfunc_dts.cc"
#endif

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// Registration code. This registers all deterministic mapping
// functions that were defined above. Once  they are registered, with a given
// name, they may be used like any decision tree.
//
// NOTE: please give all decision trees registered as such a name
// starting with "internal:" and any user defined functions
// a name starting with "user_internal:"



void
registerAllCFunctionDeterministicMappings(GMParms& gmp)
{

  ///////////////////////////////////////////////////////////////////////
  // DO NOT CHANGE ANYTHING IN THE FOLLOWING FEW LINES STARTING HERE.
  gmp.registerDeterministicCMapper("internal:copyParent",
				   COPYPARENT_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(copyParent));
  gmp.registerDeterministicCMapper("internal:alwaysZero",
				   ALWAYSZERO_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(alwaysZero));
  gmp.registerDeterministicCMapper("internal:alwaysOne",
				   ALWAYSONE_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(alwaysOne));
  gmp.registerDeterministicCMapper("internal:increment",
				   INCREMENT_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(increment));
  gmp.registerDeterministicCMapper("internal:decrement",
				   DECREMENT_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(decrement));
  gmp.registerDeterministicCMapper("internal:conditionalIncrement",
				   CONDITIONAL_INCREMENT_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(conditionalIncrement));
  gmp.registerDeterministicCMapper("internal:conditionalDecrement",
				   CONDITIONAL_DECREMENT_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(conditionalDecrement));
  gmp.registerDeterministicCMapper("internal:conditionalLimitedIncrement",
				   CONDITIONAL_LIMITED_INCREMENT_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(conditionalLimitedIncrement));
  gmp.registerDeterministicCMapper("internal:conditionalLimitedDecrement",
				   CONDITIONAL_LIMITED_DECREMENT_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(conditionalLimitedDecrement));
  // TODO: add the variable parent deterministic mappers once the DT code can accept variable numbers of parens.
  // ...
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

  /*
   * Here is an example mapping registration that one might want to define.
   * Uncomment to activate.
  gmp.registerDeterministicCMapper("user_internal:ajitMapping",
				   AJIT_MAPPING_NUM_FEATURES,
				   DETERMINISTIC_MAPPER_C_CODE_NAME(ajitMapping));
  */

  // Uncomment to register user defined DTs. You can defiine them in the below file if you like.
  // 

#ifdef USER_INTERNAL_CFUNC_DTS
#include "register_user_internal_cfunc_dts.cc"
#endif


}
