/*-
 * GMTK_DiscreteRandomVariable.cc
 *     Support code for discrete random variables.

 * Written by Jeff Bilmes<bilmes@ee.washington.edu>
 *
 * Copyright (c) 2001, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about the suitability of this software 
 * for any purpose. It is provided "as is" without express or implied warranty.
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

#include "GMTK_DiscreteRandomVariable.h"

VCID("$Header$");

DiscreteRandomVariable::DiscreteRandomVariable(string _label, vartype vt, 
int card)
    : RandomVariable(_label, vt, card) {;}

/*-
 *-----------------------------------------------------------------------
 * Function
 *      findConditionalParents()
 *     Set up conditional parents pointers and other tables.
 * 
 * Preconditions:
 *      variable must be filled in.
 *
 * Postconditions:
 *      What is true after the function is called.
 *
 * Side Effects:
 *      Changes some internal object structures such as;
 *
 * Results:
 *      What does the function return, if anything. 
 *
 *-----------------------------------------------------------------------
 */
void
DiscreteRandomVariable::findConditionalParents()
{
  cachedIntFromSwitchingState = intFromSwitchingState();
  assert ( cachedIntFromSwitchingState >= 0 && 
     cachedIntFromSwitchingState < conditionalCPTs.size() );
  curConditionalParents = & conditionalParentsList[cachedIntFromSwitchingState];
  curCPT = conditionalCPTs[cachedIntFromSwitchingState];
}


/*-
 *-----------------------------------------------------------------------
 * allocateProbabiltyTables()
 *      Allocate the internal CPT probability tables.
 * 
 * Preconditions:
 *      CPT structures shouldn't be defined. All parents
 *      both switching and otherwise should be defined.
 *
 * Postconditions:
 *      CPT structures are defined, and match the cardinalities
 *      of the current set of parents.
 *
 * Side Effects:
 *      Changes the interal CPT data structures.
 *
 * Results:
 *      What does the function return, if anything. 
 *
 *-----------------------------------------------------------------------
 */
void
DiscreteRandomVariable::allocateProbabiltyTables()
{
  

}
