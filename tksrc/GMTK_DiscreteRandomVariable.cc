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
 * Results:
 *      none
 *
 * Side Effects:
 *      Changes some internal object structures such as;
 *
 *-----------------------------------------------------------------------
 */
void
DiscreteRandomVariable::findConditionalParents()
{
  cachedIntFromSwitchingState = intFromSwitchingState();
  assert (cachedIntFromSwitchingState >= 0 && 
	  cachedIntFromSwitchingState < conditionalCPTs.size());
  curConditionalParents = & conditionalParentsList[cachedIntFromSwitchingState];
  curCPT = conditionalCPTs[cachedIntFromSwitchingState];
}



/*-
 *-----------------------------------------------------------------------
 * Function
 *     tieWith(RandomVariable *rv)
 *     Sets things up so that this random variable has all parameters
 *     that are tied with the argument.
 *  
 * Results:
 *      none
 *
 * Side Effects:
 *      Has a significant effect on the internal structures.
 *
 *-----------------------------------------------------------------------
 */
void
DiscreteRandomVariable::tieWith(RandomVariable* rv)
{

#ifndef NDEBUG
  // first do massive amounts of checking to make sure
  // that everything is kosher.
  assert ( rv->discrete );
  assert ( cardinality == rv->cardinality );
  assert ( switchingParents.size() == rv->switchingParents.size() );

#endif

  // everything checks out, now set our CPS to have same pointers
  // as rv has. Assume it is discrete and cast.
  conditionalCPTs = ((DiscreteRandomVariable*)rv)->conditionalCPTs;

  // TODO: finish this function.
  error("function tieWith not finished");

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
