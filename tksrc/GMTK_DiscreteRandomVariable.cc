/*-
 * GMTK_DiscreteRandomVariable.cc
 *     Support code for discrete random variables.
 *
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

DiscreteRandomVariable::DiscreteRandomVariable(string _label,int card)
    : RandomVariable(_label, Discrete, card) {;}

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
//  printf("DiscreteRandomVariable::findConditionalParents called\n");
  cachedIntFromSwitchingState = intFromSwitchingState();
  assert ( cachedIntFromSwitchingState >= 0 && 
     cachedIntFromSwitchingState < conditionalCPTs.size()  &&
     cachedIntFromSwitchingState < conditionalParentsList.size());
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
  error("undefined");

}



/*-
 *-----------------------------------------------------------------------
 * tieParametersWith()
 *      Ties the parameters of 'this' with whatever those of 'other' are. 
 *      'other' and 'this' must be identical structuraly.
 * 
 * Preconditions:
 *      other must be a fully instantiated RV with parameters, and 'this'
 *      and 'other' must be structurally identical.
 *
 * Postconditions:
 *      'this' has the identical _tied_ parameters with 'other'
 *
 * Side Effects:
 *      Changes the internal parameter data structures, but does not delete anything.
 *
 * Results:
 *      returns nothing.
 *
 *-----------------------------------------------------------------------
 */
void
DiscreteRandomVariable::tieParametersWith(RandomVariable*const _other)
{
  assert ( _other -> discrete );
  
  DiscreteRandomVariable * other = (DiscreteRandomVariable*) _other;

  if (!identicalStructureWith(*other))
    error("Error, trying to tie parameters of RV '%s' with RV '%s' but they have different structure.",
	  label.c_str(),other->label.c_str());

  conditionalCPTs = other->conditionalCPTs;
  curCPT = other->curCPT;
}



/*-
 *-----------------------------------------------------------------------
 * clone()
 *      Returns a clone of self.
 * 
 * Preconditions:
 *      self must be filled in.
 *
 * Postconditions:
 *      same as preconditions.
 *
 * Side Effects:
 *      No internal effects.
 *
 * Results:
 *      returns a new random variable.
 *
 *-----------------------------------------------------------------------
 */
RandomVariable*
DiscreteRandomVariable::clone()
{
  DiscreteRandomVariable* rv = 
    (DiscreteRandomVariable*) RandomVariable::clone();
  // might as well set val, although probably won't be useful.
  rv->val = val;
  rv->tieParametersWith(this);
  rv->featureElement = featureElement;
  return rv;
}


