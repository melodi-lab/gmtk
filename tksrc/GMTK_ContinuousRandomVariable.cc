/*-
 * GMTK_ContinuousRandomVariable.cc
 *     Support code for continuous random variables.

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

#include "GMTK_ContinuousRandomVariable.h"

VCID("$Header$");

ContinuousRandomVariable::ContinuousRandomVariable(string _label)
  : RandomVariable(_label,Continuous) 
{


}

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
ContinuousRandomVariable::findConditionalParents()
{
  cachedIntFromSwitchingState = intFromSwitchingState();
  assert ( cachedIntFromSwitchingState >= 0 && 
     cachedIntFromSwitchingState < conditionalParentsList.size() );
  curConditionalParents = & conditionalParentsList[cachedIntFromSwitchingState];
}

