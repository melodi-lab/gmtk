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

#include "GMTK_DiscreteRandomVariable.cc"

VCID("$Header$");


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
DiscreteRandomVariable::findConditionalParents()
{
  cachedIntFromSwitchingState = intFromSwitchingState();
  assert (cachedIntFromSwitchingState >= 0 && 
	  cachedIntFromSwitchingState < conditionalCPTs.len());
  curConditionalParents = & conditionalParentsList[cachedIntFromSwitchingState];
  curCPT = conditionalCPTs[cachedIntFromSwitchingState];
}



/*-
 *-----------------------------------------------------------------------
 * Function
 *     tieWith(randomVariable *rv)
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
DiscreteRandomVariable::tieWith(randomVariable* rv)
{

#ifndef NDEBUG
  // first do massive amounts of checking to make sure
  // that everything is kosher.
  assert ( rv->discrete );
  assert ( cardinality == rv->cardinality );
  assert ( switchingParents.len() == rv->switchingParents().len() );
  for (int i=0;i<switchingParents.len();i++) {
    assert ( switchingParents[i]->cardinality == 
	      rv->switchingParents[i]->cardinality );
  }
  assert ( conditionalParentsList.len() == rv->conditionalParentsList.len() );
  for (int i=0;i<conditionalParentsList.len();i++) {
    assert ( conditionalParentsList[i] == 
	     rv->conditionalParentsList[i] );
    for (int j=0;j<conditionalParentsList[i].len();j++) {
      assert (
	      conditionalParentsList[i].[j].len() ==
	      rv->conditionalParentsList[i].[j].len()
	      );
    }
  }
#endif

  // everything checks out, now set our CPS to have same as others.
  for (int i = 0; i < conditionalCPTs.len(); i++) {
    conditionalCPTs[i] = rv.conditionalCPTs[i];
  }

}

