/*-
 * GMTK_FileParser.cc
 *     structure file parsing
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

#include "GMTK_FileParser.h"
#include "GMTK_RandomVariable.h"
#include "GMTK_DiscreteRandomVariable.h"
#include "GMTK_ContinuousRandomVariable.h"
#include "GMTK_GM.h"
#include "GMTK_GMTemplate.h"
#include "GMTK_GMParms.h"
#include "GMTK_MDCPT.h"
#include "GMTK_MSCPT.h"
#include "GMTK_MTCPT.h"
#include "GMTK_Mixture.h"
#include "GMTK_ObservationMatrix.h"

VCID("$Header$");

////////////////////////////////////////////////////////////////////
//        General create, read, destroy routines 
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * GraphicalModel::topologicalSortRecurse
 *      Support routine for topological Sort routine.
 *
 * Preconditions:
 *      must only be called from topological sort.
 *
 * Postconditions:
 *
 * Side Effects:
 *     Position variable is modified
 *
 * Results:
 *     returns true if everything works, return false if the
 *     graph has a directed loop.
 *
 *-----------------------------------------------------------------------
 */
bool
GraphicalModel::topologicalSortRecurse(vector<RandomVariable*>& outputVarList,
				       RandomVariable* node,
				       unsigned& position)
{
  node->tag = 1;
  for (unsigned i=0;i<node->allPossibleChildren.size();i++) {
    RandomVariable*rv = node->allPossibleChildren[i];
    if (rv->tag == 0) {
      bool res = topologicalSortRecurse(outputVarList,rv,position);
      if (!res)
	// directed graph has a loop
	return false;
    } else if (rv->tag == 1)
      // directed graph has a loop
      return false;
    else
      ;
      // tag == 2, meaning we've gone down this path before and need not
      // do it again.
  }
  node->tag = 2; // done with this node
  outputVarList[--position] = node;
  return true;
}



/*-
 *-----------------------------------------------------------------------
 * GraphicalModel::topologicalSort
 *      performes a topological sort of the set of random variables
 *      that live in inputVarList, and places the result of the
 *      sort in outputVarList. If the graph has a directed loop,
 *      fail and return false otherwise return true.
 *
 * Preconditions:
 *     inputVarLIst must contain a list of random variables with
 *     the variables appropriately set up.
 *
 * Postconditions:
 *     ouputVarList contains the list of variables in order.
 *
 * Side Effects:
 *     changes the call by reference variable outputVarList. Destroys
 *     what is there before if anything.
 *
 * Results:
 *     returns true if everything works, return false if the
 *     graph has a directed loop.
 *
 *-----------------------------------------------------------------------
 */
bool
GraphicalModel::topologicalSort(vector<RandomVariable*>& inputVarList,
				vector<RandomVariable*>& outputVarList)

{
  outputVarList.clear();
  outputVarList.resize(inputVarList.size());
  for (unsigned i=0;i<inputVarList.size();i++)
    inputVarList[i]->tag = 0;
  unsigned position=inputVarList.size();
  for (unsigned i=0;i<inputVarList.size();i++) {
    if (inputVarList[i]->tag == 0)
      if (!topologicalSortRecurse(outputVarList,
				  inputVarList[i],position))
	return false;
  }
  assert (position == 0);
  return true;
}



////////////////////////////////////////////////////////////////////
//        Test Driver
////////////////////////////////////////////////////////////////////

#ifdef MAIN


#endif
