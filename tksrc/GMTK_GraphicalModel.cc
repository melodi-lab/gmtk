/*-
 * GMTK_GraphicalModel.cc
 *     various support routines for graphical models
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


////////////////////////////////////////////////////////////////////
//        Basic topological sort & support
////////////////////////////////////////////////////////////////////




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
 *     what is there before if anything. Changes the tag member of each
 *     random variable.
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



////////////////////////////////////////////////////////////////////
//        constrained topological sort & support
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * GraphicalModel::topologicalSort
 *      Another topological sort, just like the one above, but this
 *      one constrains the sort to be determined with respect only to 
 *      those variables that are in the given sortSet argument. I.e.,
 *      its as if each variables children are considered to be its
 *      real children intersected with sortSet (so some children are
 *      never considered in the sort, nor are they traversed). 
 *
 *      This is useful when we want to sort a set of variables that
 *      themselves have a partial order, but some of the variables
 *      point to children outside of the current set, and we want
 *      to ignore those children, and only sort with respect to the
 *      variables in sortSet.
 * 
 *
 * Preconditions:
 *     inputVarLIst must contain a set of random variables with
 *     the variables appropriately set up.
 *
 * Postconditions:
 *     ouputVarList contains the list of variables in order.
 *
 * Side Effects:
 *     changes the call by reference variable outputVarList. Destroys
 *     what is there before if anything. Changes the tag member of each
 *     random variable.
 *
 * Results:
 *     returns true if everything works, return false if the
 *     graph has a directed loop.
 *
 *-----------------------------------------------------------------------
 */
bool
GraphicalModel::topologicalSort(const set<RandomVariable*>& inputVarList,
				const set<RandomVariable*>& sortSet,
				vector<RandomVariable*>& outputVarList)

{
  outputVarList.clear();
  outputVarList.resize(inputVarList.size());
  set<RandomVariable*>::iterator it;
  for (it=inputVarList.begin();it != inputVarList.end();it++) {
    RandomVariable* rv = (*it);
    rv->tag = 0;
  }
  unsigned position=inputVarList.size();
  for (it=inputVarList.begin();it != inputVarList.end();it++) {
    RandomVariable* rv = (*it);
    if (rv->tag == 0)
      if (!topologicalSortRecurse(sortSet,
				  outputVarList,
				  rv,
				  position))
	return false;
  }
  assert (position == 0);
  return true;
}


/*-
 *-----------------------------------------------------------------------
 * GraphicalModel::topologicalSortRecurse
 *      Support routine for the topological Sort routine that
 *      uses the sortSet argument.
 *
 * Preconditions:
 *      must only be called from topological sort with the sortSet argument.
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
GraphicalModel::topologicalSortRecurse(const set<RandomVariable*>& sortSet,
				       vector<RandomVariable*>& outputVarList,
				       RandomVariable* node,
				       unsigned& position)
{
  node->tag = 1;
  for (unsigned i=0;i<node->allPossibleChildren.size();i++) {
    RandomVariable*rv = node->allPossibleChildren[i];
    // don't bother with children not in sort set.
    if (sortSet.find(rv) == sortSet.end())
      continue;
    if (rv->tag == 0) {
      bool res = topologicalSortRecurse(sortSet,outputVarList,rv,position);
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




////////////////////////////////////////////////////////////////////
//        constrained & Random topological sort & support
////////////////////////////////////////////////////////////////////



/*-
 *-----------------------------------------------------------------------
 * GraphicalModel::topologicalSort
 *      Another topological sort, just like the one above, but this
 *      one constrains the sort to be determined with respect only to 
 *      those variables that are in the given sortSet argument. I.e.,
 *      its as if each variables children are considered to be its
 *      real children intersected with sortSet (so some children are
 *      never considered in the sort, nor are they traversed). 
 *
 *      This version also produces a random sort, so if it is
 *      called multiple times (or with a different seed), a different
 *      resulting topological sort will result.
 *
 *      This is useful when we want to sort a set of variables that
 *      themselves have a partial order, but some of the variables
 *      point to children outside of the current set, and we want
 *      to ignore those children, and only sort with respect to the
 *      variables in sortSet.
 * 
 *
 * Preconditions:
 *     inputVarLIst must contain a set of random variables with
 *     the variables appropriately set up.
 *
 * Postconditions:
 *     ouputVarList contains the list of variables in order.
 *
 * Side Effects:
 *     changes the call by reference variable outputVarList. Destroys
 *     what is there before if anything. Changes the tag member of each
 *     random variable.
 *
 * Results:
 *     returns true if everything works, return false if the
 *     graph has a directed loop.
 *
 *-----------------------------------------------------------------------
 */
bool
GraphicalModel::topologicalSortRandom(const set<RandomVariable*>& inputVarList,
				      const set<RandomVariable*>& sortSet,
				      vector<RandomVariable*>& outputVarList)

{
  outputVarList.clear();
  outputVarList.resize(inputVarList.size());
  set<RandomVariable*>::iterator it;
  const   set<RandomVariable*>::iterator it_end = inputVarList.end();

  sArray < unsigned > permutation(inputVarList.size());
  vector< RandomVariable*> rv_vec(inputVarList.size());
  unsigned i = 0;
  for (it=inputVarList.begin();it != it_end;it++) {
    RandomVariable* rv = (*it);
    rv->tag = 0;
    permutation.ptr[i] = i;
    rv_vec[i] = rv;
    i++;
  }
  rnd.rpermute(permutation.ptr,inputVarList.size());
  unsigned position=inputVarList.size();
  for (i=0;i<inputVarList.size();i++) {
    RandomVariable* rv = rv_vec[permutation.ptr[i]];
    if (rv->tag == 0)
      if (!topologicalSortRecurseRandom(sortSet,
					outputVarList,
					rv,
					position))
	return false;
  }
  assert (position == 0);
  return true;
}


/*-
 *-----------------------------------------------------------------------
 * GraphicalModel::topologicalSortRecurseRandom
 *      Support routine for the topological Sort routine that
 *      uses the sortSet argument and is random.
 *
 * Preconditions:
 *      must only be called from topological sort with the sortSet argument.
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
GraphicalModel::topologicalSortRecurseRandom(const set<RandomVariable*>& sortSet,
					     vector<RandomVariable*>& outputVarList,
					     RandomVariable* node,
					     unsigned& position)
{
  node->tag = 1;

  const unsigned nChildren = node->allPossibleChildren.size();
  sArray <unsigned > permutation(nChildren);
  for (unsigned i=0;i<nChildren;i++) {
    permutation.ptr[i] = i;
  }
  rnd.rpermute(permutation.ptr,nChildren);
  for (unsigned i=0;i<nChildren;i++) {
    RandomVariable*rv = node->allPossibleChildren[permutation.ptr[i]];
    // don't bother with children not in sort set.
    if (sortSet.find(rv) == sortSet.end())
      continue;
    if (rv->tag == 0) {
      bool res = topologicalSortRecurseRandom(sortSet,outputVarList,rv,position);
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
 * GraphicalModel::topologicalSortContObsFirst
 *      Another topological sort, just like the ones above, but this
 *      one constrains the sort to be determined with respect only to 
 *      those variables that are in the given sortSet argument. I.e.,
 *      its as if each variables children are considered to be its
 *      real children intersected with sortSet (so some children are
 *      never considered in the sort, nor are they traversed). 
 *
 *      This version also produces a sort but this one produces sorts
 *      where continuous observations are as early in the sort as
 *      possible. This is thus useful to sort nodes in a clique
 *      to avoid calling the cont. Gaussian more often than necessary.
 *      
 *
 * Preconditions:
 *     inputVarLIst must contain a set of random variables with
 *     the variables appropriately set up.
 *
 * Postconditions:
 *     ouputVarList contains the list of variables in order.
 *
 * Side Effects:
 *     changes the call by reference variable outputVarList. Destroys
 *     what is there before if anything. Changes the tag member of each
 *     random variable.
 *
 * Results:
 *     returns true if everything works, return false if the
 *     graph has a directed loop.
 *
 *-----------------------------------------------------------------------
 */
bool
GraphicalModel::topologicalSortContFirst(const set<RandomVariable*>& inputVarList,
					    const set<RandomVariable*>& sortSet,
					    vector<RandomVariable*>& outputVarList)
{
  outputVarList.clear();
  outputVarList.resize(inputVarList.size());
  set<RandomVariable*>::iterator it;
  for (it=inputVarList.begin();it != inputVarList.end();it++) {
    RandomVariable* rv = (*it);
    rv->tag = 0;
  }
  unsigned position=0;
  // first do a pass doing just continuous variables
  for (it=inputVarList.begin();it != inputVarList.end();it++) {
    RandomVariable* rv = (*it);
    if (rv->discrete == true)
      continue;
    if (rv->tag == 0)
      if (!topologicalSortRecurseContFirst(sortSet,
					      outputVarList,
					      rv,
					      position))
	return false;
  }

  // sort using other criterion as well.

  // next do a pass for any observed variables.
  for (it=inputVarList.begin();it != inputVarList.end();it++) {
    RandomVariable* rv = (*it);
    if (rv->hidden == true)
      continue;
    if (rv->tag == 0)
      if (!topologicalSortRecurseContFirst(sortSet,
					      outputVarList,
					      rv,
					      position))
	return false;
  }
  // next do a pass for any binary variables.
  for (it=inputVarList.begin();it != inputVarList.end();it++) {
    RandomVariable* rv = (*it);
    if (!rv->discrete == true)
      continue;
    DiscreteRandomVariable* drv = (DiscreteRandomVariable*)rv;
    if (drv->cardinality > 2)
      continue;
    if (rv->tag == 0)
      if (!topologicalSortRecurseContFirst(sortSet,
					      outputVarList,
					      rv,
					      position))
	return false;
  }

  // could continue here doing ternary, etc. variables depending
  // on what exists in set of nodes.

  // last do a pass to hit any remainder.

  for (it=inputVarList.begin();it != inputVarList.end();it++) {
    RandomVariable* rv = (*it);
    if (rv->tag == 0)
      if (!topologicalSortRecurseContFirst(sortSet,
					      outputVarList,
					      rv,
					      position))
	return false;
  }
  assert (position == inputVarList.size());
  return true;
}


/*-
 *-----------------------------------------------------------------------
 * GraphicalModel::topologicalSortRecurseRandom
 *      Support routine for the topological Sort routine that
 *      uses the sortSet argument and is random.
 *
 * Preconditions:
 *      must only be called from topological sort with the sortSet argument.
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
GraphicalModel::topologicalSortRecurseContFirst(const set<RandomVariable*>& sortSet,
						   vector<RandomVariable*>& outputVarList,
						   RandomVariable* node,
						   unsigned& position)
{
  node->tag = 1;
  for (unsigned i=0;i<node->allPossibleParents.size();i++) {
    RandomVariable*rv = node->allPossibleParents[i];
    // don't bother with parents not in sort set.
    if (sortSet.find(rv) == sortSet.end())
      continue;
    if (rv->tag == 0) {
      bool res = topologicalSortRecurseContFirst(sortSet,outputVarList,rv,position);
      if (!res)
	// directed graph has a loop
	return false;
    } else if (rv->tag == 1)
      // directed graph has a loop
      return false;
    else
      ;
      // tag == 2, meaning we've gone up this path before and need not
      // do it again.
  }
  node->tag = 2; // done with this node
  outputVarList[position++] = node;
  return true;
}





////////////////////////////////////////////////////////////////////
//        Test Driver
////////////////////////////////////////////////////////////////////

#ifdef MAIN


#endif
