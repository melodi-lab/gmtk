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

#include "GMTK_RV.h"
#include "GMTK_DiscRV.h"
#include "GMTK_ContRV.h"
#include "GMTK_GraphicalModel.h"

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
GraphicalModel::topologicalSort(vector<RV*>& inputVarList,
				vector<RV*>& outputVarList)

{
  map<RV *, unsigned> tag;

  outputVarList.clear();
  outputVarList.resize(inputVarList.size());
  for (unsigned i=0;i<inputVarList.size();i++)
    tag[inputVarList[i]] = 0;
  unsigned position=inputVarList.size();
  for (unsigned i=0;i<inputVarList.size();i++) {
    if (tag[inputVarList[i]] == 0)
      if (!topologicalSortRecurse(outputVarList,
				  inputVarList[i],position,tag))
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
GraphicalModel::topologicalSortRecurse(vector<RV*>& outputVarList,
				       RV* node,
				       unsigned& position,
				       map<RV*,unsigned>& tag)
{
  tag[node] = 1;
  for (unsigned i=0;i<node->allChildren.size();i++) {
    RV*rv = node->allChildren[i];
    if (tag[rv] == 0) {
      bool res = topologicalSortRecurse(outputVarList,rv,position,tag);
      if (!res)
	// directed graph has a loop
	return false;
    } else if (tag[rv] == 1)
      // directed graph has a loop
      return false;
    else
      ;
      // tag == 2, meaning we've gone down this path before and need not
      // do it again.
  }
  tag[node] = 2; // done with this node
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
GraphicalModel::topologicalSort(const set<RV*>& inputVarList,
				const set<RV*>& sortSet,
				vector<RV*>& outputVarList)

{
  map<RV *, unsigned> tag;
  outputVarList.clear();
  outputVarList.resize(inputVarList.size());
  set<RV*>::iterator it;
  for (it=inputVarList.begin();it != inputVarList.end();it++) {
    RV* rv = (*it);
    tag[rv] = 0;
  }
  unsigned position=inputVarList.size();
  for (it=inputVarList.begin();it != inputVarList.end();it++) {
    RV* rv = (*it);
    if (tag[rv] == 0)
      if (!topologicalSortRecurse(sortSet,
				  outputVarList,
				  rv,
				  position,tag))
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
GraphicalModel::topologicalSortRecurse(const set<RV*>& sortSet,
				       vector<RV*>& outputVarList,
				       RV* node,
				       unsigned& position,
				       map<RV*,unsigned>& tag)
{
  tag[node] = 1;
  for (unsigned i=0;i<node->allChildren.size();i++) {
    RV*rv = node->allChildren[i];
    // don't bother with children not in sort set.
    if (sortSet.find(rv) == sortSet.end())
      continue;
    if (tag[rv] == 0) {
      bool res = topologicalSortRecurse(sortSet,outputVarList,rv,position,tag);
      if (!res)
	// directed graph has a loop
	return false;
    } else if (tag[rv] == 1)
      // directed graph has a loop
      return false;
    else
      ;
      // tag == 2, meaning we've gone down this path before and need not
      // do it again.
  }
  tag[node] = 2; // done with this node
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
GraphicalModel::topologicalSortRandom(const set<RV*>& inputVarList,
				      const set<RV*>& sortSet,
				      vector<RV*>& outputVarList)

{
  map<RV *, unsigned> tag;
  outputVarList.clear();
  outputVarList.resize(inputVarList.size());
  set<RV*>::iterator it;
  const   set<RV*>::iterator it_end = inputVarList.end();

  sArray < unsigned > permutation(inputVarList.size());
  vector< RV*> rv_vec(inputVarList.size());
  unsigned i = 0;
  for (it=inputVarList.begin();it != it_end;it++) {
    RV* rv = (*it);
    tag[rv] = 0;
    permutation.ptr[i] = i;
    rv_vec[i] = rv;
    i++;
  }
  rnd.rpermute(permutation.ptr,inputVarList.size());
  unsigned position=inputVarList.size();
  for (i=0;i<inputVarList.size();i++) {
    RV* rv = rv_vec[permutation.ptr[i]];
    if (tag[rv] == 0)
      if (!topologicalSortRecurseRandom(sortSet,
					outputVarList,
					rv,
					position,tag))
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
GraphicalModel::topologicalSortRecurseRandom(const set<RV*>& sortSet,
					     vector<RV*>& outputVarList,
					     RV* node,
					     unsigned& position,
					     map<RV*,unsigned>& tag)
{

  tag[node] = 1;

  const unsigned nChildren = node->allChildren.size();
  sArray <unsigned > permutation(nChildren);
  for (unsigned i=0;i<nChildren;i++) {
    permutation.ptr[i] = i;
  }
  rnd.rpermute(permutation.ptr,nChildren);
  for (unsigned i=0;i<nChildren;i++) {
    RV*rv = node->allChildren[permutation.ptr[i]];
    // don't bother with children not in sort set.
    if (sortSet.find(rv) == sortSet.end())
      continue;
    if (tag[rv] == 0) {
      bool res = topologicalSortRecurseRandom(sortSet,outputVarList,rv,position,tag);
      if (!res)
	// directed graph has a loop
	return false;
    } else if (tag[rv] == 1)
      // directed graph has a loop
      return false;
    else
      ;
      // tag == 2, meaning we've gone down this path before and need not
      // do it again.
  }
  tag[node] = 2; // done with this node
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
GraphicalModel::topologicalSortWPriority(const set<RV*>& inputVarList,
					 const set<RV*>& sortSet,
					 vector<RV*>& outputVarList,
					 const string priorityStr)
{
  map<RV *, unsigned> tag;

  outputVarList.clear();
  outputVarList.resize(inputVarList.size());
  set<RV*>::iterator it;
  for (it=inputVarList.begin();it != inputVarList.end();it++) {
    RV* rv = (*it);
    tag[rv] = 0;
  }
  unsigned position=0;

  // Different cases of sort to be prioritized by position in string:
  //    C: continuous, observed
  //    D: discrete, observed, deterministic
  //    O: observed
  //    B: binary
  //    S: switching
  //    I: ever increasing cardinality of variables.
  //    A: alphabetical (for debugging purposes)
  //    F: by frame number (for debugging purposes)
  //    N: first by alpha name and next by frame.
  // 
  // Some good ones to try for speed:  
  //   COB, 
  //   CDOI or CODI (good when no continuous vars)
  //   DOI
  //

  for (unsigned charNo=0;charNo< priorityStr.size(); charNo++) {
    const char curCase = toupper(priorityStr[charNo]);

    if (curCase == 'C') {

      // Do a pass doing just continuous observed variables We do these
      // here to avoid having to re-compute them many times (in case
      // caching is turned off).
      for (it=inputVarList.begin();it != inputVarList.end();it++) {
	RV* rv = (*it);
	if (rv->discrete() || rv->hidden())
	  continue;
	if (tag[rv] == 0)
	  if (!topologicalSortRecurseWPriorityRecurse(sortSet,
						      outputVarList,
						      rv,
						      position,tag))
	    return false;
      }
    } else if (curCase == 'D') {

      // Do a pass for discrete *deterministic* (meaning really
      // determinisitc, not sparse) observed variable.  Getting these in
      // early are likely to kill off a large chunk of unnecessary
      // computation.
      for (it=inputVarList.begin();it != inputVarList.end();it++) {
	RV* rv = (*it);
	if (!(rv->discrete() && RV2DRV(rv)->deterministic() && !rv->hidden()))
	  continue;
	if (tag[rv] == 0)
	  if (!topologicalSortRecurseWPriorityRecurse(sortSet,
						      outputVarList,
						      rv,
						      position,tag))
	    return false;
      }

    } else if (curCase == 'S') {

      // Do a pass for switching parents, namely variables whose children
      // use these as switching parents. These must be discrete
      // at the moment.

      for (it=inputVarList.begin();it != inputVarList.end();it++) {
	RV* rv = (*it);
	if (!rv->discrete() || rv->allChildren.size() == 0)
	  continue;
	// ok, so discrete, but is it used as a switching parent by at
	// least one of its children?

	bool sw_parent = false;
	for (unsigned c=0;c<rv->allChildren.size();c++) {
	  RV* chld = rv->allChildren[c];
	  if (chld->switching()) {
	    for (unsigned p=0;p<chld->switchingParentsVec().size();p++) {
	      if (chld->switchingParentsVec()[p] == rv) {
		sw_parent = true;
		goto foundSwParent;
	      }
	    }
	  }
	}
      foundSwParent:
	
	if (!sw_parent)
	  continue;

	if (tag[rv] == 0)
	  if (!topologicalSortRecurseWPriorityRecurse(sortSet,
						      outputVarList,
						      rv,
						      position,tag))
	    return false;
      }

    } else if (curCase == 'O') {

      // Next do a pass for any other observed variables.
      for (it=inputVarList.begin();it != inputVarList.end();it++) {
	RV* rv = (*it);
	if (rv->hidden())
	  continue;
	if (tag[rv] == 0)
	  if (!topologicalSortRecurseWPriorityRecurse(sortSet,
						      outputVarList,
						      rv,
						      position,tag))
	    return false;
      }
    } else if (curCase == 'I') {

      // Sort the reminaing discrete variables in increasing order of
      // cardinality and do a pass in that increasing order.
      multimap< unsigned ,RV*> cardSortedNodes;
      for (it=inputVarList.begin();it != inputVarList.end();it++) {
	RV* rv = (*it);
	if (!rv->discrete())
	  continue;
	DiscRV* drv = RV2DRV(rv);
	pair< unsigned , RV*> pr (drv->cardinality, rv );
	cardSortedNodes.insert(pr);    
      }
      for (multimap< unsigned, RV*>::iterator m = cardSortedNodes.begin();
	   m != cardSortedNodes.end(); m++) {
	// unsigned lcard = (*m).first;
	RV* rv = (*m).second;
	// DiscRV* drv = (DiscRV*)rv;
	// printf("Doing node %s(%d) with card %d\n",(*m).second->name().c_str(),(*m).second->frame(),drv->cardinality);
	if (tag[rv] == 0)
	  if (!topologicalSortRecurseWPriorityRecurse(sortSet,
						      outputVarList,
						      rv,
						      position,tag))
	    return false;
      }
    } else if (curCase == 'F') {

      // Sort the reminaing variables in increasing order of
      // frame and do a pass in that increasing order.
      multimap< unsigned ,RV*> cardSortedNodes;
      for (it=inputVarList.begin();it != inputVarList.end();it++) {
	RV* rv = (*it);
	if (!rv->discrete())
	  continue;
	DiscRV* drv = RV2DRV(rv);
	pair< unsigned , RV*> pr (drv->frame(), rv );
	cardSortedNodes.insert(pr);    
      }
      for (multimap< unsigned, RV*>::iterator m = cardSortedNodes.begin();
	   m != cardSortedNodes.end(); m++) {
	// unsigned lcard = (*m).first;
	RV* rv = (*m).second;
	// DiscRV* drv = (DiscRV*)rv;
	// printf("Doing node %s(%d) with frame %d\n",(*m).second->name().c_str(),(*m).second->frame(),drv->cardinality);
	if (tag[rv] == 0)
	  if (!topologicalSortRecurseWPriorityRecurse(sortSet,
						      outputVarList,
						      rv,
						      position,tag))
	    return false;
      }
    } else if (curCase == 'A') {
      // Sort the reminaing variables in increasing alphabetical order of
      // name and do a pass in that increasing order.
      multimap< string ,RV*> cardSortedNodes;
      for (it=inputVarList.begin();it != inputVarList.end();it++) {
	RV* rv = (*it);
	pair< string , RV*> pr (rv->name(), rv );
	cardSortedNodes.insert(pr);    
      }
      for (multimap< string, RV*>::iterator m = cardSortedNodes.begin();
	   m != cardSortedNodes.end(); m++) {
	// unsigned lcard = (*m).first;
	RV* rv = (*m).second;
	// DiscRV* drv = (DiscRV*)rv;
	// printf("Doing node %s(%d) with name\n",(*m).second->name().c_str(),(*m).second->frame());
	if (tag[rv] == 0)
	  if (!topologicalSortRecurseWPriorityRecurse(sortSet,
						      outputVarList,
						      rv,
						      position,tag))
	    return false;
      }
    } else if (curCase == 'N') {
      // Sort the reminaing variables in increasing lexicographic  order of
      // string name and frame number, and do a pass in that increasing order.
      multimap< pair<string,unsigned> ,RV*> cardSortedNodes;
      for (it=inputVarList.begin();it != inputVarList.end();it++) {
	RV* rv = (*it);
	pair< pair<string,unsigned> , RV*> pr ( pair<string,unsigned> (rv->name(),rv->frame()), rv );
	cardSortedNodes.insert(pr);    
      }
      for (multimap< pair<string,unsigned>, RV*>::iterator m = cardSortedNodes.begin();
	   m != cardSortedNodes.end(); m++) {
	// unsigned lcard = (*m).first;
	RV* rv = (*m).second;
	// DiscRV* drv = (DiscRV*)rv;
	// printf("Doing node %s(%d) with name,frame\n",(*m).second->name().c_str(),(*m).second->frame());
	if (tag[rv] == 0)
	  if (!topologicalSortRecurseWPriorityRecurse(sortSet,
						      outputVarList,
						      rv,
						      position,tag))
	    return false;
      }
    } else if (curCase == 'B') {
      // next do a pass for any binary variables.
      for (it=inputVarList.begin();it != inputVarList.end();it++) {
	RV* rv = (*it);
	if (!rv->discrete())
	  continue;
	DiscRV* drv = RV2DRV(rv);
	if (drv->cardinality > 2)
	  continue;
	if (tag[rv] == 0)
	  if (!topologicalSortRecurseWPriorityRecurse(sortSet,
						      outputVarList,
						      rv,
						      position,tag))
	    return false;
      }
    }
  }

  // Last do a pass to hit any remainder that the above might not have
  // done.
  for (it=inputVarList.begin();it != inputVarList.end();it++) {
    RV* rv = (*it);
    if (tag[rv] == 0)
      if (!topologicalSortRecurseWPriorityRecurse(sortSet,
						  outputVarList,
						  rv,
						  position,tag))
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
GraphicalModel::topologicalSortRecurseWPriorityRecurse(const set<RV*>& sortSet,
						       vector<RV*>& outputVarList,
						       RV* node,
						       unsigned& position,
						       map<RV*,unsigned>& tag)
{
  tag[node] = 1;
  for (unsigned i=0;i<node->allParents.size();i++) {
    RV*rv = node->allParents[i];
    // don't bother with parents not in sort set.
    if (sortSet.find(rv) == sortSet.end())
      continue;
    if (tag[rv] == 0) {
      bool res = topologicalSortRecurseWPriorityRecurse(sortSet,outputVarList,rv,position,tag);
      if (!res)
	// directed graph has a loop
	return false;
    } else if (tag[rv] == 1)
      // directed graph has a loop
      return false;
    else
      ;
      // tag == 2, meaning we've gone up this path before and need not
      // do it again.
  }
  tag[node] = 2; // done with this node
  outputVarList[position++] = node;
  return true;
}





////////////////////////////////////////////////////////////////////
//        Test Driver
////////////////////////////////////////////////////////////////////

#ifdef MAIN


#endif
