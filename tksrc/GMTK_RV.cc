/*
 * GMTK_RV.cc
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *  $Header$
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
 *
 * The top level GMTK random variable object for the RV class hierarchy.
 *
 *
 *
 */


#include "general.h"
VCID("$Header$")

#include <iostream>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <string.h>


#include "GMTK_RV.h"


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Static variables used by classes
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


bool RV::disconnectChildrenOfObservedParents = true;


/*-
 *-----------------------------------------------------------------------
 * printParentInfo()
 *      print information about all parents of this RV.
 *
 * Preconditions:
 *      all parents and all children member variables must already
 *      be initialized with a valid graph
 *
 * Postconditions:
 *      parent info printed.
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */

void RV::printParentInfo(FILE*f, bool nl)
{
  for (unsigned i=0;i<allParents.size();i++) {
    allParents[i]->printNameFrameValue(f,false);
    if (i < (allParents.size()-1))
      fprintf(f,",");
  }
  pnl(f,nl);
}



/*-
 *-----------------------------------------------------------------------
 * printChildrenInfo()
 *      print information about all children of this RV.
 *
 * Preconditions:
 *      all parents and all children member variables must already
 *      be initialized with a valid graph
 *
 * Postconditions:
 *      parent info printed.
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */

void RV::printChildrenInfo(FILE*f, bool nl)
{
  for (unsigned i=0;i<allChildren.size();i++) {
    allChildren[i]->printNameFrameValue(f,false);
    if (i < (allChildren.size()-1))
      fprintf(f,",");
  }
  pnl(f,nl);
}


/*-
 *-----------------------------------------------------------------------
 * createNeighborsFromParentsChildren()
 *      initializes the neighbors member set with
 *      entries from all parents and all children.
 *
 *      TODO: when disconnected networks are working, do not add
 *      neighbors if BOTH: 1) they are children of this node, *and* 2)
 *      if this node is observed. The reason is that in the directed
 *      model such neighbors are superfluous since the nodes are
 *      independent.  Note that those edges might be added elsewhere
 *      (i.e., via a moralization step). In any event, adding fewer
 *      neighbors can in some cases make the triangulation
 *      possibilities more efficient.
 *
 *
 * Preconditions:
 *      all parents and all children member variables must already
 *      be initialized with a valid graph.
 *
 * Postconditions:
 *      graph now has undirected representation, where neighbors
 *      represents the undirected edges.
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      neighbors is changed.
 *
 *-----------------------------------------------------------------------
 */

void RV::createNeighborsFromParentsChildren()
{
  // make sure it is clear first
  neighbors.clear();
  // a neighbor is a parent or a child at this point.
  for (unsigned i=0;i<allParents.size();i++) {
    if (allParents[i]->observed() && disconnectChildrenOfObservedParents) {
      // do not connect parent since we are a child of an observed node.
    } else
      neighbors.insert(allParents[i]);
  }
  if (observed() && disconnectChildrenOfObservedParents) {
    // do nothing with the children.
  } else {
    // Go ahead and add the children of this as neighbors anyway,
    // ignoring the fact in the resulting UGM that 'this' renders the
    // children independent of the parents.
    for (unsigned i=0;i<allChildren.size();i++) {
      neighbors.insert(allChildren[i]);
    }
  }
  // sanity check.
  assert ( neighbors.find(this) == neighbors.end() );
}




/*-
 *-----------------------------------------------------------------------
 * connectNeighbors()
 *      Makes it such that all neighbors of self are
 *      connected (but does not touch the varibles in the
 *      'exclude' argument. This function can therefore be used
 *      as part of a node 'elimination' process.
 *
 * Preconditions:
 *      createNeighborsFromParentsChildren() must have been
 *      called at some point before.
 *
 * Postconditions:
 *      node is now such that self and all neighbors form
 *      a complete set.
 *
 * Side Effects:
 *      will change the neighbors member variables of other variables.
 *
 * Results:
 *      none.
 *
 *-----------------------------------------------------------------------
 */

void RV::connectNeighbors(set<RV*> exclude)
{
  set<RV*> nodes;

  set_difference(neighbors.begin(),neighbors.end(),
		 exclude.begin(),exclude.end(),
		 inserter(nodes,nodes.end()));
  for (set<RV*>::iterator n = nodes.begin();
       n != nodes.end();
       n++) {
    // just union together each nodes's neighbor variable
    // with self nodes
    set<RV*> tmp;
    set_union((*n)->neighbors.begin(),(*n)->neighbors.end(),
	      nodes.begin(),nodes.end(),
	      inserter(tmp,tmp.end()));
    // make sure self is not its own neighbor
    tmp.erase((*n));
    (*n)->neighbors = tmp;
  }

}



/*-
 *-----------------------------------------------------------------------
 * moralize()
 *      moralize the node, i.e., make sure that all parents of
 *      this node are neighbors of each other.
 *
 * Preconditions:
 *      createNeighborsFromParentsChildren() *MUST* have
 *      been run for *all* parents of this node.
 *
 * Postconditions:
 *      all parents are now neighbors.
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      parent's neighbors member variable is changed.
 *
 *-----------------------------------------------------------------------
 */

void RV::moralize()
{
  for (unsigned i=0;i<allParents.size();i++) {
    for (unsigned j=i+1;j<allParents.size();j++) {
      if (!disconnectChildrenOfObservedParents) {
	// then there is no chance that we will ever wish
	// an observed parent to not be connected to the rest of
	// the parents.
	allParents[i]->neighbors.insert(allParents[j]);
	allParents[j]->neighbors.insert(allParents[i]);
      } else {
	// check that both parents are hidden, and only do this if
	// they both are. The reason is that if we are not connected
	// to a parent because it it observed, there is no need to
	// moralize with respect to that parent (meaning conect that
	// parent to the other parents) since such observed parents
	// need not exist in the same clique for the child to be
	// assigned to that clique.
	if (allParents[i]->hidden() && allParents[j]->hidden()) {
	  allParents[i]->neighbors.insert(allParents[j]);
	  allParents[j]->neighbors.insert(allParents[i]);
	}
      }
    }
  }
}



/*-
 *-----------------------------------------------------------------------
 * allParentsContainedInSet()
 *      Returns true if all parents are contained within the given set.
 *
 * Preconditions:
 *      allParents member must be created.
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      bool
 *
 *-----------------------------------------------------------------------
 */
bool RV::allParentsContainedInSet(const set <RV*> givenSet)
{
  set <RV*> res;
  set_intersection(givenSet.begin(),givenSet.end(),
		   allParents.begin(),allParents.end(),
		   inserter(res,res.end()));
  return (res.size() == allParents.size());
}




/*-
 *-----------------------------------------------------------------------
 * setParents()
 *      Set the parents to the given values. Works with variables
 *      that do not have switching. Also sets children.
 *      TODO: change name of this routine.
 *
 * Preconditions:
 *      'this' Variable must not have switching.
 *
 * Postconditions:
 *      parent are re-set.
 *
 * Side Effects:
 *      parent are re-set.
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */

void RV::setParents(vector<RV *> &sparents,vector<vector<RV *> > &cpl)
{
  // since this is for a non-switching RV, we include the following
  // assertions.
  assert ( sparents.size() == 0 );
  assert ( cpl.size() == 1 );
  
  allParents = cpl[0];

  // now set this as a child of all parents, making sure to avoid
  // duplicates by creating a temporary set. Note that this might
  // change order in child array relative to parent array, but this is
  // ok, as we never rely on order of variables in child array.
  set<RV *> parentSet;
  for (unsigned i=0;i<allParents.size();i++) {
    parentSet.insert(allParents[i]);
  }
  set<RV *>::iterator si;
  for (si = parentSet.begin(); si != parentSet.end(); si++) {
    RV* rv = (*si);
    rv->allChildren.push_back(this);
  }

}




/*-
 *-----------------------------------------------------------------------
 * printRVSet{,AndValues}()
 *      Prints out the set of random variables and perhaps their values as well.
 *
 * Preconditions:
 *      f must be open, locset a set of RVs.
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      void
 *
 *-----------------------------------------------------------------------
 */
void printRVSetAndValues(FILE*f,vector<RV*>& locset,const bool nl) 
{
  bool first = true;
  for (unsigned i=0;i<locset.size();i++) {
    RV* rv = locset[i];
    if (!first)
      fprintf(f,",");
    rv->printNameFrameValue(f,false);
    first = false;
  }
  if (nl) fprintf(f,"\n");
}
void printRVSetAndValues(FILE*f,sArray<RV*>& locset,const bool nl) 
{
  bool first = true;
  for (unsigned i=0;i<locset.size();i++) {
    RV* rv = locset[i];
    if (!first)
      fprintf(f,",");
    rv->printNameFrameValue(f,false);
    first = false;
  }
  if (nl) fprintf(f,"\n");
}
void printRVSetAndValues(FILE*f,set<RV*>& locset,bool nl)
{
  bool first = true;
  set<RV*>::iterator it;
  for (it = locset.begin();
       it != locset.end();it++) {
    RV* rv = (*it);
    if (!first)
      fprintf(f,",");
    rv->printNameFrameValue(f,false);
    first = false;
  }
  if (nl) fprintf(f,"\n");
}


void printRVSet(FILE*f,vector<RV*>& locset,const bool nl) 
{
  bool first = true;
  for (unsigned i=0;i<locset.size();i++) {
    RV* rv = locset[i];
    if (!first)
      fprintf(f,",");
    rv->printNameFrame(f,false);
    first = false;
  }
  if (nl) fprintf(f,"\n");
}
void printRVSet(FILE*f,sArray<RV*>& locset,const bool nl) 
{
  bool first = true;
  for (unsigned i=0;i<locset.size();i++) {
    RV* rv = locset[i];
    if (!first)
      fprintf(f,",");
    rv->printNameFrame(f,false);
    first = false;
  }
  if (nl) fprintf(f,"\n");
}
void printRVSet(FILE*f,const set<RV*>& locset,bool nl)
{
  bool first = true;
  set<RV*>::iterator it;
  for (it = locset.begin();
       it != locset.end();it++) {
    RV* rv = (*it);
    if (!first)
      fprintf(f,",");
    rv->printNameFrame(f,false);
    first = false;
  }
  if (nl) fprintf(f,"\n");
}



void
printRVSetPtr(FILE*f,set<RV*>& locset,bool nl)
{
  bool first = true;
  set<RV*>::iterator it;
  for (it = locset.begin();
       it != locset.end();it++) {
    RV* rv = (*it);
    if (!first)
      fprintf(f,",");
    fprintf(f,"%s(%d)=0x%X",rv->name().c_str(),rv->frame(),(unsigned)rv);
    first = false;
  }
  if (nl) fprintf(f,"\n");
}




/*
 * Take the union of A and B and place it in C.
 */
void
unionRVs(const set<RV*>& A,
	 const set<RV*>& B,
	 set<RV*>& C,
	 bool do_not_clear)
{
  if (!do_not_clear)
    C.clear();
  set_union(A.begin(),A.end(),
	    B.begin(),B.end(),
	    inserter(C,C.end()));
}


/*
 * getRV: from the set of random variables that live in (rvs,pos), get
 *        the corresponding one named by 'pp' time-shifted by 'shift'.
 *        Note, shift might be positive or negative.
 *        Optionally, make a big stink if it is not there (i.e, internal error).
 */
RV * getRV(const vector <RV*>& rvs, // a set of RVs
	   map < RVInfo::rvParent, unsigned >& pos, // mappings from name(frame) to RV ptrs
	   const RVInfo::rvParent& pp, // the variable to get
	   const int shift,  
	   const bool failIfNotExist)
{

  RVInfo::rvParent desired_pp(pp.first,pp.second+shift);

  map < RVInfo::rvParent , unsigned >::iterator it;
  if ((it = pos.find(desired_pp)) == pos.end()) {
    // this could be an assertion failure as well, but we need to set 'it'
    if (failIfNotExist) 
      coredump("INTERNAL ERROR: getRV: Can't find random variable %s(%d) in unrolled collection, asked for rv %s(%d) with offset %d.\n",
	       pp.first.c_str(),pp.second+shift,
	       pp.first.c_str(),pp.second,shift);
    else
      return NULL;
  }
  return rvs[(*it).second];
}


/*
 * getRV: get from the set of random variables that live in the
 * set (rvs,pos), a time-shifted version of the random variable 'rv'.
 * Note that this routine does not require rv to be in the rv set
 * (rvs,pos) since it only depends on 'rv' via its name and frame. This
 * means we use this routine to easily grab a corresponding rv from
 * one unrolling using a rv from another unrolling using this routine.
 */
RV * getRV(const vector <RV*>& rvs, // a set of RVs
	   map < RVInfo::rvParent, unsigned >& pos, // mappings from name(frame) to RV ptrs
	   RV* rv, // the variable to be shifted
	   const int shift, // the shift amount in frames
	   const bool failIfNotExist // abort if shifted RV doesnot exist, otherwise set empty
	   )

{
  RVInfo::rvParent p(rv->name(),rv->frame()+shift);
  map < RVInfo::rvParent , unsigned >::iterator it;      
  if ((it = pos.find(p)) == pos.end()) {
    // this could be an assertion failure as well, but we need
    // to set 'it'
    if (failIfNotExist)
      coredump("INTERNAL ERROR: Can't find random variable %s(%d) when shifted by %d frames.\n",
	       rv->name().c_str(),rv->frame(),
	       shift);
    else 
      return NULL;
  }
  return rvs[(*it).second];
}


/*
 * getRVSet: get from the set of random variables that live in the set
 * (rvs,pos), a time-shifted version of the random variables in the
 * set 'rvs_to_shift'. Note that this routine does not require the rv
 * set to be in the set (rvs,pos) since it only depends on each rv via
 * its name and frame. This means we use this routine to easily grab a
 * corresponding rv from one unrolling using a rv from another
 * unrolling using this routine.
 */
set<RV*>
getRVSet(const vector <RV*>& rvs, // a set of RVs
	 map < RVInfo::rvParent, unsigned >& pos, // mappings from name(frame) to RV ptrs
	 set<RV*>& rvs_to_shift, // the variable to be shifted
	 const int shift, // the shift amount in frames
	 const bool failIfNotExist // abort if shifted RV doesnot exist, otherwise set empty
	 )
{
  set <RV*> res;
  set<RV*>::iterator i;
  for (i=rvs_to_shift.begin(); i!= rvs_to_shift.end(); i++) {
    RV *rv = (*i);
    RV *srv = getRV(rvs,pos,rv,shift,failIfNotExist);
    if (srv != NULL)
      res.insert(srv);
  }
  return res;
}



/*
 * getRVVec: vector version of getRVSet()
 */
set<RV*>
getRVVec(const vector <RV*>& rvs, // a set of RVs
	 map < RVInfo::rvParent, unsigned >& pos, // mappings from name(frame) to RV ptrs
	 vector < RVInfo::rvParent> & pps,  // the variable to be shifted
	 const int shift, // the shift amount in frames
	 const bool failIfNotExist // abort if shifted RV doesnot exist, otherwise set empty
	 )
{
  set <RV*> res;
  for (unsigned i=0;i<pps.size();i++) {
    RVInfo::rvParent& pp = pps[i];
    RV *srv = getRV(rvs,pos,pp,shift,failIfNotExist);
    if (srv != NULL)
      res.insert(srv);
  }
  return res;
}




/*
 * getRVOVec: return ordered vector, vector version of getRVSet()
 */
vector<RV*>
getRVOVec(const vector <RV*>& rvs, // a set of RVs
	 map < RVInfo::rvParent, unsigned >& pos, // mappings from name(frame) to RV ptrs
	 vector < RVInfo::rvParent> & pps,  // the variable to be shifted
	 const int shift, // the shift amount in frames
	 const bool failIfNotExist // abort if shifted RV doesnot exist, otherwise set empty
	 )
{
  vector <RV*> res;
  res.resize(pps.size());
  for (unsigned i=0;i<pps.size();i++) {
    RVInfo::rvParent& pp = pps[i];
    RV *srv = getRV(rvs,pos,pp,shift,failIfNotExist);
    if (srv != NULL)
      res[i] = srv;
  }
  return res;
}

/*
 * simple count iterator that counts the number
 * of insertions made, but doesn't do anything else.
 */
template <typename _Container>
class count_iterator 
  : public iterator<output_iterator_tag, void, void, void, void> {
  unsigned counter;
public:

  count_iterator(_Container& __x) { counter = 0; }
  count_iterator() { counter = 0; }

  // count_iterator(const count_iterator& ci) { counter = ci.counter; }
  // count_iterator& operator=(const count_iterator& ci) { counter = ci.counter; }

  count_iterator& operator=(const typename _Container::const_reference _value) 
  { counter++; return *this; }
  count_iterator& operator*() { return *this; }
  count_iterator& operator++() {  return *this; }
  count_iterator& operator++(int) { return *this; }

  void reset() { counter = 0; }
  unsigned count() { return counter; }
};

class setrv_count_iterator: public count_iterator <set <RV*> > {
public:

};


// template<typename _Container>
// inline count_iterator<_Container>
// counter(_Container& __x)
// {
//   return count_iterator<_Container>();
// }


/*
 * returns true if the first set of RVs is (not necessarily properly) contained (<=) in the
 * second set.
 * TODO: make this generic to sets of anything.
 */
bool firstRVSetContainedInSecond(set <RV*>& firstSet,
				 set <RV*>& secondSet)
{

  // TODO: figure out how to create a count_iterator without needing
  // to create dummy object.
  // set <RV*> dummy;
  // count_iterator< set <RV*> > myit(dummy);
  setrv_count_iterator myit;

  myit = set_intersection(firstSet.begin(),firstSet.end(),
			  secondSet.begin(),secondSet.end(),
			  myit);

  return (myit.count() == firstSet.size());

}




