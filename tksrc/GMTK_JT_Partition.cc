/*-
 * GMTK_JT_Partition.cc
 *     Junction Tree
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2009, < fill in later >
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
#include <ctype.h>

#include <iterator>
#include <map>
#include <set>
#include <algorithm>
#include <new>

#include "general.h"
#include "error.h"
#include "debug.h"
#include "rand.h"

#include "GMTK_FileParser.h"
#include "GMTK_RV.h"
#include "GMTK_DiscRV.h"
#include "GMTK_GMTemplate.h"
#include "GMTK_JunctionTree.h"
#include "GMTK_GMParms.h"
#include "GMTK_JT_Partition.h"

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Static variables used by classes
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

VCID("$Header$")


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Constructors/Destructors 
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * JT_Partition constructor.
 *
 *   This one creates a JT_Parititon from a partition and is used as
 *   part of the JT unroll function. It creates a set of cliques
 *   isomorphic to the Parititon, but all RVs in the new cliques have
 *   had their frame number adjusted by a frame delta. The constructor
 *   also includes arguments for the left interface variables (liVars)
 *   and right interface variables (riVars) and appropriate frame
 *   deltas for those. The reason for having two frame deltas is that
 *   the from interfaces might actually be the same, the only
 *   difference in this partition is if it is on the left or right of
 *   the new partition, something determined only by the frame delta.
 *
 *
 * Preconditions:
 *   None.
 *
 * Postconditions:
 *   JT_Partition is created and set up.
 *
 * Side Effects:
 *   Adjusts member variables.
 *
 * Results:
 *   none.
 *
 *-----------------------------------------------------------------------
 */
JT_Partition::JT_Partition(
		       Partition& from_part,
		       const unsigned int frameDelta,
		       // the left and right interface variables for
		       // this JT partition Empty if doesn't exist
		       // (say for an P or E partition). These have
		       // their own frame deltas since they might be
		       // different.
		       // Left interface:
		       const set <RV*>& from_liVars,
		       const unsigned int liFrameDelta,
		       // Right interface:
		       const set <RV*>& from_riVars,
		       const unsigned int riFrameDelta,
		       // Information todo the mapping.
		       vector <RV*>& newRvs,
		       map < RVInfo::rvParent, unsigned >& ppf)

{

  triMethod = from_part.triMethod;
  numVEseps = 0;

  set<RV*>::iterator it;

  // clone over nodes RVs.  
  // TODO: make this next code a routine
  //  nodesClone() since it is used in several places.
  for (it = from_part.nodes.begin();
       it != from_part.nodes.end();
       it++) {
    RV* rv = (*it);
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    

    if ( ppf.find(rvp) == ppf.end() ) {
      coredump("INTERNAL ERROR: can't find rv %s(%d+%d)=%s(%d) in unrolled RV set, from from_part2\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
    }

    RV* nrv = newRvs[ppf[rvp]];
    nodes.insert(nrv);
  }
  liNodes.clear();
  for (it = from_liVars.begin();
       it != from_liVars.end();
       it++) {
    RV* rv = (*it);
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+liFrameDelta;    

    if ( ppf.find(rvp) == ppf.end() ) {
      coredump("INTERNAL ERROR: can't find rv %s(%d+%d)=%s(%d) in unrolled RV set, from from_liVars\n",
	    rv->name().c_str(),rv->frame(),liFrameDelta,
	    rvp.first.c_str(),rvp.second);
    }

    RV* nrv = newRvs[ppf[rvp]];
    liNodes.insert(nrv);
  }
  riNodes.clear();
  for (it = from_riVars.begin();
       it != from_riVars.end();
       it++) {
    RV* rv = (*it);
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+riFrameDelta;    

    if ( ppf.find(rvp) == ppf.end() ) {
      coredump("INTERNAL ERROR: can't find rv %s(%d+%d)=%s(%d) in unrolled RV set, from from_riVars\n",
	    rv->name().c_str(),rv->frame(),riFrameDelta,
	    rvp.first.c_str(),rvp.second);
    }

    RV* nrv = newRvs[ppf[rvp]];
    riNodes.insert(nrv);
  }

  // make the cliques.
  cliques.reserve(from_part.cliques.size());
  // 
  // NOTE: It is ***CRUCIAL*** for the cliques in the cloned partition
  // to be inserted in the *SAME ORDER* as in the partition being
  // cloned. If this is not done, inference will crash.
  for (unsigned i=0;i<from_part.cliques.size();i++) {
    cliques.push_back(MaxClique(from_part.cliques[i],
				newRvs,ppf,frameDelta));
  }
}


/*
 * clear all memory in the origin clique and separator objects that
 * are associated with this object.
 */
void
JT_Partition::clearCliqueAndIncommingSeparatorMemoryForClique(unsigned cliqueNo)
{
  MaxClique& clique = cliques[cliqueNo];
  // first do the separators
  for (unsigned sepNumber=0;sepNumber<clique.ceReceiveSeparators.size();sepNumber++) {
    SeparatorClique& sep = 
      separators[clique.ceReceiveSeparators[sepNumber]];
    sep.clearInferenceMemory();
  }
  // and next do self.
  clique.clearInferenceMemory();
}



/*-
 *-----------------------------------------------------------------------
 * JT_Partition constructor.
 *
 *   A bear-bones constructor that does no frame adjustment and is
 *   here for the sake of quickly being able to compute JT weight from
 *   within a triangulation routine (before we really have the full JT
 *   set up). This routine should only be called from
 *   junctionTreeWeight().
 *
 *
 * Preconditions:
 *   None.
 *
 * Postconditions:
 *   JT_Partition is created and set up, but for a specific purpose (see above comment).
 *
 * Side Effects:
 *   Adjusts member variables.
 *
 * Results:
 *   none.
 *
 *-----------------------------------------------------------------------
 */
JT_Partition::JT_Partition(
 Partition& from_part,
 // Left interface:
 const set <RV*>& from_liVars,
 // Right interface:
 const set <RV*>& from_riVars
 )
{
  nodes = from_part.nodes;
  triMethod = from_part.triMethod;
  numVEseps = 0;
  liNodes = from_liVars;
  riNodes = from_riVars;
  // make the cliques.
  cliques = from_part.cliques;
}



/*-
 *-----------------------------------------------------------------------
 * JT_Partition::findInterfaceCliques()
 *
 *   This routine looks in the current JT_Partition and finds a clique
 *   that overlaps the iNodes. It is assumed, for example, that iNodes
 *   are the partition interface nodes (either left or right) and the goal
 *   is to find a clique in the current partition that could serve as an
 *   interface clique. If there are more than one candiate interface clique
 *   available, then the one with the smallest normal weight is chosen.
 *
 *   uses command line option: interfaceCliquePriorityStr 
 *
 *
 * Preconditions:
 *   JT_partition tree must have been created and all set up.
 *
 * Postconditions:
 *   none.
 *
 *
 * Side Effects:
 *   None
 *
 * Results:
 *   - candiate interface clique number in this JT_Partition via iCliques.
 *   - if the clique returned is the exact same set of nodes as the interface
 *     nodes, then iCliqueSameAsInterface is set to true on return.
 *
 *-----------------------------------------------------------------------
 */
// first, we define two structures that are used for sorting/comparing 
// weighted cliques.
struct PriorityClique {
  unsigned clique;
  vector <double> weights;
  PriorityClique(unsigned c,vector <double> w) : clique(c), weights(w) {}
};
struct PriorityCliqueCompare {  
  // sort descending
  bool operator() (const PriorityClique& a, 
		   const PriorityClique& b) {
    return (a.weights) > (b.weights);
  }
};
void
JT_Partition::findInterfaceCliques(const set <RV*>& iNodes,
				   unsigned& iClique,
				   bool& iCliqueSameAsInterface,
				   const string priorityStr)
{
  if (iNodes.size() > 0) {
    vector < PriorityClique > pcArray;

    for (unsigned cliqueNo=0;cliqueNo<cliques.size();cliqueNo++) {
      // check that clique fully covers iNodes 
      set<RV*> res;
      set_intersection(cliques[cliqueNo].nodes.begin(),
		       cliques[cliqueNo].nodes.end(),
		       iNodes.begin(),
		       iNodes.end(),
		       inserter(res,res.end()));
      
      if (res.size() == iNodes.size()) {
	// we've found a candidate.

	// Priority for determining which clique becomes the interface
	// clique.  All attributes get sorted in decreasing order, so
	// that larger values have priority.  Larger numbers are
	// prefered.

	// Options to include by priority of priorityStr:
	//   * W: negative weight of the clique
	//   * D: number of deterministic nodes in clique
	//   * H: number of hidden nodes in clique
	//   * O: number of observed nodes in clique
	//   * I: Interface Clique same as Interface
	// 
	//  Any can be preceeded by a '-' sign to flip effect. E.g., D-E-SU
	//
	// Default case: W

	vector <double> weights;
	float mult = 1.0;
	bool loc_iCliqueSameAsInterface = false;
	for (unsigned charNo=0;charNo< priorityStr.size(); charNo++) {
	  const char curCase = toupper(priorityStr[charNo]);
	  if (curCase == '-') {
	    mult = -1.0;
	    continue;
	  }

	  if (res.size() == cliques[cliqueNo].nodes.size()) {
	    loc_iCliqueSameAsInterface = true;
	  } else {
	    loc_iCliqueSameAsInterface = false;
	  }

	  if (curCase == 'W') {	  
	    weights.push_back(-mult*(double)MaxClique::computeWeight(cliques[cliqueNo].nodes));
	  } else if (curCase == 'D') {
	    set<RV*>::iterator it;
	    set<RV*>::iterator it_end = cliques[cliqueNo].nodes.end();
	    unsigned numDeterministicNodes = 0;
	    for (it = cliques[cliqueNo].nodes.begin(); it != it_end; it++) {
	      RV* rv = (*it);
	      if (rv->discrete() && RV2DRV(rv)->deterministic())
		numDeterministicNodes++;
	    }
	    weights.push_back(mult*(double)numDeterministicNodes);
	  } else if (curCase == 'H') {	  
	    set<RV*>::iterator it;
	    set<RV*>::iterator it_end = cliques[cliqueNo].nodes.end();
	    unsigned numHiddenNodes = 0;
	    for (it = cliques[cliqueNo].nodes.begin(); it != it_end; it++) {
	      RV* rv = (*it);
	      if (rv->hidden())
		numHiddenNodes++;
	    }
	    weights.push_back(mult*(double)numHiddenNodes);
	  } else if (curCase == 'O') {
	    set<RV*>::iterator it;
	    set<RV*>::iterator it_end = cliques[cliqueNo].nodes.end();
	    unsigned numObservedNodes = 0;
	    for (it = cliques[cliqueNo].nodes.begin(); it != it_end; it++) {
	      RV* rv = (*it);
	      if (rv->observed())
		numObservedNodes++;
	    }
	    weights.push_back(mult*(double)numObservedNodes);
	  } else if (curCase == 'I') {
	    weights.push_back(mult*(double)loc_iCliqueSameAsInterface);
	  } else {
	    error("ERROR: Unrecognized interface clique determiner letter '%c' in string '%s'\n",
		  curCase,priorityStr.c_str());
	  }

	  mult = 1.0;
	}
	// always push back at end.
	weights.push_back((double)loc_iCliqueSameAsInterface);

	pcArray.push_back(PriorityClique(cliqueNo,weights));
      }
    }

    assert ( pcArray.size() > 0 );
    sort(pcArray.begin(),pcArray.end(),PriorityCliqueCompare());
    iClique = pcArray[0].clique;
    iCliqueSameAsInterface  = (bool)pcArray[0].weights[pcArray[0].weights.size()-1];

  } else {
    // This could happen either when (case 1) we are trying to find
    // the r-interface for E, the left interface for P, or (case 2)
    // for disconnected graphs or when some of the partitions (P, or
    // E) are empty.  In case 1, it doesn't matter what we return (it
    // could be an invalid value such as ~0x0). In the disconnected
    // graph case, we need to return any valid clique index. Since
    // there is always a clique, we return a 0.
    iClique = 0;
    iCliqueSameAsInterface = false;
  }
}




/*-
 *-----------------------------------------------------------------------
 * JT_Partition::find{L,R}InterfaceCliques()
 *
 *   Two quick routiens to find left/right interface cliques. See
 *   the called routine for documentation.
 *
 * Preconditions:
 *   JT_partition tree must have been created and all set up.
 *
 * Postconditions:
 *   none.
 *
 *
 * Side Effects:
 *   None
 *
 * Results:
 *
 *-----------------------------------------------------------------------
 */
void
JT_Partition::findLInterfaceClique(unsigned& liClique,bool& liCliqueSameAsInterface,const string priorityStr)
{
  findInterfaceCliques(liNodes,liClique,liCliqueSameAsInterface,priorityStr);
}
void
JT_Partition::findRInterfaceClique(unsigned& riClique,bool& riCliqueSameAsInterface,const string priorityStr)
{
  findInterfaceCliques(riNodes,riClique,riCliqueSameAsInterface,priorityStr);
}




/*-
 *-----------------------------------------------------------------------
 * JT_Partition::cliqueWith{Max,Min}Weight()
 *
 *   Return the clique in the current junction tree with maximum {minimum} 
 *   normal weight. 
 *
 * Preconditions:
 *   JT_partition tree must have been created and all set up.
 *
 * Postconditions:
 *   none.
 *
 * Side Effects:
 *   None
 *
 * Results:
 *   returns the index of the max/min weight clique in this JT Part.
 *
 *-----------------------------------------------------------------------
 */
unsigned
JT_Partition::cliqueWithMaxWeight()
{
  double weight = -10e20;
  unsigned res= ~0x0;
  // TODO: for empty P or E partitions, this next assertion will need to be removed.
  // assert ( cliques.size() > 0 );
  for (unsigned i=0;i<cliques.size();i++) {
    if (MaxClique::computeWeight(cliques[i].nodes) > weight) {
      res = i;
    }
  }
  return res;
}
unsigned
JT_Partition::cliqueWithMinWeight()
{
  double weight = DBL_MAX;
  unsigned res = ~0x0;
  // TODO: for empty P or E partitions, this next assertion will need to be removed.
  // assert ( cliques.size() > 0 );
  for (unsigned i=0;i<cliques.size();i++) {
    if (MaxClique::computeWeight(cliques[i].nodes) < weight) {
      res = i;
    }
  }
  return res;
}


