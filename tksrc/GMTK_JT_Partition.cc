/*-
 * GMTK_JT_Partition.cc
 *     Junction Tree
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (C) 2009 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
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
#include "GMTK_GMParms.h"
#include "GMTK_JT_Partition.h"

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Static variables used by classes
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

#if HAVE_CONFIG_H
#include <config.h>
#endif
#if HAVE_HG_H
#include "hgstamp.h"
#endif
VCID(HGID)



////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Constructors/Destructors 
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

static void
cloneInterface(vector< set<RV*> > &destination_interface,
	       const vector< set<RV*> > &source_interface,
	       vector <RV*> &newRvs,
	       map < RVInfo::rvParent, unsigned > &ppf,
	       unsigned frame_delta)
{
  set<RV*>::iterator it;

  destination_interface.resize(source_interface.size());
  for (unsigned i=0; i < source_interface.size(); ++i) {
    destination_interface[i].clear();
    for (it =  source_interface[i].begin();
	 it != source_interface[i].end();
	 it++) 
    {
      RV* rv = (*it);
      RVInfo::rvParent rvp;
      rvp.first = rv->name();
      rvp.second = rv->frame()+frame_delta;    

      if ( ppf.find(rvp) == ppf.end() ) {
	coredump("INTERNAL ERROR: can't find rv %s(%d+%d)=%s(%d) in unrolled RV set, from from_liVars\n",
		 rv->name().c_str(),rv->frame(),frame_delta,
		 rvp.first.c_str(),rvp.second);
      }
      RV* nrv = newRvs[ppf[rvp]];
      destination_interface[i].insert(nrv);
    }
  }
}

static void 
cloneNodes(const set<RV*> &from_nodes, set<RV*> &dest_nodes, vector <RV*> &newRvs,
	   unsigned frame_delta, map < RVInfo::rvParent, unsigned > &ppf)
{
  set<RV*>::iterator it;
  for (it = from_nodes.begin();
       it != from_nodes.end();
       it++) 
  {
    RV* rv = (*it);
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frame_delta;    

    if ( ppf.find(rvp) == ppf.end() ) {
      coredump("INTERNAL ERROR: can't find rv %s(%d+%d)=%s(%d) in unrolled RV set, from from_part2\n",
	    rv->name().c_str(),rv->frame(),frame_delta,
	    rvp.first.c_str(),rvp.second);
    }
    RV* nrv = newRvs[ppf[rvp]];
    dest_nodes.insert(nrv);
  }
}


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
		       Section& from_part,
		       const unsigned int frameDelta,
		       // the left and right interface variables for
		       // this JT partition Empty if doesn't exist
		       // (say for an P or E partition). These have
		       // their own frame deltas since they might be
		       // different.
		       // Left interface:
		       const vector< set <RV*> > &from_liVars,
		       const unsigned int liFrameDelta,
		       // Right interface:
		       const vector< set <RV*> > &from_riVars,
		       const unsigned int riFrameDelta,
		       // Information todo the mapping.
		       vector <RV*>& newRvs,
		       map < RVInfo::rvParent, unsigned >& ppf)
{
  connected_components = from_part.connected_components;
  subtree_roots = from_part.subtree_roots;

  triMethod = from_part.triMethod;
  numVEseps = 0;

  cloneNodes(from_part.nodes, nodes, newRvs, frameDelta, ppf);
  cloneInterface(liNodes, from_liVars, newRvs, ppf, liFrameDelta);
  cloneInterface(riNodes, from_riVars, newRvs, ppf, riFrameDelta);

  section_li = from_part.section_li;
  section_ri = from_part.section_ri;
  ia_message_order = from_part.ia_message_order;
  ia_rev_msg_order = from_part.ia_rev_msg_order;
  clique_name_dictionary = from_part.clique_name_dictionary;

  // make the cliques.
  cliques.reserve(from_part.cliques.size());
  // 
  // NOTE: It is ***CRUCIAL*** for the cliques in the cloned partition
  // to be inserted in the *SAME ORDER* as in the partition being
  // cloned. If this is not done, inference will crash.

  set<RV*> empty;
  for (unsigned i=0;i<from_part.cliques.size();i++) {
    cliques.push_back(MaxClique(empty)); // length the vector
    new (&(cliques[i])) MaxClique(from_part.cliques[i], newRvs,ppf,frameDelta, true);
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
 Section& from_part,
 // Left interface:
 const vector< set <RV*> > &from_liVars,
 // Right interface:
 const vector< set <RV*> > &from_riVars
 )
{
  connected_components = from_part.connected_components;
  subtree_roots = from_part.subtree_roots;
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
JT_Partition::findLInterfaceClique(vector<unsigned> &liClique,bool& liCliqueSameAsInterface,const string priorityStr)
{
  liClique.clear();
  liCliqueSameAsInterface = true;
  for (unsigned i=0; i < liNodes.size(); ++i) {
    unsigned i_clique;
    bool iCliqueSameAsInterface;
    findInterfaceCliques(liNodes[i], i_clique, iCliqueSameAsInterface, priorityStr);
    liClique.push_back(i_clique);
    liCliqueSameAsInterface = liCliqueSameAsInterface && iCliqueSameAsInterface;
  }
}
void
JT_Partition::findRInterfaceClique(vector<unsigned> &riClique,bool& riCliqueSameAsInterface,const string priorityStr)
{
  riClique.clear();
  riCliqueSameAsInterface = true;
  for (unsigned i=0; i < riNodes.size(); ++i) {
    unsigned i_clique;
    bool iCliqueSameAsInterface;
    findInterfaceCliques(riNodes[i], i_clique, iCliqueSameAsInterface, priorityStr);
    riClique.push_back(i_clique);
    riCliqueSameAsInterface = riCliqueSameAsInterface && iCliqueSameAsInterface;
  }
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
    double wt = MaxClique::computeWeight(cliques[i].nodes);
    if (wt > weight) {
      weight = wt;
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
    double wt = MaxClique::computeWeight(cliques[i].nodes);
    if (wt < weight) {
      weight = wt;
      res = i;
    }
  }
  return res;
}


// for each connected component, find the clique with maximum weight
unsigned 
JT_Partition::cliqueWithMaxWeight(set<unsigned> const &subtree) {
  double weight = DBL_MIN;
  unsigned res = ~0x0;
  for (set<unsigned>::iterator it = subtree.begin(); it != subtree.end(); ++it) {
    double wt = MaxClique::computeWeight(cliques[*it].nodes);
    if (wt > weight) {
      weight = wt;
      res = *it;
    }
  }
  return res;
}

void 
JT_Partition::findSubtreeCliquesWithMaxWeight(vector<unsigned> &heaviest_cliques) {
  heaviest_cliques.clear();
  for (unsigned component=0; component < connected_components.size(); ++component) {
    set<unsigned> &subtree = connected_components[component];
    heaviest_cliques.push_back( cliqueWithMaxWeight(subtree) );
  }
}

void sv_intersect(set<unsigned> const &s, vector<unsigned> const &v, vector<unsigned> &intersection) {
  intersection.clear();
  for (unsigned i=0; i < v.size(); ++i) {
    if (s.find(v[i]) != s.end()) intersection.push_back(v[i]);
  }
}

void
JT_Partition::findSubtreeRoots(vector<unsigned> const &interface_cliques) {
  subtree_roots.clear();
  subtree_roots.resize(connected_components.size());

  for (unsigned i=0; i < connected_components.size(); ++i) {
    set<unsigned> const &subtree = connected_components[i];
    vector<unsigned> interface_cliques_in_subtree;
    sv_intersect(subtree, interface_cliques, interface_cliques_in_subtree);
    if (interface_cliques_in_subtree.size() == 0) {
      // This subtree has no interface cliques, so just pick the heaviest clique
      subtree_roots[i] = cliqueWithMaxWeight(subtree);
    } else {
      // This subtree has at least one interface clique. Pick the heaviest 
      // interface clique as the root of the subtree.
      unsigned subtree_root = ~0x0;
      double weight = DBL_MIN;
      for (unsigned j=0; j < interface_cliques_in_subtree.size(); ++j) {
	double wt = MaxClique::computeWeight(cliques[interface_cliques_in_subtree[j]].nodes);
	if (wt > weight) {
	  subtree_root = interface_cliques_in_subtree[j];
	  weight = wt;
	}
      }
      subtree_roots[i] = subtree_root;
    }
  }
}
