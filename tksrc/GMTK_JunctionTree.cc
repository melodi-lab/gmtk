/*-
 * GMTK_JunctionTree.cc
 *     Junction Tree
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2003, < fill in later >
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


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Static variables used by classes
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

VCID("$Header$");

// turns on or off VE separators, and determins the type to use PC, or PCG. 
unsigned JunctionTree::useVESeparators = (VESEP_PC | VESEP_PCG);
// turn off VE separators by default.
unsigned JunctionTree::veSeparatorWhere = 0;
// uncomment to turn on VE separators.
// unsigned JunctionTree::veSeparatorWhere = (VESEP_WHERE_P | VESEP_WHERE_C | VESEP_WHERE_E);
bool JunctionTree::jtWeightUpperBound = false;
bool JunctionTree::jtWeightMoreConservative = false;
float JunctionTree::jtWeightPenalizeUnassignedIterated = 0.0;
float JunctionTree::jtWeightSparseNodeSepScale = 1.0;
float JunctionTree::jtWeightDenseNodeSepScale = 1.0;
char* JunctionTree::priorityStr = "DSU";
char* JunctionTree::interfaceCliquePriorityStr = "W";

bool JunctionTree::probEvidenceTimeExpired = false;
bool JunctionTree::viterbiScore = false;

// default names of the three partitions for printing/debugging messages.
char* JunctionTree::P1_n = "P'";
char* JunctionTree::Co_n = "C'";
char* JunctionTree::E1_n = "E'";


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Misc Support
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

/*
 * Take the union of A and B and place it in C.
 */
static void
union_1_2_to_3(const set<RV*>& A,
	       const set<RV*>& B,
	       set<RV*>& C,
	       bool do_not_clear = false)
{
  if (!do_not_clear)
    C.clear();
  set_union(A.begin(),A.end(),
	    B.begin(),B.end(),
	    inserter(C,C.end()));
}

struct PairUnsigned1stElementCompare {  
  // sort ascending, comparing only 1st element.
  bool operator() (const pair<unsigned,unsigned>& a, 
		   const pair<unsigned,unsigned>& b) {
    // printf("comparing %d with %d\n",a.first,b.first);
    return (a.first) < (b.first);
  }
};


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        JT Partition support
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

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
      error("INTERNAL ERROR: can't find rv %s(%d+%d)=%s(%d) in unrolled RV set, from from_part2\n",
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
      error("INTERNAL ERROR: can't find rv %s(%d+%d)=%s(%d) in unrolled RV set, from from_liVars\n",
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
      error("INTERNAL ERROR: can't find rv %s(%d+%d)=%s(%d) in unrolled RV set, from from_riVars\n",
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
  double weight = -1.0;
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



////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        JT Inference Partition support
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////



JT_InferencePartition::JT_InferencePartition(JT_Partition& from_part,
					     vector <RV*>& newRvs,
					     map < RVInfo::rvParent, unsigned >& ppf,
					     const unsigned int frameDelta)
  : origin(from_part)
{

  // first allocate space with empty (and unusable) entries
  maxCliques.resize(origin.cliques.size());
  separatorCliques.resize(origin.separators.size());

  // then actually re-construct the objects in the array appropriately.
  for (unsigned i=0;i<maxCliques.size();i++) {
    new (&maxCliques[i]) InferenceMaxClique(origin.cliques[i],
					    newRvs,ppf,frameDelta);
  }
  for (unsigned i=0;i<separatorCliques.size();i++) {
    new (&separatorCliques[i]) InferenceSeparatorClique(origin.separators[i],
					    newRvs,ppf,frameDelta);
  }

}



/*-
 *-----------------------------------------------------------------------
 * JT_InferencePartition::emIncrement()
 *
 *    Go through each clique in the partition and update the assigned probabiltiy
 *    variables for all entries in each clique, based on global probability of evidence
 *    given as the argument.
 *    If 'localNormalization' == true, then we ignore the evidence provided by probE
 *    and sum the clique first to get the local normalization constant.
 *
 * See Also:
 *    0) JunctionTree::collectEvidence()
 *    1) JunctionTree::distributeEvidence()
 *    2) JunctionTree::emIncrement()
 *       
 *
 * Preconditions:
 *    The cliques must be an instantiated table. All data structures must be set up.
 *
 * Postconditions:
 *   The accumulators are increment accordingly.
 *
 * Side Effects:
 *   This will update the accumulators of all trainable parameter objects.
 *
 * Results:
 *   None
 *
 *-----------------------------------------------------------------------
 */
void
JT_InferencePartition::emIncrement(const logpr probE,
				   const bool localCliqueNormalization,
				   const double emTrainingBeam)
{
  for (unsigned cliqueNo=0;cliqueNo < maxCliques.size(); cliqueNo++ ) {
    maxCliques[cliqueNo].emIncrement(probE,localCliqueNormalization,emTrainingBeam);
  }
}


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Support for building a tree from clique graph
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////



/*-
 *-----------------------------------------------------------------------
 * JunctionTree::setUpDataStructures()
 *
 *   Sets up all the data structures in a JT (other than preparing for
 *   unrolling), and calls all the routines in the necessary
 *   order. This is like the main routine of this class, as it calls
 *   what needs to be done to set up a JT.
 *
 * Preconditions:
 *   Junction tree must have been created. Should not have called 
 *   setUpDataStructures() before. The C partitions in the template
 *   have variables, but the other partions (P, and E) might be empty.
 *
 * Postconditions:
 *   All data structures are set up. The JT is ready to have
 *   prepareForUnrolling() called.
 *
 * Side Effects:
 *   Changes many internal data structures in this object. 
 *
 * Results:
 *    None.
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::setUpDataStructures(const char* varPartitionAssignmentPrior,
				  const char *varCliqueAssignmentPrior)
{
  // main() routine for this class.
  createPartitionJunctionTrees(priorityStr);
  computePartitionInterfaces();
  createDirectedGraphOfCliques();
  assignRVsToCliques(varPartitionAssignmentPrior,varCliqueAssignmentPrior);
  computeUnassignedCliqueNodes();
  // TODO: move this next one above by one.
  setUpMessagePassingOrders();
  // create seps and VE seps.
  createSeparators(); 
  computeSeparatorIterationOrders();
  // the next routine also computes the dispositions, and a bunch
  // of other members of each MaxClique object.
  getCumulativeUnassignedIteratedNodes();
}




/*-
 *-----------------------------------------------------------------------
 * JunctionTree::createPartitionJunctionTree()
 *   Create a mini-junction tree from the cliques in the given partition.
 *   This uses Kruskal's greedy (but optimal) algorithm for MST generation.
 *
 *   TODO: move this routine to a MaxClique class at some point.
 *   TODO: do this during triangulation time.
 *
 * Preconditions:
 *   The partition must be instantiated with cliques 
 *
 * Postconditions:
 *   The cliques in the partition are now such that they
 *   form a junction tree over cliques within that partition.
 *
 * Side Effects:
 *   Modifies all neighbors variables within the cliques within the
 *   partition.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::createPartitionJunctionTree(Partition& part,const string priorityStr)
{
  const unsigned numMaxCliques = part.cliques.size();

  infoMsg(IM::Giga,"Starting create JT\n");

  if (numMaxCliques == 0) {
    // Nothing to do.
    // This could happen if the partition is empty which might occur
    // for empty P's and E's. C should never be empty.
    return;
  } else if (numMaxCliques == 1) {
    // then nothing to do
    infoMsg(IM::Giga,"Partition has only one clique\n");
    return;
  } else if (numMaxCliques == 2) {
    // then JT is easy, just connect the two cliques.
    infoMsg(IM::Giga,"Partition has only two cliques\n");
    part.cliques[0].neighbors.push_back(1);
    part.cliques[1].neighbors.push_back(0);
  } else {

    infoMsg(IM::Giga,"Partition has only %d cliques\n",numMaxCliques);

    // Run max spanning tree to construct JT from set of cliques.
    // This is basically Krusgal's algorithm, but without using the
    // fast data structures (it doesn't need to be that fast
    // since it is run one time per partition, for all inference
    // runs of any length.

    vector < set<unsigned> >  findSet;
    findSet.resize(numMaxCliques);
    for (unsigned i=0;i<numMaxCliques;i++) {
      set<unsigned> iset;
      iset.insert(i);
      findSet[i] = iset;
    }

    vector< Edge > edges;
    edges.reserve((numMaxCliques*(numMaxCliques-1))/2);
  
    for (unsigned i=0;i<numMaxCliques;i++) {
      for (unsigned j=i+1;j<numMaxCliques;j++) {
	set<RV*> sep_set;
	set_intersection(part.cliques[i].nodes.begin(),
			 part.cliques[i].nodes.end(),
			 part.cliques[j].nodes.begin(),
			 part.cliques[j].nodes.end(),
			 inserter(sep_set,sep_set.end()));
	Edge e;
	// define the edge
	e.clique1 = i; e.clique2 = j; 
	// !!!************************!!!
	// !!!** MUST DO THIS FIRST **!!!
	// !!!************************!!!
	// First push sep set size. To get a JT, we *MUST*
	// always choose from among the cliques that
	// have the largest intersection size.
	e.weights.push_back((double)sep_set.size());
	// now that the size is there, we have many other options.

	// All gets sorted in decreasing order, so that larger values
	// have priority. Thus, the remaining items we push back in
	// the case of ties.  
	//
	// *** Larger numbers are prefered. ***
	//
	// Options to include by priority of priorityStr:
	//   * D: number of deterministic nodes in separator 
	//   * E: number of deterministic nodes in union of cliques.
	//   * S: neg. weight of separator
	//   * U: neg. weight of union of cliques
	//   * V: neg. frame number variance in separator
	//   * W: neg. frame number variance in union
	//   * H: number of hidden nodes in separator
	//   * O: number of observed nodes in separator
	//   * L: number of hidden nodes in union
	//   * Q: number of observed nodes in union
	//
	//  Any can be preceeded by a '-' sign to flip effect. E.g., D-E-SU
	//
	// Default case: DSU
	//
	// TODO: Other ideas for this:
	//        a) maximize number of variables in same frame (or near each other) (like variance)
	//        b) minimize number of neighbors in each clique (i.e., 
	//           if cliques already have neighbors, choose the ones with fewer.
	//        c) integrate with RV value assignment to minimize
	//           the number of unassigned clique nodes (since they're
	//           iterated over w/o knowledge of any parents. If this
	//           ends up being a search, make this be offline, in with gmtkTriangulate
	// 

	float mult = 1.0;
	for (unsigned charNo=0;charNo< priorityStr.size(); charNo++) {
	  const char curCase = toupper(priorityStr[charNo]);
	  if (curCase == '-') {
	    mult = -1.0;
	    continue;
	  }

	  if (curCase == 'D') {

	    // push back number of deterministic nodes in
	    // the separator
	    set<RV*>::iterator it;
	    set<RV*>::iterator it_end = sep_set.end();
	    unsigned numDeterministicNodes = 0;
	    for (it = sep_set.begin(); it != it_end; it++) {
	      RV* rv = (*it);
	      if (rv->discrete() && RV2DRV(rv)->deterministic())
		numDeterministicNodes++;
	    }
	    e.weights.push_back(mult*(double)numDeterministicNodes);

	  } else if (curCase == 'E' || curCase == 'L' || curCase == 'Q' || curCase == 'W') {
	    // push back negative weight of two cliques together.
	    set<RV*> clique_union;
	    set_union(part.cliques[i].nodes.begin(),
		      part.cliques[i].nodes.end(),
		      part.cliques[j].nodes.begin(),
		      part.cliques[j].nodes.end(),
		      inserter(clique_union,clique_union.end()));

	    if (curCase == 'E') {
	      // push back number of deterministic nodes in
	      // the union
	      set<RV*>::iterator it;
	      set<RV*>::iterator it_end = clique_union.end();
	      unsigned numDeterministicNodes = 0;
	      for (it = clique_union.begin(); it != it_end; it++) {
		RV* rv = (*it);
		if (rv->discrete() && RV2DRV(rv)->deterministic())
		  numDeterministicNodes++;
	      }
	      e.weights.push_back(mult*(double)numDeterministicNodes);
	    } else if (curCase == 'L') {
	      set<RV*>::iterator it;
	      set<RV*>::iterator it_end = clique_union.end();
	      unsigned numHidden = 0;
	      for (it = clique_union.begin(); it != it_end; it++) {
		RV* rv = (*it);
		if (rv->hidden())
		  numHidden++;
	      }
	      e.weights.push_back(mult*(double)numHidden);
	    } else if (curCase == 'Q') {
	      set<RV*>::iterator it;
	      set<RV*>::iterator it_end = clique_union.end();
	      unsigned numObserved = 0;
	      for (it = clique_union.begin(); it != it_end; it++) {
		RV* rv = (*it);
		if (rv->observed())
		  numObserved++;
	      }
	      e.weights.push_back(mult*(double)numObserved);
	    } else if (curCase == 'W') {
	      // compute frame number variance in union, push back
	      // negative to prefer smalller frame variance (i.e.,
	      // connect things that on average are close to each
	      // other in time).
	      set<RV*>::iterator it;
	      set<RV*>::iterator it_end = clique_union.end();
	      double sum = 0;
	      double sumSq = 0;
	      for (it = clique_union.begin(); it != it_end; it++) {
		RV* rv = (*it);
		sum += rv->frame();
		sumSq += rv->frame()*rv->frame();
	      }
	      double invsize = 1.0/(double)clique_union.size();
	      double variance = invsize*(sumSq - sum*sum*invsize);
	      e.weights.push_back(mult*(double)-variance);
	    }
	  } else if (curCase == 'S') {

	    // push back negative weight of separator, to prefer
	    // least negative (smallest)  weight, since larger numbers
	    // are prefered.
	    e.weights.push_back(-(double)MaxClique::computeWeight(sep_set));

	    // printf("weight of clique %d = %f, %d = %f\n",
	    // i,part.cliques[i].weight(),
	    // j,part.cliques[j].weight());

	  } else if (curCase == 'U') {

	    // push back negative weight of two cliques together.
	    set<RV*> clique_union;
	    set_union(part.cliques[i].nodes.begin(),
		      part.cliques[i].nodes.end(),
		      part.cliques[j].nodes.begin(),
		      part.cliques[j].nodes.end(),
		      inserter(clique_union,clique_union.end()));
	    e.weights.push_back(-(double)MaxClique::computeWeight(clique_union));

	  } else if (curCase == 'V') {
	    // compute frame number variance in separator, push back
	    // negative to prefer smalller frame variance (i.e.,
	    // connect things that on average are close to each
	    // other in time).
	    set<RV*>::iterator it;
	    set<RV*>::iterator it_end = sep_set.end();
	    double sum = 0;
	    double sumSq = 0;
	    for (it = sep_set.begin(); it != it_end; it++) {
	      RV* rv = (*it);
	      sum += rv->frame();
	      sumSq += rv->frame()*rv->frame();
	    }
	    if (sep_set.size() == 0) {
	      e.weights.push_back(mult*(double)-FLT_MAX);
	    } else {
	      double invsize = 1.0/(double)sep_set.size();
	      double variance = invsize*(sumSq - sum*sum*invsize);
	      e.weights.push_back(mult*(double)-variance);
	    }
	  } else if (curCase == 'H') {
	    set<RV*>::iterator it;
	    set<RV*>::iterator it_end = sep_set.end();
	    unsigned numHidden = 0;
	    for (it = sep_set.begin(); it != it_end; it++) {
	      RV* rv = (*it);
	      if (rv->hidden())
		numHidden++;
	    }
	    e.weights.push_back(mult*(double)numHidden);
	  } else if (curCase == 'O') {
	    set<RV*>::iterator it;
	    set<RV*>::iterator it_end = sep_set.end();
	    unsigned numObserved = 0;
	    for (it = sep_set.begin(); it != it_end; it++) {
	      RV* rv = (*it);
	      if (rv->observed())
		numObserved++;
	    }
	    e.weights.push_back(mult*(double)numObserved);
	  } else {
	    error("ERROR: Unrecognized junction tree clique sort order letter '%c' in string '%s'\n",curCase,priorityStr.c_str());
	  }
	  mult = 1.0;
	}

	// add the edge.
	edges.push_back(e);
	if (IM::messageGlb(IM::Giga)) {
	  infoMsg(IM::Giga,"Edge (%d,%d) has sep size %.0f, ",
		  i,j,
		  e.weights[0]);
	  for (unsigned charNo=0;charNo< priorityStr.size(); charNo++) {
	    const char curCase = toupper(priorityStr[charNo]);
	    infoMsg(IM::Giga,"%c,weight[%d] = %f, ",curCase,charNo+1,e.weights[charNo+1]);
	  }
	  infoMsg(IM::Giga,"\n");
	}
      }
    }
  
    // sort in decreasing order by edge weight which in this
    // case is the sep-set size.
    sort(edges.begin(),edges.end(),EdgeCompare());

    unsigned joinsPlusOne = 1;
    for (unsigned i=0;i<edges.size();i++) {
      infoMsg(IM::Giga,"Edge %d has sep size %.0f\n",
	      i,
	      edges[i].weights[0]);

      set<unsigned>& iset1 = findSet[edges[i].clique1];
      set<unsigned>& iset2 = findSet[edges[i].clique2];
      if (iset1 != iset2) {
	// merge the two sets
	set<unsigned> new_set;
	set_union(iset1.begin(),iset1.end(),
		  iset2.begin(),iset2.end(),
		  inserter(new_set,new_set.end()));
	// make sure that all members of the set point to the
	// new set.
	set<unsigned>::iterator ns_iter;
	for (ns_iter = new_set.begin(); ns_iter != new_set.end(); ns_iter ++) {
	  const unsigned clique = *ns_iter;
	  findSet[clique] = new_set;
	}
	infoMsg(IM::Giga,"Joining cliques %d and %d (edge %d) with intersection size %.0f\n",
		edges[i].clique1,edges[i].clique2,i,edges[i].weights[0]);

	if (edges[i].weights[0] == 0.0) {
	  if (IM::messageGlb(IM::High)) {
	    // there is no way to know the difference here if the
	    // graph is non-triangualted or is simply disconnected
	    // (which is ok). A non-triangulated graph might have
	    // resulted from the user editing the trifile, but we
	    // presume that MCS has already checked for this when
	    // reading in the trifiles. We just issue an informative
	    // message just in case.
	    printf("NOTE: junction tree creation joining two cliques (%d and %d) with size 0 set intersection. Either disconnected (which is ok) or non-triangulated (which is bad) graph.\n",
		    edges[i].clique1,edges[i].clique2);
	    // TODO: print out two cliques that are trying to be joined.
	    printf("Clique %d: ",edges[i].clique1);
	    part.cliques[edges[i].clique1].printCliqueNodes(stdout);
	    printf("Clique %d: ",edges[i].clique2);
	    part.cliques[edges[i].clique2].printCliqueNodes(stdout);
	  }
	}

 	part.cliques[edges[i].clique1].neighbors.push_back(edges[i].clique2);
	part.cliques[edges[i].clique2].neighbors.push_back(edges[i].clique1);

	if (++joinsPlusOne == numMaxCliques)
	  break;
      }
    }
  }
}


/*-
 *-----------------------------------------------------------------------
 * JunctionTree::base_unroll()
 *   sets up the internal data structures and creates a junction tree corresponding
 *   to the graph with C' unrolled k=2 times, (i.e., we have
 *   a graph that has (3) copies of C'). It does this using the initial
 *   partitions given in gm_template (P', C', and E').
 *
 *   This routine should really only be called once to set up
 *   the base partitions (P1,C1,C2,C3,E1), from which all future
 *   unrolls will take place.
 *
 *   Sets things up so that inference can take place.
 *   Unrolling is in terms of new C' partitions (so that if S=1, k corresponds to frames)
 *   but in general the network will be T(P') + k*T(C') + T(E) frames long, and in general
 *   T(C') >= S. This also depends on if the boundary algorithm was run or not. 
 *
 * Preconditions:
 *   The partitions must be validly instantiated with cliques 
 *
 * Postconditions:
 *   The partitions are such that they now know who they are connected to on
 *   the right and left.
 *
 * Side Effects:
 *   Modifies all neighbors variables within the cliques within the
 *   partition.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::base_unroll()
{

  // first create the unrolled set of random variables corresponding
  // to this JT.

  // unrolled random variables
  vector <RV*> unrolled_rvs;
  // mapping from name(frame) to integer index into unrolled_rvs.
  map < RVInfo::rvParent, unsigned > ppf;
  // number of C repetitions is M + (k+1)*S, so
  // we unroll one less than that.
  fp.unroll(gm_template.M + gm_template.S - 1,
	    unrolled_rvs,ppf);
  
  set <RV*> empty;

  // copy P partition 
  new (&P1) JT_Partition(gm_template.P,0*gm_template.S,
			 empty,0*gm_template.S,
			 gm_template.PCInterface_in_P,0*gm_template.S,
			 unrolled_rvs,ppf);

  // copy E partition
  new (&E1) JT_Partition(gm_template.E,0*gm_template.S,
			 gm_template.CEInterface_in_E,0*gm_template.S,
			 empty,0*gm_template.S,
			 unrolled_rvs,ppf);

  if (gm_template.leftInterface) {
    // left interface case
    // Template has a (P)(C)(CE) = (P')(C')(E')

    // An E' in the LI case must be E' = [ C E ], so there is never an
    // empty E here (since C is never empty).
    assert ( E1.cliques.size() > 0 );

    // there might be an empty P1, though, so we need to check for that.
    if (P1.cliques.size() > 0) {
      // neither P1 nor E1 are empty.
      new (&Co) JT_Partition(gm_template.C,0*gm_template.S,
			     gm_template.PCInterface_in_C,0*gm_template.S,			   
			     gm_template.CEInterface_in_C,0*gm_template.S,
			     unrolled_rvs,ppf);
    } else {
      // P1 is empty. For Co's left interface, we use its right
      // interface since there is no interface to P1. The reason why
      // this is valid is that E' = [C E], and E can't connect to C1
      // in [C1 C2 E], so C2's li to C1 is the same as E's li to C' in
      // [C' E']. We do need to shift E'ls li to C' to the left by S
      // though.
      new (&Co) JT_Partition(gm_template.C,0*gm_template.S,
			     gm_template.CEInterface_in_C,-1*gm_template.S,			   
			     gm_template.CEInterface_in_C,0*gm_template.S,
			     unrolled_rvs,ppf);
    }
  } else {
    // right interface case, symmetric to the left interface case above.
    // Right interface case.
    // Template has a (PC)(C)(E) = (P')(C')(E')

    // An P' in the RI case must be P' = [ P C ], so
    // there is never an empty P' here.
    assert ( P1.cliques.size() > 0 );
    
    // there might be an empty E1, though, so we need to check for that.

    if (E1.cliques.size() > 0) {
      // neither E1 nor P1 are empty.
      new (&Co) JT_Partition(gm_template.C,0*gm_template.S,
			     gm_template.PCInterface_in_C,0*gm_template.S,			   
			     gm_template.CEInterface_in_C,0*gm_template.S,
			     unrolled_rvs,ppf);
    } else {
      new (&Co) JT_Partition(gm_template.C,0*gm_template.S,
			     gm_template.PCInterface_in_C,0*gm_template.S,			   
			     gm_template.PCInterface_in_C,1*gm_template.S,
			     unrolled_rvs,ppf);
    }
  }
  
}


/*-
 *-----------------------------------------------------------------------
 * JunctionTree::createPartitionInterfaces()
 *   This version finds the vairous "interface" between the partitions
 *   in a graph, but does so for all possible unrollings by doing
 *   one unrolling. It fills in the P_ri_to_C, C_li_to_P variables
 *   and is thus necessary to be called before any true junction tree
 *   creation.
 *
 *
 * Preconditions:
 *   The gm_template partitions must be validly instantiated with cliques.
 *
 * Postconditions:
 *   The member variables P_ri_to_C etc. are filled in.
 *
 * Side Effects:
 *   chnages group of member variables P_ri_to_C, C_li_to_P, etc.
 *
 * Results:
 *   Modifies member variables as above.
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::computePartitionInterfaces()
{
  // set up the base partitions
  base_unroll();

  // Use base partitions to find the various interface cliques.
  bool P_riCliqueSameAsInterface;
  bool E_liCliqueSameAsInterface;
  P1.findRInterfaceClique(P_ri_to_C,P_riCliqueSameAsInterface,interfaceCliquePriorityStr);
  E1.findLInterfaceClique(E_li_to_C,E_liCliqueSameAsInterface,interfaceCliquePriorityStr);

  // Note that here, we do the same for both left and right interface,
  // but we break it out into two cases for clarity.
  if (gm_template.leftInterface) {
    // left interface case

    bool Co_liCliqueSameAsInterface;
    bool Co_riCliqueSameAsInterface;

    Co.findLInterfaceClique(C_li_to_C,Co_liCliqueSameAsInterface,interfaceCliquePriorityStr);
    C_li_to_P = C_li_to_C;
    Co.findRInterfaceClique(C_ri_to_E,Co_riCliqueSameAsInterface,interfaceCliquePriorityStr);
    C_ri_to_C = C_ri_to_E;

    // sanity check
    assert (C_ri_to_C == C_ri_to_E);

    P_to_C_icliques_same =
      P_riCliqueSameAsInterface && Co_liCliqueSameAsInterface;

    C_to_C_icliques_same =
      Co_riCliqueSameAsInterface && Co_liCliqueSameAsInterface;    

    C_to_E_icliques_same =
      Co_riCliqueSameAsInterface && E_liCliqueSameAsInterface;

  } else { // right interface
    // left interface case

    bool Co_liCliqueSameAsInterface;
    bool Co_riCliqueSameAsInterface;

    Co.findLInterfaceClique(C_li_to_P,Co_liCliqueSameAsInterface,interfaceCliquePriorityStr);
    Co.findRInterfaceClique(C_ri_to_C,Co_riCliqueSameAsInterface,interfaceCliquePriorityStr);
    C_li_to_C = C_li_to_P;
    C_ri_to_E = C_ri_to_C;

    // sanity check
    assert (C_li_to_C == C_li_to_P);
    
    P_to_C_icliques_same =
      P_riCliqueSameAsInterface && Co_liCliqueSameAsInterface;

    C_to_C_icliques_same =
      Co_riCliqueSameAsInterface && Co_liCliqueSameAsInterface;    

    C_to_E_icliques_same =
      Co_riCliqueSameAsInterface && E_liCliqueSameAsInterface;

  }

  // E order, clique 0 is choosen as root arbitrarily for now.  TODO:
  // see if it is possible to choose a better root for E.  Perhaps use
  // the maximum weight clique.
  E_root_clique = 0;
  // E_root_clique = E1.cliqueWithMinWeight();
  // If this is updated, need also to update in all other places
  // the code, search for string "update E_root_clique"
  
}


/*-
 *-----------------------------------------------------------------------
 * JunctionTree::createPartitionInterface()
 *   Finds the "interface" between the two partitions. The interface are
 *   the two cliques in each interface that have maximum intersection.
 *   These are the two cliques that are joined in a junction tree.
 *   It is assumed that partition 1 is just to the "left" of partition 2.
 *   This does an N*M algorithm, where N (resp. M) is the number of cliques 
 *   in one (resp. the other) partition
 *
 *
 * Preconditions:
 *   The partitions must be validly instantiated with cliques, and part 1
 *   must be just to the left of part 2 (i.e., valid combinations
 *   include (P,C), (C,E), or (C1,C2). 
 *
 * Postconditions:
 *   The partitions are such that they now know who they are connected to on
 *   the right and left.
 *
 * Side Effects:
 *   none
 *
 * Results:
 *   The left partitions right interface clique and the right partitions
 *   left interface clique are returned via the part1_ric and part2_lic variables.
 *
 *-----------------------------------------------------------------------
 */

void
JunctionTree::computePartitionInterface(JT_Partition& part1,
					// partitition 1's right interface clique
					unsigned int& part1_ric,
					JT_Partition& part2,
					// partitition 2's left interface clique
					unsigned int& part2_lic,
					// are the two interface cliques identical
					bool& icliques_same)
{
  const unsigned numMaxCliques1 = part1.cliques.size();
  const unsigned numMaxCliques2 = part2.cliques.size();

  // run a dumb n^2 algorithm to find the two cliques in each
  // partition that have maximum overlap. Because the random variable
  // set is not guaranteed to be the same at this point, we need to
  // compute intersection by (name,frame) equivalance.

  infoMsg(IM::Giga,"Starting create partition interface\n");
  
  unsigned max_sep_set_size = 0;
  unsigned part1_clique=0,part2_clique=0;
  unsigned part1_clique_size=0,part2_clique_size=0;
  for (unsigned i=0;i<numMaxCliques1;i++) {
    for (unsigned j=0;j<numMaxCliques2;j++) {

      set< RV*>::iterator it1;
      const set< RV*>::iterator it1_end = part1.cliques[i].nodes.end();
      set< RV*>::iterator it2;
      const set< RV*>::iterator it2_end = part2.cliques[j].nodes.end();

      unsigned sep_set_size = 0;
      for (it1 = part1.cliques[i].nodes.begin(); it1 != it1_end ; it1++) {
	for (it2 = part2.cliques[j].nodes.begin(); it2 != it2_end ; it2++) {
	  RV* rv1 = (*it1);
	  RV* rv2 = (*it2);
	  infoMsg(IM::Giga,"comparing %s(%d) with %s(%d)\n",
		  rv1->name().c_str(),rv1->frame(),		  
		  rv2->name().c_str(),rv2->frame());
	  if (rv1->frame() == rv2->frame() && rv1->name() == rv2->name()) {
	    sep_set_size++;
	  }
	}
      }
      infoMsg(IM::Giga,"clique %d of part1 and clique %d of part2 has overlap %d\n",i,j,sep_set_size);
      if (sep_set_size > max_sep_set_size) {
	part1_clique = i;
	part1_clique_size = part1.cliques[i].nodes.size();
	part2_clique = j;
	part2_clique_size = part2.cliques[j].nodes.size();
	max_sep_set_size = sep_set_size;
      }
    }
  }
  infoMsg(IM::Giga,"clique %d of part1 (size %d) and clique %d of part2 (size %d) has max overlap of %d\n",
	  part1_clique,part1_clique_size,
	  part2_clique,part2_clique_size,
	  max_sep_set_size);
  part1_ric = part1_clique;
  part2_lic = part2_clique;
  if ((part1_clique_size == part2_clique_size)
      &&
      (part1_clique_size == max_sep_set_size)) {
    // then the interface cliques on both sides of the
    // partition are the same, so one can be removed.
    icliques_same = true;
  } else {
    icliques_same = false;
  }
}


/*-
 *-----------------------------------------------------------------------
 * JunctionTree::createDirectedGraphOfCliques()
 *   Turns undirected graph into DAG of cliques, for all partitions.
 * 
 *
 * Preconditions:
 *   Partitions must be instantiated, and interface cliques
 *   must have been computed (i.e., computePartitionInterfaces()
 *   must have been called)
 *
 * Postconditions:
 *   The member variables P_to_C_message_order (etc.) are created
 *   with the order of message passing.
 *
 * Side Effects:
 *    Changes member variables as described above.
 *
 * Results:
 *    None.
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::createDirectedGraphOfCliques()
{

  createDirectedGraphOfCliques(P1,
			       P_ri_to_C);
  if (gm_template.leftInterface) {
    createDirectedGraphOfCliques(Co,
				 C_ri_to_E);
  } else { // right interface
    createDirectedGraphOfCliques(Co,
				 C_ri_to_C);
  }
  createDirectedGraphOfCliques(E1,
			       E_root_clique);

}


/*-
 *-----------------------------------------------------------------------
 * JunctionTree::createDirectedGraphOfCliques()
 * JunctionTree::createDirectedGraphOfCliquesRecurse()
 *   From the root that is given as an argument, turn the
 *   undirected graph of cliques JT into a directed graph
 *   of cliques rooted at root. Do this by assigning
 *   to each clique's children member variable.
 *
 * Preconditions:
 *   Partitions must be instantiated, and interface cliques
 *   must have been computed (i.e., computePartitionInterfaces()
 *   must have been called)
 *
 * Postconditions:
 *   children member variables in cliques are set.
 *
 * Side Effects:
 *    none
 *
 * Results:
 *    None, but passed back in arguments.
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::createDirectedGraphOfCliques(JT_Partition& part,
					   const unsigned root)
{
  // check for empty partitions, possible for P or E.
  if (part.cliques.size() == 0)
    return;

  // make sure none have been so far visited
  vector< bool > visited(part.cliques.size());
  for (unsigned i=0;i<part.cliques.size();i++) {
    visited[i] = false;
    part.cliques[i].children.clear();
  }
  createDirectedGraphOfCliquesRecurse(part,root,visited);
}
void
JunctionTree::createDirectedGraphOfCliquesRecurse(JT_Partition& part,
						  const unsigned root,
						  vector< bool >& visited)
{
  visited[root] = true;
  for (unsigned childNo=0;
       childNo<part.cliques[root].neighbors.size();
       childNo++) {
    const unsigned child = part.cliques[root].neighbors[childNo];
    if (visited[child])
      continue;
    part.cliques[root].children.push_back(child);
    createDirectedGraphOfCliquesRecurse(part,child,visited);
  }
}




/*-
 *-----------------------------------------------------------------------
 * JunctionTree::assignRVsToCliques()
 *    Assign to the cliques the random variables that will
 *    be used to generate the tables for each clique. This is
 *    mainly a stub routine to operate on all partitions. The
 *    real code is in the other routines below.
 *
 * Preconditions:
 *     createPartitionJunctionTrees(), computePartitionInterfaces(), and
 *     createDirectedGraphOfCliques() must have been called before 
 *     this can be called.
 *
 * Postconditions:
 *     Each of the partitions so specified in the code below
 *     are now such that their cliques have their 'assignedNodes' 
 *     member filled in.
 *
 * Side Effects:
 *     Will change the clique data member variables
 *
 * Results:
 *     None.
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::assignRVsToCliques(const char* varPartitionAssignmentPrior,
				 const char *varCliqueAssignmentPrior)
{



  if (P1.nodes.size() > 0) {
    infoMsg(IM::Giga,"assigning rvs to P1 partition\n");
    assignRVsToCliques(P1_n,P1,P_ri_to_C,varPartitionAssignmentPrior,varCliqueAssignmentPrior);
    // accumulate what occured in P1 so Co can use it. 
    union_1_2_to_3(P1.cliques[P_ri_to_C].cumulativeAssignedNodes,
		   P1.cliques[P_ri_to_C].assignedNodes,
		   Co.cliques[C_li_to_P].cumulativeAssignedNodes);
    union_1_2_to_3(P1.cliques[P_ri_to_C].cumulativeAssignedProbNodes,
		   P1.cliques[P_ri_to_C].assignedProbNodes,
		   Co.cliques[C_li_to_P].cumulativeAssignedProbNodes);

    // printf("P1's cum ass prob:");
    // printRVSet(stdout,P1.cliques[P_ri_to_C].cumulativeAssignedProbNodes);
    // printf("P1's ass prob:");
    // printRVSet(stdout,P1.cliques[P_ri_to_C].assignedProbNodes);
    // printf("Co's cum ass prob:");
    // printRVSet(stdout,Co.cliques[C_li_to_P].cumulativeAssignedProbNodes);

  }
  infoMsg(IM::Max,"assigning rvs to Co partition\n");
  assignRVsToCliques(Co_n,Co,C_ri_to_C,varPartitionAssignmentPrior,varCliqueAssignmentPrior);

  if (E1.nodes.size() > 0) {
    infoMsg(IM::Max,"assigning rvs to E partition\n");
    // accumulate what occured in P1,Co so E1 can use it. 
    union_1_2_to_3(Co.cliques[C_ri_to_E].cumulativeAssignedNodes,
		   Co.cliques[C_ri_to_E].assignedNodes,
		   E1.cliques[E_li_to_C].cumulativeAssignedNodes);
    union_1_2_to_3(Co.cliques[C_ri_to_E].cumulativeAssignedProbNodes,
		   Co.cliques[C_ri_to_E].assignedProbNodes,
		   E1.cliques[E_li_to_C].cumulativeAssignedProbNodes);
    assignRVsToCliques(E1_n,E1,E_root_clique,varPartitionAssignmentPrior,varCliqueAssignmentPrior);
  }

  // lastly, check to make sure all nodes have been assigned to to
  // give probability to one clique. This should be done in case the
  // user edited the trifile, which is one possible reason for this to
  // end in error.
  set <RV*> allNodes;
  union_1_2_to_3(P1.nodes,Co.nodes,allNodes);
  union_1_2_to_3(Co.nodes,E1.nodes,allNodes,true);
  set <RV*> allAssignedProbNodes;
  if (E1.cliques.size() > 0)
    union_1_2_to_3(E1.cliques[E_root_clique].cumulativeAssignedProbNodes,
		   E1.cliques[E_root_clique].assignedProbNodes,
		   allAssignedProbNodes);
  else 
    union_1_2_to_3(Co.cliques[C_ri_to_E].cumulativeAssignedProbNodes,
		   Co.cliques[C_ri_to_E].assignedProbNodes,
		   allAssignedProbNodes);

  set <RV*> nodesThatGiveNoProb;
  set_difference(allNodes.begin(),allNodes.end(),
		 allAssignedProbNodes.begin(),allAssignedProbNodes.end(),
		 inserter(nodesThatGiveNoProb,
			  nodesThatGiveNoProb.end()));
  
  if (nodesThatGiveNoProb.size() > 0) {
    fprintf(stderr,"INTERNAL ERROR: some nodes do not give probability to any clique, those nodes include: ");
    printRVSet(stderr,nodesThatGiveNoProb);
    error("Possibly corrupt trifile. Exiting program ...");
  }


}




/*-
 *-----------------------------------------------------------------------
 * JunctionTree::assignRVsToCliques()
 *    Assign to the cliques the random variables that will
 *    be used to generate the tables for each clique, and do
 *    it for a particular partition 'part' starting at
 *    the clique number given by rootClique.
 * 
 *    The assumption is that collectEvidence will be called with this
 *    rootClique, so the assignment of variables to cliques is an attempt
 *    to assign parents before children, and assign parents farther away
 *    from the root than children, so that by the time children are encountered,
 *    the parent values are already enumerated.
 *
 *    The basic algorithm is recursive, and uses a helper routine for
 *    that purpose. It descends down the graph in depth first order,
 *    attempting to assign a random variable to a clique. Note that
 *    it is possible for a rv to not be assigned (i.e., if no clique
 *    exist that contains all the rvs parents). If this happens, it
 *    is not an error, rather it should be the case that the rv lives with
 *    its parents in another partition.
 *
 * Preconditions:
 *     createPartitionJunctionTrees(), computePartitionInterfaces() must
 *     have been called. The partition must not be empty.
 *
 * Postconditions:
 *     the partition so specified in the code below
 *     are now such that their cliques have their 'assignedNodes' 
 *     member filled in.
 *
 * Side Effects:
 *     Will change the clique data member variables
 *
 * Results:
 *     None.
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::assignRVsToCliques(const char *const partName,
				 JT_Partition& part,
				 const unsigned rootClique,
				 const char* varPartitionAssignmentPrior,
				 const char *varCliqueAssignmentPrior)
{
  vector<RV*> sortedNodes;

  // We use a constrained topological sort based on command line
  // arguments, so that variables are considered in what hopefully
  // will be a good order (e.g., might have it such that continuous
  // variables or discrete observations come as early as possible in
  // the ordering).  This will also be beneficial when sampling hidden
  // continuous nodes.
  infoMsg(IM::Giga,"Sorting partition variables using priority order (%s)\n",
	  (varPartitionAssignmentPrior?varPartitionAssignmentPrior:"NULL"));

  GraphicalModel::topologicalSortWPriority(part.nodes,part.nodes,sortedNodes,varPartitionAssignmentPrior);

  // printf("have %d sorted nodes and %d cliques\n",sortedNodes.size(),part.cliques.size());

  // update the cumulative RV assignments before we begin.
  getCumulativeAssignedNodes(part,rootClique);

  for (unsigned n=0;n<sortedNodes.size();n++) {

    RV* rv = sortedNodes[n];

    // precompute the parents set for set intersection operations.
    // The one stored in the rv is a vector, but we need a set, so we
    // pre-compute it here. While we're at it, we compute if all
    // parents are observed since the behaviour is a bit
    // different in this case.
    set<RV*> parSet;
    bool allParentsObserved=true;
    for (unsigned p=0;p<rv->allParents.size();p++) {
      // TODO: see if there is a faster STL way to go from vector to set.
      // fprintf(stderr,"about to insert rv with address %X\n",(void*)rv->allPossibleParents[p]);
      parSet.insert(rv->allParents[p]);
      if (rv->allParents[p]->hidden())
	allParentsObserved = false;
    }

    // create a map to keep track of the scores of a given clique in
    // which RV should produce probability.  This utilizes the fact
    // that the multimap stores values in *ascending* order based on
    // key (here the score), and so the first one (i.e.,
    // scoreSet.begin()) should have the lowest weight.
    multimap < vector<double>, unsigned> scoreSet;

    unsigned numberOfTimesAssigned = 0;
    assignRVToClique(partName,
		     part,rootClique,
		     0,
		     rv,
		     numberOfTimesAssigned,
		     parSet,allParentsObserved,
		     scoreSet);
    infoMsg(IM::Max,
	    "Part %s: random variable %s(%d) with its parents contained in %d cliques\n",
	    partName,
	    rv->name().c_str(),rv->frame(),numberOfTimesAssigned);

    if (numberOfTimesAssigned > 0) {
      // then it was assigned in this partition at least once, but we
      // still need to check if it should give probability to some
      // clique within this partition. We need to check both if the
      // scoreSet has a value, and also if it was not assigned, since
      // the scoreSet might have obtained a value from a branch of the
      // JT below which didn't have this rv assigned with probability.

      if ((scoreSet.size() > 0)  &&
	  (part.cliques[rootClique].cumulativeAssignedProbNodes.find(rv) == 
	   part.cliques[rootClique].cumulativeAssignedProbNodes.end())) {

	// Then the node will be a probability contributer to one clique
	// in this partition.
	
	// We choose one of those cliques to be the one the node
	// contributes probabilty to.
	unsigned clique_num = (*(scoreSet.begin())).second;
	part.cliques[clique_num].assignedProbNodes.insert(rv);
	if (IM::messageGlb(IM::Max)) {
	  infoMsg(IM::Max,
		  "Part %s: random variable %s(%d) giving probability to clique %d.\n",
		  partName,
		  rv->name().c_str(),rv->frame(),clique_num);
	}
      }
      
      // update the cumulative RV assignments.
      getCumulativeAssignedNodes(part,rootClique);

    } else {

      // rv was not assigned to this partition, it must be the case
      // that it will be assigned to a different partition. This
      // could come from the partition before or the partition after
      // the current partition. If there are only forward-directed
      // arrows, it will come from the previous partition (and vice
      // versa if there are only backwards going edges). If we have
      // both directed edges, it could be assigned in either left or
      // right partition.

      // Note that it is impossible for a node to be assigned to give
      // probabiltiy in two partitions. The reason is the
      // following. The only nodes that live in both partitions are
      // the interface nodes.  These nodes have been special cased
      // above using the 'alreadyAProbContributer' variable, and
      // since we keep cumulative track of nodes that are to
      // provide probability, we won't do this more than once
      // across partitions.

      infoMsg(IM::Max,"Part %s: random variable %s(%d) not assigned in current partition\n",partName,
	      rv->name().c_str(),rv->frame());
      // in any event, keep track of this node, perhaps useful for jtWeight.
      part.unassignedInPartition.insert(rv);
    }
  }

  // Next, we compute the order that the assigned nodes in this clique
  // are iterated. This can have a dramatic effect on performance in
  // some cases. In general, we want to iterate over large cardinality
  // nodes last, since that results in the fewest branch mis-predicts.
  if (varCliqueAssignmentPrior && strlen(varCliqueAssignmentPrior) > 0) {
    infoMsg(IM::Giga,"Sorting cliques variables using priority order (%s)\n",varCliqueAssignmentPrior);
    
    // Lastly, we re-sort the original sorted nodes in each clique
    // according to a designated (and hopefully good) topological
    // order.
    for (unsigned clique_num=0;clique_num<part.cliques.size(); clique_num++ ) {
      GraphicalModel::topologicalSortWPriority(part.cliques[clique_num].assignedNodes,
					       part.cliques[clique_num].assignedNodes,
					       sortedNodes,
					       varCliqueAssignmentPrior);
      part.cliques[clique_num].sortedAssignedNodes = sortedNodes;
      if (IM::messageGlb(IM::Max)) {
	printf("Clique %d variables after sort:",clique_num);
	printRVSet(stdout,part.cliques[clique_num].sortedAssignedNodes,true);
      }
    }
  }

}






/*-
 *-----------------------------------------------------------------------
 * JunctionTree::assignRVToClique()
 *    Helper routine for assignRVsToCliques() above.
 *
 * Preconditions:
 *     Must be called only byy assignRVsToCliques().
 *
 * Postconditions:
 *     the random variable might be assigned (but the assignment
 *     might also be done by the caller).
 *
 * Side Effects:
 *     Might change the clique data member variables
 *
 * Results:
 *     None.
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::assignRVToClique(const char *const partName,
			       JT_Partition& part,
			       const unsigned root,
			       const unsigned depth,
			       RV* rv,
			       unsigned& numberOfTimesAssigned,
			       set<RV*>& parSet,
			       const bool allParentsObserved,
			       multimap < vector<double>, unsigned >& scoreSet)
{
  // keep a reference for easy access
  MaxClique& curClique = part.cliques[root];

  // First, check if directed_closure(rv) == (rv U parents(rv)) not in
  // the current clique, and if so then continue on down.
  bool closure_in_clique = (curClique.nodes.find(rv) != curClique.nodes.end());
  
  // printf("clique%d: from rv, closure_in_clique = %d\n",root,closure_in_clique);

  if (closure_in_clique && !allParentsObserved) {
    // need to further check that parents are in clique
    set<RV*> parentsInClique;
    set_intersection(curClique.nodes.begin(),
		     curClique.nodes.end(),
		     parSet.begin(),parSet.end(),
		     inserter(parentsInClique,parentsInClique.end()));

    // First check for equal size. If the size is equal, then *all*
    // parents are in the clique.
    closure_in_clique = (parentsInClique.size() == parSet.size());
    if (!closure_in_clique) {
      // rather than only checking for equal size, we here also check
      // if no *hidden* parents are out of the clique, since that is
      // all we need to assign this node to this clique. In other
      // words, if all parents that are out of clique are observed, we
      // can still go. In fact, assuming moralization doesn't moralize
      // w.r.t. an observed parent, this step is necessary sometimes
      // to get the child assigned.

      assert ( parentsInClique.size() < parSet.size() );

      set<RV*> parentsOutOfClique;
      set_difference(parSet.begin(),parSet.end(),
		     parentsInClique.begin(),parentsInClique.end(),
		     inserter(parentsOutOfClique,parentsOutOfClique.end()));
      
      // Check to see if all parents out of the clique are observed,
      // and if so, the closure is in the clique.
      set<RV*>::iterator it;
      closure_in_clique = true;
      for (it = parentsOutOfClique.begin(); it != parentsOutOfClique.end(); it++) {
	RV* rv = (*it);
	if (rv->hidden()) {
	  closure_in_clique = false;
	  break;
	}
      }
    }
  }

  // printf("clique%d: from par(rv), closure_not_in_clique = %d, num
  // children =
  // %d\n",root,closure_not_in_clique,part.cliques[root].children.size());

  if (closure_in_clique) {
    // So closure(rv) is in current clique, meaning rv and *all* its
    // hidden parents are members of this clique, but it is not
    // nec. the case that all of rv's parents are themselves assigned
    // to produce probabilities in this clique.

    // Since closure is in the clique, we "assign" this node to this
    // clique. Note, this does not necessarily mean that the node will
    // contribute probabilties to this clique's potential function,
    // only that we will use the RV iterator in some cases (sparse,
    // deterministic, sampling) to iterate over this RV after its
    // parents have been assigned (as long as the RV doesn't come in
    // as a separator that is).
    curClique.assignedNodes.insert(rv);
    curClique.sortedAssignedNodes.push_back(rv);
    numberOfTimesAssigned++;

    infoMsg(IM::Max,
	    "Part %s: RV and its parents in clique, assigning random variable %s(%d) to clique %d\n",
	    partName,
	    rv->name().c_str(),rv->frame(),root);

    // We don't make the probability assignment decision until after
    // we consider all possible such candidate cliques.  But first,
    // check if already set up to give probability to some other
    // clique.
    if (curClique.cumulativeAssignedProbNodes.find(rv) == curClique.cumulativeAssignedProbNodes.end()) {
      // So RV is not yet giving probability anywhere. We deal with
      // the score stuff to find a good clique for this variable to
      // assign probability.

      infoMsg(IM::Max,
	      "Part %s: considering clique %d for random variable %s(%d) probability contribution\n",
	      partName,
	      root,
	      rv->name().c_str(),rv->frame());

      // In this discussion that follows, we make a distinction
      // between node parents and children (i.e., original graph
      // parents and children) and JT parents and children which are
      // cliques in the JT. A child clique c of a parent p in the JT
      // is one such that c ~ p and that depth(c) = depth(p)+1, where
      // depth(root)=0 in the junction tree.

      // The goal here is to figure out in which of a set of possible
      // cliques the node should be producing probabilities.

      // TODO: make these various options available in priority order
      // to the command line, just like the topologically-constrainted
      // variable within-clique iteration order, currently available
      // via -vcap.


      // Since the multimap stores values in *ascending* order based
      // on key, lower numbers are prefered.

      // Compute (and ultimately push back) items in decreasing order of
      // priority.  Lower numbers is better (e.g., more negative or less
      // positive is higher priority).  First thing inserted has highest
      // priority.  At some time later, the rv will be assigned to the
      // clique with the *LOWEST* score.

      // What we do, is continue on down. It may be the the case that we
      // find a clique with rv and all of rv's node parents assigned (in
      // which case we assign rv right then), but in the mean time we
      // keep track of this clique's "score" using a set of
      // heuristics. 

      // TODO: change it so that only heuristics that are
      // used are computed.

      // Compute number of previous parents with their probability assigned.
      // The heuristic here is, have this node contribute probabiltiy
      // to this clique if many of its parents are already doing so, which
      // might produce a clique with good pruning behavior.
      set<RV*> res;
      set_intersection(curClique.assignedProbNodes.begin(),
		       curClique.assignedProbNodes.end(),
		       parSet.begin(),parSet.end(),
		       inserter(res,res.end()));
      int num_parents_with_probability = (int) res.size();
      // Previous Parents (earlier in the junction tree) with their
      // probabilities in Junction Tree.  We add this to the above.
      res.clear();
      set_intersection(curClique.cumulativeAssignedProbNodes.begin(),
		       curClique.cumulativeAssignedProbNodes.end(),
		       parSet.begin(),parSet.end(),
		       inserter(res,res.end()));
      num_parents_with_probability += (int) res.size();
      // negate so that lower is preferable.
      num_parents_with_probability *= -1;

      // Distance from root, among the higher priorities that are equal,
      // try to be as far away from the root as possible so as to prune
      // away as much zero as possible as early as possible. I.e., try
      // encourage it to have as much "influence" as possible in JT
      // parents.
      // negate so that lower (farther from the root) is preferable.
      const int distance_from_root = -(int)depth;

      // Number of children in current clique. If rv has lots of
      // children in this clique, it is hopeful that other parents of
      // those children might also be assigned to the same clique.
      int numChildren = 0;
      for (unsigned i=0;i<rv->allChildren.size();i++) {
	RV* child = rv->allChildren[i];
	if (curClique.nodes.find(child) != curClique.nodes.end())
	  numChildren++;
      }
      // negate so that lower (more children in current clique) is preferable.
      numChildren *= -1;

      // when the node is continuous, assign it to a clique that is the
      // smallest possible in terms of weight.
      // Again, negate so that larger weight is preferable.
      double weight = - curClique.weight();

      // TODO: add heuristic that assigns RV to smallest clique, in
      // terms of size. This might be useful for EM training,
      // particularly of continuous Gaussian variables, in that only
      // the parents of the Gaussian will be iterated over, rather
      // than parents of parents.

      // And so on. We can create as many heuristics as we want.
      // alternatively, perhaps take a weighted average of some of them.
      // TODO: add more heuristics here, and/or produce better
      // prioritized order above.

      // Other possible heuristics:
      //    1) a clique where it has the greatest number of node
      //    descendants (but this is only a heuristic since to gain
      //    benefit, we would need to have it be such that those node
      //    descendants have their node parents in clique as well.


      // We've now computed a bunch of heursitcs, push them
      // into an array in priority order to be scored later
      vector<double> score;
      if (rv->discrete()) {
	DiscRV* drv = (DiscRV*)rv;
	// 1st thing pushed has highest priority.
	if (drv->sparse() || drv->deterministic()) {
	  // if the RV is sparse/deterministic, use distance from root
	  // as the higest priority, since a sparse node will have
	  // lots of zeros, regardless of beam width. Goal is to prune
	  // away as early as possible.
	  score.push_back(distance_from_root);
	  score.push_back(num_parents_with_probability);
	} else {
	  // if the RV is not sparse, try to assign it to a clique
	  // with as many of its parents assigning probabiltiy, since
	  // a it will make for a clique that will prune more
	  // effectively
	  score.push_back(num_parents_with_probability);
	  score.push_back(distance_from_root);
	}
	// Last, assign based on prefering cliques where this variable
	// has lots of children (the children might then use the fact
	// that this RV gives probability to encourage it to live here
	// as well, again to improve pruning performance).
	score.push_back(numChildren);
      } else {
	// RV is continue. Pefer small weight.
	score.push_back(weight);
	// RV is continue. Then perfer lots of parents with probability.
	score.push_back(num_parents_with_probability);      
      }

      // done inserting heuristicss, now insert the score and the
      // current clique into the set.
      pair < vector<double>, unsigned> p(score,root);
      scoreSet.insert(p);
    }
  }

  // continue on down.
  for (unsigned childNo=0;
       childNo<curClique.children.size();childNo++) {
    const unsigned child = part.cliques[root].children[childNo];
    assignRVToClique(partName,part,child,depth+1,rv,numberOfTimesAssigned,parSet,allParentsObserved,scoreSet);
  }

}


/*-
 *-----------------------------------------------------------------------
 * JunctionTree::getCumulativeAssignedNodes()
 *    Decends down the directed JT graph and computes and
 *    assignes the variable cumulativeAsignedNodes.
 *
 * Preconditions:
 *     Must be called only by assignRVsToCliques()
 *
 * Postconditions:
 *     the clique's member variable cumulativeAsignedNodes has
 *     been re-computed.
 *
 * Side Effects:
 *     Might change the clique data member variables
 *
 * Results:
 *     None.
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::getCumulativeAssignedNodes(JT_Partition& part,
					 const unsigned root)
{
  MaxClique& curClique = part.cliques[root];

  for (unsigned childNo=0;
       childNo<curClique.children.size();childNo++) {

    const unsigned child = curClique.children[childNo];
    getCumulativeAssignedNodes(part,child);
    // Note: this will do a bunch of redundant work (i.e., inserting
    // elements into sets where the elements are already contained in
    // the sets), but it doesn't need to run that fast, so the
    // redundant work is not a big deal.
    union_1_2_to_3(part.cliques[child].cumulativeAssignedNodes,
		   part.cliques[child].assignedNodes,
		   curClique.cumulativeAssignedNodes,true);
    union_1_2_to_3(part.cliques[child].cumulativeAssignedProbNodes,
		   part.cliques[child].assignedProbNodes,
		   curClique.cumulativeAssignedProbNodes,true);

  }
}





/*-
 *-----------------------------------------------------------------------
 * JunctionTree::computeUnassignedCliqueNodes()
 *    Go through each clique and compute the unassigned nodes.
 *    This is only done when ceSeparatorDrivenInference == false.
 *
 * Preconditions:
 *     Must be called in appropriate order (setUpDataStructures).
 *
 * Postconditions:
 *     Cliques have this member variable filled in.
 *
 * Side Effects:
 *     Will change the clique data member variables
 *
 * Results:
 *     None.
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::computeUnassignedCliqueNodes(JT_Partition& part)
{
  if (MaxClique::ceSeparatorDrivenInference == false) {
    // these are only needed by clique driven inference.
    for (unsigned cliqueNum=0;cliqueNum<part.cliques.size();cliqueNum++) {
      part.cliques[cliqueNum].computeUnassignedCliqueNodes();
    }
  }
}
void
JunctionTree::computeUnassignedCliqueNodes()
{
  computeUnassignedCliqueNodes(P1);
  computeUnassignedCliqueNodes(Co);
  computeUnassignedCliqueNodes(E1);
}




/*-
 *-----------------------------------------------------------------------
 * JunctionTree::setUpMessagePassingOrders()
 *   sets up the message passing orders for the various
 *   basic partitions.
 * 
 *
 * Preconditions:
 *   Partitions must be instantiated, and interface cliques
 *   must have been computed (i.e., computePartitionInterfaces()
 *   must have been called)
 *
 * Postconditions:
 *   The member variables P_to_C_message_order (etc.) are created
 *   with the order of message passing.
 *
 * Side Effects:
 *    Changes member variables as described above.
 *
 * Results:
 *    None.
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::setUpMessagePassingOrders()
{

  setUpMessagePassingOrder(P1,
			   P_ri_to_C,
			   P1_message_order,
			   ~0x0,
			   P1_leaf_cliques);

  setUpMessagePassingOrder(Co,
			   C_ri_to_C,
			   Co_message_order,
			   C_li_to_P,
			   Co_leaf_cliques);

  setUpMessagePassingOrder(E1,
			   E_root_clique,
			   E1_message_order,
			   E_li_to_C,
			   E1_leaf_cliques);

}



/*-
 *-----------------------------------------------------------------------
 * JunctionTree::setUpMessagePassingOrders()
 *   sets up the message passing orders for a given
 *   specified partition, and places the result in the argument.
 *   Roots the undirected clique tree at 'root' thereby creating
 *   a directed tree of cliques. This does a recursive DFS tree traversal
 *   start at root to make the order, and places the visitation
 *   order in the array.
 * 
 * Preconditions:
 *   Same as setUpMessagePassingOrders() but for the
 *   particular arguments being passed in.
 *
 * Postconditions:
 *   order is set up.
 *.
 * Side Effects:
 *    none
 *
 * Results:
 *    None, but passed back in arguments.
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::setUpMessagePassingOrder(JT_Partition& part,
				       const unsigned root,
				       vector< pair<unsigned,unsigned> >& order,
				       const unsigned excludeFromLeafCliques,
				       vector < unsigned >& leaf_cliques)
{
  // first, handle case for potential empty P or E partition. 
  if (part.cliques.size() == 0)
    return;

  order.clear();
  // a tree of N nodes has N-1 edges.
  order.reserve(part.cliques.size() - 1);
  setUpMessagePassingOrderRecurse(part,root,order,excludeFromLeafCliques,leaf_cliques);
}
/*-
 *-----------------------------------------------------------------------
 * JunctionTree::setUpMessagePassingOrderRecurse()
 *
 *   Recursive support routine for setUpMessagePassingOrder() and
 *   should only be called from there. See that
 *   routine for documentation.
 * 
 *-----------------------------------------------------------------------
 */					
void
JunctionTree::setUpMessagePassingOrderRecurse(JT_Partition& part,
					      const unsigned root,
					      vector< pair<unsigned,unsigned> >& order,
					      const unsigned excludeFromLeafCliques,
					      vector < unsigned >& leaf_cliques)
{
  if (part.cliques[root].children.size() == 0) {
    // store leaf nodes here if children.size() == 0.
    if (root != excludeFromLeafCliques)
      leaf_cliques.push_back(root);
  } else {
    for (unsigned childNo=0;
	 childNo<part.cliques[root].children.size();
	 childNo++) {
      const unsigned child = part.cliques[root].children[childNo];
      setUpMessagePassingOrderRecurse(part,child,order,
				      excludeFromLeafCliques,leaf_cliques);
      pair<unsigned,unsigned> msg;
      msg.first = child;
      msg.second = root;
      order.push_back(msg);
    }
  }
}


/*-
 *-----------------------------------------------------------------------
 * JunctionTree::createSeperators()
 *   Creates the seperator objects based on the message passing order
 *   that has been created at this point. It creates seperators
 *   relative to the base partitions 0, 1, 2, and 3. Each
 *   seperator contains nodes consisting of the intersection
 *   of the maxCliques between which a message is sent.
 *   
 *   See comment in class JT_Partition for 'separators' member for more information.
 * 
 * Preconditions:
 *   setUpMessagePassingOrder() must have already been called.
 *
 * Postconditions:
 *   The seperators for partitions[0-3] have been created
 *   the right and left.
 *
 * Side Effects:
 *   Modifies the partitions seperators data structures.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::createSeparators()
{

  if (P1.cliques.size() > 0) { 
    createSeparators(P1,
		     P1_message_order);
    P1.cliques[P_ri_to_C].ceSendSeparator = ~0x0; // set to invalid value
    // P1 has no LI, so nothing more to do for P1 here.
    
    // insert VE seps last
    if (useVESeparators && (veSeparatorWhere & VESEP_WHERE_P)) {
      infoMsg(IM::Max,"Searching for VE seps in P1\n");
      createVESeparators(P1);
    } else
      P1.numVEseps = 0;
  }
  
  // Cs are never empty.
  assert ( Co.cliques.size() > 0 );

  createSeparators(Co,
		   Co_message_order);


  // insert VE seps before the last one
  if (useVESeparators && (veSeparatorWhere & VESEP_WHERE_C)) {
    infoMsg(IM::Max,"Searching for VE seps in Co\n");
    createVESeparators(Co);
  } else
    Co.numVEseps = 0;

  // Final one is guaranteed *always* to be the interface separator to
  // the left partition.
  if (gm_template.leftInterface) {
    // left interface case
    assert ( E1.cliques.size() > 0 );

    if (P1.cliques.size() > 0) {
      // Create separator of interface cliques, specifically the
      // separator between, P1's right-interface clique and Co's left
      // interface clqiue. 
      // 
      // Note that the interface separator always has to be the last
      // separator inserted.
      Co.separators.push_back(SeparatorClique(P1.cliques[P_ri_to_C],Co.cliques[C_li_to_P]));
      // update right partitions LI clique to include new separator
      Co.cliques[C_li_to_P].ceReceiveSeparators.push_back(Co.separators.size()-1);
    } else {
      // Then P1 is empty. We insert a dummy separator that has
      // been set up so that C'_i can go to C'_{i+1}. For i>0 this
      // works fine, for i=0 we will need to make sure during
      // inference not to use this separator since it will have
      // nothing in it.

      if (Co.liNodes.size() > 0) {
	// then there will at least be a Ci <-> C_{i+1}
	// intersection. I.e., we are guaranteed that the chunks are
	// not disconnected from each other.
	Co.separators.push_back(SeparatorClique(Co.cliques[C_li_to_P],Co.cliques[C_li_to_P]));
      } else {
	// then there will not be a Ci <-> C_{i+1} intersection.  need
	// to create an empty separator since the chunks are
	// disconnected from each other.
	set<RV*> empty;
	MaxClique mc(empty);
	Co.separators.push_back(SeparatorClique(mc,mc));
      }
      // update right partitions LI clique to include new separator
      Co.cliques[C_li_to_P].ceReceiveSeparators.push_back(Co.separators.size()-1);
    }
  } else {
    // right interface case
    assert ( P1.cliques.size() > 0 );

    // Note that the interface separator always has to be the last
    // separator inserted.
    Co.separators.push_back(SeparatorClique(P1.cliques[P_ri_to_C],Co.cliques[C_li_to_P]));
    // update right partitions LI clique to include new separator
    Co.cliques[C_li_to_P].ceReceiveSeparators.push_back(Co.separators.size()-1);

  }

  // don't update left partitions RI clique's send separator since handled explicitly
  Co.cliques[C_ri_to_C].ceSendSeparator = ~0x0; //set to invalid value

  if (E1.cliques.size() > 0) { 
    createSeparators(E1,
		     E1_message_order);

    // insert VE seps before final one.
    if (useVESeparators  && (veSeparatorWhere & VESEP_WHERE_E)) {
      infoMsg(IM::Max,"Searching for VE seps in E1\n");
      createVESeparators(E1);
    } else
      E1.numVEseps = 0;

    // Final one is guaranteed *always* to be the interface separator to
    // the left partition.
    // 
    // Create separator of interface cliques. Co is never empty. Note that the interface
    // separator always has to be the last separator inserted.
    E1.separators.push_back(SeparatorClique(Co.cliques[C_ri_to_E],E1.cliques[E_li_to_C]));
    // update right partitions LI clique to include new separator
    E1.cliques[E_li_to_C].ceReceiveSeparators.push_back(E1.separators.size()-1);
    // don't update left partitions RI clique's send separator since handeled explicitly
    E1.cliques[E_root_clique].ceSendSeparator = ~0x0; //set to invalid value
  }

}



/*-
 *-----------------------------------------------------------------------
 * JunctionTree::createSeparator()
 *   Creates the separator object based on the message passing order
 *   that has been created at this point for this partition.
 *
 * Preconditions:
 *   setUpMessagePassingOrder() must have already been called for
 *   this partition.
 *
 * Postconditions:
 *   The separators for this partition have been created
 *   the right and left.
 *
 * Side Effects:
 *   Modifies the partitions separators data structures.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::createSeparators(JT_Partition& part,
			       vector< pair<unsigned,unsigned> >& order)
{
  part.separators.clear();

  // first set up all the cliques so that they know who their separators are.
  for (unsigned p=0;p<order.size();p++) {
    const unsigned left = order[p].first;
    const unsigned right = order[p].second;    
    const unsigned sepNo = part.separators.size();
    part.separators.push_back(SeparatorClique(part.cliques[left],part.cliques[right]));
    part.cliques[left].ceSendSeparator = sepNo;
    part.cliques[right].ceReceiveSeparators.push_back(sepNo);
  }

}


/*-
 *-----------------------------------------------------------------------
 * JunctionTree::createVESeparator()
 *   Creates the virtual evidence (VE) separator objects based on any
 *   virtual evidence that might live in the current cliques.
 *
 * Preconditions:
 *   createSeparators() must have already been called.
 *
 * Postconditions:
 *   The VE separators for this partition have been created
 *   the right and left.
 *
 * Side Effects:
 *   Modifies the partitions separators data structures. Pushes
 *   each such separator back at the end of separator list.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::createVESeparators(JT_Partition& part)
{
  part.numVEseps = 0;
  if (!useVESeparators) {
    return;
  }
  for (unsigned c=0;c<part.cliques.size();c++) {
    // there might several VE separators.
    unsigned numVEsepsInClique = part.cliques[c].computeVESeparators();
    if (numVEsepsInClique > 0) {
      if (IM::messageGlb(IM::Max)) {
	printf("Found %d VE seps in clique number %d of current partition\n",numVEsepsInClique,c);
	for (unsigned i=0;i<numVEsepsInClique;i++) {
	  printf("VE sep %d,",i);
	  if (part.cliques[c].veSeparators[i].grandChild != NULL) 
	    printf("pcg:");
	  else
	    printf("pc:");
	  printRVSetAndCards(stdout,part.cliques[c].veSeparators[i].parents,false);
	  printf(",c=");
	  RV2DRV(part.cliques[c].veSeparators[i].child)->printNameFrameCard(stdout,false);
	  if (part.cliques[c].veSeparators[i].grandChild != NULL) {
	    printf(",gc=");
	    RV2DRV(part.cliques[c].veSeparators[i].grandChild)->printNameFrameCard(stdout,false);
	  }
	  printf("\n");
	}
      }
    }
    for (unsigned n=0;n<numVEsepsInClique;n++) {
      const unsigned sepNo = part.separators.size();
      part.separators.push_back(SeparatorClique(part.cliques[c].VESeparatorInfo(n)));
      part.cliques[c].ceReceiveSeparators.push_back(sepNo);
    }
    part.numVEseps += numVEsepsInClique;
  }
}



/*-
 *-----------------------------------------------------------------------
 * JunctionTree::computeSeparatorIterationOrders()
 *   computes the order in which we iterate over the separators
 *   when we do seperator driven clique potential function
 *   creation.
 *
 * Preconditions:
 *   separators both within and between partitions must have been
 *   created (meaning createSeparators() must be called).
 *
 * Postconditions:
 *   Separator orders have now been created.
 *
 *
 * Side Effects:
 *    Changes member variables within cliques of each partition.
 *
 * Results:
 *    None.
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::computeSeparatorIterationOrders()
{

  computeSeparatorIterationOrders(P1);
  computeSeparatorIterationOrders(Co);
  computeSeparatorIterationOrders(E1);
}



/*-
 *-----------------------------------------------------------------------
 * JunctionTree::computeSeparatorIterationOrders()/computeSeparatorIterationOrder()
 *
 *   1) computes the order in which we iterate over the separators
 *   when we do seperator driven clique potential function creation,
 *   for the given partition.
 *
 *   2) Also, set up the separator structures so they know the
 *   intersection of each of their nodes with the accumulated union of
 *   all nodes in separators earlier in the order.
 *
 *   3) Also, assigns unassignedIteratedNodes in this clique.
 *
 * Preconditions:
 *   separators must be created, should only be called from
 *   the general routine above w/o arguments.
 *
 * Postconditions:
 *   Separator orders have now been created within given partition.
 *
 *
 * Side Effects:
 *    Changes member variables within cliques of each partition.
 *
 * Results:
 *    None.
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::computeSeparatorIterationOrders(JT_Partition& part)
{
  for (unsigned cliqueNum=0;cliqueNum<part.cliques.size();cliqueNum++) {
    computeSeparatorIterationOrder(part.cliques[cliqueNum],part);
  }
}
void
JunctionTree::computeSeparatorIterationOrder(MaxClique& clique,
					     JT_Partition& part)
{
  // do an order based on the MST over the separators.
  const unsigned numSeparators = clique.ceReceiveSeparators.size();

  // Build up the partial and then ultimately the final union of of
  // all nodes in separators for incomming messages for this clique.
  clique.unionIncommingCESeps.clear();

  if (MaxClique::ceSeparatorDrivenInference) {
    if (numSeparators == 0) {
      // This must be a leaf-node clique relatve to root.
      // 'unionIncommingCESeps' is already empty so no need to do anything there.
    } else if (numSeparators == 1) {
      // shortcut to separator 0
      SeparatorClique& s0 = part.separators[clique.ceReceiveSeparators[0]];
      clique.unionIncommingCESeps = s0.nodes;
      s0.accumulatedIntersection.clear();
      s0.remainder = s0.nodes;
      assert ( s0.accumulatedIntersection.size() + s0.remainder.size() == s0.nodes.size() );
    } else if (numSeparators == 2) {


      // iterate through smaller weight separator first. If no
      // intersection, then this still can have a benefit since we
      // want to do the most work in the inner (rather than the outer)
      // loops. Of course, we don't know the real cost of the
      // separater state space until run time (which depends on the
      // evidence set, current parameters, and the pruning
      // thresholds), but we use the weight heuristic for now.
      if (part.separators[clique.ceReceiveSeparators[0]].weight() <=
	  part.separators[clique.ceReceiveSeparators[1]].weight()) {
	// do nothing
      } else {
	swap(clique.ceReceiveSeparators[0],clique.ceReceiveSeparators[1]);
      }

      // shortcuts to separator 0 and 1
      SeparatorClique& s0 = part.separators[clique.ceReceiveSeparators[0]];
      SeparatorClique& s1 = part.separators[clique.ceReceiveSeparators[1]];


      // intersection of the two separators
      set<RV*> sepIntersection;
      set_intersection(s0.nodes.begin(),s0.nodes.end(),
		       s1.nodes.begin(),s1.nodes.end(),
		       inserter(sepIntersection,sepIntersection.end()));

      // first one in order should be empty.
      s0.accumulatedIntersection.clear();
      // and remainder of first one should get the rest
      s0.remainder = s0.nodes;
      // 2nd one in should be intersection
      s1.accumulatedIntersection = sepIntersection;
      // and remainder gets the residual.
      set_difference(s1.nodes.begin(),
		     s1.nodes.end(),
		     sepIntersection.begin(),sepIntersection.end(),
		     inserter(s1.remainder,
			      s1.remainder.end()));
      // compute union of all separators.
      set_union(s0.nodes.begin(),s0.nodes.end(),
		s1.nodes.begin(),s1.nodes.end(),
		inserter(clique.unionIncommingCESeps,clique.unionIncommingCESeps.end()));

      assert ( s0.accumulatedIntersection.size() + s0.remainder.size() == s0.nodes.size() );
      assert ( s1.accumulatedIntersection.size() + s1.remainder.size() == s1.nodes.size() );

    } else {
      // there are 3 or more separators, determine proper order and then
      // compute running accumulated intersection relative to that order.

      // Goal: to produce an optimal separator iteration order.  In
      // general, the number of "live" surviving entries of all
      // separators when they are intersected the end of a sep
      // traversal will be the same no matter the order that is
      // chosen. What we want here, however, is something that kills
      // non-live entries ASAP, to avoid extra work.  This will depend
      // on properties of the separators (which we have here) and the
      // entries contained within the separators (which we don't have
      // until run time). Therefore, our order is a heuristic.
      // We do this by reordering the entries in clique.ceReceiveSeparators.

      // Other Ideas:
      // Probably a dynamic programing algorithm would work better
      // here. Use DP to find the order that minimizes the sum (in
      // reverse order) of the intersection of the last set with all
      // previous ones.

      //
      // What we first do: place VE seps last (since they are uneffected by pruning)
      // What we next do: sort in reverse order, by separator that has the
      // minimum intersection with others. At each step, when we have
      // found the sep with min intersection among current number of
      // seps, we place this at the *end* of the separation iteration
      // order. Once we get back down to two separators, we start with
      // the one of minimal weight (to do as much as possible in the
      // inner loops).  Also, when ties exist, do the thing that
      // causes the fewest number of branches and has the inner most
      // loop run interuppted for the greatest amount of time (i.e.,
      // min weight goes first).

      // printf("before sorting\n");
      // for (unsigned i = 0; i< clique.ceReceiveSeparators.size(); i++) {
      // printf("ceRecSep[%d] = %d\n",i,clique.ceReceiveSeparators[i]);
      // }

      // TODO: need to change this to be based only on hidden nodes,
      // since the intersection of observed nodes doesn't count
      // (observed nodes don't contribute to intersection pruning)
      // since they are never inserted into the clique.

#if 0
      for (unsigned i=0;i<clique.ceReceiveSeparators.size();i++) {
	if ((clique.ceReceiveSeparators[i] == (part.separators.size()-1))) {
	  swap(clique.ceReceiveSeparators[0],clique.ceReceiveSeparators[i]);
	  break;
	}
      }
#endif

      // First, place all non VE separators first in the order, since
      // these will be affected by the pruning parameters (while the
      // VE seps (if any) don't change with pruning.
      int swapPosition = numSeparators-1;
      for (int i=0;i<swapPosition;) {
	if (part.separators[clique.ceReceiveSeparators[i]].veSeparator) {
	  swap(clique.ceReceiveSeparators[i],clique.ceReceiveSeparators[swapPosition]);
	  swapPosition--;
	} else
	  i++;
      }
      if (!part.separators[clique.ceReceiveSeparators[swapPosition]].veSeparator)
	swapPosition++;
      const int firstVESeparator = swapPosition;
      const int lastRealSeparator = swapPosition-1;

      // printf("lastReal = %d, ceseps:",lastRealSeparator); 
      // for (unsigned i=0;i<numSeparators;i++) {
      // printf("%d ",clique.ceReceiveSeparators[i]);
      // }
      // printf("\n");

      // printf("Sorting Seps\n");

      // next, sort the real separators based on maximizing intersection.
      if (lastRealSeparator > 0) {
	// then we have at least 2 real separators.

	unsigned lastSeparator = lastRealSeparator;
	set<RV*> sep_intr_set;
	set<RV*> sep_union_set;
	set<RV*> empty;
	vector < pair<unsigned,unsigned> > sepIntersections; 
	while (lastSeparator > 1) {
	  // then at least three real separators to go.

	  // compute separator index with minimum intersection with all previous separators
	  // and swap its index with lastSeparator.
	  sepIntersections.clear();
	  sepIntersections.resize(lastSeparator+1);
	  for (unsigned i=0;i<=lastSeparator;i++) {
	    SeparatorClique& sep_i =
	      part.separators[clique.ceReceiveSeparators[i]];
	    sepIntersections[i].second = i;
	    sepIntersections[i].first = 0;
	    sep_union_set.clear();
	    for (unsigned j=0;j<=lastSeparator;j++) {
	      if (j == i)
		continue;
	      SeparatorClique& sep_j =
		part.separators[clique.ceReceiveSeparators[j]];
	      set_union(empty.begin(),empty.end(),
			sep_j.nodes.begin(),sep_j.nodes.end(),
			inserter(sep_union_set,sep_union_set.end()));
	    }
	    sep_intr_set.clear();
	    set_intersection(sep_i.nodes.begin(),sep_i.nodes.end(),
			     sep_union_set.begin(),sep_union_set.end(),
			     inserter(sep_intr_set,sep_intr_set.end()));
	    sepIntersections[i].first = sep_intr_set.size();
	  }

	  // sort in (default) ascending order
	  sort(sepIntersections.begin(),sepIntersections.end(),PairUnsigned1stElementCompare());
	
	  // among number of ties (the ones with minimal intersection
	  // with others), find the one with greatest weight to place
	  // last.
	  unsigned bestIndex = 0;
	  float bestWeight = part.separators[clique.ceReceiveSeparators[sepIntersections[bestIndex].second]].weight();
	  unsigned i = 1;
	  while (i < sepIntersections.size() && sepIntersections[bestIndex].first == sepIntersections[i].first) {
	    float curWeight = part.separators[clique.ceReceiveSeparators[sepIntersections[i].second]].weight();
	    if (curWeight > bestWeight) {
	      bestWeight = curWeight;
	      bestIndex = i;
	    }
	    i++;
	  }

	  // finally, swap the best into last place.
	  // printf("swaping sep %d with %d\n",bestIndex,lastSeparator);
	  swap(clique.ceReceiveSeparators[sepIntersections[bestIndex].second],clique.ceReceiveSeparators[lastSeparator]);
	  lastSeparator --;
	}

	// Lastly, the first two separators should be iterated
	// so that the min weight one goes first.
	if (part.separators[clique.ceReceiveSeparators[0]].weight() >
	    part.separators[clique.ceReceiveSeparators[1]].weight()) {
	  swap(clique.ceReceiveSeparators[0],clique.ceReceiveSeparators[1]);
	}
      }

      // next, sort  the VE separators to maximize overlap.
      const unsigned numVESeparators = numSeparators-firstVESeparator;
      if (numVESeparators > 1) {
	// so at least two VE separators.

	int lastSeparator = numSeparators-1;
	set<RV*> sep_intr_set;
	set<RV*> sep_union_set;
	set<RV*> empty;
	vector < pair<unsigned,unsigned> > sepIntersections; 
	while (lastSeparator > firstVESeparator + 1) {
	  // then at least three to go.

	  // compute separator index with minimum intersection with
	  // *all* previous separators (both normal and VE) and swap
	  // its index with lastSeparator.
	  sepIntersections.clear();
	  sepIntersections.resize(lastSeparator+1-firstVESeparator);
	  for (int i=firstVESeparator;i<=lastSeparator;i++) {
	    SeparatorClique& sep_i =
	      part.separators[clique.ceReceiveSeparators[i]];
	    sepIntersections[i-firstVESeparator].second = i;
	    sepIntersections[i-firstVESeparator].first = 0;
	    sep_union_set.clear();
	    for (int j=0;j<=lastSeparator;j++) {
	      if (j == i)
		continue;
	      SeparatorClique& sep_j =
		part.separators[clique.ceReceiveSeparators[j]];
	      set_union(empty.begin(),empty.end(),
			sep_j.nodes.begin(),sep_j.nodes.end(),
			inserter(sep_union_set,sep_union_set.end()));
	    
	    }
	    sep_intr_set.clear();
	    set_intersection(sep_i.nodes.begin(),sep_i.nodes.end(),
		      sep_union_set.begin(),sep_union_set.end(),
		      inserter(sep_intr_set,sep_intr_set.end()));
	    sepIntersections[i-firstVESeparator].first = sep_intr_set.size();
	  }
	  

// 	  printf("before sorting: VE seps have intersection among others:");
// 	  for (int i=firstVESeparator;i<=lastSeparator;i++) {
// 	    printf("%d has %d,",
// 		   sepIntersections[i-firstVESeparator].second,
// 		   sepIntersections[i-firstVESeparator].first);
// 	  }
// 	  printf("\n");


	  // sort in (default) ascending order
	  sort(sepIntersections.begin(),sepIntersections.end(),PairUnsigned1stElementCompare());

//  	  printf("after sorting: VE seps have intersection among others size = %d:",sepIntersections.size());
//  	  for (int i=firstVESeparator;i<=lastSeparator;i++) {
//  	    printf("Sep %d w inter %d.  ",
//  		   clique.ceReceiveSeparators[sepIntersections[i-firstVESeparator].second],
//  		   sepIntersections[i-firstVESeparator].first);
//  	  }
// 	  printf("\n");

// 	  vector < pair<unsigned,unsigned> >::iterator it;
// 	  for (it = sepIntersections.begin(); 
// 		 ((it) != sepIntersections.end()); it++) {
// 	    printf("cur=(%d,%d), next=(%d,%d), cur<next=%d\n",
// 		   (*it).first,(*it).second,
// 		   (*(it+1)).first,(*(it+1)).second,
// 		   (*it)<(*(it+1)));
// 	  }

	
	  // among number of ties (the ones with minimal intersection
	  // with others), find the one with greatest weight to place
	  // last.
	  unsigned bestIndex = 0;
	  float bestWeight = part.separators[clique.ceReceiveSeparators[sepIntersections[bestIndex].second]].weight();
	  unsigned i = 1;
	  while (i < sepIntersections.size() && sepIntersections[bestIndex].first == sepIntersections[i].first) {
	    float curWeight = part.separators[clique.ceReceiveSeparators[sepIntersections[i].second]].weight();
	    if (curWeight > bestWeight) {
	      bestWeight = curWeight;
	      bestIndex = i;
	    }
	    i++;
	  }

	  // finally, swap the best into last place.
	  // printf("swaping sep %d with %d\n",bestIndex,lastSeparator);
	  swap(clique.ceReceiveSeparators[sepIntersections[bestIndex].second],clique.ceReceiveSeparators[lastSeparator]);
	  lastSeparator --;
	}

	// Lastly, the first two separators should be iterated
	// so that the min weight one goes first.
	if (part.separators[clique.ceReceiveSeparators[firstVESeparator]].weight() >
	    part.separators[clique.ceReceiveSeparators[firstVESeparator+1]].weight()) {
	  swap(clique.ceReceiveSeparators[firstVESeparator],clique.ceReceiveSeparators[firstVESeparator+1]);
	}
      }


      ///////////////////////////////////////////////////////
      // DONE SORTING clique.ceReceiveSeparators
      // printf("after sorting\n");
      // for (unsigned i = 0; i< clique.ceReceiveSeparators.size(); i++) {
      // printf("ceRecSep[%d] = %d\n",i,clique.ceReceiveSeparators[i]);
      // }
      // Compute the cummulative intersection of the sepsets using the
      // currently assigned sepset order.
      {

	// initialize union of all previous separators

	clique.unionIncommingCESeps = 
	  part.separators[clique.ceReceiveSeparators[0]].nodes;
	part.separators[clique.ceReceiveSeparators[0]].accumulatedIntersection.clear();
	part.separators[clique.ceReceiveSeparators[0]].remainder
	  = part.separators[clique.ceReceiveSeparators[0]].nodes;

	for (unsigned sep=1;sep<numSeparators;sep++) {
      
	  // reference variables for easy access
	  set<RV*>& sepNodes
	    = part.separators[clique.ceReceiveSeparators[sep]].nodes;
	  set<RV*>& sepAccumInter
	    = part.separators[clique.ceReceiveSeparators[sep]].accumulatedIntersection;


	  // create the intersection of 1) the union of all previous nodes in
	  // the sep order, and 2) the current sep nodes.
	  sepAccumInter.clear();
	  set_intersection(clique.unionIncommingCESeps.begin(),clique.unionIncommingCESeps.end(),
			   sepNodes.begin(),sepNodes.end(),
			   inserter(sepAccumInter,sepAccumInter.end()));

	  // compute the separator remainder while we're at it.
	  // specifically: remainder = sepNodes - sepAccumInter.
	  part.separators[clique.ceReceiveSeparators[sep]].remainder.clear();
	  set_difference(sepNodes.begin(),sepNodes.end(),
			 sepAccumInter.begin(),sepAccumInter.end(),
			 inserter(part.separators[clique.ceReceiveSeparators[sep]].remainder,
				  part.separators[clique.ceReceiveSeparators[sep]].remainder.end()));


	  // update the accumulated (union) of all previous sep nodes.
	  set<RV*> res;
	  set_union(sepNodes.begin(),sepNodes.end(),
		    clique.unionIncommingCESeps.begin(),clique.unionIncommingCESeps.end(),
		    inserter(res,res.end()));
	  clique.unionIncommingCESeps = res;	
	}
      }

    }
  } else {
    // we are doing clique driven iteration 

    // In the case of CLIQUE-DRIVEN inference: while we could use the
    // same code as above (producing both accumulated intersection and
    // remainder sets), we optimize for this case by making sure that
    // all accumualted intersection sets should be of zero size, so
    // that we only do one hash lookup for the remainder (which in
    // this case will be everything). Note also that the order of the
    // separators no longer matters.

    const set < RV* > empty;
    for (unsigned i=0;i<numSeparators;i++) {
      SeparatorClique& sep = part.separators[clique.ceReceiveSeparators[i]];
      sep.accumulatedIntersection.clear();
      sep.remainder = sep.nodes;
      // compute union of all separators.
      set_union(empty.begin(),empty.end(),
		sep.nodes.begin(),sep.nodes.end(),
		inserter(clique.unionIncommingCESeps,clique.unionIncommingCESeps.end()));

    }
  }

  // lastly, assign unassignedIteratedNodes in this clique
  {
    // unassigned nodes, unassigned iterated nodes, 
    // compute: unassignedIteratedNodes  = nodes - (assignedNodes U unionIncommingCESeps)
    // first: compute res = nodes - assignedNodes
    set<RV*> res;
    set_difference(clique.nodes.begin(),clique.nodes.end(),
		   clique.assignedNodes.begin(),clique.assignedNodes.end(),
		   inserter(res,res.end()));
    // next: compute unassignedIteratedNodes = res - unionIncommingCESeps
    // note at this point unionIncommingCESeps contains the union of all nodes in all separators
    set_difference(res.begin(),res.end(),
		   clique.unionIncommingCESeps.begin(),clique.unionIncommingCESeps.end(),
		   inserter(clique.unassignedIteratedNodes,
			    clique.unassignedIteratedNodes.end()));
  }
  // @@@ perhaps subtract out unassignedInPartition since they are different??

}




/*-
 *-----------------------------------------------------------------------
 * JunctionTree::getCumulativeUnassignedIteratedNodes()
 *
 *   Computes for each clique the union of the set of nodes which have
 *   been iterated unassigned in previous cliques in the JT. This is
 *   used to determine which, in each clique, nodes should be iterated
 *   over and which are already assigned by a separator driven
 *   iteration.
 *
 * Preconditions:
 *   computeSeparatorIterationOrder() must have been called. Partition
 *   in argument 'part' must not be empty.
 *
 * Postconditions:
 *   cliques know which nodes to iterate over or not.
 *
 * Side Effects:
 *    Changes member variables within cliques of this partition.
 *
 * Results:
 *    None.
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::getCumulativeUnassignedIteratedNodes(JT_Partition& part,
						   const unsigned root)
{
  MaxClique& curClique = part.cliques[root];
  set<RV*>& res = curClique.cumulativeUnassignedIteratedNodes;
  for (unsigned childNo=0;
       childNo<curClique.children.size();childNo++) {

    const unsigned child = curClique.children[childNo];
    getCumulativeUnassignedIteratedNodes(part,child);

    set_union(part.cliques[child].cumulativeUnassignedIteratedNodes.begin(),
	      part.cliques[child].cumulativeUnassignedIteratedNodes.end(),
	      part.cliques[child].unassignedIteratedNodes.begin(),
	      part.cliques[child].unassignedIteratedNodes.end(),
	      inserter(res,res.end()));
  }
  curClique.computeAssignedNodesDispositions();
}
void
JunctionTree::getCumulativeUnassignedIteratedNodes()
{
  set<RV*> res;

  if (P1.cliques.size() > 0) {
    // TODO: this should be a member function in P1.
    getCumulativeUnassignedIteratedNodes(P1,P_ri_to_C);
    res.clear();
    set_union(P1.cliques[P_ri_to_C].unassignedIteratedNodes.begin(),
	      P1.cliques[P_ri_to_C].unassignedIteratedNodes.end(),
	      P1.cliques[P_ri_to_C].cumulativeUnassignedIteratedNodes.begin(),
	      P1.cliques[P_ri_to_C].cumulativeUnassignedIteratedNodes.end(),	    
	      inserter(res,res.end()));
  }

  // Co is never empty.
  Co.cliques[C_li_to_P].cumulativeUnassignedIteratedNodes = res;
  getCumulativeUnassignedIteratedNodes(Co,C_ri_to_C);

  if (E1.cliques.size() > 0) {
    res.clear();
    set_union(Co.cliques[C_ri_to_C].unassignedIteratedNodes.begin(),
	      Co.cliques[C_ri_to_C].unassignedIteratedNodes.begin(),
	      Co.cliques[C_ri_to_C].cumulativeUnassignedIteratedNodes.begin(),
	      Co.cliques[C_ri_to_C].cumulativeUnassignedIteratedNodes.begin(),	    
	      inserter(res,res.end()));
    E1.cliques[E_li_to_C].cumulativeUnassignedIteratedNodes = res;
    getCumulativeUnassignedIteratedNodes(E1,E_root_clique);
  }

}



/*-
 *-----------------------------------------------------------------------
 * JunctionTree::junctionTreeWeight()
 *
 * Compute the 'junction tree weight' (roughly, the log10(cost of
 * doing inference)) for the set of cliques given in cliques. Note,
 * cliques *must* be a valid set of maxcliques of a junction tree --
 * if they are not, unexpected results are returned.
 *
 * Note also, this returns an upper bound on the cost of a JT in this
 * partition, rather than the actual cost. It returns an upper bound
 * since it assumes that the separator driven iteration will be at
 * full cardinality (it can't tell the difference between previous
 * assigned nodes which have been heavily cut back and previous
 * unassigned nodes which are really at full cardinality). Therefore,
 * this routine returns the conservative upper bound. The goal
 * is to minimize this upper bound.
 *
 *
 * Preconditions:
 *   The partition must have been fully instantiated. I.e.,
 *   we must have that assignRVsToCliques have been called.
 *
 * Postconditions:
 *   The weight, as inference would do it, of this clique is 
 *   computed.
 *
 * Side Effects:
 *   none.
 *
 * Results:
 *   the weight
 *
 *-----------------------------------------------------------------------
 */
double
JunctionTree::junctionTreeWeight(JT_Partition& part,
				 const unsigned rootClique,
				 set<RV*>* lp_nodes,
				 set<RV*>* rp_nodes)
{
  // check for case of empty E or P partition.
  if (part.cliques.size() == 0)
    return 0;

  MaxClique& curClique = part.cliques[rootClique];

  set <RV*> empty;
  // @@ add in unassigned in partition information to next call
  double weight = curClique.weightInJunctionTree(part.unassignedInPartition,
						 jtWeightUpperBound,
						 jtWeightMoreConservative,
						 true,
						 lp_nodes,rp_nodes);
  // double weight = curClique.weight();
  for (unsigned childNo=0;
       childNo<curClique.children.size();childNo++) {
    double child_weight = junctionTreeWeight(part,
					     curClique.children[childNo],lp_nodes,rp_nodes);
    // i.e., log addition for weight = weight + child_weight
    weight = log10add(weight,child_weight);
  }
  return weight;
}



/*-
 *-----------------------------------------------------------------------
 * JunctionTree::junctionTreeWeight()
 *
 * Given a set of maxcliques for a partition, and an interface for
 * this (can be left right, or any set including empty, the only
 * condition is that it must be covered by at least one of the
 * cliques), compute the junction tree for this set and return the
 * estimated JT cost. This is a static routine so can be called from
 * anywhere.
 *
 * Preconditions:
 *   one of the cliques must cover interface nodes. I.e.,
 *   there must be an i such that cliques[i].nodes >= interfaceNodes.
 *
 * Postconditions:
 *   An estimate of the weight, as in inference, of this clique set is 
 *   computed and returned.
 *
 * Side Effects:
 *   none.
 *
 * Results:
 *   the weight
 *
 *-----------------------------------------------------------------------
 */
double
JunctionTree::junctionTreeWeight(vector<MaxClique>& cliques,
				 const set<RV*>& interfaceNodes,
				 set<RV*>* lp_nodes,
				 set<RV*>* rp_nodes)

{

  // might be a partition that is empty.
  if (cliques.size() == 0)
    return 0.0;

  Partition part;
  const set <RV*> emptySet;
  part.cliques = cliques;


  // set up the nodes into part
  for (unsigned i=0;i<cliques.size();i++) {
    // while we're at it, clear up the JT structures we're
    // about to create.
    cliques[i].clearJTStructures();
    set_union(emptySet.begin(),emptySet.end(),
	      cliques[i].nodes.begin(),cliques[i].nodes.end(),
	      inserter(part.nodes,part.nodes.end()));
  }
  createPartitionJunctionTree(part);

  // a JT version of this partition.
  JT_Partition jt_part(part,emptySet,interfaceNodes);

  bool tmp;
  unsigned root;
  if (interfaceNodes.size() > 0)
    jt_part.findRInterfaceClique(root,tmp,interfaceCliquePriorityStr);
  else {
    // 'E_root_clique case'
    // Presumably, this is an E partition, so the root should be done
    // same as E_root_clique computed above.
    // "update E_root_clique"
    root = 0; 
  }
    
  createDirectedGraphOfCliques(jt_part,root);
  assignRVsToCliques("candidate partition",jt_part,root,"","");

  vector< pair<unsigned,unsigned> > message_order;
  vector< unsigned > leaf_cliques;
  setUpMessagePassingOrder(jt_part,
			   root,
			   message_order,
			   ~0x0,
			   leaf_cliques);

  // Note that since we're calling createSeparators(jt_part,msg_order)
  // here directly, this means that the left interface clique will not
  // have its separator to the left partition filled in.
  createSeparators(jt_part,message_order);

  computeSeparatorIterationOrders(jt_part);
  getCumulativeUnassignedIteratedNodes(jt_part,root);
  
  // return jt_part.cliques[root].cumulativeUnassignedIteratedNodes.size();
  double weight = junctionTreeWeight(jt_part,root,lp_nodes,rp_nodes);

  if (jtWeightPenalizeUnassignedIterated > 0.0) {
    unsigned badness_count=0;
    set <RV*>::iterator it;
    for (it = jt_part.cliques[root].cumulativeUnassignedIteratedNodes.begin();
	 it != jt_part.cliques[root].cumulativeUnassignedIteratedNodes.end();
	 it++) 
      {
	RV* rv = (*it);
	
	assert ( rv-> discrete() );
	DiscRV* drv = (DiscRV*)rv;
	if (!drv->sparse())
	  continue;

	// if it has not at all been assigned in this partition, we
	// assume it has been assigned in previous partition, and dont'
	// count it's cost.
	if (jt_part.unassignedInPartition.find(rv) ==
	    jt_part.unassignedInPartition.end())
	  continue;
	
	badness_count ++;
      }
    weight = badness_count*jtWeightPenalizeUnassignedIterated + weight;
  }
  return weight;

  // return badness_count*1000 + weight;
  // return badness_count;

}


/*-
 *-----------------------------------------------------------------------
 * JunctionTree::printAllJTInfo()
 *
 *   Prints all information to a file that has been computed
 *   for this junction tree.
 *
 * Preconditions:
 *   All JT creation functions, ending in
 *   computeSeparatorIterationOrders() must have been
 *   called. If not, this function will print garbage (or will crash).
 *
 * Postconditions:
 *   All stuff printed to given file.
 *
 * Side Effects:
 *    none
 *
 * Results:
 *    None.
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::printAllJTInfo(char *fileName) 
{
  FILE* f;
  if ((fileName == NULL)
      ||
      ((f = ::fopen(fileName,"w")) == NULL))
    return;


  // print partition (clique,separator) information

  fprintf(f,"===============================\n");
  fprintf(f,"   P1 partition information: JT_weight = %f\n",
	  junctionTreeWeight(P1,P_ri_to_C,
			     NULL,&Co.nodes));
  printAllJTInfo(f,P1,P_ri_to_C,NULL,&Co.nodes);
  fprintf(f,"\n\n");

  fprintf(f,"===============================\n");
  fprintf(f,"   Co partition information: JT_weight = %f\n",
	  junctionTreeWeight(Co,C_ri_to_E,&P1.nodes,&E1.nodes));
  printAllJTInfo(f,Co,C_ri_to_C,&P1.nodes,&E1.nodes);
  fprintf(f,"\n\n");

  fprintf(f,"===============================\n");
  fprintf(f,"   E1 partition information: JT_weight = %f\n",
	  junctionTreeWeight(E1,E_root_clique,
			     &Co.nodes,NULL));
  printAllJTInfo(f,E1,E_root_clique,&Co.nodes,NULL);
  fprintf(f,"\n\n");

  // print message order information
  fprintf(f,"===============================\n\n");    

  fprintf(f,"===============================\n");  
  fprintf(f,"   P1 message order\n");
  printMessageOrder(f,P1_message_order);
  fprintf(f,"\n\n");

  fprintf(f,"===============================\n");  
  fprintf(f,"   Co message order\n");
  printMessageOrder(f,Co_message_order);
  fprintf(f,"\n\n");

  fprintf(f,"===============================\n");  
  fprintf(f,"   E1 message order\n");
  printMessageOrder(f,E1_message_order);
  fprintf(f,"\n\n");


  fclose(f);
}


/*-
 *-----------------------------------------------------------------------
 * JunctionTree::printAllJTInfo()
 *
 *   Prints all information to a file that has been computed
 *   for this junction tree for a particular given partition and given root clique.
 *
 * Preconditions:
 *   All JT creation functions, ending in
 *   computeSeparatorIterationOrders() must have been
 *   called. If not, this function will print garbage (or will crash).
 *
 * Postconditions:
 *   All stuff printed to given file.
 *
 * Side Effects:
 *    none
 *
 * Results:
 *    None.
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::printAllJTInfo(FILE* f,
			     JT_Partition& part,
			     const unsigned root,
			     set <RV*>* lp_nodes,
			     set <RV*>* rp_nodes)
			     
{
  // print cliques information
  fprintf(f,"=== Clique Information ===\n");
  fprintf(f,"Number of cliques = %d\n",part.cliques.size());
  if (part.cliques.size() > 0)
    printAllJTInfoCliques(f,part,root,0,lp_nodes,rp_nodes);

  // print separator information
  fprintf(f,"\n=== Separator Information ===\n");
  fprintf(f,"Number of separators = %d\n",part.separators.size());
  for (unsigned sepNo=0;sepNo<part.separators.size();sepNo++) {
    fprintf(f,"== Separator number: %d\n",sepNo);
    part.separators[sepNo].printAllJTInfo(f);
  }

}

/*-
 *-----------------------------------------------------------------------
 * JunctionTree::printAllJTInfo()
 *
 *   Recursive version that prints all junction tree information to a
 *   file that has been computed for this junction tree for a
 *   particular given partition and given root clique. It is recursive
 *   so that it can print the JT info to a file in tree-traversal order.
 *
 * Preconditions:
 *   All JT creation functions, ending in
 *   computeSeparatorIterationOrders() must have been
 *   called. If not, this function will print garbage (or will crash).
 *
 * Postconditions:
 *   All stuff printed to given file.
 *
 * Side Effects:
 *    none
 *
 * Results:
 *    None.
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::printAllJTInfoCliques(FILE* f,
				    JT_Partition& part,
				    const unsigned root,
				    const unsigned treeLevel,
				    set <RV*>* lp_nodes,
				    set <RV*>* rp_nodes)
{
  // print cliques information
  for (unsigned i=0;i<treeLevel;i++) fprintf(f,"  ");
  fprintf(f,"== Clique number: %d",root);
  if (treeLevel == 0)
    fprintf(f,", root/right-interface clique");
  if (part.cliques[root].ceReceiveSeparators.size() == 
      part.cliques[root].numVESeparators() + part.cliques[root].children.size() + 1) {
    if (part.cliques[root].ceReceiveSeparators.size() > 1 + part.cliques[root].numVESeparators()) 
      fprintf(f,", left-interface clique");
    else
      fprintf(f,", leaf/left-interface clique");
  } else {
    assert ( part.cliques[root].ceReceiveSeparators.size() == 
	     part.cliques[root].children.size() + part.cliques[root].numVESeparators());
    if (part.cliques[root].ceReceiveSeparators.size() == part.cliques[root].numVESeparators()) 
      fprintf(f,", leaf");
  }
  fprintf(f,"\n");
  part.cliques[root].printAllJTInfo(f,treeLevel,part.unassignedInPartition,
				    jtWeightUpperBound,
				    jtWeightMoreConservative,
				    true,
				    lp_nodes,rp_nodes);
  for (unsigned childNo=0;
       childNo<part.cliques[root].children.size();childNo++) {
    unsigned child = part.cliques[root].children[childNo];
    printAllJTInfoCliques(f,part,child,treeLevel+1,lp_nodes,rp_nodes);
  }
}



/*-
 *-----------------------------------------------------------------------
 * JunctionTree::printAllJTInfo()
 *
 *   Print the given message_order to the given file.
 *
 * Preconditions:
 *   All JT creation functions, ending in
 *   computeSeparatorIterationOrders() must have been
 *   called. If not, this function will print garbage (or will crash).
 *
 * Postconditions:
 *   All stuff printed to given file.
 *
 * Side Effects:
 *    none
 *
 * Results:
 *    None.
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::printMessageOrder(FILE *f,
				vector< pair<unsigned,unsigned> >& message_order)
{
  fprintf(f,"Number of messages: %d\n",message_order.size());
  for (unsigned m=0;m<message_order.size();m++) {
    const unsigned from = message_order[m].first;
    const unsigned to = message_order[m].second;
    fprintf(f,"  %d: %d --> %d\n",m,from,to);
  }
}




/*-
 *-----------------------------------------------------------------------
 * JunctionTree::prepareForUnrolling()
 *
 *  Sets up a number of data structures that must be ready before
 *  we can do general unrolling. This includes:
 *  
 *   1) computes the number of hidden variables in each clique/separator
 *      (i.e., the number of discrete hidden and non-continous variables)
 *
 * Preconditions:
 *   computeSeparatorIterationOrders() must have been called.
 *
 * Postconditions:
 *   everything is now ready for general unrolling.
 *
 *
 * Side Effects:
 *    Changes member variables within cliques/separators within each partition.
 *
 * Results:
 *    None.
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::prepareForUnrolling()
{

  const unsigned totalNumVESeps =
    P1.numVEseps + Co.numVEseps + E1.numVEseps;
  if (useVESeparators && totalNumVESeps > 0) {
    // possibly set VE sep file

    SeparatorClique::veSeparatorFile = NULL;

    if (SeparatorClique::recomputeVESeparatorTables || 
	((SeparatorClique::veSeparatorFileName != NULL) &&
	 (::fsize(SeparatorClique::veSeparatorFileName) == 0))) {
      // open file for writing since we're re-generating the information.

      if (SeparatorClique::veSeparatorFileName != NULL) {
	SeparatorClique::veSeparatorFile =
	  ::fopen(SeparatorClique::veSeparatorFileName,"w");
	if (SeparatorClique::veSeparatorFile == NULL) {
	  error("ERROR: cannot open VE separator file (%s) for writing.",SeparatorClique::veSeparatorFileName);
	}
      }
      SeparatorClique::generatingVESeparatorTables = true;
      infoMsg(IM::Default,"Computing information for %d total VE separators\n",totalNumVESeps);
    } else if (SeparatorClique::veSeparatorFileName != NULL) {
      // assume that the current ve sep file is valid.
      SeparatorClique::veSeparatorFile =
	::fopen(SeparatorClique::veSeparatorFileName,"r");
      if (SeparatorClique::veSeparatorFile == NULL) {
	error("ERROR: cannot open VE separator file (%s) for reading.",SeparatorClique::veSeparatorFileName);
      }
      SeparatorClique::generatingVESeparatorTables = false;
      infoMsg(IM::Default,"Reading information for %d total VE separators\n",totalNumVESeps);
    }


  }

  prepareForUnrolling(P1);
  prepareForUnrolling(Co);
  prepareForUnrolling(E1);

  if (useVESeparators && totalNumVESeps > 0) {
    // close file
    if (SeparatorClique::veSeparatorFile != NULL)
      fclose(SeparatorClique::veSeparatorFile);
    if (SeparatorClique::generatingVESeparatorTables)
      infoMsg(IM::Default,"Done computing information for %d total VE separators\n",totalNumVESeps);
  }
}

void
JunctionTree::prepareForUnrolling(JT_Partition& part)
{
  for (unsigned i=0;i<part.cliques.size();i++) {
    part.cliques[i].prepareForUnrolling();
  }
  for (unsigned i=0;i<part.separators.size();i++) {
    part.separators[i].prepareForUnrolling();
  }
}


/*-
 *-----------------------------------------------------------------------
 * JunctionTree::unroll()
 *   sets up the internal data structures and creates a junction tree corresponding
 *   to the graph with C' unrolled k times, k >= 0. (i.e., we have
 *   a graph that has (k+1) copies of C'). 
 *
 *   This also set up data structures for inference. Some of what this routine
 *   does is copy data out of STL structures into faster customized structures
 *   (such as hash tables, arrays, etc.).
 *
 *   Sets things up so that inference can take place.
 *   Unrolling is in terms of new C' partitions (so that if S=1, k corresponds to frames)
 *   but in general the network will be T(P') + k*T(C') + T(E) frames long, and in general
 *   T(C') >= S. This also depends on if the boundary algorithm was run or not. 
 *
 * Preconditions:
 *   The partitions must be validly instantiated with cliques, and
 *   the routine assignRVsToCliques() must have been called.
 *
 * Postconditions:
 *   The partitions are such that they now know who they are connected to on
 *   the right and left.
 *   Also, unrolled data structures are set up and ready for inference for this length.
 *
 * Side Effects:
 *   Modifies all neighbors variables within the cliques within the
 *   partition.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
unsigned
JunctionTree::unroll(const unsigned int numFrames)
{

  // first create the unrolled set of random variables corresponding
  // to this JT.

  unsigned basicTemplateUnrollAmount;
  int modifiedTemplateUnrollAmount;
  unsigned numUsableFrames;
  unsigned frameStart;

  if (!gm_template.computeUnrollParameters(numFrames,
					   basicTemplateUnrollAmount,
					   modifiedTemplateUnrollAmount,
					   numUsableFrames,
					   frameStart))
    error("Segment of %d frames too short with current GMTK template of length [P=%d,C=%d,E=%d] %d frames, and M=%d,S=%d boundary parameters. Use longer utterances, different template, or decrease M,S if >1.\n",
	  numFrames,
	  fp.numFramesInP(),fp.numFramesInC(),fp.numFramesInE(),
	  fp.numFrames(),
	  gm_template.M,gm_template.S);

  infoMsg(IM::Info,"Number of Current Frames = %d, Number of Currently Usable Frames = %d\n",numFrames,numUsableFrames);
  infoMsg(IM::Tiny,"Number Of Frames = %d, Unrolling Basic Template %d times, Modified Template %d times\n",
	  numFrames,
	  basicTemplateUnrollAmount,
	  modifiedTemplateUnrollAmount);

  // unrolled random variables
  // vector <RV*> unrolled_rvs;
  // mapping from 'name+frame' to integer index into unrolled_rvs.
  // map < RVInfo::rvParent, unsigned > ppf;

  fp.unroll(basicTemplateUnrollAmount,cur_unrolled_rvs,cur_ppf);
  setObservedRVs(cur_unrolled_rvs);

  // TODO: clear out the old and pre-allocate for new size.
  jtIPartitions.clear();
  // re-allocate.
  jtIPartitions.resize(modifiedTemplateUnrollAmount+3);
  // this clears the shared caches. 
  clearCliqueSepValueCache();

  unsigned partNo = 0;
  const int numCoPartitions = modifiedTemplateUnrollAmount+1;

  new (&jtIPartitions[partNo++]) JT_InferencePartition(P1,cur_unrolled_rvs,cur_ppf,0*gm_template.S);
  for (int p=0;p<numCoPartitions;p++) {
    new (&jtIPartitions[partNo]) JT_InferencePartition(Co,cur_unrolled_rvs,cur_ppf,p*gm_template.S);
    if (p == 0 && P1.cliques.size() == 0) {
      // Alternatively to the below, we might set this inference
      // partition to not use the LI separator, since it is will
      // always be empty at the start. At the moment, we handle this
      // in the inference CE/DE code below.
    }
    partNo++;
  }

  new (&jtIPartitions[partNo++]) 
    JT_InferencePartition(E1,cur_unrolled_rvs,cur_ppf,
			  (modifiedTemplateUnrollAmount)*gm_template.S);

  assert (partNo == jtIPartitions.size());

  return numUsableFrames;
}



/*-
 *-----------------------------------------------------------------------
 * JunctionTree::setObservedRVs()
 *   sets the observed RVs to their values, either taking values from
 *   the global observation matrix, or taking the values from the files.
 *
 * Preconditions:
 *   The given RVs must come from the result of unroll, and the observation matrix
 *   *must* be set up and ready to be used.
 *
 * Postconditions:
 *   All discrete observed random variables have a value that is their appropriate observed
 *   value.
 *
 * Side Effects:
 *   Modifies values of discrete observed random variables.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::setObservedRVs(vector <RV*>& rvs)
{
  // Set all discrete observed variables to their values here
  vector <RV*>::iterator it;
  for (it = rvs.begin(); it!= rvs.end(); it++) {
    RV *rv = (*it);
    if (rv->discrete() && !rv->hidden()) {
      DiscRV* drv = (DiscRV*)rv;
      drv->setToObservedValue();
    }
  }
}




/*-
 *-----------------------------------------------------------------------
 * JunctionTree::ceGatherIntoRoot
 *   
 *   Collect Evidence Gather Into Root: This routine does a collect
 *   evidence pass for this partition, and gathers all messages into
 *   the root within the current partition. It does so using the
 *   message order given in the argument 'message_order', and gathers
 *   into the provided root clique.  
 *
 * See Also:
 *   Dual routine: JunctionTree::deScatterOutofRoot()
 *
 *
 * Preconditions:
 *   It is assumed that either:
 *     1) this is the left-most partition
 *  or 2) that the left interface clique within this partition has had a message
 *        sent to it from the left neighbor partition.
 *
 * Postconditions:
 *     All cliques in the partition have all messages but one sent to it.
 *
 * Side Effects:
 *     all partitions will have been instantiated to the extent that the messages (with
 *     the current pruning ratios)  have been created.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::ceGatherIntoRoot(// partition
			       JT_InferencePartition& part,
			       // index of root clique in the partition
			       const unsigned root,
			       // message order of the JT in this partition
			       vector< pair<unsigned,unsigned> >& message_order,
			       // the name of the partition (for debugging/status msgs)
			       const char*const part_type_name,
			       // number of the partition in unrolled graph 
			       // (for debugging/status msgs)
			       const unsigned part_num)
{
  // first check that this is not an empty partition.
  if (part.maxCliques.size() == 0)
    return;

  // Now, do partition messages.
  for (unsigned msgNo=0;msgNo < message_order.size(); msgNo ++) {
    const unsigned from = message_order[msgNo].first;
    const unsigned to = message_order[msgNo].second;
    // TODO: add another info message here
    part.maxCliques[from].
      ceGatherFromIncommingSeparators(part);
    infoMsg(IM::Mod,
	    "CE: message %s,part[%d]: clique %d --> clique %d\n",
	    part_type_name,part_num,from,to);
    part.maxCliques[from].
      ceSendToOutgoingSeparator(part);
  }
  // collect to partition's root clique
  // TODO: add another info message here
  part.maxCliques[root].
    ceGatherFromIncommingSeparators(part);
}



/*-
 *-----------------------------------------------------------------------
 * JunctionTree::ceSendToNextPartition
 *   
 *   Collect Evidence Send To Next Partition: This routine sends a
 *   message from the right interface clique of a left (or previous)
 *   partition to the left interface clique of a right (or next)
 *   partition in the partition series. It is assumed that the right
 *   interface clique has had all its incomming messages sent to it.
 *
 * See Also:
 *   Dual routine: JunctionTree::deReceiveToPreviousPartition()
 *
 *
 * Preconditions:
 *   It is assumed that:
 *     1) the right interface of the previous partition must have had
 *        all messages sent to it.
 * 
 * Postconditions:
 *     the left interface of the next partition is now set up.
 *
 * Side Effects:
 *     all partitions will have been instantiated to the extent that the messages (with
 *     the current pruning ratios)  have been created.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::ceSendToNextPartition(// previous partition
				    JT_InferencePartition& previous_part,
				    // root clique of the previous partition (i.e., the
				    // right interface clique) 
				    const unsigned previous_part_root,
				    // name of previous partition (for debugging/status msgs)
				    const char*const previous_part_type_name,
				    // sequence number (in unrolling) of previous partition
				    // (for debugging/status msgs)
				    const unsigned previous_part_num,
				    // next partition
				    JT_InferencePartition& next_part,
				    // leaf clique of next partition (i.e., index number
				    // of the left interface clique of next partition)
				    const unsigned next_part_leaf,
				    // name (debugging/status msgs)
				    const char*const next_part_type_name,
				    // partitiiton number (debugging/status msgs)
				    const unsigned next_part_num)
{
  // check for empty partitions.
  if (previous_part.maxCliques.size() == 0 || next_part.maxCliques.size() == 0)
    return;

  infoMsg(IM::Mod,"CE: message %s,part[%d],clique(%d) --> %s,part[%d],clique(%d)\n",
	   previous_part_type_name,
	   previous_part_num,
	   previous_part_root,
	   next_part_type_name,
	   next_part_num,
	   next_part_leaf);
  previous_part.maxCliques[previous_part_root].
    ceSendToOutgoingSeparator(previous_part,
			      next_part.
			      separatorCliques[next_part.separatorCliques.size()-1]);
}



/*-
 *-----------------------------------------------------------------------
 * JunctionTree::collectEvidence()
 *   
 *   Collect Evidence: This routine performes a complete collect
 *   evidence pass for this series of partitions that have been
 *   unrolled a certain amount, given by the jtIPartitions array. It
 *   sets the appropriate names, etc. of the partitions depending on
 *   if this is a left-interface or right-interface form of inference.
 *   the root within the current partition. 
 *  
 *   This routine demonstrates a simple version of linear space
 *   collect evidence inference. It keeps everything in memory at the
 *   same time when doing going forward, so it is suitable for use in
 *   a collect evidence/distribute evidence (forward/backward)
 *   framework.  A companion routine (aptly named
 *   'distributeEvidence') will, assuming collect evidence has been
 *   called, do the distribute evidence stage, thereby leaving all
 *   cliques in all partitions locally and therefore globally
 *   consistent, ready to be used say in EM training.
 *
 * See Also:
 *    1) contant memory (not dept. on time T) version of collect evidence
 *       JunctionTree::probEvidence()
 *    2) log space version of collect/distribute evidence  
 *       JunctionTree::collectDistributeIsland()
 *
 * Preconditions:
 *   The parititons must have been created and placed in the array jtIPartitions.
 *
 * Postconditions:
 *     All cliques in the partition have all messages but one sent to it.
 *
 * Side Effects:
 *     all partitions will have been instantiated to the extent that the messages (with
 *     the current pruning ratios)  have been created.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::collectEvidence()
{
  // This routine handles all of:
  // 
  //    unrolled 0 times: (so there is a single P1, and E1)  
  //    unrolled 1 time: so there is a P1, C1, C3, E1
  //    unrolled 2 or more times: so there is a P1 C1 [C2 ...] C3, E1
  //    etc.

  const unsigned numCoPartitions = jtIPartitions.size()-2;
  unsigned partNo = 0;
  // set up appropriate name for debugging output.
  const char* prv_nm;
  prv_nm = P1_n;
  unsigned prv_ri = P_ri_to_C;
  ceGatherIntoRoot(jtIPartitions[partNo],
		   P_ri_to_C,
		   P1_message_order,
		   prv_nm,partNo);

  for (unsigned p = 0; p < numCoPartitions; p++ ) {
    ceSendToNextPartition(jtIPartitions[partNo],prv_ri,prv_nm,partNo,
			  jtIPartitions[partNo+1],C_li_to_C,Co_n,partNo+1);
    partNo++;
    prv_nm = Co_n;
    prv_ri = C_ri_to_C;

    // we skip the first Co's LI separator if there is no P1
    // partition, since otherwise we'll get zero probability.
    if (p == 0 && P1.cliques.size() == 0)
      Co.skipLISeparator();
    ceGatherIntoRoot(jtIPartitions[partNo],
		     C_ri_to_C,
		     Co_message_order,
		     prv_nm,partNo);
    // if the LI separator was turned off, we need to turn it back on.
    if (p == 0 && P1.cliques.size() == 0)
      Co.useLISeparator();
  }

  ceSendToNextPartition(jtIPartitions[partNo],prv_ri,prv_nm,partNo,
			jtIPartitions[partNo+1],E_li_to_C,E1_n,partNo+1);
  partNo++;
  ceGatherIntoRoot(jtIPartitions[partNo],
		   E_root_clique,
		   E1_message_order,
		   E1_n,partNo);


  assert ( partNo == jtIPartitions.size()-1 );

  // root clique of last partition did not do partition, since it
  // never sent to next separator (since there is none). We explicitly
  // call pruning on the root clique of the last partition.
  if (jtIPartitions[partNo].maxCliques.size() > 0)
    jtIPartitions[partNo].maxCliques[E_root_clique].ceDoAllPruning();

}




/*-
 *-----------------------------------------------------------------------
 * JunctionTree::deScatterOutofRoot
 *   
 *   Distribute Evidence Scatter Outof Root: This routine does a
 *   distribute evidence pass for this partition, and scatters all
 *   messages outof the root within the current partition. It does so
 *   using the message order given in the argument 'message_order',
 *   and scatters out of the provided root clique (which is the right
 *   interface clique of this partition). By "Scatter", I mean it
 *   sends messages from the root clique distributing everything
 *   ultimately to all leaf cliques in this partition.
 *
 * See Also:
 *   Dual routine: JunctionTree::ceGatherIntoRoot()
 *
 *
 * Preconditions:
 *   It is assumed that either:
 *     1) this is the ritht-most partition
 *  or 2) that the right interface clique within this partition has had a message
 *        sent to it from the right neighbor partition.
 *
 * Postconditions:
 *     All cliques in the partition have all messages but one sent to it.
 *
 * Side Effects:
 *     all partitions will have been instantiated to the extent that the messages (with
 *     the current pruning ratios)  have been created.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::deScatterOutofRoot(// the partition
				 JT_InferencePartition& part,
				 // root (right interface clique) of this partition
				 const unsigned root,
				 // message order
				 vector< pair<unsigned,unsigned> >& message_order,
				 // name (debugging/status msgs)
				 const char*const part_type_name,
				 // partition number (debugging/status msgs)
				 const unsigned part_num)
{
  if (part.maxCliques.size() == 0)
    return;

  const unsigned stopVal = (unsigned)(-1);
  part.maxCliques[root].
    deScatterToOutgoingSeparators(part);
  for (unsigned msgNo=(message_order.size()-1);msgNo != stopVal; msgNo --) {
    const unsigned to = message_order[msgNo].first;
    const unsigned from = message_order[msgNo].second;
    infoMsg(IM::Mod,"DE: message %s,part[%d]: clique %d <-- clique %d\n",
	    part_type_name,part_num,to,from);
    part.maxCliques[to].
      deReceiveFromIncommingSeparator(part);
    part.maxCliques[to].
      deScatterToOutgoingSeparators(part);
  }
}




/*-
 *-----------------------------------------------------------------------
 * JunctionTree::deReceiveToPreviousPartition()
 *   
 *   Distribute Evidence Receive To Previous Partition: This routine
 *   sends a message from the left interface clique of a right (or
 *   next) partition to the right interface clique of a left (or
 *   previous) partition in the partition series. It is assumed that
 *   the left interface clique has had all its messages sent to it.
 *
 * See Also:
 *   Dual routine: JunctionTree::ceSendToNextPartition()
 *
 * Preconditions:
 *   It is assumed that:
 *     1) the left interface of the next partition must have had
 *        all messages sent to it.
 * 
 * Postconditions:
 *     the right interface of the previous partition is now set up.
 *
 * Side Effects:
 *     all partitions will have been instantiated to the extent that the messages (with
 *     the current pruning ratios)  have been created.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::deReceiveToPreviousPartition(// the next partition
					   JT_InferencePartition& next_part,
					   // leaf (left interface cliuqe) of next partition
					   const unsigned next_part_leaf,
					   // name of next part (debugging/status msgs)
					   const char*const next_part_type_name,
					   // number of next part (debugging/status msgs)
					   const unsigned next_part_num,
					   // previous partition
					   JT_InferencePartition& previous_part,
					   // root (right interface clique) of previous partition 
					   const unsigned previous_part_root,
					   // name of prev part (debugging/status msgs)
					   const char*const previous_part_type_name,
					   // number of prev part (debugging/status msgs)
					   const unsigned previous_part_num)
{
  // check for empty partitions.
  if (previous_part.maxCliques.size() == 0 || next_part.maxCliques.size() == 0)
    return;

  infoMsg(IM::Mod,"DE: message %s,part[%d],clique(%d) <-- %s,part[%d],clique(%d)\n",
	  previous_part_type_name,previous_part_num,previous_part_root,
	  next_part_type_name,next_part_num,next_part_leaf);
  previous_part.maxCliques[previous_part_root].
    deReceiveFromIncommingSeparator(previous_part,
				    next_part.
				    separatorCliques[next_part.
						     separatorCliques.size()-1]);
}



/*-
 *-----------------------------------------------------------------------
 * JunctionTree::distributeEvidence()
 *   
 *   Distribute Evidence: This routine performs a complete distribute
 *   evidence pass for this series of partitions that have been
 *   unrolled a certain amount, given by the jtIPartitions array. It
 *   sets the appropriate names, etc. of the partitions depending on
 *   if this is a left-interface or right-interface form of inference.
 *   the root within the current partition.
 *  
 *   This routine demonstrates a simple version of linear space
 *   distribute evidence inference. It uses data kept in memory when
 *   doing going backwards, and it is assumed that it will be used in
 *   tandem with a previous collect evidence call.  A companion
 *   routine ('collectEvidence') will, do collect evidence appropriately
 *   leaving all data structures set up for distributeEvidence() to be called.
 *
 *   After distributeEvidence() is called, all cliques in all
 *   partitions will be locally (& globally if it is a JT) consistant, and so will
 *   be ready for EM training, etc.
 *
 * See Also:
 *    0) collectEvidence()
 *    1) contant memory (not dept. on time T) version of collect evidence
 *       JunctionTree::probEvidence()
 *    2) log space version of collect/distribute evidence  
 *       JunctionTree::collectDistributeIsland()
 *
 * Preconditions:
 *   - collectEvidence() must have been called right before this.
 *   - The parititons must have been created and placed in the array jtIPartitions.
 * 
 *
 * Postconditions:
 *     All cliques are now locally consistant, and will be globally consistant
 *     if we have a junction tree. 
 *
 * Side Effects:
 *     all partitions will have been instantiated to the extent that the messages (with
 *     the current pruning ratios)  have been created.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::distributeEvidence()
{
  // this routine handles all of:
  // unrolled 0 times: (so there is a single P1, and E1)  
  // unrolled 1 time: so there is a P1, C1, C3, E1
  // unrolled 2 or more times: so there is a P1 C1 [C2 ...] C3, E1

  const unsigned numCoPartitions = jtIPartitions.size()-2;
  unsigned partNo = jtIPartitions.size()-1;
  // set up appropriate name for debugging output.
  const char* prv_nm;

  deScatterOutofRoot(jtIPartitions[partNo],
		     E_root_clique,
		     E1_message_order,
		     E1_n,partNo);

  partNo--;

  prv_nm = Co_n;
  deReceiveToPreviousPartition(jtIPartitions[partNo+1],E_li_to_C,E1_n,partNo+1,
			       jtIPartitions[partNo],C_ri_to_E,Co_n,partNo);

  for (unsigned p=numCoPartitions; p > 0; p--) {

    if (p == 1 && P1.cliques.size() == 0)
      Co.skipLISeparator();
    deScatterOutofRoot(jtIPartitions[partNo],
		       C_ri_to_C,
		       Co_message_order,
		       Co_n,partNo);
    if (p == 1 && P1.cliques.size() == 0)
      Co.useLISeparator();

    partNo--;
    prv_nm = ((p == 1)?P1_n:Co_n);
    unsigned prv_ri;
    if (p > 1) 
      prv_ri = C_ri_to_C;
    else
      prv_ri = P_ri_to_C;
    deReceiveToPreviousPartition(jtIPartitions[partNo+1],C_li_to_C,Co_n,partNo+1,
				 jtIPartitions[partNo],prv_ri,prv_nm,partNo);
  }

  deScatterOutofRoot(jtIPartitions[partNo],
		     P_ri_to_C,
		     P1_message_order,
		     P1_n,partNo);

  assert ( partNo == 0 );
}

/*
 * Must be called after collectEvidence() or distributeEvidence() has been called.
 *
 */
void
JunctionTree::printAllCliques(FILE* f,const bool normalize)
{
  unsigned partNo = 0;
  char buff[2048];
  if (pPartCliquePrintRange != NULL) {
    BP_Range::iterator it = pPartCliquePrintRange->begin();
    while (!it.at_end()) {
      const unsigned cliqueNum = (unsigned)(*it);
      if (cliqueNum < jtIPartitions[partNo].maxCliques.size()) {
	sprintf(buff,"Partition %d (P), Clique %d:",partNo,cliqueNum); 
	jtIPartitions[partNo].maxCliques[cliqueNum].printCliqueEntries(f,buff,normalize);
      }
      it++;
    }
  }
  partNo++;

  if (cPartCliquePrintRange != NULL) {
    for (;partNo<jtIPartitions.size()-1;partNo++) {
      BP_Range::iterator it = cPartCliquePrintRange->begin();
      while (!it.at_end()) {
	const unsigned cliqueNum = (unsigned)(*it);
	if (cliqueNum < jtIPartitions[partNo].maxCliques.size()) {
	  sprintf(buff,"Partition %d (C), Clique %d:",partNo,cliqueNum); 
	  jtIPartitions[partNo].maxCliques[cliqueNum].printCliqueEntries(f,buff,normalize);
	}
	it++;
      }
    }
  } else {
    partNo = jtIPartitions.size()-1;
  }
    
  if (ePartCliquePrintRange != NULL) {
    BP_Range::iterator it = ePartCliquePrintRange->begin();
    while (!it.at_end()) {
      const unsigned cliqueNum = (unsigned)(*it);
      if (cliqueNum < jtIPartitions[partNo].maxCliques.size()) {
	sprintf(buff,"Partition %d (E), Clique %d:",partNo,cliqueNum); 
	jtIPartitions[partNo].maxCliques[cliqueNum].printCliqueEntries(f,buff,normalize);
      }
      it++;
    }
  }
}



/*-
 *-----------------------------------------------------------------------
 * JunctionTree::emIncrement()
 *
 *    A version of emIncrement that works with the data structures set up by 
 *    collectEvidence() and distributeEvidence() above. Note, EM training
 *    can also be performed with the island algorithm.
 *
 * See Also:
 *    0) collectEvidence()
 *    1) distributeEvidence()
 *    2) log space version of collect/distribute evidence  
 *       JunctionTree::collectDistributeIsland()
 *
 * Preconditions:
 *    collectEvidence() AND distributeEvidence() must have just been
 *    called setting up the data structures. All cliques must exist and
 *    must have their clique tables filled out and ready.
 *    Also, all parametr accumulators must be set up and ready to go.
 *
 * Postconditions:
 *   The accumulators are increment accordingly.
 *
 * Side Effects:
 *   This will update the accumulators of all trainable parameter objects.
 *
 * Results:
 *   None
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::emIncrement(const logpr probE,
			  const bool localCliqueNormalization,
			  const double emTrainingBeam)
{
  // Quite simply, just iterate through all partitions and call emIncrement
  // therein.
  for (unsigned partNo = 0; partNo < jtIPartitions.size() ; partNo ++ ) {
    jtIPartitions[partNo].emIncrement(probE,localCliqueNormalization,emTrainingBeam);
  }
}


/*-
 *-----------------------------------------------------------------------
 * JunctionTree::probEvidence()
 *
 *    A constant memory (i.e., indep. of T), combination of unroll and
 *    collectEvidence() above. It is constnat memory in that
 *    it keeps no more than two partitions in memory simultaneously.
 *
 * See Also:
 *    0) collectEvidence()
 *    1) distributeEvidence()
 *    2) log space version of collect/distribute evidence  
 *       JunctionTree::collectDistributeIsland()
 *
 * Preconditions:
 *    same as unroll() above.
 *
 * Postconditions:
 *   The probability of the evidence is returned to the caller.
 *
 * Side Effects:
 *   This will update the space managers averages and statistics.
 *
 *   TODO: hash table usage, sharingetc. so that hash tables are not re-used 
 *         for the entire length.
 *
 * Results:
 *
 *-----------------------------------------------------------------------
 */
logpr 
JunctionTree::probEvidence(const unsigned int numFrames,
			   unsigned& numUsableFrames)
{

  // first create the unrolled set of random variables corresponding
  // to this JT.

  unsigned basicTemplateUnrollAmount;
  int modifiedTemplateUnrollAmount;
  unsigned frameStart;
  if (!gm_template.computeUnrollParameters(numFrames,
					   basicTemplateUnrollAmount,
					   modifiedTemplateUnrollAmount,
					   numUsableFrames,
					   frameStart))
    error("Segment of %d frames too short with current GMTK template of length [P=%d,C=%d,E=%d] %d frames, and M=%d,S=%d boundary parameters. Use longer utterances, different template, or decrease M,S if >1.\n",
	  numFrames,
	  fp.numFramesInP(),fp.numFramesInC(),fp.numFramesInE(),
	  fp.numFrames(),
	  gm_template.M,gm_template.S);


  infoMsg(IM::Info,"numFrames = %d, numUsableFrames = %d\n",numFrames,numUsableFrames);
  infoMsg(IM::Tiny,"numFrames = %d, unrolling BT %d times, MT %d times\n",
	  numFrames,
	  basicTemplateUnrollAmount,
	  modifiedTemplateUnrollAmount);

  // unrolled random variables
  // vector <RV*> unrolled_rvs;
  // mapping from 'name+frame' to integer index into unrolled_rvs.
  // map < RVInfo::rvParent, unsigned > ppf;

  fp.unroll(basicTemplateUnrollAmount,cur_unrolled_rvs,cur_ppf);
  setObservedRVs(cur_unrolled_rvs);

  // this clears the shared caches. 
  clearCliqueSepValueCache();


  // actual absolute part numbers
  unsigned partNo;
  const int numCoPartitions = modifiedTemplateUnrollAmount+1;
  // never use more than two partitions at a time.
  JT_InferencePartition *curPart, *prevPart;
  // set up appropriate name for debugging output.
  const char* prv_nm;
  unsigned prv_ri;

  partNo = 0;
  curPart = new JT_InferencePartition(P1,cur_unrolled_rvs,cur_ppf,0*gm_template.S);
  prevPart = NULL;
  prv_nm = P1_n;
  prv_ri = P_ri_to_C;
  ceGatherIntoRoot(*curPart,
		   P_ri_to_C,
		   P1_message_order,
		   prv_nm,partNo);
  swap(curPart,prevPart);

  for (int p = 0; p < numCoPartitions; p++ ) {
    delete curPart;
    curPart = new JT_InferencePartition(Co,cur_unrolled_rvs,cur_ppf,p*gm_template.S);
    ceSendToNextPartition(*prevPart,prv_ri,prv_nm,partNo,
			  *curPart,C_li_to_C,Co_n,partNo+1);
    partNo++;
    prv_nm = Co_n;
    prv_ri = C_ri_to_C;

    if (p == 0 && P1.cliques.size() == 0)
      Co.skipLISeparator();
    ceGatherIntoRoot(*curPart,
		     C_ri_to_C,
		     Co_message_order,
		     prv_nm,partNo);
    if (p == 0 && P1.cliques.size() == 0)
      Co.useLISeparator();

    swap(curPart,prevPart);
  }

  delete curPart;
  curPart = new JT_InferencePartition(E1,cur_unrolled_rvs,cur_ppf,
				      modifiedTemplateUnrollAmount*gm_template.S);
  ceSendToNextPartition(*prevPart,prv_ri,prv_nm,partNo,
			*curPart,E_li_to_C,E1_n,partNo+1);
  partNo++;
  ceGatherIntoRoot(*curPart,
		   E_root_clique,
		   E1_message_order,
		   E1_n,partNo);

  // root clique of last partition did not do partition, since it
  // never sent to next separator (since there is none). We explicitly
  // call pruning on the root clique of the last partition.
  curPart->maxCliques[E_root_clique].ceDoAllPruning();

  logpr rc = curPart->maxCliques[E_root_clique].sumProbabilities();
  
  delete curPart;
  delete prevPart;
  
  return rc;

}



/*-
 *-----------------------------------------------------------------------
 * JunctionTree::probEvidenceTime()
 *    Just like probEvidence() but keeps track of how many messages between partitions
 *    have been done, and returns (even if not complete) when timer expires.  
 *
 * See Also:
 *    0) probEvidence()
 *
 * Preconditions:
 *    same as probEvidence() above.
 *
 * Postconditions:
 *    same as probEvidence() above.
 *
 * Side Effects:
 *    same as probEvidence() above.
 *
 * Results:
 *    same as probEvidence() above.
 *-----------------------------------------------------------------------
 */
logpr 
JunctionTree::probEvidenceTime(const unsigned int numFrames,
			       unsigned& numUsableFrames,
			       unsigned& numPartitionsDone)
{

  // first create the unrolled set of random variables corresponding
  // to this JT.
  numPartitionsDone = 0;
  logpr res;

  unsigned basicTemplateUnrollAmount;
  int modifiedTemplateUnrollAmount;
  unsigned frameStart;
  if (!gm_template.computeUnrollParameters(numFrames,
					   basicTemplateUnrollAmount,
					   modifiedTemplateUnrollAmount,
					   numUsableFrames,
					   frameStart))
    error("Segment of %d frames too short with current GMTK template of length [P=%d,C=%d,E=%d] %d frames, and M=%d,S=%d boundary parameters. Use longer utterances, different template, or decrease M,S if >1.\n",
	  numFrames,
	  fp.numFramesInP(),fp.numFramesInC(),fp.numFramesInE(),
	  fp.numFrames(),
	  gm_template.M,gm_template.S);


  infoMsg(IM::Info,"numFrames = %d, numUsableFrames = %d\n",numFrames,numUsableFrames);
  infoMsg(IM::Tiny,"numFrames = %d, unrolling BT %d times, MT %d times\n",
	  numFrames,
	  basicTemplateUnrollAmount,
	  modifiedTemplateUnrollAmount);

  // unrolled random variables
  // vector <RV*> unrolled_rvs;
  // mapping from 'name+frame' to integer index into unrolled_rvs.
  // map < RVInfo::rvParent, unsigned > ppf;

  fp.unroll(basicTemplateUnrollAmount,cur_unrolled_rvs,cur_ppf);
  setObservedRVs(cur_unrolled_rvs);

  // this clears the shared caches. 
  clearCliqueSepValueCache();

  // actual absolute part numbers
  unsigned partNo;
  const unsigned numCoPartitions = modifiedTemplateUnrollAmount+1;
  // never use more than two partitions at a time.
  JT_InferencePartition *curPart, *prevPart;
  // set up appropriate name for debugging output.
  const char* prv_nm;
  unsigned prv_ri;

  partNo = 0;
  curPart = new JT_InferencePartition(P1,cur_unrolled_rvs,cur_ppf,0*gm_template.S);
  prevPart = NULL;
  prv_nm = P1_n;
  prv_ri = P_ri_to_C;
  ceGatherIntoRoot(*curPart,
		   P_ri_to_C,
		   P1_message_order,
		   prv_nm,partNo);
  swap(curPart,prevPart);
  if (probEvidenceTimeExpired)
    goto finished;

  for (unsigned p = 0; p < numCoPartitions; p++ ) {
    delete curPart;
    curPart = new JT_InferencePartition(Co,cur_unrolled_rvs,cur_ppf,p*gm_template.S);
    ceSendToNextPartition(*prevPart,prv_ri,prv_nm,partNo,
			  *curPart,C_li_to_C,Co_n,partNo+1);
    partNo++;
    prv_nm = Co_n;
    prv_ri = C_ri_to_C;

    if (p == 0 && P1.cliques.size() == 0)
      Co.skipLISeparator();
    ceGatherIntoRoot(*curPart,
		     C_ri_to_C,
		     Co_message_order,
		     prv_nm,partNo);
    if (p == 0 && P1.cliques.size() == 0)
      Co.useLISeparator();

    swap(curPart,prevPart);
    if (probEvidenceTimeExpired)
      goto finished;
  }

  delete curPart;
  curPart = new JT_InferencePartition(E1,cur_unrolled_rvs,cur_ppf,
				      modifiedTemplateUnrollAmount*gm_template.S);
  ceSendToNextPartition(*prevPart,prv_ri,prv_nm,partNo,
			*curPart,E_li_to_C,E1_n,partNo+1);
  partNo++;
  ceGatherIntoRoot(*curPart,
		   E_root_clique,
		   E1_message_order,
		   E1_n,partNo);

  // root clique of last partition did not do partition, since it
  // never sent to next separator (since there is none). We explicitly
  // call pruning on the root clique of the last partition.
  curPart->maxCliques[E_root_clique].ceDoAllPruning();

  res = curPart->maxCliques[E_root_clique].sumProbabilities();
 finished:
  numPartitionsDone = partNo;
  
  delete curPart;
  delete prevPart;

  return res;

}




/*-
 *-----------------------------------------------------------------------
 * JunctionTree:: interface routines for island algorithm.
 *
 *    The various routines are interface routines for island algorithm
 *    to the lower level collect/distribute evidence routines defined
 *    above. These routines each only know the current partition (and
 *    possibly the previous or next partition) number, take the actual
 *    partition information from a pre-allocated array, and call
 *    the appropriate real routine. The routines included here
 *    are:
 *        ceGatherIntoRoot - call CE into the root (RI) of a partition
 *        createPartition - create the partition at location given
 *        ceSendToNextPartition - send from RI of left partition to LI of next right partition.
 *        cePruneRootCliqueOfPartition - explicitly prune the root cliuqe of the partition 
 *        deReceiveToPreviousPartition - send from LI of right partition to RI of previous left partition
 *        deletePartition - free memory associated with partition, and set array ptr to NULL
 *        deScatterOutofRoot - call DE from root (RI) of a partition
 *        probEvidenceRoot - return the prob of evidence from the root clique of the partition
 *        emIncrementIsland - name says it all
 *        printAllCliques - print out all clique entries according to clique ranges.
 *
 *    Note that these routines only do real work if the current prob evidence is
 *    something other than zero. As if it is zero, that means that
 *    we've already gotten to the end and found it is zero, and we're
 *    in the mode where we're just freeing up memory.
 *
 * Preconditions:
 *   For the non create/destroy routines, each routine assumes that
 *   the array location for the corresponding partition has been
 *   assigned appropriately. The other pre-conditions are the
 *   same as the routines (above) that are called. 
 *   The create routine has the same preconditions.
 *   The destroy routine assumes that the partition at the current location
 *   has been created.
 *
 * Postconditions:
 *   Since these routiens are mainly stubs, see the corresponding function that is called.
 *
 * Side Effects:
 *   These routines will affect internal structures, space managers, hash tables
 *   within the partitions that they refer to.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void JunctionTree::
ceGatherIntoRoot(const unsigned part)
{
  infoMsg(IM::Mod,"==> ceGatherIntoRoot: part = %d, nm = %s\n",
	  part,partPArray[part].nm);

  if (part == 1 && P1.cliques.size() == 0)
    Co.skipLISeparator();
  if (cur_prob_evidence.not_essentially_zero())
    ceGatherIntoRoot(*partPArray[part].p,
		     partPArray[part].ri,
		     *partPArray[part].mo,
		     partPArray[part].nm,
		     part);
  if (part == 1 && P1.cliques.size() == 0)
    Co.useLISeparator();

}
void JunctionTree::
createPartition(const unsigned part)
{
  infoMsg(IM::Mod,"$$$ createPartition: part = %d, nm = %s\n",
	  part,partPArray[part].nm);
  if (cur_prob_evidence.not_essentially_zero())
    partPArray[part].p = new JT_InferencePartition(*partPArray[part].JT,
						   cur_unrolled_rvs,cur_ppf,
						   partPArray[part].offset);
}
void JunctionTree::
ceSendToNextPartition(const unsigned part,const unsigned nextPart)
{
  infoMsg(IM::Mod,"--> ceSendToNextPartition: part = %d (%s), nextPart = %d (%s)\n",
	  part,partPArray[part].nm,
	  nextPart,partPArray[nextPart].nm);
  if (cur_prob_evidence.not_essentially_zero())
    ceSendToNextPartition(*partPArray[part].p,partPArray[part].ri,partPArray[part].nm,part,
			  *partPArray[nextPart].p,partPArray[nextPart].li,partPArray[nextPart].nm,nextPart);
}
void JunctionTree::
cePruneRootCliqueOfPartition(const unsigned part)
{
  infoMsg(IM::Mod,"\\/- cePruneRootCliqueOfPartition: part = %d, nm = %s\n",
	  part,partPArray[part].nm);
  if (cur_prob_evidence.not_essentially_zero())
    partPArray[part].p->maxCliques[partPArray[part].ri].ceDoAllPruning();
}
void JunctionTree::
deReceiveToPreviousPartition(const unsigned part,const unsigned prevPart)
{
  infoMsg(IM::Mod,"<-- deReceiveToPreviousPartition: part = %d (%s), prevPart = %d (%s)\n",
	  part,partPArray[part].nm,
	  prevPart,partPArray[prevPart].nm);
  if (cur_prob_evidence.not_essentially_zero())
    deReceiveToPreviousPartition(*partPArray[part].p,partPArray[part].li,partPArray[part].nm,part,
				 *partPArray[prevPart].p,partPArray[prevPart].ri,partPArray[prevPart].nm,prevPart);
}
void JunctionTree::
deletePartition(const unsigned part)
{
  infoMsg(IM::Mod,"*** deletePartition: part = %d (%s)\n",
	  part,partPArray[part].nm);
  delete partPArray[part].p;
  partPArray[part].p = NULL;
}
void JunctionTree::
deScatterOutofRoot(const unsigned part)
{
  infoMsg(IM::Mod,"<== deScatterOutofRoot: part = %d (%s)\n",
	  part,partPArray[part].nm);

  if (part == 2 && P1.cliques.size() == 0)
    Co.skipLISeparator();
  if (cur_prob_evidence.not_essentially_zero())
    deScatterOutofRoot(*partPArray[part].p,
		       partPArray[part].ri,
		       *partPArray[part].mo,
		       partPArray[part].nm,
		       part);
  if (part == 2 && P1.cliques.size() == 0)
    Co.useLISeparator();


}
logpr
JunctionTree::probEvidenceRoot(const unsigned part)
{
  // return the sum of probs for the root (right interface) clique of the given partition.
  infoMsg(IM::Mod,"^^^ computing evidence for JT root: part = %d (%s)\n",
	  part,partPArray[part].nm);
  return partPArray[part].p->maxCliques[partPArray[part].ri].sumProbabilities();
}
logpr
JunctionTree::setRootToMaxCliqueValue(const unsigned part)
{
  // return the sum of probs for the root (right interface) clique of the given partition.
  infoMsg(IM::Mod,"^^^ setting JT root to max clique value: part = %d (%s)\n",
	  part,partPArray[part].nm);
  return partPArray[part].p->maxCliques[partPArray[part].ri].setCliqueToMaxCliqueValue();
}
void
JunctionTree::emIncrementIsland(const unsigned part,
				const logpr cur_prob_evidence,
				const bool localCliqueNormalization)
{
  // increment for this partition.
  infoMsg(IM::Mod,"^^^ incrementing EM: part = %d (%s)\n",
	  part,partPArray[part].nm);
  return partPArray[part].p->emIncrement(cur_prob_evidence,localCliqueNormalization,curEMTrainingBeam);
}
void
JunctionTree::printAllCliques(const unsigned part,FILE* f,const bool normalize)
{
  BP_Range* rng;
  if (part == 0) 
    rng = pPartCliquePrintRange;
  else if (part == partPArray.size()-1)
    rng = ePartCliquePrintRange;
  else 
    rng = cPartCliquePrintRange;
  char buff[2048];
  if (rng != NULL) {
    BP_Range::iterator it = rng->begin();
    while (!it.at_end()) {
      const unsigned cliqueNum = (unsigned)(*it);
      sprintf(buff,"Partition %d (%s), Clique %d:",part,partPArray[part].nm,cliqueNum); 
      partPArray[part].p->maxCliques[cliqueNum].printCliqueEntries(f,buff,normalize);
      it++;
    }
  }
}





/*-
 *-----------------------------------------------------------------------
 * JunctionTree::collectDistributeIslandBase()
 *
 * Preconditions:
 *   1) Guaranteed that partition[start] has
 *   been allocated. If start > 0, then we are also
 *   guaranteed that ceSendNext(start-1,start) has
 *   been called. In either case, we are also guaranteed
 *   taht ceGatherIntoRoot(start) has been called.
 *
 *   2) Also, we are guaranteed that either end == final, or that
 *   partition[end+1] has also already been created and that we've
 *   already sent a forward message to partition[end+1]. Therefore, if
 *   end < final, we are ready to do a deReceiveToPrev(end+1,end) to
 *   the partition at position end.
 *
 *   This function should only be called by collectDistributeIslandRecurse()
 *
 * Postconditions:
 *
 *  All partitions between start and end will exist and will be
 *  complete, meaning they will have had *all* messages sent, and so
 *  they could be used for EM updating. All partitions
 *  from [start+1,end] will be deallocated. 
 *
 *
 * Side Effects:
 *
 * Results:
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::collectDistributeIslandBase(const unsigned start,
					  const unsigned end,
					  const bool runEMalgorithm,
					  const bool runViterbiAlgorithm,
					  const bool localCliqueNormalization)
{
  for (unsigned part = start; part <= end; part ++) {
    if (part > start) {
      // first one is done already (see pre-conditions) , so don't do.
      ceGatherIntoRoot(part);
    }
    if (part != end) {
      // then we create the next partition and send starting message
      // to it.
      createPartition(part+1);
      ceSendToNextPartition(part,part+1);
    } else {
      // We are at the end, so we don't create and send to the next
      // partition. This is because either 1)
      // end=(partPArray.size()-1) (we are really at the last
      // partition on the right), so there is no right neighboring
      // partition, or 2) the right neighboring partition has already
      // been sent a message from the previous
      // incarnation/instantiation of the current partition.  Note,
      // that if we don't send a message to the next partition, we do
      // need to prune last clique of this partition here since that
      // is what would have happened ordinarily. 
      // All the partition's root cliques need to be pruned, as
      // otherwise there might be entries in that clique not contained
      // in its previously-created outgoing separator. Therefore, we
      // explicitly prune here in that case.
      cePruneRootCliqueOfPartition(part);
    }
  }
  for (unsigned part = end; (1);) { 
    if (part < (partPArray.size()-1)) {
      // then there is something on the right to receive a message
      // from. 
      deReceiveToPreviousPartition(part+1,part);
    } else {
      // this is the true end and we can get the probability of evidence here.
      // only initialize this if we are not doing localCliqueNormalization
      if (runEMalgorithm) {
	cur_prob_evidence = probEvidenceRoot(part);
	if (cur_prob_evidence.essentially_zero()) {
	  infoMsg(IM::Default,"Island not training segment since probability is essentially zero\n");
	  // note that we can't just jump out as we have to free up all the
	  // memory that we allocated. We thus have to check a bunch of conditions on 
	  // the way out and do EM training only when appropriate, but always delete.
	}
      } else if (runViterbiAlgorithm) {
	cur_prob_evidence = probEvidenceRoot(part);
	if (cur_prob_evidence.essentially_zero()) {
	  infoMsg(IM::Default,"Island not decoding segment since probability is essentially zero\n");
	  // note that we can't just jump out as we have to free up
	  // all the memory that we allocated. We thus have to check a
	  // bunch of conditions on the way out and do decoding only
	  // when appropriate, but always delete.
	} else
	  setRootToMaxCliqueValue(part);
      }
      if (IM::messageGlb(IM::Low)) {
	infoMsg(IM::Low,"XXX Island Finished Inference: part = %d (%s): log probE = %f\n",
		part,partPArray[part].nm,probEvidenceRoot(part).valref());
      }
    }
    // 
    if (part != end) {
      // delete it if we created it here.
      deletePartition(part+1);
    }

    if ( ((runEMalgorithm == runViterbiAlgorithm) && (runEMalgorithm == false))
	 || cur_prob_evidence.not_essentially_zero()) {
      // scatter into all cliques within this separator (at the very least)
      deScatterOutofRoot(part);
      // We now have a completed partition that we can use
      // for EM, viterbi decoding, scoring, etc.
      if (IM::messageGlb(IM::Mod)) {
	infoMsg(IM::Mod,"!!! finished partition: part = %d (%s)\n",
		part,partPArray[part].nm);
	for (unsigned cliqueNo=0;cliqueNo<partPArray[part].p->maxCliques.size();cliqueNo++) {
	  printf("XXX Island: Part no %d: clique no %d: log probE = %f\n",
		 part,cliqueNo,partPArray[part].p->maxCliques[cliqueNo].sumProbabilities().valref());
	}
      }
    }
    
    // and do em updating if appropriate.
    if (runEMalgorithm && cur_prob_evidence.not_essentially_zero()) {
      emIncrementIsland(part,cur_prob_evidence,localCliqueNormalization);
    }
    printAllCliques(part,stdout,true);

    if (part == start)
      break;
    part--;
  }
}



/*-
 *-----------------------------------------------------------------------
 * JunctionTree::collectDistributeIslandRecurse()
 *
 *  Run a collect distribute evidence stage
 *  between the partitions that are at locations [start,end] INCLUSIVE.
 *  Note if this is the very first call, it should be called
 *  with [start=0,end=(nparts-1)].
 *
 * Preconditions:
 *    For this to work, we must have that:
 *     1) partitions[start] exists (i.e., allocated), and
 *        ceGather(start) has already been called.
 *     2) Either end == partitions[final], 
 *        Or,  end != partitions[final], and
 *             partitions[end+1] has been allocated, and
 *             it is all ready to have deReceiveToPrev(end+1,end)
 *             be called.
 *    This function should only be called by collectDistributeIsland()
 *
 * Postconditions:
 *  All partitions between start and end will have at one time
 *  been in a state where all cliques will have had all
 *  messages, and so could be used for EM updating.
 *
 * Side Effects:
 *
 * Results:
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::collectDistributeIslandRecurse(const unsigned start,
					     const unsigned end,
					     const unsigned base,
					     const unsigned linear_section_threshold,
					     const bool runEMalgorithm,
					     const bool runViterbiAlgorithm,
					     const bool localCliqueNormalization)
{
  // We're doing from [start,end] inclusive, so compute length
  // accordingly
  const unsigned len = end-start + 1;
  if (len <= linear_section_threshold) {
    // do base case.
    collectDistributeIslandBase(start,end,runEMalgorithm,
				runViterbiAlgorithm,
				localCliqueNormalization);
  } else { 
    const unsigned section_size = len/base;

    if (section_size <= 1) {
      infoMsg(IM::High,"Island collect/distribute inference, log base (%d) too large for current sect length (%d) & linear sect threshold (%d), it would result in sub sect len (%d). Backing off to linear case.\n",base,len,linear_section_threshold,section_size);
      return collectDistributeIslandBase(start,end,runEMalgorithm,
					 runViterbiAlgorithm,
					 localCliqueNormalization);
    }
    // so we are assured there that section_size is at least two here. 

    // compute what is left over at end.
    const unsigned remainder = len - section_size*base;

    unsigned section_start = start;
    // number of sections that get an extra partition.
    unsigned num_to_distribute = remainder;
    // The size of section is stored in this variable.
    unsigned cur_section_size;
    while (1) {
      // Update the current section size based on a uniform
      // allocation of the remaining partitions over the
      // first 'remainder' sections.
      cur_section_size = section_size;
      if (num_to_distribute > 0) {
	// distribute the remainder evenly over the sections.
	cur_section_size ++;
	num_to_distribute--;
      }
      // don't do last section, as is handled by recursive call.
      if (section_start + cur_section_size == (end+1))
	break;

      // the last partition included within section.
      const unsigned section_end = section_start + cur_section_size-1;
      // At this point, we are Guaranteed that:
      //  - partition[section_start] has been allocated and ahs
      //    already had ceGatherIntoRoot(section) called.
      // We now run a constant space collect evidence stage, where we:
      //    a) leave partition[section_start] in place
      //    b) delete all partitions starting from and including (section_start+1)
      //       up to and including section_end.
      //    c) the last section we allocate (section_end+1) gets
      //       a ceSendTo message, and is also gathered (since it is
      //       going to be the start for the recursive call), and we leave
      //       that one allocated in place.
      createPartition(section_start+1);
      ceSendToNextPartition(section_start,section_start+1);
      ceGatherIntoRoot(section_start+1);
      for (unsigned part = section_start+1; part <= section_end; part ++) {
	createPartition(part+1);
	ceSendToNextPartition(part,part+1);
	ceGatherIntoRoot(part+1);
	deletePartition(part);
      }
      // move to point to next partition.
      section_start += cur_section_size;
    }

    unsigned num_not_to_distribute = base - remainder;
    section_start = (end+1);
    bool first_iteration = true;
    while (1) {
      cur_section_size = section_size;
      if (num_not_to_distribute > 0) 
	num_not_to_distribute--;
      else {
	// distribute remainder over the first sections
	cur_section_size++;
      }
      section_start -= cur_section_size;
      // the last partition included within section.
      const unsigned section_end = section_start + cur_section_size-1;
      // recurse
      collectDistributeIslandRecurse(section_start,section_end,
				     base,linear_section_threshold,
				     runEMalgorithm,
				     runViterbiAlgorithm,
				     localCliqueNormalization);
      // need to delete partition at location section_start+cur_section_size
      // if it is one that we created.
      if (!first_iteration) {
	assert (section_start+cur_section_size != (end+1));
	deletePartition(section_start+cur_section_size);
      }
      first_iteration = false;
      // if this was the first section, we end now.
      if (section_start == start)
	break;
    }
  }
}



/*-
 *-----------------------------------------------------------------------
 * JunctionTree::collectDistributeIsland()
 *
 *  Do a collect evidence/distribute evidence pass using the island
 *  algorithm, a log-space version of collect/distribute evidence.
 *  The original idea was presented in the paper:
 *
 *      J. Binder, K. Murphy, and S. Russel. Space-efficient inference in
 *      dynamic probabilistic networks. In IJCAI, 1997.
 *      http://citeseer.nj.nec.com/30635.html
 *
 *  By word of mouth, I (bilmes) heard that the idea arose from an
 *  algorithm that is used in genetic sequencing, and the connection
 *  was originally suggested by Paul Horton, but I am not sure of the
 *  original source in genetics (as of 2004. Please let me know if you
 *  are reading this and you know.).
 *
 *  Here, we adopt this idea to GMTK's notion of graph partitions.
 *
 *  Most simply, in the linear case we pay O(T) memory and O(T)
 *  compute. In the island case, we pay O(b*log_b(T)) memory and
 *  O(Tlog_b(T)) compute, where b is the base of the logarithm (note,
 *  all costs are of course multiplied by the average within-partition
 *  cost, we explain things only interms of the time/space complexity
 *  in how it relates to the length of the segment T). This can make
 *  the difference between being able to use collect/distribute
 *  evidence if in going from O(T) to O(logT) we go from *not* being
 *  able to fit in main memory to *being* able to fit in memory, but
 *  the additional time cost of logT can still be quite large. If T =
 *  1024, say, and base = 2, then we pay an extra time cost of a
 *  factor of 10!! (e.g., 1 day to 10 days).
 *
 *  More specifically, for a segment of length T, this algorithm never
 *  keeps simultaneously in memory more than about M partitions, where
 *
 *            M = T/base^k + k*(base-1)
 *
 *  and where k is the *smallest* integer value such that
 *
 *            T/base^k <= lst
 *
 *  and where lst = linear_section_threshold (the threshold at which
 *  we drop down to linear inference).  Therefore, we keep
 *  in memory no more than about
 *
 *            M = lst + k*(base-1)
 *
 *  partitions, where k is defined as above, log_b(T/lst) = k. For
 *  example, if lst == 3, and base == 3, then we keep simultaneously
 *  in memory no more than about 3+(base-1)log3(T/3) partitions,
 *  having logarithmic growth in T.
 *
 *  As mentioned above, this doesn't come for free, however, as we pay
 *  a cost in time complexity to save memory. Specifically, we do
 *  multiple collect evidence stages (unfortunately, the more time
 *  costly of the two, collect vs. distribute). Specifically, we need
 *  to do T' = R(T) collect evidence stages between partitions, where
 * 
 *        R(T) = T + (base-1)*R(T/base), 
 *
 *  This is the recurance relationship corresponding to what is
 *  implemented here in GMTK. To solve it, however, we can simplify by
 *  saying
 *  
 *        R(T) < T + base*R(T/base) 
 *
 *  meaning R(T) < T + T + ... + T, and there are log_{base}(T) terms
 *  in the sum. In actuality, there are k terms in the sum, where k is
 *  defined as above.
 * 
 *  Therefore, we can decrease running time and increase space
 *  requirements either by increasing base (from 2 on up) or
 *  alternatively (or also) by increasing lst (from 2 on up). Note
 *  that we increase mem in increasing b since O(b*log_b(T)) =
 *  O((b/ln(b))*ln(T)), so mem is growing as b/ln(b). Note also that
 *  one should set base and lst such that everything just fits in main
 *  memory (i.e., so we don't start swapping to disk), but it is not
 *  worth it to set these so that it takes any less than main memory,
 *  as that will slow things down further than necessary. Also note
 *  that argmin_{ b \in 2,3,4,... } b/ln(b) = 3, so a base of 3 is
 *  a good starting point. 
 *
 * See Also:
 *    0) collectEvidence()
 *    1) distributeEvidence()
 *    2) probEvidence()
 *
 * Preconditions:
 *   Same as unroll()
 *
 * Postconditions:
 *   Forward/backward has been run.
 *
 * Side Effects:
 *  Updates space manager statitics in cliques. 
 *  TODO: add side effects here as routine evolves.
 *
 * Results:
 *
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::collectDistributeIsland(// number of frames in this segment.
				      const unsigned int numFrames,
				      // return value, number of frames
				      // that are actually used.
				      unsigned& numUsableFrames,
				      // the base of the logarithm
				      const unsigned base,
				      // the threshold at which we drop
				      // down to the linear collect/distribute
				      // evidence stage.
				      const unsigned linear_section_threshold,
				      const bool runEMalgorithm,
				      const bool runViterbiAlgorithm,
				      const bool localCliqueNormalization)
{

  // cant run both EM and viterbi at the same time.
  assert (!runEMalgorithm || !runViterbiAlgorithm);

  // must have a linear_section_threshold of at least two partitions.
  if (linear_section_threshold < 2)
    error("ERROR: Island algorithm collect/distribute inference. linear section threshold value (%d) is too small.\n",linear_section_threshold);


  // the log base must be a number that actually causes a split.
  if (base <= 1)
    error("ERROR: Island algorithm collect/distribute inference. base of log (%d) is too small.\n",base);


  // first create the unrolled set of random variables corresponding
  // to this JT.
  unsigned basicTemplateUnrollAmount;
  int modifiedTemplateUnrollAmount;
  unsigned frameStart;
  if (!gm_template.computeUnrollParameters(numFrames,
					   basicTemplateUnrollAmount,
					   modifiedTemplateUnrollAmount,
					   numUsableFrames,
					   frameStart))
    error("Segment of %d frames too short with current GMTK template of length [P=%d,C=%d,E=%d] %d frames, and M=%d,S=%d boundary parameters. Use longer utterances, different template, or decrease M,S if >1.\n",
	  numFrames,
	  fp.numFramesInP(),fp.numFramesInC(),fp.numFramesInE(),
	  fp.numFrames(),
	  gm_template.M,gm_template.S);

  infoMsg(IM::Info,"numFrames = %d, numUsableFrames = %d\n",numFrames,numUsableFrames);
  infoMsg(IM::Tiny,"numFrames = %d, unrolling BT %d times, MT %d times\n",
	  numFrames,
	  basicTemplateUnrollAmount,
	  modifiedTemplateUnrollAmount);

  fp.unroll(basicTemplateUnrollAmount,cur_unrolled_rvs,cur_ppf);
  setObservedRVs(cur_unrolled_rvs);

  const unsigned numPartitions = modifiedTemplateUnrollAmount+3;
  const unsigned numCoPartitions = modifiedTemplateUnrollAmount+1;
  unsigned partNo = 0;  

  // this clears the shared caches. 
  clearCliqueSepValueCache();

  // TODO: turn this array into a function so we don't alloate the entire array length.
  // re-allocate.
  partPArray.resize(numPartitions);
  // Redundantly store partition information in an array both for
  // debug messages and easy access to *significantly* simplify the
  // island code above.
  partPArray[partNo].JT = &P1;
  partPArray[partNo].offset = 0*gm_template.S;
  partPArray[partNo].p = NULL;
  partPArray[partNo].mo = &P1_message_order;
  partPArray[partNo].nm = P1_n;
  partPArray[partNo].ri = P_ri_to_C;
  partPArray[partNo].li = ~0x0;
  partNo++;
  for (unsigned p=0;p<numCoPartitions;p++) {
    partPArray[partNo].JT = &Co;
    partPArray[partNo].offset = p*gm_template.S;
    partPArray[partNo].p = NULL;
    partPArray[partNo].mo = &Co_message_order;
    partPArray[partNo].nm = Co_n;
    partPArray[partNo].ri = C_ri_to_C;
    partPArray[partNo].li = C_li_to_C;
    partNo++;
  }
  partPArray[partNo].JT = &E1;
  partPArray[partNo].offset = (modifiedTemplateUnrollAmount)*gm_template.S;
  partPArray[partNo].p = NULL;
  partPArray[partNo].mo = &E1_message_order;
  partPArray[partNo].nm = E1_n;
  partPArray[partNo].ri = E_root_clique;
  partPArray[partNo].li = E_li_to_C;
  partNo++;

  assert (partNo == partPArray.size());

  // Start off with unity probability of evidence (we can set it to be
  // anything other than zero actually). This will signal to the
  // island algorithm to keep going and do real work. If we reach the
  // end and if the segment decodes with zero probability, this
  // variable will be set to such, and which will cause the island
  // algorithm to, during its backward phase, just free up the islands
  // of memory that have been allocated rather than doing anything
  // else.
  cur_prob_evidence.set_to_one();

  // The recursion assumes that its first partition is already
  // allocated, so we make sure to do that here.
  partPArray[0].p = new JT_InferencePartition(P1,cur_unrolled_rvs,cur_ppf,0*gm_template.S);
  ceGatherIntoRoot(0);
  collectDistributeIslandRecurse(0,partPArray.size()-1,base,linear_section_threshold,
				 runEMalgorithm,
				 runViterbiAlgorithm,
				 localCliqueNormalization);
  delete partPArray[0].p;
  partPArray[0].p = NULL;

  partPArray.clear();
}


/*-
 *-----------------------------------------------------------------------
 * JunctionTree::probEvidence()
 *   
 *   This version of probEvidence() merely sums up the entries in E's
 *   root clique and return the result. This, if collect evidence has
 *   been called all the way to E's root clique, this function will
 *   return the probabilty of the evidence prob(E).
 *
 * See Also:
 *   probEvidence() above, the other (overloaded) version which is
 *   a const. mem version of collectEvidence().
 *
 * Preconditions:
 *   - collectEvidence() (or the other probEvidence()) must have been called
 *   - The E partition must exist in the jtIPartitions array
 *     and must have been called for this to work right.
 *
 * Postconditions:
 *   none
 *
 * Side Effects:
 *   none 
 *
 * Results:
 *   The probability of the evidence
 *
 *-----------------------------------------------------------------------
 */
logpr
JunctionTree::probEvidence()
{

  if (jtIPartitions[jtIPartitions.size()-1].maxCliques.size() > 0)
    return jtIPartitions[jtIPartitions.size()-1].maxCliques[E_root_clique].
      sumProbabilities();
  else {
    // this means that E is empty, we use the root clique of the last
    // Co partition.
    return jtIPartitions[jtIPartitions.size()-2].maxCliques[C_ri_to_E].
      sumProbabilities();
  }
}


/*-
 *-----------------------------------------------------------------------
 * JunctionTree::setRootToMaxCliqueValue()
 *   
 *   This version of probEvidence() merely sets E's root clique to its
 *   current max value, and return the max value.
 *
 * See Also:
 *   probEvidence() above, the other (overloaded) version which is
 *   a const. mem version of collectEvidence().
 *
 * Preconditions:
 *   - collectEvidence() (or the other probEvidence()) must have been called
 *     for this to work.
 *   - The E partition must exist in the jtIPartitions array
 *     and must have been called for this to work right.
 *
 * Postconditions:
 *   none
 *
 * Side Effects:
 *   none 
 *
 * Results:
 *   The probability of the evidence
 *
 *-----------------------------------------------------------------------
 */
logpr
JunctionTree::setRootToMaxCliqueValue()
{
  if (jtIPartitions[jtIPartitions.size()-1].maxCliques.size() > 0)
    return jtIPartitions[jtIPartitions.size()-1].maxCliques[E_root_clique].
      setCliqueToMaxCliqueValue();
  else 
    return jtIPartitions[jtIPartitions.size()-2].maxCliques[C_ri_to_E].
      setCliqueToMaxCliqueValue();
}


/*-
 *-----------------------------------------------------------------------
 * JunctionTree::printAllCliquesProbEvidence()
 *   
 *   This routine will cycle through all cliques of all partitions and will
 *   print out the sums of the probs. of each cliuqe. Therefore, if 
 *   collect/distribute evidence has been called (and it is working)
 *   this routine should print out exactly the same value for all
 *   cliques (to within numerical precision and logp.h table error roundoff).
 *
 *   This routine really is only used for debugging/status messages.
 *
 * Preconditions:
 *   - collectEvidence() and distributeEvidence() should have been called
 *     and the jtIPartitions array left set up.
 *
 * Postconditions:
 *   none
 *
 * Side Effects:
 *   none 
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::printAllCliquesProbEvidence()
{
  for (unsigned part=0;part<jtIPartitions.size();part++) {
    for (unsigned cliqueNo=0;cliqueNo<jtIPartitions[part].maxCliques.size();cliqueNo++) {
      printf("Part no %d: clique no %d: log probE = %f\n",
	     part,cliqueNo,jtIPartitions[part].maxCliques[cliqueNo].sumProbabilities().valref());
    }
  }
}




/*-
 *-----------------------------------------------------------------------
 * JunctionTree::printCurrentRVValues()
 *   
 *   Prints out all random variables in the current unrolled set.
 *
 * Preconditions:
 *   - cur_unrolled_rvs must be defined and unrolled
 *
 * Postconditions:
 *   cur_unrolled_rvs is printed.
 *
 * Side Effects:
 *   none 
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::printCurrentRVValues(FILE *f)
{
  for (unsigned i=0;i<cur_unrolled_rvs.size();i++) {
    cur_unrolled_rvs[i]->printSelf(f);
  }
}

void
JunctionTree::setCliquePrintRanges(char *p,char*c,char*e)
{
  if (p != NULL) {
    BP_Range* tmp = new BP_Range(p,0,gm_template.P.cliques.size());
    if (!tmp || tmp->length() <= 0)
      delete tmp;
    else
      pPartCliquePrintRange = tmp;
  }
  if (c != NULL) {
    BP_Range* tmp = new BP_Range(c,0,gm_template.C.cliques.size());
    if (!tmp || tmp->length() <= 0)
      delete tmp;
    else
      cPartCliquePrintRange = tmp;
  }
  if (e != NULL) {
    BP_Range* tmp = new BP_Range(e,0,gm_template.E.cliques.size());
    if (!tmp || tmp->length() <= 0)
      delete tmp;
    else
      ePartCliquePrintRange = tmp;
  }
}
