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
#include "GMTK_RandomVariable.h"
#include "GMTK_DiscreteRandomVariable.h"
#include "GMTK_ContinuousRandomVariable.h"
#include "GMTK_GM.h"
#include "GMTK_GMTemplate.h"
#include "GMTK_JunctionTree.h"
#include "GMTK_GMParms.h"



////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Static variables used by classes
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

VCID("$Header$");

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Misc Support
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

#if 0
// TODO: put this function somewhere more generally available.
static void
printRVSet(set<RandomVariable*>& locset)
{
  set<RandomVariable*>::iterator it;
  bool first = true;
  for (it=locset.begin();it!=locset.end();it++) {
    RandomVariable* rv = (*it);
    if (!first)
      printf(", ");
    printf("%s(%d)",rv->name().c_str(),rv->frame());
    first = false;
  }
  printf("\n");
}
#endif

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

JT_Partition::JT_Partition(
		       Partition& from_part,
		       vector <RandomVariable*>& newRvs,
		       map < RVInfo::rvParent, unsigned >& ppf,
		       const unsigned int frameDelta)
{

  triMethod = from_part.triMethod;

  set<RandomVariable*>::iterator it;

  // clone over nodes RVs.  
  // TODO: make this next code a routine
  //  nodesClone() since it is used in several places.
  for (it = from_part.nodes.begin();
       it != from_part.nodes.end();
       it++) {
    RandomVariable* rv = (*it);
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    

    if ( ppf.find(rvp) == ppf.end() ) {
      error("INTERNAL ERROR: can't find rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
    }

    RandomVariable* nrv = newRvs[ppf[rvp]];
    nodes.insert(nrv);
  }
  cliques.reserve(from_part.cliques.size());
  // 
  // NOTE: It is ***CRUCIAL*** for the cliques in the cloned partition
  // to be inserted in the *SAME ORDER* as in the partition being
  // cloned.
  for (unsigned i=0;i<from_part.cliques.size();i++) {
    cliques.push_back(MaxClique(from_part.cliques[i],
				newRvs,ppf,frameDelta));
  }
}


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        JT Inference Partition support
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

JT_InferencePartition::JT_InferencePartition(JT_Partition& from_part,
					     vector <RandomVariable*>& newRvs,
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

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Support for building a tree from clique graph
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * JunctionTree::createPartitionJunctionTree()
 *   Create a mini-junction tree from the cliques in the given partition.
 *   This uses Kruskal's greedy (but optimal) algorithm for MST generation.
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
JunctionTree::createPartitionJunctionTree(Partition& part)
{
  const unsigned numMaxCliques = part.cliques.size();

  infoMsg(IM::Med,"Starting create JT\n");

  if (numMaxCliques == 0) {
    // this shouldn't happen
    assert(0);
    infoMsg(IM::Med,"Partition is empty\n");
    return;
  } else if (numMaxCliques == 1) {
    // then nothing to do
    infoMsg(IM::Med,"Partition has only one clique\n");
    return;
  } else if (numMaxCliques == 2) {
    // then JT is easy, just connect the two cliques.
    infoMsg(IM::Med,"Partition has only two cliques\n");
    part.cliques[0].neighbors.push_back(1);
    part.cliques[1].neighbors.push_back(0);
  } else {

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
	set<RandomVariable*> sep_set;
	set_intersection(part.cliques[i].nodes.begin(),
			 part.cliques[i].nodes.end(),
			 part.cliques[j].nodes.begin(),
			 part.cliques[j].nodes.end(),
			 inserter(sep_set,sep_set.end()));
	Edge e;
	e.clique1 = i; e.clique2 = j; e.weight = sep_set.size();
	edges.push_back(e);
	infoMsg(IM::Huge,"Edge (%d,%d) has weight %d\n",i,j,e.weight);
      }
    }
  
    // sort in decreasing order by edge weight which in this
    // case is the sep-set size.
    sort(edges.begin(),edges.end(),EdgeCompare());

    unsigned joinsPlusOne = 1;
    for (unsigned i=0;i<edges.size();i++) {
      infoMsg(IM::Huge,"edge %d has weight %d\n",i,edges[i].weight);

      // TODO: optimize this to deal with ties to make message
      // passing cheaper.

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
	infoMsg(IM::Med,"Joining cliques %d and %d\n",
		edges[i].clique1,edges[i].clique2);

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
 *   the base partitions (P1,C1,C2,C3,E1,Cu0), from which all future
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

  const unsigned k = 2;

  // first create the unrolled set of random variables corresponding
  // to this JT.

  // unrolled random variables
  vector <RandomVariable*> unrolled_rvs;
  // mapping from name(frame) to integer index into unrolled_rvs.
  map < RVInfo::rvParent, unsigned > ppf;
  // number of C repetitions is M + (k+1)*S, so
  // we unroll one less than that.
  fp.unroll(gm_template.M + (k+1)*gm_template.S - 1,
	    unrolled_rvs,ppf);
  
  // copy P partition 
  new (&P1) JT_Partition(gm_template.P,unrolled_rvs,ppf,0);

  // copy C partition
  new (&C1) JT_Partition(gm_template.C,unrolled_rvs,ppf,0*gm_template.S);
  new (&Cu0) JT_Partition(gm_template.C,unrolled_rvs,ppf,0*gm_template.S);
  new (&C2) JT_Partition(gm_template.C,unrolled_rvs,ppf,1*gm_template.S);
  new (&C3) JT_Partition(gm_template.C,unrolled_rvs,ppf,2*gm_template.S);

  // copy E partition
  new (&E1) JT_Partition(gm_template.E,unrolled_rvs,ppf,2*gm_template.S);
  
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

  infoMsg(IM::Huge,"computing interface for partitions P to C");
  computePartitionInterface(P1,
			    P_ri_to_C,
			    C1,
			    C_li_to_P,
			    P_to_C_icliques_same);

  infoMsg(IM::Huge,"computing interface for partitions C to C");
  computePartitionInterface(C1,
			    C_ri_to_C,
			    C2,
			    C_li_to_C,
			    C_to_C_icliques_same);

  infoMsg(IM::Huge,"computing interface for partitions C to E");
  computePartitionInterface(C3,
			    C_ri_to_E,
			    E1,
			    E_li_to_C,
			    C_to_E_icliques_same);

  // E order, clique 0 is choosen as root arbitrarily for now.
  // TODO: see if it is possible to choose a better root for E.
  E_root_clique = 0;
  
}



/*-
 *-----------------------------------------------------------------------
 * JunctionTree::createPartitionInterface()
 *   finds the "interface" between the two partitions. The interface are
 *   the two cliques in each interface that have maximum intersection.
 *   These are the two cliques that are joined in a junction tree.
 *   It is assumed that partition 1 is just to the "left" of partition 2.
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

  infoMsg(IM::Med,"Starting create partition interface\n");
  
  unsigned max_sep_set_size = 0;
  unsigned part1_clique=0,part2_clique=0;
  unsigned part1_clique_size=0,part2_clique_size=0;
  for (unsigned i=0;i<numMaxCliques1;i++) {
    for (unsigned j=0;j<numMaxCliques2;j++) {

      set< RandomVariable*>::iterator it1;
      const set< RandomVariable*>::iterator it1_end = part1.cliques[i].nodes.end();
      set< RandomVariable*>::iterator it2;
      const set< RandomVariable*>::iterator it2_end = part2.cliques[j].nodes.end();

      unsigned sep_set_size = 0;
      for (it1 = part1.cliques[i].nodes.begin(); it1 != it1_end ; it1++) {
	for (it2 = part2.cliques[j].nodes.begin(); it2 != it2_end ; it2++) {
	  RandomVariable* rv1 = (*it1);
	  RandomVariable* rv2 = (*it2);
	  infoMsg(IM::Huge,"comparing %s(%d) with %s(%d)\n",
		  rv1->name().c_str(),rv1->frame(),		  
		  rv2->name().c_str(),rv2->frame());
	  if (rv1->frame() == rv2->frame() && rv1->name() == rv2->name()) {
	    sep_set_size++;
	  }
	}
      }
      infoMsg(IM::Huge,"clique %d of part1 and clique %d of part2 has overlap %d\n",i,j,sep_set_size);
      if (sep_set_size > max_sep_set_size) {
	part1_clique = i;
	part1_clique_size = part1.cliques[i].nodes.size();
	part2_clique = j;
	part2_clique_size = part2.cliques[i].nodes.size();
	max_sep_set_size = sep_set_size;
      }
    }
  }
  infoMsg(IM::Med,"clique %d of part1 (size %d) and clique %d of part2 (size %d) has max overlap of %d\n",
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
  createDirectedGraphOfCliques(C1,
			       C_ri_to_C);
  createDirectedGraphOfCliques(C2,
			       C_ri_to_C);
  createDirectedGraphOfCliques(C3,
			       C_ri_to_E);
  createDirectedGraphOfCliques(E1,
			       E_root_clique);
  createDirectedGraphOfCliques(Cu0,
			       C_ri_to_E);

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
JunctionTree::assignRVsToCliques()
{

  infoMsg(IM::Med,"assigning rvs to P1 partition\n");
  assignRVsToCliques("P1",P1,P_ri_to_C);

  infoMsg(IM::Med,"assigning rvs to C1 partition\n");
  C1.cliques[C_li_to_P].cumulativeAssignedNodes = 
    P1.cliques[P_ri_to_C].cumulativeAssignedNodes;
  // printf("C1's li to P cum nodes are:\n");
  // printRVSet(C1.cliques[C_li_to_P].cumulativeAssignedNodes);
  assignRVsToCliques("C1",C1,C_ri_to_C);

  infoMsg(IM::Med,"assigning rvs to C2 partition\n");
  C2.cliques[C_li_to_C].cumulativeAssignedNodes = 
    C1.cliques[C_ri_to_C].cumulativeAssignedNodes;
  assignRVsToCliques("C2",C2,C_ri_to_C);

  infoMsg(IM::Med,"assigning rvs to C3 partition\n");
  C3.cliques[C_li_to_C].cumulativeAssignedNodes = 
    C2.cliques[C_ri_to_C].cumulativeAssignedNodes;
  assignRVsToCliques("C3",C3,C_ri_to_E);

  infoMsg(IM::Med,"assigning rvs to last Cu0 partition\n");
  Cu0.cliques[C_li_to_P].cumulativeAssignedNodes =
    P1.cliques[P_ri_to_C].cumulativeAssignedNodes;
  assignRVsToCliques("Cu0",Cu0,C_ri_to_E);

  infoMsg(IM::Med,"assigning rvs to E partition\n");
  E1.cliques[E_li_to_C].cumulativeAssignedNodes =
    C3.cliques[C_ri_to_E].cumulativeAssignedNodes;
  assignRVsToCliques("E1",E1,E_root_clique);


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
 *     have been called.
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
				 const unsigned rootClique)
{
  vector<RandomVariable*> sortedNodes;
  GraphicalModel::topologicalSort(part.nodes,part.nodes,sortedNodes);

  // printf("have %d sorted nodes and %d cliques\n",sortedNodes.size(),part.cliques.size());

  for (unsigned n=0;n<sortedNodes.size();n++) {

    RandomVariable* rv = sortedNodes[n];

    // precompute the parents set for set intersection operations.
    // The one stored in the rv is a vector, but we need a set, so we
    // pre-compute it here.
    set<RandomVariable*> parSet;
    for (unsigned p=0;p<rv->allPossibleParents.size();p++) {
      parSet.insert(rv->allPossibleParents[p]);
    }
    bool assigned = false;
    multimap < vector<float>, unsigned> scoreSet;
    assignRVToClique(partName,
		     part,rootClique,
		     0,
		     rv,
		     parSet,assigned,
		     scoreSet);
    if (!assigned) {
      if (scoreSet.size() > 0) {
	// choose the first one to assign.
	unsigned clique_num = (*(scoreSet.begin())).second;
	part.cliques[clique_num].assignedNodes.insert(rv);
	part.cliques[clique_num].sortedAssignedNodes.push_back(rv);
	assigned = true;
	infoMsg(IM::Med,
		"Part %s: random variable %s(%d) assigned to clique %d\n",
		partName,
		rv->name().c_str(),rv->frame(),clique_num);
      } else {
	// rv was not assigned to this partition, it must
	// be the case that it will be assigned to a different
	// partition.
	infoMsg(IM::Med,"Part %s: random variable %s(%d) not assigned in current partition\n",partName,
		rv->name().c_str(),rv->frame());
      }
    }
    // update the cumulative RV assignments.
    if (assigned)
      getCumulativeAssignedNodes(part,rootClique);
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
			       unsigned depth,
			       RandomVariable* rv,
			       set<RandomVariable*>& parSet,
			       bool& assigned,
			       multimap < vector<float>, unsigned >& scoreSet)
{
  depth++;

  // keep a reference for easy access
  MaxClique& curClique = part.cliques[root];

  // First, check if directed_closure(rv) == (rv U parents(rv)) not in
  // the current clique, and if so then continue on down.
  bool closure_not_in_clique = (curClique.nodes.find(rv) == curClique.nodes.end());
  
  // printf("clique%d: from rv, closure_not_in_clique = %d\n",root,closure_not_in_clique);

  if (!closure_not_in_clique) {
    // need to further check that parents are in clique
    set<RandomVariable*> res;
    set_intersection(curClique.nodes.begin(),
		     curClique.nodes.end(),
		     parSet.begin(),parSet.end(),
		     inserter(res,res.end()));
    closure_not_in_clique = (res.size() != parSet.size());
  }

  // printf("clique%d: from par(rv), closure_not_in_clique = %d, num children = %d\n",root,closure_not_in_clique,part.cliques[root].children.size());

  if (!closure_not_in_clique) {
    // So closure(rv) is in current clique, meaning rv and all its
    // parents are members of this clique, but it is not nec. the case
    // that all of rv's parents are themselves assigned to this clique
    // (i.e., the parents conditional probabilities might not be
    // assigned to this clique).

    infoMsg(IM::Huge,
	    "Part %s: found random variable %s(%d) in clique %d along with all its parents\n",
	    partName,
	    rv->name().c_str(),rv->frame(),root);


    // If parents(rv) are actually assigned to this clique already,
    // then assign rv to this clique right now.
    set<RandomVariable*> res;
    set_intersection(curClique.assignedNodes.begin(),
		     curClique.assignedNodes.end(),
		     parSet.begin(),parSet.end(),
		     inserter(res,res.end()));
    const bool parents_assigned_to_clique
      = (res.size() == parSet.size());
    if (parents_assigned_to_clique) {
      // we've got all our parents assigned, so assign
      // child to this clique right now.
      curClique.assignedNodes.insert(rv);
      curClique.sortedAssignedNodes.push_back(rv);
      assigned = true;
      infoMsg(IM::Med,
	      "Part %s: random variable %s(%d) assigned to clique %d with all its parents\n",partName,
		rv->name().c_str(),rv->frame(),root);
      // all done with this routine since rv has been assigned.
      return;
    } else {

      // Note: in this discussion, we make a distinction between node
      // parents and children (i.e., original graph parents and
      // children) and JT parents and children which are cliques in
      // the JT. A child clique c of a parent p in the JT is one such
      // that c ~ p and that depth(c) = depth(p)+1, where
      // depth(root)=0 in the junction tree.

      // While we have a clique that contains this rv and all of its
      // node parents, it is not the case that all of rv's node
      // parents are assigned to this clique. Moreover, it will never
      // be the case that rv's node parents will be assigned to *this*
      // clique since we are presumably assigning rv nodes to cliques
      // in topological order, and at this point, we've already
      // assigned rv's node parents somewhere else.

      // Therefore, we have to do the best we can (i.e., we're now in
      // a situation where rv will live in a clique without all of its
      // node parents assigned).

      // What we do, is continue on down. It may be the the case that
      // we find a clique with rv and all of rv's node parents
      // assigned (in which case we assign rv right then), but in the
      // mean time we keep track of this clique's "score" using a set
      // of heuristics. At some time later, the rv will be assigned to
      // the clique with the *LOWEST* score.

      // Heuristics: want to assign rv to a clique where it:
      // 
      //    A) has all of its node parents in clique (already satisfied here)
      //
      //    B) has its node parents actually assigned (not satisfied here)
      //
      //    C) be in a clique CL close enough to the root clique so
      //    that all rv's node parents are assigned to cliques that
      //    are JT decendants of CL in the JT rooted at root. This way
      //    the separator driven iterations will preserve any
      //    sparsity driven zeros in the CPTs.
      //
      //    D) has lots of node children whose other node parents are
      //    already (or will be) assigned.
      // 
      //    E) is assign as far away from root as possible, to try
      //    encourage it to have as much "influence" as possible in
      //    JT parents (but this is only a poor heuristic)
      //
      //    F) (NOT DONE BELOW) is assigned to a clique where it has the greatest
      //    number of node descendants (but this is only a heuristic
      //    since to gain benefit, we would need to have it be such
      //    that those node descendants have their node parents in
      //    clique as well.
      // 

      // TODO: The variable being in a clique with its parents
      // *assigned* is really only important for sparse or
      // deterministic nodes. Ideally, this bit of info should also
      // affect the assignment process.

      // TODO: think this through as this code is only
      // heuristic. There is most likely some theoretically best thing
      // to do here.
      
      vector<float> score;

      // Push back items in decreasing order of priority.  Lower
      // numbers is better (e.g., more negative or less positive is
      // higher priority).  First thing inserted has highest priority.

      // First (top priority), 
      // insert value:
      //   0, if all parents have been assigned in the cummulative set, or
      //   1, if not.
      {
	set<RandomVariable*> res;
	set_intersection(curClique.cumulativeAssignedNodes.begin(),
			 curClique.cumulativeAssignedNodes.end(),
			 parSet.begin(),parSet.end(),
			 inserter(res,res.end()));
	bool parents_not_assigned = (res.size() != parSet.size());
	score.push_back(parents_not_assigned);
	if (!parents_not_assigned) {
	  // then this is good! we've found at least one.
	  infoMsg(IM::Med,
		  "Part %s: found random variable %s(%d) in clique %d with all parents in children cliques\n",partName,
		  rv->name().c_str(),rv->frame(),root);
	}
      }

      // Next, add number of children in current clique. If rv
      // has lots of children in this clique, it is hopeful that
      // other parents of those children might also be assigned
      // to the same clique.
      unsigned numChildren = 0;
      for (unsigned i=0;i<rv->allPossibleChildren.size();i++) {
	RandomVariable* child = rv->allPossibleChildren[i];
	if (curClique.nodes.find(child) != curClique.nodes.end())
	  numChildren++;
      }
      score.push_back(-numChildren);

      // Next, distance from root, among the higher priorities that
      // are equal, try to be as far away from the root as possible so
      // as to prune away as much zero as possible as early as
      // possible.
      score.push_back(-depth);

      // And so on. We can push back as many heuristics as we want.
      // alternatively, perhaps take a weighted average of some of
      // them??

      // TODO: add more heuristics here, and/or produce better prioritized
      // order above.

      // ...
      
      // done inserting heuristicss, now insert the score and the
      // current clique into the set.
      pair < vector<float>, unsigned> p(score,root);
      scoreSet.insert(p);
    }
  }

  // continue on down.
  for (unsigned childNo=0;
       childNo<curClique.children.size();childNo++) {
    const unsigned child = part.cliques[root].children[childNo];
    assignRVToClique(partName,part,child,depth,rv,parSet,assigned,scoreSet);
    if (assigned)
      break;
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

  set<RandomVariable*> res;
  const set<RandomVariable*> empty;
  for (unsigned childNo=0;
       childNo<curClique.children.size();childNo++) {

    const unsigned child = part.cliques[root].children[childNo];

    getCumulativeAssignedNodes(part,child);

    // Note: this will do a bunch of redundant work (i.e., inserting
    // elements into sets where the elements are already contained in
    // the sets), but it doesn't need to run that fast, so the
    // redundant work is not a big deal.
    set_union(empty.begin(),empty.end(),
	      part.cliques[child].cumulativeAssignedNodes.begin(),
	      part.cliques[child].cumulativeAssignedNodes.end(),
	      inserter(res,res.end()));
  }
  set_union(curClique.assignedNodes.begin(),
	    curClique.assignedNodes.end(),
	    res.begin(),res.end(),
	    inserter(curClique.cumulativeAssignedNodes,
		     curClique.cumulativeAssignedNodes.end()));

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

  setUpMessagePassingOrder(C1,
			   C_ri_to_C,
			   C1_message_order,
			   C_li_to_P,
			   C1_leaf_cliques);

  setUpMessagePassingOrder(C3,
			   C_ri_to_E,
			   C3_message_order,
			   C_li_to_C,
			   C3_leaf_cliques);

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
  order.clear();
  // a tree of N nodes has N-1 edges.
  order.reserve(part.cliques.size() - 1);
  setUpMessagePassingOrderRecurse(part,root,order,excludeFromLeafCliques,leaf_cliques);
}
/*-
 *-----------------------------------------------------------------------
 * JunctionTree::setUpMessagePassingOrderRecurse()
 *   Recursive support routine for setUpMessagePassingOrder() and
 *   should only be called from there. See that
 *   routine for documentation.
d * 
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

  createSeparators(P1,
		   P1_message_order);
  P1.cliques[P_ri_to_C].ceSendSeparator = ~0x0; // set to invalid value
  // P1 has no LI, so nothing more to do for P1 here.

  createSeparators(C1,
		   C1_message_order);
  // Create separator of interface cliques
  C1.separators.push_back(SeparatorClique(P1.cliques[P_ri_to_C],C1.cliques[C_li_to_P]));
  // update right partitions LI clique to include new separator
  C1.cliques[C_li_to_P].ceReceiveSeparators.push_back(C1.separators.size()-1);
  // don't update left partitions RI clique's send separator since handeled explicitly
  C1.cliques[C_ri_to_C].ceSendSeparator = ~0x0; //set to invalid value

  createSeparators(C2,
		   C1_message_order);
  // Create separator of interface cliques
  C2.separators.push_back(SeparatorClique(C1.cliques[C_ri_to_C],C2.cliques[C_li_to_C]));
  // update right partitions LI clique to include new separator
  C2.cliques[C_li_to_C].ceReceiveSeparators.push_back(C2.separators.size()-1);
  // don't update left partitions RI clique's send separator since handeled explicitly
  C2.cliques[C_ri_to_C].ceSendSeparator = ~0x0; //set to invalid value

  createSeparators(C3,
		   C3_message_order);
  // Create separator of interface cliques
  C3.separators.push_back(SeparatorClique(C2.cliques[C_ri_to_C],C3.cliques[C_li_to_C]));
  // update right partitions LI clique to include new separator
  C3.cliques[C_li_to_C].ceReceiveSeparators.push_back(C3.separators.size()-1);
  // don't update left partitions RI clique's send separator since handeled explicitly
  C3.cliques[C_ri_to_E].ceSendSeparator = ~0x0; //set to invalid value


  createSeparators(E1,
		   E1_message_order);
  // Create separator of interface cliques
  E1.separators.push_back(SeparatorClique(C3.cliques[C_ri_to_E],E1.cliques[E_li_to_C]));
  // update right partitions LI clique to include new separator
  E1.cliques[E_li_to_C].ceReceiveSeparators.push_back(E1.separators.size()-1);
  // don't update left partitions RI clique's send separator since handeled explicitly
  C3.cliques[E_root_clique].ceSendSeparator = ~0x0; //set to invalid value


  createSeparators(Cu0,
		   C3_message_order);
  // Create separator of interface cliques
  Cu0.separators.push_back(SeparatorClique(P1.cliques[P_ri_to_C],Cu0.cliques[C_li_to_P]));
  // update right partitions LI clique to include new separator
  Cu0.cliques[C_li_to_P].ceReceiveSeparators.push_back(Cu0.separators.size()-1);
  // don't update left partitions RI clique's send separator since handeled explicitly
  Cu0.cliques[C_ri_to_E].ceSendSeparator = ~0x0; //set to invalid value

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
  computeSeparatorIterationOrders(C1);
  computeSeparatorIterationOrders(C2);
  computeSeparatorIterationOrders(C3);
  computeSeparatorIterationOrders(E1);
  computeSeparatorIterationOrders(Cu0);
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

  // partial and then ultimately the final union of of all nodes in
  // separators for incomming messages for this clique.
  set<RandomVariable*> accumSeps;


  if (numSeparators == 0) {
    // This must be a leaf-node clique relatve to root.
    // 'accumSeps' is already empty so no need to do anything there.
  } else if (numSeparators == 1) {
    accumSeps = part.separators[clique.ceReceiveSeparators[0]].nodes;
    part.separators[clique.ceReceiveSeparators[0]].accumulatedIntersection.clear();
    part.separators[clique.ceReceiveSeparators[0]].remainder =
      part.separators[clique.ceReceiveSeparators[0]].nodes;
  } else if (numSeparators == 2) {
    

    // shortcuts to separator 0 and 1
    SeparatorClique& s0 = part.separators[clique.ceReceiveSeparators[0]];
    SeparatorClique& s1 = part.separators[clique.ceReceiveSeparators[1]];

    // intersection of the two separators
    set<RandomVariable*> sepIntersection;
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
	      inserter(accumSeps,accumSeps.end()));

    assert ( s0.accumulatedIntersection.size() + s0.remainder.size() == s0.nodes.size() );
    assert ( s1.accumulatedIntersection.size() + s1.remainder.size() == s1.nodes.size() );

  } else {
    // there are 3 or more separators, determine proper order and then
    // compute running accumulated intersection relative to that order.

    // TODO: need to determine if there is an optimum order or not.
    // note: code to compute 'an' order is commented out and taged
    // with ABCDEFGHIJK at end of this file.

    // Compute the cummulative intersection of the sepsets
    // using the current sepset order.

    {

      // initialize union of all previous separators
      accumSeps = 
	part.separators[clique.ceReceiveSeparators[0]].nodes;

      part.separators[clique.ceReceiveSeparators[0]].accumulatedIntersection.clear();
      for (unsigned sep=1;sep<numSeparators;sep++) {
      
	// reference variables for easy access
	set<RandomVariable*>& sepNodes
	  = part.separators[clique.ceReceiveSeparators[sep]].nodes;
	set<RandomVariable*>& sepAccumInter
	  = part.separators[clique.ceReceiveSeparators[sep]].accumulatedIntersection;


	// create the intersection of 1) the union of all previous nodes in
	// the sep order, and 2) the current sep nodes.
	sepAccumInter.clear();
	set_intersection(accumSeps.begin(),accumSeps.end(),
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
	set<RandomVariable*> res;
	set_union(sepNodes.begin(),sepNodes.end(),
		  accumSeps.begin(),accumSeps.end(),
		  inserter(res,res.end()));
	accumSeps = res;	
      }
    }

  }

  // lastly, assign unassignedIteratedNodes in this clique
  {
    // compute: unassignedIteratedNodes  = nodes - (assignedNodes U accumSeps)
    // first: compute res = nodes - assignedNodes
    set<RandomVariable*> res;
    set_difference(clique.nodes.begin(),clique.nodes.end(),
		   clique.assignedNodes.begin(),clique.assignedNodes.end(),
		   inserter(res,res.end()));
    // next: compute unassignedIteratedNodes = res - accumSeps
    // note at this point accumSeps contains the union of all nodes in all separators
    set_difference(res.begin(),res.end(),
		   accumSeps.begin(),accumSeps.end(),
		   inserter(clique.unassignedIteratedNodes,
			    clique.unassignedIteratedNodes.end()));
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
  prepareForUnrolling(P1);
  prepareForUnrolling(C1);
  prepareForUnrolling(C2);
  prepareForUnrolling(C3);
  prepareForUnrolling(E1);
  prepareForUnrolling(Cu0);
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
void
JunctionTree::unroll(const unsigned int k)
{

  // first create the unrolled set of random variables corresponding
  // to this JT.

  unsigned basicTemplateUnrollAmount;
  unsigned modifiedTemplateUnrollAmount;
  unsigned numUsableFrames;
  unsigned frameStart;
  if (!gm_template.computeUnrollParamaters(k,
					   basicTemplateUnrollAmount,
					   modifiedTemplateUnrollAmount,
					   numUsableFrames,
					   frameStart))
    error("Can't unroll\n"); // TODO: fix this error.

  fprintf(stderr,"numFrames = %d, unrolling BT %d times, MT %d times\n",
	  k,
	  basicTemplateUnrollAmount,
	  modifiedTemplateUnrollAmount);

  // unrolled random variables
  vector <RandomVariable*> unrolled_rvs;
  // mapping from 'name+frame' to integer index into unrolled_rvs.
  map < RVInfo::rvParent, unsigned > ppf;

  fp.unroll(basicTemplateUnrollAmount,unrolled_rvs,ppf);


  // TODO: clear out the old and pre-allocate for new size.

  // preallocate
  jtIPartitions.resize(modifiedTemplateUnrollAmount+3);

  // printf("unroll: doing P\n");
  // copy P partition into partitions[0]
  new (&jtIPartitions[0]) JT_InferencePartition(P1,unrolled_rvs,ppf,0);


  if (modifiedTemplateUnrollAmount == 0) {
    new (&jtIPartitions[1]) JT_InferencePartition(Cu0,unrolled_rvs,ppf,0*gm_template.S);
    new (&jtIPartitions[2]) JT_InferencePartition(E1,unrolled_rvs,ppf,-2*gm_template.S);
  } else if (modifiedTemplateUnrollAmount == 1) {
    new (&jtIPartitions[1]) JT_InferencePartition(C1,unrolled_rvs,ppf, 0*gm_template.S);
    new (&jtIPartitions[2]) JT_InferencePartition(C3,unrolled_rvs,ppf,-1*gm_template.S);
    new (&jtIPartitions[3]) JT_InferencePartition(E1,unrolled_rvs,ppf,-1*gm_template.S);
  } else {
    new (&jtIPartitions[1]) JT_InferencePartition(C1,unrolled_rvs,ppf, 0*gm_template.S);

    for (unsigned i = 1; i < modifiedTemplateUnrollAmount; i++) {
      new (&jtIPartitions[i+1]) JT_InferencePartition(C2,unrolled_rvs,ppf, (i-1)*gm_template.S);
    }

    new (&jtIPartitions[modifiedTemplateUnrollAmount+1]) 
      JT_InferencePartition(C3,unrolled_rvs,ppf, 
			    ((int)modifiedTemplateUnrollAmount-2)*gm_template.S);
    new (&jtIPartitions[modifiedTemplateUnrollAmount+2]) 
      JT_InferencePartition(E1,unrolled_rvs,ppf, 
			    ((int)modifiedTemplateUnrollAmount-2)*gm_template.S);

  }
}




/*-
 *-----------------------------------------------------------------------
 * JunctionTree::collectEvidence()
 *
 * Preconditions:
 *
 * Postconditions:
 *
 * Side Effects:

 *
 * Results:

 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::collectEvidence()
{

  // this routine handles all of:
  // unrolled 0 times: (so there is a single P1,Cu0, and E1)  
  // unrolled 1 time: so there is a P1, C1, C3, E1
  // unrolled 2 or more times: so there is a P1 C1 [C2 ...] C3, E1


  // do P1 messages
  for (unsigned msgNo=0;msgNo < P1_message_order.size(); msgNo ++) {
    const unsigned from = P1_message_order[msgNo].first;
    const unsigned to = P1_message_order[msgNo].second;
    infoMsg(IM::Med,"message in P1,part[0] from clique %d --> clique %d\n",
	    from,to);
    jtIPartitions[0].maxCliques[from].
      collectEvidenceFromSeparators(jtIPartitions[0]);
    jtIPartitions[0].maxCliques[from].
      ceCollectToSeparator(jtIPartitions[0]);
  }
  // do P1-C1 interface 
  jtIPartitions[0].maxCliques[P_ri_to_C].
    collectEvidenceFromSeparators(jtIPartitions[0]);
  jtIPartitions[0].maxCliques[P_ri_to_C].
    ceCollectToSeparator(jtIPartitions[0],
			 jtIPartitions[1].
			 separatorCliques[jtIPartitions[1].
					  separatorCliques.size()-1]);

  // then do C1 [C2 ...] messages
  unsigned partNo;
  for (partNo = 1; partNo < (jtIPartitions.size()-2); partNo ++ ) {
    for (unsigned msgNo=0;msgNo < C1_message_order.size(); msgNo ++) {
      const unsigned from = C1_message_order[msgNo].first;
      const unsigned to = C1_message_order[msgNo].second;
      infoMsg(IM::Med,"message in C1,part[%d] from clique %d --> clique %d\n",
	      partNo,from,to);
      jtIPartitions[partNo].maxCliques[from].
	collectEvidenceFromSeparators(jtIPartitions[partNo]);
      jtIPartitions[partNo].maxCliques[from].
	ceCollectToSeparator(jtIPartitions[partNo]);
    }

    // do C1-nextC interface
    jtIPartitions[partNo].maxCliques[C_ri_to_C].
      collectEvidenceFromSeparators(jtIPartitions[partNo]);
    jtIPartitions[partNo].maxCliques[C_ri_to_C].
      ceCollectToSeparator(jtIPartitions[partNo],
			   jtIPartitions[partNo+1].
			   separatorCliques[jtIPartitions[partNo+1].
					    separatorCliques.size()-1]);
  }

  // then do C3 messages
  for (unsigned msgNo=0;msgNo < C3_message_order.size(); msgNo ++) {
    const unsigned from = C3_message_order[msgNo].first;
    const unsigned to = C3_message_order[msgNo].second;
    infoMsg(IM::Med,"message in C3,part[%d] from clique %d --> clique %d\n",
	    partNo,from,to);
    jtIPartitions[partNo].maxCliques[from].
      collectEvidenceFromSeparators(jtIPartitions[partNo]);
    jtIPartitions[partNo].maxCliques[from].
      ceCollectToSeparator(jtIPartitions[partNo]);
  }
  // do C3-E1 interface
  jtIPartitions[partNo].maxCliques[C_ri_to_E].
    collectEvidenceFromSeparators(jtIPartitions[partNo]);
  jtIPartitions[partNo].maxCliques[C_ri_to_E].
    ceCollectToSeparator(jtIPartitions[partNo],
			 jtIPartitions[partNo+1].
			 separatorCliques[jtIPartitions[partNo+1].
					  separatorCliques.size()-1]);



  // finally E1 messages
  partNo++;
  for (unsigned msgNo=0;msgNo < E1_message_order.size(); msgNo ++) {
    const unsigned from = E1_message_order[msgNo].first;
    const unsigned to = E1_message_order[msgNo].second;
    infoMsg(IM::Med,"message in E1,part[%d] from clique %d --> clique %d\n",
	    partNo,from,to);
    jtIPartitions[partNo].maxCliques[from].
      collectEvidenceFromSeparators(jtIPartitions[partNo]);
    jtIPartitions[partNo].maxCliques[from].
      ceCollectToSeparator(jtIPartitions[partNo]);
  }
  jtIPartitions[partNo].maxCliques[E_root_clique].
    collectEvidenceFromSeparators(jtIPartitions[partNo]);

}







/*-
 *-----------------------------------------------------------------------
 * JunctionTree::probEvidence()
 *
 * Preconditions:
 *    collectEvidence must have been called.
 * Postconditions:
 *
 * Side Effects:

 *
 * Results:

 *
 *-----------------------------------------------------------------------
 */
logpr
JunctionTree::probEvidence()
{
  return jtIPartitions[jtIPartitions.size()-1].maxCliques[E_root_clique].
    sumProbabilities();
}



#if 0


  // code from above that is commented out but not yet deleted. 
  // it is kept here for now to clarify development.
  // tag: ABCDEFGHIJK at end of this file.

    {
    // TODO: need to determine if there is an optimum order or not. If
    // not, no need to sort here. In general, the number of "live"
    // entries at the end of a sep traversal will be the same no
    // matter the order that is chosen. What we want here, however,
    // is something that kills non-live entries ASAP, to avoid extra
    // work. The question is if this is depending on properties
    // of just the separators independent of the particular set of
    // sparse entries these separators currently contain.

    // current code is set up to do a max spanning tree on sep sets
    // but it is commented out for now.
    
    // 3 or more separators, need to compute max spanning tree

    vector < set<unsigned> >  findSet;
    findSet.resize(numMaxCliques);
    for (unsigned i=0;i<numSeparators;i++) {
      set<unsigned> iset;
      iset.insert(i);
      findSet[i] = iset;
    }

    vector< Edge > edges;
    edges.reserve((numSeparators*(numSeparators-1))/2);

    set<RandomVariable*> sep_intr_set;
    for (unsigned i=0;i<numSeparators;i++) {
      for (unsigned j=i+1;j<numSeparators;j++) {

	const set<RandomVariable*>& sep_i =
	  part.separators[clique.ceReceiveSeparators[i]];
	const set<RandomVariable*>& sep_j =
	  part.separators[clique.ceReceiveSeparators[j]];
	sep_intr_set.clear();
	set_intersection(sep_i.begin(),sep_i.end(),
			 sep_j.begin(),sep_j.end(),
			 inserter(sep_intr_set,sep_intr_set.end()));
	Edge e;
	e.clique1 = i; e.clique2 = j; 
	e.weight = sep_intr_set.size();
	edges.push_back(e);
	infoMsg(IM::Huge,"Edge (%d,%d) has weight %d\n",i,j,e.weight);
      }
    }
  
    // sort in decreasing order by edge weight which in this
    // case is the sep-set size.
    sort(edges.begin(),edges.end(),EdgeCompare());

    unsigned joinsPlusOne = 1;
    for (unsigned i=0;i<edges.size();i++) {
      infoMsg(IM::Huge,"edge %d has weight %d\n",i,edges[i].weight);

      // TODO: optimize this to deal with ties to make message
      // passing cheaper.

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
	infoMsg(IM::Med,"Joining cliques %d and %d\n",
		edges[i].clique1,edges[i].clique2);

	// make edge between edges[i].clique1 and edges[i].clique2

	// part.cliques[edges[i].clique1].neighbors.push_back(edges[i].clique2);
	// part.cliques[edges[i].clique2].neighbors.push_back(edges[i].clique1);

	if (++joinsPlusOne == numSeparators)
	  break;
      }
    }
    
    // since these are separators, and all for a given clique A
    // with some other clique B_i, and since all separators are
    // subsets of A, but it is not the case that all separators have
    // an intersection with each other. Goal is to find an order to
    // traverse the sep sets so that we maximize the amount of overlap
    // from the so-far-traversed to the to-be-traversed.

    }

#endif

