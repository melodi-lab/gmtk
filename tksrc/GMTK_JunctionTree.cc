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
//        Inference Partition support
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Constructors/Destructors 
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

InferencePartition::InferencePartition(
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

    // TODO: ultimately turn this just into an assert
    if ( ppf.find(rvp) == ppf.end() ) {
      warning("ERROR: can't find rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
      assert ( ppf.find(rvp) != ppf.end() );
    }

    RandomVariable* nrv = newRvs[ppf[rvp]];
    nodes.insert(nrv);
  }
  cliques.reserve(from_part.cliques.size());
  // 
  // NOTE: It is Crucial for the cliques in the cloned partition to be
  // inserted in the *SAME ORDER* as in the partition being cloned.
  for (unsigned i=0;i<from_part.cliques.size();i++) {
    cliques.push_back(MaxClique(from_part.cliques[i],
				newRvs,ppf,frameDelta));
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
    // then nothing to do
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
    // need to do some work.

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
 *   Modifies member variables, as above.
 *
 * Results:
 *   Modifies member variables as above.
 *
 *-----------------------------------------------------------------------
 */
void
JunctionTree::computePartitionInterfaces()
{
  // unroll internal data structures one time, setting
  // up the base_partition structures.
  base_unroll(1);

  // now we have that
  //  base_partition[0] is a P partition
  //  base_partition[1] is a C partition (C1)
  //  base_partition[2] is a C partition (C2)
  //  base_partition[3] is an E 
  // Use these new partitions to find the various interface
  // cliques.

  infoMsg(IM::Huge,"computing interface for partitions P to C");
  computePartitionInterface(base_partitions[0],
			    P_ri_to_C,
			    base_partitions[1],
			    C_li_to_P,
			    P_to_C_icliques_same);

  infoMsg(IM::Huge,"computing interface for partitions C to C");
  computePartitionInterface(base_partitions[1],
			    C_ri_to_C,
			    base_partitions[2],
			    C_li_to_C,
			    C_to_C_icliques_same);

  infoMsg(IM::Huge,"computing interface for partitions C to E");
  computePartitionInterface(base_partitions[2],
			    C_ri_to_E,
			    base_partitions[3],
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
 *   Modifies all neighbors variables within the cliques within the
 *   partition.
 *
 * Results:
 *   The left partitions right interface clique and the right partitions
 *   left interface clique are returned via the part1_ric and part2_lic variables.
 *
 *-----------------------------------------------------------------------
 */

void
JunctionTree::computePartitionInterface(InferencePartition& part1,
					// partitition 1's right interface clique
					unsigned int& part1_ric,
					InferencePartition& part2,
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
 * JunctionTree::unroll()
 *   sets up the internal data structures and creates a junction tree corresponding
 *   to the graph with C' unrolled k times, k >= 0. (i.e., we have
 *   a graph that has (k+1) copies of C'). It does this using the initial
 *   partitions given in gm_template (P', C', and E').
 *
 *   This routine should really only be called once to set up
 *   the base_partitions data structure, from which all future
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
JunctionTree::base_unroll(const unsigned int k)
{

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

  // clear out the old and pre-allocate for new size.
  base_partitions.clear();
  base_partitions.reserve(k+3);
  
  // copy P partition into base_partitions[0]
  base_partitions.push_back(InferencePartition(gm_template.P,unrolled_rvs,ppf,0));

  for (unsigned i = 0; i <= k; i++) {
    base_partitions.push_back(InferencePartition(gm_template.C,unrolled_rvs,ppf,i*gm_template.S));
  }

  // copy E partition into base_partitions[k+1]
  base_partitions.push_back(InferencePartition(gm_template.E,unrolled_rvs,ppf,k*gm_template.S));
  
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
  assert (base_partitions.size() == 4);

  // do P to C
  createDirectedGraphOfCliques(base_partitions[0],
			       P_ri_to_C);

  // do C to C
  createDirectedGraphOfCliques(base_partitions[1],
			       C_ri_to_C);

  // do C to E
  createDirectedGraphOfCliques(base_partitions[2],
			       C_ri_to_E);

  // do E
  createDirectedGraphOfCliques(base_partitions[3],
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
JunctionTree::createDirectedGraphOfCliques(InferencePartition& part,
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
JunctionTree::createDirectedGraphOfCliquesRecurse(InferencePartition& part,
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
  assert (base_partitions.size() == 4);

  // infoMsg(IM::Med,"assigning rvs to orig P partition\n");
  // assignRVsToCliques(gm_template.P,P_ri_to_C);

  infoMsg(IM::Med,"assigning rvs to P partition\n");
  assignRVsToCliques(0,base_partitions[0],P_ri_to_C);


  infoMsg(IM::Med,"assigning rvs to C partition\n");
  assignRVsToCliques(1,base_partitions[1],C_ri_to_C);


  infoMsg(IM::Med,"assigning rvs to last C partition\n");
  assignRVsToCliques(2,base_partitions[2],C_ri_to_E);


  infoMsg(IM::Med,"assigning rvs to E partition\n");
  assignRVsToCliques(3,base_partitions[3],E_root_clique);



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
JunctionTree::assignRVsToCliques(const unsigned partNo,
				 InferencePartition& part,
				 const unsigned rootClique)
{
  vector<RandomVariable*> sortedNodes;
  GraphicalModel::topologicalSort(part.nodes,part.nodes,sortedNodes);

  // printf("have %d sorted nodes and %d cliques\n",sortedNodes.size(),part.cliques.size());

  for (unsigned n=0;n<sortedNodes.size();n++) {

    RandomVariable* rv = sortedNodes[n];

    // precompute the parents set for set intersection operations.
    set<RandomVariable*> parSet;
    for (unsigned p=0;p<rv->allPossibleParents.size();p++) {
      parSet.insert(rv->allPossibleParents[p]);
    }
    bool assigned = false;
    multimap < vector<float>, unsigned> scoreSet;
    assignRVToClique(partNo,
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
		"Part %d: random variable %s(%d) assigned to clique %d\n",
		partNo,
		rv->name().c_str(),rv->frame(),clique_num);
      } else {
	// rv was not assigned to this partition, it must
	// be the case that it will be assigned to a different
	// partition.
	infoMsg(IM::Med,"Part %d: random variable %s(%d) not assigned in current partition\n",partNo,
		rv->name().c_str(),rv->frame());
      }
    }
    // re-compute the cumulative RV assignments.
    if (assigned)
      getCumulativeAssignedNodes(part,
				 rootClique,
				 part.cliques[rootClique].cumulativeAssignedNodes);
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
JunctionTree::assignRVToClique(const unsigned partNo,
			       InferencePartition& part,
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
	    "Part %d: found random variable %s(%d) in clique %d along with all its parents\n",
	    partNo,
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
	      "Part %d: random variable %s(%d) assigned to clique %d with all its parents\n",partNo,
		rv->name().c_str(),rv->frame(),root);
      // all done with this routine since rv has been assigned.
      return;
    } else {

      // While we have a clique that contains this rv and all of its
      // parents, it is not the case that all of rv's parents are
      // assigned to this clique. Moreover, it will never be the case
      // that rv's parents will be assigned to *this* clique since we
      // are presumably assigning rv's to cliques in topological
      // order, and at this point, we've already assigned rv's parents
      // somewhere else.

      // Therefore, we have to do the best we can (i.e., we're now in
      // a situation where rv will live in a clique without all of its
      // parents assigned).

      // What we do, is continue on down. It may be the the case that
      // we find a clique with rv and all of rv's parents assigned (in
      // which case we assign rv right then), but in the mean time we
      // keep track of this clique's "score" using a set of
      // heuristics. At some time later, the variable will be assigned
      // to the clique with the *LOWEST* score.

      // Possible Heuristics: want to assign rv to a clique where it 
      //    A) has all of its parents  in clique (already satisfied here)
      //    B) has its parents actually assigned (not satisfied here)

      //    C) be in a clique close enough to the root (say at depth
      //       d) so that all rv's parents are assigned to cliques that
      //       JT decendants in the JT rooted at root.

      //    D) lots of children whose other parents are already
      //       (or will be) assigned.
      //    E) assign as far away from root as possible.

      // Note. The variable being in a clique with its parents
      // *assigned* is really only important for sparse or deterministic
      // nodes. Ideally, this bit of info should also affect the
      // assignment process.
      
      vector<float> score;

      // Push back items in decreasing order of priority.
      // Lower is better (e.g., more negative or less positive is better).
      // First thing inserted has highest priority.


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
		  "Part %d: found random variable %s(%d) in clique %d with all parents in children cliques\n",partNo,
		  rv->name().c_str(),rv->frame(),root);
	}
      }

      // Next, distance from root, among the higher priorities that
      // are equal, try to be as far away from the root as possible so
      // as to prune away as much zero as possible as early as
      // possible.
      score.push_back(-depth);

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

      // And so on. We can push back as many heuristics as we want.
      // alternatively, perhaps take a weighted average of some of
      // them??

      // insert the score and the current clique into the set.
      pair < vector<float>, unsigned> p(score,root);
      scoreSet.insert(p);
    }
  }

  // continue on down.
  for (unsigned childNo=0;
       childNo<curClique.children.size();childNo++) {
    const unsigned child = part.cliques[root].children[childNo];
    assignRVToClique(partNo,part,child,depth,rv,parSet,assigned,scoreSet);
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
JunctionTree::getCumulativeAssignedNodes(InferencePartition& part,
					 const unsigned root,
					 set<RandomVariable*> &res)
{
  MaxClique& curClique = part.cliques[root];
  // res.clear(); @@@ no need to have, do this and do cumulative clique from before.
  for (unsigned childNo=0;
       childNo<curClique.children.size();childNo++) {

    const unsigned child = part.cliques[root].children[childNo];
    getCumulativeAssignedNodes(part,child,part.cliques[child].cumulativeAssignedNodes);
    set_union(part.cliques[child].assignedNodes.begin(),
	      part.cliques[child].assignedNodes.end(),
	      part.cliques[child].cumulativeAssignedNodes.begin(),
	      part.cliques[child].cumulativeAssignedNodes.end(),
	      inserter(res,res.end()));

  }

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
  // do P to C order
  setUpMessagePassingOrder(base_partitions[0],
			   P_ri_to_C,
			   P_to_C_message_order);
  // do C to C order
  setUpMessagePassingOrder(base_partitions[1],
			   C_ri_to_C,
			   C_to_C_message_order);
  // do C to E order
  setUpMessagePassingOrder(base_partitions[2],
			   C_ri_to_E,
			   C_to_E_message_order);

  // do E order
  setUpMessagePassingOrder(base_partitions[3],
			   E_root_clique,
			   E_message_order);

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
JunctionTree::setUpMessagePassingOrder(InferencePartition& part,
				       const unsigned root,
				       vector< pair<unsigned,unsigned> >& order)
{
  order.clear();
  // a tree of N nodes has N-1 edges.
  order.reserve(part.cliques.size() - 1);
  setUpMessagePassingOrderRecurse(part,root,order);
}
/*-
 *-----------------------------------------------------------------------
 * JunctionTree::setUpMessagePassingOrderRecurse()
 *   Recursive support routine for setUpMessagePassingOrder() and
 *   should only be called from there. See that
 *   routine for documentation.
 * 
 *-----------------------------------------------------------------------
 */					
void
JunctionTree::setUpMessagePassingOrderRecurse(InferencePartition& part,
					     const unsigned root,
					     vector< pair<unsigned,unsigned> >& order)
{

  for (unsigned childNo=0;
       childNo<part.cliques[root].children.size();
       childNo++) {
    const unsigned child = part.cliques[root].children[childNo];
    setUpMessagePassingOrderRecurse(part,child,order);
    pair<unsigned,unsigned> msg;
    msg.first = child;
    msg.second = root;
    order.push_back(msg);
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
  // do P to C 
  createSeparators(base_partitions[0],
		   P_to_C_message_order);

  // do C to C
  createSeparators(base_partitions[1],
		   C_to_C_message_order);

  // do C to E
  createSeparators(base_partitions[2],C_to_E_message_order);

  // do E
  createSeparators(base_partitions[3],E_message_order);

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
JunctionTree::createSeparators(InferencePartition& part,
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
 * JunctionTree::unroll()
 *   sets up the internal data structures and creates a junction tree corresponding
 *   to the graph with C' unrolled k times, k >= 0. (i.e., we have
 *   a graph that has (k+1) copies of C'). 
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

  assert ( base_partitions.size() == 4 );

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

  // printf("unrolled variables %d times\n",gm_template.M + (k+1)*gm_template.S - 1);

  // clear out the old and pre-allocate for new size.
  partitions.clear();
  partitions.reserve(k+3);

  // printf("unroll: doing P\n");
  // copy P partition into partitions[0]
  partitions.push_back(InferencePartition(base_partitions[0],unrolled_rvs,ppf,0));

  // copy the C partitions into partitions[1:k-1]
  // this does nothing if k == 0.
  for (unsigned i = 0; i < k; i++) {
    // printf("unroll: doing C[%d]\n",i);
    partitions.push_back(InferencePartition(base_partitions[1],unrolled_rvs,ppf,i*gm_template.S));

  }

  // copy the last C partition into partitions[k]
  // printf("unroll: doing C[%d]\n",k);
  partitions.push_back(InferencePartition(base_partitions[2],unrolled_rvs,ppf,((int)k-1)*gm_template.S));

  // printf("unroll: doing E\n");
  // copy E partition into partitions[k+1]
  partitions.push_back(InferencePartition(base_partitions[3],unrolled_rvs,ppf,((int)k-1)*gm_template.S));
  
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

  // first do P messages
  for (unsigned msgNo=0;msgNo < P_to_C_message_order.size(); msgNo ++){
    infoMsg(IM::Med,"message in P from clique %d --> clique %d\n",
	   P_to_C_message_order[msgNo].first,
	   P_to_C_message_order[msgNo].second);
  }

  infoMsg(IM::Med,"message from P clique %d --> C clique %d\n",
	 P_ri_to_C,C_li_to_P);

  // then do C messsages
  for (unsigned partNo = 1; partNo <= (partitions.size()-2); partNo ++ ) {
    
    for (unsigned msgNo=0;msgNo < C_to_C_message_order.size(); msgNo ++){
      infoMsg(IM::Med,"message in C%d from clique %d --> clique %d\n",
	     partNo,
	     C_to_C_message_order[msgNo].first,
	     C_to_C_message_order[msgNo].second);
    }

    if (partNo < (partitions.size()-2)) {
      infoMsg(IM::Med,"message from C%d clique %d --> C%d clique %d\n",
	     partNo,C_ri_to_C,
	     partNo+1,C_li_to_C);
    }
  }

  infoMsg(IM::Med,"message from C clique %d --> E clique %d\n",
	 C_ri_to_E,E_li_to_C);

  // and lastly do E messages
  for (unsigned msgNo=0;msgNo < E_message_order.size(); msgNo ++) {
    infoMsg(IM::Med,"message in E from clique %d --> clique %d\n",
	   E_message_order[msgNo].first,
	   E_message_order[msgNo].second);
  }
}


