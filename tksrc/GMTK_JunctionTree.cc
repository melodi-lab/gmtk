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
		       const unsigned int frameDelta,
		       // the left and right interface variables for
		       // this JT partition Empty if doesn't exist
		       // (say for an P or E partition). These have
		       // their own frame deltas since they might be
		       // different.
		       // Left interface:
		       const set <RandomVariable*>& from_liVars,
		       const unsigned int liFrameDelta,
		       // Right interface:
		       const set <RandomVariable*>& from_riVars,
		       const unsigned int riFrameDelta,
		       // Information todo the mapping.
		       vector <RandomVariable*>& newRvs,
		       map < RVInfo::rvParent, unsigned >& ppf)

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
      error("INTERNAL ERROR: can't find rv %s(%d+%d)=%s(%d) in unrolled RV set, from from_part\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
    }

    RandomVariable* nrv = newRvs[ppf[rvp]];
    nodes.insert(nrv);
  }
  liNodes.clear();
  for (it = from_liVars.begin();
       it != from_liVars.end();
       it++) {
    RandomVariable* rv = (*it);
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+liFrameDelta;    

    if ( ppf.find(rvp) == ppf.end() ) {
      error("INTERNAL ERROR: can't find rv %s(%d+%d)=%s(%d) in unrolled RV set, from from_liVars\n",
	    rv->name().c_str(),rv->frame(),liFrameDelta,
	    rvp.first.c_str(),rvp.second);
    }

    RandomVariable* nrv = newRvs[ppf[rvp]];
    liNodes.insert(nrv);
  }
  riNodes.clear();
  for (it = from_riVars.begin();
       it != from_riVars.end();
       it++) {
    RandomVariable* rv = (*it);
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+riFrameDelta;    

    if ( ppf.find(rvp) == ppf.end() ) {
      error("INTERNAL ERROR: can't find rv %s(%d+%d)=%s(%d) in unrolled RV set, from from_riVars\n",
	    rv->name().c_str(),rv->frame(),riFrameDelta,
	    rvp.first.c_str(),rvp.second);
    }

    RandomVariable* nrv = newRvs[ppf[rvp]];
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


// bear-bones constructor to quickly compute
// JT weight.
JT_Partition::JT_Partition(
 Partition& from_part,
 // Left interface:
 const set <RandomVariable*>& from_liVars,
 // Right interface:
 const set <RandomVariable*>& from_riVars
 )
{
  nodes = from_part.nodes;
  liNodes = from_liVars;
  riNodes = from_riVars;
  // make the cliques.
  cliques = from_part.cliques;
}



void
JT_Partition::findInterfaceCliques(const set <RandomVariable*>& iNodes,
				   unsigned& iClique,
				   bool& iCliqueSameAsInterface)
{
  if (iNodes.size() > 0) {
    iCliqueSameAsInterface = false;
    // starting invalid clique.
    iClique = ~0x0u; 
    // the weight of the interface clique, to break ties.
    float interfaceCliqueWeight=HUGE;

    for (unsigned cliqueNo=0;cliqueNo<cliques.size();cliqueNo++) {
      // check that clique fully covers iNodes 
      set<RandomVariable*> res;
      set_intersection(cliques[cliqueNo].nodes.begin(),
		       cliques[cliqueNo].nodes.end(),
		       iNodes.begin(),
		       iNodes.end(),
		       inserter(res,res.end()));
      
      if (res.size() == iNodes.size()) {
	// we've found a candidate.
	if (iClique == ~0x0u) {
	  // first time
	  iClique = cliqueNo;
	  interfaceCliqueWeight = MaxClique::computeWeight(cliques[cliqueNo].nodes);
	  if (res.size() == cliques[cliqueNo].nodes.size()) {
	    iCliqueSameAsInterface = true;
	  } else {
	    iCliqueSameAsInterface = false;
	  }
	} else {
	  float new_weight = MaxClique::computeWeight(cliques[cliqueNo].nodes);
	  if (new_weight < interfaceCliqueWeight) {
	    iClique = cliqueNo;
	    interfaceCliqueWeight = new_weight;
	    if (res.size() == cliques[cliqueNo].nodes.size()) {
	      iCliqueSameAsInterface = true;
	    } else {
	      iCliqueSameAsInterface = false;
	    }
	  }
	}
      }
    }
  } else {
    iClique = ~0x0;
    iCliqueSameAsInterface = false;
  }
}


void
JT_Partition::findLInterfaceClique(unsigned& liClique,bool& liCliqueSameAsInterface)
{
  findInterfaceCliques(liNodes,liClique,liCliqueSameAsInterface);
}


void
JT_Partition::findRInterfaceClique(unsigned& riClique,bool& riCliqueSameAsInterface)
{
  findInterfaceCliques(riNodes,riClique,riCliqueSameAsInterface);
}

unsigned
JT_Partition::cliqueWithMaxWeight()
{
  float weight = -1.0;
  unsigned res= ~0x0;
  assert ( cliques.size() > 0 );
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
  float weight = HUGE;
  unsigned res = ~0x0;
  assert ( cliques.size() > 0 );
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
					     vector <RandomVariable*>& newRvs,
					     map < RVInfo::rvParent, unsigned >& ppf,
					     const unsigned int frameDelta)
  : origin(from_part)
{
  unsigned i;
  // first allocate space with empty (and unusable) entries
  maxCliques.resize(origin.cliques.size());
  separatorCliques.resize(origin.separators.size());

  // then actually re-construct the objects in the array appropriately.
  for ( i=0;i<maxCliques.size();i++) {
    new (&maxCliques[i]) InferenceMaxClique(origin.cliques[i],
					    newRvs,ppf,frameDelta);
  }
  for ( i=0;i<separatorCliques.size();i++) {
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
	// define the edge
	e.clique1 = i; e.clique2 = j; 
	// first push sep set size. To get a JT, we must
	// always choose from among the cliques that
	// have the largest intersection size.
	e.weights.push_back((float)sep_set.size());

	// for ties, we next push back negative weight of separator
	e.weights.push_back(-(float)MaxClique::computeWeight(sep_set));

	// printf("weight of clique %d = %f, %d = %f\n",
	// i,part.cliques[i].weight(),
	// j,part.cliques[j].weight());

	// if ties still, we next push back negative weight of two
	// cliques together.
	set<RandomVariable*> clique_union;
	set_union(part.cliques[i].nodes.begin(),
		  part.cliques[i].nodes.end(),
		  part.cliques[j].nodes.begin(),
		  part.cliques[j].nodes.end(),
		  inserter(clique_union,clique_union.end()));
	e.weights.push_back(-(float)MaxClique::computeWeight(clique_union));

	// add the edge.
	edges.push_back(e);
	infoMsg(IM::Huge,"Edge (%d,%d) has sep size %.0f, log(sep state) %f, log(union state) %f \n",
		i,j,
		e.weights[0],
		-e.weights[1],
		-e.weights[2]);
      }
    }
  
    // sort in decreasing order by edge weight which in this
    // case is the sep-set size.
    sort(edges.begin(),edges.end(),EdgeCompare());

    unsigned joinsPlusOne = 1;
    for (unsigned i=0;i<edges.size();i++) {
	infoMsg(IM::Huge,"Edge %d has sep size %.0f, log(sep state) %f, log(union state) %f \n",
		i,
		edges[i].weights[0],
		-edges[i].weights[1],
		-edges[i].weights[2]);

      // TODO: optimize this to deal with ties to make message
      // passing cheaper.
      //    Ideas: among all ties
      //        a) maximize number of variables in same frame (or near each other)
      //        b) minimize cardinality or state space of separator
      //        c) minimize number of neighbors in each clique (i.e., 
      //           if cliques already have neighbors, choose the ones with fewer.
      //        d) integrate with RV value assignment to minimize
      //           the number of unassigned clique nodes (since they're
      //           iterated over w/o knowledge of any parents. If this
      //           ends up being a search, make this be offline, in with gmtkTriangulate


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
	infoMsg(IM::Med,"Joining cliques %d and %d (edge %d) with intersection size %.0f\n",
		edges[i].clique1,edges[i].clique2,i,edges[i].weights[0]);

	if (edges[i].weights[0] == 0.0) {
	  warning("ERROR: junction tree creation trying to join two cliques (%d and %d) with size 0 set intersection. Possible non-triangulated graph.",
		edges[i].clique1,edges[i].clique2);
	  // TODO: print out two cliques that are trying to be joined.
	  fprintf(stderr,"Clique %d: ",edges[i].clique1);
	  part.cliques[edges[i].clique1].printCliqueNodes(stderr);
	  fprintf(stderr,"Clique %d: ",edges[i].clique2);
	  part.cliques[edges[i].clique2].printCliqueNodes(stderr);
	  error("exiting");
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

  const unsigned k = 1;

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
  

  set <RandomVariable*> empty;
  // copy P partition 
  new (&P1) JT_Partition(gm_template.P,0*gm_template.S,
			 empty,0*gm_template.S,
			 gm_template.PCInterface_in_P,0*gm_template.S,
			 unrolled_rvs,ppf);

  if (gm_template.leftInterface) { 
    // left interface case.
    // Template has a (P)(C)(CE) = (P')(C')(E')
    // Unroll to a P Cu0 Co CE
    new (&Cu0) JT_Partition(gm_template.C,0*gm_template.S,
			    gm_template.PCInterface_in_C,0*gm_template.S,
			    gm_template.CEInterface_in_C,0*gm_template.S,
			    unrolled_rvs,ppf);
    new (&Co) JT_Partition(gm_template.C,1*gm_template.S,
			   gm_template.CEInterface_in_C,0*gm_template.S,
			   gm_template.CEInterface_in_C,1*gm_template.S,			   			   unrolled_rvs,ppf);
  } else { // right interface
    // Right interface case.
    // Template has a (PC)(C)(E) = (P')(C')(E')
    // Unroll to a P Co Cu0 CE

    new (&Co) JT_Partition(gm_template.C,0*gm_template.S,
			   gm_template.PCInterface_in_C,0*gm_template.S,			   
			   gm_template.PCInterface_in_C,1*gm_template.S,
			   unrolled_rvs,ppf);

    new (&Cu0) JT_Partition(gm_template.C,1*gm_template.S,
			    gm_template.PCInterface_in_C,1*gm_template.S,
			    gm_template.CEInterface_in_C,1*gm_template.S,
			    unrolled_rvs,ppf);

  }

  // copy E partition
  new (&E1) JT_Partition(gm_template.E,1*gm_template.S,
			 gm_template.CEInterface_in_E,1*gm_template.S,
			 empty,0*gm_template.S,
			 unrolled_rvs,ppf);
  
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
  P1.findRInterfaceClique(P_ri_to_C,P_riCliqueSameAsInterface);
  E1.findLInterfaceClique(E_li_to_C,E_liCliqueSameAsInterface);

  if (gm_template.leftInterface) {
    bool Cu0_liCliqueSameAsInterface;
    bool Cu0_riCliqueSameAsInterface;
    bool Co_liCliqueSameAsInterface;
    bool Co_riCliqueSameAsInterface;

    Cu0.findLInterfaceClique(C_li_to_P,Cu0_liCliqueSameAsInterface);
    Cu0.findRInterfaceClique(C_ri_to_C,Cu0_riCliqueSameAsInterface);
    Co.findLInterfaceClique(C_li_to_C,Co_liCliqueSameAsInterface);
    Co.findRInterfaceClique(C_ri_to_E,Co_riCliqueSameAsInterface);
    
    // sanity check
    assert (C_ri_to_C == C_ri_to_E);

    P_to_C_icliques_same =
      P_riCliqueSameAsInterface && Cu0_liCliqueSameAsInterface;

    C_to_C_icliques_same =
      Cu0_riCliqueSameAsInterface && Co_liCliqueSameAsInterface;    

    C_to_E_icliques_same =
      Co_riCliqueSameAsInterface && E_liCliqueSameAsInterface;

  } else { // right interface

    bool Cu0_liCliqueSameAsInterface;
    bool Cu0_riCliqueSameAsInterface;
    bool Co_liCliqueSameAsInterface;
    bool Co_riCliqueSameAsInterface;

    Co.findLInterfaceClique(C_li_to_P,Co_liCliqueSameAsInterface);
    Co.findRInterfaceClique(C_ri_to_C,Co_riCliqueSameAsInterface);
    Cu0.findLInterfaceClique(C_li_to_C,Cu0_liCliqueSameAsInterface);
    Cu0.findRInterfaceClique(C_ri_to_E,Cu0_riCliqueSameAsInterface);

    // sanity check
    assert (C_li_to_C == C_li_to_P);
    
    P_to_C_icliques_same =
      P_riCliqueSameAsInterface && Co_liCliqueSameAsInterface;

    C_to_C_icliques_same =
      Co_riCliqueSameAsInterface && Cu0_liCliqueSameAsInterface;    

    C_to_E_icliques_same =
      Cu0_riCliqueSameAsInterface && E_liCliqueSameAsInterface;

  }

  // E order, clique 0 is choosen as root arbitrarily for now.
  // TODO: see if it is possible to choose a better root for E.
  // Perhaps use the maximum weight clique???
  E_root_clique = 0;
  // E_root_clique = E1.cliqueWithMinWeight();
  
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
	part2_clique_size = part2.cliques[j].nodes.size();
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
  if (gm_template.leftInterface) {
    createDirectedGraphOfCliques(Cu0,
				 C_ri_to_C);
    createDirectedGraphOfCliques(Co,
				 C_ri_to_E);
  } else { // right interface
    createDirectedGraphOfCliques(Co,
				 C_ri_to_C);
    createDirectedGraphOfCliques(Cu0,
				 C_ri_to_E);
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

  if (gm_template.leftInterface) {
    infoMsg(IM::Med,"assigning rvs to Cu0 partition\n");
    Cu0.cliques[C_li_to_P].cumulativeAssignedNodes = 
      P1.cliques[P_ri_to_C].cumulativeAssignedNodes;
    assignRVsToCliques("Cu0",Cu0,C_ri_to_C);

    infoMsg(IM::Med,"assigning rvs to Co partition\n");
    Co.cliques[C_li_to_C].cumulativeAssignedNodes = 
      Cu0.cliques[C_ri_to_C].cumulativeAssignedNodes;
    assignRVsToCliques("Co",Co,C_ri_to_C);

    infoMsg(IM::Med,"assigning rvs to E partition\n");
    E1.cliques[E_li_to_C].cumulativeAssignedNodes =
      Co.cliques[C_ri_to_E].cumulativeAssignedNodes;
    assignRVsToCliques("E1",E1,E_root_clique);

  } else { // right interface 
    infoMsg(IM::Med,"assigning rvs to Co partition\n");
    Co.cliques[C_li_to_P].cumulativeAssignedNodes = 
      P1.cliques[P_ri_to_C].cumulativeAssignedNodes;
    assignRVsToCliques("Co",Co,C_ri_to_C);

    infoMsg(IM::Med,"assigning rvs to Cu0 partition\n");
    Cu0.cliques[C_li_to_C].cumulativeAssignedNodes = 
      Co.cliques[C_ri_to_C].cumulativeAssignedNodes;
    assignRVsToCliques("Cu0",Cu0,C_ri_to_E);

    infoMsg(IM::Med,"assigning rvs to E partition\n");
    E1.cliques[E_li_to_C].cumulativeAssignedNodes =
      Cu0.cliques[C_ri_to_E].cumulativeAssignedNodes;
    assignRVsToCliques("E1",E1,E_root_clique);
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
      // fprintf(stderr,"about to insert rv with address %X\n",(void*)rv->allPossibleParents[p]);
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
	// rv was not assigned to this partition, it must be the case
	// that it will be assigned to a different partition. This
	// could come from the partition before or the partition after
	// the current partition. If there are only forward-directed
	// arrows, it will come from the previous partition (and vice
	// versa if there are only backwards going edges). If we have
	// both directed edges, it could be assigned in either left or
	// right partition.
	infoMsg(IM::Med,"Part %s: random variable %s(%d) not assigned in current partition\n",partName,
		rv->name().c_str(),rv->frame());
	// in any event, keep track of this node.
	part.unassignedInPartition.insert(rv);
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

      // Previous Parents in Junction Tree.
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


      // Distance from root, among the higher priorities that
      // are equal, try to be as far away from the root as possible so
      // as to prune away as much zero as possible as early as
      // possible.
      score.push_back(-depth);


      // Number of children in current clique. If rv has lots of
      // children in this clique, it is hopeful that other parents of
      // those children might also be assigned to the same clique.
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

    const unsigned child = curClique.children[childNo];

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

  if (gm_template.leftInterface) {
    setUpMessagePassingOrder(Cu0,
			     C_ri_to_C,
			     Cu0_message_order,
			     C_li_to_P,
			     Cu0_leaf_cliques);

    setUpMessagePassingOrder(Co,
			     C_ri_to_C,
			     Co_message_order,
			     C_li_to_C,
			     Co_leaf_cliques);

  } else { // right interface
    setUpMessagePassingOrder(Co,
			     C_ri_to_C,
			     Co_message_order,
			     C_li_to_P,
			     Co_leaf_cliques);
    setUpMessagePassingOrder(Cu0,
			     C_ri_to_E,
			     Cu0_message_order,
			     C_li_to_C,
			     Cu0_leaf_cliques);
  }

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

  createSeparators(Cu0,
		   Cu0_message_order);
  createSeparators(Co,
		   Co_message_order);
  createSeparators(E1,
		   E1_message_order);

  if (gm_template.leftInterface) {

    // Create separator of interface cliques
    Cu0.separators.push_back(SeparatorClique(P1.cliques[P_ri_to_C],Cu0.cliques[C_li_to_P]));
    // update right partitions LI clique to include new separator
    Cu0.cliques[C_li_to_P].ceReceiveSeparators.push_back(Cu0.separators.size()-1);
    // don't update left partitions RI clique's send separator since handeled explicitly
    Cu0.cliques[C_ri_to_C].ceSendSeparator = ~0x0; //set to invalid value

    // Create separator of interface cliques
    Co.separators.push_back(SeparatorClique(Cu0.cliques[C_ri_to_C],Co.cliques[C_li_to_C]));
    // update right partitions LI clique to include new separator
    Co.cliques[C_li_to_C].ceReceiveSeparators.push_back(Co.separators.size()-1);
    // don't update left partitions RI clique's send separator since handeled explicitly
    Co.cliques[C_ri_to_E].ceSendSeparator = ~0x0; //set to invalid value


    // Create separator of interface cliques
    E1.separators.push_back(SeparatorClique(Co.cliques[C_ri_to_E],E1.cliques[E_li_to_C]));
    // update right partitions LI clique to include new separator
    E1.cliques[E_li_to_C].ceReceiveSeparators.push_back(E1.separators.size()-1);
    // don't update left partitions RI clique's send separator since handeled explicitly
    E1.cliques[E_root_clique].ceSendSeparator = ~0x0; //set to invalid value


  } else { // right interface 

    // Create separator of interface cliques
    Co.separators.push_back(SeparatorClique(P1.cliques[P_ri_to_C],Co.cliques[C_li_to_P]));
    // update right partitions LI clique to include new separator
    Co.cliques[C_li_to_P].ceReceiveSeparators.push_back(Co.separators.size()-1);
    // don't update left partitions RI clique's send separator since handeled explicitly
    Co.cliques[C_ri_to_C].ceSendSeparator = ~0x0; //set to invalid value

    // Create separator of interface cliques
    Cu0.separators.push_back(SeparatorClique(Co.cliques[C_ri_to_C],Cu0.cliques[C_li_to_C]));
    // update right partitions LI clique to include new separator
    Cu0.cliques[C_li_to_C].ceReceiveSeparators.push_back(Cu0.separators.size()-1);
    // don't update left partitions RI clique's send separator since handeled explicitly
    Cu0.cliques[C_ri_to_E].ceSendSeparator = ~0x0; //set to invalid value

    // Create separator of interface cliques
    E1.separators.push_back(SeparatorClique(Cu0.cliques[C_ri_to_E],E1.cliques[E_li_to_C]));
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
  computeSeparatorIterationOrders(Cu0);
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
  clique.accumSeps.clear();

  if (numSeparators == 0) {
    // This must be a leaf-node clique relatve to root.
    // 'accumSeps' is already empty so no need to do anything there.
  } else if (numSeparators == 1) {
    // shortcut to separator 0
    SeparatorClique& s0 = part.separators[clique.ceReceiveSeparators[0]];
    clique.accumSeps = s0.nodes;
    s0.accumulatedIntersection.clear();
    s0.remainder = s0.nodes;
    assert ( s0.accumulatedIntersection.size() + s0.remainder.size() == s0.nodes.size() );
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
	      inserter(clique.accumSeps,clique.accumSeps.end()));

    assert ( s0.accumulatedIntersection.size() + s0.remainder.size() == s0.nodes.size() );
    assert ( s1.accumulatedIntersection.size() + s1.remainder.size() == s1.nodes.size() );

  } else {
    // there are 3 or more separators, determine proper order and then
    // compute running accumulated intersection relative to that order.

    // TODO: need to determine if there is an optimum order or not.
    // note: code to compute 'an' order is commented out and taged
    // with ABCDEFGHIJK at end of this file.
    // 
    // In otherwords, rerder clique.ceReceiveSeparators for maximal overlap.
    // 

    // Compute the cummulative intersection of the sepsets
    // using the current sepset order.

    {

      // initialize union of all previous separators

      clique.accumSeps = 
	part.separators[clique.ceReceiveSeparators[0]].nodes;
      part.separators[clique.ceReceiveSeparators[0]].accumulatedIntersection.clear();
      part.separators[clique.ceReceiveSeparators[0]].remainder
	= part.separators[clique.ceReceiveSeparators[0]].nodes;

      for (unsigned sep=1;sep<numSeparators;sep++) {
      
	// reference variables for easy access
	set<RandomVariable*>& sepNodes
	  = part.separators[clique.ceReceiveSeparators[sep]].nodes;
	set<RandomVariable*>& sepAccumInter
	  = part.separators[clique.ceReceiveSeparators[sep]].accumulatedIntersection;


	// create the intersection of 1) the union of all previous nodes in
	// the sep order, and 2) the current sep nodes.
	sepAccumInter.clear();
	set_intersection(clique.accumSeps.begin(),clique.accumSeps.end(),
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
		  clique.accumSeps.begin(),clique.accumSeps.end(),
		  inserter(res,res.end()));
	clique.accumSeps = res;	
      }
    }

  }

  // lastly, assign unassignedIteratedNodes in this clique
  {
    // unassigned nodes, unassigned iterated nodes, 
    // compute: unassignedIteratedNodes  = nodes - (assignedNodes U accumSeps)
    // first: compute res = nodes - assignedNodes
    set<RandomVariable*> res;
    set_difference(clique.nodes.begin(),clique.nodes.end(),
		   clique.assignedNodes.begin(),clique.assignedNodes.end(),
		   inserter(res,res.end()));
    // next: compute unassignedIteratedNodes = res - accumSeps
    // note at this point accumSeps contains the union of all nodes in all separators
    set_difference(res.begin(),res.end(),
		   clique.accumSeps.begin(),clique.accumSeps.end(),
		   inserter(clique.unassignedIteratedNodes,
			    clique.unassignedIteratedNodes.end()));
  }

}




/*-
 *-----------------------------------------------------------------------
 * JunctionTree::getPrecedingIteratedUnassignedNodes()
 *
 *   Computes for each clique the union of the set of nodes which have
 *   been iterated unassigned in previous cliques in the JT. This is
 *   used to determine which, in each clique, nodes should be iterated
 *   over and which are already assigned by a separator driven
 *   iteration.
 *
 * Preconditions:
 *   computeSeparatorIterationOrder() must have been called.
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
JunctionTree::getPrecedingIteratedUnassignedNodes(JT_Partition& part,
						  const unsigned root)
{

  MaxClique& curClique = part.cliques[root];
  set<RandomVariable*>& res = curClique.precedingUnassignedIteratedNodes;
  for (unsigned childNo=0;
       childNo<curClique.children.size();childNo++) {

    const unsigned child = curClique.children[childNo];
    getPrecedingIteratedUnassignedNodes(part,child);

    set_union(part.cliques[child].precedingUnassignedIteratedNodes.begin(),
	      part.cliques[child].precedingUnassignedIteratedNodes.end(),
	      part.cliques[child].unassignedIteratedNodes.begin(),
	      part.cliques[child].unassignedIteratedNodes.end(),
	      inserter(res,res.end()));
  }
  curClique.computeAssignedNodesToIterate();
}
void
JunctionTree::getPrecedingIteratedUnassignedNodes()
{
  set<RandomVariable*> res;

  // TODO: this should be a member function in P1.
  getPrecedingIteratedUnassignedNodes(P1,P_ri_to_C);
  res.clear();
  set_union(P1.cliques[P_ri_to_C].unassignedIteratedNodes.begin(),
	    P1.cliques[P_ri_to_C].unassignedIteratedNodes.end(),
	    P1.cliques[P_ri_to_C].precedingUnassignedIteratedNodes.begin(),
	    P1.cliques[P_ri_to_C].precedingUnassignedIteratedNodes.end(),	    
	    inserter(res,res.end()));

  if (gm_template.leftInterface) {
    
    Cu0.cliques[C_li_to_P].precedingUnassignedIteratedNodes = res;
    getPrecedingIteratedUnassignedNodes(Cu0,C_ri_to_C);

    res.clear();
    set_union(Cu0.cliques[C_ri_to_C].unassignedIteratedNodes.begin(),
	      Cu0.cliques[C_ri_to_C].unassignedIteratedNodes.begin(),
	      Cu0.cliques[C_ri_to_C].precedingUnassignedIteratedNodes.begin(),
	      Cu0.cliques[C_ri_to_C].precedingUnassignedIteratedNodes.begin(),	    
	      inserter(res,res.end()));
    Co.cliques[C_li_to_C].precedingUnassignedIteratedNodes = res;
    getPrecedingIteratedUnassignedNodes(Co,C_ri_to_E);

    res.clear();
    set_union(Co.cliques[C_ri_to_E].unassignedIteratedNodes.begin(),
	      Co.cliques[C_ri_to_E].unassignedIteratedNodes.end(),
	      Co.cliques[C_ri_to_E].precedingUnassignedIteratedNodes.begin(),
	      Co.cliques[C_ri_to_E].precedingUnassignedIteratedNodes.end(),
	      inserter(res,res.end()));
    E1.cliques[E_li_to_C].precedingUnassignedIteratedNodes = res;
    getPrecedingIteratedUnassignedNodes(E1,E_root_clique);

  } else { // right interface 

    Co.cliques[C_li_to_P].precedingUnassignedIteratedNodes = res;
    getPrecedingIteratedUnassignedNodes(Co,C_ri_to_C);

    res.clear();
    set_union(Co.cliques[C_ri_to_C].unassignedIteratedNodes.begin(),
	      Co.cliques[C_ri_to_C].unassignedIteratedNodes.begin(),
	      Co.cliques[C_ri_to_C].precedingUnassignedIteratedNodes.begin(),
	      Co.cliques[C_ri_to_C].precedingUnassignedIteratedNodes.begin(),	    
	      inserter(res,res.end()));
    Cu0.cliques[C_li_to_C].precedingUnassignedIteratedNodes = res;
    getPrecedingIteratedUnassignedNodes(Cu0,C_ri_to_E);

    res.clear();
    set_union(Cu0.cliques[C_ri_to_E].unassignedIteratedNodes.begin(),
	      Cu0.cliques[C_ri_to_E].unassignedIteratedNodes.end(),
	      Cu0.cliques[C_ri_to_E].precedingUnassignedIteratedNodes.begin(),
	      Cu0.cliques[C_ri_to_E].precedingUnassignedIteratedNodes.end(),
	      inserter(res,res.end()));
    E1.cliques[E_li_to_C].precedingUnassignedIteratedNodes = res;
    getPrecedingIteratedUnassignedNodes(E1,E_root_clique);

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
				 const unsigned rootClique)
{
  MaxClique& curClique = part.cliques[rootClique];

  set <RandomVariable*> empty;
  // @@ add in unassigned in partition information to next call
  double weight = curClique.weightInJunctionTree();
  // double weight = curClique.weight();
  for (unsigned childNo=0;
       childNo<curClique.children.size();childNo++) {
    double child_weight = junctionTreeWeight(part,
					     curClique.children[childNo]);
    // i.e., log addition for weight = weight + child_weight
    weight = weight + log10(1+pow(10,child_weight - weight));
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
				 const set<RandomVariable*>& interfaceNodes)
{
  Partition part;
  const set <RandomVariable*> emptySet;
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
    jt_part.findRInterfaceClique(root,tmp);
  else {
    // Presumably, this is an E partition, so the root should be done
    // same as E_root_clique computed above.
    root = 0; 
  }
    
  createDirectedGraphOfCliques(jt_part,root);
  assignRVsToCliques("candidate partition",jt_part,root);

  vector< pair<unsigned,unsigned> > message_order;
  vector< unsigned > leaf_cliques;
  setUpMessagePassingOrder(jt_part,
			   root,
			   message_order,
			   ~0x0,
			   leaf_cliques);
  createSeparators(jt_part,message_order);
  computeSeparatorIterationOrders(jt_part);
  getPrecedingIteratedUnassignedNodes(jt_part,root);
  
  // return jt_part.cliques[root].precedingUnassignedIteratedNodes.size();
  double weight = junctionTreeWeight(jt_part,root);

#if 0
  unsigned badness_count=0;
  set <RandomVariable*>::iterator it;
  for (it = jt_part.cliques[root].precedingUnassignedIteratedNodes.begin();
       it != jt_part.cliques[root].precedingUnassignedIteratedNodes.end();
       it++) 
    {
      RandomVariable* rv = (*it);
      assert ( rv-> discrete );
      DiscreteRandomVariable* drv = (DiscreteRandomVariable*)rv;

      if (drv->sparse())
	badness_count ++;
    }
#endif

  return weight;

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
	  junctionTreeWeight(P1,P_ri_to_C));
  printAllJTInfo(f,P1,P_ri_to_C);
  fprintf(f,"\n\n");

  if (gm_template.leftInterface) {

    fprintf(f,"===============================\n");
    fprintf(f,"   Cu0 partition information: JT_weight = %f\n",
	    junctionTreeWeight(Cu0,C_ri_to_C));
    printAllJTInfo(f,Cu0,C_ri_to_C);
    fprintf(f,"\n\n");


    fprintf(f,"===============================\n");
    fprintf(f,"   Co partition information: JT_weight = %f\n",
	    junctionTreeWeight(Co,C_ri_to_E));
    printAllJTInfo(f,Co,C_ri_to_E);
    fprintf(f,"\n\n");

  } else { // right interface

    fprintf(f,"===============================\n");
    fprintf(f,"   Co partition information: JT_weight = %f\n",
	    junctionTreeWeight(Co,C_ri_to_C));
    printAllJTInfo(f,Co,C_ri_to_C);
    fprintf(f,"\n\n");

    fprintf(f,"===============================\n");
    fprintf(f,"   Cu0 partition information: JT_weight = %f\n",
	    junctionTreeWeight(Cu0,C_ri_to_E));
    printAllJTInfo(f,Cu0,C_ri_to_E);
    fprintf(f,"\n\n");

  }

  fprintf(f,"===============================\n");
  fprintf(f,"   E1 partition information: JT_weight = %f\n",
	  junctionTreeWeight(E1,E_root_clique));
  printAllJTInfo(f,E1,E_root_clique);
  fprintf(f,"\n\n");

  // print message order information
  fprintf(f,"===============================\n\n");    

  fprintf(f,"===============================\n");  
  fprintf(f,"   P1 message order\n");
  printMessageOrder(f,P1_message_order);
  fprintf(f,"\n\n");

  if (gm_template.leftInterface) {
    fprintf(f,"===============================\n");  
    fprintf(f,"   Cu0 message order\n");
    printMessageOrder(f,Cu0_message_order);
    fprintf(f,"\n\n");

    fprintf(f,"===============================\n");  
    fprintf(f,"   Co message order\n");
    printMessageOrder(f,Co_message_order);
    fprintf(f,"\n\n");
  } else {

    fprintf(f,"===============================\n");  
    fprintf(f,"   Co message order\n");
    printMessageOrder(f,Co_message_order);
    fprintf(f,"\n\n");

    fprintf(f,"===============================\n");  
    fprintf(f,"   Cu0 message order\n");
    printMessageOrder(f,Cu0_message_order);
    fprintf(f,"\n\n");
  }

  fprintf(f,"===============================\n");  
  fprintf(f,"   E1 message order\n");
  printMessageOrder(f,E1_message_order);
  fprintf(f,"\n\n");


  fclose(f);
}

void
JunctionTree::printAllJTInfo(FILE* f,
			     JT_Partition& part,
			     const unsigned root)
{
  // print cliques information
  fprintf(f,"=== Clique Information ===\n");
  fprintf(f,"Number of cliques = %d\n",part.cliques.size());
  printAllJTInfoCliques(f,part,root,0);

  // print separator information
  fprintf(f,"\n=== Separator Information ===\n");
  fprintf(f,"Number of separators = %d\n",part.separators.size());
  for (unsigned sepNo=0;sepNo<part.separators.size();sepNo++) {
    fprintf(f,"== Separator number: %d\n",sepNo);
    part.separators[sepNo].printAllJTInfo(f);
  }

}

void
JunctionTree::printAllJTInfoCliques(FILE* f,
				    JT_Partition& part,
				    const unsigned root,
				    const unsigned treeLevel)
{
  // print cliques information
  for (unsigned i=0;i<treeLevel;i++) fprintf(f,"  ");
  fprintf(f,"== Clique number: %d\n",root);
  part.cliques[root].printAllJTInfo(f,treeLevel);
  for (unsigned childNo=0;
       childNo<part.cliques[root].children.size();childNo++) {
    unsigned child = part.cliques[root].children[childNo];
    printAllJTInfoCliques(f,part,child,treeLevel+1);
  }
}

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
  prepareForUnrolling(P1);
  prepareForUnrolling(Cu0);
  prepareForUnrolling(Co);
  prepareForUnrolling(E1);
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
  unsigned modifiedTemplateUnrollAmount;
  unsigned numUsableFrames;
  unsigned frameStart;
  if (!gm_template.computeUnrollParamaters(numFrames,
					   basicTemplateUnrollAmount,
					   modifiedTemplateUnrollAmount,
					   numUsableFrames,
					   frameStart))
    error("Can't unroll\n"); // TODO: fix this error.
    // return 0

  infoMsg(IM::Default,"numFrames = %d, numUsableFrames = %d\n",numFrames,numUsableFrames);
  infoMsg(IM::Default,"numFrames = %d, unrolling BT %d times, MT %d times\n",
	  numFrames,
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

  unsigned partNo = 0;
  const unsigned numCoPartitions = modifiedTemplateUnrollAmount;

  new (&jtIPartitions[partNo++]) JT_InferencePartition(P1,unrolled_rvs,ppf,0*gm_template.S);
  if (gm_template.leftInterface) 
    new (&jtIPartitions[partNo++]) JT_InferencePartition(Cu0,unrolled_rvs,ppf,0*gm_template.S);
  for (unsigned p=0;p<numCoPartitions;p++) 
    new (&jtIPartitions[partNo++]) JT_InferencePartition(Co,unrolled_rvs,ppf,p*gm_template.S);
  if (!gm_template.leftInterface) 
    new (&jtIPartitions[partNo++]) 
      JT_InferencePartition(Cu0,unrolled_rvs,ppf,
			    ((int)modifiedTemplateUnrollAmount-1)*gm_template.S);
  new (&jtIPartitions[partNo++]) 
    JT_InferencePartition(E1,unrolled_rvs,ppf,
			  ((int)modifiedTemplateUnrollAmount-1)*gm_template.S);

  assert (partNo == jtIPartitions.size());

  return numUsableFrames;
}


void
JunctionTree::ceGatherIntoRoot(JT_InferencePartition& part,
			       const unsigned root,
			       vector< pair<unsigned,unsigned> >& message_order,
			       const char*const part_type_name,
			       const unsigned part_num)
{
  // do partition messages
  for (unsigned msgNo=0;msgNo < message_order.size(); msgNo ++) {
    const unsigned from = message_order[msgNo].first;
    const unsigned to = message_order[msgNo].second;
    part.maxCliques[from].
      ceGatherFromIncommingSeparators(part);
    infoMsg(IM::Mod,
	    "CE: message %s,part[%d]: clique %d --> clique %d\n",
	    part_type_name,part_num,from,to);
    part.maxCliques[from].
      ceSendToOutgoingSeparator(part);
  }
  // collect to partition's root clique
  part.maxCliques[root].
    ceGatherFromIncommingSeparators(part);
}

void
JunctionTree::ceSendToNextPartition(JT_InferencePartition& previous_part,
				    const unsigned previous_part_root,
				    const char*const previous_part_type_name,
				    const unsigned previous_part_num,
				    JT_InferencePartition& next_part,
				    const unsigned next_part_leaf,
				    const char*const next_part_type_name,
				    const unsigned next_part_num)
{
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

  const unsigned numCoPartitions = jtIPartitions.size()-3;
  unsigned partNo = 0;
  // set up appropriate name for debugging output.
  const char* prv_nm;

  prv_nm = "P1";
  ceGatherIntoRoot(jtIPartitions[partNo],
		   P_ri_to_C,
		   P1_message_order,
		   prv_nm,partNo);


  if (gm_template.leftInterface) {
    ceSendToNextPartition(jtIPartitions[partNo],P_ri_to_C,"P1",partNo,
			  jtIPartitions[partNo+1],C_li_to_P,"Cu0",partNo+1);
    partNo++;
    prv_nm = "Cu0";
    ceGatherIntoRoot(jtIPartitions[partNo],
		     C_ri_to_C,
		     Cu0_message_order,
		     prv_nm,partNo);
  }

  for (unsigned p = 0; p < numCoPartitions; p++ ) {
    ceSendToNextPartition(jtIPartitions[partNo],C_ri_to_C,prv_nm,partNo,
			  jtIPartitions[partNo+1],C_li_to_C,"Co",partNo+1);
    partNo++;
    prv_nm = "Co";
    ceGatherIntoRoot(jtIPartitions[partNo],
		     C_ri_to_C,
		     Co_message_order,
		     prv_nm,partNo);
  }

  if (!gm_template.leftInterface) {
    ceSendToNextPartition(jtIPartitions[partNo],C_ri_to_C,prv_nm,partNo,
			  jtIPartitions[partNo+1],C_li_to_C,"Cu0",partNo+1);
    partNo++;
    prv_nm = "Cu0";
    ceGatherIntoRoot(jtIPartitions[partNo],
		     C_ri_to_E,
		     Cu0_message_order,
		     prv_nm,partNo);
  }

  ceSendToNextPartition(jtIPartitions[partNo],C_ri_to_E,prv_nm,partNo,
			jtIPartitions[partNo+1],E_li_to_C,"E1",partNo+1);
  partNo++;
  ceGatherIntoRoot(jtIPartitions[partNo],
		   E_root_clique,
		   E1_message_order,
		   "E1",partNo);

}


void
JunctionTree::deScatterOutofRoot(JT_InferencePartition& part,
				 const unsigned root,
				 vector< pair<unsigned,unsigned> >& message_order,
				 const char*const part_type_name,
				 const unsigned part_num)
{
  const unsigned stopVal = (unsigned)(-1);
  part.maxCliques[root].
    deScatterToOutgoingSeparators(part);
  for (unsigned msgNo=(message_order.size()-1);msgNo != stopVal; msgNo --) {
    const unsigned to = message_order[msgNo].first;
    const unsigned from = message_order[msgNo].second;
    infoMsg(IM::Mod,"DE: message %s,part[%d]: clique %d --> clique %d\n",
	    part_type_name,part_num,from,to);
    part.maxCliques[to].
      deReceiveFromIncommingSeparator(part);
    part.maxCliques[to].
      deScatterToOutgoingSeparators(part);
  }
}


void
JunctionTree::deReceiveToPreviousPartition(JT_InferencePartition& next_part,
					   const unsigned next_part_leaf,
					   const char*const next_part_type_name,
					   const unsigned next_part_num,
					   JT_InferencePartition& previous_part,
					   const unsigned previous_part_root,
					   const char*const previous_part_type_name,
					   const unsigned previous_part_num)
{
  infoMsg(IM::Mod,"DE: message %s,part[%d],clique(%d) -> %s,part[%d],clique(%d)\n",
	  next_part_type_name,next_part_num,next_part_leaf,
	  previous_part_type_name,previous_part_num,previous_part_root);
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
JunctionTree::distributeEvidence()
{
  // this routine handles all of:
  // unrolled 0 times: (so there is a single P1,Cu0, and E1)  
  // unrolled 1 time: so there is a P1, C1, C3, E1
  // unrolled 2 or more times: so there is a P1 C1 [C2 ...] C3, E1

  const unsigned numCoPartitions = jtIPartitions.size()-3;
  unsigned partNo = jtIPartitions.size()-1;
  // set up appropriate name for debugging output.
  const char* prv_nm;
  const char *frst_name = "Cu0";

  deScatterOutofRoot(jtIPartitions[partNo],
		     E_root_clique,
		     E1_message_order,
		     "E1",partNo);
  partNo--;

  prv_nm = ((!gm_template.leftInterface)?"Cu0":"Co");
  deReceiveToPreviousPartition(jtIPartitions[partNo+1],E_li_to_C,"E1",partNo+1,
			       jtIPartitions[partNo],C_ri_to_E,prv_nm,partNo);

  if (!gm_template.leftInterface) {
    deScatterOutofRoot(jtIPartitions[partNo],
		       C_ri_to_E,
		       Cu0_message_order,
		       "Cu0",partNo);
    partNo--;
    prv_nm = ((numCoPartitions>0)?"Co":"P1");
    // if there's a Cu0 here, it doesn't occur at the beginning.
    frst_name = "P1";
    deReceiveToPreviousPartition(jtIPartitions[partNo+1],C_li_to_C,"Cu0",partNo+1,
				 jtIPartitions[partNo],C_ri_to_C,prv_nm,partNo);

  }

  for (unsigned p=numCoPartitions; p > 0; p--) {
    deScatterOutofRoot(jtIPartitions[partNo],
		       C_ri_to_C,
		       Co_message_order,
		       "Co",partNo);
    partNo--;
    prv_nm = ((p == 1)?frst_name:"Co");
    deReceiveToPreviousPartition(jtIPartitions[partNo+1],C_li_to_C,"Co",partNo+1,
				 jtIPartitions[partNo],C_ri_to_C,prv_nm,partNo);
  }


  if (gm_template.leftInterface) {
    deScatterOutofRoot(jtIPartitions[partNo],
		       C_ri_to_C,
		       Cu0_message_order,
		       "Cu0",partNo);
    partNo--;
    deReceiveToPreviousPartition(jtIPartitions[partNo+1],C_li_to_P,"Cu0",partNo+1,
				 jtIPartitions[partNo],P_ri_to_C,"P1",partNo);

  }

  deScatterOutofRoot(jtIPartitions[partNo],
		     P_ri_to_C,
		     P1_message_order,
		     "P1",partNo);

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


#if 0

spare code:

  // build a JT
  // assign RVs.

  
  Partition part;

  part.cliques = cliques;

  for (unsigned clique = 0; clique <= part.clique.size(); clique ++) {
    part.cliques[i].neighbors.clear();
  }
  createPartitionJunctionTrees(part);

  JT_Partition jt_part(part);



  // code from above that is commented out but not yet deleted. 
  // it is kept here for now to clarify development.
  // tag: ABCDEFGHIJK at end of this file.

   // Idea: sort the separators in decreasing order, where
   // sorted by size of intersectio of this separator with all
   // others. That way, we first iterate over the sep that has
   // the max overlap with all other seps, next iterate over 
   // sep that has next max.
   // actually, first compute max, and then remove it,
   // then repeat (so this is an N^2 algorithm for "sorting").
   

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

