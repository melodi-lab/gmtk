/*
 * GMTK_MAXCLIQUE.h
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *  $Header$
 *
 *
 *   A clique class.
 *   Note: some texts define a 'clique' as any complete set
 *   while other texts define a 'clique' as a maximally
 *   complete set with respect to the subset operator (i.e., a
 *   clique is one such that no proper superset of the set
 *   of nodes is a clique). In order to avoid confusion,
 *   I adopt here the term 'maxclique' which corresponds
 *   to a maximally complete set. Note, however, that in this
 *   program, the concepts are such that 
 *
 *               'clique == maxclique != complete set'
 *
 *   meaning that cliques are taken to be max cliques.
 *
 * Copyright (C) 2001 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 *
 */

// TODO: perhaps create a subclass or member of maxClique at some point, rather than
// adding everything for exact inference to the base class.


#ifndef GMTK_MAXCLIQUE_H
#define GMTK_MAXCLIQUE_H

#if HAVE_CONFIG_H
#include <config.h>
#endif

#if 0

// USE_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL is now set by configure

////////////////////////////////////////////////////////////////////////
// Comment/Uncomment to optimize for speed/reducing memory usage using
// another trick that only works when pruning is turned on.
#ifndef USE_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL
#define USE_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL
#endif
////////////////////////////////////////////////////////////////////////
#endif

#include "general.h"
#include "vhash_set.h"
#include "vhash_map.h"
#include "logp.h"
#include "cArray.h"
#include "sArray.h"
#include "debug.h"
#include "fixed_filter.h"
#include "lms_filter.h"
#include "rls_filter.h"
#include "counted_ptr.h"

#include "GMTK_RV.h"
#include "GMTK_DiscRV.h"
#include "GMTK_PackCliqueValue.h"
#include "GMTK_SpaceManager.h"
#include "GMTK_FactorInfo.h"
#include "GMTK_ObservationFile.h"

#include <stdio.h>
#include <stdlib.h>

#include <vector>
#include <string>
#include <map>

class PartitionStructures;
class PartitionTables;
class SeparatorClique;
class ConditionalSeparatorTable;
// class ConditionalSeparatorTable::SharedLocalStructure;

class MaxCliqueTable;


#include "GMTK_CliqueValueHolder.h"
#include "GMTK_VHashMapUnsignedUnsignedKeyUpdatable.h"
#include "GMTK_MaxCliqueBase.h"


class MaxClique : public MaxCliqueBase {

  friend class FileParser;
  friend class GraphicalModel;
  friend class GMTemplate;
  friend class SeparatorClique;

 public:

  // basic constructor with a set of nodes
 MaxClique(set<RV*> arg) : MaxCliqueBase(arg) {}
  // TODO: Figure out polymorphicness of copy ctor

  // Clone constructor from another MaxClique, but that uses a new set
  // of random variables, and adjusts the frame of each new set of
  // random variable with offset
  MaxClique(MaxClique& from_clique,
	    vector <RV*>& newRvs,
	    map < RVInfo::rvParent, unsigned >& ppf,
	    const unsigned int frameDelta = 0);

  // clone more of clique than the above...
  MaxClique(MaxClique& from_clique,
	    vector <RV*>& newRvs,
	    map < RVInfo::rvParent, unsigned >& ppf,
	    unsigned int frameDelta,
	    bool dummy);

  void checkClique(MaxClique const &target);

  virtual ~MaxClique() { 
    // TODO: do this right so that it works with tmp values in these objects.
    //   right now, not deleting these causes a memory leak.
    // if (maxCEValuePredictor != NULL) 
    // delete maxCEValuePredictor; 
  }

   float weightInJunctionTree(const set<RV*>& unassignedInPartition,
			     const bool upperBound,
			     const bool moreConservative,
			     const bool useDeterminism,
			     vector< set<RV*> > *lp_nodes,
			     vector< set<RV*> > *rp_nodes) const { 
    // We pass in all the arguments needed to the static routine that
    // actually does the work.
    return computeWeightInJunctionTree(nodes,
				       assignedNodes,
				       cumulativeAssignedNodes,
				       unassignedIteratedNodes,
				       cumulativeUnassignedIteratedNodes,
				       unionIncommingCESeps,
				       unassignedInPartition,
				       lp_nodes,rp_nodes,
				       upperBound,
				       moreConservative,
				       useDeterminism);
  }

  // memory reporting
   void reportMemoryUsageTo(FILE *f);


   // TODO: write a destructor for this (not a problem for now since
  // we only create a fixed set of maxcliques per GMTK run). 

  // USED ONLY IN JUNCTION TREE INFERENCE
  // The set of nodes that are *assigned* to this maxClique. "Assigned"
  // means that the node and its parents exist in this clique, but
  // assigned does not necessarily mean that such a node will be used
  // to contribute a probability to this clique potential. Assigned
  // nodes will have the option of using CPT-based iteration in certain
  // cases, unless it can't (i.e., the node is also in the incomming 
  // CE separator) or it is not worth it (the node is not sparse).
  // Computed in JunctionTree::assignRVsToCliques().
  // Used to:
  //   1) compute cumulativeAssignedNodes
  //   2) compute unassignedIteratedNodes (and its cumulative verison)
  //   3) compute jt weight
  set<RV*> assignedNodes;

  // USED ONLY IN JUNCTION TREE INFERENCE
  // A topologically sorted vector version of assignedNodes except for
  // those nodes with disposition AN_CONTINUE (since we don't want to
  // iterate those as it would be wasteful).  The topological sort is
  // relative only to the variables in this clique and relative to
  // those that are not disposition AN_CONTINUE. The topological sort
  // is used during CE clique creation, so that CPT iteration can be
  // used when available.  
  // Computed in: 
  //     JunctionTree::assignRVsToCliques() or JunctionTree::sortCliqueAssignedNodesAndComputeDispositions()
  // Used to: 
  //    1) compute assigned nodes dispositions (in appropriate order) 
  //    2) compute fSortedAssignedNodes in an inference clique
  vector<RV*> sortedAssignedNodes;

  // USED ONLY IN JUNCTION TREE INFERENCE
  // continuation scores to use for clique build pruning.
  // this array is one larger than sortedAssignedNodes, and for
  // node at position i (up to and including the length) gives
  // the best possible continuation score.
  sArray<logpr> sortedAssignedContinuationScores;


  // USED ONLY IN JUNCTION TREE INFERENCE
  // These are the nodes in this clique that are used to
  // contribute a probability to the clique potential. Necessarily, it is the
  // case that assignedProbNodes <= assignedNodes.
  // Computed in JunctionTree::assignRVsToCliques().
  // Used to:
  //   1) compute cumulativeAssignedProbNodes
  //   2) determine nodes to assign to clique
  //   3) compute assigned nodes dispositions 
  //   4) which nodes in cliques during EM get sent posterior probability.
  set<RV*> assignedProbNodes;


  // USED ONLY IN JUNCTION TREE INFERENCE
  // The set of assigned nodes cummulative in the JT relative to root
  // in the JT for the current partition and clique.  Note that
  // cumulativeAssignedNodes *DOES NOT* include the assignedNodes in the current
  // clique.  In other words, we have that:
  // intersect(cumulativeAssignedNodes,assignedNodes) = EMPTY
  // Also note, that during inference, this includes all
  // assigned nodes in the previous partition, but during jt-weight
  // evaluation in triangulation, it does not include those previous
  // partition nodes. Instead, those nodes are a subset of
  // JT_Partition::unassignedInPartition.
  // Computed in JunctionTree::getCumulativeAssignedNodes
  //    via JunctionTree::assignRVsToCliques()
  // Used to:
  //   1) compute JT weight
  //   2) compute assigned nodes dispositions 
  set<RV*> cumulativeAssignedNodes;

  // USED ONLY IN JUNCTION TREE INFERENCE
  // The set of assigned probability nodes cummulative in the JT relative to root
  // in the JT for the current partition. This *DOES NOT* include
  // the set of nodes in the current clique that contribute
  // probability. In other words, we have that:
  // intersect(cumulativeAssignedProbNodes,assignedProbNodes) = EMPTY
  // Computed in JunctionTree::getCumulativeAssignedNodes
  //    via JunctionTree::assignRVsToCliques()
  // Used to:
  //     1) select clique for later nodes to be assigned to.
  set<RV*> cumulativeAssignedProbNodes;


  // USED ONLY IN JUNCTION TREE INFERENCE
  // The union of of all nodes in separators for incomming messages
  // for this clique during collect evidence stage. By 'incomming',
  // what is meant is incomming during the collect evidence stage.
  // Computed in JunctionTree::computeSeparatorIterationOrder()
  // Used to:
  //   1) when computing the sep sets division between previous
  //      separator iteration, and remainder, this is
  //      built up accumulatively. When done, it contains
  //      the union of all CE incoming separators, but 
  //      is not currently used outside the routine it is created in
  //      (other than printing).
  //   2) computing  unassignedIteratedNodes
  //   3) compute assigned nodes dispositions 
  set<RV*> unionIncommingCESeps;

  // USED ONLY IN JUNCTION TREE INFERENCE
  // The set of nodes that are not assigned to this clique but that
  // also are *NOT* iterated over by the incomming
  // separators. Therefore, these nodes in the clique must be iterated
  // over [0,card-1] unfortunately.  The hope is that assignment and setup can be
  // done so there are very few or zero of such nodes. Note however
  // that once a clique iterates over such a node, any later cliques
  // (closer to the JT root) will not have to re-do this as these
  // assignments here will be represented in this cliques out-going
  // sepset --- this is guaranteed by the running intersection property.
  // Computed in JunctionTree::computeSeparatorIterationOrders()
  // Used to:
  //  1) compute JT weight in clique
  //  2) compute fUnassignedIteratedNodes in inference clique (which says
  //     which nodes in a clique are unassigend so iterated [0,card-1]
  //  3) in CE, choose which routine to call first.
  set<RV*> unassignedIteratedNodes;


  // USED ONLY IN JUNCTION TREE INFERENCE 
  // The set of nodes that are simply not assigned. In other words,
  // these nodes can be defined as:
  // 
  //    unassignedNodes = nodes - assignedNodes;
  //
  // Computed in JunctionTree::assignRVsToCliques().
  // Used to:
  //  1) These nodes are used for iteration when doing
  //     clique driven inference.
  set<RV*> unassignedNodes;


  // USED ONLY IN JUNCTION TREE INFERENCE 
  // This enum defines the different things that could
  // happen to an assigned node in a clique during CE assigned node iteration.
  // TODO: come up with better (and only one name) for each of the below.
  enum AssignedNodeDisposition {
    // In the discussion below, say v is an assignedNode (meaning it
    // exists with its parents in the clique).
    // We have 6 cases, each case has 1 or more name.

    //
    // 0)
    //   v is not in sepNodes, and v is a prob node as well, v either
    //   sparse or dense:  meaning we multiply this clique potential by
    //   prob(v|parents(v)). We use CPT iteration to iterate over v
    //   and compute and apply probabilities.
    AN_NOTSEP_PROB_SPARSEDENSE=0,
    AN_CPT_ITERATION_COMPUTE_AND_APPLY_PROB=0,

    //
    // 1) 
    //   v is not in sepNodes, and v is not a prob node, and V is
    //   sparse. We still use CPT iteration to iterate over v given
    //   parents, but do not multiply clique potential by
    //   probabilties. By using CPT iteration, we automatically remove any zeros
    //   now.
    AN_NOTSEP_NOTPROB_SPARSE=1,
    AN_CPT_ITERATION_COMPUTE_PROB_REMOVE_ZEROS=1,

    //
    // 2)
    //   v is not in sepNodes, and v is not a prob node, and V is NOT
    //   sparse.  We iterate using [0,card-1] and continue, and do not
    //   change any probabilities. Note, we could still
    //   do CPT iteration in this case, as there might be zeros
    //   in a dense CPT that are being missed due to the fact
    //   that we're doing card iteration. I.e., this case
    //   could be set to value 1. TODO: test this.
    AN_NOTSEP_NOTPROB_DENSE=2,
    AN_CARD_ITERATION=2,

    //
    // 3)
    //   v is in sepNodes, and v is a probability node, v either
    //   sparse or dense. Then we compute the probability for v, apply
    //   it, and continue only if the probabilty is non zero.
    AN_SEP_PROB_SPARSEDENSE=3,
    AN_COMPUTE_AND_APPLY_PROB=3,

    //
    // 4)
    //   A) v is in sepNodes, v is not a probability node, v either
    //   sparse or dense, but v *WAS* assigned in some previous JT
    //   clique. In this case, we just continue on since we know at
    //   least one of the incomming separators killed off any zero
    //   entries, i.e., no need to use the CPT, even when this node is
    //   sparse since the zero was used in a previous clique.
    AN_SEP_NOTPROB_SPARSEDENSE_PRVASSIGNED=4,
    //   -
    //   B) v is in sepNodes, v is not a probability node, v not sparse,
    //   but v was *NOT* assigned in a previous JT clique. We just
    //   continue, since the CPT won't tell us anything.
    AN_SEP_NOTPROB_DENSE_NOTPRVASSIGNED=4,
    AN_CONTINUE=4,

    //
    // 5)
    //   v is in sepNodes, v is not a probability node, v sparse, but
    //   v was *NOT* assigned in a previous JT clique. Then, since v
    //   is sparse, we check here to make sure that the probability of
    //   v|parents is not zero. If it is zero, we backout, but if it
    //   is not zero we continue, but *do not* multiply the
    //   probability into this clique potential.
    AN_SEP_NOTPROB_SPARSE_NOTPRVASSIGNED=5,
    AN_CONTINUE_COMPUTE_PROB_REMOVE_ZEROS=5
  };


  static char const *disp2str(AssignedNodeDisposition d) {
    switch (d) {
    case 0: return "AN_NOTSEP_PROB_SPARSEDENSE/AN_CPT_ITERATION_COMPUTE_AND_APPLY_PROB";
    case 1: return "AN_NOTSEP_NOTPROB_SPARSE/AN_CPT_ITERATION_COMPUTE_PROB_REMOVE_ZEROS";
    case 2: return "AN_NOTSEP_NOTPROB_DENSE/AN_CARD_ITERATION";
    case 3: return "AN_SEP_PROB_SPARSEDENSE/AN_COMPUTE_AND_APPLY_PROB";
    case 4: return "AN_SEP_NOTPROB_SPARSEDENSE_PRVASSIGNED/AN_SEP_NOTPROB_DENSE_NOTPRVASSIGNED/AN_CONTINUE";
    case 5: return "AN_SEP_NOTPROB_SPARSE_NOTPRVASSIGNED/AN_CONTINUE_COMPUTE_PROB_REMOVE_ZEROS";
    default: assert(false);
    }
  }

  // USED ONLY IN JUNCTION TREE INFERENCE 
  // Predicate for each each sorted assigned node, used to say what
  // should be done to that node during assigned node iteration.
  // There are a number of different options, and they are documented
  // in the routine. The options are specified by the corresponding
  // enum also defined above. Note that if this
  // array is of size 0, then the default case is assumed
  // to be AN_CPT_ITERATION_COMPUTE_AND_APPLY_PROB
  // computed in MaxClique::computeSortedAssignedNodesDispositions()
  //     which is called in JunctionTree::getCumulativeUnassignedIteratedNodes()
  // Used to:
  //   1) determine way in which fSortedAssignedNodes in inf clique should be iterated.
  sArray < AssignedNodeDisposition > dispositionSortedAssignedNodes;


  // TODO: TRY TO REMOVE THIS AS IT IS REDUNDANT!!
  // USED ONLY IN JUNCTION TREE INFERENCE
  // the preceding cumulative set of unassigned and iterated nodes
  // relative to the root in the JT. During JT weight scoring,
  // this only includes the current partition, but during
  // actual inference (and jt_info.txt file creation), this will also include 
  // previous partitions as well.
  // 
  // This var is used to determine which of the assigned nodes
  // need actually be iterated within a clique, and which are already
  // set by the separator iterations. Note that
  // precedingUnassignedIteratedNodes does *NOT* include the
  // unassignedIteratedNodes in the current clique (thus the word
  // preceding above).  
  // Computed in JunctionTree::getCumulativeUnassignedIteratedNodes()
  // Used to:
  //    1) comptue jt weight
  // Note: should be replaced by unionIncommingCESeps, perhaps remove.
  set<RV*> cumulativeUnassignedIteratedNodes;


  // USED ONLY IN JUNCTION TREE INFERENCE
  // These are the nodes that are hidden (i.e., they are
  // the non-continous hidden variables). These are
  // the nodes whose values could be hashed.
  // Computed in MaxClique::prepareForUnrolling()
  // Used to:
  //   1) compute hashableNodes
  vector<RV*> hiddenNodes;


  // USED ONLY IN JUNCTION TREE INFERENCE
  // These are the hidden nodes that are hashed (i.e., they are
  // the nodes that when unpacked, and along with
  // the observations, are such that any remaining
  // nodes in the clique can be determined uniquely (e.g.,
  // deterministic nodes)
  // Computed in MaxClique::prepareForUnrolling()
  // Used to:
  //   1) set the size of the clique packer
  //   2) set the size of discreteValuePtrs in inference clique
  //      to get quick access to the set of RV vals.
  vector<RV*> hashableNodes;

  // USED ONLY IN JUNCTION TREE INFERENCE
  // These are the nodes that given the hashable nodes
  // and any observations can be determinied deterministically.
  // Note that if this vector is empty, then hashableNodes contains
  // everything to be hashed. The nodes are so ordered that
  // these  nodes should be determined in order (i.e., later ones
  // might use earlier ones as parents).
  // Computed in MaxClique::prepareForUnrolling()
  // Used to:
  //   1) when unpacking clique, compute the rest of the clique values
  //      (by determining these variables deterministic values)
  vector<RV*> determinableNodes;


  // USED ONLY IN JUNCTION TREE INFERENCE
  // Clique neighbors in a junction tree, integer indexes into a table
  // of cliques contained in the Partition class containing this
  // maxclique. By storing ints rather than actual pointers to
  // maxcliques, cloned partitions (with same corresponding clique
  // array ordering) can use the same structure (ints offsets 
  // in an array) on different cliques in unrolled partitions.
  // Computed in JunctionTree::createPartitionJunctionTree()
  vector<unsigned> neighbors;

  // USED ONLY IN JUNCTION TREE INFERENCE
  // Clique children in a junction tree relative to the JT that has
  // been rooted at the root of the partition. The variable contains
  // integer indexes into a table of cliques contained in the
  // Partition class containing this maxclique. By storing ints rather
  // than actual pointers to maxcliques, cloned partitions (with same
  // corresponding clique array ordering) can use the same structure
  // on different cliques. This variable is used only by the code that
  // assigns CPTs to cliques, so the variable DOES NOT SURIVE during a
  // clique cloning.
  // Computed in JunctionTree::createDirectedGraphOfCliques()
  vector<unsigned> children;

  // USED ONLY IN JUNCTION TREE INFERENCE
  // set of separators that we receive from in the collect evidence
  // stage.  (equivalently, the set of separators we send to in
  // distribute evidence stage) Again, ints indexing into parent
  // partition. Note that if this clique is the LI of the partition,
  // then there will be one extra separator at the end of this list
  // that corresponds to the separator that links us to the left
  // partition.  This is of course only true for the C and E
  // partitions.
  // 
  // Created in JunctionTree::createSeparators().
  // Modified/Sorted in JunctionTree::computeSeparatorIterationOrder()
  vector<unsigned> ceReceiveSeparators;

  // USED ONLY IN JUNCTION TREE INFERENCE
  // The separator that we send to in the collect collect evidence
  // stage (equivalently, the separator that we receive from in the
  // distribute evidence stage).  Again, ints indexing into parent
  // partition. Set to ~0x0 when the appropriate separator lives in
  // another partition (and so needs to be explicitly given).
  vector<unsigned> ceSendSeparators;

  // USED ONLY IN JUNCTION TREE INFERENCE
  // structure used to pack and unpack clique values for this
  // clique. Shared by all infernece instances of this clique.
  PackCliqueValue packer;

  // USED ONLY IN JUNCTION TREE INFERENCE 
  // Structure to allocate clique values from.  These things are
  // allocated and deleted in bulk (rather than individually) and that
  // is handled by this object, thereby acting effectively as a
  // customized memory management unit.
  // 
  // Note that if a packed clique value can be held in less than
  // IMC_NWWOH*sizeof(unsigned)*8 bits, then this isn't used since the
  // direct storage (overlaped with a pointer) is used instead.
  CliqueValueHolder valueHolder;

  // USED ONLY IN JUNCTION TREE INFERENCE
  // hash table that holds the clique values (a packed vector of RV
  // values). The actual values are pointers to within
  // valueHolder. This hash table is used to look up a clique value
  // quickly and see if it is already there (so that storage for
  // clique the same values over all instances of a clique in an
  // unrolling are stored only once (so that memory for clique value
  // storage stays roughly constant for any amount of unrolling).
  //
  // Note that if a packed clique value can be held in < IMC_NWWOH*
  // sizeof(unsigned) bytes, then this isn't used since the direct
  // storage (overlaped with pointer) can be used instead.
  vhash_set< unsigned > cliqueValueHashSet;


#ifdef USE_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL
  // USED ONLY IN JUNCTION TREE INFERENCE
  // this next array is used when creating a clique. It is used when
  // pruning is turned on, and when a clique insert has occured,
  // the clique values are stored here first, before putting in the global
  // shared pool. Then after pruning, only those clique values that correspond
  // to clique value scores that survive pruning are added to the globally shared
  // poll, thus saving memory. The array is stored here so that we don't need
  // to reallocate it each time we create a clique if we don't want to.
  // CliqueValueHolder localValueHolder;
  sArray < unsigned > temporaryCliqueValuePool;
#endif

  // USED ONLY IN JUNCTION TREE INFERENCE
  // Manages and memorizes the size and space requests made
  // by all corresponding MaxCliqueTables. This way,
  // the next time an MaxCliqueTables asks for an initial
  // amount of memory, we'll be able to give it something
  // closer to what it previously asked for rather than
  // something too small.
  SpaceManager cliqueValueSpaceManager;

  // USED ONLY IN JUNCTION TREE INFERENCE An array containing the
  // partial/cumulative probabilities of partial clique values and is
  // used by the non-recursive version of the clique value iterationt
  // through clique assigned nodes. In particular, used in:
  //    void ceIterateAssignedNodesNoRecurse()
  // Should be one greater than the size of: sortedAssignedNodes.size()
  sArray< logpr > probArrayStorage;

  // USED ONLY IN JUNCTION TREE INFERENCE
  // structure to store the VE separators for this clique, if there are 
  // any.
  struct VESepInfo {
    vector< RV* > parents;
    // if this is a case of just parents and child,
    // then grandChild is NULL. If this is parents, child, and child/grandChild relationship,
    // then grandChild is non NULL.
    RV* child;
    RV* grandChild;
    VESepInfo() {
      child = grandChild = NULL;
    }
  };
  vector < VESepInfo > veSeparators;

  // USED ONLY IN JUNCTION TREE INFERENCE
  // The set of factors/constraints that are assigned to this
  // clique. In this case 'assigned' means that the factor is covered
  // by this max clique (meaning the max clique variables are a
  // superset of the factor variables). Any assigned factor might help
  // to remove zeros from, or its score might be used to affect
  // pruning in the clique, but any score provided by the factor is
  // not necessarily applied.  Again, ints indexing into parent
  // partition's factor array.
  // 
  // Created in JunctionTree::assignFactorsToCliques().
  vector<unsigned> assignedFactors;

  // USED ONLY IN JUNCTION TREE INFERENCE
  // The set of factors/constraints that are assigned and scored in
  // this clique. In this case 'assigned' and 'scored' means that the
  // factor is covered by this max clique (meaning the max clique
  // variables are a superset of the factor variables), and that
  // assigned factor might not only help to remove zeros, pruning, but
  // it will also provide a numeric score (e.g., log probability) that
  // is applied to any clique entry.  Again, ints indexing into parent
  // partition's factor array.
  // 
  // Created in JunctionTree::assignFactorsToCliques().
  vector<unsigned> assignedScoredFactors;


  // USED ONLY IN JUNCTION TREE INFERENCE
  // a value that is roughly equal to the order of the state space of
  // the clique. It should be smaller than the clique state space,
  // however, since pruning will make the state space smaller.
  // This will be used for initial memory allocation units, etc.
  // by this clique object and clones of this clique object.
  // It is only usable after we have prepared for unrolling.
  // TODO: add proper allocation statistics object and use that.
  // unsigned allocationUnitChunkSize;


  // Clear up the things that are just to create and hold information
  // about this maxclique being used in a junction tree.
  void clearJTStructures() {
    assignedNodes.clear();
    sortedAssignedNodes.clear();
    assignedProbNodes.clear();
    cumulativeAssignedNodes.clear();
    unionIncommingCESeps.clear();
    unassignedIteratedNodes.clear();
    dispositionSortedAssignedNodes.clear();
    cumulativeUnassignedIteratedNodes.clear();
    hiddenNodes.clear();
    neighbors.clear();
    children.clear();
    ceReceiveSeparators.clear();
  }


  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////

  // USED ONLY IN JUNCTION TREE INFERENCE
  // Prepare the last set of data structures so that JT inference
  // 'clones' of this object can be unrolled and inference can occur.
  void prepareForUnrolling();

  // USED ONLY IN JUNCTION TREE INFERENCE
  // compute the way in which assigned nodes are iterated (or not).
  bool computeSortedAssignedNodesDispositions();

  // USED ONLY IN JUNCTION TREE INFERENCE
  // sort the node and assign the dispositions.
  void sortAndAssignDispositions(); // new

  void sortAndAssignDispositions(const char *varCliqueAssignmentPrior); // old


  // USED ONLY IN JUNCTION TREE INFERENCE
  // print out everything in this junction tree clique to a file.
  void printAllJTInfo(FILE* f,const unsigned indent,const set<RV*>& unassignedInPartition,
		      const bool upperBound,
		      const bool moreConservative,
		      const bool useDeterminism,
		      vector< set<RV*> > *lp_nodes,
		      vector< set<RV*> > *rp_nodes);
  void printAllJTInfo(FILE* f,const unsigned indent,const set<RV*>& unassignedInPartition,
		      const bool upperBound,
		      const bool moreConservative,
		      const bool useDeterminism,
		      set<RV*> *lp_nodes,
		      set<RV*> *rp_nodes);

  // USED ONLY IN JUNCTION TREE INFERENCE
  // used to clear out hash table memory between segments
  void clearCliqueValueCache(bool force = false);

  // USED ONLY IN JUNCTION TREE INFERENCE
  // used to compute the unassigned nodes in this clique.
  void computeUnassignedCliqueNodes();

  // USED ONLY IN JUNCTION TREE INFERENCE
  // computes initial information about any VE separators that
  // might exist in this clique, and sets up any data structures
  // needed to produce these separators. Returns the number of VE separators.
  // See the routine code for the definition of a VE separator.
  unsigned computeVESeparators();

  // returns number of VE seps, but only after computeVESeparators() is called.
  unsigned numVESeparators() { return veSeparators.size(); }


  // USED ONLY IN JUNCTION TREE INFERENCE
  // resets variables that might be used in inference for a "round" (meaning
  // a full run of collect or collect/distribute evidence.
  void prepareForNextInferenceRound() {
    prevMaxCEValue.valref() = (-LZERO);
    prevPrevMaxCEValue.valref() = (-LZERO);
    if (maxCEValuePredictor.ptr() != NULL)
      maxCEValuePredictor->init();
  }

  // USED ONLY IN JUNCTION TREE INFERENCE
  // Clears or resets to initial status all memory that was created
  // for this clque during an inference run.
  void clearInferenceMemory() {
    valueHolder.prepare();
    cliqueValueHashSet.clear();
#ifdef USE_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL
    temporaryCliqueValuePool.clear();
#endif
  }

  void clearCliqueAndIncommingSeparatorMemory(SeparatorClique& sep);

  // return the s'th VE separator. Must call computeVESeparators()
  // first. Caller assumed to copy this out since we return a
  // reference.
  VESepInfo& VESeparatorInfo(unsigned s) { 
    assert ( s < veSeparators.size() );
    return veSeparators[s];
  }


};


#endif
