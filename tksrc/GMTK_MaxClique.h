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
 * Copyright (c) 2001, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about the suitability of this software
 * for any purpose. It is provided "as is" without express or implied warranty.
 *
 */

// TODO: perhaps create a subclass or member of maxClique at some point, rather than
// adding everything for exact inference to the base class.


#ifndef GMTK_MAXCLIQUE_H
#define GMTK_MAXCLIQUE_H


////////////////////////////////////////////////////////////////////////
// Comment/Uncomment to optimize for speed/reducing memory usage using
// another trick that only works when pruning is turned on.
#ifndef USE_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL
#define USE_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL
#endif
////////////////////////////////////////////////////////////////////////


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
class MaxClique;
class MaxCliqueTable;

/////////////////////////////
// TODO: think of a way of subclassing or something so that members
// are only in the objects when they are used. (perhaps one base class
// with two children will work).
/////////////////////////////


class CliqueValueHolder  {

  // The packed clique value size (in unsigned).
  // @@@ sould be const
  // const unsigned cliqueValueSize;
  unsigned cliqueValueSize;

  // The amount that successive allocation chunks grow in size.  Must
  // be >= 1.0.
  // Ideally, should be a const.
  float growthFactor;

  // The size of the initial allocation unit. A chunk is an
  // array of packed clique values. When we allocate the n'th
  // chunk, we allocate k(n) new packed clique values, 
  // where k(n) = allocationUnitChunkSize*growthFactor^(n-1)
  // Ideally, should be a const.
  unsigned allocationUnitChunkSize;

  // manage the space here.
  // SpaceManager spaceManager;

  // Current capacity, total number of unsigned values that this
  // object can currently potentially hold max without a resize.
  unsigned capacity;

  // a chunk, i.e., matrix of unsigned numbers constituting
  // allocationUnitChunkSize*growthFactor^(n-1) packed clique values
  // for some n.
  typedef sArray< unsigned > AllocationChunk;

  // Array of chunks.
  cArray< AllocationChunk >  values;

  // The the pointernext position in the current chunk to obtain a clique value
  // to use.
  unsigned* curAllocationPosition;

  // the end position in the current chunk, meaning
  // that we need to reallocate
  unsigned* curAllocationEnd;

  // total number of values currently used (not really
  // needed but kept anyway for debugging).
  unsigned numAllocated;




public:

  // create an empty object to re-construct later
  CliqueValueHolder() {}
  
  // real constructor
  CliqueValueHolder(unsigned cliqueValueSize,
		    unsigned allocationUnitChunkSize,
		    float growthFactor=1.25);

  ~CliqueValueHolder() { makeEmpty();  }

  // clear out all existing memory, and get ready for next use.
  void prepare();

  // Empty out and free up all memory, and reset to having
  // nothing added. Make invalid as well.
  void makeEmpty();

  // return a pointer to the next unused clique value
  // for scratch (etc.) without actually allocating it.
  unsigned* curCliqueValuePtr() { return curAllocationPosition; }

  // Actually "allocate" the space pointed to by curCliqueValuePtr()
  // and then advance curCliqueValuePtr() to the next position.
  void allocateCurCliqueValue();

  // return the total number of bytes requested to the OS memory system by this structure.
  unsigned long bytesRequested() { return (unsigned long)capacity*sizeof(unsigned); }

};


// a special vhash class, for mapping from keys consisting of
// compressed sets of RV values, to items consisting of array indices.
typedef vhash_map < unsigned, unsigned > VHashMapUnsignedUnsigned;
class VHashMapUnsignedUnsignedKeyUpdatable : public VHashMapUnsignedUnsigned {
public:
  //////////////////////
  // constructor for empty invalid object. 
  // WARNING: this will create an invalid object. It is assumed
  // that this object will re-reconstructed later.
  VHashMapUnsignedUnsignedKeyUpdatable() {}
  // constructor
  VHashMapUnsignedUnsignedKeyUpdatable(const unsigned arg_vsize,
				       unsigned approximateStartingSize = 
				       HashTableDefaultApproxStartingSize);

  //
  // Direct access to tables and keys (made available for speed).  Use
  // sparingly, and only if you know what you are doing about the
  // internals of a hash table. Note also that this breaks
  // encapsulation, meaning that if the implementation of the
  // internals of the parent hash table change, this code might break.
  unsigned*& tableKey(const unsigned i) { return table.ptr[i].key; }
  unsigned& tableItem(const unsigned i) { return table.ptr[i].item; }
  bool tableEmpty(const unsigned i) { return table.ptr[i].empty(); }
  unsigned tableSize() { return table.size(); }

};


//////////////////////////////////////////////////////////////////////////////
// Number of words to exist in a packed clique without a hash. These
// may be set to zero to turn on hashing even when packed clique
// values take up less than 1 machine word. To turn off hashing
// always, set to something larger than the largest packed clique
// value in machine words (e.g., 2, 3, etc. but this might consume
// lots of memory). The number of bits required for a packed clique
// value is equal to:
//
//    num_bits_required = \sum_{v \in h(C)} ceil(log2(card(v)))
//
// where C is a clique, h(C) is all *hidden* variables in the clique C
// (i.e., vars with card > 1), and card(v) is the cardinality of the
// variable v (the card is taken from the structure file). Note, we
// use h(C) so that observed variables (or vars with card = 1) are not
// needlessly stored in the clique value.
// --
// MaxCliqueTable Number Words WithOut a Hash (IMC_NWWOH): Namely,
// the number of words that can be stored directly as a packed clique
// value before we resort to using a shared hash table for all
// instances of this origin clique.
#define IMC_NWWOH (1)
// ConditionalSeparatorTable Number Words WithOut a Hash (ISC_NWWOH):
// Namely, the number of words that can be stored directly as a packed
// clique value before we resort to using a shared hash table for all
// instances of this origin clique.  One for the accumulated
// Intersection packed values
#define ISC_NWWOH_AI (1)
// And the same for the remainder in a Separator.
#define ISC_NWWOH_RM (1)
// -- 
//////////////////////////////////////////////////////////////////////////////


class MaxClique : public IM {

  friend class FileParser;
  friend class GraphicalModel;
  friend class GMTemplate;
  friend class SeparatorClique;


public:
  // Thresholds for -cpbeam pruning. We keep previous clique max value
  // and previous previous clique max value around.
  logpr prevMaxCEValue;
  logpr prevPrevMaxCEValue;
  counted_ptr<AdaptiveFilter> maxCEValuePredictor;

  double prevMaxCEValPrediction;


  // beam width for clique-based beam pruning.
  static double cliqueBeam;
  // beam width to use while building cliques for partial clique pruning.
  static double cliqueBeamBuildBeam;
  // properties of the filter used for partial clique pruning.
  static char* cliqueBeamBuildFilter;
  // use clique cont heuristic
  static bool cliqueBeamContinuationHeuristic;
  // clique build beam expansion factor 
  static double cliqueBeamBuildExpansionFactor;
  // clique build beam maximum number of expansions
  static unsigned cliqueBeamBuildMaxExpansions;
  // clique beam cluster pruning number of clsuters, 0 or 1 to turn off.
  static unsigned cliqueBeamClusterPruningNumClusters;
  // beam width for cluster-based beam pruning.
  static double cliqueBeamClusterBeam;
  // state beam width for state-based cluster pruning
  static unsigned cliqueBeamClusterMaxNumStates;

  // TODO: @@@@ remove this variable, it is here just to gather a few stats for
  // Wed Dec 10 10:47:27 2008 Friday's talk. This will *NOT* work
  // with multiple cliques per partition.
  double prevFixedPrediction;



  // forced max number of states in a clique. Set to 0 to turn it off.
  static unsigned cliqueBeamMaxNumStates;
  // fraction of clique to retain, forcibly pruning away everything else. Must be

  // between 0 and 1 (i.e., 0 < v <= 1).
  static float cliqueBeamRetainFraction;
  // a version for clustered states
  static float cliqueBeamClusterRetainFraction;

  // between 0 and 1 (i.e., 0 < v <= 1).
  static double cliqueBeamMassRetainFraction;
  // exponentiate the clique scores before pruning by this amount.
  static double cliqueBeamMassExponentiate;
  // min possible resulting size.
  static unsigned cliqueBeamMassMinSize;
  // additional normal beam to include once mass is covered.
  static double cliqueBeamMassFurtherBeam;

  // a version of the above for, but for cluster/diversity pruning
  static double cliqueBeamClusterMassRetainFraction;
  static double cliqueBeamClusterMassExponentiate;
  static unsigned cliqueBeamClusterMassMinSize;
  static double cliqueBeamClusterMassFurtherBeam;


  // amount to re-sample from pruned clique. =0.0 to turn off.
  static double cliqueBeamUniformSampleAmount;


  // set to true to store all deterministic children that exist in
  // this clique with its parents in the clique sorage. Otherwise, set
  // to false so that the det. children are not stored (which saves
  // space and might but does not necessarily slow things down).
  static bool storeDeterministicChildrenInClique;


  // if this is turned on, GMTK will normalize the score of each clique so
  //  that we do not run out of numeric dynamic range. This will
  //  render prob(evidence) meaningless, but will give the same
  //  results for training or decoding. The special values
  //  f the variabale are:
  //      val == 1 ==> don't do any score normalize (default)
  //      val == 0 ==> normalize (divide) by the maximum clique value
  //      val > 1  ==> normalize (divide) by the given value. A good
  //                value to use woudl be exp(- log(prob(evidence)) / T ) 
  //                and the per-frame log-likelyhood can be estimaed
  //                by a short-length run of the program.
  static double normalizeScoreEachClique;

  // @@@ need to take out, here for now to satisify STL call of vector.clear().
#if 0
  MaxClique& operator=(const MaxClique& f) {
    return *this;
  }
#endif

  // the set of nodes which form a max clique, in arbitrary order.
  set<RV*> nodes;

  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////

  // basic constructor with a set of nodes
  MaxClique(set<RV*> arg) {
    nodes = arg;
  }

  // Clone constructor from another MaxClique, but that uses a new set
  // of random variables, and adjusts the frame of each new set of
  // random variable with offset
  MaxClique(MaxClique& from_clique,
	    vector <RV*>& newRvs,
	    map < RVInfo::rvParent, unsigned >& ppf,
	    const unsigned int frameDelta = 0);


  ~MaxClique() { 
    // TODO: do this right so that it works with tmp values in these objects.
    //   right now, not deleting these causes a memory leak.
    // if (maxCEValuePredictor != NULL) 
    // delete maxCEValuePredictor; 
  }


  // Complete the max clique by adjusting the random variables
  // themselves.
  void makeComplete() { makeComplete(nodes); }
  // Static version of variable set completion, to complete set of
  // random variables passed in.
  static void makeComplete(const set<RV*> &rvs);



  // The penalty per feature that we pay to include a continous
  // observation in a clique. We include this since if we always
  // charge penalty of a factor of one, the continuous observation
  // might end up in a clique with many more than its own parents,
  // something that might possibly slow things down. The
  // default value is 0, meaning no additional penalty.
  static double continuousObservationPerFeaturePenalty;

  // Static version of various compute weight routines, useful in
  // certain places outside this class where we just have a collection
  // of nodes to compute the weight of.
  static float computeWeight(const set<RV*>& nodes,
			     const RV* node = NULL,
			     const bool useDeterminism = true);
  static float computeWeightWithExclusion(const set<RV*>& nodes,
					  const set<RV*>& unassignedIteratedNodes,
					  const set<RV*>& unionSepNodes,
					  const bool useDeterminism = true);
  static float computeWeightInJunctionTree(const set<RV*>& nodes,
					   const set<RV*>& assignedNodes,
					   const set<RV*>& cumulativeAssignedNodes,
					   const set<RV*>& unassignedIteratedNodes,
					   const set<RV*>& cumulativeUnassignedIteratedNodes,
					   const set<RV*>& separatorNodes,
					   const set<RV*>& unassignedInPartition,
					   set<RV*>* lp_nodes,
					   set<RV*>* rp_nodes,
					   const bool upperBound,
					   const bool moreConservative,
					   const bool useDeterminism);

  // compute the weight (log10 state space) of this clique.
  float weight(const bool useDeterminism = true) const { 
    return computeWeight(nodes,NULL,useDeterminism); 
  }
  float weightInJunctionTree(const set<RV*>& unassignedInPartition,
			     const bool upperBound,
			     const bool moreConservative,
			     const bool useDeterminism,
			     set<RV*>* lp_nodes,
			     set<RV*>* rp_nodes) const { 
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


  // print just the clique nodes
  void printCliqueNodes(FILE* f);


  // memory reporting
  void reportMemoryUsageTo(FILE *f);

  // score reporting
  void reportScoreStats();

  //////////////////////////////////////////////
  // TODO: figure out a way so that the member variables below exist only in
  // a subclass since all of the below is not needed for basic
  // triangulation.


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
  unsigned ceSendSeparator;

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
  void sortAndAssignDispositions(const char *varCliqueAssignmentPrior);


  // USED ONLY IN JUNCTION TREE INFERENCE
  // print out everything in this junction tree clique to a file.
  void printAllJTInfo(FILE* f,const unsigned indent,const set<RV*>& unassignedInPartition,
		      const bool upperBound,
		      const bool moreConservative,
		      const bool useDeterminism,
		      set<RV*>* lp_nodes,set<RV*>* rp_nodes);

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

  void numParentsSatisfyingChild(unsigned& num,unsigned par,vector <RV*> & parents, RV* child);

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


  // could start adding stuff for v3 inference here, and then clean up code later??

  

};



class SeparatorClique : public IM
{
  friend class FileParser;
  friend class GraphicalModel;
  friend class GMTemplate;
  friend class ConditionalSeparatorTable;
  friend class JunctionTree;

public:

  // beam width for separator-based beam pruning.
  static double separatorBeam;

  // the set of nodes which forms the separator
  set<RV*> nodes;

  //////////////////////////////////////////////////////////////////////
  // For a given iteration order, the intersection of 'nodes' with all
  // nodes in previous seperators according to the iteration.
  set<RV*> accumulatedIntersection;
  // hidden variables of the above
  vector<RV*> hAccumulatedIntersection;
  // structure used to pack and unpack clique values
  PackCliqueValue accPacker;
  // structure to allocate clique values from.  These things are
  // allocated and deleted in bulk (rather than individually) and that
  // is handled by this object, thereby acting effectively as a
  // customized memory management unit.
  CliqueValueHolder accValueHolder;
  // hash table that holds the clique values (a packed vector of RV
  // values). The actual values are pointers to within valueHolder.
  vhash_set< unsigned > accSepValHashSet;


  //////////////////////////////////////////////////////////////////////
  // remainder = nodes - accumulatedIntersection which is precomputed
  // in the non-unrolled version of a separator.
  set<RV*> remainder;
  // hidden variables of the above
  vector<RV*> hRemainder;
  // structure used to pack and unpack clique values
  PackCliqueValue remPacker;
  // structure to allocate clique values from.  These things are
  // allocated and deleted in bulk (rather than individually) and that
  // is handled by this object, thereby acting effectively as a
  // customized memory management unit.
  CliqueValueHolder remValueHolder;
  // hash table that holds the clique values (a packed vector of RV
  // values). The actual values are pointers to within valueHolder.
  vhash_set< unsigned > remSepValHashSet;

  ////////////////////////////////////////////////////////////
  // Data structures for when this is a VE separator. 
  // set to true if this is a virtual evidence (VE) separator.
  bool veSeparator;
  // information about this VE separator
  MaxClique::VESepInfo veSepInfo;
  // a pointer to the ve sep clique that contains the actual probabily table.
  // This is computed in prepareForUnrolling(). Once it is computed,
  // we are not allowed to make any additional copies of this object,
  // or otherwise this pointer might get deleted twice.
  ConditionalSeparatorTable* veSepClique;
  ///////////////////////////////////////////////
  // VE separator files information.
  ///////////////////////////////////////////////
  // Command line: recompute the VE separator tables and save to disk in all cases.
  static bool recomputeVESeparatorTables;
  // File name to read/write VE separator table.
  static const char* veSeparatorFileName;
  // set to true if we are (re-)generating the VE tables,
  // or set to false if we are just reading them in from disk.
  static bool generatingVESeparatorTables;
  // actual file to get ve sep stuff.
  static FILE* veSeparatorFile;
  // The log (base 10) upper limit on a VE sep variable cardinality
  // product. I.e., if the number of parents that need to be iterated
  // over to produce the VE sep table has a prod. of cardinalties
  // greater than this (in log base 10), we don't use this VE sep.
  static float veSeparatorLogProdCardLimit;



  // A boolean flag that the inference code uses to determine if it
  // should skip this separator. This is used when a P partition is
  // empty and the C partition needs to skip its incomming
  // interface separator (which is empty and if it wasn't
  // skipped would lead to a zero probability).
  bool skipMe;

  // copy constructor 
  SeparatorClique(const SeparatorClique& sep)
    : veSeparator(sep.veSeparator)
  { 
    // this constructor only copies the non-filled out information
    // (nodes and veSep status and information) since the other stuff
    // isn't needed for this type of constructor. Note that if other
    // code changes, we might need to add more here.
    nodes = sep.nodes; 
    veSepInfo = sep.veSepInfo;
    // make sure this is NULL as if we copy this in, it will get
    // deleted twice.
    assert ( sep.veSepClique == NULL );;
    veSepClique = NULL;
    skipMe = sep.skipMe;
  }

  // constructor for VE separators.
  SeparatorClique(const MaxClique::VESepInfo& _veSepInfo)
    : veSeparator(true)
  { 
    veSepInfo = _veSepInfo;
    // need nodes to reflect union, to sort, etc.
    for (unsigned i=0;i<veSepInfo.parents.size();i++) {
      nodes.insert(veSepInfo.parents[i]);
    }
    // child is guaranteed not to be NULL.
    assert ( veSepInfo.child != NULL );
    if (veSepInfo.grandChild == NULL) {
      // then this is a PC case.
      // Here, we do not insert the child since
      // it is not officially part of the separator.
    } else {
      // the PCG case. Here we insert the child, since the separator
      // is relative to the grandchild, so we iterate over both the
      // parents and the child. 
      nodes.insert(veSepInfo.child);
      // We do not insert the grandchild since it is not officially
      // part of the separator.
    }
    veSepClique = NULL;
    skipMe = false;
  }

  // construct a separator between two cliques
  SeparatorClique(MaxClique& c1, MaxClique& c2);

#if 0
  // not used any longer. Keep here in case we 
  // need it again in future.
  SeparatorClique(SeparatorClique& from_sep,
		  vector <RV*>& newRvs,
		  map < RVInfo::rvParent, unsigned >& ppf,
		  const unsigned int frameDelta = 0);
#endif
  
  // Create an empty separator, used for re-construction later.
  // Hopefully STL doesn't allocate anything when using default
  // constructor, as if it does, it will be lost. Note that we do not
  // repeatedly construct one of these objects, so while we might have
  // a bit of lost memory as a result of this, it won't constitute an
  // ever-growing memory leak.
  SeparatorClique() : veSeparator(false),veSepClique(NULL) {}

  ~SeparatorClique();

  // compute the weight (log10 state space) of this separator clique.
  float weight(const bool useDeterminism = true) const { 
    return MaxClique::computeWeight(nodes,NULL,useDeterminism); 
  }

  // prepare the last set of data structures so that clones of this
  // can be unrolled and inference can occur.
  void prepareForUnrolling();

  // print out everything in this clique to a file.
  void printAllJTInfo(FILE* f);

  // memory reporting
  void reportMemoryUsageTo(FILE *f);

  // Manages and memorizes the size and space requests made by all
  // corresponding ConditionalSeparatorTables regarding the size of
  // 'separatorValues' array. I.e., the next time a
  // ConditionalSeparatorTable asks for an initial amount of memory,
  // we'll be able to give it something closer to what it previously
  // asked for rather than something too small.
  SpaceManager separatorValueSpaceManager;

  // TODO: Manages and memorizes space for the remainder portion of a
  // separator. Note that this is used quite differently as the
  // 'separatorValueSpaceManager' variable, as a potentially different
  // sized remainder exists for all accumulated intersection
  // size. This therefore only keeps track of the maximum size.
  SpaceManager remainderValueSpaceManager; 

  // used to clear out hash table memory between segments
  void clearSeparatorValueCache(bool force=false);

  void clearInferenceMemory() {
    accSepValHashSet.clear();
    remSepValHashSet.clear();
    // could change to makeEmpty if running on a static graph only 
    if (accPacker.packedLen() > ISC_NWWOH_AI)
      accValueHolder.prepare();
    if (remPacker.packedLen() > ISC_NWWOH_RM)
      remValueHolder.prepare();
  }


  // @@@ need to take out, here for now to satisify STL call of vector.clear().
#if 0
  SeparatorClique& operator=(const SeparatorClique& f) {
    return *this;
  }
#endif

};


// A version of separatorclique that:
//   1) has no STL and uses only fast datastructures with custom
//      memory managment
//   2) keeps a pointer back to it's original clique for hash tables,
//      etc.
// Note that this is not a subclass of SeparatorClique since we do not want
// these objects to have to have all SeparatorClique's member variables.
class ConditionalSeparatorTable : public IM
{
  friend class MaxCliqueTable;
  friend class PartitionStructures;
  friend class PartitionTables;
  friend class SeparatorClique;

  // Member functions that are shared accross multiple instances of
  // this table, and thus are kept separate and passed in via
  // arguments as needed, rather than wasting storage to keep these
  // member functions around. We manage the data structure here though
  // since here is where it is defined what is and isn't needed.
  struct SharedLocalStructure {


    // the original separator clique from which this object has been cloned.
    SeparatorClique* origin;

    // Non-STL "f=fast" versions of arrays that are instantiated
    // only in the final unrolled versions of the separator cliques.
    sArray< RV*> fNodes;
    sArray< RV*> fAccumulatedIntersection; // remove
    sArray< RV*> fRemainder;

    // Direct fast access pointers to the values within the discrete
    // *HIDDEN* RVs within this sep clique.  Note, the observed discrete
    // and continuous variables are not contained here.  1) one for the
    // accumulated intersection
    sArray < DiscRVType*> accDiscreteValuePtrs;
    // 2) and one for the remainder variables.
    sArray < DiscRVType*> remDiscreteValuePtrs;


    // @@@@ NOW TODO: include a full function that takes a
    // maxClique/separator, offset, and new rv set and fills in all of
    // the member slots here (like what hte old constructor of the
    // containing object used to do). (the consructure soudl do it, no 1/15)

    set <RV*> returnRVsAsSet();

    SharedLocalStructure(SeparatorClique& _origin,
			 vector <RV*>& newRvs,
			 map < RVInfo::rvParent, unsigned >& ppf,
			 const unsigned int frameDelta);
    // empty constructor
    SharedLocalStructure() : origin(NULL) {}

  };

  // Two indices to get at the veterbi values for current
  // separator. In other words, these indices give the separator entry
  // corresponding to the variable assigmnets that give the max value
  // for the clique closer to the root relative to this separator. We
  // store these indices here so that that clique need not retain its
  // values in the set of rvs, while the current separator can still
  // store the values it needs. For VE separators, we never use these
  // for their normal use, so instead we encode an id to indicate that
  // this is a VE separator (so that certain other data items do not
  // get deleted on destruction). Actualy, this can be done with
  // a single pointer directly to the appropriate entry!!
  struct ForwardPointer {
    unsigned viterbiAccIndex;
    unsigned viterbiRemIndex;
    ForwardPointer() 
    {
      viterbiAccIndex = viterbiRemIndex = 0;
    }
    void setToVeSeparatorId() {
      viterbiAccIndex = viterbiRemIndex = (unsigned)~0x0;
    }
    bool veSeparator() {
      // we use special values here to indicate if the current
      // separator is a VE (virtual evidence) separator which is
      // treated very differently in the inference code.
      return (viterbiAccIndex == (unsigned)~0x0 && viterbiRemIndex == (unsigned)~0x0);
    }
  };


  class RemainderValue {
  public:
    // shared space.
    union {
      // When a packed clique value is > ISC_NWWOH_RM words
      // (unsigned), then we keep a pointer to the packed clique value
      // (list of variables) where the actuall clique value is
      // obtained. This points to a data structure maintained by
      // origin.
      unsigned *ptr;
      // When a packed clique value is only ISC_NWWOH_RM words
      // (unsigned) or less we don't bother with a hash table and just
      // store the packed clique value right here, thereby saving
      // needing to look things up in a hash table and also saving
      // memory.
      // unsigned val[ISC_NWWOH_RM];
      unsigned val[((ISC_NWWOH_RM>1)?ISC_NWWOH_RM:1)];
    };
    // forward probability
    logpr p;

    // make a union to share space between doing
    // backward inference and backward viterbi pass
    union {
      char _bpo[sizeof(logpr)];
      // TODO: make this an unsigned* to be able
      // todo n-best.
      unsigned backPointer;
      // pointer to array of nbest backpointers.
      // use sign bit for 'used' status.
      // unsigned *nBestList;
    };

    // probability for distribute evidence pass
    // easy access to different fields
    inline logpr& bp() { return (*((logpr*)(&_bpo[0]))); }

    RemainderValue() { 
      bp().set_to_zero(); 
    }
  };


  // Accumulated intersection separator values.
  class AISeparatorValue {
  public:
    // shared space.
    union {
      // When a packed clique value is > ISC_NWWOH_AI words
      // (unsigned vals), then we keep a pointer to the packed clique value
      // (list of variables) where the actuall clique value is
      // obtained. This points to a data structure maintained by
      // origin.
      unsigned *ptr;
      // When a packed clique value is only ISC_NWWOH_AI words
      // (unsigned vals) or less we don't bother with a hash table and
      // just store the packed clique value right here, thereby saving
      // needing to look things up in a hash table and also saving
      // memory.
      // unsigned val[ISC_NWWOH_AI];
      unsigned val[((ISC_NWWOH_AI>1)?ISC_NWWOH_AI:1)];
    };
    // Array of remainder separator values for this accumulated
    // intersection value. 
    cArray <RemainderValue> remValues;
    // number of currently used remainder values
    unsigned numRemValuesUsed;

    // Hash table into remainder. This is used during clique iteration
    // to project down into the outgoing separator in a CE stage.
    // When we iterate over a clique to project into the outgoing
    // separator during CE, we need to index (hash) into both the AI,
    // and for the AI found, see if the remainder exists. This hash
    // table is used for that purpose.  Note, by default, these hash
    // tables are constructed empty and invalid, and they need to be
    // explicitly re-constructed.
    VHashMapUnsignedUnsignedKeyUpdatable iRemHashMap;

    // ensure that we start with nothing inserted.
    AISeparatorValue() { 
      numRemValuesUsed = 0;
      // fprintf(stderr,"in aisv init\n");
    }
  };


  // Forward pointer for viterbi, and for encoding VE separator identity status.
  // TODO: figure out a way to avoid needing to store these here.
  // TODO: make an array for N-best decoding.
  ForwardPointer forwPointer;
  bool veSeparator() { return forwPointer.veSeparator(); }
  void setToVeSeparatorId() { forwPointer.setToVeSeparatorId(); }

  // The collection of AI separator values.
  cArray< AISeparatorValue >* separatorValues;
  // number of currently used clique values
  unsigned numSeparatorValuesUsed;
  // The hash table from accum intersection rv values to integer index
  // into separatorValues. This hash table is used to quickly lookup
  // the values in this separator that are common (intersected) with
  // the accumulated union of all previous separators that we have
  // already iterated over. During seprator driven clique iterations,
  // if this hash table is looked up and we get a hit, it means that
  // there is a common value between this separator and all previous
  // separators that we have so far seen, and we must continue
  // iteration using the remaining rvs in this separator.
  VHashMapUnsignedUnsignedKeyUpdatable* iAccHashMap;

public:


  // WARNING: constructor hack to create a very broken object with
  // non-functional reference objects (in order to create an array of
  // these objects and then initialize them later with appropriate
  // references). Do not use until after proper re-constructor.
  ConditionalSeparatorTable()  
  { iAccHashMap = NULL; separatorValues = NULL; }
  // normal (or re-)constructor.
  ConditionalSeparatorTable(SeparatorClique& origin);
  void init(SeparatorClique& origin);

  ///////////////////////////////////////////////////////////////////////
  // version to create a VE separators that lives in a SeparatorClique
  ConditionalSeparatorTable(SeparatorClique& origin,
			    ConditionalSeparatorTable::SharedLocalStructure& sepSharedStructure);

  // destructor
  ~ConditionalSeparatorTable() 
  {
    clearInferenceMemory();
  }

  void clearInferenceMemory() {
    // only delete when not a VE separator, since when it is these
    // guys are shared accross multiple ConditionalSeparatorTables.
    // Note: the 'mother' VE ConditionalSeparatorTable is actually a
    // placeholder VE separator, but it is deleted only when the
    // containing SeparatorClique is deleted.
    if (veSeparator())
      return;
    delete iAccHashMap;
    iAccHashMap = NULL;
    delete separatorValues;
    separatorValues = NULL;
  }


  // insert current value of RVs into separator
  void insert(SeparatorClique& origin,
	      ConditionalSeparatorTable::SharedLocalStructure& sepSharedStructure);

  // separator based pruning
  void ceSeparatorPrune(SeparatorClique& origin);

  // memory reporting
  void reportMemoryUsageTo(SeparatorClique&,FILE *f);

};



// The table entries for a maxclique that:
//   1) has no STL and uses only fast datastructures with custom memory managment
//   2) keeps a pointer back to it's original clique for hash tables, etc.
// 
// Note that this is not a subclass of MaxClique since we want this
// object to be as absolutely small as possible to save memory on very
// long segments (e.g., genomes), and we definitely do not want these
// objects to have to have all MaxClique's member variables which
// would be highly redundant.
class MaxCliqueTable  : public IM
{
  friend struct CliqueValueDescendingProbCompare;
  friend class PartitionStructures;
  friend class PartitionTables;

  // integer value to keep track of indenting when running in trace
  // mode. TODO: reentrant issues.
  static int traceIndent;
  static const unsigned spi;


  // Member functions that are shared accross multiple instances of
  // this table, and thus are kept separate and passed in via
  // arguments as needed, rather than wasting storage to keep these
  // member functions around. We manage the data structure here though
  // since here is where it is defined what is and isn't needed.
  struct SharedLocalStructure {
    // The original maxclique from which this object has been created
    // and where we get access to some data structures that are common
    // to all MaxCliqueTable's of this particular base clique.
    MaxClique* origin;

    // Non-STL "fast" versions of arrays that exist in the final
    // unrolled versions of the cliques. These give rapid access to the
    // random variables involved in this clique and used in the unrolled
    // random variables during inference (i.e., the variables that are
    // shared amongst cliques and their separators both within and
    // accross cliques and partitions).
    sArray< RV*> fNodes;
    sArray< RV*> fSortedAssignedNodes;
    sArray< RV*> fUnassignedIteratedNodes;
    sArray< RV*> fDeterminableNodes;
    // Direct pointers to the values within the discrete hidden RVs
    // within this clique.  Note, the observed discrete and continuous
    // variables are not contained here.
    sArray < DiscRVType* > discreteValuePtrs;

    // pointers to RVs having the max frame value and the min frame value
    RV* rv_w_max_frame_num;
    RV* rv_w_min_frame_num;


    // initialize the above members.
    SharedLocalStructure(MaxClique& _origin,
			 vector <RV*>& newRvs,
			 map < RVInfo::rvParent, unsigned >& ppf,
			 const unsigned int frameDelta);
    // empty constructor
    SharedLocalStructure() : origin(NULL) {}

    set <RV*> returnRVsAsSet();

    // return as as set the random variables and any of their observed
    // parents as a set. The main use of this routine is to return the
    // random variables that might get their frame number changed.  The
    // reason for also returning the observed parents is that when a
    // frame number changes, the observed values might also need to
    // change, and observed parents might live in a different clique or
    // even a different partition (see the variable
    // disconnectChildrenOfObservedParents in RV.h).
    set <RV*> returnRVsAndTheirObservedParentsAsSet();

  };


  // Data structure holding a clique value, i.e., a set of random
  // variable values and a probability for that set of values.
  // -
  // TODO: potentially create a global table of these clique values and
  // in this clique, store rather than actual clique values, store
  // instead integer indices into this global table. When clique
  // values get pruned away, they can be restored into global table.
  // Advantages of such a scheme:
  //    a) pruned clique values get removed right away, reclaiming
  //       storage for additional clique values later in inference
  //    b) in distribute evidence (backward pass), might be useful
  //       and faster not to have to keep checking forward pass pruning
  //       threshold (i.e., pruned away clique values will need
  //       to be ignored in some way or another)
  //    c) wasted storage might be less since unused array
  //       in each clique holds only one word (index) rather
  //       than entire contents of a clique value. (but
  //       global clique value pool needs to hold unused values as well
  //       but how many?? not clear.
  // 
  // Dis-advantages of such a scheme:
  //    a) uses more storage per clique value (an integer index +
  //     the entries in the clique value
  //    b) still more indirection toget to data we need.
  //    c) more bookeeping to keep
  // 
  // Note: in any event, this class needs to be as small as possible!!!
  class CliqueValue {
  public:

    // shared space.
    union {
      // When a packed clique value is > ISC_NWWOH words (unsigned),
      // then we keep a pointer to the packed clique value (list of
      // variables) where the actuall clique value is obtained. This
      // points to a data structure maintained by origin.
      union {
	// Keep both a ptr and an unsigned long version of the ptr.
	// The ptr is used when pointing to a shared global value pool
	// and the ival is used to point to an index entry in a
	// temporary clique value pool which is used to store values
	// as the clique table is being created but before pruning
	// happens. Values that are in the global clique value pool
	// are never removed, but this option allows the insertion of
	// values only that have not been pruned away, thereby saving
	// much memory.
	unsigned *ptr;
	unsigned long ival;
      };
      // When a packed clique value is only ISC_NWWOH words (unsigned)
      // or less we don't bother with a hash table and just store the
      // packed clique value right here, thereby saving needing to
      // look things up in a hash table and also saving memory.
      // unsigned val[IMC_NWWOH];
      unsigned val[((IMC_NWWOH>1)?IMC_NWWOH:1)];
    };

    // The probability p. Note that we could keep a collect evidence
    // and distribute evidence probability here (and thereby avoid
    // doing the Hugin-style divide on the distribute evidence stage)
    // but we only keep one value 1) to save space, as adding an extra
    // probability will increase storage requirements (especially if a
    // logpr is a 64-bit fp number), and 2) since everything is done
    // in log arithmetic, a divide is really a floating point
    // subtraction which is cheap.
    logpr p;

  };

  // the collection of clique values in the clique table for this
  // clique.
  sArray< CliqueValue > cliqueValues;
  // Number of currently used clique values. Note that this value
  // might be different than the size contained in the sArray since
  // the sArray indicates the amount allocated and/or expanded, and
  // differences in allocation vs. total number that have been
  // generated, and also pruning, will mean that these values can be
  // different.
  unsigned numCliqueValuesUsed;


  // extra storage to store special clique values for things like
  // n-best, etc.
  union {
    unsigned back_max_cvn;
    // for n-best, a length n array, should create a n-chunk array
    // class that allocates these in fixed units using operator new.
    // unsigned* back_max_cvn_arr;
  };

#ifdef TRACK_NUM_CLIQUE_VALS_SHARED
  // Number of times that the clique value was shared from a time before
  unsigned numCliqueValuesShared;
#endif


public:

  // WARNING: constructor hack to create a very broken object with
  // non-functional reference objects (in order to create an array of
  // these objects and then initialize them later with appropriate
  // references). Do not use until after proper re-constructor.
  MaxCliqueTable() {}
  // normal constructor (i.e., re-constructor).
  MaxCliqueTable(MaxClique& origin);
  void init(MaxClique& origin);

  // destructor
  ~MaxCliqueTable() {}

  ///////////////////////////////////////////////////////////////////////////////
  // collect evidence functions.

  void ceGatherFromIncommingSeparators(MaxCliqueTable::SharedLocalStructure& sharedStructure,
				       ConditionalSeparatorTable* separatorTableArray,
				       ConditionalSeparatorTable::SharedLocalStructure* sepSharedStructureArray);
				       

  // special case when the clique is fully observed.
  void ceGatherFromIncommingSeparatorsCliqueObserved(MaxCliqueTable::SharedLocalStructure& sharedStructure,
						     ConditionalSeparatorTable* separatorCliqueArray,
						     ConditionalSeparatorTable::SharedLocalStructure* sepSharedStructureArray,
						     logpr& maxCEValue);


  ////////////////////////////////////////////////////////////////////////////
  // iterate unassigned nodes
  void inline ceIterateUnassignedIteratedNodes(MaxCliqueTable::SharedLocalStructure& sharedStructure,
					       logpr cliqueBeamThresholdEstimate,
					       logpr& maxCEValue,
					       const unsigned nodeNumber,
					       const logpr p)
  {
    if (nodeNumber == sharedStructure.fUnassignedIteratedNodes.size()) {
      ceIterateAssignedNodes(sharedStructure,cliqueBeamThresholdEstimate,maxCEValue,0,p);
      return;
    }
    ceIterateUnassignedIteratedNodesRecurse(sharedStructure,
					    cliqueBeamThresholdEstimate,
					    maxCEValue,
					    nodeNumber,
					    p);
  }
  void ceIterateUnassignedIteratedNodesRecurse(MaxCliqueTable::SharedLocalStructure& sharedStructure,
					       logpr cliqueBeamThresholdEstimate,
					       logpr& maxCEValue,
					       const unsigned nodeNumber,
					       const logpr p);

  ///////////////////////////////////////
  // iterate separators

  void inline ceIterateSeparators(MaxCliqueTable::SharedLocalStructure& sharedStructure,
				  ConditionalSeparatorTable* separatorTableArray,
				  ConditionalSeparatorTable::SharedLocalStructure* sepSharedStructureArray,
				  logpr cliqueBeamThresholdEstimate,
				  logpr& maxCEValue,
				  const unsigned sepNumber,
				  const logpr p) 
  {
    if (sepNumber == sharedStructure.origin->ceReceiveSeparators.size()) {
      // move on to the iterated nodes.
      ceIterateUnassignedIteratedNodes(sharedStructure,
				       cliqueBeamThresholdEstimate,
				       maxCEValue,
				       0,p);
      return;
    }
    ceIterateSeparatorsRecurse(sharedStructure,separatorTableArray,sepSharedStructureArray,
			       cliqueBeamThresholdEstimate,
			       maxCEValue,
			       sepNumber,p);
  }
  void ceIterateSeparatorsRecurse(MaxCliqueTable::SharedLocalStructure& sharedStructure,
				  ConditionalSeparatorTable* separatorTableArray,
				  ConditionalSeparatorTable::SharedLocalStructure* sepSharedStructureArray,
				  logpr cliqueBeamThresholdEstimate,
				  logpr& maxCEValue,
				  const unsigned sepNumber,
				  const logpr p);

  ///////////////////////////////////////
  // iterate assigned nodes.


  void ceIterateAssignedNodesRecurse(MaxCliqueTable::SharedLocalStructure& sharedStructure,
				     logpr cliqueBeamThresholdEstimate,
				     logpr& maxCEValue,
				     const unsigned nodeNumber,
				     const logpr p);
  void ceIterateAssignedNodesNoRecurse(MaxCliqueTable::SharedLocalStructure& sharedStructure,
				       logpr cliqueBeamThresholdEstimate,
				       logpr& maxCEValue,
				       const logpr p);

  void inline ceIterateAssignedNodes(MaxCliqueTable::SharedLocalStructure& sharedStructure,
				     logpr cliqueBeamThresholdEstimate,
				     logpr& maxCEValue,
				     const unsigned nodeNumber,
				     const logpr p)
  {
    // apply factors so far from separators or unassigned nodes.

    // if (true || sharedStructure.fSortedAssignedNodes.size() == 0 || message(High)) {
    if (sharedStructure.fSortedAssignedNodes.size() == 0 || message(High)) {
      // let recursive version handle degenerate or message full case
      ceIterateAssignedNodesRecurse(sharedStructure,
				    cliqueBeamThresholdEstimate,
				    maxCEValue,
				    0,
				    p);
    } else {
      // ceIterateAssignedNodesRecurse(part,0,p);
      ceIterateAssignedNodesNoRecurse(sharedStructure,
				      cliqueBeamThresholdEstimate,
				      maxCEValue,
				      p);
    }
  }


  /////////////////////////////////////////
  // a version that automatically selects which separator to use within partition from the clique.
  void ceSendToOutgoingSeparator(MaxCliqueTable::SharedLocalStructure& sharedStructure,
				 ConditionalSeparatorTable* separatorTableArray,
				 ConditionalSeparatorTable::SharedLocalStructure* sepSharedStructureArray);
  // a version with an explicit outgoing separator.
  void ceSendToOutgoingSeparator(MaxCliqueTable::SharedLocalStructure& sharedStructure,
				 ConditionalSeparatorTable& sep,
				 ConditionalSeparatorTable::SharedLocalStructure&);




  /////////////////////////////////////////
  // memory clearing.
  void clearInferenceMemory() {
    // clear out all memory used by this inference clique.
    cliqueValues.clear();
  }
  void clearCliqueAndIncommingSeparatorMemory(MaxCliqueTable::SharedLocalStructure& sharedStructure,
					      ConditionalSeparatorTable*,
					      ConditionalSeparatorTable::SharedLocalStructure*);


  /////////////////////////////////////////
  // Pruning
  /////////////////////////////////////////

  void ceDoAllPruning(MaxClique& origin,logpr maxCEValue);
  void ceCliqueBeamPrune(MaxClique& origin,logpr maxCEValue);
  unsigned ceCliqueStatePrune(const unsigned k,
			      CliqueValue*,
			      const unsigned);
  unsigned ceCliqueMassPrune(const double removeFraction,
			     const double exponentiate,
			     const double furtherBeam,
			     const unsigned minSize,
			     CliqueValue*,
			     const unsigned);
  void ceCliqueDiversityPrune(MaxClique& origin,const unsigned numClusters);

  void ceCliqueUniformSamplePrunedCliquePortion(MaxClique& origin,
						const unsigned origNumCliqueValuesUsed);

  void ceDoCliqueScoreNormalization(MaxCliqueTable::SharedLocalStructure& sharedStructure);


#ifdef USE_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL
  void insertLocalCliqueValuesIntoSharedPool(MaxClique& origin);
#endif


  /////////////////////////////////////////
  // distribute evidence functions.
  /////////////////////////////////////////

  void deScatterToOutgoingSeparators(MaxCliqueTable::SharedLocalStructure& sharedStructure,
				     ConditionalSeparatorTable* separatorTableArray,
				     ConditionalSeparatorTable::SharedLocalStructure* sepSharedStructureArray);


  void deScatterToOutgoingSeparatorsViterbi(MaxCliqueTable::SharedLocalStructure& sharedStructure,
					    ConditionalSeparatorTable* separatorTableArray,
					    ConditionalSeparatorTable::SharedLocalStructure* sepSharedStructureArray);

  // a version that automatically selects which separator to use within partition from the clique.
  void deReceiveFromIncommingSeparator(MaxCliqueTable::SharedLocalStructure& sharedStructure,
				       ConditionalSeparatorTable* separatorTableArray,
				       ConditionalSeparatorTable::SharedLocalStructure* sepSharedStructureArray);
  // a version with an explicit incomming separator
  void deReceiveFromIncommingSeparator(MaxCliqueTable::SharedLocalStructure& sharedStructure,
				       ConditionalSeparatorTable&,
				       ConditionalSeparatorTable::SharedLocalStructure&);
  // a version specific to viterbi decoding.
  void deReceiveFromIncommingSeparatorViterbi(MaxCliqueTable::SharedLocalStructure& sharedStructure,
					      ConditionalSeparatorTable&,
					      ConditionalSeparatorTable::SharedLocalStructure&);


  ////////////////////////////////////////////////////////////////////////////////
  // MISC.

  // sum up the probabilities in the current clique table and return
  // their value.
  logpr sumProbabilities();
  logpr sumExponentiatedProbabilities(double exponent,
				      CliqueValue* curCliqueVals,
				      const unsigned curNumCliqueValuesUsed);

  // compute the clique entropy
  double cliqueEntropy();
  // compute a form of clique diversity
  double cliqueDiversity(MaxClique& origin);

  // memory reporting.
  void reportMemoryUsageTo(MaxClique& origin,FILE *f);


  // compute the max probability and return its value, and also
  // optionally sets the rvs to its max value.
  logpr maxProbability(MaxCliqueTable::SharedLocalStructure&,
		       bool setCliqueToMaxValue = true);
  // a faster version that just operates on the table.
  logpr maxProb();

  // print all clique values and prob to given file.
  void printCliqueEntries(MaxCliqueTable::SharedLocalStructure&,
			  FILE*f,const char*str=NULL,
			  const bool normalize = false,
			  const bool justPrintEntropy = false);
  
  // EM accumulation support.
  void emIncrement(MaxCliqueTable::SharedLocalStructure&,
		   const logpr probE, 
		   const bool localCliqueNormalization = false,
		   const double emTrainingBeam = -LZERO);


};




// A factor clique is a (not-necessarily max) clique that implements
// some form of directed "factor" or "constraint" within a regular max
// clique. This coresponds to the 'factor' construct in the .str file.
// A factor clique's nodes will, like a separator, necessarily be a
// subset of some MaxClique, but a factor is implemented quite
// differently. Factors mixed with DAG cpts allow the specificatio of
// true hybrid directed/undirected graphical models in GMTK.
class FactorClique : public IM
{
  friend class FileParser;
  friend class GraphicalModel;
  friend class GMTemplate;
  friend class ConditionalSeparatorTable;
  friend class JunctionTree;

public:


  // information about the factor corresponding
  // to this factor
  FactorInfo* factorInfo;

  // the set of nodes that are involved in the factor/constraint.
  set<RV*> nodes;
  // vector of the same nodes, in order that they
  // appear in the .str file, and the order used to index into
  // the factor functions.
  vector<RV*> orderedNodes;

  // copy constructor 
  FactorClique(const FactorClique& factor)
  { 
    factorInfo = factor.factorInfo;
    nodes = factor.nodes; 
    orderedNodes = factor.orderedNodes;
  }
  FactorClique(FactorInfo& factorInfo,
	       vector <RV*>& unrolled_rvs,
	       map < RVInfo::rvParent, unsigned > ppf,
	       const unsigned offset);


  ~FactorClique() {}

  // prepare the last set of data structures so that clones of this
  // can be unrolled and inference can occur.
  void prepareForUnrolling();

  // print out everything in this clique to a file.
  void printAllJTInfo(FILE* f);

};


class InferenceFactorClique : public IM
{

  // Non-STL "f=fast" versions of arrays that are instantiated
  // only in the final unrolled versions of the separator cliques.
  sArray< RV*> fOrderedNodes;

  // the original factor clique from which this object has been cloned.
  FactorClique& origin;


public:

  // WARNING: constructor hack to create a very broken object with
  // non-functional reference objects (in order to create an array of
  // these objects and then initialize them later with appropriate
  // references). Do not use until after proper re-constructor.
  InferenceFactorClique() : origin(*((FactorClique*)NULL)) {}
  // normal (or re-)constructor
  InferenceFactorClique(FactorClique& _origin,
			vector <RV*>& newRvs,
			map < RVInfo::rvParent, unsigned >& ppf,
			const unsigned int frameDelta);
  // version for VE separators
  InferenceFactorClique(FactorClique& _origin);

};





#endif
