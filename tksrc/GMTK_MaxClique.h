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

#ifndef GMTK_MAXCLIQUE_H
#define GMTK_MAXCLIQUE_H

#include "general.h"
#include "vhash_set.h"
#include "vhash_map.h"
#include "logp.h"
#include "cArray.h"
#include "sArray.h"
#include "debug.h"


#include "GMTK_RandomVariable.h"
#include "GMTK_PackCliqueValue.h"

#include <stdio.h>
#include <stdlib.h>

#include <vector>
#include <string>
#include <map>

class JT_InferencePartition;
class SeperatorClique;
class InferenceSeparatorClique;
class MaxClique;
class InferenceMaxClique;

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
  // @@@ sould be const
  float growthFactor;

  // The size of the initial allocation unit. A chunk is an
  // array of packed clique values. When we allocate the n'th
  // chunk, we allocate k(n) new packed clique values, 
  // where k(n) = allocationUnitChunkSize*growthFactor^(n-1)
  // @@@ sould be const
  unsigned allocationUnitChunkSize;

  // Current capacity, total number of packed clique values that this
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

  ~CliqueValueHolder() { makeEmpty(); }

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


};


// a (someday to be) special vhash class
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
  unsigned*& tableKey(const unsigned i) { return table[i].key; }
  unsigned& tableItem(const unsigned i) { return table[i].item; }
  bool tableEmpty(const unsigned i) { return empty(table[i]); }
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
// InferenceMaxClique Number Words WithOut a Hash (IMC_NWWOH): Namely,
// the number of words that can be stored directly as a packed clique
// value before we resort to using a shared hash table for all
// instances of this origin clique.
#define IMC_NWWOH (1)
// InferenceSeparatorClique Number Words WithOut a Hash (ISC_NWWOH):
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
  friend class SeperatorClique;

public:


  // bool which if true means we globally do separator driven (rather
  // than clique driven) inference. Separator driven is good when
  // pruning is turned on. Clique driven is good when you want
  // inference cost to directly correspond to normal weight.
  static bool ceSeparatorDrivenInference;

  // beam width for clique-based beam pruning.
  static double cliqueBeam;


  // @@@ need to take out, here for now to satisify STL call of vector.clear().
#if 0
  MaxClique& operator=(const MaxClique& f) {
    return *this;
  }
#endif

  // the set of nodes which form a max clique, in arbitrary order.
  set<RandomVariable*> nodes;

  // weight of this clique, keep pre-computed here.
  float cliqueWeight;

  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////

  // basic constructor with a set of nodes
  MaxClique(set<RandomVariable*> arg) {
    nodes = arg;
  }

  // Clone constructor from another MaxClique, but that uses a new set
  // of random variables, and adjusts the frame of each new set of
  // random variable with offset
  MaxClique(MaxClique& from_clique,
	    vector <RandomVariable*>& newRvs,
	    map < RVInfo::rvParent, unsigned >& ppf,
	    const unsigned int frameDelta = 0);


  ~MaxClique() {}


  // Complete the max clique by adjusting the random variables
  // themselves.
  void makeComplete() { makeComplete(nodes); }
  // Static version of variable set completion, to complete set of
  // random variables passed in.
  static void makeComplete(const set<RandomVariable*> &rvs);



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
  static float computeWeight(const set<RandomVariable*>& nodes,
			     const RandomVariable* node = NULL,
			     const bool useDeterminism = true);
  static float computeWeightWithExclusion(const set<RandomVariable*>& nodes,
					  const set<RandomVariable*>& unassignedIteratedNodes,
					  const set<RandomVariable*>& unionSepNodes,
					  const bool useDeterminism = true);
  static float computeWeightInJunctionTree(const set<RandomVariable*>& nodes,
					   const set<RandomVariable*>& assignedNodes,
					   const set<RandomVariable*>& cumulativeAssignedNodes,
					   const set<RandomVariable*>& unassignedIteratedNodes,
					   const set<RandomVariable*>& cumulativeUnassignedIteratedNodes,
					   const set<RandomVariable*>& separatorNodes,
					   const set<RandomVariable*>& unassignedInPartition,
					   const bool upperBound,
					   const bool useDeterminism);

  // compute the weight (log10 state space) of this clique.
  float weight(const bool useDeterminism = true) const { 
    return computeWeight(nodes,NULL,useDeterminism); 
  }
  float weightInJunctionTree(const set<RandomVariable*>& unassignedInPartition,
			     const bool upperBound = false,
			     const bool useDeterminism = true) const { 
    return computeWeightInJunctionTree(nodes,
				       assignedNodes,
				       cumulativeAssignedNodes,
				       unassignedIteratedNodes,
				       cumulativeUnassignedIteratedNodes,
				       unionIncommingCESeps,
				       unassignedInPartition,
				       upperBound,
				       useDeterminism);
  }



  // print just the clique nodes
  void printCliqueNodes(FILE* f);

  
  //////////////////////////////////////////////
  // TODO: figure out a way so that the member variables below exist only in
  // a subclass since all of the below is not needed for basic
  // triangulation.

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
  set<RandomVariable*> assignedNodes;

  // USED ONLY IN JUNCTION TREE INFERENCE
  // A topologically sorted vector version of assignedNodes.
  // The topological sort is relative only to the variables
  // in this clique. The topological sort is used during
  // CE clique creation, so that CPT iteration can be used when available.
  // Computed in JunctionTree::assignRVsToCliques().
  // Used to:
  //   1) compute assigned nodes dispositions (in appropriate order)
  //   2) compute fSortedAssignedNodes in an inference clique
  vector<RandomVariable*> sortedAssignedNodes;

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
  set<RandomVariable*> assignedProbNodes;


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
  set<RandomVariable*> cumulativeAssignedNodes;

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
  set<RandomVariable*> cumulativeAssignedProbNodes;


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
  set<RandomVariable*> unionIncommingCESeps;

  // USED ONLY IN JUNCTION TREE INFERENCE
  // The set of nodes that are not assigned to this clique but that
  // also are *NOT* iterated over by the incomming
  // separators. Therefore, these nodes in the clique must be iterated
  // over [0,card-1] unfortunately.  The hope is that assignment and setup can be
  // done so there are very few or zero of such nodes. Note however
  // that once a clique iterates over such a node, any later cliques
  // (closer to the JT root) will not have to re-do this as these
  // assignments here will be represented in this cliques out-going
  // sepset.
  // Computed in JunctionTree::computeSeparatorIterationOrders()
  // Used to:
  //  1) compute JT weight in clique
  //  2) compute fUnassignedIteratedNodes in inference clique (which says
  //     which nodes in a clique are unassigend so iterated [0,card-1]
  //  3) in CE, choose which routine to call first.
  set<RandomVariable*> unassignedIteratedNodes;


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
  set<RandomVariable*> unassignedNodes;


  // USED ONLY IN JUNCTION TREE INFERENCE 
  // This enum defines the different things that could
  // happen to an assigned node in a clique during CE assigned node iteration.
  enum AssignedNodeDisposition {
    // In the discussion below, say v is an assignedNode (meaning it
    // exists with its parents in the clique).
    // We have 6 cases, each case has 1 or more name.

    //
    // 1)
    //   v is not in sepNodes, and v is a prob node as well, v either
    //   sparse or dense:  meaning we multiply this clique potential by
    //   prob(v|parents(v)). We use CPT iteration to iterate over v
    //   and compute and apply probabilities.
    AN_NOTSEP_PROB_SPARSEDENSE=0,
    AN_CPT_ITERATION_COMPUTE_AND_APPLY_PROB=0,

    //
    // 2) 
    //   v is not in sepNodes, and v is not a prob node, and V is
    //   sparse. We still use CPT iteration to iterate over v given
    //   parents, but do not multiply clique potential by
    //   probabilties. By using CPT iteration, we automatically remove any zeros
    //   now.
    AN_NOTSEP_NOTPROB_SPARSE=1,
    AN_CPT_ITERATION_COMPUTE_PROB_REMOVE_ZEROS=1,

    //
    // 3)
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
    // 4)
    //   v is in sepNodes, and v is a probability node, v either
    //   sparse or dense. Then we compute the probability for v, apply
    //   it, and continue only if the probabilty is non zero.
    AN_SEP_PROB_SPARSEDENSE=3,
    AN_COMPUTE_AND_APPLY_PROB=3,

    //
    // 5) 
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
    // 6)
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
  // computed in MaxClique::computeAssignedNodesToIterate()
  // Used to:
  //   1) determine way in which fSortedAssignedNodes in inf clique should be iterated.
  sArray < AssignedNodeDisposition > dispositionSortedAssignedNodes;


  // TODO: TRY TO REMOVE THIS AS IT IS REDUNDANT!!
  // USED ONLY IN JUNCTION TREE INFERENCE
  // the preceding cumulative set of unassigned and iterated nodes
  // relative to the root in the JT for the current partition and
  // clique. This is used to determine which of the assigned nodes
  // need actually be iterated within a clique, and which are already
  // set by the separator iterations. Note that
  // precedingUnassignedIteratedNodes does *NOT* include the
  // unassignedIteratedNodes in the current clique (thus the name
  // preceding).  
  // Computed in JunctionTree::getCumulativeUnassignedIteratedNodes()
  // Used to:
  //    1) comptue jt weight
  // Note: should be replaced by unionIncommingCESeps, perhaps remove.
  set<RandomVariable*> cumulativeUnassignedIteratedNodes;


  // USED ONLY IN JUNCTION TREE INFERENCE
  // These are the nodes that are hidden (i.e., they are
  // the non-continous hidden variables). These are
  // the nodes whose values will be hashed.
  // Computed in MaxClique::prepareForUnrolling()
  // Used to:
  //   1) set the size of the clique packer
  //   2) set the size of discreteValuePtrs in inference clique
  //      to get quick access to the set of RV vals.
  vector<RandomVariable*> hiddenNodes;

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
  // partition.
  // Created in JunctionTree::createSeparators().
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


  // USED ONLY IN JUNCTION TREE INFERENCE An array containing the
  // partial/cumulative probabilities of partial clique values and is
  // used by the non-recursive version of the clique value iterationt
  // through clique assigned nodes. In particular, used in:
  //    void ceIterateAssignedNodesNoRecurse()
  // Should be one greater than the size of: sortedAssignedNodes.size()
  sArray< logpr > probArrayStorage;


  // USED ONLY IN JUNCTION TREE INFERENCE
  // a value that is roughly equal to the order of the state space of
  // the clique. It should be smaller than the clique state space,
  // however, since pruning will make the state space smaller.
  // This will be used for initial memory allocation units, etc.
  // by this clique object and clones of this clique object.
  // It is only usable after we have prepared for unrolling.
  // TODO: add proper allocation statistics object and use that.
  unsigned allocationUnitChunkSize;


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
  // compute which assigned nodes to iterate over or not.
  void computeAssignedNodesDispositions();

  // USED ONLY IN JUNCTION TREE INFERENCE
  // print out everything in this junction tree clique to a file.
  void printAllJTInfo(FILE* f,const unsigned indent,const set<RandomVariable*>& unassignedInPartition);

  // USED ONLY IN JUNCTION TREE INFERENCE
  // used to clear out hash table memory between segments
  void clearDataMemory() {
    if (packer.packedLen() > IMC_NWWOH) {
      valueHolder.prepare();
      cliqueValueHashSet.clear();    
    }
  }

  // USED ONLY IN JUNCTION TREE INFERENCE
  // used to compute the unassigned nodes in this clique.
  void computeUnassignedCliqueNodes();

};


// A version of maxclique that:
//   1) has no STL and uses only fast datastructures with custom memory managment
//   2) keeps a pointer back to it's original clique for hash tables, etc.
// Note that this is not a subclass of MaxClique since we do not want
// these objects to have to have all MaxClique's member variables.
class InferenceMaxClique  : public IM
{

  // the original maxclique from which this object has been cloned and
  // where we get access to some data structures that are common to
  // all InferenceMaxClique of this particular base clique.
  MaxClique& origin;

  // Non-STL "fast" versions of arrays that exist in the final
  // unrolled versions of the cliques. These give rapid access
  // to the random variables involved in this clique.
  sArray< RandomVariable*> fNodes;
  sArray< RandomVariable*> fSortedAssignedNodes;
  sArray< RandomVariable*> fUnassignedIteratedNodes;
  sArray< RandomVariable*> fUnassignedNodes;
  // Direct pointers to the values within the discrete hidden RVs
  // within this clique.  Note, the observed discrete and continuous
  // variables are not contained here.
  sArray < RandomVariable::DiscreteVariableType* > discreteValuePtrs;

  // Data structure holding a clique value, i.e., a set of random
  // variable values and a probability for that set of values.
  // -
  // TODO: potentially create a gobal table of these clique values and
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

  class CliqueValue {
  public:

    // shared space.
    union {
      // When a packed clique value is > ISC_NWWOH words (unsigned),
      // then we keep a pointer to the packed clique value (list of
      // variables) where the actuall clique value is obtained. This
      // points to a data structure maintained by origin.
      unsigned *ptr;
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

  // the collection of clique values for this clique.
  sArray< CliqueValue > cliqueValues;
  // Number of currently used clique values
  unsigned numCliqueValuesUsed;

  // Max collect-evidence probability for this clique. Used for beam
  // pruning.
  logpr maxCEValue;

public:

  // WARNING: constructor hack to create a very broken object with
  // non-functional reference objects (in order to create an array of
  // these objects and then initialize them later with appropriate
  // references). Do not use until after proper re-constructor.
  InferenceMaxClique() : origin(*((MaxClique*)NULL)) {}
  // normal constructor (i.e., re-constructor).
  InferenceMaxClique(MaxClique& _origin,
		     vector <RandomVariable*>& newRvs,
		     map < RVInfo::rvParent, unsigned >& ppf,
		     const unsigned int frameDelta);
  // destructor
  ~InferenceMaxClique() {}

  // collect evidence functions.
  void ceGatherFromIncommingSeparators(JT_InferencePartition& part);
  void ceIterateSeparators(JT_InferencePartition& part,
			   const unsigned sepNumber,
			   const logpr p);

  void ceIterateAssignedNodesRecurse(JT_InferencePartition& part,
				     const unsigned nodeNumber,
				     const logpr p);
  void ceIterateAssignedNodesNoRecurse(JT_InferencePartition& part,
				       const logpr p);

  inline void ceIterateAssignedNodes(JT_InferencePartition& part,
				     const unsigned nodeNumber,
				     const logpr p)
  {
    // ceIterateAssignedNodesRecurse(part,nodeNumber,p);
    ceIterateAssignedNodesNoRecurse(part,p);
  }


  void ceIterateUnassignedIteratedNodes(JT_InferencePartition& part,
					const unsigned nodeNumber,
					const logpr p);
  void ceSendToOutgoingSeparator(JT_InferencePartition& part,
			    InferenceSeparatorClique& sep); 
  void ceSendToOutgoingSeparator(JT_InferencePartition& part);
  void ceCliquePrune();

  // support for collect evidence clique driven operations.
  void ceGatherFromIncommingSeparatorsCliqueDriven(JT_InferencePartition& part);
  void ceIterateAssignedNodesCliqueDriven(JT_InferencePartition& part,
					  const unsigned nodeNumber,
					  logpr p);
  void ceIterateUnassignedNodesCliqueDriven(JT_InferencePartition& part,
					    const unsigned nodeNumber,
					    const logpr p);
  void ceGatherFromIncommingSeparatorsCliqueObserved(JT_InferencePartition& part);


  // distribute evidence functions.
  void deScatterToOutgoingSeparators(JT_InferencePartition& part);
  void deReceiveFromIncommingSeparator(JT_InferencePartition& part,
				       InferenceSeparatorClique& sep);
  void deReceiveFromIncommingSeparator(JT_InferencePartition& part);


  // sum up the probabilities in the current clique and return their value.
  logpr sumProbabilities();

  // EM accumulation support.
  void emIncrement(const logpr probE, 
		   const bool localCliqueNormalization = false);

};



class SeparatorClique : public IM
{
  friend class FileParser;
  friend class GraphicalModel;
  friend class GMTemplate;

public:

  // beam width for separator-based beam pruning.
  static double separatorBeam;

  // the set of nodes which forms the separator
  set<RandomVariable*> nodes;

  //////////////////////////////////////////////////////////////////////
  // For a given iteration order, the intersection of 'nodes' with all
  // nodes in previous seperators according to the iteration.
  set<RandomVariable*> accumulatedIntersection;
  // hidden version of above
  vector<RandomVariable*> hAccumulatedIntersection;
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
  set<RandomVariable*> remainder;
  // hidden version of above
  vector<RandomVariable*> hRemainder;
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


  SeparatorClique(const SeparatorClique& sep)
  { nodes = sep.nodes; }

  // construct a separator between two cliques
  SeparatorClique(MaxClique& c1, MaxClique& c2);

  SeparatorClique(SeparatorClique& from_sep,
		  vector <RandomVariable*>& newRvs,
		  map < RVInfo::rvParent, unsigned >& ppf,
		  const unsigned int frameDelta = 0);
  
  // create an empty separator
  SeparatorClique() {}

  // compute the weight (log10 state space) of this separator clique.
  float weight(const bool useDeterminism = true) const { 
    return MaxClique::computeWeight(nodes,NULL,useDeterminism); 
  }

  // prepare the last set of data structures so that clones of this
  // can be unrolled and inference can occur.
  void prepareForUnrolling();

  // print out everything in this clique to a file.
  void printAllJTInfo(FILE* f);


  // used to clear out hash table memory between segments
  void clearDataMemory() {
    if (accPacker.packedLen() > ISC_NWWOH_AI) {
      accValueHolder.prepare();
      accSepValHashSet.clear();
    }
    if (remPacker.packedLen() > ISC_NWWOH_RM) { 
      remValueHolder.prepare();
      remSepValHashSet.clear();
    }
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
class InferenceSeparatorClique : public IM
{
  friend class InferenceMaxClique;

  // Non-STL "f=fast" versions of arrays that are instantiated
  // only in the final unrolled versions of the separator cliques.
  sArray< RandomVariable*> fNodes;
  sArray< RandomVariable*> fAccumulatedIntersection;
  sArray< RandomVariable*> fRemainder;

  // Direct pointers to the values within the discrete *HIDDEN* RVs
  // within this sep clique.  Note, the observed discrete and
  // continuous variables are not contained here.
  // 1) one for the accumulated intersection
  sArray < RandomVariable::DiscreteVariableType*> accDiscreteValuePtrs;
  // 2) and one for the remainder variables.
  sArray < RandomVariable::DiscreteVariableType*> remDiscreteValuePtrs;

  // the original separator clique from which this object has been cloned.
  SeparatorClique& origin;

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
    // probability
    logpr p;
    // 
    // union upi {
    //   char bp[sizeof(logpr)];
    //   unsigned backPointer;
    // };
    // 
    // probability for distribute evidence pass
    logpr bp;
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
    // Note, by default, these hash tables are constructed empty and
    // invalid, and they need to be explicitly re-constructed.
    VHashMapUnsignedUnsignedKeyUpdatable iRemHashMap;

    // ensure that we start with nothing inserted.
    AISeparatorValue() { 
      numRemValuesUsed = 0;
      // fprintf(stderr,"in aisv init\n");
    }
  };

  // The collection of AI separator values.
  cArray< AISeparatorValue > separatorValues;
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
  VHashMapUnsignedUnsignedKeyUpdatable iAccHashMap;

public:


  // WARNING: constructor hack to create a very broken object with
  // non-functional reference objects (in order to create an array of
  // these objects and then initialize them later with appropriate
  // references). Do not use until after proper re-constructor.
  InferenceSeparatorClique() : origin(*((SeparatorClique*)NULL)) {}
  // normal (or re-)constructor
  InferenceSeparatorClique(SeparatorClique& _origin,
			   vector <RandomVariable*>& newRvs,
			   map < RVInfo::rvParent, unsigned >& ppf,
			   const unsigned int frameDelta);
  // destructor
  ~InferenceSeparatorClique() {}

  // separator based pruning
  void ceSeparatorPrune();

};



#endif
