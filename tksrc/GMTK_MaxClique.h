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
  const unsigned cliqueValueSize;

  // The amount that successive allocation chunks grow in size.  Must
  // be >= 1.0.
  const float growthFactor;

  // The size of the initial allocation unit. A chunk is an
  // array of packed clique values. When we allocate the n'th
  // chunk, we allocate k(n) new packed clique values, 
  // where k(n) = allocationUnitChunkSize*growthFactor^(n-1)
  const unsigned allocationUnitChunkSize;

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
};

class MaxClique : public IM {

  friend class FileParser;
  friend class GraphicalModel;
  friend class GMTemplate;
  friend class SeperatorClique;

public:


  // @@@ need to take out, here for now to satisify STL call of vector.clear().
  MaxClique& operator=(const MaxClique& f) {
    return *this;
  }

  // the set of nodes which form a max clique, in arbitrary order.
  set<RandomVariable*> nodes;


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


  // compute the weight (log10 state space) of this clique.
  float weight(const bool useDeterminism = true) const { 
    return computeWeight(nodes,NULL,useDeterminism); 
  }
  // Static version of routine, useful in certain places outside this
  // class where we just have a collection of nodes to compute the
  // weight of.
  static float computeWeight(const set<RandomVariable*>& nodes,
			     const RandomVariable* node = NULL,
			     const bool useDeterminism = true);

  // print just the clique nodes
  void printCliqueNodes(FILE* f);


  
  //////////////////////////////////////////////
  // TODO: figure out a way so that the member variables below exist only in
  // a subclass since all of the below is not needed for basic
  // triangulation.

  // USED ONLY IN JUNCTION TREE INFERENCE
  // The set of nodes that are assigned to this
  // maxClique from which probabilities are extracted.
  // Necessarily it is the case
  // that assignedNodes <= nodes (i.e., set containment)
  // TODO: needs to be an array in topological partial order.
  set<RandomVariable*> assignedNodes;

  // USED ONLY IN JUNCTION TREE INFERENCE
  // a topologically sorted vector version of assignedNodes
  // The topological sort is relative only to the variables
  // in this clique.
  vector<RandomVariable*> sortedAssignedNodes;

  // USED ONLY IN JUNCTION TREE INFERENCE 
  // predicate for each each sorted assigned node, used to say if that
  // node should be iterated, or if its value has already been set by
  // a separator iteration. If this has zero length, then we iterate
  // everything. If it does not have zero length, it has same length
  // as sortedAssignedNodes array.
  sArray < bool > iterateSortedAssignedNodesP;


  // USED ONLY IN JUNCTION TREE INFERENCE
  // The set of nodes that are not
  // assigned to this clique but that also are not iterated over by
  // the incomming separators. Therefore, these nodes in the clique
  // must be iterated over from scratch.  The hope is that assignment
  // and setup can be done so there are very few or zero of such
  // nodes. Note however that once a clique iterates over such a node,
  // any later cliques (closer to the JT root) will not have to re-do
  // this as these assignments here will be represented in this
  // cliques out-going sepset.
  set<RandomVariable*> unassignedIteratedNodes;

  // USED ONLY IN JUNCTION TREE INFERENCE
  // The set of assigned nodes cummulative in the JT relative to root
  // in the JT for the current partition.  This variable is used only
  // by the code that assigns CPTs to cliques. Note that
  // cumulativeAssignedNodes includes the assignedNodes
  // in the current clique.
  set<RandomVariable*> cumulativeAssignedNodes;

  // USED ONLY IN JUNCTION TREE INFERENCE
  // the preceding cumulative set of unassigned and iterated nodes
  // relative to the root in the JT for the current partition adn
  // clique. This is used to determine which of the assigned nodes
  // need actually be iterated within a clique, and which are already
  // set by the separator iterations. Note that
  // cumulativeUnassignedIteratedNodes does *NOT* include the
  // assignedNodes in the current clique (thus the name preceding).
  set<RandomVariable*> precedingUnassignedIteratedNodes;


  // USED ONLY IN JUNCTION TREE INFERENCE
  // These are the nodes that are hidden (i.e., they are
  // the non-continous hidden variables). These are
  // the nodes whose values will be hashed.
  vector<RandomVariable*> hiddenNodes;

  // USED ONLY IN JUNCTION TREE INFERENCE
  // Clique neighbors in a junction tree, integer indexes into a table
  // of cliques contained in the Partition class containing this
  // maxclique. By storing ints rather than actual pointers to
  // maxcliques, cloned partitions (with same corresponding clique
  // array ordering) can use the same structure (ints offsets 
  // in an array) on different cliques in unrolled partitions.
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
  vector<unsigned> children;

  // USED ONLY IN JUNCTION TREE INFERENCE
  // set of separators that we receive from in the collect evidence stage.
  // Again, ints indexing into parent partition.
  vector<unsigned> ceReceiveSeparators;

  // USED ONLY IN JUNCTION TREE INFERENCE
  // separator that we send to in the collect collect evidence stage
  // Again, ints indexing into parent partition.
  unsigned ceSendSeparator;

  // USED ONLY IN JUNCTION TREE INFERENCE
  // structure used to pack and unpack clique values
  PackCliqueValue packer;

  // USED ONLY IN JUNCTION TREE INFERENCE 
  // structure to allocate clique values from.  These things are
  // allocated and deleted in bulk (rather than individually) and that
  // is handled by this object, thereby acting effectively as a
  // customized memory management unit.
  // 
  // Note that if a packed clique value can be held in <
  // sizeof(unsigned)*8, then this isn't used since the pointer
  // storage can be used instead.
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
  // sizeof(unsigned) bytes, then this isn't used since the pointer
  // storage can be used instead.
  vhash_set< unsigned > cliqueValueHashSet;

  // USED ONLY IN JUNCTION TREE INFERENCE
  // a value that is roughly equal to the order of the state space of
  // the clique. It should be smaller than the clique state space,
  // however, since pruning will make the state space smaller.
  // This will be used for initial memory allocation units, etc.
  // by this clique object and clones of this clique object.
  // It is only usable after we have prepared for unrolling.
  unsigned allocationUnitChunkSize;

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
  void computeAssignedNodesToIterate();

  // USED ONLY IN JUNCTION TREE INFERENCE
  // print out everything in this junction tree clique to a file.
  void printAllJTInfo(FILE* f,const unsigned indent);

};



//////////////////////////////////////////////////////////////////////////////
// Number of words to exist in a pack clique without a hash. These may
// be set to zero to turn on hashing even when packed clique values
// fall into less than 1 machine word.
// --
// InferenceMaxClique Number Words WithOut a Hash: Namely,
// the number of words that can be stored directly as
// a packed clique value before we resort to using
// a shared hash table for all instances of this origin clique.
#define IMC_NWWOH (0)
// InferenceSeparatorClique Number Words WithOut a Hash: Namely,
// the number of words that can be stored directly as
// a packed clique value before we resort to using
// a shared hash table for all instances of this origin clique.
// One for the accumulated Intersection packed values
#define ISC_NWWOH_AI (0)
// And for the remainder
#define ISC_NWWOH_RM (0)
// -- 
//////////////////////////////////////////////////////////////////////////////


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

    // probability
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
  void collectEvidenceFromSeparators(JT_InferencePartition& part);
  void ceIterateSeparators(JT_InferencePartition& part,
			   const unsigned sepNumber,
			   const logpr p);
  void ceIterateAssignedNodes(JT_InferencePartition& part,
			      const unsigned nodeNumber,
			      const logpr p);
  void ceIterateUnassignedIteratedNodes(JT_InferencePartition& part,
					const unsigned nodeNumber,
					const logpr p);
  void ceCollectToSeparator(JT_InferencePartition& part,
			    InferenceSeparatorClique& sep); 
  void ceCollectToSeparator(JT_InferencePartition& part);

  // sum up the probabilities in the current clique and return their value.
  logpr sumProbabilities();

};




class SeparatorClique : public IM
{
  friend class FileParser;
  friend class GraphicalModel;
  friend class GMTemplate;

public:

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

  // prepare the last set of data structures so that clones of this
  // can be unrolled and inference can occur.
  void prepareForUnrolling();

  // print out everything in this clique to a file.
  void printAllJTInfo(FILE* f);

  // @@@ need to take out, here for now to satisify STL call of vector.clear().
  SeparatorClique& operator=(const SeparatorClique& f) {
    return *this;
  }


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


};



#endif
