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
 * Copyright (C) 2015 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 *
 */

// TODO: perhaps create a subclass or member of maxClique at some point, rather than
// adding everything for exact inference to the base class.


#ifndef GMTK_PEDAGOGICALCLIQUETABLE_H
#define GMTK_PEDAGOGICALCLIQUETABLE_H

#if HAVE_CONFIG_H
#  include <config.h>
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
#include "GMTK_ConditionalSeparatorTable.h"

#include <stdio.h>
#include <stdlib.h>

#include <vector>
#include <string>
#include <map>

class PartitionStructures;
class SeparatorClique;
class ConditionalSeparatorTable;
//class ConditionalSeparatorTable::SharedLocalStructure;
class MaxClique;


// The table entries for a maxclique that:
//   1) has no STL and uses only fast datastructures with custom memory managment
//   2) keeps a pointer back to it's original clique for hash tables, etc.
// 
// Note that this is not a subclass of MaxClique since we want this
// object to be as absolutely small as possible to save memory on very
// long segments (e.g., genomes), and we definitely do not want these
// objects to have to have all MaxClique's member variables which
// would be highly redundant.
class PedagogicalCliqueTable  : public IM
{
  friend struct CliqueValueDescendingProbCompare;
  friend class PartitionStructures;
  friend class SectionTablesBase;
  friend class SparseJoinSectionTables;

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
    // to all PedagogicalCliqueTable's of this particular base clique.
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
    sArray< RV*> fUnassignedNodes;
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
    vector<RV*> returnRVsAsVector();

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

  // Memory management options set by -memoryGrowth
  static float valuePoolGrowthRate;

  // WARNING: constructor hack to create a very broken object with
  // non-functional reference objects (in order to create an array of
  // these objects and then initialize them later with appropriate
  // references). Do not use until after proper re-constructor.
  PedagogicalCliqueTable() {}
  // normal constructor (i.e., re-constructor).
  PedagogicalCliqueTable(MaxClique& origin);
  void init(MaxClique& origin);

  // destructor
  ~PedagogicalCliqueTable() {}

  ///////////////////////////////////////////////////////////////////////////////
  // collect evidence functions.

  void ceGatherFromIncommingSeparators(PedagogicalCliqueTable::SharedLocalStructure& sharedStructure,
				       ConditionalSeparatorTable* separatorTableArray,
				       ConditionalSeparatorTable::SharedLocalStructure* sepSharedStructureArray);
				       
  // clique driven pedagogical inference
  void ceGatherFromIncommingSeparatorsCliqueDriven(PedagogicalCliqueTable::SharedLocalStructure& sharedStructure,
						   ConditionalSeparatorTable* separatorCliqueArray,
						   ConditionalSeparatorTable::SharedLocalStructure* sepSharedStructureArray,
						   logpr cliqueBeamThresholdEstimate,
						   logpr& maxCEValue);

  // special case when the clique is fully observed.
  void ceGatherFromIncommingSeparatorsCliqueObserved(PedagogicalCliqueTable::SharedLocalStructure& sharedStructure,
						     ConditionalSeparatorTable* separatorCliqueArray,
						     ConditionalSeparatorTable::SharedLocalStructure* sepSharedStructureArray,
						     logpr& maxCEValue);


  ////////////////////////////////////////////////////////////////////////////
  // iterate unassigned nodes

  void ceIterateUnassignedNodesCliqueDriven(PedagogicalCliqueTable::SharedLocalStructure& sharedStructure,
					    ConditionalSeparatorTable *separatorTableArray,
					    ConditionalSeparatorTable::SharedLocalStructure* sepSharedStructureArray,
					    logpr cliqueBeamThresholdEstimate,
					    logpr& maxCEValue,
					    const unsigned nodeNumber,
					    const logpr p);

  ///////////////////////////////////////
  // iterate assigned nodes.

  void ceIterateAssignedNodesCliqueDriven(PedagogicalCliqueTable::SharedLocalStructure& sharedStructure,
					  ConditionalSeparatorTable *separatorTableArray,
					  ConditionalSeparatorTable::SharedLocalStructure* sepSharedStructureArray,
					  logpr cliqueBeamThresholdEstimate,
					  logpr& maxCEValue,
					  const unsigned nodeNumber,
					  logpr p);


  /////////////////////////////////////////
  // a version that automatically selects which separator to use within partition from the clique.
  void ceSendToOutgoingSeparators(PedagogicalCliqueTable::SharedLocalStructure& sharedStructure,
				  ConditionalSeparatorTable* separatorTableArray,
				  ConditionalSeparatorTable::SharedLocalStructure* sepSharedStructureArray);
  // a version with an explicit outgoing separator.
  void ceSendToOutgoingSeparator(PedagogicalCliqueTable::SharedLocalStructure& sharedStructure,
				 ConditionalSeparatorTable& sep,
				 ConditionalSeparatorTable::SharedLocalStructure&);




  /////////////////////////////////////////
  // memory clearing.
  void clearInferenceMemory() {
    // clear out all memory used by this inference clique.
    cliqueValues.clear();
  }
  void clearCliqueAndIncommingSeparatorMemory(PedagogicalCliqueTable::SharedLocalStructure& sharedStructure,
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

  void ceDoCliqueScoreNormalization(PedagogicalCliqueTable::SharedLocalStructure& sharedStructure);


#ifdef USE_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL
  void insertLocalCliqueValuesIntoSharedPool(MaxClique& origin);
#endif


  /////////////////////////////////////////
  // distribute evidence functions.
  /////////////////////////////////////////

  void deScatterToOutgoingSeparators(PedagogicalCliqueTable::SharedLocalStructure& sharedStructure,
				     ConditionalSeparatorTable* separatorTableArray,
				     ConditionalSeparatorTable::SharedLocalStructure* sepSharedStructureArray);


  void deScatterToOutgoingSeparatorsViterbi(PedagogicalCliqueTable::SharedLocalStructure& sharedStructure,
					    ConditionalSeparatorTable* separatorTableArray,
					    ConditionalSeparatorTable::SharedLocalStructure* sepSharedStructureArray);

  // a version that automatically selects which separator to use within partition from the clique.
  void deReceiveFromIncommingSeparator(PedagogicalCliqueTable::SharedLocalStructure& sharedStructure,
				       ConditionalSeparatorTable* separatorTableArray,
				       ConditionalSeparatorTable::SharedLocalStructure* sepSharedStructureArray);
  // a version with an explicit incomming separator
  void deReceiveFromIncommingSeparator(PedagogicalCliqueTable::SharedLocalStructure& sharedStructure,
				       ConditionalSeparatorTable&,
				       ConditionalSeparatorTable::SharedLocalStructure&);
  // a version specific to viterbi decoding.
  void deReceiveFromIncommingSeparatorViterbi(PedagogicalCliqueTable::SharedLocalStructure& sharedStructure,
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
  logpr maxProbability(PedagogicalCliqueTable::SharedLocalStructure&,
		       bool setCliqueToMaxValue = true);
  // a faster version that just operates on the table.
  logpr maxProb();

  // print all clique values and prob to given file.
  void printCliqueEntries(PedagogicalCliqueTable::SharedLocalStructure&,
			  FILE*f,const char*str=NULL,
			  const bool normalize = true, const bool unlog = true,
			  const bool justPrintEntropy = false);
  
  void printCliqueEntries(PedagogicalCliqueTable::SharedLocalStructure&,
			  ObservationFile *f, const bool normalize = true, 
			  const bool unlog = true);
  
  int cliqueValueDistance(SharedLocalStructure& sharedStructure, 
			  unsigned a, unsigned b);

  static unsigned cliqueDomainSize(SharedLocalStructure& sharedStructure);

  static void printCliqueOrder(FILE *f, SharedLocalStructure& sharedStructure, int frameDelta=0);

  unsigned cliqueValueMagnitude(SharedLocalStructure& sharedStructure, unsigned cliqueIndex);


  class CliqueValueIndex {

    SharedLocalStructure *sharedStructure;
    PedagogicalCliqueTable       *table;

  public:
    unsigned index;

    CliqueValueIndex(SharedLocalStructure *sharedStructure, 
		     PedagogicalCliqueTable       *table,
		     unsigned index)
      : sharedStructure(sharedStructure), table(table), index(index)
    {}

    CliqueValueIndex() 
      : sharedStructure(NULL), table(NULL), index(0)
    {}
 
    bool operator<(const CliqueValueIndex& rhs) const {
      return table->cliqueValueDistance(*sharedStructure, index, rhs.index) < 0;
    }

    bool operator>(const CliqueValueIndex& rhs) const {
      return table->cliqueValueDistance(*sharedStructure, index, rhs.index) > 0;
    }

    bool operator==(const CliqueValueIndex& rhs) const {
      return table->cliqueValueDistance(*sharedStructure, index, rhs.index) == 0;
    }

    CliqueValueIndex& operator=(CliqueValueIndex rhs) {
      this->sharedStructure = rhs.sharedStructure;
      this->table = rhs.table;
      this->index = rhs.index;
      return *this;
    }
    
  };


  // EM accumulation support.
  void emIncrement(PedagogicalCliqueTable::SharedLocalStructure&,
		   const logpr probE, 
		   const bool localCliqueNormalization = false,
		   const double emTrainingBeam = -LZERO);


};



#endif
