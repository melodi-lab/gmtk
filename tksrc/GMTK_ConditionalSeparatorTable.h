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


#ifndef GMTK_CONDITIONALSEPARATORCLIQUE_H
#define GMTK_CONDITIONALSEPARATORCLIQUE_H

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
//class ConditionalSeparatorTable::SharedLocalStructure;
class MaxClique;
class MaxCliqueTable;

/////////////////////////////
// TODO: think of a way of subclassing or something so that members
// are only in the objects when they are used. (perhaps one base class
// with two children will work).
/////////////////////////////


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

    // note that this casts through void * to avoid a type
    // punning warning. this is potentially unsafe, but we
    // believe it works OK on our target platforms
    inline logpr& bp() { 
	void *t = (void*) (&_bpo[0]);
	logpr *p = reinterpret_cast<logpr*>(t);
	return *p ;
    }

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

  // Preserve separator for later use in linear space inference
  bool preserve;

  // Memory management options set by -memoryGrowth
  static unsigned remHashMapStartingSize;

  // WARNING: constructor hack to create a very broken object with
  // non-functional reference objects (in order to create an array of
  // these objects and then initialize them later with appropriate
  // references). Do not use until after proper re-constructor.
  ConditionalSeparatorTable() : preserve(false)
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
    if (!preserve) {
      delete iAccHashMap;
      iAccHashMap = NULL;
      delete separatorValues;
      separatorValues = NULL;
    }
  }


  // insert current value of RVs into separator
  void insert(SeparatorClique& origin,
	      ConditionalSeparatorTable::SharedLocalStructure& sepSharedStructure);

  // separator based pruning
  void ceSeparatorPrune(SeparatorClique& origin);

  // memory reporting
  void reportMemoryUsageTo(SeparatorClique&,FILE *f);

};



#endif
