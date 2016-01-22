/*-
 * GMTK_MaxClique.cc
 *
 *     maxClique support and JT probabilistic inference. Includes the
 *     implementation of a message in the message passing algorithm.
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (C) 2003 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 *
 */




/*
 * TODO: turn this into multiple files
 *   mc, mctable CE, mctable DE, mctable prune, csctable
 *
 */



/*
 * TODO:
 *   1) Maybe: In some of the structures the choice of when
 *   to use a hash table vs. an inlined value is done
 *   by either using an unsigned* or a unsigned. Create
 *   options where we do 'unsigned val[LEN]' and where
 *   only if packed value is > LEN words do we resort to
 *   the hash table.
 *
 *   2) Figure out way to remove the (unsigned**) type casts
 *   in front of call to packer object (see below).   
 *
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif


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
#include "psp.h"

#include "GMTK_FileParser.h"
#include "GMTK_RV.h"
#include "GMTK_DiscRV.h"
#include "GMTK_HidDiscRV.h"
#include "GMTK_ObsDiscRV.h"
#include "GMTK_ContRV.h"
#include "GMTK_GMTemplate.h"
#include "GMTK_GMParms.h"
#include "GMTK_MaxClique.h"
#include "GMTK_MDCPT.h"
#include "GMTK_MSCPT.h"
#include "GMTK_MTCPT.h"
#include "GMTK_Mixture.h"
#include "GMTK_ObservationSource.h"
#include "GMTK_SectionScheduler.h"

#if HAVE_CONFIG_H
#include <config.h>
#endif
#if HAVE_HG_H
#include "hgstamp.h"
#endif
VCID(HGID)

#if 0


#if 0
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Constants
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////
// Comment/Uncomment to optimize for speed/reducing memory usage.
#ifndef OPTIMIZE_FOR_MEMORY_USAGE
#define OPTIMIZE_FOR_MEMORY_USAGE
#endif
////////////////////////////////////////////////////////////////////////

#ifdef OPTIMIZE_FOR_MEMORY_USAGE

#ifdef USE_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL
#define TEMPORARY_LOCAL_CLIQUE_VALUE_POOL_GROWTH_RATE 1.15
#endif

// reasonable value is 1.25, but could go as low as 1.05 or so.
#define MEM_OPT_GROWTH_RATE 1.25

#define CLIQUE_VALUE_HOLDER_STARTING_SIZE 23
#define CLIQUE_VALUE_HOLDER_GROWTH_RATE   MEM_OPT_GROWTH_RATE

#define AI_SEP_VALUE_HOLDER_STARTING_SIZE 1
#define AI_SEP_VALUE_HOLDER_GROWTH_RATE   MEM_OPT_GROWTH_RATE

#define REM_SEP_VALUE_HOLDER_STARTING_SIZE 1
#define REM_SEP_VALUE_HOLDER_GROWTH_RATE   MEM_OPT_GROWTH_RATE

#define REM_HASH_MAP_STARTING_SIZE 1

#define CLIQUE_VALUE_SPACE_MANAGER_GROWTH_RATE   MEM_OPT_GROWTH_RATE
#define CLIQUE_VALUE_SPACE_MANAGER_DECAY_RATE    0.0

#define SEPARATOR_VALUE_SPACE_MANAGER_GROWTH_RATE  MEM_OPT_GROWTH_RATE
#define SEPARATOR_VALUE_SPACE_MANAGER_DECAY_RATE   0.0

#define REMAINDER_VALUE_SPACE_MANAGER_GROWTH_RATE  MEM_OPT_GROWTH_RATE
#define REMAINDER_VALUE_SPACE_MANAGER_DECAY_RATE   0.0

#else
// Then we optimize more for speed at the expense of memory usage.
// Note, with these settings, information about memory usage will be
// used from one segment to pre-allocate structures for the next
// segment processed (which means that if segment i uses less memory
// than segment i+1, we might use more memory than segment i+1
// needs).


#ifdef USE_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL
#define TEMPORARY_LOCAL_CLIQUE_VALUE_POOL_GROWTH_RATE 2.0
#endif


#define CLIQUE_VALUE_HOLDER_STARTING_SIZE 23
#define CLIQUE_VALUE_HOLDER_GROWTH_RATE   2.0

#define AI_SEP_VALUE_HOLDER_STARTING_SIZE 23
#define AI_SEP_VALUE_HOLDER_GROWTH_RATE   2.0

#define REM_SEP_VALUE_HOLDER_STARTING_SIZE 23
#define REM_SEP_VALUE_HOLDER_GROWTH_RATE   2.0

#define REM_HASH_MAP_STARTING_SIZE 2

#define CLIQUE_VALUE_SPACE_MANAGER_GROWTH_RATE   2.0
#define CLIQUE_VALUE_SPACE_MANAGER_DECAY_RATE    0.9

#define SEPARATOR_VALUE_SPACE_MANAGER_GROWTH_RATE  2.0
#define SEPARATOR_VALUE_SPACE_MANAGER_DECAY_RATE   0.9

#define REMAINDER_VALUE_SPACE_MANAGER_GROWTH_RATE  2.0
#define REMAINDER_VALUE_SPACE_MANAGER_DECAY_RATE   0.9

#endif
#endif




#endif


// for sorting an array of RVs ascending based on increasing cardinality
struct ParentCardinalityCompare 
{  
  bool operator() (const RV* rv1,
		   const RV* rv2) 
  {
    // place observed ones first.
    if (RV2DRV(rv1)->discreteObservedImmediate())
      return false;
    else if (RV2DRV(rv2)->discreteObservedImmediate())
      return true;
    else return (RV2DRV(rv1)->cardinality < RV2DRV(rv2)->cardinality);
  }
};







////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        MaxClique support
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * MaxClique::MaxClique()
 *    Clone constructor with frame delta to create a clone but under an unrolling.
 *    I.e., this isn't really a normal constructor, this is a contructor that
 *    sets up a clone of the MaxClique given by the argument from_clique. The
 *    clone is adjusted so that:
 *        1) it is shifted in time by frameDelta
 *        2) not all data structures are retained, only the ones 
 *           necessary to do exact inference.
 *
 * Preconditions:
 *
 * Postconditions:
 *
 * Side Effects:
 *
 * Results:
 *     nothing
 *
 *-----------------------------------------------------------------------
 */
MaxClique::MaxClique(MaxClique& from_clique,
		     vector <RV*>& newRvs,
		     map < RVInfo::rvParent, unsigned >& ppf,
		     const unsigned int frameDelta)

  : MaxCliqueBase(), // MaxCliqueBase(from_clique.nodes), 
     cliqueValueSpaceManager(1,     // starting size
			     spaceMgrGrowthRate,   // growth rate
			     1,     // growth addition
			     spaceMgrDecayRate)    // decay rate 
{
  set<RV*>::iterator it;


  // clone over nodes RVs.
  for (it = from_clique.nodes.begin();
       it != from_clique.nodes.end();
       it++) {
    RV* rv = (*it);
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    

    if ( ppf.find(rvp) == ppf.end() ) {
      coredump("ERROR: can't find rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	       rv->name().c_str(),rv->frame(),frameDelta,
	       rvp.first.c_str(),rvp.second);
    }

    RV* nrv = newRvs[ppf[rvp]];
    nodes.insert(nrv);
  }

  // and clone over assigned nodes and sorted assigned nodes
  sortedAssignedNodes.reserve(from_clique.sortedAssignedNodes.size());
  for (unsigned i=0;i<from_clique.sortedAssignedNodes.size();i++) {
    RV* rv = from_clique.sortedAssignedNodes[i];
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    

    if ( ppf.find(rvp) == ppf.end() ) {
      coredump("ERROR: can't find assigned rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	       rv->name().c_str(),rv->frame(),frameDelta,
	       rvp.first.c_str(),rvp.second);
    }

    RV* nrv = newRvs[ppf[rvp]];
    assignedNodes.insert(nrv);
    sortedAssignedNodes.push_back(nrv);
  }

  // do unassignedIteratedNodes
  for (it = from_clique.unassignedIteratedNodes.begin();
       it != from_clique.unassignedIteratedNodes.end();
       it++) {
    RV* rv = (*it);
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    

    if ( ppf.find(rvp) == ppf.end() ) {
      coredump("ERROR: can't find assigned rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	       rv->name().c_str(),rv->frame(),frameDelta,
	       rvp.first.c_str(),rvp.second);
    }
    RV* nrv = newRvs[ppf[rvp]];
    unassignedIteratedNodes.insert(nrv);
  }

  // do cumulativeAssignedNodes
  for (it = from_clique.cumulativeAssignedNodes.begin();
       it != from_clique.cumulativeAssignedNodes.end();
       it++) {
    RV* rv = (*it);
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    

    if ( ppf.find(rvp) == ppf.end() ) {
      coredump("ERROR: can't find assigned rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	       rv->name().c_str(),rv->frame(),frameDelta,
	       rvp.first.c_str(),rvp.second);
    }
    RV* nrv = newRvs[ppf[rvp]];
    cumulativeAssignedNodes.insert(nrv);
  }

  // these are just integer indices, so we can copy them.
  neighbors = from_clique.neighbors;
  children = from_clique.children;
  ceReceiveSeparators = from_clique.ceReceiveSeparators;
  ceSendSeparator = from_clique.ceSendSeparator;

  if (cliqueBeamBuildBeam != (-LZERO)) {
    // then we're doing clique build pruning.
    if (cliqueBeamBuildFilter == NULL || strlen(cliqueBeamBuildFilter) == 0 || strncmp(cliqueBeamBuildFilter,"fixed",5) == 0 ) {
      maxCEValuePredictor = 
	counted_ptr<AdaptiveFilter>(new FixedFilter());
    } else if ( strncmp(cliqueBeamBuildFilter,"lms",3) == 0 ) {
      // TODO: make a routine to parse arguments like this.

      // expecting syntax of the form "lms,3,0.9"
      unsigned order = 3;
      double lr = 0.9;
      unsigned slen = strlen(cliqueBeamBuildFilter);
      if (slen > 4) {
	char *startp = &cliqueBeamBuildFilter[4];
	char *endp;
	unsigned tmp = strtol(startp, &endp, 10);
	if (endp != startp) {
	  order = tmp;
	  if (*endp == ',') {
	    startp = endp+1;
	    double d = strtod(startp, &endp);
	    if (endp != startp) {
	      lr = d;
	    }
	  }
	}
      }
      // fprintf(stderr,"lms order = %d, lr = %f\n",order,lr);
      maxCEValuePredictor = 
	counted_ptr<AdaptiveFilter>(new LMSFilter(order,lr));
    } else if ( strncmp(cliqueBeamBuildFilter,"rls",3) == 0 ) {
      // expecting syntax of the form "rls,3,1.0"
      unsigned order = 3;
      double fc = 1.0;
      unsigned slen = strlen(cliqueBeamBuildFilter);
      if (slen > 4) {
	char *startp = &cliqueBeamBuildFilter[4];
	char *endp;
	unsigned tmp = strtol(startp, &endp, 10);
	if (endp != startp) {
	  order = tmp;
	  if (*endp == ',') {
	    startp = endp+1;
	    double d = strtod(startp, &endp);
	    if (endp != startp) {
	      fc = d;
	    }
	  }
	}
      }
      // fprintf(stderr,"rls order = %d, fc = %f\n",order,fc);
      maxCEValuePredictor =
	counted_ptr<AdaptiveFilter>(new RLSFilter(order,fc));
    } else {
      error("Error: unknown clique build beam pruning filter type '%s'\n",
	    cliqueBeamBuildFilter);
    }
    prevMaxCEValPrediction = 0.0;
  } else {
    // nothing to do.
    // maxCEValuePredictor = NULL;
  }

}


/*-
 *-----------------------------------------------------------------------
 * MaxClique::clearCliqueValueCache()
 *   
 *   clear out hash and clique value memory for use between segments
 *   (but not yet working between partitions/frames/slices/chunks).
 *
 * Preconditions:
 *   Data structures must be set up.
 *
 * Postconditions:
 *   all hash and clique value memory is released.
 *
 * Side Effects:
 *     changes internal data structures.
 *
 * Results:
 *     none
 *
 *
 *-----------------------------------------------------------------------
 */
void
MaxClique::clearCliqueValueCache(bool force)
{
  if (force && packer.packedLen() > IMC_NWWOH) {
    valueHolder.prepare();
    cliqueValueHashSet.clear(spaceMgrStartingSize);
#ifdef USE_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL
    temporaryCliqueValuePool.clear();
#endif
  }
  // shrink space asked for by clique values. 
  cliqueValueSpaceManager.decay();
}




/*-
 *-----------------------------------------------------------------------
 * MaxClique::reportMemoryUsageTo()
 *
 *    Report current memory usage of this inference clique *and* the origin clique to the file in
 *    units of MBs.
 *
 * Preconditions:
 *      Clique data structures must be created.
 *
 * Postconditions:
 *      results printed.
 *
 * Side Effects:
 *      None
 *
 * Results:
 *     none
 *
 *-----------------------------------------------------------------------
 */
void
MaxClique::
reportMemoryUsageTo(FILE *f)
{

  if (packer.packedLen() > IMC_NWWOH) {
#ifdef USE_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL
    fprintf(f,"OC(VH=%luMB,HS=%luMB,TP=%luMB)",
	    1+valueHolder.bytesRequested()/(1024*1024),
	    1+cliqueValueHashSet.bytesRequested()/(1024*1024),
	    1+(sizeof(unsigned)*((unsigned long)temporaryCliqueValuePool.size()))/(1024*1024));
#else
    fprintf(f,"OC(VH=%luMB,HS=%luMB)",
	    1+valueHolder.bytesRequested()/(1024*1024),
	    1+cliqueValueHashSet.bytesRequested()/(1024*1024));
#endif
  } else {
    // in this case, we report zero usage since the clique values are
    // all being stored in the inference clique.
#ifdef USE_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL
    fprintf(f,"OC(VH=0MB,HS=0MB,TP=0MB)");
#else
    fprintf(f,"OC(VH=0MB,HS=0MB)");
#endif
  }
}



// TODO: this should probably be a pure virtual method in MaxCliqueBase
/*-
 *-----------------------------------------------------------------------
 * MaxClique::prepareForUnrolling()
 *   
 *   Sets a few last data structures that are needed for both unrolling
 *   and inference to work. 
 *
 * Preconditions:
 *   Set of nodes must be valid meaning that it has valid neighbors,
 *   parents, and children member variables.
 *
 * Postconditions:
 *   computed weight is provided.
 *
 * Side Effects:
 *     none
 *
 * Results:
 *     none
 *
 *
 *-----------------------------------------------------------------------
 */
void
MaxClique::prepareForUnrolling()
{

  set<RV*> hiddenNodesSet;

  // set up number of hidden ndoes
  set<RV*>::iterator it;
  for (it = nodes.begin();
       it != nodes.end();
       it++) {
    RV* rv = (*it);
    if (rv->hidden()) {
      hiddenNodes.push_back(rv);
      hiddenNodesSet.insert(rv);
    }
  }

  if (storeDeterministicChildrenInClique) {
    hashableNodes = hiddenNodes;
    determinableNodes.clear();
  } else {

    vector<RV*> hiddenNodesSorted;
    GraphicalModel::topologicalSort(hiddenNodesSet,hiddenNodesSet,hiddenNodesSorted);

    determinableNodes.clear();
    hashableNodes.clear();
    for (unsigned n=hiddenNodesSorted.size();n>0;n--) {
      RV* rv = hiddenNodesSorted[n-1];
      if (   rv->discrete() 
	     && rv->hidden()    // we know the first two are true but include them for clarity.
	     && RV2DRV(rv)->deterministic() ) {
	// printf("found deterministic node, in assignedNodes = %d :",(assignedNodes.find(rv) != assignedNodes.end())); 
	// rv->printNameFrame(stdout,true);

	// Now check to see if this node hidden parents live in earlier
	// nodes in clique. Because we are checking the nodes reverse
	// topologically, we need not worry about the nodes that have
	// already been removed. Also, we don't care about this nodes
	// observed parents, since any observed parents can live
	// anywhere.

	// create a set version of all the parents for this rv.
	set<RV*> parSet;
	for (unsigned j=0;j<rv->allParents.size();j++) {
	  parSet.insert(rv->allParents[j]);
	}
	set<RV*> parentsInClique;
	set_intersection(parSet.begin(),parSet.end(),
			 nodes.begin(),nodes.end(),
			 inserter(parentsInClique,parentsInClique.end()));
	// printf("number of parents in clique = %d, num parents = %d\n",parentsInClique.size(),rv->allParents.size());
	if (parentsInClique.size() == rv->allParents.size()) {
	  // then rv is discrete, hidden, deterministic, and its parents live
	  // in the earlier nodes in the clique. 
	  determinableNodes.push_back(rv);
	  continue; // onto next loop
	} else {
	  // check if all the parents out of the clique are all observed, and
	  // if so, we're still fine (i.e., the node rv can be removed).
	  set<RV*> parentsOutOfClique;
	  set_difference(parSet.begin(),parSet.end(),
			 parentsInClique.begin(),parentsInClique.end(),
			 inserter(parentsOutOfClique,parentsOutOfClique.end()));
	  set<RV*>::iterator it;
	  bool parentsOutOfCliqueObserved = true;
	  // printf("num parents out of clique = %d: ",parentsOutOfClique.size());
	  for (it = parentsOutOfClique.begin(); it != parentsOutOfClique.end(); it++) {
	    RV* par = (*it);
	    if (par->hidden()) {
	      parentsOutOfCliqueObserved = false;
	      // par->printNameFrame(stdout,false);
	      // break;
	    }
	  }
	  // printf("\n");

	  if (parentsOutOfCliqueObserved) {
	    determinableNodes.push_back(rv);
	    continue; // onto next loop
	  }
	}
      }
      // we failed, so this must be a hashable node.
      hashableNodes.push_back(rv);
    }

    assert( hiddenNodes.size() == hashableNodes.size() + determinableNodes.size() );

    // Reverse nodes so that increasing order is the order that they
    // should be set (i.e., so that later ones might depend on earlier
    // ones, but not the reverse).
    reverse(determinableNodes.begin(),determinableNodes.end());

    // printf("Number of words required for original %d hidden variables = %d\n",
    // hiddenNodes.size(),PackCliqueValue::numWordsRequiredFor(hiddenNodes));
    // printf("Number of words required for %d hidden hashable variables w/o determinism = %d\n",
    // hashableNodes.size(),PackCliqueValue::numWordsRequiredFor(hashableNodes));

    const unsigned hash_pack_size = PackCliqueValue::numWordsRequiredFor(hashableNodes);
    const unsigned orig_pack_size = PackCliqueValue::numWordsRequiredFor(hiddenNodes);
    assert ( hash_pack_size <= orig_pack_size );

    infoMsg(High,"Packed clique size %d words with %d all hidden variables, and %d words without %d deterministic vars\n",
	    orig_pack_size,hiddenNodes.size(),hash_pack_size,determinableNodes.size());

    if (hash_pack_size == orig_pack_size) {
      // then don't bother with not storing the deterministic nodes
      // since we wouldn't save any memory, and it is faster to not have to 
      // recompute teh determinable nodes.
      hashableNodes = hiddenNodes;
      determinableNodes.clear();
    } else {
      // So hash_pack_size < orig_pack_size,
      // Last, we take some of the nodes in determinable nodes and add them
      // back to hashable nodes as long as the number of words does not
      // increase. Do this in increasing order of node cardinality.

      // TODO: option to do this in increasing order of complexity of evaluating
      // a deterministic node (i.e., store the nodes if they involve a complicated
      // deterministic evaluator).

      vector<RV*> determinableNodesSortedByCard = determinableNodes;
      sort(determinableNodesSortedByCard.begin(),determinableNodesSortedByCard.end(),ParentCardinalityCompare());

      // this is not the most efficient way of doing this, but it is
      // only run 1X so probably not that crucial to optimize.
      set <RV*> determinableNodesToRemove;
      for (unsigned j=0; 
	   j < determinableNodesSortedByCard.size() &&
	     PackCliqueValue::numWordsRequiredFor(hashableNodes,determinableNodesSortedByCard[j]) == hash_pack_size;
	   j++) {
	hashableNodes.push_back(determinableNodesSortedByCard[j]);
	determinableNodesToRemove.insert(determinableNodesSortedByCard[j]);
      }

      if (determinableNodesToRemove.size() > 0) {
	// need to create a new determinableNodes without the ones we removed, but otherwise in the same order.
	infoMsg(High,"Placing %d deterministic variables back into hashable set since word size (%d) the same\n",
		determinableNodesToRemove.size(),hash_pack_size);
	vector <RV*> updatedDeterminableNodes;
	for (unsigned j=0;j<determinableNodes.size();j++) {
	  if (determinableNodesToRemove.find(determinableNodes[j]) == determinableNodesToRemove.end())
	    updatedDeterminableNodes.push_back(determinableNodes[j]);
	}
	determinableNodes = updatedDeterminableNodes;
      }
    }
  }

  // setup and re-construct packer
  assert (packer.unPackedLen() == 0); // make sure it is empty.

  // If we have no hidden variables in this clique, than it means that
  // the entire clique consists of observed values. Since we do not
  // store observed values in clique tables, we need to do a bit of
  // extra checking here.

  if (hashableNodes.size() > 0) {
    new (&packer) PackCliqueValue(hashableNodes);
    // ensure that we have something to store.
    assert (packer.packedLen() > 0);
  } else {
    // fully observed clique, or at least a clique where all nodes are determinable.
    // Note, this might also include a node that is "hidden", but that has
    // a DeterminsticCPT with no parents that always returns one particular value.
  }

  // TODO: optimize initial size and growth factor.  Compute an
  // estimate of the state space of this clique for a starting
  // allocation of the value holder.  Take 1/4 of the estimated weight
  // for now.
  // TODO: note that this might FPE if weight() is too large, since
  // we get larger than pow() range. Fix this.
  //   allocationUnitChunkSize =
  //     (unsigned)(::pow(10,weight() - ::log10(4.0)));
  // lower bound.
  // if (allocationUnitChunkSize < 16)
  // allocationUnitChunkSize = 16;

  // set to small size now to test out re-allocation schemes.
  // allocationUnitChunkSize = 1;
  // allocationUnitChunkSize = 10000;

  if (packer.packedLen() > IMC_NWWOH) {
    // setup value hodler
    new (&valueHolder) CliqueValueHolder(packer.packedLen());

    // set up common clique hash tables 
    // TODO: add appropriate default staring hash sizes.
    // new (&cliqueValueHashSet) vhash_set< unsigned > (packer.packedLen(),2);
    new (&cliqueValueHashSet) vhash_set< unsigned > (packer.packedLen(),spaceMgrStartingSize);
#ifdef USE_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL
    temporaryCliqueValuePool.resize(spaceMgrStartingSize*packer.packedLen());
#endif

  } else {
    // then no need to do a hash table at all, either
    // we just store the packed values in a local integer
    // (in the case when there are hidden variables in the
    // clique) or we store nothing at all (in the case when
    // there are only observed variables).
  }

  // allocate storage for the non-recursive clique iteration code.
  probArrayStorage.resize( sortedAssignedNodes.size() + 1 );

}


/*-
 *-----------------------------------------------------------------------
 * MaxClique::computeSortedAssignedNodesDispositions()
 *   
 *   computes the disposition of the set of nodes that are sorted and
 *   assigned to this clique (i.e., that have the chance to use cpt
 *   iteration). Note that 'assigned' nodes are the ones that exist in
 *   the clique with their parents, but this does not mean that those
 *   nodes contribute probability to this clique. Even if they don't
 *   contribute probability, if the nodes are sparse, it can be
 *   beneficial to use that sparsity to consider child values rather
 *   than iterating over child cardinality (as would be done in
 *   ceIterateUnassignedIteratedNodes()). Note that even some assigned
 *   nodes, however, will be cardinality-based iterated, namely those
 *   that are both not probability contributers to this clique and
 *   also that are dense. Moreover, some nodes are not sorted assigned
 *   even if they are assigned, namely those that are contained in the
 *   separator and that are also not probability nodes (namely
 *   disposition AN_CONTINUE). In other words, we should never see
 *   disposition AN_CONTINUE nodes in this routine since we assume
 *   they have already been removed.  In any event, this routine
 *   decides how each node in the clique is iterated.
 *
 *
 * Preconditions:
 *   assignedNodes, sortedAssignedNodes, and cumulativeUnassignedIteratedNodes members 
 *   must have already been computed.
 *
 * Postconditions:
 *   dispositionSortedAssignedNodes is now set appropriately according to the
 *   nodes to be iterated in the clique as specified by sortedAssignedNodes.
 *
 * Side Effects:
 *     potentially modifies one member variable
 *
 * Results:
 *     Return true if there is at least one node that has a disposition that
 *     corresponds to doing no useful work (so it should ideally be removed).
 *
 *
 *-----------------------------------------------------------------------
 */
bool
MaxClique::computeSortedAssignedNodesDispositions()
{

  bool res = false;

  dispositionSortedAssignedNodes.resize(sortedAssignedNodes.size());
  for (unsigned i=0;i<sortedAssignedNodes.size();i++) {
    RV*rv = sortedAssignedNodes[i];
    if (!rv->hidden()) {
      // observed
      if (assignedProbNodes.find(rv) != assignedProbNodes.end()) {
	dispositionSortedAssignedNodes[i] = AN_CPT_ITERATION_COMPUTE_AND_APPLY_PROB;
      } else {
	if (rv->discrete()) {
	  // discrete observed variable. Prune here if sparse.
	  DiscRV* drv = (DiscRV*)rv;	  
	  if (drv->sparse())
	    dispositionSortedAssignedNodes[i] = AN_CPT_ITERATION_COMPUTE_PROB_REMOVE_ZEROS;
	  else {
	    dispositionSortedAssignedNodes[i] = AN_CONTINUE;
	    res = true;
	  }
	} else {
	  // continuous observed variable, but we get prob. in another clique, and
	  // we don't want to do this again, so just continue here.
	  dispositionSortedAssignedNodes[i] = AN_CONTINUE;
	  res = true;
	}
      }      
    } else {
      // hidden node, must be discrete
      DiscRV* drv = (DiscRV*)rv;
      if (unionIncommingCESeps.find(rv) == unionIncommingCESeps.end()) {
	// not contained in any CE incomming separator
	if (assignedProbNodes.find(rv) != assignedProbNodes.end()) {
	  // probabilities from this node contribute to clique potential
	  dispositionSortedAssignedNodes[i] = AN_CPT_ITERATION_COMPUTE_AND_APPLY_PROB;
	} else {
	  if (drv->sparse()) {
	    dispositionSortedAssignedNodes[i] = AN_CPT_ITERATION_COMPUTE_PROB_REMOVE_ZEROS;
	  } else {
	    dispositionSortedAssignedNodes[i] = AN_CARD_ITERATION;
	  }
	}
      } else {
	// contained in an CE incomming separator
	if (assignedProbNodes.find(rv) != assignedProbNodes.end()) {
	  // probabilities from this node contribute to clique potential
	  dispositionSortedAssignedNodes[i] = AN_COMPUTE_AND_APPLY_PROB;
	} else {
	  if (cumulativeAssignedNodes.find(rv) != cumulativeAssignedNodes.end()) {
	    dispositionSortedAssignedNodes[i] = AN_CONTINUE;
	    res = true;
	  } else if (!drv->sparse()) {
	    dispositionSortedAssignedNodes[i] = AN_CONTINUE;
	    res = true;
	  } else {
	    dispositionSortedAssignedNodes[i] = AN_CONTINUE_COMPUTE_PROB_REMOVE_ZEROS;
	  }
	}
      }
    }
  }
  return res;
}




/*-
 *-----------------------------------------------------------------------
 * MaxClique::sortAndAssignDispositions()
 *   
 *   sorts the nodes and assignes dispositions. Also, optinally removes
 *   nodes that would have disposition AN_CONTINUE since such nodes
 *   have no reason to be considered by the assigned nodes iteration (i.e., 
 *   they are set by the incomming separators and either have their probabilities
 *   assigned already, or do not need to be).
 *
 *   Note that clique driven inference must have all nodes sorted
 *   (including AN_CONTINUE nodes) since that form is not separator driven.
 *
 *
 * Preconditions:
 *   assignedNodes, sortedAssignedNodes, and cumulativeUnassignedIteratedNodes members 
 *   must have already been computed. Note, sortedAssignedNodes needs to be set
 *   but it doesn't need to be set to any "good" value.
 *
 * Postconditions:
 *   sortedAssignedNodes is probably re-set according to the desired sort order.
 *   dispositionSortedAssignedNodes is now set appropriately according to the
 *   nodes to be iterated in the clique.
 *
 * Side Effects:
 *     potentially modifies member variables
 *
 * Results:
 *     none
 *
 *
 *-----------------------------------------------------------------------
 */
void
MaxClique::sortAndAssignDispositions(const char *varCliqueAssignmentPrior)
{

  if (!varCliqueAssignmentPrior || strlen(varCliqueAssignmentPrior) == 0) {
    // then in this case, we don't try to remove any AN_CONTINUE nodes
    // nor do we do any sorting.  TODO: change things so that there
    // are never any AN_CONTINUE dispositions regardless of if we sort
    // or not.
    computeSortedAssignedNodesDispositions();
  } else {

    // the first thing we do is compute the node disposistions. We don't
    // actually use the disposistions here, as we use this routine just
    // to tell us if there are any AN_CONTINUE nodes in the clique, and
    // if there are, we remove them (since they are wasteful). The
    // variable nodesToRemove is a boolean telling us if there
    // are any such nodes to remove.
    bool nodesToRemove =  computeSortedAssignedNodesDispositions();


    if (nodesToRemove == false) {
      // just sort and leave.
      GraphicalModel::topologicalSortWPriority(assignedNodes,
					       assignedNodes,
					       sortedAssignedNodes,
					       varCliqueAssignmentPrior);
      // potentially need to re-compute disposistions if sorting order
      // changed anything. These are the final and real dispositions
      // used.
      nodesToRemove = computeSortedAssignedNodesDispositions();
      assert (!nodesToRemove);
    } else {

      // so there are nodes to remove, we go through remove the
      // nodes, re-sort, and then recompute the dispositions (yes, a bit
      // wasteful but this happens only once).

      set<RV*> toSortNodes;
      for (unsigned i=0;i<dispositionSortedAssignedNodes.size();i++) {
	if (dispositionSortedAssignedNodes[i] != AN_CONTINUE) {
	  toSortNodes.insert(sortedAssignedNodes[i]);
	}
      }
      GraphicalModel::topologicalSortWPriority(toSortNodes,
					       toSortNodes,
					       sortedAssignedNodes,
					       varCliqueAssignmentPrior);
      // now re-compute dispositions.
      nodesToRemove = computeSortedAssignedNodesDispositions();
      assert (!nodesToRemove);
    }
  }

  // next, we compute the continuation scores. We need one more (+1)
  // continuation score than there are sorted assigned nodes because
  // we may want to do pruning at the very last node.
  sortedAssignedContinuationScores.resize(sortedAssignedNodes.size()+1);
  if (cliqueBeamContinuationHeuristic && cliqueBeamBuildBeam != (-LZERO) /* ie, default -cpbeam meaning no pruning */ ) {
    // TODO: do something better than just using global max value, like local max value
    // of the current rv.
    sortedAssignedContinuationScores.ptr[sortedAssignedNodes.size()] = 1.0;
    for (int i = (sortedAssignedContinuationScores.size()-2); i>=0; i--) {
      sortedAssignedContinuationScores.ptr[i]
	= sortedAssignedContinuationScores.ptr[i+1] * 
	sortedAssignedNodes[i]->maxValue();
    }
  } else {
    for (unsigned i=0;i<sortedAssignedContinuationScores.size();i++)
      sortedAssignedContinuationScores.ptr[i] = 1.0;
  }

  //   for (unsigned i=0;i<sortedAssignedContinuationScores.size();i++) {
  //     printf("sortedAssignedContinuationScores[%d] = %f\n",i,
  // 	   sortedAssignedContinuationScores.ptr[i].valref());
  //   }

}




/*-
 *-----------------------------------------------------------------------
 * MaxClique::printAllJTInfo()
 *   
 *   prints everything JT (not not inference) about this clique to file.
 *
 * Preconditions:
 *   all variables must have been set up.  prepareForUnrolling() must have
 *   been called.
 *
 * Postconditions:
 *   none
 *
 * Side Effects:
 *     none
 *
 * Results:
 *     none
 *
 *
 *-----------------------------------------------------------------------
 */
void
MaxClique::printAllJTInfo(FILE*f,const unsigned indent,const set<RV*>& unassignedInPartition,
			  const bool upperBound,
			  const bool moreConservative,
			  const bool useDeterminism,
			  vector< set<RV*> > *lp_nodes,
			  vector< set<RV*> > *rp_nodes)
{

  // TODO: also print out nubmer of bits for acc and rem.

  psp(f,indent*2);
  fprintf(f,"Clique information: %d packed bits, %d unsigned words (%d splits), weight = %f, jt_weight = %f\n",
	  packer.packedLenBits(),packer.packedLen(),packer.numSplits(),weight(),
	  weightInJunctionTree(unassignedInPartition,
			       upperBound,
			       moreConservative,
			       useDeterminism,
			       lp_nodes,rp_nodes));


  psp(f,indent*2);
  fprintf(f,"%ld Nodes: ",(unsigned long)nodes.size()); printRVSetAndCards(f,nodes);

  psp(f,indent*2);
  fprintf(f,"%ld Assigned: ",(unsigned long)assignedNodes.size()); printRVSet(f,assignedNodes);

  psp(f,indent*2);
  fprintf(f,"%ld Assigned Sorted: ",(unsigned long)sortedAssignedNodes.size()); printRVSetAndCards(f,sortedAssignedNodes);

  psp(f,indent*2);
  fprintf(f,"%d Dispositions:",dispositionSortedAssignedNodes.size());
  for (unsigned i=0;i<dispositionSortedAssignedNodes.size();i++)
    fprintf(f," %d",dispositionSortedAssignedNodes[i]);
  fprintf(f,"\n");


  psp(f,indent*2);
  fprintf(f,"%ld Assigned Prob: ",(unsigned long)assignedProbNodes.size()); printRVSet(f,assignedProbNodes);  

  psp(f,indent*2);
  fprintf(f,"%ld Cum Assigned Prob: ",(unsigned long)cumulativeAssignedProbNodes.size()); printRVSet(f,cumulativeAssignedProbNodes);  

  psp(f,indent*2);
  fprintf(f,"%ld Union Incomming Seps: ",(unsigned long)unionIncommingCESeps.size()); printRVSet(f,unionIncommingCESeps);

  psp(f,indent*2);
  fprintf(f,"%ld Unassigned Iterated: ",(unsigned long)unassignedIteratedNodes.size()); printRVSet(f,unassignedIteratedNodes);


  psp(f,indent*2);
  fprintf(f,"%ld Cumulative Unassigned: ",(unsigned long)cumulativeUnassignedIteratedNodes.size()); printRVSet(f,cumulativeUnassignedIteratedNodes);

  if (hiddenNodes.size() == hashableNodes.size()) {
    psp(f,indent*2);
    fprintf(f,"%ld Hidden/Hashable: ",(unsigned long)hiddenNodes.size()); printRVSetAndCards(f,hiddenNodes);
  } else {
    psp(f,indent*2);
    fprintf(f,"%ld Hidden: ",(unsigned long)hiddenNodes.size()); printRVSetAndCards(f,hiddenNodes);

    psp(f,indent*2);
    fprintf(f,"%ld Hashable: ",(unsigned long)hashableNodes.size()); printRVSetAndCards(f,hashableNodes);
  }


  psp(f,indent*2);
  fprintf(f,"%ld Clique Neighbors: ",(unsigned long)neighbors.size());
  for (unsigned i=0;i<neighbors.size();i++) fprintf(f,"%d,",neighbors[i]); fprintf(f,"\n");

  psp(f,indent*2);
  fprintf(f,"%ld Clique Children: ",(unsigned long)children.size());
  for (unsigned i=0;i<children.size();i++) fprintf(f,"%d,",children[i]); fprintf(f,"\n");

  psp(f,indent*2);
  fprintf(f,"%ld Receive Seps: ",(unsigned long)ceReceiveSeparators.size());
  for (unsigned i=0;i<ceReceiveSeparators.size();i++) fprintf(f,"%d,",ceReceiveSeparators[i]); fprintf(f,"\n");

  psp(f,indent*2);
  fprintf(f,"%ld incomming VE Separators\n",(unsigned long)veSeparators.size());

  psp(f,indent*2);
  fprintf(f,"Send Sep: %d\n",ceSendSeparator);

}



/*-
 *-----------------------------------------------------------------------
 * MaxClique::computeUnassignedCliqueNodes()
 *   Compute the unassignedNodes variable
 *
 * Preconditions:
 *   Nodes and assignedNodes must be filled in.
 *   Note that this routine is a nop if we're doing separator driven inference.
 *
 * Postconditions:
 *   Unassigned nodes udpated.
 *
 * Side Effects:
 *   changes member variable unasssignedNodes
 *
 * Results:
 *     nothing
 *
 *-----------------------------------------------------------------------
 */
void 
MaxClique::computeUnassignedCliqueNodes()
{
  unassignedNodes.clear();
  set_difference(nodes.begin(),nodes.end(),
		 assignedNodes.begin(),assignedNodes.end(),
		 inserter(unassignedNodes,unassignedNodes.end()));
}


/*-
 *-----------------------------------------------------------------------
 * MaxClique::computeVESeparators()
 *   
 * computes initial information about any VE separators that
 * might exist in this clique, and sets up any data structures
 * needed to produce these separators. Returns the number of VE separators.
 * A VE separator is basically either:
 *   1) TYPE PC: When a child c exists such that:
 *        A: c is immediate observed (meaning it is const observed in the structure file)
 *        B: c a deterministic function of its parents, and the deterministic function
 *           is static (e.g., not an iterable DeterministicCPT or DT), and c is
 *           not switching (i.e., there are no switching parents, so only one CPT is involved)   
 *        C: all c's parents (and there are more than one) live in the current clique
 *        D: c gives probability to the current clique (so ideally is as
 *           far away from the current JT root as possible).
 *      - in this case the VE sep consists of the parents of c, and c.
 *   2) TYPE PCG: when a child c exists such that:
 *        A: c is immediate observed (meaning it is const observed in the structure file)
 *        B: c is a random or deterministic function of a c's single parent p,
 *           and c is not switching (i.e., there are no switching parents).
 *        C: p lives in the clique and is hidden
 *        D: p is a deterministic function of its parents, and the deterministic function
 *           is static (e.g., not an iterable DeterministicCPT or DT), and again there
 *           is no switching.
 *        C: all p's parents live in the current clique
 *        D: both c and p give probability to the current clique (so ideally as
 *           far away from the current JT root as possible).
 *      - in this case the VE sep consists of the parents of p, p, and c.
 *
 * TODO: switching() is possible, but we no longer can use just the DT to generate
 *       the table in those cases.
 *
 *
 * Preconditions:
 *   Nodes and assignedNodes must be filled in.
 *   Note that this routine is a nop if we're doing separator driven inference.
 *
 * Postconditions:
 *   Unassigned nodes udpated.
 *
 * Side Effects:
 *   changes some member variables relative to the VE seps.
 *
 * Results:
 *     returns number of possible VE seps.
 *
 *-----------------------------------------------------------------------
 */
unsigned
MaxClique::computeVESeparators()
{
  veSeparators.clear();

  // search through all assigned nodes looking for immediate observed nodes.

  // TODO: figure out a way to not store child value redundantly in separator.

  set<RV*>::iterator apn_it;
  for (apn_it = assignedProbNodes.begin(); apn_it != assignedProbNodes.end(); apn_it ++) {
    RV* rv = (*apn_it);
    if (rv->discreteObservedImmediate() && !rv->switching() && !rv->iterable()) {
      DiscRV* c = RV2DRV(rv);
      if (c->allParents.size() == 1 && 
	  c->allParents[0]->hidden() && 
	  c->allParents[0]->discrete() && 
	  RV2DRV(c->allParents[0])->deterministic() &&
	  !c->allParents[0]->iterable() &&
	  !c->allParents[0]->switching()) {
	// continue checking condition PCG
	RV * p = c->allParents[0];
	if ((SectionScheduler::useVESeparators & SectionScheduler::VESEP_PCG)  && 
	    (assignedProbNodes.find(p) !=  assignedProbNodes.end())) {

	  float logProdCard = 
	    log10((double)RV2DRV(p->allParents[0])->cardinality);
	  for (unsigned i=1;i<p->allParents.size();i++) {
	    if (!p->allParents[i]->discreteObservedImmediate())
	      logProdCard += log10((double)RV2DRV(p->allParents[i])->cardinality);
	  }

	  if (logProdCard > SeparatorClique::veSeparatorLogProdCardLimit)
	    continue;

	  infoMsg(Inference, Max,"Found VE sep of type PCG (currently not using it).\n");
	  // uncomment when rest of code below is finished.

	  // since p is assigned prob, all it parents live in the
	  // current clique. Iterate through all parents finding combinations
	  // that explain the child (with score given by Pr(c|p)).
	  VESepInfo vesep;
	  vesep.parents = p->allParents;
	  vesep.child = p;
	  vesep.grandChild = c;
	  veSeparators.push_back(vesep);
	}
      } else if ((SectionScheduler::useVESeparators & SectionScheduler::VESEP_PC) 
		 &&
		 (c->deterministic() && c->allParents.size() > 1)) {
	// continue checking condition PC.

	// NOTE: we allow the parents to be observed even if they are
	// not immediate observed, since in that case we generate the
	// table for all possible parent values.

	// check that the product of cardinalities of parents does not exceed threshold.


	float logProdCard = 
	  log10((double)RV2DRV(c->allParents[0])->cardinality);
	for (unsigned i=1;i<c->allParents.size();i++) {
	  if (!c->allParents[i]->discreteObservedImmediate())
	    logProdCard += log10((double)RV2DRV(c->allParents[i])->cardinality);
	}

	if (logProdCard > SeparatorClique::veSeparatorLogProdCardLimit)
	  continue;

	infoMsg(Inference, Max,"Found VE sep of type PC. iterable = %d\n",rv->iterable());

	// since c is an assigned prob node, we know its parents live
	// in current clique as well. We build a table of all parent
	// values that explain this child being observed this value,
	// and the table becomes a new separator used for all
	// inference instantiations of this clique.
	VESepInfo vesep;
	vesep.parents = c->allParents;
	vesep.child = c; 
	vesep.grandChild = NULL;
	veSeparators.push_back(vesep);

      }
    }
  }
  return veSeparators.size();
}



/////////////////////////////////////////////	
/// END OF FILE
/////////////////////////////////////////////
