/*-
 * GMTK_MaxCliqueBase.cc
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
#if HAVE_HG_H
#  include "hgstamp.h"
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

#include "GMTK_FileParser.h"
#include "GMTK_RV.h"
#include "GMTK_DiscRV.h"
#include "GMTK_HidDiscRV.h"
#include "GMTK_ObsDiscRV.h"
#include "GMTK_ContRV.h"
#include "GMTK_GMTemplate.h"
#include "GMTK_GMParms.h"
#include "GMTK_MaxCliqueBase.h"
#include "GMTK_MDCPT.h"
#include "GMTK_MSCPT.h"
#include "GMTK_MTCPT.h"
#include "GMTK_Mixture.h"
#include "GMTK_ObservationSource.h"
#include "GMTK_SectionScheduler.h"

VCID(HGID)


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


#if 0
    // TODO: remove this; it's in subclass?

// for sorting an array of CliqueValue descending based on the contained logpr 
struct CliqueValueDescendingProbCompare
{  
  bool operator() (const MaxCliqueTable::CliqueValue& cv1,
		   const MaxCliqueTable::CliqueValue& cv2)
  {
    return (cv1.p > cv2.p);
  }
};
#endif



/*
 *
 * Static memory for use for packed clique values. We do this here
 * rather than placing things on the stack since some routines might
 * do it repeatly and recursively. Also, making this global will make
 * the inference code below non-reentrant. TODO: change this when
 * getting working with POSIX threads.
 *
 */

namespace CliqueBuffer {
  // TODO: change this for multi-threading.

  // allocate some temporary storage for packed separator values.
  // 128 words is *much* bigger than any possible packed clique value
  // will take on, but it is easy/fast to allocate on the stack right now.
  unsigned packedVal[128];
}




/*
 *
 * Continuous observation per-feature penalty, default value defined here.
 *
 */
double 
MaxCliqueBase::continuousObservationPerFeaturePenalty = 0.0;


unsigned MaxCliqueBase::spaceMgrStartingSize = 1;
float    MaxCliqueBase::spaceMgrGrowthRate   = 1.05;
float    MaxCliqueBase::spaceMgrDecayRate    = 0.0;

bool MaxCliqueBase::storeDeterministicChildrenInClique = true;

double MaxCliqueBase::normalizeScoreEachClique = 1.0;

bool MaxCliqueBase::failOnZeroClique = true;


/*
 *
 * clique beam width, for clique-based beam pruning.  Default value is
 * very large (1.0/0.0 = est. of infty) meaning that we do no beam
 * pruning.
 *
 */
double
MaxCliqueBase::cliqueBeam=(-LZERO);

/*
 *
 * clique beam width, for use when building clique.  Default
 * value is very large (1.0/0.0 = est. of infty) meaning that we do no
 * beam pruning.
 *
 */
double
MaxCliqueBase::cliqueBeamBuildBeam=(-LZERO);

char* MaxCliqueBase::cliqueBeamBuildFilter = (char*)"";

bool MaxCliqueBase::cliqueBeamContinuationHeuristic = true;

double MaxCliqueBase::cliqueBeamBuildExpansionFactor = 1.0;
// could also just have a maximum beam value rather than max expansions.
unsigned MaxCliqueBase::cliqueBeamBuildMaxExpansions = 1;


unsigned MaxCliqueBase::cliqueBeamClusterPruningNumClusters = 0;
double MaxCliqueBase::cliqueBeamClusterBeam = (-LZERO);

unsigned MaxCliqueBase::cliqueBeamClusterMaxNumStates = NO_PRUNING_CLIQUEBEAMCLUSTERPRUNINGMAXSTATES;


/*
 *
 * The max number of states in a clique (or set to 0 to make ineffectual).
 * Basically, only the top cliqueBeamMaxNumStates cliuqe entries (ranked by their
 * probability entry) will be kept, everything below will be pruned away.
 *
 */
unsigned
MaxCliqueBase::cliqueBeamMaxNumStates = 0;

/*
 * Fraction of clique to retain. Default (1.0) means prune nothing.
 *
 */
float
MaxCliqueBase::cliqueBeamRetainFraction = 1.0;
// a version for cluster/diversity pruning
float
MaxCliqueBase::cliqueBeamClusterRetainFraction = 1.0;


/*
 * Fraction of clique mass to relinquish. Default (1.0) means prune nothing.
 *
 */
double
MaxCliqueBase::cliqueBeamMassRetainFraction = 1.0;

/*
 * When using cliqueBeamMassRetainFraction, this option determins the lower bound of the
 * min min clique state space size.
 *
 */
unsigned 
MaxCliqueBase::cliqueBeamMassMinSize = 1;



/*
 * When using mass beam pruning, we use this to further allow additional states through below
 * the one below first one after the mass has been acounted for. 0.0 to turn off.
 *
 */
double
MaxCliqueBase::cliqueBeamMassFurtherBeam = 0.0;



/*
 * When using mass beam pruning, we can exponentiate the clique scores 
 * to make them more uniform if need be, so as to hopefully prune more effectively.
 *
 */
double
MaxCliqueBase::cliqueBeamMassExponentiate = 1.0;

// A version of the above four variables that are used for cluster/diversity pruning.

double MaxCliqueBase::cliqueBeamClusterMassRetainFraction = 1.0;
unsigned MaxCliqueBase::cliqueBeamClusterMassMinSize = 1;
double MaxCliqueBase::cliqueBeamClusterMassFurtherBeam = 0.0;
double MaxCliqueBase::cliqueBeamClusterMassExponentiate = 1.0;


/*
 * If 0 <= val <= 1, then fraction of pruned clique table to uniformly re-sample from (without replacement)
 * if > 1, then number of times to re-sample (without replacement) from pruned clique table.
 * Set to 0.0 to turn off.
 */
double
MaxCliqueBase::cliqueBeamUniformSampleAmount = 0.0;


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        MaxClique support
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////



/*-
 *-----------------------------------------------------------------------
 * GMTemplate::makeComplete()
 *      make complete the set of random variables given
 *
 * Preconditions:
 *      - The set of random variables should be given.
 *
 * Postconditions:
 *      - set of random variables are made complete (via
 *        their neighbors variables)
 *
 * Side Effects:
 *      - random variables in rvs are changed.
 *
 * Results:
 *     nothing
 *
 *-----------------------------------------------------------------------
 */
void
MaxCliqueBase::makeComplete(const set<RV*> &rvs)
{
  // just go through each rv and union its neighbors
  // set with all of rvs.

  for (set<RV*>::iterator i=rvs.begin();
       i != rvs.end(); i++) {
    set<RV*> res;
    set_union(rvs.begin(),rvs.end(),
	      (*i)->neighbors.begin(),
	      (*i)->neighbors.end(),
	      inserter(res,res.end()));
    // make sure self is not its own neighbor
    res.erase((*i));
    (*i)->neighbors = res;
  }

  return;
}





/*-
 *-----------------------------------------------------------------------
 * MaxCliqueBase::computeWeight()
 *   Computes the log base 10 weight of a set of nodes (i.e.,
 *   the union of 'node' and 'nodes', ignores 'node' if 'node == NULL').
 *   
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
 *     the weight
 *
 *
 *-----------------------------------------------------------------------
 */
float
MaxCliqueBase::
computeWeight(const set<RV*>& nodes,
	      const RV* node,
	      const bool useDeterminism)
{
  // compute weight in log base 10 so that
  //   1) we don't overflow
  //   2) base 10 is an easy to understand magnitude rating of state space.

  float tmp_weight = 0;
  // First get cardinality of 'node', but if
  // it is continuous or observed, it does not change the weight.
  // TODO: The assumption here (for now) is that all continuous variables
  // are observed. This will change in a future version.
  if (node != NULL) {
    if (node->discrete() && node->hidden()) {
      DiscRV *const drv = (DiscRV*)node;
      // weight changes only if node is not deterministic (Lauritzen CG inference).
      if (useDeterminism && drv->sparse()) {
	// then there is a possibility that this node does not affect
	// the state space, as long as all of this nodes parents are
	// in the clique.  The variable 'truly_sparse' indicates that.
	bool truly_sparse = true;
	for (unsigned i=0;i<drv->allParents.size();i++) {
	  if (nodes.find(drv->allParents[i]) == nodes.end()) {
	    // found a parent of drv that is not in 'nodes' set so the
	    // node would not truly be sparse/deterministic here.
	    truly_sparse = false;
	    break;
	  }
	}
	if (truly_sparse)
	  tmp_weight += log10((double)drv->useCardinality());	
	else 
	  tmp_weight += log10((double)drv->cardinality);
      } else
	tmp_weight += log10((double)drv->cardinality);
    } else if (!node->discrete()) {
      // node is continuous observed.
      ContRV *crv = (ContRV*)node;
      tmp_weight += crv->dimensionality()*continuousObservationPerFeaturePenalty;
    }
  }
  // Next, get weight of all 'nodes'
  for (set<RV*>::iterator j=nodes.begin();
       j != nodes.end();
       j++) {
    RV *const rv = (*j);
    // First get cardinality of 'node', but if
    // it is continuous or observed, it does not change the weight.
    // TODO: The assumption here (for now) is that all continuous variables
    // are observed. This will change in a future version (Lauritzen CG inference).
    if (rv->discrete() && rv->hidden()) {
      DiscRV *const drv = (DiscRV*)rv;
      if (useDeterminism && drv->sparse()) {
	// then there is a possibility that this node does not affect
	// the state space, as long as all of this nodes parents are
	// in the clique.  The variable 'truly_sparse' indicates that.
	bool truly_sparse = true;
	for (unsigned i=0;i<drv->allParents.size();i++) {
	  if ((nodes.find(drv->allParents[i]) == nodes.end())
	      &&
	      (drv->allParents[i] != node)) {
	    // found a parent that is not in 'node' set so the node
	    // would not truly be sparse/deterministic here.
	    truly_sparse = false;
	    break;
	  }
	}
	if (truly_sparse)
	  tmp_weight += log10((double)drv->useCardinality());	
	else 
	  tmp_weight += log10((double)drv->cardinality);
      } else
	tmp_weight += log10((double)drv->cardinality);
    } else if (!rv->discrete()) {
      // node is continuous observed.
      ContRV *crv = (ContRV*)rv;
      tmp_weight += crv->dimensionality()*continuousObservationPerFeaturePenalty;
    } 
  }
  return tmp_weight;
}



/*-
 *-----------------------------------------------------------------------
 * TODO: remove this routine since it is not being used.
 *
 * MaxCliqueBase::computeWeightWithExclusion()
 *   Computes the log base 10 weight of a set of nodes (i.e.,
 *   the union of 'node' and 'nodes', ignores 'node' if 'node == NULL').
 *
 *   Note, that all nodes in the exclusion set (excludeSet) are counted
 *   as their full cardinality, regardless of if their parents are
 *   in the current node set or not and regardless of if the rv
 *   is sparse.
 *
 * Preconditions:
 *   Set of nodes must be valid meaning that it has valid neighbors,
 *   parents, and children member variables.
 *
 * Postconditions:
 *   computed weight (with exclusion) is provided.
 *
 * Side Effects:
 *     none
 *
 * Results:
 *     the weight
 *
 *
 *-----------------------------------------------------------------------
 */
float
MaxCliqueBase::
computeWeightWithExclusion(const set<RV*>& nodes,
			   const set<RV*>& unassignedIteratedNodes,
			   const set<RV*>& unionSepNodes,
			   const bool useDeterminism)
{
  // compute weight in log base 10 so that
  //   1) we don't overflow
  //   2) base 10 is an easy to understand magnitude rating of state space.

  float tmp_weight = 0;
  // Next, get weight of all 'nodes'
  for (set<RV*>::iterator j=nodes.begin();
       j != nodes.end();
       j++) {
    RV *const rv = (*j);
    // First get cardinality of 'node', but if
    // it is continuous or observed, it does not change the weight.
    // TODO: The assumption here (for now) is that all continuous variables
    // are observed. This will change in a future version (Lauritzen CG inference).
    if (rv->discrete() && rv->hidden()) {
      DiscRV *const drv = (DiscRV*)rv;
      bool truly_sparse = true;
      if (!useDeterminism || !drv->sparse()) {
	truly_sparse = false;
      } else { 
	// the potential exists for finding sparsity
	if (unassignedIteratedNodes.find(rv) == unassignedIteratedNodes.end() && 
	    unionSepNodes.find(rv) == unionSepNodes.end()) {
	  // then there is a possibility that this node does not affect
	  // the state space, as long as all of this nodes parents are
	  // in the clique.  The variable 'truly_sparse' indicates that.
	  for (unsigned i=0;i<drv->allParents.size();i++) {
	    if (nodes.find(drv->allParents[i]) == nodes.end()) {
	      // found a parent that is not in 'node' set so the node
	      // would not truly be sparse/deterministic here.
	      truly_sparse = false;
	      break;
	    }
	  }	  
	} else
	  truly_sparse = false;
      }
      // Finally, multiply in the weight depending on if it is "truly
      // sparse" or not.
      if (truly_sparse)
	tmp_weight += log10((double)drv->useCardinality());	
      else 
	tmp_weight += log10((double)drv->cardinality);
    } else if (!rv->discrete()) {
      // node is continuous observed.
      ContRV *crv = (ContRV*)rv;
      tmp_weight += crv->dimensionality()*continuousObservationPerFeaturePenalty;
    } 
  }
  return tmp_weight;
}





/*-
 *-----------------------------------------------------------------------
 * TODO: update this routine once new inference is in place.
 *
 * MaxCliqueBase::computeWeightInJunctionTree()
 *   Computes an ESTIMATE of the log base 10 weight of a set of nodes (i.e.,
 *   the union of 'node' and 'nodes', ignores 'node' if 'node == NULL').
 *
 *   This routine approximates the weight as the clique appears
 *   in a junction tree given a set of arguments. The arguments
 are as follows:

 nodes: The set of all nodes in clique.
 assigned_nodes: The set of nodes that are assigned (i.e., meaning
 node and its parents live in this clique).
 unassigned_iterated_nodes: nodes in this clique that   
 are not assigned and are not incomming
 separator nodes (so they need to be iterated in
 full, based in their cardinality).
 separator_nodes: Nodes that are part of an incomming separator
 during the collect evidence stage of inference. 
 The cost of these nodes (if sparse) depends on two things.

 If the separator nodes are also assigned here, we
 pay full cardinality (since we're doing separator
 driven clique instantiation). This is true even
 if the node are assigned (have parents in clique)
 since we can't be sure that the node in a
 separator has all its parents in the same separator.

 If, on the other hand, the separator nodes are
 previously unassigned, we also pay full card.

 If, the separator nodes were previously assigned
 we pay only use_card (i.e., if they are deterministic
 we pay nothing). This could produce a lower bound for
 weight, but we really want to bias in favor of such
 JTs, since in the DT case, this will significantly
 prune away zero probabilities.

 cumulativeAssignedNodes: all nodes lower (farther away from root) in JT that 
 are assigned. NOTE: These nodes INCLUDE the
 assigned nodes in this clique, meaning:
 assigned_nodes <= cumulativeAssignedNodes

 Here is the basic algorithm (acting as comments to the below):

 for each node v
 if hidden(v) 
 if (!sparse(v))
 multiply by card since not sparse.
 else if unassigned_iterated(v)
 multiply by card since we iterate over all values.
 else if separator_node(v)
 if assigned_node(v)
 multiply by card since wasn't assigned before and
 we need to iterate over all separator values.
 else
 if cumulativeAssignedNodes(v)
 if parents(v) <= separator_nodes
 multiply by use_card, since all of v's 
 parents are in the separator, and v will 
 cost nothing.
 else
 multiply by use_card (assigned earlier in JT).
 While it might be that this node costs more than
 use_card here, we are not penalizing for this case
 since we assume that the node will only come into
 this clique with parent values (if any) such that
 p(node|parents) > 0. This will in some cases
 give us a lower bound on the weight though.
 else
 multiply by card, since in this case it must
 be comming in as a previous unassignediterated node.
 This is a bad case.
 else // must be assigned non sep node
 multiply by use_card, since we know parents live in cur 
 clique since assigned here.
 *
 * Preconditions:
 *   All sets nodes must be valid meaning that it has valid neighbors,
 *   parents, and children member variables. All input arguments must
 *   refer to the same overall set of unrolled RVs.
 *
 * Postconditions:
 *   computed weight (with exclusion) is provided.
 *
 * Side Effects:
 *     none
 *
 * Results:
 *     the weight
 *
 *
 *-----------------------------------------------------------------------
 */

float
MaxCliqueBase
::computeWeightInJunctionTree(const set<RV*>& nodes,
			      const set<RV*>& assignedNodes,
			      const set<RV*>& cumulativeAssignedNodes,
			      const set<RV*>& unassignedIteratedNodes,
			      const set<RV*>& cumulativeUnassignedIteratedNodes,
			      const set<RV*>& separatorNodes,
			      const set<RV*>& unassignedInPartition,
			      vector< set<RV*> > *lp_nodes,
			      vector< set<RV*> > *rp_nodes,
			      const bool upperBound,
			      const bool moreConservative,
			      const bool useDeterminism)
{

  // Tue Mar 21 19:30:58 2006: comment this out for now until new
  // inference code is done.
  return 1.0;

#if 0

  // compute weight in log base 10 so that
  //   1) we don't overflow
  //   2) base 10 is an easy to understand magnitude rating of state space.

  // printf("separatorNodes.size() = %d, lp_nodes = 0x%X, rp_nodes = 0x%X\n",separatorNodes.size(),
  // 	 lp_nodes,rp_nodes);
  //   if (lp_nodes) {
  //     printf("lp_nodes = ");
  //     printRVSetPtr(stdout,(*lp_nodes));
  //     printf("\n");
  //   }
  //   if (rp_nodes) {
  //     printf("rp_nodes = ");
  //     printRVSetPtr(stdout,(*rp_nodes));
  //     printf("\n");
  //   }

  // weight for the sparse separator nodes
  float weight_sep_sparse = 0;
  // weight for the dense separator nodes
  float weight_sep_dense = 0;
  // weight for all the rest.
  float weight_remainder = 0;

  // Next, get weight of all 'nodes'
  for (set<RV*>::iterator j=nodes.begin();
       j != nodes.end();
       j++) {
    RV *const rv = (*j);

    //     printf("computing charge for RV %s(%d)=0x%X\n",
    // 	   rv->name().c_str(),
    // 	   rv->frame(),
    // 	   (unsigned)rv);

    // First get cardinality of 'node', but if
    // it is continuous or observed, it does not change the weight.
    // TODO: The assumption here (for now) is that all continuous variables
    // are observed. This will change in a future version (Lauritzen CG inference).
    if (rv->discrete() && rv->hidden()) {
      DiscRV *const drv = (DiscRV*)rv;
      if (!drv->sparse()) {
	// printf("   RV %s(%d) is dense\n",rv->name().c_str(),rv->frame());
	// node is dense.
	if (separatorNodes.find(rv) != separatorNodes.end()) {
	  // then node lives in separator.
	  weight_sep_dense += log10((double)drv->cardinality);
	} else {
	  weight_remainder += log10((double)drv->cardinality);
	}
      } else {
	// printf("   RV %s(%d) is sparse\n",rv->name().c_str(),rv->frame());
	// node is sparse
	if (!useDeterminism) {
	  // we're not using determinism/sparsity for node charging
	  if (separatorNodes.find(rv) != separatorNodes.end()) {
	    // printf("   RV %s(%d) charged to weight_sep_sparse, full since not using det.\n",rv->name().c_str(),rv->frame());
	    // then node lives in separator.
	    weight_sep_sparse += log10((double)drv->cardinality);
	  } else {
	    // printf("   RV %s(%d) charged to weight_remainder, full since not using det.\n",rv->name().c_str(),rv->frame());
	    weight_remainder += log10((double)drv->cardinality);
	  }
	} else {
	  // we are using determinism/sparsity for node charging.
	  if ((separatorNodes.find(rv) != separatorNodes.end())
	      ||
	      // if this is not a separator node, but if we are given the the nodes in
	      // the left partition and rv is in the left partition, then it is still a
	      // separator node.
	      ((lp_nodes != NULL) && ((*lp_nodes).find(rv) != (*lp_nodes).end()))) {
	    // printf("   RV %s(%d) is a separator\n",rv->name().c_str(),rv->frame());
	    // node in separator case.
	    if (unassignedInPartition.find(rv) != unassignedInPartition.end()) {
	      // node in separator, UNassigned in the current partition (so either
	      // assigned in previous or next partition, depending on direction
	      // of edges):
	      if (upperBound)
		weight_sep_sparse += log10((double)drv->cardinality);
	      else {
		// Then variable is unasigned in this partition. We assume
		// that the node has been assigned in another (say
		// previous) partition, but it could have been assigned in
		// the next partition (if there are backwards time links).

		if (moreConservative) {
		  // now you know that more conservative really means try a cludge/hack.
		  // (min(card,unavailble_parents_prod_card) + use_card)/2 (but in log domain).
		  // Even so, this could either be a lower or upper bound.
		  weight_sep_sparse += (min(drv->log10ProductCardOfParentsNotContainedInSet(separatorNodes),
					    log10((double)drv->cardinality)) + log10((double)drv->useCardinality()))/2.0;
		} else {
		  // What we do: we don't charge full amount. Note that this
		  // could cause the estimate to be LOWER than the true
		  // weight.
		  weight_sep_sparse += log10((double)drv->useCardinality());
		}
	      }
	    } else if (cumulativeUnassignedIteratedNodes.find(rv) !=
		       cumulativeUnassignedIteratedNodes.end()) {
	      // node in separator, assigned in this partition, but NOT assigned
	      // in any previous clique.
	      if (assignedNodes.find(rv) == assignedNodes.end()) {
		// node in separator, assigned in this partition, unassigned previously, not assigned 
		// in current clique either:

		// Charge full amount since we do separator iteration over
		// something that is not assigned in any previous cliques
		// and nor in this clique.

		// add to weight_remainder since there is no
		// way this will be pruned coming into this clique.
		weight_remainder += log10((double)drv->cardinality);
	      } else {
		// node in separator, assigned in this partition, unassigned
		// previously, assigned here in this clique:

		// This is the case we would like to avoid since for
		// separator driven iteration, we are iterating over all
		// values of var, and don't remove the zeros until this
		// clique. If there are many of these cases, we might
		// consider doing clique driven rather than separator
		// driven clique potential instantiation.

		// While it will come into this clique without zeros being
		// removed, this clique will remove them (since it is
		// assigned), so from a memory point of view, we could
		// charge useCard. For now, we are conservative here,
		// however, and charge full card (which is the computational
		// but not the memory cost).

		// add to weight_remainder since there is no
		// way this will be pruned coming into this clique.
		weight_remainder += log10((double)drv->cardinality);
	      }
	    } else {
	      // node in separator, assigned in this partition, and assigned in a
	      // previous clique.

	      if (upperBound) 
		weight_sep_sparse += log10((double)drv->cardinality);
	      else {
		// This is an important one since this is quite likely to occur
		// when determinism abounds.
		if (moreConservative) {
		  // now you know that more conservative really means try a cludge/hack.
		  // (min(card,unavailble_parents_prod_card) + use_card)/2 (but in log domain).
		  // Even so, this could either be a lower or upper bound.
		  weight_sep_sparse += (min(drv->log10ProductCardOfParentsNotContainedInSet(separatorNodes),
					    log10((double)drv->cardinality)) + log10((double)drv->useCardinality()))/2.0;
		} else {
		  // Charge low amount since it has been assigned in some
		  // previous clique, and at least one of the separators will
		  // kill off the zero prob entries. This could cause the
		  // estimate to be LOWER than the true weight.
		  weight_sep_sparse += log10((double)drv->useCardinality());
		}
	      }
	    }
	  } else {
	    // printf("   RV %s(%d) is NOT a separator\n",rv->name().c_str(),rv->frame());
	    // node NOT in separator case.
	    if (unassignedInPartition.find(rv) != unassignedInPartition.end()) {
	      // node NOT in separator, unassigned in partition (probably assigned
	      // in right partition).

	      if (upperBound) 
		weight_remainder += log10((double)drv->cardinality);
	      else {
		if (moreConservative) {
		  // (min(card,unavailble_parents_prod_card) + use_card)/2 (but in log domain).
		  // Even so, this could either be a lower or upper bound.
		  weight_remainder += (min(drv->log10ProductCardOfParentsNotContainedInSet(nodes),
					   log10((double)drv->cardinality)) + log10((double)drv->useCardinality()))/2.0;	      
		} else {
		  // Then unasigned in this partition. We assume that the
		  // node has been assigned in another (say previous) partition
		  // and we don't charge full amount. Note that
		  // this could cause the estimate to be LOWER than
		  // the true weight.
		  weight_remainder += log10((double)drv->useCardinality());
		}
	      }
	    } else if (assignedNodes.find(rv) != assignedNodes.end()) {
	      // node NOT in separator, assigned here in this clique:

	      // Then assigned in this clique. Charge correct amount.
	      weight_remainder += log10((double)drv->useCardinality());
	    } else {
	      // node NOT in separator, not assigned here in this clique:

	      // Not assigned in this clique. We know it can't be in
	      // cumulativeAssignedNodes since it is not a sep node.
	      assert ( cumulativeAssignedNodes.find(rv) == cumulativeAssignedNodes.end());
	      // charge full amount.
	      weight_remainder += log10((double)drv->cardinality);
	    }
	  }
	}
      }
    } else if (!rv->discrete()) {
      // node is continuous observed.
      ContRV *crv = (ContRV*)rv;
      weight_remainder += crv->dimensionality()*continuousObservationPerFeaturePenalty;
    } else {
      // node is discrete observed, charge nothing.
    }
  }

  //   printf("weight_sep_sparse = %f, weight_sep_dense = %f, weight_remainder = %f\n",
  // 	  weight_sep_sparse,
  // 	  weight_sep_dense,
  // 	  weight_remainder);

  return SectionScheduler::jtWeightSparseNodeSepScale*weight_sep_sparse + 
    SectionScheduler::jtWeightDenseNodeSepScale*weight_sep_dense 
    + weight_remainder;
#endif
}






/*-
 *-----------------------------------------------------------------------
 * MaxCliqueBase::reportScoreStats()
 *
 *    Report score stats for this clique.
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
MaxCliqueBase::
reportScoreStats()
{
  for (set<RV*>::iterator j=nodes.begin();
       j != nodes.end();
       j++) {
    RV *const rv = (*j);
    logpr p = rv->maxValue();
    printf("max val of rv %s(%d) = %f\n",
	   rv->name().c_str(),rv->frame(),p.val());
  }
}




/*-
 *-----------------------------------------------------------------------
 * MaxCliqueBase::printCliqueNodes()
 *   
 *   prints the names and frames of the clique nodes.
 *
 * Preconditions:
 *   all variables must have been set up.
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
MaxCliqueBase::printCliqueNodes(FILE*f)
{
  fprintf(f,"%ld Nodes: ",(unsigned long)nodes.size()); printRVSet(f,nodes);
}



void
MaxCliqueBase::numParentsSatisfyingChild(unsigned& num,unsigned par,vector <RV*> & parents, RV* child)
{
  if (par == parents.size()) {
    logpr cur_p;
    child->probGivenParents(cur_p);
    if (!cur_p.zero())
      num++;
#if 0
    printf("Pr[%s(%d)=%d|",child->name().c_str(),child->frame(),
	   RV2DRV(child)->val);
    printRVSetAndValues(stdout,parents,false);
    printf("]=%d\n",!cur_p.zero());
#endif
  } else {
    DiscRV*drv = RV2DRV(parents[par]);
    if (drv->hidden()) {
      HidDiscRV*hdrv = (HidDiscRV*)drv;
      for (hdrv->val = 0; hdrv->val < hdrv->cardinality; hdrv->val ++) {
	numParentsSatisfyingChild(num,par+1,parents,child);
      }
    } else {
      ObsDiscRV*odrv = (ObsDiscRV*)drv;
      odrv->setToObservedValue();
      numParentsSatisfyingChild(num,par+1,parents,child);
    }
  }
}



/*
  A note on trace printing and debugging values to control it for the
  following few routines: 
  Rough range of debug values:

  Med  = 50,
  High = 60,
  Huge = 70,
  Mega = 80,
  Giga = 90,

  We have three routines for iterating a clique.
  ceIterateSeparators
  ceIterateUnassignedIteratedNodes
  ceIterateAssignedNodesRecurse
  - High: print just final clique insertions, nothign else. No indentation. 
  done just in ceIterateAssignedNodesRecurse
  - High+5: add all starts (Separator, Unassigned, and RV iteration starts)
  - Huge: Print all iters (separator, unassigned, & RV iters), but not parent values in RV case.
  - Mega: Also print all parent values at all iterations.
  - Mega+5: also prints continuous observation values (rather than just "=C").

*/

