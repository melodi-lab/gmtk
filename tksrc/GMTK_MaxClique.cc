/*-
 * GMTK_MaxClique.cc
 *
 *     maxClique support and JT probabilistic inference. Includes the
 *     implementation of a message in the message passing algorithm.
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
#include "GMTK_ContRV.h"
#include "GMTK_GMTemplate.h"
#include "GMTK_GMParms.h"
#include "GMTK_MDCPT.h"
#include "GMTK_MSCPT.h"
#include "GMTK_MTCPT.h"
#include "GMTK_Mixture.h"
#include "GMTK_ObservationMatrix.h"
#include "GMTK_JunctionTree.h"

VCID("$Header$");


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Constants
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


// set to 1 to test
#define CLIQUE_VALUE_HOLDER_STARTING_SIZE 23
#define CLIQUE_VALUE_HOLDER_GROWTH_RATE   2.0

#define AI_SEP_VALUE_HOLDER_STARTING_SIZE 23
#define AI_SEP_VALUE_HOLDER_GROWTH_RATE   2.0

#define REM_SEP_VALUE_HOLDER_STARTING_SIZE 23
#define REM_SEP_VALUE_HOLDER_GROWTH_RATE   2.0


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Static variables and functions
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////





/*
 *
 * Continuous observation per-feature penalty, default value defined here.
 *
 */
double 
MaxClique::continuousObservationPerFeaturePenalty = 0.0;

/*
 *
 * true to do separator driven inference, false means do clique driven. 
 *
 */
bool
MaxClique::ceSeparatorDrivenInference = true;

bool MaxClique::perSegmentClearCliqueValueCache = true;

/*
 *
 * clique beam width, for clique-based beam pruning.  Default value is
 * very large (1.0/0.0 = est. of infty) meaning that we do no beam
 * pruning.
 *
 */
double
MaxClique::cliqueBeam=(-LZERO);

/*
 *
 * The max number of states in a clique (or set to 0 to make ineffectual).
 * Basically, only the top cliqueBeamMaxNumStates cliuqe entries (ranked by their
 * probability entry) will be kept, everything below will be pruned away.
 *
 */
unsigned
MaxClique::cliqueBeamMaxNumStates = 0;

/*
 * Fraction of clique to retain. Default (1.0) means prune nothing.
 *
 */
float
MaxClique::cliqueBeamRetainFraction = 1.0;



/*
 *
 * separator beam width, for separator-based beam pruning.  Default value is
 * very large (1.0/0.0 = est. of infty) meaning that we do no beam
 * pruning.
 *
 */
double
SeparatorClique::separatorBeam=(-LZERO);



// TODO: put this in misc support
static void
psp(FILE*f,const int numSpaceChars,const char c = ' ')
{
  int tmp = numSpaceChars;
  while (tmp--)
    fprintf(f,"%c",c);
}



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

  :  cliqueValueSpaceManager(1,     // starting size
			     2.0,   // growth rate
			     1,     // growth addition
			     0.90)  // decay rate 
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

}





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
MaxClique::makeComplete(const set<RV*> &rvs)
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
 * MaxClique::computeWeight()
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
MaxClique::
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
 * MaxClique::computeWeightWithExclusion()
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
MaxClique::
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
 * MaxClique::computeWeightInJunctionTree()
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
MaxClique
::computeWeightInJunctionTree(const set<RV*>& nodes,
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
			      const bool useDeterminism)
{
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
  
  return JunctionTree::jtWeightSparseNodeSepScale*weight_sep_sparse + 
    JunctionTree::jtWeightDenseNodeSepScale*weight_sep_dense 
    + weight_remainder;
}



#if 0
      // see comments above for description and rational of this algorithm
      if (!useDeterminism || !drv->sparse()) {
	// charge full amount since not sparse.
	tmp_weight += log10((double)drv->cardinality);
      } else if (separatorNodes.find(rv) != separatorNodes.end()) {
	// separator node case.
	if (unassignedInPartition.find(rv) != unassignedInPartition.end()) {
	  // separator, UNassigned in the current partition (so either
	  // assigned in previous or next partition, depending on direction
	  // of edges):

	  if (upperBound)
	    tmp_weight += log10((double)drv->cardinality);
	  else {
	    // Then variable is unasigned in this partition. We assume
	    // that the node has been assigned in another (say
	    // previous) partition, but it could have been assigned in
	    // the next partition (if there are backwards time links).

	    if (moreConservative) {
	      // now you know that more conservative really means try a cludge/hack.
	      // (min(card,unavailble_parents_prod_card) + use_card)/2 (but in log domain).
	      // Even so, this could either be a lower or upper bound.
	      tmp_weight += (min(drv->log10ProductCardOfParentsNotContainedInSet(separatorNodes),
				 log10((double)drv->cardinality)) + log10((double)drv->useCardinality()))/2.0;
	    } else {
	      // What we do: we don't charge full amount. Note that this
	      // could cause the estimate to be LOWER than the true
	      // weight.
	      tmp_weight += log10((double)drv->useCardinality());
	    }
	  }
	} else if (cumulativeUnassignedIteratedNodes.find(rv) !=
		   cumulativeUnassignedIteratedNodes.end()) {
	  // separator, assigned in this partition, but NOT assigned
	  // in any previous clique.
	  if (assignedNodes.find(rv) == assignedNodes.end()) {
	    // separator, assigned in this partition, unassigned previously, not assigned 
	    // in current clique either:

	    // Charge full amount since we do separator iteration over
	    // something that is not assigned in any previous cliques
	    // and nor in this clique.
	    tmp_weight += log10((double)drv->cardinality);
	  } else {
	    // separator, assigned in this partition, unassigned
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
	    tmp_weight += log10((double)drv->cardinality);
	  }
	} else {
	  // separator, assigned in this partition, and assigned in a
	  // previous clique.

	  if (upperBound) 
	    tmp_weight += log10((double)drv->cardinality);
	  else {
	    // This is an important one since this is quite likely to occur
	    // when determinism abounds.
	    if (moreConservative) {
	      // now you know that more conservative really means try a cludge/hack.
	      // (min(card,unavailble_parents_prod_card) + use_card)/2 (but in log domain).
	      // Even so, this could either be a lower or upper bound.
	      tmp_weight += (min(drv->log10ProductCardOfParentsNotContainedInSet(separatorNodes),
				 log10((double)drv->cardinality)) + log10((double)drv->useCardinality()))/2.0;
	    } else {
	      // Charge low amount since it has been assigned in some
	      // previous clique, and at least one of the separators will
	      // kill off the zero prob entries. This could cause the
	      // estimate to be LOWER than the true weight.
	      tmp_weight += log10((double)drv->useCardinality());
	    }
	  }
	}
      } else {
	// non separator node case.
	if (unassignedInPartition.find(rv) != unassignedInPartition.end()) {
	  // non separator, unassigned in partition.

	  if (upperBound) 
	    tmp_weight += log10((double)drv->cardinality);
	  else {
	    if (moreConservative) {
	      // (min(card,unavailble_parents_prod_card) + use_card)/2 (but in log domain).
	      // Even so, this could either be a lower or upper bound.
	      tmp_weight += (min(drv->log10ProductCardOfParentsNotContainedInSet(nodes),
				 log10((double)drv->cardinality)) + log10((double)drv->useCardinality()))/2.0;	      
	    } else {
	      // Then unasigned in this partition. We assume that the
	      // node has been assigned in another (say previous) partition
	      // and we don't charge full amount. Note that
	      // this could cause the estimate to be LOWER than
	      // the true weight.
	      tmp_weight += log10((double)drv->useCardinality());
	    }
	  }
	} else if (assignedNodes.find(rv) != assignedNodes.end()) {
	  // non separator, assigned here:

	  // Then assigned in this clique. Charge correct amount.
	  tmp_weight += log10((double)drv->useCardinality());
	} else {
	  // non separator, not assigned here:

	  // Not assigned in this clique. We know it can't be in
	  // cumulativeAssignedNodes since it is not a sep node.
	  assert ( cumulativeAssignedNodes.find(rv) == cumulativeAssignedNodes.end());
	  // charge full amount.
	  tmp_weight += log10((double)drv->cardinality);
	}
      }
    } else if (!rv->discrete()) {
      // node is continuous observed.
      ContRV *crv = (ContRV*)rv;
      tmp_weight += crv->dimensionality()*continuousObservationPerFeaturePenalty;
    } else {
      // node is discrete observed.
      
    } 
  }
#endif


#if 0
float
MaxClique::
computeWeightInJunctionTree(const set<RV*>& nodes,
			    const set<RV*>& assignedNodes,
			    const set<RV*>& unassignedIteratedNodes,
			    const set<RV*>& separatorNodes,
			    const set<RV*>& cumulativeAssignedNodes,
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
      // see comments above for description and rational of this algorithm
      if (!useDeterminism || !drv->sparse()) {
	truly_sparse = false;
      } else if (unassignedIteratedNodes.find(rv) != unassignedIteratedNodes.end()) {
	truly_sparse = false;
      } else if (separatorNodes.find(rv) != separatorNodes.end()) {
	if (assignedNodes.find(rv) != assignedNodes.end()) {
	  truly_sparse = false;	  
	} else {
	  if (cumulativeAssignedNodes.find(rv) != cumulativeAssignedNodes.end()) {
	    if (rv->allParentsContainedInSet(separatorNodes)) {
	      ; // do nothing, we can multiply by use_card.
	    } else {
	      // TODO: This is the problem case here. we need to estimate the
	      // cost.
	      truly_sparse = true;;
	    }
	  } else {
	    truly_sparse = false;
	  }
	}
      } else {
	; // do nothing.
      }
      // Finally, multiply in the weight depending on if it is "truly
      // sparse" or not.
      if (truly_sparse)
	tmp_weight += log10((double)drv->useCardinality());	
      else 
	tmp_weight += log10((double)drv->cardinality);
    }
  }
  return tmp_weight;
}
#endif







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

  // set up number of hidden ndoes
  set<RV*>::iterator it;
  for (it = nodes.begin();
       it != nodes.end();
       it++) {
    RV* rv = (*it);
    if (rv->hidden())
      hiddenNodes.push_back(rv);
  }

  // setup and re-construct packer
  assert (packer.unPackedLen() == 0); // make sure it is empty.

  // If we have no hidden variables in this clique, than it means that
  // the entire clique consists of observed values. Since we do not
  // store observed values in clique tables, we need to do a bit of
  // extra checking here.

  if (hiddenNodes.size() > 0) {
    new (&packer) PackCliqueValue(hiddenNodes);
    // ensure that we have something to store.
    assert (packer.packedLen() > 0);
  } else {
    // fully observed clique.
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
    new (&valueHolder) CliqueValueHolder(packer.packedLen(),
					 CLIQUE_VALUE_HOLDER_STARTING_SIZE, // set to 1 to test.
					 CLIQUE_VALUE_HOLDER_GROWTH_RATE); // 1.25
    // set up common clique hash tables 
    // TODO: add appropriate default staring hash sizes.
    // new (&cliqueValueHashSet) vhash_set< unsigned > (packer.packedLen(),2);
    new (&cliqueValueHashSet) vhash_set< unsigned > (packer.packedLen(),CLIQUE_VALUE_HOLDER_STARTING_SIZE);
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
 * MaxClique::computeAssignedNodesDispositions()
 *   
 *   computes the number and set of nodes that are assigned to this clique
 *   that we are actually to iterate over (rather than have them be iterated
 *   by a separator driven iteration).
 *
 * Preconditions:
 *   assignedNodes, sortedAssignedNodes, and cumulativeUnassignedIteratedNodes members 
 *   must have already been computed.
 *
 * Postconditions:
 *   If we iterate over all nodes, then iterateSortedAssignedNodesP is of size 0.
 *   If we iterate over some nodes, then iterateSortedAssignedNodesP is of same size
 *        as assignedNodes and indicates which nodes to iterate over.
 *
 * Side Effects:
 *     potentially modifies one member variable
 *
 * Results:
 *     none
 *
 *
 *-----------------------------------------------------------------------
 */
void
MaxClique::computeAssignedNodesDispositions()
{


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
	  else 
	    dispositionSortedAssignedNodes[i] = AN_CONTINUE;
	} else {
	  // continuous observed variable, but we get prob. in another clique, and
	  // we don't want to do this again, so just continue here.
	  dispositionSortedAssignedNodes[i] = AN_CONTINUE;
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
	  } else if (!drv->sparse()) {
	    dispositionSortedAssignedNodes[i] = AN_CONTINUE;
	  } else {
	    dispositionSortedAssignedNodes[i] = AN_CONTINUE_COMPUTE_PROB_REMOVE_ZEROS;
	  }
	}
      }
    }
  }
  
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
			  set<RV*>* lp_nodes,set<RV*>* rp_nodes)
{

  // TODO: also print out nubmer of bits for acc and rem.

  psp(f,indent*2);
  fprintf(f,"Clique information: %d packed bits, %d unsigned words, weight = %f, jt_weight = %f\n",
	  packer.packedLenBits(),packer.packedLen(),weight(),
	  weightInJunctionTree(unassignedInPartition,
			       upperBound,
			       moreConservative,
			       useDeterminism,
			       lp_nodes,rp_nodes));


  psp(f,indent*2);
  fprintf(f,"%d Nodes: ",nodes.size()); printRVSet(f,nodes);

  psp(f,indent*2);
  fprintf(f,"%d Assigned: ",assignedNodes.size()); printRVSet(f,assignedNodes);

  psp(f,indent*2);
  fprintf(f,"%d Assigned Sorted: ",sortedAssignedNodes.size()); printRVSet(f,sortedAssignedNodes);

  psp(f,indent*2);
  fprintf(f,"%d Dispositions:",dispositionSortedAssignedNodes.size());
  for (unsigned i=0;i<dispositionSortedAssignedNodes.size();i++)
    fprintf(f," %d",dispositionSortedAssignedNodes[i]);
  fprintf(f,"\n");


  psp(f,indent*2);
  fprintf(f,"%d Assigned Prob: ",assignedProbNodes.size()); printRVSet(f,assignedProbNodes);  

  psp(f,indent*2);
  fprintf(f,"%d Cum Assigned Prob: ",cumulativeAssignedProbNodes.size()); printRVSet(f,cumulativeAssignedProbNodes);  

  psp(f,indent*2);
  fprintf(f,"%d Union Incomming Seps: ",unionIncommingCESeps.size()); printRVSet(f,unionIncommingCESeps);

  psp(f,indent*2);
  fprintf(f,"%d Unassigned Iterated: ",unassignedIteratedNodes.size()); printRVSet(f,unassignedIteratedNodes);



  psp(f,indent*2);
  fprintf(f,"%d Cumulative Unassigned: ",cumulativeUnassignedIteratedNodes.size()); printRVSet(f,cumulativeUnassignedIteratedNodes);


  psp(f,indent*2);
  fprintf(f,"%d Hidden: ",hiddenNodes.size()); printRVSet(f,hiddenNodes);


  psp(f,indent*2);
  fprintf(f,"%d Clique Neighbors: ",neighbors.size());
  for (unsigned i=0;i<neighbors.size();i++) fprintf(f,"%d,",neighbors[i]); fprintf(f,"\n");

  psp(f,indent*2);
  fprintf(f,"%d Clique Children: ",children.size());
  for (unsigned i=0;i<children.size();i++) fprintf(f,"%d,",children[i]); fprintf(f,"\n");

  psp(f,indent*2);
  fprintf(f,"%d Receive Seps: ",ceReceiveSeparators.size());
  for (unsigned i=0;i<ceReceiveSeparators.size();i++) fprintf(f,"%d,",ceReceiveSeparators[i]); fprintf(f,"\n");

  psp(f,indent*2);
  fprintf(f,"Send Sep: %d\n",ceSendSeparator);

}



/*-
 *-----------------------------------------------------------------------
 * MaxClique::printCliqueNodes()
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
MaxClique::printCliqueNodes(FILE*f)
{
  fprintf(f,"%d Nodes: ",nodes.size()); printRVSet(f,nodes);
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
  if (ceSeparatorDrivenInference)
    return; // don't compute if we're doing separator driven.
  unassignedNodes.clear();
  set_difference(nodes.begin(),nodes.end(),
		 assignedNodes.begin(),assignedNodes.end(),
		 inserter(unassignedNodes,unassignedNodes.end()));
}


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        InferenceMaxClique support
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////




/*-
 *-----------------------------------------------------------------------
 * InferenceMaxClique::InferenceMaxClique()
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
InferenceMaxClique::InferenceMaxClique(MaxClique& from_clique,
				       vector <RV*>& newRvs,
				       map < RVInfo::rvParent, unsigned >& ppf,
				       const unsigned int frameDelta)
  : origin(from_clique)

{

  set<RV*>::iterator it;

  // clone over nodes RVs.
  fNodes.resize(origin.nodes.size());
  unsigned i=0;
  for (it = origin.nodes.begin();
       it != origin.nodes.end();
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
    fNodes[i++] = nrv;
  }

  // and clone over assigned nodes and sorted assigned nodes
  fSortedAssignedNodes.resize(origin.sortedAssignedNodes.size());
  for (i=0;i<origin.sortedAssignedNodes.size();i++) {
    RV* rv = origin.sortedAssignedNodes[i];
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    

    if ( ppf.find(rvp) == ppf.end() ) {
      coredump("ERROR: can't find assigned rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
    }

    RV* nrv = newRvs[ppf[rvp]];
    fSortedAssignedNodes[i] = nrv;
  }

  // do unassignedIteratedNodes
  i=0;
  fUnassignedIteratedNodes.resize(origin.unassignedIteratedNodes.size());
  for (it = origin.unassignedIteratedNodes.begin();
       it != origin.unassignedIteratedNodes.end();
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
    fUnassignedIteratedNodes[i++] = nrv;
  }


  // do unassignedNodes
  if (MaxClique::ceSeparatorDrivenInference == false) {
    i=0;
    fUnassignedNodes.resize(origin.unassignedNodes.size());
    for (it = origin.unassignedNodes.begin();
	 it != origin.unassignedNodes.end();
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
      fUnassignedNodes[i++] = nrv;
    }
  }


  // Clique values only store/hash values of hidden (thus necessarily
  // discrete) variables since they are the only thing that change.
  discreteValuePtrs.resize(origin.hiddenNodes.size());
  for (i=0;i<discreteValuePtrs.size();i++) {
    // get the unrolled rv for this hidden node
    RV* rv = origin.hiddenNodes[i];
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    
    if ( ppf.find(rvp) == ppf.end() ) {
      warning("ERROR: can't find assigned rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
      assert ( ppf.find(rvp) != ppf.end() );
    }
    RV* nrv = newRvs[ppf[rvp]];

    // hidden nodes are always discrete (in this version).
    DiscRV* drv = 
      (DiscRV*)nrv;
    // grab a pointer directly to its value for easy access later.
    discreteValuePtrs[i] = &(drv->val);
  }

  numCliqueValuesUsed = 0;
  maxCEValue.set_to_zero();

  // TODO: optimize this and make depend on if clique is all hidden, has observed, etc.
  // NOTE: This must be set to something greater than 0.
  // cliqueValues.resize(3); // 10000
  cliqueValues.resize(origin.cliqueValueSpaceManager.currentSize());

}



//
// TODO: make proper comments to all of the functions below.
//

/*-
 *-----------------------------------------------------------------------
 * MaxClique::ceGatherFromIncommingSeparators()
 *
 *    Collect Evidence, Gather from Incomming Separators: This routine
 *    is the main driver for creating a clique table during the
 *    collect evidence stage of inference. It assumes that all
 *    separators have been created, based on that the clique will be created.
 *    There are two forms of clique table creation here, the one that
 *    is used is decided at the top of the routine.
 *    
 *        1) clique driven clique creation: Here, we iterate through
 *        all clique values of a cliuqe (based on the unassigned and
 *        assigned nodes in that clique), and for each one, check all
 *        incomming separators to make sure that the separators have
 *        entries for the corresponding intersected variables at the
 *        current clique value. The separator scores are then
 *        multiplied together with the scoring assigned variables,
 *        all of which becomes the clique value score. The clique
 *        value and its score are assigned to the clique.
 *
 *        2) separator driven clique instantiation: Here, the
 *        separators are iterated directly, and later separators are
 *        iterated over only those values which are "compatible" with
 *        earlier separators (meaning any zero probability entries are
 *        skipped entirely, this is done using the separator data
 *        structures which were set up prior to this routine). Once an
 *        entry survives each separator (meaning the separator
 *        intersection is non zero for that entry), it is subjected to
 *        assigned clique nodes, and if it survies them as well
 *        (meaning not zero prob), then the score is multiplied and
 *        added to the clique.
 *
 *        Note that this version first iterates through separators,
 *        then through unassigned clique nodes, and finally assigned
 *        clique nodes.
 *
 *    Inference method 1 is good because normal weight tells you
 *    exactly how computationally costly it is. Inference method 1 is
 *    bad because it is apparentlly *much* slower then method 2 when
 *    we have either many deterministic variables, or 2 when there is
 *    much pruning. Unfortunately, the computational cost of method 2
 *    is much more difficult to predict a-priori since it depends on
 *    pruning, and the degree of sparse/deterministic CPTs. Junction
 *    Tree weight (jtweight) is an attempt to predict this cost.
 *
 * Preconditions:
 *
 *   All the cliques incomming separators *must* have been created and are ready
 *   to be used to produce the clique table.
 *     
 *
 * Postconditions:
 *    Clique table has been created.
 *
 * Side Effects:
 *    Will significantly affect member variables in cliques.
 *
 * Results:
 *     nothing
 *
 *-----------------------------------------------------------------------
 */
void
InferenceMaxClique::ceGatherFromIncommingSeparators(JT_InferencePartition& part)
{
  if (!origin.ceSeparatorDrivenInference)
    return ceGatherFromIncommingSeparatorsCliqueDriven(part);

  if (origin.hiddenNodes.size() == 0) {
    ceGatherFromIncommingSeparatorsCliqueObserved(part);
  } else {
    // if we're still here, we do regular separator driven inference.
    logpr p = 1.0;
    if (origin.ceReceiveSeparators.size() == 0) {
      if (origin.unassignedIteratedNodes.size() == 0) {
	ceIterateAssignedNodes(part,0,p);
      } else {
	ceIterateUnassignedIteratedNodes(part,0,p);
      }
    } else {
      ceIterateSeparators(part,0,p);
    }
  }
}



/*-
 *-----------------------------------------------------------------------
 * MaxClique::ceIterateSeparators()
 *
 *    Collect Evidence, Iterate Separators: This routine
 *    is part of separator driven clique instantiation. It 
 *    iterates through each separator and moves on to the next
 *    separator checking it based on the accumulated set of variables
 *    that have been assigned in previous separators. This routine
 *    might therefore be called a 'sparse join', since it does a join
 *    of a bunch of separator tables each of which might be very sparse. 
 *
 *    If we have reached the last separator in this clique, we move
 *    directly on to the unassigned nodes (if any). 
 *
 * Preconditions:
 *
 *   All the cliques incomming separators *must* have been created and
 *   are ready to be used to produce the clique table. 
 *     
 *
 * Postconditions:
 *    We move on to the unasssigned nodes with a probabilty value that
 *    includes the multiplication of all compatible separator entries for this var
 *    value set.
 *
 * Side Effects:
 *    none
 *
 * Results:
 *     nothing
 *
 *-----------------------------------------------------------------------
 */
void
InferenceMaxClique::ceIterateSeparators(JT_InferencePartition& part,
					const unsigned sepNumber,
					const logpr p)
{
  if (sepNumber == origin.ceReceiveSeparators.size()) {
    ceIterateUnassignedIteratedNodes(part,0,p);
    return;
  }
  
  // get a handy reference to the current separator
  InferenceSeparatorClique& sep = 
    part.separatorCliques[origin.ceReceiveSeparators[sepNumber]];

  if (message(Giga)) {
    fprintf(stdout,"Starting separator iteration, sepNumber =%d, part sepNo = %d,p = %f, nodes:",
	    sepNumber,origin.ceReceiveSeparators[sepNumber],p.val());
    printRVSet(stdout,sep.fNodes);
  }


  unsigned sepValueNumber;  
  if (sep.origin.hAccumulatedIntersection.size() > 0) {
    // look up existing intersected values to see if we have a match
    // and only proceed if we do.
    
    // find index if it exists.
    // unsigned tmp[ISC_NWWOH_AI];
    unsigned tmp[((ISC_NWWOH_AI>1)?ISC_NWWOH_AI:1)];
    unsigned *key;
    // TODO: optimize this check out of loop (perhaps also assign to const local variable).
    if (sep.origin.accPacker.packedLen() <= ISC_NWWOH_AI) {
      key = &tmp[0];
    } else {
      // use since it always holds at least one extra slot ready
      key = sep.origin.accValueHolder.curCliqueValuePtr();
    }
    sep.origin.accPacker.pack((unsigned**)sep.accDiscreteValuePtrs.ptr,
			      (unsigned*)key);
    // key = &tmp;
    // const unsigned** foo = (const unsigned**)&key;
    // int** foo = sep.accDiscreteValuePtrs.ptr;
    // sep.origin.accPacker.pack(foo,&tmp);

    unsigned* indexp = sep.iAccHashMap.find(key);
    if (indexp == NULL) {
      // Then not found in this separator, so it must have zero
      // probability. We continue with the next value of the previous
      // separator.
      return;
    } else {
      // need to further iterate.
      sepValueNumber = *indexp;
    }

  } else {
    // This condition would fail if we've completely pruned away
    // parent clique, but that can only happen when
    // we've got a negative beam. non-negative beam means
    // that we will always have at least one entry.
    // TODO: probably ok to remove this assertion.
    assert ( sep.separatorValues.size() == 1);
    sepValueNumber = 0;
  }

  // Iterate through remainder of separator. 
  // NOTE: we could do some online pruning here as well, but instead
  // we do it in a special separator prune routine, called ceSeparatorPrune().
  if (sep.origin.hRemainder.size() == 0) {
    // Only one remainder entry (in position 0) and also no need to
    // unpack since all has been covered by accumulated intersection
    // set above in a previous separator. Just continue on with single
    // value. 
    // - 
    // Note that this case also includes the case when the entire
    // separator is observed, i.e., when we'll get an accum inter size
    // of 0 and a hreminder size of 0 -- nothing to unpack here either
    // (no hash tables even exist), so we just continue along.

    if (message(Giga)) {
      fprintf(stdout,"Separator iteration no-unpack, sepNumber =%d, part sepNo = %d,p=%f, sp=%f, nodes:",
	      sepNumber,origin.ceReceiveSeparators[sepNumber],
	      p.val(),
	      sep.separatorValues.ptr[sepValueNumber].remValues.ptr[0].p.val());
      printRVSetAndValues(stdout,sep.fNodes);
    }

    // We should have either one or zero entries. The only way zero
    // entries could arrise is if we have done severe separator pruning.
    assert ( (sep.separatorValues.ptr[sepValueNumber].remValues.size() & ~0x1)
	     == 0x0 );

    // assert ( sep.separatorValues.ptr[sepValueNumber].remValues.size() == 1 );
    if (sep.separatorValues.ptr[sepValueNumber].numRemValuesUsed == 1) {
      // Continue down with new probability value.
      // Search for tag 'ALLOCATE_REMVALUES_OPTION' in this file for
      // more info why remValues.ptr[0] exists.
      // Note: could do more separator pruning here.
      ceIterateSeparators(part,sepNumber+1,
			  p*
			  sep.separatorValues.ptr[sepValueNumber].remValues.ptr[0].p);
    }

  } else {

    // TODO: this assertion should be redundant (check above)
    assert ( sep.origin.remPacker.packedLen() > 0 );

    if (sep.origin.remPacker.packedLen() <= ISC_NWWOH_RM) {
      for (unsigned i=0;i< sep.separatorValues.ptr[sepValueNumber].numRemValuesUsed; i++) {

	// if (i == 0 && sep.fNodes[0]->timeIndex == 5 && sep.fNodes[0]->label == "state") {
	// printf("here1");
	// } 

	// TODO: optimize this, pre-compute base array outside of loop.
	sep.origin.remPacker.unpack(
		  (unsigned*)&(sep.separatorValues.ptr[sepValueNumber].remValues.ptr[i].val[0]),
		  (unsigned**)sep.remDiscreteValuePtrs.ptr);

	// if (i == 0 && sep.fNodes[0]->timeIndex == 5 && sep.fNodes[0]->label == "state") {
	// printf("here2");
	// } 

	if (message(Giga)) {
	  fprintf(stdout,"Separator iteration %d, sepNumber =%d, part sepNo = %d,p=%f, sp=%f, nodes:",
		  i,
		  sepNumber,origin.ceReceiveSeparators[sepNumber],
		  p.val(),
		  sep.separatorValues.ptr[sepValueNumber].remValues.ptr[i].p.val());
	  printRVSetAndValues(stdout,sep.fNodes);
	}

	// continue down with new probability value.
	// TODO: separator prune here.
	ceIterateSeparators(part,sepNumber+1,
		    p*
		    sep.separatorValues.ptr[sepValueNumber].remValues.ptr[i].p);
      }
    } else {
      for (unsigned i=0;i< sep.separatorValues.ptr[sepValueNumber].numRemValuesUsed; i++) {

	// TODO: optimize this, pre-compute base array outside of loop.
	sep.origin.remPacker.unpack(
		  (unsigned*)sep.separatorValues.ptr[sepValueNumber].remValues.ptr[i].ptr,
		  (unsigned**)sep.remDiscreteValuePtrs.ptr);

	if (message(Giga)) {
	  fprintf(stdout,"pSeparator iteration %d, sepNumber =%d, part sepNo = %d,p=%f, sp=%f, nodes:",
		  i,
		  sepNumber,origin.ceReceiveSeparators[sepNumber],
		  p.val(),
		  sep.separatorValues.ptr[sepValueNumber].remValues.ptr[i].p.val());
	  printRVSetAndValues(stdout,sep.fNodes);
	}

	// continue down with new probability value.
	// TODO: separator prune here.
	ceIterateSeparators(part,sepNumber+1,
			    p*
			    sep.separatorValues.ptr[sepValueNumber].remValues.ptr[i].p);
      }
    }
  }

}



/*-
 *-----------------------------------------------------------------------
 * MaxClique::ceIterateAssignedNodesRecurse()
 *
 *    Collect Evidence, Iterate Separators: This routine
 *    is part of separator driven clique instantiation. 
 *  
 *    This routine is a main workhorse of inference. Once we have a
 *    partial clique value that has survied the separators, and has
 *    (possibly) unassigned variable values, the next thing we need to
 *    do is iterate through the assigned nodes in this cliuqe.
 *
 *
 *    If we have reached the last unassigned node in this clique, we
 *    go ahead and assign the value to this clique if it is non-zero.
 *
 *    If we have not reached the last, we 'iterate' over the assigned
 *    node depending on its type (see MaxClique.h for a definition of
 *    the different types of assigned nodes).
 *    
 *    Note: this routine will *never* knowingly add a zero-probability
 *    clique value to the clique. I.e., all current zeros are not
 *    added. This way, during backward pass if anyway, we never need
 *    to check for divide by zero.
 *
 *
 * Preconditions:
 *
 *   All the cliques incomming separators *must* have been iterated, and
 *   same for the unassigned nodes (if any). The clique is assumed
 *   to have at least one hidden node in it.
 *
 * Postconditions:
 *    We move on to the asssigned nodes with an probabilty. The prob. has
 *    not been updated from what is passed in since these nodes are unassigned
 *    and so never contribute score to this clique.
 *
 * Side Effects:
 *    none
 *
 * Results:
 *     nothing
 *
 *-----------------------------------------------------------------------
 */
void
InferenceMaxClique::ceIterateAssignedNodesRecurse(JT_InferencePartition& part,
						  const unsigned nodeNumber,
						  const logpr p)
{

  if (nodeNumber == fSortedAssignedNodes.size()) {
    // time to store clique value and total probability, p is
    // current clique probability.

    // printf("ceIterateAssignedNodesRecurse: nodeNumber = %d, p = %f,",nodeNumber,p.val());
    // printRVSet(fNodes);

    // keep track of the max clique probability right here.
    if (p > maxCEValue)
      maxCEValue = p;

    if (numCliqueValuesUsed >= cliqueValues.size()) {
      // TODO: optimize this.
      if (numCliqueValuesUsed >= origin.cliqueValueSpaceManager.currentSize())
	origin.cliqueValueSpaceManager.advanceToNextSize();
      // cliqueValues.resizeAndCopy(cliqueValues.size()*2);
      cliqueValues.resizeAndCopy(origin.cliqueValueSpaceManager.currentSize());
    }

    // TODO: figure out if it is possible to get around doing this
    // check (or to help branch prediction predict it, since it will
    // be different for differnet cliques). Answer: We can do this
    // check once at beginning of iteration of assigned nodes, and
    // have two versions of this code.
    if (origin.packer.packedLen() <= IMC_NWWOH) {
      // pack the clique values directly into place
      origin.packer.pack(
			 (unsigned**)discreteValuePtrs.ptr,
			 (unsigned*)&(cliqueValues.ptr[numCliqueValuesUsed].val[0]));
    } else {
      // Deal with the hash table to re-use clique values.
      // First, grab pointer to storge where next clique value would
      // be stored if it ends up being used.
      unsigned *pcv = origin.valueHolder.curCliqueValuePtr();
      // Next, pack the clique values into this position.
      origin.packer.pack((unsigned**)discreteValuePtrs.ptr,(unsigned*)pcv);
      // Look it up in the hash table.
      bool foundp;
      unsigned *key;
      key = origin.cliqueValueHashSet.insert(pcv,foundp);
      if (!foundp) {
	// if it was not found, need to claim this storage that we
	// just used.
	origin.valueHolder.allocateCurCliqueValue();
      }
      // Save the pointer to whatever the hash table decided to use.
      cliqueValues.ptr[numCliqueValuesUsed].ptr = key;
    }
    // save the probability
    cliqueValues.ptr[numCliqueValuesUsed].p = p;
    numCliqueValuesUsed++;

    if (message(Mega)) {
      psp(stdout,2*nodeNumber);
      infoMsg(Mega,"%d:Inserting New Clique Value. prob = %f, sum = %f: ",
	      nodeNumber,
	      cliqueValues.ptr[numCliqueValuesUsed-1].p.val(),sumProbabilities().val());
      printRVSetAndValues(stdout,fNodes);
    }



    return;
  }
  RV* rv = fSortedAssignedNodes[nodeNumber];
  // do the loop right here

  if (message(Giga)) {
    psp(stdout,2*nodeNumber);
    infoMsg(Giga,"%d:Starting assigned iteration of rv %s(%d), p = %f\n",
	    nodeNumber,
	    rv->name().c_str(),rv->frame(),p.val());
  }

  switch (origin.dispositionSortedAssignedNodes[nodeNumber]) {
  case MaxClique::AN_CPT_ITERATION_COMPUTE_AND_APPLY_PROB: 
    {
      logpr cur_p;
      rv->begin(cur_p);
      do {
	if (message(Giga)) {
	  if (!rv->discrete()) {
	    psp(stdout,2*nodeNumber);
	    infoMsg(Giga,"%d:assigned iter/prob app rv %s(%d)=C,p=%f,cur_p=%f\n",nodeNumber,
		    rv->name().c_str(),rv->frame(),p.val(),cur_p.val());

	  } else {
	    DiscRV* drv = (DiscRV*)rv;
	    psp(stdout,2*nodeNumber);
	    infoMsg(Giga,"%d:assigned CPT iter/prob app rv %s(%d)=%d,p=%f,cur_p=%f\n",
		    nodeNumber,
		    rv->name().c_str(),rv->frame(),drv->val,p.val(),cur_p.val());
	  }
	  if (message(Max)) {
	    psp(stdout,2*nodeNumber);
	    infoMsg(Max,"%d:RV %s(%d)'s parents:",nodeNumber,rv->name().c_str(),rv->frame());
	    printRVSetAndValues(stdout,rv->allParents);
	  }
	}
	// if at any step, we get zero, then back out.
	if (!cur_p.essentially_zero()) {
	  // Continue, updating probability by cur_p, contributing this
	  // probability to the clique potential.
	  ceIterateAssignedNodesRecurse(part,nodeNumber+1,p*cur_p);
	}
      } while (rv->next(cur_p));
    }
    break;

  case MaxClique::AN_CPT_ITERATION_COMPUTE_PROB_REMOVE_ZEROS:
    {
      logpr cur_p;
      rv->begin(cur_p);
      do {
	// At each step, we compute probability
	if (message(Giga)) {
	  if (!rv->discrete()) {
	    psp(stdout,2*nodeNumber);
	    infoMsg(Giga,"%d:assigned iter rv %s(%d)=C,p=%f,cur_p=%f\n",
		    nodeNumber,
		    rv->name().c_str(),rv->frame(),p.val(),cur_p.val());
	  } else {
	    DiscRV* drv = (DiscRV*)rv;
	    psp(stdout,2*nodeNumber);
	    infoMsg(Giga,"%d:assigned CPT iter/zero removal rv %s(%d)=%d,p=%f,cur_p=%f\n",
		    nodeNumber,
		    rv->name().c_str(),rv->frame(),drv->val,p.val(),cur_p.val());
	  }
	  if (message(Max)) {
	    psp(stdout,2*nodeNumber);
	    infoMsg(Max,"%d:RV %s(%d)'s parents:",nodeNumber,rv->name().c_str(),rv->frame());
	    printRVSetAndValues(stdout,rv->allParents);
	  }
	}
	// if at any step, we get zero, then back out.
	if (!cur_p.essentially_zero()) {
	  // Continue, do not update probability!!
	  ceIterateAssignedNodesRecurse(part,nodeNumber+1,p);
	}
      } while (rv->next(cur_p));
    }
    break;


  case MaxClique::AN_CARD_ITERATION:
    {
      DiscRV* drv = (DiscRV*)rv;
      // do the loop right here
      drv->val = 0;
      do {
	if (message(Giga)) {
	  psp(stdout,2*nodeNumber);
	  infoMsg(Giga,"%d:assigned card iter rv %s(%d)=%d,p=%f\n",
		  nodeNumber,
		  rv->name().c_str(),rv->frame(),drv->val,p.val());
	  if (message(Max)) {
	    psp(stdout,2*nodeNumber);
	    infoMsg(Max,"%d:RV %s(%d)'s parents:",nodeNumber,rv->name().c_str(),rv->frame());
	    printRVSetAndValues(stdout,rv->allParents);
	  }
	}
	// Continue, do not update probability!!
	ceIterateAssignedNodesRecurse(part,nodeNumber+1,p);
      } while (++drv->val < drv->cardinality);
    }
    break;

  case MaxClique::AN_COMPUTE_AND_APPLY_PROB:
    {
      // TODO: Make more efficient version of this, based on the type of
      // RV.
      logpr cur_p = rv->probGivenParents();
      // if at any step, we get zero, then back out.
      if (message(Giga)) {
	psp(stdout,2*nodeNumber);
	if (!rv->discrete()) {
	  infoMsg(Giga,"%d:assigned compute/prob app rv %s(%d)=C,p=%f,cur_p=%f\n",
		  nodeNumber,
		  rv->name().c_str(),rv->frame(),p.val(),cur_p.val());
	} else {
	  DiscRV* drv = (DiscRV*)rv;
	  infoMsg(Giga,"%d:assigned compute/prob app rv %s(%d)=%d,p=%f,cur_p=%f\n",
		  nodeNumber,
		  rv->name().c_str(),rv->frame(),drv->val,p.val(),cur_p.val());
	}
	if (message(Max)) {
	  psp(stdout,2*nodeNumber);
	  infoMsg(Max,"%d:RV %s(%d)'s parents:",nodeNumber,rv->name().c_str(),rv->frame());
	  printRVSetAndValues(stdout,rv->allParents);
	}
      }
      if (!cur_p.essentially_zero()) {
	// Continue, updating probability by cur_p.
	ceIterateAssignedNodesRecurse(part,nodeNumber+1,p*cur_p);
      }
    }
    break;

  case MaxClique::AN_CONTINUE:
    ceIterateAssignedNodesRecurse(part,nodeNumber+1,p);
    break;

  case MaxClique::AN_CONTINUE_COMPUTE_PROB_REMOVE_ZEROS:
    {
      // TODO: Make more efficient version of this, based on the type of
      // RV.
      logpr cur_p = rv->probGivenParents();
      if (message(Giga)) {
	psp(stdout,2*nodeNumber);
	if (rv->discrete()) {
	  infoMsg(Giga,"%d:assigned compute continue rv %s(%d)=%d,p=%f,cur_p=%f\n",
		  nodeNumber,
		  rv->name().c_str(),rv->frame(),RV2DRV(rv)->val,p.val(),cur_p.val());
	} else {
	  infoMsg(Giga,"%d:assigned compute continue rv %s(%d)=C,p=%f,cur_p=%f\n",
		  nodeNumber,
		  rv->name().c_str(),rv->frame(),p.val(),cur_p.val());
	}
	if (message(Max)) {
	  psp(stdout,2*nodeNumber);
	  infoMsg(Max,"%d:RV %s(%d)'s parents:",nodeNumber,rv->name().c_str(),rv->frame());
	  printRVSetAndValues(stdout,rv->allParents);
	}
      }

      if (!cur_p.essentially_zero()) {
	// Continue, do not update probability!!
	ceIterateAssignedNodesRecurse(part,nodeNumber+1,p);
      }
    }
    break;

  default:
    assert(0);
    break;
  }

}



/*-
 *-----------------------------------------------------------------------
 * MaxClique::ceIterateAssignedNodesNoRecurse()
 *
 *    This is a non-recursive version of
 *    ceIterateAssignedNodesRecurse(). The reason for this routine is
 *    to try to produce code that has better branch prediction
 *    behavior and that also creates fewer temporary variables.
 *
 *    This routine is a bit more involved than the recursive version, as
 *    it uses customize loops (with some goto statements). Note that
 *    the recursive version is called as a backup when there are no
 *    assigned nodes in a clique. For that reason, and also for
 *    pedagogical reasons, the recursive version is still in place.
 *
 * Preconditions:
 *
 *   Same as ceIterateAssignedNodesRecurse()
 *
 * Postconditions:
 *
 *   Same as ceIterateAssignedNodesRecurse()
 *
 * Side Effects:
 *
 *   Same as ceIterateAssignedNodesRecurse()
 *
 * Results:
 *
 *   Same as ceIterateAssignedNodesRecurse()
 *
 *-----------------------------------------------------------------------
 */
void
InferenceMaxClique::ceIterateAssignedNodesNoRecurse(JT_InferencePartition& part,
						    const logpr p)
{

  if (fSortedAssignedNodes.size() == 0 || message(Giga)) {
    // let recursive version handle degenerate or message full case
    return ceIterateAssignedNodesRecurse(part,0,p);
  }

  // parray has to be 1 offset, storing p in entry -1
  logpr* parray = origin.probArrayStorage.ptr + 1;
  parray[-1] = p;

  int nodeNumber;
  logpr cur_p;
  for (nodeNumber=0;nodeNumber<(int)fSortedAssignedNodes.size();nodeNumber++) {
    RV* rv = fSortedAssignedNodes[nodeNumber];
    switch (origin.dispositionSortedAssignedNodes[nodeNumber]) {
    case MaxClique::AN_CPT_ITERATION_COMPUTE_AND_APPLY_PROB: 
      {
	rv->begin(cur_p);
	// check for possible zero (could occur with
	// zero score or observations).
	if (cur_p.essentially_zero()) {
	  // Since we have zero here, we cancel iterations of all
	  // subsequent variables right now, rather than iterate them
	  // with what will end up being zero probability.
	  parray[nodeNumber].set_to_zero();
	  nodeNumber++;
	  goto end_of_initial_clique_value;
	} else 
	  parray[nodeNumber] = parray[nodeNumber-1]*cur_p;
      }
      break;

    case MaxClique::AN_CPT_ITERATION_COMPUTE_PROB_REMOVE_ZEROS:
      {
	rv->begin(cur_p);
	// check for possible zero (could occur with
	// zero score or observations).
	if (cur_p.essentially_zero()) {
	  // Since we have zero here, we cancel iterations of all
	  // subsequent variables right now, rather than iterate them
	  // with what will end up being zero probability.
	  parray[nodeNumber].set_to_zero();
	  nodeNumber++;
	  goto end_of_initial_clique_value;
	} else 
	  parray[nodeNumber] = parray[nodeNumber-1];
      }
      break;


    case MaxClique::AN_CARD_ITERATION:
      {
	DiscRV* drv = (DiscRV*)rv;
	drv->val = 0;
	parray[nodeNumber] = parray[nodeNumber-1];
      }
      break;

    case MaxClique::AN_COMPUTE_AND_APPLY_PROB:
      {
	// TODO: Make more efficient version of this, based on the type of
	// RV.
	rv->probGivenParents(cur_p);
	// might get a zero probability, check condition here.
	if (cur_p.essentially_zero()) {
	  // Since we have zero here, we cancel iterations of all
	  // subsequent variables right now, rather than iterate them
	  // with what will end up being zero probability.
	  parray[nodeNumber].set_to_zero();
	  nodeNumber++;
	  goto end_of_initial_clique_value;
	} else 
	  parray[nodeNumber] = parray[nodeNumber-1]*cur_p;
      }
      break;

    case MaxClique::AN_CONTINUE:
      parray[nodeNumber] = parray[nodeNumber-1];
      break;

    case MaxClique::AN_CONTINUE_COMPUTE_PROB_REMOVE_ZEROS:
      {
	// TODO: Make more efficient version of this, based on the type of
	// RV.
	rv->probGivenParents(cur_p);
	// might get a zero probability check condition here.
	if (cur_p.essentially_zero()) {
	  // Since we have zero here, we cancel iterations of all
	  // subsequent variables right now, rather than iterate them
	  // with what will end up being zero probability.
	  parray[nodeNumber].set_to_zero();
	  nodeNumber++;
	  goto end_of_initial_clique_value;
	} else 
	  parray[nodeNumber] = parray[nodeNumber-1];
      }
      break;
    }
  }
  
 end_of_initial_clique_value:
  
  // get ready for main loop
  const bool imc_nwwoh_p = (origin.packer.packedLen() <= IMC_NWWOH);
  const int maxNodeNumber = (int)fSortedAssignedNodes.size()-1;
  nodeNumber--;

  MaxClique::AssignedNodeDisposition cur_disp = origin.dispositionSortedAssignedNodes.ptr[nodeNumber];
  RV* cur_rv = fSortedAssignedNodes.ptr[nodeNumber];

  // check if zero probability, and if so, skip the first one and
  // continue on.
  if (parray[nodeNumber].essentially_zero())
    goto next_iteration;


  // main loop, iterate through all assigned nodes in this clique.
  do {

    // add a clique value to the clique.
    {
      const logpr final_p = parray[nodeNumber];

      // time to store clique value and total probability, p is current
      // clique probability.
      // keep track of the max clique probability right here.

      if (final_p > maxCEValue)
	maxCEValue = final_p;

      if (numCliqueValuesUsed >= cliqueValues.size()) {
	// TODO: optimize this.
	// cliqueValues.resizeAndCopy(cliqueValues.size()*2);
	if (numCliqueValuesUsed >= origin.cliqueValueSpaceManager.currentSize())
	  origin.cliqueValueSpaceManager.advanceToNextSize();
	cliqueValues.resizeAndCopy(origin.cliqueValueSpaceManager.currentSize());
      }

      // Possibly remove this check if possible. It's probably ok though
      // since as long as this loop runs several times, branch
      // predictions should handle it (so 1 cycle latency).
      if (imc_nwwoh_p) {
	// pack the clique values directly into place
	origin.packer.pack(
			   (unsigned**)discreteValuePtrs.ptr,
			   (unsigned*)&(cliqueValues.ptr[numCliqueValuesUsed].val[0]));
      } else {
	// Deal with the hash table to re-use clique values.
	// First, grab pointer to storge where next clique value would
	// be stored if it ends up being used.
	unsigned *pcv = origin.valueHolder.curCliqueValuePtr();
	// Next, pack the clique values into this position.
	origin.packer.pack((unsigned**)discreteValuePtrs.ptr,(unsigned*)pcv);
	// Look it up in the hash table.
	bool foundp;
	unsigned *key;
	key = origin.cliqueValueHashSet.insert(pcv,foundp);
	if (!foundp) {
	  // if it was not found, need to claim this storage that we
	  // just used.
	  origin.valueHolder.allocateCurCliqueValue();
	}
	// Save the pointer to whatever the hash table decided to use.
	cliqueValues.ptr[numCliqueValuesUsed].ptr = key;
      }
      // save the probability
      cliqueValues.ptr[numCliqueValuesUsed].p = final_p;
      numCliqueValuesUsed++;

      /*
       * uncomment this next code to produce lots of messages.
       * note above, if high verbosity is on, recursive version of this
       * routine will print this message.
       */
      /*
      if (message(Mega)) {
	// psp(stdout,2*nodeNumber);
	infoMsg(Mega,"Inserting New Clique Value. prob = %f, sum = %f: ",
		cliqueValues.ptr[numCliqueValuesUsed-1].p.val(),sumProbabilities().val());
	printRVSetAndValues(stdout,fNodes);
      }
      */

    }

  next_iteration:
    // Now we need to move to next clique value. This bunch of next
    // code is like an inlined iterator over clique values. It is done
    // as "inline" in an attempt to avoid branch mis-predicts and since
    // this code appears only one time.
    do {

      bool unfinished;
      switch (cur_disp) {
      case MaxClique::AN_CPT_ITERATION_COMPUTE_AND_APPLY_PROB: 
	{
	  unfinished = cur_rv->next(cur_p);
	  // assert ( !unfinished || cur_p.not_essentially_zero() );
	  parray[nodeNumber] = parray[nodeNumber-1]*cur_p;
	}
	break;

      case MaxClique::AN_CPT_ITERATION_COMPUTE_PROB_REMOVE_ZEROS:
	{
	  // guaranteed not to get zero prob here. Note that
	  // we get the probability into cur_p, but we do not
	  // use it in this case.
	  unfinished = cur_rv->next(cur_p);
	  // we could include:
	  // assert ( !cur_p.essentially_zero() );
	}
	break;

      case MaxClique::AN_CARD_ITERATION:
	{
	  DiscRV* drv = (DiscRV*)cur_rv;
	  unfinished = (++drv->val < drv->cardinality);
	}
	break;

      default:
	// these cases are handeled by 'default':
	// case MaxClique::AN_COMPUTE_AND_APPLY_PROB:
	// case MaxClique::AN_CONTINUE:
	// case MaxClique::AN_CONTINUE_COMPUTE_PROB_REMOVE_ZEROS:
	{
	  unfinished = false;
	}
	break;
      }
      
      if (unfinished) {
	// Best case, we continue on to next outer iteration filling
	// in next clique value.
	break;
      } else {
	// we're finished with the current variable.
	if (nodeNumber == 0) {
	  // then we're really done.
	  goto done_with_clique_iteration;
	} else {
	  nodeNumber--;
	  cur_disp = origin.dispositionSortedAssignedNodes.ptr[nodeNumber];
	  cur_rv = fSortedAssignedNodes.ptr[nodeNumber];
	}
      }

      // still here? Continue on incrementing previous variable.
    } while (1);

    while (nodeNumber < maxNodeNumber) {
      nodeNumber++;
      cur_disp = origin.dispositionSortedAssignedNodes.ptr[nodeNumber];
      cur_rv = fSortedAssignedNodes.ptr[nodeNumber];


      // need to fill up the rest of the table.
      switch (cur_disp) {
      case MaxClique::AN_CPT_ITERATION_COMPUTE_AND_APPLY_PROB: 
	{
	  cur_rv->begin(cur_p);
	  // might get a zero probability, check condition here.
	  if (cur_p.essentially_zero()) {
	    // Since we have zero here, we cancel iterations of all
	    // subsequent variables right now, rather than iterate
	    // them with what will end up being zero probability.
	    nodeNumber--;
	    cur_disp = origin.dispositionSortedAssignedNodes.ptr[nodeNumber];
	    cur_rv = fSortedAssignedNodes.ptr[nodeNumber];
	    goto next_iteration;
	  } else 
	    parray[nodeNumber] = parray[nodeNumber-1]*cur_p;
	}
	break;

      case MaxClique::AN_CPT_ITERATION_COMPUTE_PROB_REMOVE_ZEROS:
	{
	  cur_rv->begin(cur_p);
	  // might get a zero probability, check condition here.
	  if (cur_p.essentially_zero()) {
	    // Since we have zero here, we cancel iterations of all
	    // subsequent variables right now, rather than iterate them
	    // with what will end up being zero probability.
	    nodeNumber--;
	    cur_disp = origin.dispositionSortedAssignedNodes.ptr[nodeNumber];
	    cur_rv = fSortedAssignedNodes.ptr[nodeNumber];
	    goto next_iteration;
	  } else {
	    // TODO: update post-condition comments in RV and CPT iterators as
	    // to when a zero RV can occur and when not.
	    parray[nodeNumber] = parray[nodeNumber-1];
	  }
	}
	break;

      case MaxClique::AN_CARD_ITERATION:
	{
	  DiscRV* drv = (DiscRV*)cur_rv;
	  drv->val = 0;
	  parray[nodeNumber] = parray[nodeNumber-1];
	}
	break;

      case MaxClique::AN_COMPUTE_AND_APPLY_PROB:
	{
	  // TODO: Make more efficient version of this, based on the type of
	  // RV.
	  cur_rv->probGivenParents(cur_p);
	  // might get a zero probability, check condition here.
	  if (cur_p.essentially_zero()) {
	    // Since we have zero here, we cancel iterations of all
	    // subsequent variables right now, rather than iterate them
	    // with what will end up being zero probability.
	    nodeNumber--;
	    cur_disp = origin.dispositionSortedAssignedNodes.ptr[nodeNumber];
	    cur_rv = fSortedAssignedNodes.ptr[nodeNumber];
	    goto next_iteration;
	  } else 
	    parray[nodeNumber] = parray[nodeNumber-1]*cur_p;
	}
	break;

      case MaxClique::AN_CONTINUE:
	parray[nodeNumber] = parray[nodeNumber-1];
	break;

      case MaxClique::AN_CONTINUE_COMPUTE_PROB_REMOVE_ZEROS:
	{
	  // TODO: Make more efficient version of this, based on the type of
	  // RV.
	  cur_rv->probGivenParents(cur_p);
	  // might get a zero probability check condition here.
	  if (cur_p.essentially_zero()) {
	    // Since we have zero here, we cancel iterations of all
	    // subsequent variables right now, rather than iterate them
	    // with what will end up being zero probability.
	    nodeNumber--;
	    cur_disp = origin.dispositionSortedAssignedNodes.ptr[nodeNumber];
	    cur_rv = fSortedAssignedNodes.ptr[nodeNumber];
	    goto next_iteration;
	  } else 
	    parray[nodeNumber] = parray[nodeNumber-1];
	}
	break;
      }
    }

  } while (1);

 done_with_clique_iteration:
  ;
}


/*-
 *-----------------------------------------------------------------------
 * MaxClique::ceIterateUnassignedIteratedNodes()
 *
 *    Collect Evidence, Iterate Separators: This routine
 *    is part of separator driven clique instantiation. Once
 *    we have a partial clique value that has survied the separators,
 *    the next thing we need to do is iterate through any unassigned
 *    nodes in this cliuqe (if any, hopefully not since these can be costly
 *    and might indicate a poor triangulation, depending on pruning and/or
 *    the sparse/deterministic variables).
 *
 *    If we have reached the last unassigned node in this clique, we
 *    move directly on to the assigned nodes.
 *
 * Preconditions:
 *
 *   All the cliques incomming separators *must* have been created and are ready
 *   to be used to produce the clique table.
 *     
 *
 * Postconditions:
 *    We move on to the asssigned nodes with an probabilty. The prob. has
 *    not been updated from what is passed in since these nodes are unassigned
 *    and so never contribute score to this clique.
 *
 * Side Effects:
 *    none
 *
 * Results:
 *     nothing
 *
 *-----------------------------------------------------------------------
 */

void
InferenceMaxClique::ceIterateUnassignedIteratedNodes(JT_InferencePartition& part,
						     const unsigned nodeNumber,
						     const logpr p)
{
  if (nodeNumber == fUnassignedIteratedNodes.size()) {
    ceIterateAssignedNodes(part,0,p);
    return;
  }
  RV* rv = fUnassignedIteratedNodes[nodeNumber];
  infoMsg(Giga,"Starting Unassigned iteration of rv %s(%d), nodeNumber = %d, p = %f\n",
	  rv->name().c_str(),rv->frame(),nodeNumber,p.val());

  if (rv->hidden()) {
    // only discrete RVs can be hidden for now.
    DiscRV* drv = (DiscRV*)rv;
    // do the loop right here
    drv->val = 0;
    do {
      infoMsg(Giga,"Unassigned iteration of rv %s(%d)=%d, nodeNumber = %d, p = %f\n",
	      rv->name().c_str(),rv->frame(),drv->val,nodeNumber,p.val());
      // continue on, effectively multiplying p by unity.
      ceIterateUnassignedIteratedNodes(part,nodeNumber+1,p);
    } while (++drv->val < drv->cardinality);
  } else {
    // TODO: Perhaps unassignedIteratedNodes should contain
    // no observed nodes at all since we are not updating 
    // probability here anyway.

    // observed, either discrete or continuous
    if (rv->discrete()) {
      // DiscRV* drv = (DiscRV*)rv;
      // TODO: for observed variables, do this once at the begining
      // before any looping here.
      // drv->setToObservedValue();
      infoMsg(Giga,"Unassigned pass through of observed rv %s(%d)=%d, nodeNumber = %d, p = %f\n",
	      rv->name().c_str(),rv->frame(),RV2DRV(rv)->val,nodeNumber,p.val());
    } else {
      // nothing to do since we get continuous observed
      // value indirectly
      infoMsg(Giga,"Unassigned pass through of observed rv %s(%d)=C, nodeNumber = %d, p = %f\n",
	      rv->name().c_str(),rv->frame(),nodeNumber,p.val());
    }
    // continue on, effectively multiplying p by unity.
    ceIterateUnassignedIteratedNodes(part,nodeNumber+1,p);
  }
}


/*-
 *-----------------------------------------------------------------------
 * MaxClique::ceSendToOutgoingSeparator()
 *
 *    Collect Evidence, Send to outgoing separator.
 *
 *    We have now a fully instantiated clique. This routine Iterates
 *    through the values that are above beam in the clique table, and
 *    instantiates the outgoing separator with those values. If a separator
 *    value has multiple clique values, then values are accumulated into
 *    the separator (in a Veterbi approach, we would take the max, see
 *    veterbi code). 
 *
 *    This code also combine clique pruning right here, rather than needing
 *    to do a separate pruning stage by calling ceCliquePrune(). It also
 * 
 *
 * Preconditions:
 *
 *   The clique table must be fully instantiated, and contain entries
 *   for this clique.  The outgoing separator must have been created,
 *   setup, and ready to accept a message (projection down) from a
 *   clique
 *
 * Postconditions:
 *    The outgoing separator has been instantiated.
 *
 * Side Effects:
 *    Changes 
 *      1) the clique table, since with pruning, it will potentially shrink the clique table
 *      2) the outgoing separator table, since it is being created.
 *
 * Results:
 *     nothing
 *
 *-----------------------------------------------------------------------
 */
void 
InferenceMaxClique::
ceSendToOutgoingSeparator(JT_InferencePartition& part)
{
  ceSendToOutgoingSeparator(part,
		       part.separatorCliques[origin.ceSendSeparator]);
}
void 
InferenceMaxClique::
ceSendToOutgoingSeparator(JT_InferencePartition& part,
			  InferenceSeparatorClique& sep)
{

  // first check if this is an all observed clique.
  if (origin.hiddenNodes.size() == 0) {
    // Everything in this clique is observed.  Therefore, there should
    // be one and only one clique value in this clique which is the
    // probability of the assigned probability nodes in this clique.

    // Since this clique is a superset of the separator, it means also
    // that the seperator will consist only of observed nodes. We
    // therefore create one accumulated entry and one rem entry and
    // insert our current probability.

    // first do some sanity checks.
    assert (sep.origin.hAccumulatedIntersection.size() == 0);
    assert (sep.origin.hRemainder.size() == 0);
    assert (sep.remDiscreteValuePtrs.size() == 0);
    
    const unsigned accIndex = 0;

    // We here keep handy reference for readability.
    // To find out where has this been allocated, 
    // search for tag: ALLOCATE_REMVALUES_ALL_OBSERVED.
    InferenceSeparatorClique::AISeparatorValue& sv
      = sep.separatorValues.ptr[accIndex];

    // This must be first time for this entry.  We never need more
    // than one rem value.  Search for tag
    // 'ALLOCATE_REMVALUES_ALL_OBSERVED' in this file for where else
    // this could be done, but we allocate it here.
    if (sv.remValues.size() == 0)
      sv.remValues.resize(1); 
    // in all cases
    sv.numRemValuesUsed = 1;	  

    sv.remValues.ptr[0].p = cliqueValues.ptr[0].p;

    // and we're done already. This was easy!
    return;
  } 

  // do k-pruning (ideally we would do this after
  // beam pruning, but since beam pruning is integrated
  // into the code below, we do max state pruning first).
  // Prune the minimum of the fixed K size and the percentage size.
  unsigned k;
  k = 2 + (unsigned)((origin.cliqueBeamRetainFraction)*(double)numCliqueValuesUsed);
  if (origin.cliqueBeamMaxNumStates > 0) {
    k = min(k,origin.cliqueBeamMaxNumStates);
  }
  //   printf("nms = %d, pf = %f, ncv = %d, k = %d\n",origin.cliqueBeamMaxNumStates,
  // origin.cliqueBeamRetainFraction,numCliqueValuesUsed,k);
  ceCliquePrune(k);

  // create an ininitialized variable using dummy argument
  logpr beamThreshold((void*)0);
  if (origin.cliqueBeam != (-LZERO)) {
    // then we do clique table pruning right here rather
    // than a separate call to ceCliquePrune().
    // break into the logp to avoid unnecessary zero checking.
    beamThreshold.valref() = maxCEValue.valref() - origin.cliqueBeam;
  } else {
    // set beam threshold to a value that will never cause pruning.
    beamThreshold.set_to_zero();
  }


  const unsigned origNumCliqueValuesUsed = numCliqueValuesUsed;  

  // next check if the outgoing separator has only obseved values.
  if (sep.origin.hAccumulatedIntersection.size() == 0 && sep.origin.hRemainder.size() == 0) {

    // Then indeed, outgoing separator is all observed values. We do
    // this special separately from the general case since that case
    // is already getting a bit unwieldy.

    InferenceSeparatorClique::AISeparatorValue& sv
      = sep.separatorValues.ptr[0];

    // This must be first time for this entry.  We never need more
    // than one rem value.  Search for tag
    // 'ALLOCATE_REMVALUES_ALL_OBSERVED' in this file for where else
    // this could be done, but we allocate it here.
    if (sv.remValues.size() == 0)
      sv.remValues.resize(1); 

    // Just sum up the clique entries projecting down into the
    // observed separator, but only do the ones that pass the beam
    // threshold. Arguably, we might not want or need to do beam
    // pruning here since when the separator is all observed, any
    // sparseness that we introduce by clique pruning isn't going to
    // change things later on (since everything is projected to the
    // same point), but we do it here anyway for numerical consistency
    // with the general case.
    for (unsigned cvn=0;cvn<numCliqueValuesUsed;) {
      if (cliqueValues.ptr[cvn].p < beamThreshold) {
	// swap with last entry, and decrease numCliqueValuesUsed by
	// one.
	swap(cliqueValues.ptr[cvn],cliqueValues.ptr[numCliqueValuesUsed-1]);
	numCliqueValuesUsed--;
	// continue on to next iteration without incrementing cvn
	continue;
      } else {
	if (JunctionTree::viterbiScore) {
	  // sv.remValues.ptr[0].p.assign_if_greater(cliqueValues.ptr[cvn].p);
	  if (cliqueValues.ptr[cvn].p > sv.remValues.ptr[0].p) {
	    sv.remValues.ptr[0].p = cliqueValues.ptr[cvn].p;
	    sv.remValues.ptr[0].backPointer = cvn;
	  }
	} else {
	  sv.remValues.ptr[0].p += cliqueValues.ptr[cvn].p;
	}
	cvn++;
      }
    }

  } else {
    // We are guaranteed that either we have an accumulated
    // intersection, a reminder, or both (but not neither).  Go
    // through clique values and accumulate into appropriate place
    // within separator.

    // pre-compute a few values that the compiler might spend more
    // time than necessary to re-check.
    const bool imc_nwwoh_p = (origin.packer.packedLen() <= IMC_NWWOH);
    const bool isc_nwwoh_ai_p = (sep.origin.accPacker.packedLen() <= ISC_NWWOH_AI);
    const bool isc_nwwoh_rm_p = (sep.origin.remPacker.packedLen() <= ISC_NWWOH_RM);
    const bool sep_origin_hAccumulatedIntersection_exists_p =
      (sep.origin.hAccumulatedIntersection.size() > 0);
    const bool sep_remDiscreteValuePtrs_exists_p = 
      (sep.remDiscreteValuePtrs.size() > 0);

    for (unsigned cvn=0;cvn<numCliqueValuesUsed;) {{

      if (cliqueValues.ptr[cvn].p < beamThreshold) {
	// swap with last entry, and decrease numCliqueValuesUsed by
	// one.
	swap(cliqueValues.ptr[cvn],cliqueValues.ptr[numCliqueValuesUsed-1]);
	numCliqueValuesUsed--;
	// continue on to next iteration without incrementing cvn
	continue;
      }

      // TODO: optimize away this conditional check. (and/or use const
      // local variable to indicate it wont change)
      if (imc_nwwoh_p) {
	origin.packer.unpack((unsigned*)&(cliqueValues.ptr[cvn].val[0]),
			     (unsigned**)discreteValuePtrs.ptr);
      } else {
	origin.packer.unpack((unsigned*)cliqueValues.ptr[cvn].ptr,
			     (unsigned**)discreteValuePtrs.ptr);
      }

      // All hidden random variables now have their discrete
      // value. Accumulate this probability into the given separator.

      /*
       * There are 3 cases.
       * 1) AI exists and REM exist
       * 2) AI exists and REM doesnt exist
       * 3) AI does not exist, but REM exists
       * AI not exist and REM not exist can't occur.
       */

      unsigned accIndex;
      // TODO: optimize this check away out of loop.
      if (sep_origin_hAccumulatedIntersection_exists_p) { 
	// an accumulated intersection exists.

	// make sure there is at least one available accumulated intersection entry
	assert ( sep.numSeparatorValuesUsed <= sep.separatorValues.size());
	if (sep.numSeparatorValuesUsed >= sep.separatorValues.size()) {
	  const unsigned old_size = sep.separatorValues.size();
	  // TODO: optimize this size re-allocation.
	  // sep.separatorValues.resizeAndCopy(sep.separatorValues.size()*2);
	  if (sep.numSeparatorValuesUsed >= sep.origin.separatorValueSpaceManager.currentSize()) 
	    sep.origin.separatorValueSpaceManager.advanceToNextSize();
	  sep.separatorValues.resizeAndCopy(sep.origin.separatorValueSpaceManager.currentSize()); 
	  if (isc_nwwoh_ai_p) {
	    // Then the above resize just invalided all our pointers to keys,
	    // but it did not invalidate the array indices. Go through
	    // and correct the keys within the hash table.
	    // TODO: think of a better way to do this that also looses no efficiency.
	    for (unsigned i=0;i<sep.iAccHashMap.tableSize();i++) {
	      if (!sep.iAccHashMap.tableEmpty(i)) {
		sep.iAccHashMap.tableKey(i)
		  = &(sep.separatorValues.ptr[sep.iAccHashMap.tableItem(i)].val[0]);
	      }
	    }
	  }
	  const unsigned new_size = sep.separatorValues.size();
	  // if (sep.remDiscreteValuePtrs.size() > 0) {
	  if (sep_remDiscreteValuePtrs_exists_p) {
	    for (unsigned i=old_size;i<new_size;i++) {
	      // re-construct hash tables only for new entries.
	      new (&sep.separatorValues.ptr[i].iRemHashMap)
		VHashMapUnsignedUnsignedKeyUpdatable
		(sep.origin.remPacker.packedLen(),2);
	      // TODO: potentially preallocate default size of  separatorValues.ptr[i].remValues.resize(default);
	    }
	  }
	}
      
	unsigned *accKey;
	// TODO: optimize this check out of loop.
	if (isc_nwwoh_ai_p) {
	  accKey = &(sep.separatorValues.ptr[sep.numSeparatorValuesUsed].val[0]);
	  sep.origin.accPacker.pack((unsigned**)sep.accDiscreteValuePtrs.ptr,
				    accKey);
	} else {
	  accKey = sep.origin.accValueHolder.curCliqueValuePtr();
	  sep.origin.accPacker.pack((unsigned**)sep.accDiscreteValuePtrs.ptr,
				    accKey);
	  // check if this value combination already lives in
	  // origin's value holder hash table and if so, use that.
	  bool foundp;
	  accKey = sep.origin.accSepValHashSet.insert(accKey,foundp);
	  if (!foundp) {
	    // only allocate a new value if it was inserted.
	    sep.origin.accValueHolder.allocateCurCliqueValue();
	  }
	  // store the pointer in case we use it.
	  sep.separatorValues.ptr[sep.numSeparatorValuesUsed].ptr = accKey;
	}
      
	bool foundp;
	unsigned* accIndexp =
	  sep.iAccHashMap.insert(accKey,
				 sep.numSeparatorValuesUsed,
				 foundp);

	if (!foundp) {
	  //  add the values we just used. 
	  sep.numSeparatorValuesUsed++;
	}
	accIndex = *accIndexp;

	// TODO: optimize this check out of loop.
	// if (sep.remDiscreteValuePtrs.size() == 0) {
	if (!sep_remDiscreteValuePtrs_exists_p) {
	  // 2) AI exists and REM doesnt exist
	  // Then this separator is entirely covered by one or 
	  // more other separators earlier in the order.

	  // go ahead and insert it here to the 1st entry (entry 0).

	  // handy reference for readability.
	  InferenceSeparatorClique::AISeparatorValue& sv
	    = sep.separatorValues.ptr[accIndex];

	  // Accumulate the clique's
	  // probability into this separator's probability.
	  if (sv.remValues.size() < 1) {
	    // This must be first time for this entry.
	    // Search for tag 'ALLOCATE_REMVALUES_OPTION' in this file for where else this could
	    // be done.
	    sv.remValues.resize(1);
	    sv.numRemValuesUsed = 1;	  
	    // initialize and assign.
	    sv.remValues.ptr[0].p = cliqueValues.ptr[cvn].p;
	    if (JunctionTree::viterbiScore)
	      sv.remValues.ptr[0].backPointer = cvn;
	  } else {
	    // already there so must have hit before.
	    // we thus accumulate.
	    if (JunctionTree::viterbiScore) {
	      // sv.remValues.ptr[0].p.assign_if_greater(cliqueValues.ptr[cvn].p);
	      if (cliqueValues.ptr[cvn].p > sv.remValues.ptr[0].p) {
		sv.remValues.ptr[0].p = cliqueValues.ptr[cvn].p;
		sv.remValues.ptr[0].backPointer = cvn;
	      }
	    } else {
	      sv.remValues.ptr[0].p += cliqueValues.ptr[cvn].p;
	    }
	  }

	  goto next_iteration;
	}

      } else {
	accIndex = 0;
      }

      // If we're here, then we are guaranteed must have some remainder
      // pointers, i.e., we could do:
      //    assert (sep.remDiscreteValuePtrs.size() > 0);
      // So, the remainder exists in this separator.
      // either:
      //   1) AI exists and REM exist
      //     or
      //   3) AI does not exist (accIndex == 0), but REM exists
      // 

      // keep handy reference for readability.
      InferenceSeparatorClique::AISeparatorValue& sv
	= sep.separatorValues.ptr[accIndex];
    
      // make sure there is at least one available entry
      assert (sv.numRemValuesUsed <= sv.remValues.size());
      if (sv.numRemValuesUsed >= sv.remValues.size()) {
	// TODO: optimize this growth rate.
	// start small but grow fast.
	// sv.remValues.resizeAndCopy(1+sv.remValues.size()*2); // *3
	sv.remValues.resizeAndCopy(sep.origin.remainderValueSpaceManager.nextSizeFrom(sv.remValues.size()));
	sep.origin.remainderValueSpaceManager.setCurrentAllocationSizeIfLarger(sv.remValues.size());
	if (isc_nwwoh_rm_p) {
	  // Then the above resize just invalided all sv.iRemHashMap's pointers to keys,
	  // but it did not invalidate its array indices. Go through
	  // and correct the keys within the hash table.
	  // TODO: think of a better way to do this that looses no efficiency.
	  for (unsigned i=0;i<sv.iRemHashMap.tableSize();i++) {
	    if (!sv.iRemHashMap.tableEmpty(i)) {
	      sv.iRemHashMap.tableKey(i)
		= &(sv.remValues.ptr[sv.iRemHashMap.tableItem(i)].val[0]);
	    }
	  }
	}
      }

      unsigned *remKey;
      // pack relevant variable values
      // TODO: optimize away this check.
      if (isc_nwwoh_rm_p) {
	// grab pointer to next location to be used in this case.
	remKey = &(sv.remValues.ptr[sv.numRemValuesUsed].val[0]);
	// pack the remainder pointers
	sep.origin.remPacker.pack((unsigned**)sep.remDiscreteValuePtrs.ptr,
				  remKey);
      } else {
	// grab pointer to next packed clique value to be used.
	remKey = sep.origin.remValueHolder.curCliqueValuePtr();
	sep.origin.remPacker.pack((unsigned**)sep.remDiscreteValuePtrs.ptr,
				  remKey);
	// check if this value combination already lives in
	// origin's value holder hash table and if so, use that.
	bool foundp;
	remKey = sep.origin.remSepValHashSet.insert(remKey,foundp);
	if (!foundp) {
	  // only allocate a new value if it was inserted.
	  sep.origin.remValueHolder.allocateCurCliqueValue();
	}
	// store the pointer in case we use it.
	sv.remValues.ptr[sv.numRemValuesUsed].ptr = remKey;
      }

      bool foundp;
      unsigned* remIndexp =
	sv.iRemHashMap.insert(remKey,
			      sv.numRemValuesUsed,
			      foundp);
      if (!foundp) {
	// add the values we just used. 
	sv.numRemValuesUsed++;
      }

      // We've finally got the entry, so accumulate the clique's
      // probability into this separator's probability.
      if (JunctionTree::viterbiScore) {
	// sv.remValues.ptr[*remIndexp].p.assign_if_greater(cliqueValues.ptr[cvn].p);
	if (cliqueValues.ptr[cvn].p > sv.remValues.ptr[*remIndexp].p) {
	  sv.remValues.ptr[*remIndexp].p = cliqueValues.ptr[cvn].p;
	  sv.remValues.ptr[*remIndexp].backPointer = cvn;
	}
      } else {
	sv.remValues.ptr[*remIndexp].p += cliqueValues.ptr[cvn].p;
      }
    }
    next_iteration:
    cvn++;
    }
  }
  
  // To reallocate or not to reallocate, that is the question.  here,
  // we just reallocate for now.
  // 
  // TODO: reallocate only if change is > some percentage (say 5%),
  // and export to command line.
  // TODO: possibly tell origin.cliqueValueSpaceManager about this pruning event.
  // 
  if (numCliqueValuesUsed < origNumCliqueValuesUsed)
    cliqueValues.resizeAndCopy(numCliqueValuesUsed);

  // Keep this out of previous check so we can also use large beams to
  // take a look at what state space size is with large verbosity
  // value.
  if (origin.cliqueBeam != (-LZERO)) {
    infoMsg(IM::Huge,"Clique beam pruning, Max cv = %f, thres = %f. Original clique state space = %d, new clique state space = %d\n",
	    maxCEValue.valref(),
	    beamThreshold.valref(),
	    origNumCliqueValuesUsed,
	    numCliqueValuesUsed);
  }


  /////////////////////////////////////
  // And prune the separator as well.
  sep.ceSeparatorPrune();

}


/*-
 *-----------------------------------------------------------------------
 * InferenceMaxClique::ceCliquePrune()
 *
 *    Collect Evidence, Clique Prune: This routine will prune away
 *    part of a previously instantiated clique based on the current
 *    clique beam width.
 *    
 *    Note that InferenceMaxClique::ceSendToOutgoingSeparator() does
 *    its own pruning, so when using ceSendToOutgoingSeparator(), this
 *    pruning routine does not need to be called (at least with the
 *    same beam width). This routine is kept here for pedagogical
 *    purposes however.
 *
 * Preconditions:
 *   1) the value of the max clique 'maxCEValue' must have been
 *      computed already.
 *
 *   2) clique table must be created, meaning that either:
 *
 *        InferenceMaxClique::ceIterateAssignedNodesRecurse()
 *   or
 *        InferenceMaxClique::ceIterateAssignedNodesCliqueDriven
 *   must have just been called.
 *
 * Postconditions:
 *    Clique table has been pruned, and memory for it has been re-allocated to
 *    fit the smaller size.
 *
 * Side Effects:
 *    changes the clique size
 *
 * Results:
 *     nothing
 *
 *-----------------------------------------------------------------------
 */

void 
InferenceMaxClique::ceCliquePrune()
{
  // return immediately if beam pruning is turned off.
  if (origin.cliqueBeam == (-LZERO))
    return;
  
  // create an ininitialized variable
  logpr beamThreshold((void*)0);
  // break into the logp to avoid unnecessary zero checking.
  beamThreshold.valref() = maxCEValue.valref() - origin.cliqueBeam;

  const unsigned origNumCliqueValuesUsed = numCliqueValuesUsed;
  for (unsigned cvn=0;cvn<numCliqueValuesUsed;) {
    if (cliqueValues.ptr[cvn].p < beamThreshold) {
      // swap with last entry, and decrease numCliqueValuesUsed by one.
      swap(cliqueValues.ptr[cvn],cliqueValues.ptr[numCliqueValuesUsed-1]);
      numCliqueValuesUsed--;
    } else {
      cvn++;
    }
  }

  infoMsg(IM::Huge,"Clique beam pruning: Max cv = %f, thres = %f. Original clique state space = %d, new clique state space = %d\n",
	  maxCEValue.valref(),
	  beamThreshold.valref(),
	  origNumCliqueValuesUsed,
	  numCliqueValuesUsed);
  
  // To reallocate or not to reallocate, that is the question.  here,
  // we just reallocate for now.
  // 
  // TODO: reallocate only if change is > some percentage (say 5%),
  // and export to command line.
  // e.g., if ((origNumCliqueValuesUsed - numCliqueValuesUsed) > 0.05*origNumCliqueValuesUsed)
  // 
  if (numCliqueValuesUsed < origNumCliqueValuesUsed)
    cliqueValues.resizeAndCopy(numCliqueValuesUsed);

}



/*-
 *-----------------------------------------------------------------------
 * InferenceMaxClique::ceCliquePrune(unsigned k)
 *
 *    Collect Evidence, Clique Prune: This routine will prune away
 *    part of a previously instantiated clique so that it has only
 *    k entries left. It does this by chosing the k-best (maximum partial
 *    probability) clique entries.
 *
 * Preconditions:
 *   1) the value of the max clique 'maxCEValue' must have been
 *      computed already.
 *
 *   2) clique table must be created, meaning that either:
 *
 *        InferenceMaxClique::ceIterateAssignedNodesRecurse()
 *   or
 *        InferenceMaxClique::ceIterateAssignedNodesCliqueDriven
 *   must have just been called.
 *
 * Postconditions:
 *    Clique table has been pruned, and memory for it has been re-allocated to
 *    fit the smaller size.
 *
 * Side Effects:
 *    changes the clique size.
 *
 * Results:
 *     nothing
 *
 *-----------------------------------------------------------------------
 */

void 
InferenceMaxClique::ceCliquePrune(const unsigned k)
{
  // k can't be larger than the number of clique entries.
  if (k == 0 || k >= numCliqueValuesUsed)
    return;

  // inclusive range to process.
  unsigned lower = 0;
  unsigned upper = numCliqueValuesUsed-1;

  // We find the k max elements using a quick sort-like algorithm, and
  // start by processing the elements between lower and upper
  // inclusive. Note that unlike quicksort (which is O(NlogN) average
  // case and O(N^2) worse case), this algorithm is O(N) average case
  // (yes, not dep. on k) since we're only finding the top k values,
  // and those top k values need not be sorted.

  do {
    // k'th largest lives within [lower,upper] inclusive.
    // printf("lower = %d, upper = %d\n",lower,upper);

    if (lower == upper)
      break;

    // Choose a random pivot and then divide the array so that
    // elements >= pivot are on the left, and elements < pivot are on
    // the right. We choose a pivot by taking the median of 3 randomly
    // choosen values.
    unsigned pl1 = rnd.uniform(lower,upper);
    unsigned pl2 = rnd.uniform(lower,upper);
    unsigned pl3 = rnd.uniform(lower,upper);
    unsigned pl;
    // find rough median
    if (cliqueValues[pl1].p < cliqueValues[pl2].p) {
      if (cliqueValues[pl2].p < cliqueValues[pl3].p) {
	pl = pl2;
      } else {
	pl = pl3;
      }
    } else {
      if (cliqueValues[pl1].p < cliqueValues[pl3].p) {
	pl = pl1;
      } else {
	pl = pl3;
      }
    }

    logpr pivot = cliqueValues[pl].p;
    // printf("pivot location = %d, pivot = %f\n",pl,pivot);

    // swap pivot into first position, which is a valid position for
    // pivot, since pivot >= pivot.
    swap(cliqueValues[pl],cliqueValues[lower]);

    unsigned l = lower+1;
    unsigned u = upper;
    while (l < u) {
      if (cliqueValues[l].p >= pivot) {
	l++;
      } else {
	swap(cliqueValues[l],cliqueValues[u]);
	u--;
      }
    }

    // 1) Now all entries with index < l are >= pivot,
    //   which means at this point we know we have
    //   the top l entries at indices [0,l]
    // 2) All entries with index > u are < pivot
    // 3) we have l == u
    // 4) The entry at index u==l is unknown however.


    if (cliqueValues[l].p >= pivot) {
      l++;
      u++;
    }
    // 1) Now all entries with index < l are >= pivot,
    //   which means at this point we know we have
    //   the top l entries at indices [0,l]
    // 2) All entries with index >= u are < pivot
    // 3) we still have l == u

    // put pivot in its appropriate place.
    l--;
    swap(cliqueValues[lower],cliqueValues[l]);
    // 1) Now all entries with index <= l are >= pivot,
    //   which means at this point we know we have
    //   the top l entries at indices [0,l]
    // 2) All entries with index >= u are < pivot
    // 3) we have l +1 == u

    // Now, deal with potential ties.
    // Move any other entries that are equal to pivot to the center.
    // -- Note that this next step helps to significantly speed up the
    // -- algorithm in the case when there are lots of ties --- this
    // -- is done here since one of the main reasons for having this
    // -- form of pruning is that normal beam pruning doesn't work
    // -- well on very "uniform"-like unnormalized distributions,
    // -- i.e., ones for which many ties might be present.
    unsigned ll;
    if (l > lower) {
      ll = l-1; // left of left-most known pivot value.
      unsigned i = lower;
      while (i < ll) {
	if (cliqueValues[i].p == pivot) {
	  swap(cliqueValues[i],cliqueValues[ll]);
	  ll--;
	} else {
	  i++;
	}
      }
    } else
      ll = l;

    // 1) Now all entries with index <= l are >= pivot,
    //   which means at this point we know we have
    //   the top l entries at indices [0,l]
    // 2) All entries with index >= u are < pivot
    // 3) we have l +1 == u
    // 4) either
    //    a) all entries with index <= ll are > pivot, or
    //         (meaning that ll+1 is the number of entries
    //          that are strictly greater than pivot)
    //    b) ll = lower and arr[ll] == pivot
    //       (in which case, ll is the number of entries
    //          that are strictly greater than pivot).

    //     printf("after swapping, l = %d, u = %d,array now:\n",l,u);
    //     for (unsigned j=0;j<len;j++)
    //       printf("arr[%d] = %f\n",j,arr[j]);

    if (ll > k) {
      upper = ll;
    } else if (l > k) {
      // then we're done since we must
      // have a bunch of ties, and pivot value
      // must be correct k'th value.
      break;
    } else if (l == k) {
      // then we have l+1 entries
      // printf("finished, we have l == k == %d\n",l);
      break;
    } else { // l < k, need to search to the right.
      // printf("setting lower to %d\n",u);
      lower = u;
    }

  } while(1);

  infoMsg(IM::Huge,"Clique k-beam pruning: Original clique state space = %d, new clique state space = %d\n",
	  numCliqueValuesUsed,k);
  
  // To reallocate or not to reallocate, that is the question.  here,
  // we just reallocate for now.
  // 
  // TODO: reallocate only if change is > some percentage (say 5%),
  // and export to command line.
  // e.g., if ((origNumCliqueValuesUsed - numCliqueValuesUsed) > 0.05*origNumCliqueValuesUsed)
  // 
  cliqueValues.resizeAndCopy(k);
  numCliqueValuesUsed = k;

}




//////////////
// Clique driven version of gather from incomming separators
// TODO: add comments
/////////////

void
InferenceMaxClique::ceGatherFromIncommingSeparatorsCliqueDriven(JT_InferencePartition& part)
{
  assert (MaxClique::ceSeparatorDrivenInference == false); 

  if (origin.hiddenNodes.size() == 0) {
    ceGatherFromIncommingSeparatorsCliqueObserved(part);
  } else {
    logpr p = 1.0;
    if (origin.unassignedNodes.size() == 0) {
      ceIterateAssignedNodesCliqueDriven(part,0,p);
    } else {
      ceIterateUnassignedNodesCliqueDriven(part,0,p);
    }
  }

}


//
// a version that is specifically for cliques that are all observed which
// means that all surrounding separators are also all observed.
// observed clique observedclique
void
InferenceMaxClique::ceGatherFromIncommingSeparatorsCliqueObserved(JT_InferencePartition& part)
{

  // just do the assigned nodes as the unassigned observed nodes have no effect here.
  unsigned nodeNumber;

  logpr p = 1.0;
  bool zero = false;
  for (nodeNumber = 0; nodeNumber < fSortedAssignedNodes.size(); nodeNumber ++ ) {
    RV* rv = fSortedAssignedNodes[nodeNumber];
    if (origin.dispositionSortedAssignedNodes[nodeNumber] == MaxClique::AN_NOTSEP_PROB_SPARSEDENSE || 
	origin.dispositionSortedAssignedNodes[nodeNumber] == MaxClique::AN_NOTSEP_PROB_SPARSEDENSE) {
      // apply the probabiltiy
      logpr cur_p = rv->probGivenParents();

      if (message(Giga)) {
	if (!rv->discrete()) {
	  infoMsg(Giga,"Observed Clique prob application of rv %s(%d)=C, nodeNumber =%d, cur_p = %f, p = %f\n",
		  rv->name().c_str(),rv->frame(),nodeNumber,cur_p.val(),p.val());
	} else {
	  DiscRV* drv = (DiscRV*)rv;
	  infoMsg(Giga,"Observed Clique prob application of rv %s(%d)=%d, nodeNumber =%d, cur_p = %f, p = %f\n",
		  drv->name().c_str(),drv->frame(),drv->val,nodeNumber,cur_p.val(),p.val());
	}
      }

      // if at any step, we get zero, then back out.
      if (cur_p.essentially_zero()) {
	// we've got a zero, so might as well stop here now. 
	zero = true;
	break;
      } else {
	p *= cur_p;
      }
    } else { 
      // in none of the other cases do we apply the probability. We
      // still check for zeros though.
      logpr cur_p = rv->probGivenParents();

      if (message(Giga)) {
	if (!rv->discrete()) {
	  infoMsg(Giga,"Observed Clique prob non-application of rv %s(%d)=C, nodeNumber =%d, cur_p = %f, p = %f\n",
		  rv->name().c_str(),rv->frame(),nodeNumber,cur_p.val(),p.val());
	} else {
	  DiscRV* drv = (DiscRV*)rv;
	  infoMsg(Giga,"Observed Clique prob non-application of rv %s(%d)=%d, nodeNumber =%d, cur_p = %f, p = %f\n",
		  drv->name().c_str(),drv->frame(),drv->val,nodeNumber,cur_p.val(),p.val());
	}
      }

      if (cur_p.essentially_zero()) {
	zero = true;
	// we've got a zero, so might as well stop here now. 
	break;
      }
    }
  }

  // make sure we have an entry.
  if (cliqueValues.size() == 0) {
    cliqueValues.resize(1);
  }
  // we will only be using one value.
  numCliqueValuesUsed = 1;  

  if (zero) {
    // don't bother checking separators.
    maxCEValue.set_to_zero();
    cliqueValues.ptr[0].p.set_to_zero();
  } else {

    // Now, we check all incoming CE separators, make sure the entry for
    // the current clique value it exists, and if it does, multiply by
    // separator probability checking for zeros. Since the clique
    // is all observed, we know that all separators are observed as well.

    for (unsigned sepNumber=0;sepNumber<origin.ceReceiveSeparators.size();sepNumber++) {
      // get a handy reference to the current separator
      InferenceSeparatorClique& sep = 
	part.separatorCliques[origin.ceReceiveSeparators[sepNumber]];

      unsigned accIndex = 0;
      // separator consists of all observed values. We just multiply
      // in the one entry the separator, that is guaranteed to be at
      // postiion 0,0. The key thing is that we must not do ANY hash
      // lookup as there are no hash tables in a separator that
      // has no hidden variables.
      
      InferenceSeparatorClique::AISeparatorValue& sv
	= sep.separatorValues.ptr[accIndex];
    
      // Where was this allocated?  Search in file for key string
      // "ALLOCATE_REMVALUES_ALL_OBSERVED'

      // if anyone is almost zero, stop right now.a 
      if (sv.remValues.ptr[0].p.essentially_zero()) {
	maxCEValue.set_to_zero();
	cliqueValues.ptr[0].p.set_to_zero();
	return;
      }

      p *= sv.remValues.ptr[0].p;
    }

    // store max (and only) value.
    maxCEValue = p;

    // finally, save the probability
    cliqueValues.ptr[0].p = p;
  }

  if (message(Giga)) {
    infoMsg(Giga,"Inserting New (Observed) Clique Value. prob = %f, sum = %f: ",
	    cliqueValues.ptr[0].p.val(),sumProbabilities().val());
    printRVSetAndValues(stdout,fNodes);
  }

}






void
InferenceMaxClique::ceIterateAssignedNodesCliqueDriven(JT_InferencePartition& part,
						       const unsigned nodeNumber,
						       logpr p)
{

  if (nodeNumber == fSortedAssignedNodes.size()) {
    // time to store maxclique value and total probability, p is
    // current clique probability.

    // Now, we iterate through all incoming CE separators, make sure the entry 
    // for the current clique value it exists, and if it does, multiply by separator probability.
    for (unsigned sepNumber=0;sepNumber<origin.ceReceiveSeparators.size();sepNumber++) {
      // get a handy reference to the current separator
      InferenceSeparatorClique& sep = 
	part.separatorCliques[origin.ceReceiveSeparators[sepNumber]];

      unsigned packedVal[128];
      // If these assertions fail (at some time in the future, probably in
      // the year 2150), then it is fine to increase 128 to something larger.
      // In fact, 128 is so large, lets not even do the assert.
      // assert ( sep.origin.accPacker.packedLen() < 128 );
      // assert ( sep.origin.remPacker.packedLen() < 128 );

      /*
       * There are 3 cases.
       * 1) AI exists and REM exist
       * 2) AI exists and REM doesnt exist
       * 3) AI does not exist, but REM exists
       * AI not exist and REM not exist can't occur.
       */

      unsigned accIndex;
      // TODO: optimize this check away out of loop.
      if (sep.origin.hAccumulatedIntersection.size() > 0) {
	// an accumulated intersection exists.

	sep.origin.accPacker.pack((unsigned**)sep.accDiscreteValuePtrs.ptr,
				  &packedVal[0]);
	unsigned* accIndexp =
	  sep.iAccHashMap.find(&packedVal[0]);

	
	// If it doesn't exist in this separator, then it must have
	// zero probability some where. We therefore do not insert it
	// into this clique, and continue on with next cliuqe value.
	if ( accIndexp == NULL )
	  return;

	accIndex = *accIndexp;

	// TODO: optimize this check out of loop.
	if (sep.remDiscreteValuePtrs.size() == 0) {
	  // 2) AI exists and REM doesnt exist
	  // Then this separator is entirely covered by one or 
	  // more other separators earlier in the order.

	  // go ahead and insert it here to the 1st entry (entry 0).

	  // handy reference for readability.
	  InferenceSeparatorClique::AISeparatorValue& sv
	    = sep.separatorValues.ptr[accIndex];

	  // Multiply in the separator values probability into the clique value's current entry.
	  p *= sv.remValues.ptr[0].p;
	  // done, move on to next separator.
	  continue; 
	}
      } else {
	// no accumulated intersection exists, everything
	// is in the remainder.
	accIndex = 0;
      }


      if (sep.remDiscreteValuePtrs.size() > 0) {
	// if we're here, then we must have some remainder pointers,
	// meaning something in the separator was hidden.


	// Do the remainder exists in this separator.
	// 
	// either:
	//   1) AI exists and REM exist
	//     or
	//   3) AI does not exist (accIndex == 0), but REM exists
	// 

	// keep handy reference for readability.
	InferenceSeparatorClique::AISeparatorValue& sv
	  = sep.separatorValues.ptr[accIndex];

	sep.origin.remPacker.pack((unsigned**)sep.remDiscreteValuePtrs.ptr,
				  &packedVal[0]);

	unsigned* remIndexp =
	  sv.iRemHashMap.find(&packedVal[0]);

	// If it doesn't exist in this separator, then it must have
	// zero probability some where. We therefore do not insert it
	// into this clique, and continue on with next cliuqe value.
	if ( remIndexp == NULL )
	  return;

	// We've finally got the sep entry.  Multiply in the separator
	// values probability into the clique value's current entry.
	p *= sv.remValues.ptr[*remIndexp].p;
      } else {
	// separator consists of all observed values. We just multiply
	// in the one entry the separator, that is guaranteed to be at
	// postiion 0,0. The key thing is that we must not do ANY hash
	// lookup as there are no hash tables in a separator without
	// hidden variables.

	InferenceSeparatorClique::AISeparatorValue& sv
	  = sep.separatorValues.ptr[accIndex];

	// Where was this allocated?  Search in file for key string
	// "ALLOCATE_REMVALUES_ALL_OBSERVED'
	p *= sv.remValues.ptr[0].p;
      }
    }


    // keep track of the max clique probability right here.
    if (p > maxCEValue)
      maxCEValue = p;

    if (numCliqueValuesUsed >= cliqueValues.size()) {
      // TODO: optimize this.
      // cliqueValues.resizeAndCopy(cliqueValues.size()*2);
      if (numCliqueValuesUsed >= origin.cliqueValueSpaceManager.currentSize())
	origin.cliqueValueSpaceManager.advanceToNextSize();
      cliqueValues.resizeAndCopy(origin.cliqueValueSpaceManager.currentSize());
    }

    // TODO: figure out if it is possible to get around doing this
    // check (or to help branch prediction predict it, since it will
    // be different for differnet cliques). Answer: We can do this
    // check once at beginning of iteration of assigned nodes, and
    // have two versions of this code.
    if (origin.packer.packedLen() <= IMC_NWWOH) {
      // pack the clique values directly into place
      origin.packer.pack(
			 (unsigned**)discreteValuePtrs.ptr,
			 (unsigned*)&(cliqueValues.ptr[numCliqueValuesUsed].val[0]));
    } else {
      // Deal with the hash table to re-use clique values.
      // First, grab pointer to storge where next clique value would
      // be stored if it ends up being used.
      unsigned *pcv = origin.valueHolder.curCliqueValuePtr();
      // Next, pack the clique values into this position.
      origin.packer.pack((unsigned**)discreteValuePtrs.ptr,(unsigned*)pcv);
      // Look it up in the hash table.
      bool foundp;
      unsigned *key;
      key = origin.cliqueValueHashSet.insert(pcv,foundp);
      if (!foundp) {
	// if it was not found, need to claim this storage that we
	// just used.
	origin.valueHolder.allocateCurCliqueValue();
      }
      // Save the pointer to whatever the hash table decided to use.
      cliqueValues.ptr[numCliqueValuesUsed].ptr = key;
    }
    // save the probability
    cliqueValues.ptr[numCliqueValuesUsed].p = p;
    numCliqueValuesUsed++;

    // we are now done with this maxclique value.
    return;
  }

  RV* rv = fSortedAssignedNodes[nodeNumber];
  // do the loop right here

  infoMsg(Giga,"Starting assigned iteration of rv %s(%d), nodeNumber=%d, p = %f\n",
	  rv->name().c_str(),rv->frame(),nodeNumber,p.val());

  switch (origin.dispositionSortedAssignedNodes[nodeNumber]) {
  case MaxClique::AN_NOTSEP_PROB_SPARSEDENSE:
  case MaxClique::AN_SEP_PROB_SPARSEDENSE:
    {
      logpr cur_p;
      rv->begin(cur_p);
      do {
	// At each step, we compute probability
	if (message(Giga)) {
	  if (!rv->discrete()) {
	    infoMsg(Giga,"Assigned iteration and prob application of rv %s(%d)=C, nodeNumber =%d, cur_p = %f, p = %f\n",
		    rv->name().c_str(),rv->frame(),nodeNumber,cur_p.val(),p.val());
	  } else {
	    DiscRV* drv = (DiscRV*)rv;
	    infoMsg(Giga,"Assigned CPT iteration and prob application of rv %s(%d)=%d, nodeNumber =%d, cur_p = %f, p = %f\n",
		    rv->name().c_str(),rv->frame(),drv->val,nodeNumber,cur_p.val(),p.val());
	  }
	}
	// if at any step, we get zero, then back out.
	if (!cur_p.essentially_zero()) {
	  // Continue, updating probability by cur_p, contributing this
	  // probability to the clique potential.
	  ceIterateAssignedNodesCliqueDriven(part,nodeNumber+1,p*cur_p);
	}
      } while (rv->next(cur_p));
    }
    break;


  case MaxClique::AN_NOTSEP_NOTPROB_DENSE:
    {
      DiscRV* drv = (DiscRV*)rv;
      // do the loop right here
      drv->val = 0;
      do {
	infoMsg(Giga,"Assigned card iteration of rv %s(%d)=%d, nodeNumber = %d, p = %f\n",
		rv->name().c_str(),rv->frame(),drv->val,nodeNumber,p.val());
	// Continue, do not update probability!!
	ceIterateAssignedNodesCliqueDriven(part,nodeNumber+1,p);
      } while (++drv->val < drv->cardinality);
    }
    break;


  default: // all other cases, we do CPT iteration removing zeros without updating probabilities.
    {
      logpr cur_p;
      rv->begin(cur_p);
      do {
	// At each step, we compute probability
	if (message(Giga)) {
	  if (!rv->discrete()) {
	    infoMsg(Giga,"Assigned iteration of rv %s(%d)=C, nodeNumber =%d, p = %f\n",
		    rv->name().c_str(),rv->frame(),nodeNumber,p.val());
	  } else {
	    DiscRV* drv = (DiscRV*)rv;
	    infoMsg(Giga,"Assigned CPT iteration and zero removal of rv %s(%d)=%d, nodeNumber =%d, cur_p = %f, p = %f\n",
		    rv->name().c_str(),rv->frame(),drv->val,nodeNumber,cur_p.val(),p.val());
	  }
	}
	// if at any step, we get zero, then back out.
	if (!cur_p.essentially_zero()) {
	  // Continue, do not update probability!!
	  ceIterateAssignedNodesCliqueDriven(part,nodeNumber+1,p);
	}
      } while (rv->next(cur_p));
    }
    break;
  }

}




void
InferenceMaxClique::ceIterateUnassignedNodesCliqueDriven(JT_InferencePartition& part,
							 const unsigned nodeNumber,
							 const logpr p)
{
  if (nodeNumber == fUnassignedNodes.size()) {
    ceIterateAssignedNodesCliqueDriven(part,0,p);
    return;
  }
  RV* rv = fUnassignedNodes[nodeNumber];
  infoMsg(Giga,"Starting Unassigned iteration of rv %s(%d), nodeNumber = %d, p = %f\n",
	  rv->name().c_str(),rv->frame(),nodeNumber,p.val());

  if (rv->hidden()) {
    // only discrete RVs can be hidden for now.
    DiscRV* drv = (DiscRV*)rv;
    // do the loop right here
    drv->val = 0;
    do {
      infoMsg(Giga,"Unassigned iteration of rv %s(%d)=%d, nodeNumber = %d, p = %f\n",
	      rv->name().c_str(),rv->frame(),drv->val,nodeNumber,p.val());
      // continue on, effectively multiplying p by unity.
      ceIterateUnassignedNodesCliqueDriven(part,nodeNumber+1,p);
    } while (++drv->val < drv->cardinality);
  } else {
    // observed, either discrete or continuous
    if (rv->discrete()) {
      // DiscRV* drv = (DiscRV*)rv;
      // TODO: for observed variables, do this once at the begining
      // before any looping here.
      // drv->setToObservedValue();
      infoMsg(Giga,"Unassigned iteration of rv %s(%d)=%d, nodeNumber = %d, p = %f\n",
	      rv->name().c_str(),rv->frame(),RV2DRV(rv)->val,nodeNumber,p.val());
    } else {
      // nothing to do since we get continuous observed
      // value indirectly
      infoMsg(Giga,"Unassigned iteration of rv %s(%d)=C, nodeNumber = %d, p = %f\n",
	      rv->name().c_str(),rv->frame(),nodeNumber,p.val());
    }
    // continue on, effectively multiplying p by unity.
    ceIterateUnassignedNodesCliqueDriven(part,nodeNumber+1,p);
  }
}






logpr
InferenceMaxClique::
sumProbabilities()
{
  logpr p;
  if (numCliqueValuesUsed > 0) {
    // We directly assign first one rather than adding to initialized
    // zero so that logpr's log(0) floating point value is preserved.
    p = cliqueValues.ptr[0].p;
    for (unsigned i=1;i<numCliqueValuesUsed;i++)
      p += cliqueValues.ptr[i].p;
  }
  return p;
}


void
InferenceMaxClique::
emIncrement(const logpr probE,
	    const bool localCliqueNormalization,
	    const double emTrainingBeam)
{
  // recompute here each time, shouldn't take too long
  // TODO: re-compute this once for each inference clique.
  sArray< RV*> fAssignedProbNodes;
  unsigned numAssignedProbNodes = 0;
  for (unsigned nodeNumber = 0; nodeNumber < fSortedAssignedNodes.size(); nodeNumber ++ ) {
    if (origin.dispositionSortedAssignedNodes[nodeNumber] == MaxClique::AN_CPT_ITERATION_COMPUTE_AND_APPLY_PROB ||
	origin.dispositionSortedAssignedNodes[nodeNumber] == MaxClique::AN_COMPUTE_AND_APPLY_PROB) {
      numAssignedProbNodes++;
    }
  }
  if (numAssignedProbNodes == 0)
    return; // nothing to do for this clique

  infoMsg(Huge,"EM Incrementing for clique\n");

  fAssignedProbNodes.resize(numAssignedProbNodes);
  numAssignedProbNodes = 0;
  for (unsigned nodeNumber = 0; nodeNumber < fSortedAssignedNodes.size(); nodeNumber ++ ) {
    if (origin.dispositionSortedAssignedNodes[nodeNumber] == MaxClique::AN_CPT_ITERATION_COMPUTE_AND_APPLY_PROB ||
	origin.dispositionSortedAssignedNodes[nodeNumber] == MaxClique::AN_COMPUTE_AND_APPLY_PROB) {
      RV* rv = fSortedAssignedNodes[nodeNumber];
      fAssignedProbNodes[numAssignedProbNodes++] = rv;
    }
  }

  logpr locProbE((void*)0);
  if (localCliqueNormalization) {
    locProbE = sumProbabilities();
  } else {
    locProbE = probE;
  }

  const bool imc_nwwoh_p = (origin.packer.packedLen() <= IMC_NWWOH);
  const bool clique_has_hidden_vars = (origin.hiddenNodes.size() > 0);
  
  logpr beamThreshold((void*)0);
  if (origin.cliqueBeam != (-LZERO)) {
    // then we do clique table pruning right here rather
    // than a separate call to ceCliquePrune().
    // break into the logp to avoid unnecessary zero checking.
    beamThreshold.valref() = - emTrainingBeam;
  } else {
    // set beam threshold to a value that will never cause pruning.
    beamThreshold.set_to_zero();
  }
  // create unnormalized beam threshold.
  beamThreshold *= locProbE;

  // now go through updating each thing
  for (unsigned cvn=0;cvn<numCliqueValuesUsed;cvn++) {

    // EM pruning here based on unnormalized posterior. Don't bother
    // with things that are below threshold.
    if (cliqueValues.ptr[cvn].p < beamThreshold)
      continue;

    // if still here, then create the posterior to update the
    // parameters.
    logpr posterior = cliqueValues.ptr[cvn].p/locProbE;

    if (clique_has_hidden_vars) {
      if (imc_nwwoh_p) {
	origin.packer.unpack((unsigned*)&(cliqueValues.ptr[cvn].val[0]),
			     (unsigned**)discreteValuePtrs.ptr);
      } else {
	origin.packer.unpack((unsigned*)cliqueValues.ptr[cvn].ptr,
			     (unsigned**)discreteValuePtrs.ptr);
      }
    }
    // Increment all assigned probability nodes.  
    // 
    // TODO: integerate out all but cont. variables parents so we
    // don't multiply increment those varialbes for the same parent
    // values (waisting time).
    for (unsigned nodeNumber = 0; nodeNumber < fAssignedProbNodes.size(); nodeNumber ++ ) {
      RV* rv = fAssignedProbNodes[nodeNumber];
      rv->emIncrement(posterior);
    }
  }

}



/*
 * we have now a fully instantiated clique and are ready for backwards
 * pass. Iterate through the values that are above beam and
 * instantiate the outgoing separator with those values.
 *
 */

void 
InferenceMaxClique::
deReceiveFromIncommingSeparator(JT_InferencePartition& part)
{
  deReceiveFromIncommingSeparator(part,
				  part.separatorCliques[origin.ceSendSeparator]);
}

// For each clique value, we need to look up appropriate value in the
// separator and multiply it into the current clique probability.
void 
InferenceMaxClique::
deReceiveFromIncommingSeparator(JT_InferencePartition& part,
				InferenceSeparatorClique& sep)
{
  if (JunctionTree::viterbiScore) {
    return deReceiveFromIncommingSeparatorViterbi(part,sep);
  }

  if (origin.hiddenNodes.size() == 0) {
    // do the observed clique case up front right here so we don't
    // need to keep checking below.
    InferenceSeparatorClique::AISeparatorValue& sv
      = sep.separatorValues.ptr[0];
    cliqueValues.ptr[0].p *= sv.remValues.ptr[0].bp();
    return;
  }

  // allocate some temporary storage for packed separator values.
  // 128 words is *much* bigger than any possible packed clique value
  // will take on, but it is easy/fast to allocate on the stack right now.
  unsigned packedVal[128];
  // but just in case, we assert.
  assert ((sep.origin.hAccumulatedIntersection.size() == 0)
	  ||
	  (sep.origin.accPacker.packedLen() < 128)
	  );
  assert ((sep.origin.hRemainder.size() == 0) 
	  ||
	  (sep.origin.remPacker.packedLen() < 128 )
	  );
  // If this assertion fails (at some time in the future, probably in
  // the year 2150), then it is fine to increase 128 to something larger.

  // cache check here.
  const bool imc_nwwoh_p = (origin.packer.packedLen() <= IMC_NWWOH);
  for (unsigned cvn=0;cvn<numCliqueValuesUsed;cvn++) {{

    // TODO: optimize away this conditional check.
    if (imc_nwwoh_p) {
      origin.packer.unpack((unsigned*)&(cliqueValues.ptr[cvn].val[0]),
			   (unsigned**)discreteValuePtrs.ptr);
    } else {
      origin.packer.unpack((unsigned*)cliqueValues.ptr[cvn].ptr,
			   (unsigned**)discreteValuePtrs.ptr);
    }

    // All hidden random variables now have their discrete
    // value. Get appropriate entry in the separator
    // for this clique entry.

    /*
     * There are 3 cases.
     * 1) AI exists and REM exist
     * 2) AI exists and REM doesnt exist
     * 3) AI does not exist, but REM exists
     * AI not exist and REM not exist can't occur.
    */

    unsigned accIndex;
    // TODO: optimize this check away out of loop.
    if (sep.origin.hAccumulatedIntersection.size() > 0) {
      // an accumulated intersection exists.

      sep.origin.accPacker.pack((unsigned**)sep.accDiscreteValuePtrs.ptr,
				&packedVal[0]);
      unsigned* accIndexp =
	sep.iAccHashMap.find(&packedVal[0]);

      // we should always find something or else something is wrong.
      assert ( accIndexp != NULL ); 
      accIndex = *accIndexp;

      // TODO: optimize this check out of loop.
      if (sep.remDiscreteValuePtrs.size() == 0) {
	// 2) AI exists and REM doesnt exist
	// Then this separator is entirely covered by one or 
	// more other separators earlier in the order.

	// go ahead and insert it here to the 1st entry (entry 0).

	// handy reference for readability.
	InferenceSeparatorClique::AISeparatorValue& sv
	  = sep.separatorValues.ptr[accIndex];

	// Multiply in this separator value's probability.
	cliqueValues.ptr[cvn].p *= sv.remValues.ptr[0].bp();
	// done
	goto next_iteration;
      }
    } else {
      // no accumulated intersection exists, everything
      // is in the remainder.
      accIndex = 0;
    }

    if (sep.remDiscreteValuePtrs.size() > 0) {
      // if we're here, then we must have some remainder pointers.
      // Do the remainder exists in this separator.
      // either:
      //   1) AI exists and REM exist
      //     or
      //   3) AI does not exist (accIndex == 0), but REM exists
      // 

      // keep handy reference for readability.
      InferenceSeparatorClique::AISeparatorValue& sv
	= sep.separatorValues.ptr[accIndex];

      sep.origin.remPacker.pack((unsigned**)sep.remDiscreteValuePtrs.ptr,
				&packedVal[0]);

      unsigned* remIndexp =
	sv.iRemHashMap.find(&packedVal[0]);

      // it must exist
      // assert ( remIndexp != NULL );
      if ( remIndexp == NULL ) {
	// NOTE: we could prune this clique entry away, but for now we just issue the error.  
	error("ERROR: During distribute evidence, found clique entry without corresponding incomming separator entry. Try increasing -sbeam beam, or just use -cbeam without the -sbeam option.\n");
      }

      // We've finally got the sep entry. Multiply it it into the
      // current clique value.
      cliqueValues.ptr[cvn].p *= sv.remValues.ptr[*remIndexp].bp();
    } else {
      // Either separator is all observed, or the separator
      // is completely contained in the accumulated intersection.
      // In either case, we multiply in its one value.

      InferenceSeparatorClique::AISeparatorValue& sv
	= sep.separatorValues.ptr[accIndex];

      if (sv.numRemValuesUsed == 1) {
	// We've finally got the sep entry. Multiply it it into the
	// current clique value.
	cliqueValues.ptr[cvn].p *= sv.remValues.ptr[0].bp();
      } else {
	// NOTE: we could prune this clique entry away, but for now we just issue the error.  
	error("ERROR: During distribute evidence, found clique entry without corresponding incomming separator entry. Try increasing -sbeam beam, or just use -cbeam without the -sbeam option.\n");
      }
    }

  }
  next_iteration:
  ;    
  }

  // TODO: backwards/distribute evidence beam pruning here.
}


// Viterbi version of above.
// For each clique value, we need to look up appropriate value in the
// separator and multiply it into the current clique probability.
void 
InferenceMaxClique::
deReceiveFromIncommingSeparatorViterbi(JT_InferencePartition& part,
				       InferenceSeparatorClique& sep)
{

  if (origin.hiddenNodes.size() == 0) {
    // All clique values already set to their max (and only) settings,
    // so no need to do anything.
    return;
  }


  // cache check here.
  const bool imc_nwwoh_p = (origin.packer.packedLen() <= IMC_NWWOH);
  // grab backpointer from seps forward pointer entry directly.
  InferenceSeparatorClique::AISeparatorValue& sv
    = sep.separatorValues.ptr[sep.forwPointer.viterbiAccIndex];
  unsigned cvn = sv.remValues.ptr[sep.forwPointer.viterbiRemIndex].backPointer;

  if (imc_nwwoh_p) {
    origin.packer.unpack((unsigned*)&(cliqueValues.ptr[cvn].val[0]),
			 (unsigned**)discreteValuePtrs.ptr);
  } else {
    origin.packer.unpack((unsigned*)cliqueValues.ptr[cvn].ptr,
			 (unsigned**)discreteValuePtrs.ptr);
  }

#if 0

  // allocate some temporary storage for packed separator values.
  // 128 words is *much* bigger than any possible packed clique value
  // will take on, but it is easy/fast to allocate on the stack right now.
  unsigned packedVal[128];
  // but just in case, we assert.
  assert ((sep.origin.hAccumulatedIntersection.size() == 0)
	  ||
	  (sep.origin.accPacker.packedLen() < 128)
	  );
  assert ((sep.origin.hRemainder.size() == 0) 
	  ||
	  (sep.origin.remPacker.packedLen() < 128 )
	  );
  // If this assertion fails (at some time in the future, probably in
  // the year 2150), then it is fine to increase 128 to something larger.



  // We assume here that the RVs according to the separator are
  // already set to the appropriate value. Get appropriate sep entry
  // and fill in rest of clique RVs with that value.

  /*
   * There are 3 cases.
   * 1) AI exists and REM exist
   * 2) AI exists and REM doesnt exist
   * 3) AI does not exist, but REM exists
   * AI not exist and REM not exist can't occur.
   */

  unsigned accIndex,cvn;
  // TODO: optimize this check away out of loop.
  if (sep.origin.hAccumulatedIntersection.size() > 0) {
    // an accumulated intersection exists.

    sep.origin.accPacker.pack((unsigned**)sep.accDiscreteValuePtrs.ptr,
			      &packedVal[0]);
    unsigned* accIndexp =
      sep.iAccHashMap.find(&packedVal[0]);

    // we should always find something or else something is wrong.
    assert ( accIndexp != NULL ); 
    accIndex = *accIndexp;

    // TODO: optimize this check out of loop.
    if (sep.remDiscreteValuePtrs.size() == 0) {
      // 2) AI exists and REM doesnt exist
      // Then this separator is entirely covered by accumulated
      // intersection (one or more other separators earlier in the
      // order).

      // go ahead and insert it here to the 1st entry (entry 0).

      // handy reference for readability.
      InferenceSeparatorClique::AISeparatorValue& sv
	= sep.separatorValues.ptr[accIndex];

      // Multiply in this separator value's probability.

      cvn = sv.remValues.ptr[0].backPointer;

      // unpack the rest of the clique.
      if (imc_nwwoh_p) {
	origin.packer.unpack((unsigned*)&(cliqueValues.ptr[cvn].val[0]),
			     (unsigned**)discreteValuePtrs.ptr);
      } else {
	origin.packer.unpack((unsigned*)cliqueValues.ptr[cvn].ptr,
			     (unsigned**)discreteValuePtrs.ptr);
      }

      // done
      goto done;
    }
  } else {
    // no accumulated intersection exists, everything
    // is in the remainder.
    accIndex = 0;
  }

  if (sep.remDiscreteValuePtrs.size() > 0) {
    // if we're here, then we must have some remainder pointers.
    // Do the remainder exists in this separator.
    // either:
    //   1) AI exists and REM exist
    //     or
    //   3) AI does not exist (accIndex == 0), but REM exists
    // 

    // keep handy reference for readability.
    InferenceSeparatorClique::AISeparatorValue& sv
      = sep.separatorValues.ptr[accIndex];

    sep.origin.remPacker.pack((unsigned**)sep.remDiscreteValuePtrs.ptr,
			      &packedVal[0]);

    unsigned* remIndexp =
      sv.iRemHashMap.find(&packedVal[0]);

    // it must exist
    assert ( remIndexp != NULL );

    // We've finally got the sep entry. Unpack the pointed-to clique
    // entry.

    cvn = sv.remValues.ptr[*remIndexp].backPointer;
    // unpack the rest of the clique.
    if (imc_nwwoh_p) {
      origin.packer.unpack((unsigned*)&(cliqueValues.ptr[cvn].val[0]),
			   (unsigned**)discreteValuePtrs.ptr);
    } else {
      origin.packer.unpack((unsigned*)cliqueValues.ptr[cvn].ptr,
			   (unsigned**)discreteValuePtrs.ptr);
    }



  } else {
    // Either separator is all observed, or the separator
    // is completely contained in the accumulated intersection.
    // In either case, we use its one value.

    InferenceSeparatorClique::AISeparatorValue& sv
      = sep.separatorValues.ptr[accIndex];

    // We've finally got the sep entry. Unpack the pointed-to clique
    // entry.
    cvn = sv.remValues.ptr[0].backPointer;
    if (imc_nwwoh_p) {
      origin.packer.unpack((unsigned*)&(cliqueValues.ptr[cvn].val[0]),
			   (unsigned**)discreteValuePtrs.ptr);
    } else {
      origin.packer.unpack((unsigned*)cliqueValues.ptr[cvn].ptr,
			   (unsigned**)discreteValuePtrs.ptr);
    }
  }
 done:
  ;
#endif
}



// Scatter out to the outgoing separators which are the same as the
// receive separators in the collect evidence stage, so we use that
// array here.
void 
InferenceMaxClique::
deScatterToOutgoingSeparators(JT_InferencePartition& part)
{

  if (JunctionTree::viterbiScore) {
    // While we might think that there is nothing to do in this case
    // since if the RVs associated with the current clique have been
    // set to the appropriate clique table entry, the associated
    // separator RVs have also been set. But in some cases (namely the
    // island algorithm) it is possible for the separator to get
    // changed by an additional forward pass in which case we would be
    // using the wrong value.

    // Specifically, to get island decoding working, we need to have
    // separator keep track of which entry is current max rather than
    // having it assume that its values are set from right most
    // clique.
    // Consider the following:
    //  
    //    C0 -- s0 -- C1 -- s1 -- C2 -- s2 -- C3 --
    //             1>    2>    3>    <4          <b  
    // 
    //  1> and 3> are ceGatherIntoRoot msgs
    //  2> is a ceSendToNextPartition msg
    //  <4 is a deReceveFromIncommingSep msg
    //      and deScatterOutofRoot is (currently) a nop
    //  <b is a deReceveFromIncommingSep msg from before.
    // 
    // The potential problem is that <4 might assume that it's 's2' is currently
    // set to appropriate values (max of C3 done by <b), but message 3>
    // changes clique C2 and thus changes s2 so that it no longer
    // holds max value.
    // 
    // One solution is to add two variables in InferenceSep, i.e., 
    // two indices to get at the veterbi values for current separator.
    //      unsigned viterbiAccIndex;
    //      unsigned viterbiRemIndex;
    // to keep track of which entries in seps (such as s2) are max.
    // Then <4 would restore those values for each separator.  and
    // deScatterOutofRoot would compute them (since deScatterOutofRoot
    // is called at the time C3 is at its true max value). Note that
    // this works between partitions since the separator between
    // two partitions is contained in the right partition.
    //
    // Note also that this can probably be easily extended to the
    // N-best list case.
    deScatterToOutgoingSeparatorsViterbi(part);
    return; 
  }

  if (origin.ceReceiveSeparators.size() == 0)
    return;


  // Note. All separator .bp values have already been initialized to
  // zero when the structure containing them 'a RemainderValue' was
  // constructed. All memory reallocations will have preserved these
  // initializations, so there is no need to scan through initializing
  // bp to zero here.


  if (origin.hiddenNodes.size() == 0) {
    // Do the observed clique case up front right here so we don't
    // need to keep checking below. Here, the clique is observed which
    // means that all connecting separators are also observed. We just
    // sweep through all separators updating the values.
    for (unsigned sepNumber=0;sepNumber<origin.ceReceiveSeparators.size();sepNumber++) {
      InferenceSeparatorClique& sep = 
	part.separatorCliques[origin.ceReceiveSeparators[sepNumber]];
      InferenceSeparatorClique::AISeparatorValue& sv
	= sep.separatorValues.ptr[0];
      // can use assignment rather than += here since there is only one value.
      sv.remValues.ptr[0].bp() = cliqueValues.ptr[0].p;      
    }
  } else {

    // allocate some temporary storage for packed separator values.
    // 128 words is *much* bigger than any possible packed clique value
    // will take on, but it is easy/fast to allocate on the stack right now.
    unsigned packedVal[128];

    // cache check here.
    const bool imc_nwwoh_p = (origin.packer.packedLen() <= IMC_NWWOH); 
    for (unsigned cvn=0;cvn<numCliqueValuesUsed;cvn++) {

      // TODO: optimize away this conditional check.
      if (imc_nwwoh_p) {
	origin.packer.unpack((unsigned*)&(cliqueValues.ptr[cvn].val[0]),
			     (unsigned**)discreteValuePtrs.ptr);
      } else {
	origin.packer.unpack((unsigned*)cliqueValues.ptr[cvn].ptr,
			     (unsigned**)discreteValuePtrs.ptr);
      }

      // now we iterate through all the separators.
      for (unsigned sepNumber=0;sepNumber<origin.ceReceiveSeparators.size();sepNumber++) {
	// get a handy reference to the current separator
	InferenceSeparatorClique& sep = 
	  part.separatorCliques[origin.ceReceiveSeparators[sepNumber]];

	// If these assertions fail (at some time in the future, probably in
	// the year 2150), then it is fine to increase 128 to something larger.
	// In fact, 128 is so large, lets not even do the assert.
	// assert ( sep.origin.accPacker.packedLen() < 128 );
	// assert ( sep.origin.remPacker.packedLen() < 128 );


	/*
	 * There are 3 cases.
	 * 1) AI exists and REM exist
	 * 2) AI exists and REM doesnt exist
	 * 3) AI does not exist, but REM exists
	 * AI not exist and REM not exist can't occur.
	 */

	unsigned accIndex;
	// TODO: optimize this check away out of loop.
	if (sep.origin.hAccumulatedIntersection.size() > 0) {
	  // an accumulated intersection exists.

	  sep.origin.accPacker.pack((unsigned**)sep.accDiscreteValuePtrs.ptr,
				    &packedVal[0]);
	  unsigned* accIndexp =
	    sep.iAccHashMap.find(&packedVal[0]);

	  // we should always find something or else something is wrong.
	  assert ( accIndexp != NULL ); 
	  accIndex = *accIndexp;

	  // TODO: optimize this check out of loop.
	  if (sep.remDiscreteValuePtrs.size() == 0) {
	    // 2) AI exists and REM doesnt exist
	    // Then this separator is entirely covered by one or 
	    // more other separators earlier in the order.

	    // go ahead and insert it here to the 1st entry (entry 0).

	    // handy reference for readability.
	    InferenceSeparatorClique::AISeparatorValue& sv
	      = sep.separatorValues.ptr[accIndex];

	    // Add in this clique value's probability.  Note that bp was
	    // initialized during forward pass.
	    sv.remValues.ptr[0].bp() += cliqueValues.ptr[cvn].p;
	    // done, move on to next separator.
	    continue; 
	  }
	} else {
	  // no accumulated intersection exists, everything
	  // is in the remainder.
	  accIndex = 0;
	}

	if (sep.remDiscreteValuePtrs.size() > 0) {
	  // if we're here, then we must have some remainder
	  // pointers.

	  // Do the remainder exists in this separator.
	  // 
	  // either:
	  //   1) AI exists and REM exist
	  //     or
	  //   3) AI does not exist (accIndex == 0), but REM exists
	  // 
	
	  // keep handy reference for readability.
	  InferenceSeparatorClique::AISeparatorValue& sv
	    = sep.separatorValues.ptr[accIndex];
	
	  sep.origin.remPacker.pack((unsigned**)sep.remDiscreteValuePtrs.ptr,
				    &packedVal[0]);

	  unsigned* remIndexp =
	    sv.iRemHashMap.find(&packedVal[0]);

	  // it must exist
	  assert ( remIndexp != NULL );
	
	  // We've finally got the sep entry.  Add in this clique value's
	  // probability.  Note that bp was initialized during forward
	  // pass.
	  sv.remValues.ptr[*remIndexp].bp() += cliqueValues.ptr[cvn].p;
	} else {
	  // Either separator is all observed, or the separator
	  // is completely contained in the accumulated intersection.
	
	  // keep handy reference for readability.
	  InferenceSeparatorClique::AISeparatorValue& sv
	    = sep.separatorValues.ptr[accIndex];

	  // We've finally got the sep entry.  Add in this clique value's
	  // probability.  Note that bp was initialized during forward
	  // pass.
	  sv.remValues.ptr[0].bp() += cliqueValues.ptr[cvn].p;
	}
      }
    }
  }

  // lastly iterate through all separators, and all entries in
  // each separator and do the "divide" (subtraction)
  for (unsigned sepNumber=0;sepNumber<origin.ceReceiveSeparators.size();sepNumber++) {
    // get a handy reference to the current separator
    InferenceSeparatorClique& sep = 
      part.separatorCliques[origin.ceReceiveSeparators[sepNumber]];

    for (unsigned aiNo=0;aiNo < sep.numSeparatorValuesUsed; aiNo ++) {
      InferenceSeparatorClique::AISeparatorValue* aisep = &(sep.separatorValues.ptr[aiNo]);
      for (unsigned remNo=0; remNo < aisep->numRemValuesUsed; remNo++) {
	InferenceSeparatorClique::RemainderValue* rv = &(aisep->remValues.ptr[remNo]);
	// We remove p from bp since bp will already have a factor of
	// p in it. We do this by dividing it out.
	// -
	// We must make sure that if CE stage is zero, we do not run
	// DE stage, as in that case it might be the case that rv->p
	// == 0.
	// -
	// In the normal case (CE != 0), we could do direct value
	// reference subtraction in log domain (corresponding to
	// divison in original domain) to ensure that compiler creates
	// no temporaries. In other words, this operation could be
	// "rv->bp = rv->bp / rv->p;" rv->bp.valref() =
	// rv->bp.valref() - rv->p.valref(); Do slower version for now
	// until debugged: We assume here that (!rv->p.zero()) is true
	// since we pruned all zero p's above. If we didn't prune,
	// then rv->p == zero would imply that rv->bp == zero, and we
	// would need to do a check. Note that this pruning always
	// occurs, regardless of beam.
	rv->bp() /= rv->p;
      }
    }
  }
}


// Viterbi version. 
// Since all sep value are currently set to appropriate max
// value from cur clique, through each sep and find the index to
// re-set the sep to that value so it can be restored later. We
// do this since the seps otherwise might change their value.
void 
InferenceMaxClique::
deScatterToOutgoingSeparatorsViterbi(JT_InferencePartition& part)
{

  // allocate some temporary storage for packed separator values.
  // 128 words is *much* bigger than any possible packed clique value
  // will take on, but it is easy/fast to allocate on the stack right now.
  unsigned packedVal[128];

  // now we iterate through all the separators.
  for (unsigned sepNumber=0;sepNumber<origin.ceReceiveSeparators.size();sepNumber++) {
    // get a handy reference to the current separator
    InferenceSeparatorClique& sep = 
      part.separatorCliques[origin.ceReceiveSeparators[sepNumber]];
    
    /*
     * There are 3 cases.
     * 1) AI exists and REM exist
     * 2) AI exists and REM doesnt exist
     * 3) AI does not exist, but REM exists
     * AI not exist and REM not exist can't occur.
     */

    unsigned accIndex;
    // TODO: optimize this check away out of loop.
    if (sep.origin.hAccumulatedIntersection.size() > 0) {
      // an accumulated intersection exists.

      sep.origin.accPacker.pack((unsigned**)sep.accDiscreteValuePtrs.ptr,
				&packedVal[0]);
      unsigned* accIndexp =
	sep.iAccHashMap.find(&packedVal[0]);

      // we should always find something or else something is wrong.
      assert ( accIndexp != NULL ); 
      accIndex = *accIndexp;

    } else {
      // no accumulated intersection exists, everything
      // is in the remainder.
      accIndex = 0;
    }

    if (sep.remDiscreteValuePtrs.size() == 0) {
      // 2) AI exists and REM doesnt exist
      // Then this separator is entirely covered by one or 
      // more other separators earlier in the order.
      sep.forwPointer.viterbiAccIndex = accIndex;
      sep.forwPointer.viterbiRemIndex = 0;
    } else {
      // if we're here, then we must have some remainder
      // pointers.

      // Do the remainder exists in this separator.
      // 
      // either:
      //   1) AI exists and REM exist
      //     or
      //   3) AI does not exist (accIndex == 0), but REM exists
      // 
	
      // keep handy reference for readability.
      InferenceSeparatorClique::AISeparatorValue& sv
	= sep.separatorValues.ptr[accIndex];
	
      sep.origin.remPacker.pack((unsigned**)sep.remDiscreteValuePtrs.ptr,
				&packedVal[0]);

      unsigned* remIndexp =
	sv.iRemHashMap.find(&packedVal[0]);

      // it must exist
      assert ( remIndexp != NULL );
	
      // We've finally got the sep entry.  Store the sep entry's id.
      sep.forwPointer.viterbiAccIndex = accIndex;
      sep.forwPointer.viterbiRemIndex = *remIndexp;
    }
  }
}



/*-
 *-----------------------------------------------------------------------
 * InferenceMaxClique::setCliqueToMaxCliqueValue()
 *
 *    This routine just sets the clique to the clique value to the one
 *    that has maximum score. It sets all RVs to the values associated
 *    with the clique value that has maximum score.
 *    
 * Preconditions:
 *
 *     Clique table must be at least partially instantiated.
 *
 * Postconditions:
 *
 *     Hidden RVs now have values of clique value having maximum score
 *     in clique table.
 *
 * Side Effects:
 *
 *    Changes values of RVs associated with this clique
 *
 * Results:
 *    nothing
 *
 *-----------------------------------------------------------------------
 */
logpr
InferenceMaxClique::
setCliqueToMaxCliqueValue()
{

  if (origin.hiddenNodes.size() == 0) {
    // The observed clique case requires no action since this
    // means that the cliuqe (and therefore all its separators)
    // are all observed and already set to their max prob (and only) values.
    return cliqueValues.ptr[0].p;
  } else {
    unsigned max_cvn = 0;
    logpr max_cvn_score = cliqueValues.ptr[0].p;
    
    // find the max score clique entry
    for (unsigned cvn=1;cvn<numCliqueValuesUsed;cvn++) {
      if (cliqueValues.ptr[cvn].p > max_cvn_score) {
	max_cvn_score = cliqueValues.ptr[cvn].p;
	max_cvn = cvn;
      }
    }

    const bool imc_nwwoh_p = (origin.packer.packedLen() <= IMC_NWWOH); 
    if (imc_nwwoh_p) {
      origin.packer.unpack((unsigned*)&(cliqueValues.ptr[max_cvn].val[0]),
			   (unsigned**)discreteValuePtrs.ptr);
    } else {
      origin.packer.unpack((unsigned*)cliqueValues.ptr[max_cvn].ptr,
			   (unsigned**)discreteValuePtrs.ptr);
    }

    return max_cvn_score;
  }

}



////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        SeparatorClique support
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Constructors
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

SeparatorClique::SeparatorClique(MaxClique& c1, MaxClique& c2)
  :  separatorValueSpaceManager(1,     // starting size
				2.0,   // growth rate
				1,     // growth addition
				0.90),  // decay rate 
     remainderValueSpaceManager(1,     // starting size
				2.0,   // growth rate
				1,     // growth addition
				0.90)  // decay rate
{
  nodes.clear();

  // create a set of nodes that is the intersection of the two
  set_intersection(c1.nodes.begin(),c1.nodes.end(),
		   c2.nodes.begin(),c2.nodes.end(),
		   inserter(nodes,nodes.end()));
  assert (nodes.size() > 0);
  
}


SeparatorClique::SeparatorClique(SeparatorClique& from_sep,
				 vector <RV*>& newRvs,
				 map < RVInfo::rvParent, unsigned >& ppf,
				 const unsigned int frameDelta)
  :  separatorValueSpaceManager(1,     // starting size
				2.0,   // growth rate
				1,     // growth addition
				0.90),  // decay rate 
     remainderValueSpaceManager(1,     // starting size
				2.0,   // growth rate
				1,     // growth addition
				0.90)  // decay rate 
{
  set<RV*>::iterator it;

  // clone over nodes RVs and accumulated intersection.
  for (it = from_sep.nodes.begin();
       it != from_sep.nodes.end();
       it++) {
    RV* rv = (*it);
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    

    if ( ppf.find(rvp) == ppf.end() ) {
      coredump("INTERNAL ERROR: can't find rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
    }

    RV* nrv = newRvs[ppf[rvp]];
    nodes.insert(nrv);
  }

  // copy over accumulated intersection
  for (it=from_sep.accumulatedIntersection.begin();
       it != from_sep.accumulatedIntersection.end();
       it++) {
    RV* rv = (*it);
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    


    if ( ppf.find(rvp) == ppf.end() ) {
      coredump("INTERNAL ERROR: can't find rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
    }

    RV* nrv = newRvs[ppf[rvp]];
    accumulatedIntersection.insert(nrv);
  }

  // and 'remainder'
  for (it=from_sep.remainder.begin();
       it != from_sep.remainder.end();
       it++) {
    RV* rv = (*it);
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    

    if ( ppf.find(rvp) == ppf.end() ) {
      coredump("INTERNAL ERROR: can't find rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
    }

    RV* nrv = newRvs[ppf[rvp]];
    remainder.insert(nrv);
  }


}


/*-
 *-----------------------------------------------------------------------
 * SeparatorClique::prepareForUnrolling()
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
SeparatorClique::prepareForUnrolling()
{

  // set up number of hidden ndoes within accumulated intersection.

  // create a vector form of the variables.
  set<RV*>::iterator it;
  for (it = accumulatedIntersection.begin();
       it != accumulatedIntersection.end();
       it++) {
    RV* rv = (*it);
    if (rv->hidden())
      hAccumulatedIntersection.push_back(rv);
  }

  if (hAccumulatedIntersection.size() > 0) {
    // we only use these structures if there is an intersection.
    new (&accPacker) PackCliqueValue(hAccumulatedIntersection);
    assert( accPacker.packedLen() > 0);
    if (accPacker.packedLen() > ISC_NWWOH_AI) {
      // only setup hash table if the packed accumulated insersection
      // set is larger than one machine word (unsigned).
      new (&accValueHolder) CliqueValueHolder(accPacker.packedLen(),
					      // TODO: optimize this 1000 value.
					      AI_SEP_VALUE_HOLDER_STARTING_SIZE,
					      AI_SEP_VALUE_HOLDER_GROWTH_RATE); // 1.25
      // TODO: optimize starting size.
      new (&accSepValHashSet) vhash_set< unsigned > (accPacker.packedLen(),AI_SEP_VALUE_HOLDER_STARTING_SIZE);
    }
  }

  // create a vector form of the variables.
  for (it = remainder.begin();
       it != remainder.end();
       it++) {
    RV* rv = (*it);
    if (rv->hidden())
      hRemainder.push_back(rv);
  }

  if (hRemainder.size() > 0) {
    new (&remPacker) PackCliqueValue(hRemainder);
    assert (remPacker.packedLen() > 0);
    if (remPacker.packedLen() > ISC_NWWOH_RM) { 
      // Only setup hash table if the packed remainder set is larger
      // than one machine word (unsigned).
      new (&remValueHolder) CliqueValueHolder(remPacker.packedLen(),
					      // TODO: optimize this starting sizse
					      REM_SEP_VALUE_HOLDER_STARTING_SIZE, // 2
					      REM_SEP_VALUE_HOLDER_GROWTH_RATE); // 1.25
      new (&remSepValHashSet) vhash_set< unsigned > (remPacker.packedLen(),REM_SEP_VALUE_HOLDER_STARTING_SIZE);
    }
  }

}






/*-
 *-----------------------------------------------------------------------
 * SeparatorClique::printAllJTInfo()
 *   
 *   prints everything JT (not not inference) about this clique to file.
 *
 * Preconditions:
 *   all variables must have been set up. prepareForUnrolling() must have
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
SeparatorClique::printAllJTInfo(FILE*f)
{
  // TODO: also print out nubmer of bits for acc and rem.
  fprintf(f,"Separator information: %d acc packed bits (%d words), %d rem packed bits (%d words)\n",
	  accPacker.packedLenBits(),accPacker.packedLen(),
	  remPacker.packedLenBits(),remPacker.packedLen());


  fprintf(f,"%d Nodes: ",nodes.size()); printRVSet(f,nodes);
  fprintf(f,"%d Acc Inter: ",accumulatedIntersection.size()); printRVSet(f,accumulatedIntersection);  
  fprintf(f,"%d Hid Acc Inter: ",hAccumulatedIntersection.size()); printRVSet(f,hAccumulatedIntersection);  
  fprintf(f,"%d remainder: ",remainder.size()); printRVSet(f,remainder);  
  fprintf(f,"%d hRemainder: ",hRemainder.size()); printRVSet(f,hRemainder);  
}




////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        InferenceSeparatorClique support
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////



/*-
 *-----------------------------------------------------------------------
 * InferenceSeparatorClique::InferenceSeparatorClique()
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
InferenceSeparatorClique::InferenceSeparatorClique(SeparatorClique& from_clique,
						   vector <RV*>& newRvs,
						   map < RVInfo::rvParent, unsigned >& ppf,
						   const unsigned int frameDelta)
  : origin(from_clique)
{

  set<RV*>::iterator it;

  // clone over nodes RVs.
  fNodes.resize(origin.nodes.size());
  unsigned i=0;
  for (it = origin.nodes.begin();
       it != origin.nodes.end();
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
    fNodes[i++] = nrv;
  }

  i=0;
  fAccumulatedIntersection.resize(origin.accumulatedIntersection.size());
  for (it = origin.accumulatedIntersection.begin();
       it != origin.accumulatedIntersection.end();
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
    fAccumulatedIntersection[i] = nrv;
  }

  i=0;
  fRemainder.resize(origin.remainder.size());
  for (it = origin.remainder.begin();
       it != origin.remainder.end();
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
    fRemainder[i++] = nrv;
  }

  // Separator accumulated intersection values only store/hash values
  // of hidden (thus necessarily discrete) variables since they are
  // the only thing that change.
  accDiscreteValuePtrs.resize(origin.hAccumulatedIntersection.size());
  for (i=0;i<accDiscreteValuePtrs.size();i++) {
    // get the hidden rv for this location
    RV* rv = origin.hAccumulatedIntersection[i];;
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    

    if ( ppf.find(rvp) == ppf.end() ) {
      coredump("ERROR: can't find assigned rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
    }
    RV* nrv = newRvs[ppf[rvp]];

    // hidden nodes are always discrete (in this version).
    // TODO: add hidden continuous variable
    DiscRV* drv = 
      (DiscRV*)nrv;

    // grab a pointer directly to its value for easy access later.
    accDiscreteValuePtrs[i] = &(drv->val);
  }

  // Separator remainder values only store/hash values of hidden (thus
  // necessarily discrete) variables since they are the only thing
  // that change.
  remDiscreteValuePtrs.resize(origin.hRemainder.size());
  for (i=0;i<remDiscreteValuePtrs.size();i++) {
    // get the hidden rv for this location
    RV* rv = origin.hRemainder[i];
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    

    if ( ppf.find(rvp) == ppf.end() ) {
      coredump("ERROR: can't find assigned rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
    }
    RV* nrv = newRvs[ppf[rvp]];

    // hidden nodes are always discrete (in this version).
    DiscRV* drv = 
      (DiscRV*)nrv;

    // grab a pointer directly to its value for easy access later.
    remDiscreteValuePtrs[i] = &(drv->val);
  }

  // allocate at one value for now.
  if (origin.hAccumulatedIntersection.size() == 0) {
    // in this case, we'll only need one and never more.
    separatorValues.resize(1);
    // there will always be one used value here.
    numSeparatorValuesUsed = 1;
    if (origin.hRemainder.size() > 0) {
      // So we have no accumulated intersection, and a remainder which
      // means we are in a good position to predict the size of the
      // (necessarily single) remainder vectors from the previous
      // times we used it. Therefore, we do just that, but only in
      // this case.
      separatorValues.ptr[0].remValues.resize(origin.remainderValueSpaceManager.currentSize());
      new (&separatorValues.ptr[0].iRemHashMap)VHashMapUnsignedUnsignedKeyUpdatable
	(origin.remPacker.packedLen(),origin.remainderValueSpaceManager.currentSize());
      // new (&separatorValues.ptr[0].iRemHashMap)VHashMapUnsignedUnsignedKeyUpdatable(origin.remPacker.packedLen(),2);
    } else {
      // The separator consists of all observed nodes. 
      // Search in file for key string
      // "ALLOCATE_REMVALUES_ALL_OBSERVED'
      // to find where the nec. single entry is allocated.
    }
  } else {
    // start with something a bit larger
    // TODO: optimize this.
    const unsigned starting_size = origin.separatorValueSpaceManager.currentSize(); // 3,2000;
    separatorValues.resize(starting_size);
    if (origin.hRemainder.size() > 0) {
      for (unsigned i=0;i<starting_size;i++) {
	// need to re-construct individual hash tables.
	new (&separatorValues.ptr[i].iRemHashMap)VHashMapUnsignedUnsignedKeyUpdatable
	  (origin.remPacker.packedLen(),2);
	// TODO: while we potentially could preallocate default size
	// of separatorValues.ptr[i].remValues.resize(default); here,
	// we don't really know what it should be. Since there are
	// multiple remainders here, and each might be drastically
	// different in size (some even being zero length), it is not
	// a good idea to allocate anything at all and we let it
	// lazily be sized as needed. The TODO is to come up with a
	// better scheme (e.g., keep a counter array to keep track of
	// the number at each allocation level, and then allocate this
	// to the min size of the previous time that that a non-zero
	// number of cases.
      }
    } else {
      // things such as array separatorValues.ptr[i].remValues will be sized as needed later.
      // Search for tag 'ALLOCATE_REMVALUES_OPTION' in this file for where it is allocated.
    }
    // need to re-construct the hash table.
    new (&iAccHashMap) VHashMapUnsignedUnsignedKeyUpdatable
      (origin.accPacker.packedLen(),starting_size); // 2
    numSeparatorValuesUsed = 0;
  }

}



/*-
 *-----------------------------------------------------------------------
 * InferenceSeparatorClique::ceSeparatorPrune()
 *
 *    Collect Evidence, Separator Prune: This routine will prune away
 *    part of a previously instantiated separator based on the current
 *    separator beam width.
 *
 * Preconditions:
 *   1) separator table must be created, meaning that either:
 *
 *        InferenceMaxClique::ceSendToOutgoingSeparator()
 *      must have been called sending a message (projection downto) this separator.
 *
 * Postconditions:
 *    Separator table has been pruned, and memory for it has been re-allocated to
 *    fit the smaller size. All hash tables adjusted.
 *
 * Side Effects:
 *    potential memory allocations and hash table adjustments.
 *
 * Results:
 *     nothing
 *
 *-----------------------------------------------------------------------
 */

void 
InferenceSeparatorClique::ceSeparatorPrune()
{
  // return immediately if separator beam pruning is turned off.
  if (origin.separatorBeam == (-LZERO))
    return;

  // compute max and current state space.
  logpr maxCEsepValue;
  unsigned originalTotalStateSpace = 0;
  for (unsigned asv=0;asv<numSeparatorValuesUsed;asv++) {
    originalTotalStateSpace += separatorValues.ptr[asv].numRemValuesUsed;
    for (unsigned rsv=0;rsv<separatorValues.ptr[asv].numRemValuesUsed;rsv++) {
      if (separatorValues.ptr[asv].remValues.ptr[rsv].p > maxCEsepValue)
	maxCEsepValue = separatorValues.ptr[asv].remValues.ptr[rsv].p;
    }
  }

  // Check for all observed case, or a case where there is only one
  // entry, in which case we never do anything.
  if (originalTotalStateSpace == 1)
    return;

  // create an ininitialized variable
  logpr beamThreshold((void*)0);
  // break into the logp to avoid unnecessary zero checking.
  beamThreshold.valref() = maxCEsepValue.valref() - origin.separatorBeam;

  // go through and shrink guys less than maximum.
  unsigned newTotalStateSpace = 0;  
  for (unsigned asv=0;asv<numSeparatorValuesUsed;asv++) {
    const unsigned origNumRemValuesUsed = separatorValues.ptr[asv].numRemValuesUsed;
    for (unsigned rsv=0;rsv<separatorValues.ptr[asv].numRemValuesUsed;) {
      if (separatorValues.ptr[asv].remValues.ptr[rsv].p < beamThreshold) {

	if (separatorValues.ptr[asv].numRemValuesUsed > 1) {


	  // We prune away entry for rsv, by swapping it in last
	  // position. Here, however, it is not as easy as with a clique
	  // separator as we have also to deal with the hash
	  // table. Specifically, we need to swap index entries in hash
	  // table as well. Note that we can not remove the hash table
	  // entry for the one that got pruned away without re-hashing
	  // the entire hash table. The reason is that if the entry that
	  // got removed was a collision for another entry that is in
	  // the table, then removing the collision will make the other
	  // entry inaccessible.
	  // TODO: test if it is better to just prune here and just
	  // rehash everything.

	  // the index of the entry being swapped with the
	  // one that is being pruned.
	  const unsigned swap_index = separatorValues.ptr[asv].numRemValuesUsed-1;

	  // First, get pointers to hash table index values for the two
	  // entries corresponding to the one we are prunning
	  // and the one ware swapping it with.
	  unsigned* prune_index_p;
	  unsigned* swap_index_p;

	  // the keys for the two entries.
	  unsigned* prune_key_p;
	  unsigned* swap_key_p;

	  // pointers to the ht keys for the two entries.
	  unsigned** ht_prune_key_p;
	  unsigned** ht_swap_key_p;

	  if (origin.remPacker.packedLen() <= ISC_NWWOH_RM) {
	    prune_key_p = &(separatorValues.ptr[asv].remValues.ptr[rsv].val[0]);
	    swap_key_p = &(separatorValues.ptr[asv].remValues.ptr[swap_index].val[0]);
	  } else {
	    prune_key_p = separatorValues.ptr[asv].remValues.ptr[rsv].ptr;
	    swap_key_p = separatorValues.ptr[asv].remValues.ptr[swap_index].ptr;
	  }
	
	  prune_index_p =  separatorValues.ptr[asv].iRemHashMap.find(prune_key_p,ht_prune_key_p);
	  // it must exist
	  assert ( prune_index_p != NULL );
	  swap_index_p =  separatorValues.ptr[asv].iRemHashMap.find(swap_key_p,ht_swap_key_p);
	  // it must exist
	  assert ( swap_index_p != NULL );

	  // swap the entries
	  swap(separatorValues.ptr[asv].remValues.ptr[rsv],
	       separatorValues.ptr[asv].remValues.ptr[swap_index]);
	  // and swap the hash table pointers
	  swap((*prune_index_p),(*swap_index_p));
	  // and swap the hash table keys if they are pointers to the arrays which
	  // just got swapped.
	  if (origin.remPacker.packedLen() <= ISC_NWWOH_RM) {
	    // printf("foobarbaz");
	    swap((*ht_prune_key_p),(*ht_swap_key_p));
	  }

	  // decrease values
	  separatorValues.ptr[asv].numRemValuesUsed--;
	} else {
	  // then separatorValues.ptr[asv].numRemValuesUsed == 1. This
	  // will happen under two conditiosn.
	  //  1) There is no remainder, meaning (sep.origin.hRemainder.size() == 0), 
	  //     and the entire separator is
	  //     in the accumulated intersection with other separators.
	  //     In this case there is no hash table at all.
	  //  2) All other entries in the remainder for this particular
	  //     accumulated intersection slot have been pruned away. In this
	  //     case, the hash table does exist, but it'll still be deleted
	  //     upon calling the destructor.

	  // so, what we do is just set the the number of values to zero.
	  // Other code will need to thus check for empty separator acc. intersection
	  // values.
	  separatorValues.ptr[asv].numRemValuesUsed = 0;

	}

      } else {
	rsv++;
      }
    }
    newTotalStateSpace += separatorValues.ptr[asv].numRemValuesUsed;
    if (separatorValues.ptr[asv].numRemValuesUsed < origNumRemValuesUsed) {
      if (separatorValues.ptr[asv].numRemValuesUsed == 0 ) {
	// should/could remove accumulator entry here as well. 
      }
      // - re-allocate memory & adjust hash table.
      // - separatorValues.ptr[asv].remValues
      // - possibly re-hash hash tables if necessary.
    }
  }

  infoMsg(IM::Huge,"Separator beam pruning, Max cv = %f, thres = %f. Original sep state space = %d, new sep state space = %d\n",
	  maxCEsepValue.valref(),
	  beamThreshold.valref(),
	  originalTotalStateSpace,newTotalStateSpace);

}



////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        CliqueValueHolder support
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


CliqueValueHolder::CliqueValueHolder(unsigned _cliqueValueSize,
				     unsigned _allocationUnitChunkSize,
				     float _growthFactor)
  : cliqueValueSize(_cliqueValueSize),
    growthFactor(_growthFactor),
    allocationUnitChunkSize(_allocationUnitChunkSize)
{
  assert ( allocationUnitChunkSize > 0 );
  prepare();
}


// make sure there is rom for one element.
void
CliqueValueHolder::prepare()
{
  makeEmpty();
  values.resize(1);
  // newSize *MUST* be a multiple of 'cliqueValueSize' or else
  // this code will fail.
  unsigned newSize = cliqueValueSize*allocationUnitChunkSize;
  values[values.size()-1].resize(newSize);
  capacity = newSize;
  numAllocated = 0;
  curAllocationPosition = values[values.size()-1].ptr;
  curAllocationEnd = values[values.size()-1].ptr + newSize;

}


// free up all memory, postcondition: make an invalid object until
// next prepare is called.
void
CliqueValueHolder::makeEmpty()
{
  for (unsigned i=0;i<values.size();i++) {
    values[i].clear();
  }
  values.clear();
  capacity = 0;
  numAllocated = 0;
}



void
CliqueValueHolder::allocateCurCliqueValue()
{

  curAllocationPosition += cliqueValueSize;
  numAllocated++;

  // first to a cheap and fast allocation of new clique value storage
  if (curAllocationPosition != curAllocationEnd) {
    return;
  } 

  // if here, we need to allocate another chunk add a new chunk so we
  // don't need to re-copy all the existing ones already.
  values.resizeAndCopy(values.size()+1);

  // newSize *MUST* be a multiple of 'cliqueValueSize' or else
  // this code will fail.
  // TODO: optimize this re-sizing.
  unsigned newSize = cliqueValueSize*
    unsigned(1+allocationUnitChunkSize*
	     ::pow(growthFactor,values.size()-1));

  values[values.size()-1].resize(newSize);

  capacity += newSize;
  curAllocationPosition = values[values.size()-1].ptr;
  curAllocationEnd = values[values.size()-1].ptr + newSize;

}



VHashMapUnsignedUnsignedKeyUpdatable::
VHashMapUnsignedUnsignedKeyUpdatable(const unsigned arg_vsize,
				     unsigned approximateStartingSize)
  : VHashMapUnsignedUnsigned(arg_vsize,approximateStartingSize)
{
  // do nothing
}
	

