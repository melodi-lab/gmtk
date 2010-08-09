/*-
 * GMTK_MaxClique.cc
 *
 *     maxClique table pruning for Maximal cliques
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2009, < fill in later >
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
 * TODO: turn this into multiple files
 *   mc, mctable CE, mctable DE, mctable prune, csctable
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
#include "GMTK_HidDiscRV.h"
#include "GMTK_ObsDiscRV.h"
#include "GMTK_ContRV.h"
#include "GMTK_GMTemplate.h"
#include "GMTK_GMParms.h"
#include "GMTK_MDCPT.h"
#include "GMTK_MSCPT.h"
#include "GMTK_MTCPT.h"
#include "GMTK_Mixture.h"
#include "GMTK_ObservationMatrix.h"
#include "GMTK_JunctionTree.h"

VCID("$Header$")


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Constants
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


// for sorting an array of CliqueValue descending based on the contained logpr 
struct CliqueValueDescendingProbCompare
{  
  bool operator() (const MaxCliqueTable::CliqueValue& cv1,
		   const MaxCliqueTable::CliqueValue& cv2)
  {
    return (cv1.p > cv2.p);
  }
};


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
 * clique beam width, for use when building clique.  Default
 * value is very large (1.0/0.0 = est. of infty) meaning that we do no
 * beam pruning.
 *
 */
double
MaxClique::cliqueBeamBuildBeam=(-LZERO);

char* MaxClique::cliqueBeamBuildFilter = (char*)"";

bool MaxClique::cliqueBeamContinuationHeuristic = true;

double MaxClique::cliqueBeamBuildExpansionFactor = 1.0;
// could also just have a maximum beam value rather than max expansions.
unsigned MaxClique::cliqueBeamBuildMaxExpansions = 1;


unsigned MaxClique::cliqueBeamClusterPruningNumClusters = 0;
double MaxClique::cliqueBeamClusterBeam = (-LZERO);

#define NO_PRUNING_CLIQUEBEAMCLUSTERPRUNINGMAXSTATES (UINT_MAX)
unsigned MaxClique::cliqueBeamClusterPruningMaxStates = NO_PRUNING_CLIQUEBEAMCLUSTERPRUNINGMAXSTATES;


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
 * Fraction of clique mass to relinquish. Default (0.0) means prune nothing.
 *
 */
double
MaxClique::cliqueBeamMassRelinquishFraction = 0.0;

/*
 * When using cliqueBeamMassRelinquishFraction, this option determins the lower bound of the
 * min min clique state space size.
 *
 */
unsigned 
MaxClique::cliqueBeamMassMinSize = 1;



/*
 * When using mass beam pruning, we use this to further allow additional states through below
 * the one below first one after the mass has been acounted for. 0.0 to turn off.
 *
 */
double
MaxClique::cliqueBeamMassFurtherBeam = 0.0;



/*
 * When using mass beam pruning, we can exponentiate the clique scores 
 * to make them more uniform if need be, so as to hopefully prune more effectively.
 *
 */
double
MaxClique::cliqueBeamMassExponentiate = 1.0;


/*
 * If 0 <= val <= 1, then fraction of pruned clique table to uniformly re-sample from (without replacement)
 * if > 1, then number of times to re-sample (without replacement) from pruned clique table.
 * Set to 0.0 to turn off.
 */
double
MaxClique::cliqueBeamUniformSampleAmount = 0.0;


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
psp2(FILE*f,const int numSpaceChars,const char c = ' ')
{
  int tmp = numSpaceChars;
  // prints a '|' every 4 chars. 
  if (!tmp) return;
  fprintf(f,"|"); tmp--;
  while (tmp > 4) {
    fprintf(f,"   |");
    tmp -= 4;
  }
  while (tmp--)
    fprintf(f,"%c",c);
}

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
			     CLIQUE_VALUE_SPACE_MANAGER_GROWTH_RATE,   // growth rate
			     1,     // growth addition
			     CLIQUE_VALUE_SPACE_MANAGER_DECAY_RATE)    // decay rate 
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
 * TODO: remove this routine since it is not being used.
 *
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
 * TODO: update this routine once new inference is in place.
 *
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

  return JunctionTree::jtWeightSparseNodeSepScale*weight_sep_sparse + 
    JunctionTree::jtWeightDenseNodeSepScale*weight_sep_dense 
    + weight_remainder;
#endif
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
    cliqueValueHashSet.clear(CLIQUE_VALUE_HOLDER_STARTING_SIZE);
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




/*-
 *-----------------------------------------------------------------------
 * MaxClique::reportScoreStats()
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
MaxClique::
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
    new (&valueHolder) CliqueValueHolder(packer.packedLen(),
					 CLIQUE_VALUE_HOLDER_STARTING_SIZE, // set to 1 to test.
					 CLIQUE_VALUE_HOLDER_GROWTH_RATE); // 1.25
    // set up common clique hash tables 
    // TODO: add appropriate default staring hash sizes.
    // new (&cliqueValueHashSet) vhash_set< unsigned > (packer.packedLen(),2);
    new (&cliqueValueHashSet) vhash_set< unsigned > (packer.packedLen(),CLIQUE_VALUE_HOLDER_STARTING_SIZE);
#ifdef USE_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL
    temporaryCliqueValuePool.resize(CLIQUE_VALUE_HOLDER_STARTING_SIZE*packer.packedLen());
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
  if (cliqueBeamContinuationHeuristic) {
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
			  set<RV*>* lp_nodes,set<RV*>* rp_nodes)
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
  fprintf(f,"%ld Nodes: ",(unsigned long)nodes.size()); printRVSet(f,nodes);
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
	if ((JunctionTree::useVESeparators & JunctionTree::VESEP_PCG)  && 
	    (assignedProbNodes.find(p) !=  assignedProbNodes.end())) {

	  float logProdCard = 
	    log10((double)RV2DRV(p->allParents[0])->cardinality);
	  for (unsigned i=1;i<p->allParents.size();i++) {
	    if (!p->allParents[i]->discreteObservedImmediate())
	      logProdCard += log10((double)RV2DRV(p->allParents[i])->cardinality);
	  }

	  if (logProdCard > SeparatorClique::veSeparatorLogProdCardLimit)
	    continue;

	  infoMsg(Max,"Found VE sep of type PCG (currently not using it).\n");
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
      } else if ((JunctionTree::useVESeparators & JunctionTree::VESEP_PC) 
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

	infoMsg(Max,"Found VE sep of type PC. iterable = %d\n",rv->iterable());

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

void
MaxClique::numParentsSatisfyingChild(unsigned& num,unsigned par,vector <RV*> & parents, RV* child)
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



////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        MaxCliqueTable support
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////




/*-
 *-----------------------------------------------------------------------
 * MaxCliqueTable::MaxCliqueTable()
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
MaxCliqueTable::MaxCliqueTable(MaxClique& origin)
{
  init(origin);
}

// get clique ready for use.
void MaxCliqueTable::init(MaxClique& origin)
{
  numCliqueValuesUsed = 0;
#ifdef TRACK_NUM_CLIQUE_VALS_SHARED
  numCliqueValuesShared = 0;
#endif

  // TODO: optimize this and make depend on if clique is all hidden, has observed, etc.
  // NOTE: This must be set to something greater than 0.
  // cliqueValues.resize(3); // 10000
  cliqueValues.resize(origin.cliqueValueSpaceManager.currentSize());

}


/*-
 *-----------------------------------------------------------------------
 * MaxCliqueTable::SharedLocalStructure::MaxCliqueTable() constructor.
 *    INitialize the shared part of a MaxCliqueTable
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
MaxCliqueTable::SharedLocalStructure::
SharedLocalStructure(MaxClique& _origin,
		     vector <RV*>& newRvs,
		     map < RVInfo::rvParent, unsigned >& ppf,
		     const unsigned int frameDelta)
{


  origin = &_origin;

  set<RV*>::iterator it;

  // clone over nodes RVs.
  fNodes.resize(origin->nodes.size());
  unsigned i=0;
  for (it = origin->nodes.begin();
       it != origin->nodes.end();
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
  fSortedAssignedNodes.resize(origin->sortedAssignedNodes.size());
  for (i=0;i<origin->sortedAssignedNodes.size();i++) {
    RV* rv = origin->sortedAssignedNodes[i];
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
  fUnassignedIteratedNodes.resize(origin->unassignedIteratedNodes.size());
  for (it = origin->unassignedIteratedNodes.begin();
       it != origin->unassignedIteratedNodes.end();
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

  // Clique values only store/hash values of hidden (thus necessarily
  // discrete) variables since they are the only thing that change.
  discreteValuePtrs.resize(origin->hashableNodes.size());
  for (i=0;i<discreteValuePtrs.size();i++) {
    // get the unrolled rv for this hidden node
    RV* rv = origin->hashableNodes[i];
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

  fDeterminableNodes.resize(origin->determinableNodes.size());
  for (i=0;i<fDeterminableNodes.size();i++) {
    // get the unrolled rv for this hidden node
    RV* rv = origin->determinableNodes[i];
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
    fDeterminableNodes[i] = RV2DRV(nrv);
  }

}


set <RV*> 
MaxCliqueTable::SharedLocalStructure::returnRVsAndTheirObservedParentsAsSet()
{
  set<RV*> rc;
  for (unsigned i=0;i<fNodes.size();i++) {
    rc.insert(fNodes[i]);
    set <RV*> tmp = fNodes[i]->observedParents();
    // add in observed parents
    std::copy(tmp.begin(), tmp.end(),
	      inserter(rc,rc.end()));
  }
  return rc;
}

set <RV*> 
MaxCliqueTable::SharedLocalStructure::returnRVsAsSet()
{
  set<RV*> rc;
  for (unsigned i=0;i<fNodes.size();i++) {
    rc.insert(fNodes[i]);
  }
  return rc;
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
MaxCliqueTable::
ceGatherFromIncommingSeparators(MaxCliqueTable::SharedLocalStructure& sharedStructure,
				ConditionalSeparatorTable* separatorTableArray,
				ConditionalSeparatorTable::SharedLocalStructure* sepSharedStructureArray)
{

  // syntactic convenience variables.
  MaxClique& origin = *(sharedStructure.origin);

  // we should never try more than 1x in this case.
  assert ( origin.cliqueBeamBuildBeam != (-LZERO) || origin.cliqueBeamBuildMaxExpansions == 1 );

  unsigned cliqueExpansionTry = 0;

  // Max collect-evidence probability for this clique. Used for beam
  // pruning.
  logpr maxCEValue;

  // An estimate of the beam threshold (i.e., est ~ maxCeValue -
  // threshold) for the current clique based on predictions from a
  // history of previous clique's maxCE values. If clique values fall
  // below this estimated threshold, they are pruned.
  logpr cliqueBeamThresholdEstimate;

  while (cliqueExpansionTry < origin.cliqueBeamBuildMaxExpansions) {

    traceIndent=-1; 
    // this is like the sub-main() for collect evidence.
    if (origin.cliqueBeamBuildBeam != (-LZERO)
	&& origin.maxCEValuePredictor.ptr() != NULL
	&& origin.maxCEValuePredictor->readyToMakePrediction()) {

      double currentCliqueBeamBuildBeam = 
	origin.cliqueBeamBuildBeam * (::pow(origin.cliqueBeamBuildExpansionFactor,cliqueExpansionTry));
      double maxCEValPrediction = origin.maxCEValuePredictor->makePrediction();
      double fixedPrediction = 2*origin.prevMaxCEValue.valref() - origin.prevPrevMaxCEValue.valref();

      cliqueBeamThresholdEstimate.valref() = maxCEValPrediction - 
	currentCliqueBeamBuildBeam;
      if (cliqueBeamThresholdEstimate.essentially_zero())
	cliqueBeamThresholdEstimate.set_to_almost_zero();
      // report the various values, 'p' = previous or predicted.

      double denom = origin.prevMaxCEValue.valref();
      // report absolute error if we can't compute relative error.
      if (denom == 0.0) denom = 1.0;
      infoMsg(IM::Med,"Partial clique beam pruning, ppmax= %f, pmax= %f, ppred= %f, pfpred= %f, p_rel_%%err = %f, p_frel_%%err= %f, pred= %f, pthres= %f.\n",
	      origin.prevPrevMaxCEValue.valref(),
	      origin.prevMaxCEValue.valref(),
	      origin.prevMaxCEValPrediction,
	      origin.prevFixedPrediction,
	      100*::fabs(origin.prevMaxCEValPrediction - origin.prevMaxCEValue.valref())/
	      ::fabs(denom),
	      100*::fabs(origin.prevFixedPrediction - origin.prevMaxCEValue.valref())/
	      ::fabs(denom),
	      maxCEValPrediction,
	      cliqueBeamThresholdEstimate.valref());
      if (cliqueExpansionTry == 0) {
	origin.prevMaxCEValPrediction = maxCEValPrediction;
	origin.prevFixedPrediction = fixedPrediction;
      }
    } else {

      // always prune if we fall below or equal to almost zero.
      cliqueBeamThresholdEstimate.set_to_almost_zero(); 
      origin.prevFixedPrediction = 0;
      origin.prevMaxCEValPrediction = 0;
    }

    maxCEValue.set_to_zero();
    // next, do the actual collect message.
    if (origin.hashableNodes.size() == 0) {
      ceGatherFromIncommingSeparatorsCliqueObserved(sharedStructure,
						    separatorTableArray,
						    sepSharedStructureArray,
						    maxCEValue);
    } else {
      // if we're still here, we do regular separator driven inference.
      logpr p = 1.0;
      if (origin.ceReceiveSeparators.size() == 0) {
	if (origin.unassignedIteratedNodes.size() == 0) {
	  ceIterateAssignedNodes(sharedStructure,
				 cliqueBeamThresholdEstimate,
				 maxCEValue, // max value that is returned
				 0, // initial node number
				 p);
	} else {
	  ceIterateUnassignedIteratedNodes(sharedStructure,
					   cliqueBeamThresholdEstimate,
					   maxCEValue,
					   0,
					   p);
	}
      } else {
	ceIterateSeparators(sharedStructure,
			    separatorTableArray,
			    sepSharedStructureArray,
			    cliqueBeamThresholdEstimate,
			    maxCEValue,
			    0,
			    p);
      }
    }

    if (numCliqueValuesUsed == 0) {
      // if we have a zero clique, print message and possibly continue
      // with expanded clique.
      if (message(IM::Med)) {
	printf("WARNING: ZERO CLIQUE: clique with no entries, try %d out of %d.\n",
	       cliqueExpansionTry+1,origin.cliqueBeamBuildMaxExpansions);
      }
    } else {
      // current pruning level worked.
      break;
    }
    cliqueExpansionTry ++;
  }

  // check if we have a zero clique, and if we do, print message and exit.
  // TODO: rather than exit, pop back to the top and allow continuation and/or
  // beam expansion.
  if (numCliqueValuesUsed == 0)
    error("ERROR: ZERO CLIQUE: clique with no entries. Final probability will be zero.\n");

  // We have some clique entries, so we store new previous max CE
  // values, before any pruning.
  if (origin.cliqueBeamBuildBeam != (-LZERO)) {
    origin.prevPrevMaxCEValue.valref() = origin.prevMaxCEValue.valref();
    origin.prevMaxCEValue.valref() = maxCEValue.valref();
    if (origin.maxCEValuePredictor.ptr() != NULL) {
      origin.maxCEValuePredictor->addNextSampleAndUpdate(maxCEValue.valref());
    }
  }

  // if (numCliqueValuesUsed == 278664)
  // printCliqueEntries(sharedStructure,stdout,NULL,false,false);

  // we prune here right away.
  ceDoAllPruning(origin,maxCEValue);

}


/*-
 *-----------------------------------------------------------------------
 * MaxCliqueTable::ceGatherFromIncommingSeparatorsCliqueObserved()
 *
 *    A version of ceGatherFromIncommingSeparators that is is
 *    specifically for cliques that are all observed which means that
 *    all surrounding separators are also all observed.  observed
 *    clique observedclique
 *
 * Preconditions:
 *      Same as the separator driven case.
 *
 * Postconditions:
 *      Same as the separator driven case.
 *
 * Side Effects:
 *      Same as the separator driven case.
 *
 * Results:
 *     nothing
 *
 *-----------------------------------------------------------------------
 */

void
MaxCliqueTable::
ceGatherFromIncommingSeparatorsCliqueObserved(MaxCliqueTable::SharedLocalStructure& sharedStructure,
					      ConditionalSeparatorTable* separatorTableArray,
					      ConditionalSeparatorTable::SharedLocalStructure* sepSharedStructureArray,
					      logpr& maxCEValue)
{

  MaxClique& origin = *(sharedStructure.origin);

  // just do the assigned nodes as the unassigned observed nodes have no effect here.
  unsigned nodeNumber;

  logpr p = 1.0;
  bool zero = false;
  for (nodeNumber = 0; nodeNumber < sharedStructure.fSortedAssignedNodes.size(); nodeNumber ++ ) {
    RV* rv = sharedStructure.fSortedAssignedNodes[nodeNumber];
    if (origin.dispositionSortedAssignedNodes[nodeNumber] == MaxClique::AN_NOTSEP_PROB_SPARSEDENSE || 
	origin.dispositionSortedAssignedNodes[nodeNumber] == MaxClique::AN_NOTSEP_PROB_SPARSEDENSE) {
      // apply the probabiltiy
      logpr cur_p = rv->probGivenParents();

      if (message(Huge)) {
	psp2(stdout,spi*(traceIndent+1+nodeNumber));
	printf("%d:assigned obs/prob app, Pr[",nodeNumber);
	rv->printNameFrameValue(stdout,false);
	if (message(Mega)) {
	  if (rv->allParents.size() > 0) {
	    printf("|");
	    printRVSetAndValues(stdout,rv->allParents,false);
	  }
	} else {
	  if (rv->allParents.size() > 0) {
	    printf("|parents");
	  }
	}
	printf("]=%f,crClqPr=%f\n",cur_p.val(),p.val());
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

      if (message(Huge)) {
	psp2(stdout,spi*(traceIndent+1+nodeNumber));
	printf("%d:assigned obs/zero rmv, Pr[",nodeNumber);
	rv->printNameFrameValue(stdout,false);
	if (message(Mega)) {
	  if (rv->allParents.size() > 0) {
	    printf("|");
	    printRVSetAndValues(stdout,rv->allParents,false);
	  }
	} else {
	  if (rv->allParents.size() > 0) {
	    printf("|parents");
	  } 
	}
	printf("]=%f,crClqPr=%f\n",cur_p.val(),p.val());
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
      // get a handy reference to the current separator table and its origin
      ConditionalSeparatorTable& sep = 
	separatorTableArray[origin.ceReceiveSeparators[sepNumber]];
      SeparatorClique& sepOrigin = 
	*(sepSharedStructureArray[origin.ceReceiveSeparators[sepNumber]].origin);

      // if the separator is currently not active, we should not be
      // using it for the additional probability here.
      if (sepOrigin.skipMe)
	continue;

      // keep a local variable copy of this around to avoid potential dereferencing.
      ConditionalSeparatorTable::AISeparatorValue * const
	sepSeparatorValuesPtr = sep.separatorValues->ptr; 


      unsigned accIndex = 0;
      // separator consists of all observed values. We just multiply
      // in the one entry the separator, that is guaranteed to be at
      // postiion 0,0. The key thing is that we must not do ANY hash
      // lookup as there are no hash tables in a separator that
      // has no hidden variables.

      ConditionalSeparatorTable::AISeparatorValue& sv
	= sepSeparatorValuesPtr[accIndex];

      // If anyone is almost zero, stop right now. 
      // Where was this allocated?  Search in file for key string
      // "ALLOCATE_REMVALUES_ALL_OBSERVED'
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

  if (message(High)) {
    psp2(stdout,spi*(traceIndent+1+sharedStructure.fSortedAssignedNodes.size()));
    infoMsg(High,"CI:Inserting Observed %d-clique ent #0,pr=%f,sm=%f:",
	    sharedStructure.fNodes.size(),
	    cliqueValues.ptr[0].p.val(),sumProbabilities().val());
    printRVSetAndValues(stdout,sharedStructure.fNodes);
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
MaxCliqueTable::ceIterateSeparatorsRecurse(MaxCliqueTable::SharedLocalStructure& sharedStructure,
					   ConditionalSeparatorTable* separatorTableArray,
					   ConditionalSeparatorTable::SharedLocalStructure* sepSharedStructureArray,
					   logpr cliqueBeamThresholdEstimate,
					   logpr& maxCEValue,
					   const unsigned sepNumber,
					   const logpr p)
{

  if (p <= cliqueBeamThresholdEstimate)
    return;

  // syntactic simplicity and cached variables ...
  MaxClique& origin = *(sharedStructure.origin);
  ConditionalSeparatorTable& sep = 
    separatorTableArray[origin.ceReceiveSeparators[sepNumber]];
  ConditionalSeparatorTable::SharedLocalStructure& sepSharedStructure = 
    sepSharedStructureArray[origin.ceReceiveSeparators[sepNumber]];

  DiscRVType** accDiscreteValuePtrs = 
    sepSharedStructure.accDiscreteValuePtrs.ptr;
  DiscRVType** remDiscreteValuePtrs = 
    sepSharedStructure.remDiscreteValuePtrs.ptr;

  SeparatorClique& sepOrigin = 
    *(sepSharedStructure.origin);
  ConditionalSeparatorTable::AISeparatorValue * const
    sepSeparatorValuesPtr = sep.separatorValues->ptr; 


  if (message(High+5)) {
    traceIndent++;
    psp2(stdout,spi*traceIndent);    
    infoMsg(High+5,"S%d:Starting separator iter,partSepNo=%d,p=%f,nodes:",
	    sepNumber,origin.ceReceiveSeparators[sepNumber],p.val());
    printRVSet(stdout,sepSharedStructure.fNodes);
  }

  if (sepOrigin.skipMe) {
    // then we completely skip this separator passing onto the next
    // one.
    ceIterateSeparators(sharedStructure,
			separatorTableArray,
			sepSharedStructureArray,
			cliqueBeamThresholdEstimate,
			maxCEValue,
			sepNumber+1,
			p);
    goto ceIterateSeparatorsFinished;
  }


  unsigned sepValueNumber;  
  if (sepOrigin.hAccumulatedIntersection.size() > 0) {
    // look up existing intersected values to see if we have a match
    // and only proceed if we do.

    unsigned *key = &CliqueBuffer::packedVal[0];
    sepOrigin.accPacker.pack((unsigned**)accDiscreteValuePtrs,
			     (unsigned*)key);
    unsigned* indexp = sep.iAccHashMap->find(key);
    if (indexp == NULL) {
      // Then not found in this separator, so it must have (pruned) zero
      // probability. We continue with the next value of the previous
      // separator.
      if (message(Huge)) {
	psp2(stdout,spi*traceIndent);
	infoMsg(Huge,"S%d:Separator iter accumulated intersection prune\n",
		sepNumber);
	// TODO: @@@ figure out why we can't do: 
	// printRVSet(stdout,sep.fAccumulatedIntersection);
      }
      goto ceIterateSeparatorsFinished;
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
    assert ( sep.separatorValues->size() == 1);
    sepValueNumber = 0;
  }

  // Iterate through remainder of separator. 
  // NOTE: we could do some online pruning here as well, but instead
  // we do it in a special separator prune routine, called ceSeparatorPrune().
  if (sepOrigin.hRemainder.size() == 0) {
    // Only one remainder entry (in position 0) and also no need to
    // unpack since all has been covered by accumulated intersection
    // set above in a previous separator. Just continue on with single
    // value. 
    // - 
    // Note that this case also includes the case when the entire
    // separator is observed, i.e., when we'll get an accum inter size
    // of 0 and a hreminder size of 0 -- nothing to unpack here either
    // (no hash tables even exist), so we just continue along.

    if (message(Huge)) {
      psp2(stdout,spi*traceIndent);
      infoMsg(Huge,"S%d:Separator iter no-unpack %d,%d,partSepNo=%d,p=%f,sp=%f,nodes:",
	      sepNumber,
	      sepSeparatorValuesPtr[sepValueNumber].remValues.size(),
	      sepSeparatorValuesPtr[sepValueNumber].numRemValuesUsed,
	      origin.ceReceiveSeparators[sepNumber],
	      p.val(),
	      sepSeparatorValuesPtr[sepValueNumber].remValues.ptr[0].p.val());
      printRVSetAndValues(stdout,sepSharedStructure.fNodes);
    }

    // We should have either one or zero entries. The only way zero
    // entries could arrise is if we have done severe separator pruning.
    assert ( (sepSeparatorValuesPtr[sepValueNumber].remValues.size() & ~0x1)
	     == 0x0 );

    // assert ( sepSeparatorValuesPtr[sepValueNumber].remValues.size() == 1 );
    if (sepSeparatorValuesPtr[sepValueNumber].numRemValuesUsed == 1) {
      // Continue down with new probability value.
      // Search for tag 'ALLOCATE_REMVALUES_OPTION' in this file for
      // more info why remValues.ptr[0] exists.
      // Note: could do more separator pruning here.
      ceIterateSeparators(sharedStructure,
			  separatorTableArray,
			  sepSharedStructureArray,
			  cliqueBeamThresholdEstimate,
			  maxCEValue,
			  sepNumber+1,
			  p*
			  sepSeparatorValuesPtr[sepValueNumber].remValues.ptr[0].p);
    }
  } else {

    // TODO: this assertion should be redundant (check above)
    assert ( sepOrigin.remPacker.packedLen() > 0 );

    // TODO: perhaps special case for VE seps, since all probs are == 1, so no need to multiply.

    if (sepOrigin.remPacker.packedLen() <= ISC_NWWOH_RM) {
      for (unsigned i=0;i< sepSeparatorValuesPtr[sepValueNumber].numRemValuesUsed; i++) {

	// TODO: optimize this, pre-compute base array outside of loop.
	sepOrigin.remPacker.unpack(
				   (unsigned*)&(sepSeparatorValuesPtr[sepValueNumber].remValues.ptr[i].val[0]),
				   (unsigned**)remDiscreteValuePtrs);


	if (message(Huge)) {
	  psp2(stdout,spi*traceIndent);
	  infoMsg(Huge,"S%d:Separator iter %d of %d,partSepNo=%d,p=%f,sp=%f,nodes:",
		  sepNumber,
		  i,sepSeparatorValuesPtr[sepValueNumber].numRemValuesUsed,
		  origin.ceReceiveSeparators[sepNumber],
		  p.val(),
		  sepSeparatorValuesPtr[sepValueNumber].remValues.ptr[i].p.val());
	  printRVSetAndValues(stdout,sepSharedStructure.fNodes);
	}

	// continue down with new probability value.
	ceIterateSeparators(sharedStructure,
			    separatorTableArray,
			    sepSharedStructureArray,			    
			    cliqueBeamThresholdEstimate,
			    maxCEValue,
			    sepNumber+1,
			    p*sepSeparatorValuesPtr[sepValueNumber].remValues.ptr[i].p);
      }
    } else {
      for (unsigned i=0;i< sepSeparatorValuesPtr[sepValueNumber].numRemValuesUsed; i++) {

	// TODO: optimize this, pre-compute base array outside of loop.
	sepOrigin.remPacker.unpack(
				   (unsigned*)sepSeparatorValuesPtr[sepValueNumber].remValues.ptr[i].ptr,
				   (unsigned**)remDiscreteValuePtrs);

	if (message(Huge)) {
	  psp2(stdout,spi*traceIndent);
	  infoMsg(Huge,"S%d:Separator iter %d of %d,partSepNo=%d,p=%f,sp=%f,nodes:",
		  sepNumber,
		  i,sepSeparatorValuesPtr[sepValueNumber].numRemValuesUsed,
		  origin.ceReceiveSeparators[sepNumber],
		  p.val(),
		  sepSeparatorValuesPtr[sepValueNumber].remValues.ptr[i].p.val());
	  printRVSetAndValues(stdout,sepSharedStructure.fNodes);
	}

	// continue down with new probability value.
	ceIterateSeparators(sharedStructure,
			    separatorTableArray,
			    sepSharedStructureArray,
			    cliqueBeamThresholdEstimate,
			    maxCEValue,
			    sepNumber+1,
			    p*sepSeparatorValuesPtr[sepValueNumber].remValues.ptr[i].p);
      }
    }

  }

 ceIterateSeparatorsFinished:
  if (message(High+5))
    traceIndent--;
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
 * TODO:
 *     should not force iterating unassigned at the top of the recurssion since it might 
 *     no be relevant or useful for the clique expansion (e.g., switching parents might
 *     mean that this variable isn't needed for most cases).
 *
 *-----------------------------------------------------------------------
 */

void
MaxCliqueTable::ceIterateUnassignedIteratedNodesRecurse(MaxCliqueTable::SharedLocalStructure& sharedStructure,
							logpr cliqueBeamThresholdEstimate,
							logpr& maxCEValue,
							const unsigned nodeNumber,
							const logpr p)
{

  RV* rv = sharedStructure.fUnassignedIteratedNodes[nodeNumber];
  // TODO: update comments here to match others.
  if (message(High+5)) {
    traceIndent++;
    psp2(stdout,spi*traceIndent);
    infoMsg(High+5,"U%d:Starting Unassigned iteration of rv %s(%d),p=%f\n",
	    nodeNumber,
	    rv->name().c_str(),rv->frame(),p.val());
  }
  if (rv->hidden()) {
    // only discrete RVs can be hidden for now.
    DiscRV* drv = (DiscRV*)rv;
    // do the loop right here
    drv->val = 0;
    do {
      if (message(Huge)) {
	psp2(stdout,spi*traceIndent);
	infoMsg(Huge,"U%d:Unassigned iter of rv %s(%d)=%d,p=%f\n",
		nodeNumber,
		rv->name().c_str(),rv->frame(),drv->val,p.val());
      }
      // continue on, effectively multiplying p by unity.
      ceIterateUnassignedIteratedNodes(sharedStructure,
				       cliqueBeamThresholdEstimate,
				       maxCEValue,
				       nodeNumber+1,p);
    } while (++drv->val < drv->cardinality);
  } else {
    // TODO: Perhaps unassignedIteratedNodes should contain
    // no observed nodes at all since we are not updating 
    // probability here anyway.

    if (message(Huge)) {
      psp2(stdout,spi*traceIndent);
      // observed, either discrete or continuous
      if (rv->discrete()) {
	infoMsg(Huge,"U%d:Unassigned pass through observed rv %s(%d)=%d,p=%f\n",
		nodeNumber,
		rv->name().c_str(),rv->frame(),RV2DRV(rv)->val,p.val());
      } else {
	// nothing to do since we get continuous observed value
	// indirectly
	infoMsg(Huge,"U%d:Unassigned pass through of observed rv %s(%d)=C,p=%f\n",
		nodeNumber,
		rv->name().c_str(),rv->frame(),p.val());
      }
    }
    // continue on, effectively multiplying p by unity.
    ceIterateUnassignedIteratedNodes(sharedStructure,
				     cliqueBeamThresholdEstimate,
				     maxCEValue,
				     nodeNumber+1,p);
  }
  if (message(High+5))
    traceIndent--;
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
MaxCliqueTable::ceIterateAssignedNodesRecurse(MaxCliqueTable::SharedLocalStructure& sharedStructure,
					      logpr cliqueBeamThresholdEstimate,
					      logpr& maxCEValue,
					      const unsigned nodeNumber,
					      const logpr p)
{

  // syntactic convenience variables.
  MaxClique& origin = *(sharedStructure.origin);

  // First, do (potentially) partial clique beam pruning here.
  // 
  // Note, we may wish, if (nodeNumber == fSortedAssignedNodes.size())
  // to just go ahead and add the clique entry since it is already
  // here, but we still prune since the very last variable that was
  // assigned might have given us a value that falls below threshold.
  if (p*origin.sortedAssignedContinuationScores[nodeNumber] <= cliqueBeamThresholdEstimate)
    return;

  // No pruning, so we go ahead with the continued expansion.

  if (nodeNumber == sharedStructure.fSortedAssignedNodes.size()) {
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
			 (unsigned**)sharedStructure.discreteValuePtrs.ptr,
			 (unsigned*)&(cliqueValues.ptr[numCliqueValuesUsed].val[0]));
    } else {

#ifdef USE_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL      
      const unsigned long lindex = numCliqueValuesUsed*origin.packer.packedLen();
      if (lindex >= origin.temporaryCliqueValuePool.size()) {
	// use aggressive growth factor for now to avoid expensive copies.
	origin.temporaryCliqueValuePool.resizeAndCopy(
						      origin.packer.packedLen()*
						      int(1.5+(double)origin.temporaryCliqueValuePool.size()*TEMPORARY_LOCAL_CLIQUE_VALUE_POOL_GROWTH_RATE));
      }
      unsigned *pcv = 
	&origin.temporaryCliqueValuePool.ptr[lindex];
      origin.packer.pack((unsigned**)sharedStructure.discreteValuePtrs.ptr,(unsigned*)pcv);
      // store integer value of the location.
      cliqueValues.ptr[numCliqueValuesUsed].ival = lindex;
#else

      // Deal with the hash table to re-use clique values.
      // First, grab pointer to storge where next clique value would
      // be stored if it ends up being used.
      unsigned *pcv = origin.valueHolder.curCliqueValuePtr();
      // Next, pack the clique values into this position.
      origin.packer.pack((unsigned**)sharedStructure.discreteValuePtrs.ptr,(unsigned*)pcv);
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
#endif


    }
    // save the probability
    cliqueValues.ptr[numCliqueValuesUsed].p = p;
    numCliqueValuesUsed++;

    if (message(High)) {
      psp2(stdout,spi*(traceIndent+1));
      infoMsg(High,"CI:Inserting %d-clique ent #%d,pr=%f,sm=%f:",
	      sharedStructure.fNodes.size(),
	      (numCliqueValuesUsed-1),
	      cliqueValues.ptr[numCliqueValuesUsed-1].p.val(),sumProbabilities().val());
      printRVSetAndValues(stdout,sharedStructure.fNodes);
    }
    return;
  }
  RV* rv = sharedStructure.fSortedAssignedNodes[nodeNumber];
  // do the loop right here

  if (message(High+5)) {
    traceIndent++;
    psp2(stdout,spi*traceIndent);
    infoMsg(High+5,"A%d:Starting assigned iteration of rv %s(%d),crClqPr=%f\n",
	    nodeNumber,
	    rv->name().c_str(),rv->frame(),p.val());
  }

  switch (origin.dispositionSortedAssignedNodes[nodeNumber]) {
  case MaxClique::AN_CPT_ITERATION_COMPUTE_AND_APPLY_PROB: 
    {
      logpr cur_p;
      rv->begin(cur_p);
      do {
	if (message(Huge)) {
	  psp2(stdout,spi*traceIndent);
	  printf("A%d:assigned iter/prob app, Pr[",nodeNumber);
	  rv->printNameFrameValue(stdout,false);
	  if (message(Mega)) {
	    if (rv->allParents.size() > 0) {
	      printf("|");
	      printRVSetAndValues(stdout,rv->allParents,false);
	    }
	  } else {
	    if (rv->allParents.size() > 0) {
	      printf("|parents");
	    }
	  }
	  printf("]=%f,crClqPr=%f\n",cur_p.val(),p.val());
	}
	// if at any step, we get zero, then back out.
	if (!cur_p.essentially_zero()) {
	  // Continue, updating probability by cur_p, contributing this
	  // probability to the clique potential.
	  ceIterateAssignedNodesRecurse(sharedStructure,
					cliqueBeamThresholdEstimate,maxCEValue,
					nodeNumber+1,p*cur_p);
	}
	// TODO: backjump
      } while (rv->next(cur_p));
    }
    break;

  case MaxClique::AN_CPT_ITERATION_COMPUTE_PROB_REMOVE_ZEROS:
    {
      logpr cur_p;
      rv->begin(cur_p);
      do {
	// At each step, we compute probability
	if (message(Huge)) {
	  psp2(stdout,spi*traceIndent);
	  printf("A%d:assigned iter/zero rmv, Pr[",nodeNumber);
	  rv->printNameFrameValue(stdout,false);
	  if (message(Mega)) {
	    if (rv->allParents.size() > 0) {
	      printf("|");
	      printRVSetAndValues(stdout,rv->allParents,false);
	    }
	  } else {
	    if (rv->allParents.size() > 0) {
	      printf("|parents");
	    } 
	  }
	  printf("]=%f,crClqPr=%f\n",cur_p.val(),p.val());
	}
	// if at any step, we get zero, then back out.
	if (!cur_p.essentially_zero()) {
	  // Continue, do not update probability!!
	  ceIterateAssignedNodesRecurse(sharedStructure,
					cliqueBeamThresholdEstimate,maxCEValue,
					nodeNumber+1,p);
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
	if (message(Huge)) {
	  psp2(stdout,spi*traceIndent);
	  printf("A%d:assigned card iter, Pr[",nodeNumber);
	  rv->printNameFrameValue(stdout,false);
	  if (message(Mega)) {
	    if (rv->allParents.size() > 0) {
	      printf("|");
	      printRVSetAndValues(stdout,rv->allParents,false);
	    }
	  } else {
	    if (rv->allParents.size() > 0) {
	      printf("|parents");
	    }
	  }
	  printf("]=???,crClqPr=%f\n",p.val());
	}
	// Continue, do not update probability!!
	ceIterateAssignedNodesRecurse(sharedStructure,
				      cliqueBeamThresholdEstimate,maxCEValue,
				      nodeNumber+1,p);
      } while (++drv->val < drv->cardinality);
    }
    break;

  case MaxClique::AN_COMPUTE_AND_APPLY_PROB:
    {
      // TODO: Make more efficient version of this, based on the type of
      // RV.
      logpr cur_p = rv->probGivenParents();
      // if at any step, we get zero, then back out.
      if (message(Huge)) {
	psp2(stdout,spi*traceIndent);
	printf("A%d:assigned compute appl prob, Pr[",nodeNumber);
	rv->printNameFrameValue(stdout,false);
	if (message(Mega)) {
	  if (rv->allParents.size() > 0) {
	    printf("|");
	    printRVSetAndValues(stdout,rv->allParents,false);
	  }
	} else {
	  if (rv->allParents.size() > 0) {
	    printf("|parents");
	  }
	}
	printf("]=%f,crClqPr=%f\n",cur_p.val(),p.val());
      }
      if (!cur_p.essentially_zero()) {
	// Continue, updating probability by cur_p.
	ceIterateAssignedNodesRecurse(sharedStructure,
				      cliqueBeamThresholdEstimate,maxCEValue,				      
				      nodeNumber+1,p*cur_p);
      }
    }
    break;

  case MaxClique::AN_CONTINUE:
    if (message(Huge)) {
      psp2(stdout,spi*traceIndent);
      printf("A%d:sep cont, non prob, Pr[",nodeNumber);
      rv->printNameFrameValue(stdout,false);
      if (message(Mega)) {
	if (rv->allParents.size() > 0) {
	  printf("|");
	  printRVSetAndValues(stdout,rv->allParents,false);
	}
      } else {
	if (rv->allParents.size() > 0) {
	  printf("|parents");
	}
      }
      printf("]=???,crClqPr=%f\n",p.val());
    }
    ceIterateAssignedNodesRecurse(sharedStructure,
				  cliqueBeamThresholdEstimate,maxCEValue,
				  nodeNumber+1,p);
    break;

  case MaxClique::AN_CONTINUE_COMPUTE_PROB_REMOVE_ZEROS:
    {
      // TODO: Make more efficient version of this, based on the type of
      // RV.
      logpr cur_p = rv->probGivenParents();
      if (message(Huge)) {
	psp2(stdout,spi*traceIndent);
	printf("A%d:assigned compute continue, Pr[",nodeNumber);
	rv->printNameFrameValue(stdout,false);
	if (message(Mega)) {
	  if (rv->allParents.size() > 0) {
	    printf("|");
	    printRVSetAndValues(stdout,rv->allParents,false);
	  }
	} else {
	  if (rv->allParents.size() > 0) {
	    printf("|parents");
	  }
	}
	printf("]=%f,crClqPr=%f\n",cur_p.val(),p.val());
      }
      if (!cur_p.essentially_zero()) {
	// Continue, do not update probability!!
	ceIterateAssignedNodesRecurse(sharedStructure,
				      cliqueBeamThresholdEstimate,
				      maxCEValue,
				      nodeNumber+1,p);
      }
    }
    break;

  default:
    assert(0);
    break;
  }
  if (message(High+5))
    traceIndent--;

}



/*-
 *-----------------------------------------------------------------------
 * MaxClique::ceIterateAssignedNodesNoRecurse()
 *
 *    This is a non-recursive version of
 *    ceIterateAssignedNodesRecurse(). The reason for this routine is
 *    to try to produce code that has better branch prediction
 *    behavior and that also creates fewer temporary variables on
 *    routine calls. Timings show that GMTK is about 10-15% faster
 *    using this routine.
 *
 *    This routine is a bit more involved than the recursive version
 *    above, as it uses customize loops (with much use of the evil
 *    goto statements, so the code is a bit ugly). It is essentially
 *    assembly code but expressed in C. Note that the recursive
 *    version is called as a backup when there are no assigned nodes
 *    in a clique or when verbose printing is turned on. For that
 *    reason, and also for pedagogical reasons, the recursive version
 *    is still in place. Any modifications not associated with verbose
 *    printing will need to be done in both places.
 *
 *
 *  TODO: when we find that all values of a variable cause a prune, don't back up to the previous variable (which 
 *        might not be a parent),
 *        rather back up to the latest parent of that variable and change its value. We might even
 *        back out to the separators.
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
MaxCliqueTable::ceIterateAssignedNodesNoRecurse(MaxCliqueTable::SharedLocalStructure& sharedStructure,
						logpr cliqueBeamThresholdEstimate,
						logpr& maxCEValue,
						const logpr p)
{
  // syntactic convenience variables.
  MaxClique& origin = *(sharedStructure.origin);

  if (p*origin.sortedAssignedContinuationScores[0]  <= cliqueBeamThresholdEstimate)
    return;


  // parray has to be 1 offset, storing p in entry -1
  logpr* parray = origin.probArrayStorage.ptr + 1;
  parray[-1] = p;

  // Another option: compute a local clique beam estimate that can be used
  //    to compare directly against cur_p.
  // set as follows:
  //    cliqueBeamThresholdEstimate = cliqueBeamThresholdEstimate/p;
  //    cliqueBeamThresholdEstimate = cliqueBeamThresholdEstimate/parray[nodeNumber];
  //
  // or could subtract off continuation scores from an array estimate.

  int nodeNumber;
  logpr cur_p;
  for (nodeNumber=0;nodeNumber<(int)sharedStructure.fSortedAssignedNodes.size();nodeNumber++) {
    RV* rv = sharedStructure.fSortedAssignedNodes.ptr[nodeNumber];
    switch (origin.dispositionSortedAssignedNodes.ptr[nodeNumber]) {

    case MaxClique::AN_CPT_ITERATION_COMPUTE_AND_APPLY_PROB: 
      rv->begin(cur_p);
      goto applyProbTag;
    case MaxClique::AN_COMPUTE_AND_APPLY_PROB:
      rv->probGivenParents(cur_p);
      {
      applyProbTag:
	// check for possible zero (could occur with
	// zero score or observations).
	parray[nodeNumber] = parray[nodeNumber-1]*cur_p;
	if (parray[nodeNumber]*origin.sortedAssignedContinuationScores[nodeNumber] 
	    <= cliqueBeamThresholdEstimate) {
	  // Since we have small here, we cancel iterations of all
	  // subsequent clique variables right now, rather than
	  // iterate them with what will end up being below threshold
	  // probability (assuming the Gaussians are never > 1).
	  parray[nodeNumber].set_to_zero(); // @@@ REMOVE
	  nodeNumber++;
	  goto end_of_initial_clique_value;
	}
      }
      break;

    case MaxClique::AN_CPT_ITERATION_COMPUTE_PROB_REMOVE_ZEROS:
      rv->begin(cur_p);
      goto removeZeroTag;
    case MaxClique::AN_CONTINUE_COMPUTE_PROB_REMOVE_ZEROS:
      rv->probGivenParents(cur_p);
      {
      removeZeroTag:
	// check for possible zero (could occur with
	// zero score or observations).
	// TODO: could do cpbeam pruning here and see if 
	//       parray[nodeNumber-1]*cur_p falls below threshold pruning if it does.
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

    case MaxClique::AN_CONTINUE:
      parray[nodeNumber] = parray[nodeNumber-1];
      break;

    }
  }

 end_of_initial_clique_value:

  // get ready for main loop
  const bool imc_nwwoh_p = (origin.packer.packedLen() <= IMC_NWWOH);
  const int maxNodeNumber = (int)sharedStructure.fSortedAssignedNodes.size()-1;
  nodeNumber--;

  MaxClique::AssignedNodeDisposition cur_disp = origin.dispositionSortedAssignedNodes.ptr[nodeNumber];
  RV* cur_rv = sharedStructure.fSortedAssignedNodes.ptr[nodeNumber];

  // check if zero probability, and if so, skip the first one and
  // continue on.
  if (parray[nodeNumber]*origin.sortedAssignedContinuationScores[nodeNumber]  
      <= cliqueBeamThresholdEstimate)
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
			   (unsigned**)sharedStructure.discreteValuePtrs.ptr,
			   (unsigned*)&(cliqueValues.ptr[numCliqueValuesUsed].val[0]));
      } else {


#ifdef USE_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL
	const unsigned long lindex = numCliqueValuesUsed*origin.packer.packedLen();
	if (lindex >= origin.temporaryCliqueValuePool.size()) {
	  // use aggressive growth factor for now to avoid expensive copies.
	  origin.temporaryCliqueValuePool.resizeAndCopy(
							origin.packer.packedLen()*
							int(1.5+(double)origin.temporaryCliqueValuePool.size()*TEMPORARY_LOCAL_CLIQUE_VALUE_POOL_GROWTH_RATE));
	}
	unsigned *pcv = 
	  &origin.temporaryCliqueValuePool.ptr[lindex];
	origin.packer.pack((unsigned**)sharedStructure.discreteValuePtrs.ptr,(unsigned*)pcv);
	// store integer value of the location.
	cliqueValues.ptr[numCliqueValuesUsed].ival = lindex;

#if 0
	// grab pointer to storage.
	unsigned *pcv = origin.localValueHolder.curCliqueValuePtr();
	origin.packer.pack((unsigned**)discreteValuePtrs.ptr,(unsigned*)pcv);
	// claim stored value.
	origin.localValueHolder.allocateCurCliqueValue();
	// store pointer to appropriate location.
	cliqueValues.ptr[numCliqueValuesUsed].ptr = pcv;
#endif

#else

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
#ifdef TRACK_NUM_CLIQUE_VALS_SHARED	
	else
	  numCliqueValuesShared++;
#endif
	// Save the pointer to whatever the hash table decided to use.
	cliqueValues.ptr[numCliqueValuesUsed].ptr = key;


#endif

      }
      // save the probability
      cliqueValues.ptr[numCliqueValuesUsed].p = final_p;
      numCliqueValuesUsed++;

      /*
       * uncomment this next code to produce lots of messages.
       * note above, if high verbosity is on, recursive version of this
       * routine will print this information.
       */
      /*
	if (message(Mega)) {
	// psp2(stdout,spi*traceIndent);
	infoMsg(Mega,"Inserting New Clique Val,pr=%f,sm=%f: ",
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
	  do {
	    // continue going until we get one that is unfinished and within beam.
	    unfinished = cur_rv->next(cur_p);
	    if (unfinished) {
	      register logpr tmp = parray[nodeNumber-1]*cur_p;
	      if (tmp*origin.sortedAssignedContinuationScores[nodeNumber]  
		  > cliqueBeamThresholdEstimate) {
		parray[nodeNumber] = tmp;
		goto next_node_number;
	      }
	    } else
	      break;
	  } while (1);
	}
	break;

      case MaxClique::AN_CPT_ITERATION_COMPUTE_PROB_REMOVE_ZEROS:
	{
	  // guaranteed not to get zero prob here. Note that
	  // we get the probability into cur_p, but we do not
	  // use it in this case. Also, parray[nodeNumber] is
	  // already up to date from the "begin()" at the beginning
	  // of this variables iteration.
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
	  // then we're really done since we're finished with the
	  // first clique variable.
	  return;
	} else {
	  // we're finished with the current variable, need to
	  // continue on with previous variable.
	  nodeNumber--;
	  cur_disp = origin.dispositionSortedAssignedNodes.ptr[nodeNumber];
	  cur_rv = sharedStructure.fSortedAssignedNodes.ptr[nodeNumber];
	}
      }

      // still here? Continue on incrementing previous variable.
    } while (1);

  next_node_number:
    while (nodeNumber < maxNodeNumber) {
      nodeNumber++;
      cur_disp = origin.dispositionSortedAssignedNodes.ptr[nodeNumber];
      cur_rv = sharedStructure.fSortedAssignedNodes.ptr[nodeNumber];

      // need to fill up the rest of the table.
      switch (cur_disp) {
      case MaxClique::AN_CPT_ITERATION_COMPUTE_AND_APPLY_PROB: 
	{
	  cur_rv->begin(cur_p);
	  // might get a zero probability, check condition here.
	  // TODO: keep result in temp and only write back out if need be. @@@
	  register logpr tmp =  parray[nodeNumber-1]*cur_p;
	  if (tmp*origin.sortedAssignedContinuationScores[nodeNumber] 
	      <= cliqueBeamThresholdEstimate) {
	    // We just did a begin and got zero on the first try. 
	    // we need to continue on with this variable until we finish.
	    goto next_iteration;
	  } else {
	    // only write when we have to.
	    parray[nodeNumber] = tmp;
	  }
	}
	break;

      case MaxClique::AN_COMPUTE_AND_APPLY_PROB:
	{
	  cur_rv->probGivenParents(cur_p);
	  // might get a zero probability, check condition here.
	  register logpr tmp = parray[nodeNumber-1]*cur_p;
	  if (tmp*origin.sortedAssignedContinuationScores[nodeNumber] 
	      <= cliqueBeamThresholdEstimate) {
	    // Since we have zero here, we cancel iterations of all
	    // subsequent variables right now, rather than iterate
	    // them with what will end up being zero probability.
	    nodeNumber--;
	    cur_disp = origin.dispositionSortedAssignedNodes.ptr[nodeNumber];
	    cur_rv = sharedStructure.fSortedAssignedNodes.ptr[nodeNumber];
	    goto next_iteration;
	  } else
	    parray[nodeNumber] = tmp;
	}
	break;

      case MaxClique::AN_CPT_ITERATION_COMPUTE_PROB_REMOVE_ZEROS:
	{
	  cur_rv->begin(cur_p);
	  // might get a zero probability, check condition here.
	  if (cur_p.essentially_zero()) {
	    goto next_iteration;
	  } else {
	    // TODO: update post-condition comments in RV and CPT iterators as
	    // to when a zero RV can occur and when not.
	    parray[nodeNumber] = parray[nodeNumber-1];
	  }
	}
	break;

      case MaxClique::AN_CONTINUE_COMPUTE_PROB_REMOVE_ZEROS:
	{
	  cur_rv->probGivenParents(cur_p);
	  // might get a zero probability, check condition here.
	  if (cur_p.essentially_zero()) {
	    // Since we have zero here, we cancel iterations of all
	    // subsequent variables right now, rather than iterate them
	    // with what will end up being zero probability.
	    nodeNumber--;
	    cur_disp = origin.dispositionSortedAssignedNodes.ptr[nodeNumber];
	    cur_rv = sharedStructure.fSortedAssignedNodes.ptr[nodeNumber];
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

      case MaxClique::AN_CONTINUE:
	parray[nodeNumber] = parray[nodeNumber-1];
	break;

      }
    }

  } while (1);

}




/*-
 *-----------------------------------------------------------------------
 * MaxCliqueTable::clearCliqueAndIncommingSeparatorMemory()
 *
 *    memory clearing routine, this routine clears all significant memory
 *    associated with this clique and all its incomming separators.
 *    It should only be called when using this clique in a collect-evidence
 *    form, where the clique memory is never going to be used again (such as
 *    during prob(Evidence) form of inference, where the goal is just to
 *    compute the probability of evidence.
 *
 *
 * Preconditions:
 *
 *   The basic data structures should be set up.
 *
 * Postconditions:
 *    The memory has been completely freed. It should be possible to reconstruct the
 *    separator again though. Also, if the 'alsoClearOrigins=true', then
 *    we can no longer use an instance of this clique again without constructing
 *    it from scratch (i.e., it would always do a hash insert).
 *    
 *
 * Side Effects:
 *    Changes 
 *      1) all of the memory associated with this clique and its incomming separators
 *      2) if 'alsoClearOrigins' is true, it will delete all of the memory associated with
 *         the origin of the clique.
 *
 * Results:
 *     nothing
 *
 *-----------------------------------------------------------------------
 */
void 
MaxCliqueTable::clearCliqueAndIncommingSeparatorMemory(MaxCliqueTable::SharedLocalStructure& sharedStructure,
						       ConditionalSeparatorTable* separatorTableArray,
						       ConditionalSeparatorTable::SharedLocalStructure* sepSharedStructureArray)
{
  MaxClique& origin = *(sharedStructure.origin);
  // first do the separators
  for (unsigned sepNumber=0;sepNumber<origin.ceReceiveSeparators.size();sepNumber++) {
    ConditionalSeparatorTable& sep = 
      separatorTableArray[origin.ceReceiveSeparators[sepNumber]];
    sep.clearInferenceMemory();
  }
  // and next do self.
  clearInferenceMemory();
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
 *    to do a separate pruning stage by calling ceCliqueBeamPrune().
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
MaxCliqueTable::
ceSendToOutgoingSeparator(MaxCliqueTable::SharedLocalStructure& sharedStructure,
			  ConditionalSeparatorTable* separatorTableArray,
			  ConditionalSeparatorTable::SharedLocalStructure* sepSharedStructureArray)
{
  // syntactic convenience variables.
  MaxClique& origin = *(sharedStructure.origin);
  ceSendToOutgoingSeparator(sharedStructure,
			    separatorTableArray[origin.ceSendSeparator],
			    sepSharedStructureArray[origin.ceSendSeparator]
			    );
}
void 
MaxCliqueTable::
ceSendToOutgoingSeparator(MaxCliqueTable::SharedLocalStructure& sharedStructure,
			  ConditionalSeparatorTable& sep,
			  ConditionalSeparatorTable::SharedLocalStructure& sepSharedStructure)
{

  // printf("At start of sending to outgoing separator\n");

  // keep a local variable copy of this around to avoid potential
  // dereferencing.  This one cannot be const since it might change
  // during a resize, in which case we need to reassign this variable.
  ConditionalSeparatorTable::AISeparatorValue * 
    sepSeparatorValuesPtr = sep.separatorValues->ptr; 

#ifdef TRACK_NUM_CLIQUE_VALS_SHARED  
  infoMsg(IM::High-2,"ceSendToOutgoingSep: from clique w state space = %d. NumShared = %d, %2.2f percent\n",
	  numCliqueValuesUsed,numCliqueValuesShared, 100*(float)numCliqueValuesShared/(float)numCliqueValuesUsed);
#else
  infoMsg(IM::High-2,"ceSendToOutgoingSep: from clique w state space = %d.\n",numCliqueValuesUsed);
#endif

  // syntactic convenience variables.
  MaxClique& origin = *(sharedStructure.origin);
  SeparatorClique& sepOrigin = 
    *(sepSharedStructure.origin);  

  // first check if this is an all "observed" clique
  if (origin.hashableNodes.size() == 0) {
    // Everything in this clique is observed or deterministic.  Therefore, there should
    // be one and only one clique value in this clique which is the
    // probability of the assigned probability nodes in this clique.

    // Since this clique is a superset of the separator, it means also
    // that the seperator will consist only of observed nodes. We
    // therefore create one accumulated entry and one rem entry and
    // insert our current probability.

    // printf("=====>>>> Observed clique in forward pass\n");

    // first do some sanity checks.
    assert (sepOrigin.hAccumulatedIntersection.size() == 0);
    assert (sepOrigin.hRemainder.size() == 0);
    assert (sepSharedStructure.remDiscreteValuePtrs.size() == 0);

    const unsigned accIndex = 0;

    // We here keep handy reference for readability.
    // To find out where has this been allocated, 
    // search for tag: ALLOCATE_REMVALUES_ALL_OBSERVED.
    ConditionalSeparatorTable::AISeparatorValue& sv
      = sepSeparatorValuesPtr[accIndex];

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

  // Note, we could optionally do pruning here rather than where we do
  // it now, namely when we construct the clique. There may be some
  // utility in doing it here (delaying the pruning) since the values
  // at the end of the pruned clique table could be used again before
  // they are deallocated.

  // TODO: do sampling code here, i.e., optionally sample from
  //       remaining portion of what otherwise would be pruned away
  //       portion of the clique.

  // Now, we send the clique to the outgoing separator.

  // next check if the outgoing separator has only obseved values.
  if (sepOrigin.hAccumulatedIntersection.size() == 0 
      && sepOrigin.hRemainder.size() == 0) {

    // printf("Observed separator\n");

    // Then indeed, outgoing separator is all observed values. We do
    // this special separately from the general case since that case
    // is already getting a bit unwieldy.

    // TODO: this can also probably handle the case of disconnected networks.
    // 

    ConditionalSeparatorTable::AISeparatorValue& sv
      = sepSeparatorValuesPtr[0];

    // This must be first time for this entry.  We never need more
    // than one rem value.  Search for tag
    // 'ALLOCATE_REMVALUES_ALL_OBSERVED' in this file for where else
    // this could be done, but we allocate it here.
    if (sv.remValues.size() == 0)
      sv.remValues.resize(1); 
    // in all cases, we use only one entry when the separator is observed.
    sv.numRemValuesUsed = 1;

    // Just sum up the clique entries projecting down into the
    // observed separator, but only do the ones that pass the beam
    // threshold. Arguably, we might not want or need to do beam
    // pruning here since when the separator is all observed, any
    // sparseness that we introduce by clique pruning isn't going to
    // change things later on (since everything is projected to the
    // same point), but we do it here anyway for numerical consistency
    // with the general case.
    for (unsigned cvn=0;cvn<numCliqueValuesUsed;) {
      if (JunctionTree::viterbiScore) {
	// sv.remValues.ptr[0].p.assign_if_greater(cliqueValues.ptr[cvn].p);
	// TODO: add k-best
	if (cliqueValues.ptr[cvn].p > sv.remValues.ptr[0].p) {
	  sv.remValues.ptr[0].p = cliqueValues.ptr[cvn].p;
	  sv.remValues.ptr[0].backPointer = cvn;
	}
      } else {
	sv.remValues.ptr[0].p += cliqueValues.ptr[cvn].p;
      }
      cvn++;
    }

  } else {
    // We are guaranteed that either we have an accumulated
    // intersection, a reminder, or both (but not neither).  Go
    // through clique values and accumulate into appropriate place
    // within separator.

    // pre-compute a few values that the compiler might spend more
    // time than necessary to re-check.
    const bool imc_nwwoh_p = (origin.packer.packedLen() <= IMC_NWWOH);
    const bool isc_nwwoh_ai_p = (sepOrigin.accPacker.packedLen() <= ISC_NWWOH_AI);
    const bool isc_nwwoh_rm_p = (sepOrigin.remPacker.packedLen() <= ISC_NWWOH_RM);
    const bool sep_origin_hAccumulatedIntersection_exists_p =
      (sepOrigin.hAccumulatedIntersection.size() > 0);
    const bool sep_remDiscreteValuePtrs_exists_p = 
      (sepSharedStructure.remDiscreteValuePtrs.size() > 0);

    for (unsigned cvn=0;cvn<numCliqueValuesUsed;) {{

	// printf("Iteration through clique, iter = %d\n",cvn);

	// TODO: optimize away this conditional check. (and/or use const
	// local variable to indicate it wont change)
	if (imc_nwwoh_p) {
	  origin.packer.unpack((unsigned*)&(cliqueValues.ptr[cvn].val[0]),
			       (unsigned**)sharedStructure.discreteValuePtrs.ptr);
	} else {
	  origin.packer.unpack((unsigned*)cliqueValues.ptr[cvn].ptr,
			       (unsigned**)sharedStructure.discreteValuePtrs.ptr);
	}
	for (unsigned j=0;j<sharedStructure.fDeterminableNodes.size();j++) {
	  RV* rv = sharedStructure.fDeterminableNodes.ptr[j];
	  RV2DRV(rv)->assignDeterministicChild();
	}


	// All hidden random variables now have their discrete
	// value. Accumulate this probability into the given separator.

	/*
	 * There are 3 cases.
	 * 1) Both AI exists and REM exist.
	 * 2) AI exists but REM doesnt exist.
	 * 3) AI does not exist, but REM exists.
	 * x) AI not exist and REM not exist can't occur (covered above).
	 */

	unsigned accIndex;
	// TODO: optimize this check away out of loop.
	if (sep_origin_hAccumulatedIntersection_exists_p) { 
	  // an accumulated intersection exists.

	  // make sure there is at least one available accumulated intersection entry
	  assert ( sep.numSeparatorValuesUsed <= sep.separatorValues->size());
	  if (sep.numSeparatorValuesUsed >= sep.separatorValues->size()) {
	    const unsigned old_size = sep.separatorValues->size();
	    // TODO: optimize this size re-allocation.
	    if (sep.numSeparatorValuesUsed >= sepOrigin.separatorValueSpaceManager.currentSize()) 
	      sepOrigin.separatorValueSpaceManager.advanceToNextSize();
	    sep.separatorValues->resizeAndCopy(sepOrigin.separatorValueSpaceManager.currentSize()); 
	    sepSeparatorValuesPtr = sep.separatorValues->ptr;
	    if (isc_nwwoh_ai_p) {
	      // Then the above resize just invalided all our pointers to keys,
	      // but it did not invalidate the array indices. Go through
	      // and correct the keys within the hash table.
	      // TODO: think of a better way to do this that also looses no efficiency.
	      for (unsigned i=0;i<sep.iAccHashMap->tableSize();i++) {
		if (!sep.iAccHashMap->tableEmpty(i)) {
		  sep.iAccHashMap->tableKey(i)
		    = &(sepSeparatorValuesPtr[sep.iAccHashMap->tableItem(i)].val[0]);
		}
	      }
	    }
	    const unsigned new_size = sep.separatorValues->size();
	    // if (sep.remDiscreteValuePtrs.size() > 0) {
	    if (sep_remDiscreteValuePtrs_exists_p) {
	      for (unsigned i=old_size;i<new_size;i++) {
		// re-construct hash tables only for new entries.
		new (&sepSeparatorValuesPtr[i].iRemHashMap)
		  VHashMapUnsignedUnsignedKeyUpdatable
		  (sepOrigin.remPacker.packedLen(),REM_HASH_MAP_STARTING_SIZE);
		// TODO: potentially preallocate default size of  
		// separatorValues->ptr[i].remValues.resize(default);
		// TODO: potentially create zero size here, and only
		//       grow bigger when we start adding things.
	      }
	    }
	  }

	  unsigned *accKey;
	  // TODO: optimize this check out of loop.
	  if (isc_nwwoh_ai_p) {
	    accKey = &(sepSeparatorValuesPtr[sep.numSeparatorValuesUsed].val[0]);
	    sepOrigin.accPacker.pack((unsigned**)sepSharedStructure.accDiscreteValuePtrs.ptr,
				     accKey);
	  } else {
	    accKey = sepOrigin.accValueHolder.curCliqueValuePtr();
	    sepOrigin.accPacker.pack((unsigned**)sepSharedStructure.accDiscreteValuePtrs.ptr,
				     accKey);
	    // check if this value combination already lives in
	    // origin's value holder hash table and if so, use that.
	    bool foundp;
	    accKey = sepOrigin.accSepValHashSet.insert(accKey,foundp);
	    if (!foundp) {
	      // only allocate a new value if it was inserted.
	      sepOrigin.accValueHolder.allocateCurCliqueValue();
	    }
	    // store the pointer in case we use it.
	    sepSeparatorValuesPtr[sep.numSeparatorValuesUsed].ptr = accKey;
	  }

	  bool foundp;
	  unsigned* accIndexp =
	    sep.iAccHashMap->insert(accKey,
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
	    ConditionalSeparatorTable::AISeparatorValue& sv
	      = sepSeparatorValuesPtr[accIndex];

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
	ConditionalSeparatorTable::AISeparatorValue& sv
	  = sepSeparatorValuesPtr[accIndex];

	// make sure there is at least one available entry
	assert (sv.numRemValuesUsed <= sv.remValues.size());
	if (sv.numRemValuesUsed >= sv.remValues.size()) {
	  // TODO: optimize this growth rate.
	  // start small but grow fast.
	  // sv.remValues.resizeAndCopy(1+sv.remValues.size()*2); // *3
	  sv.remValues.resizeAndCopy(sepOrigin.remainderValueSpaceManager.nextSizeFrom(sv.remValues.size()));
	  sepOrigin.remainderValueSpaceManager.setCurrentAllocationSizeIfLarger(sv.remValues.size());
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
	  sepOrigin.remPacker.pack((unsigned**)sepSharedStructure.remDiscreteValuePtrs.ptr,
				   remKey);
	} else {
	  // grab pointer to next packed clique value to be used.
	  remKey = sepOrigin.remValueHolder.curCliqueValuePtr();
	  sepOrigin.remPacker.pack((unsigned**)sepSharedStructure.remDiscreteValuePtrs.ptr,
				   remKey);
	  // check if this value combination already lives in
	  // origin's value holder hash table and if so, use that.
	  bool foundp;
	  remKey = sepOrigin.remSepValHashSet.insert(remKey,foundp);
	  if (!foundp) {
	    // only allocate a new value if it was inserted.
	    sepOrigin.remValueHolder.allocateCurCliqueValue();
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

	// printf("Inserted sep value, iter = %d\n",cvn);

      }
    next_iteration:
      cvn++;
    }
  }

  // To reallocate or not to reallocate, that is the question.  here,
  // we just reallocate for now. The reason why we might want to reallocate is that
  // pruning might have adjusted numCliqueValuesUsed in this routine above.
  // 
  // TODO: reallocate only if change is > some percentage (say 5%),
  // and export to command line.
  // TODO: possibly tell origin.cliqueValueSpaceManager about this pruning event.
  // TODO: if -probE option or gmtkDecode is running, no need to re-allocate here since this
  //       will soon be deleted anyway.
  // TODO: if we do a backward pass, we might not want to resize, and if an entry ends up
  //       becoming probable, re-insert the entry as a valid clique entry.
  if (numCliqueValuesUsed < cliqueValues.size())
    cliqueValues.resizeAndCopy(numCliqueValuesUsed);

  // printCliqueEntries(stdout);

  /////////////////////////////////////
  // And prune the separator as well.
  sep.ceSeparatorPrune(sepOrigin);

}


/*-
 *-----------------------------------------------------------------------
 * MaxCliqueTable::ceCliqueBeamPrune()
 *
 *    Collect Evidence, Clique Prune: This routine will prune away
 *    part of a previously instantiated clique based on the current
 *    clique beam width.
 *    
 *    Note that MaxCliqueTable::ceSendToOutgoingSeparator() does
 *    its own pruning, so when using ceSendToOutgoingSeparator(), this
 *    pruning routine does not need to be called (at least with the
 *    same beam width). This routine, however, is kept here since 1)
 *    at the very right clique of the right most partition, we need to
 *    explicitely prune, and 2) the island algorithm sometimes also
 *    needs to explicitly call pruning.
 *
 * Preconditions:
 *   1) the value of the max clique 'maxCEValue' must have been
 *      computed already.
 *
 *   2) clique table must be created, meaning that either:
 *
 *        MaxCliqueTable::ceIterateAssignedNodesRecurse()
 *
 * Postconditions:
 *    Clique table has been pruned, and memory for it has NOT been re-allocated to
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
MaxCliqueTable::ceCliqueBeamPrune(MaxClique& origin,
				  logpr maxCEValue)
{
  // return immediately if beam pruning is turned off.
  if (origin.cliqueBeam == (-LZERO))
    return;

  // create an ininitialized variable
  logpr beamThreshold((void*)0);
  if (origin.cliqueBeam != (-LZERO)) {
    // then we do clique table pruning right here rather
    // than a separate call to ceCliqueBeamPrune().
    // break into the logp to avoid unnecessary zero checking.
    beamThreshold.valref() = maxCEValue.valref() - origin.cliqueBeam;
  } else {
    // set beam threshold to a value that will never cause pruning.
    beamThreshold.set_to_zero();
  }

  const unsigned origNumCliqueValuesUsed = numCliqueValuesUsed;
  for (unsigned cvn=0;cvn<numCliqueValuesUsed;) {
    if (cliqueValues.ptr[cvn].p < beamThreshold) {
      // swap with last entry, and decrease numCliqueValuesUsed by one. We
      // swap so that entries at the end can be added back in by a future stage.
      swap(cliqueValues.ptr[cvn],cliqueValues.ptr[numCliqueValuesUsed-1]);
      numCliqueValuesUsed--;
    } else {
      cvn++;
    }
  }

  infoMsg(IM::Med,"Clique beam pruning: Max cv = %f, thres = %f. Original clique state space = %d, new clique state space = %d, %2.2f%% reduction\n",
	  maxCEValue.valref(),
	  beamThreshold.valref(),
	  origNumCliqueValuesUsed,
	  numCliqueValuesUsed,
	  100*(1.0 - (double)numCliqueValuesUsed/(double)(origNumCliqueValuesUsed>0?origNumCliqueValuesUsed:1)) );

#if 0
  // A version with a bit more information printed.
  infoMsg(IM::Med,"Clique beam pruning, Max cv = %f, thres = %f. Original clique state space = %d, orig sum = %f, new clique state space = %d, new sum = %f\n",
	  maxCEValue.valref(),
	  beamThreshold.valref(),
	  origNumCliqueValuesUsed,
	  numCliqueValuesUsed,
	  sumProbabilities().valref());
#endif


}



/*-
 *-----------------------------------------------------------------------
 * MaxCliqueTable::ceDoAllPruning()
 *
 *    Clique Prune: This routine will do both k-pruning and beam pruning just
 *    like the routine that projects down to the next separator, but this routine does
 *    just the pruning part.
 *
 * Preconditions:
 *   same as other pruning routines
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
MaxCliqueTable::ceDoAllPruning(MaxClique& origin,
			       logpr maxCEValue)
{


  // if all observed and/or deterministic clique, then only one state,
  // so nothing to prune.
  if (origin.hashableNodes.size() == 0)
    return;

  // We seed the random number generator specifically for this clique
  // since it might get called again during island algorithm. 
  // Search for string K-BEAM-SEED elsewhere in this file for further information.
  rnd.seed(&(cliqueValues[0].p.valref()));

  const unsigned long origNumCliqueValuesUsed = numCliqueValuesUsed;
  // printf("DEBUG: orig = %lu, allo = %lu\n",origNumCliqueValuesUsed,cliqueValues.size());

  // TODO: give a command line option to change the order of these pruning
  // options.

  // First, do k-pruning (ideally we would do this after
  // beam pruning, but since beam pruning is integrated
  // into the code above, we do max state pruning first).
  // Prune the minimum of the fixed K size and the percentage size.
  unsigned k;
  k = 2 + (unsigned)((origin.cliqueBeamRetainFraction)*(double)numCliqueValuesUsed);
  if (origin.cliqueBeamMaxNumStates > 0) {
    k = min(k,origin.cliqueBeamMaxNumStates);
  }
  //   printf("nms = %d, pf = %f, ncv = %d, k = %d\n",origin.cliqueBeamMaxNumStates,
  // origin.cliqueBeamRetainFraction,numCliqueValuesUsed,k);
  // printf("starting k pruning with state space %d\n",numCliqueValuesUsed); fflush(stdout);

  if (k < numCliqueValuesUsed) {
    infoMsg(IM::Med,"Clique k-beam pruning with k=%d: Original clique state space = %d\n",k,
	    numCliqueValuesUsed);
    numCliqueValuesUsed = ceCliqueStatePrune(k,cliqueValues.ptr,numCliqueValuesUsed);
  }

  // printf("ending k pruning\n"); fflush(stdout);

  // next do mass pruning.
  numCliqueValuesUsed = ceCliqueMassPrune(origin.cliqueBeamMassRelinquishFraction,
					  origin.cliqueBeamMassExponentiate,
					  origin.cliqueBeamMassFurtherBeam,
					  origin.cliqueBeamMassMinSize,
					  cliqueValues.ptr,
					  numCliqueValuesUsed);

  // next, do normal beam pruning.
  ceCliqueBeamPrune(origin,maxCEValue);

  // do diversity pruning.
  // printf("starting diversity pruning with state space %d\n",numCliqueValuesUsed); fflush(stdout);
  ceCliqueDiversityPrune(origin,origin.cliqueBeamClusterPruningNumClusters);
  // printf("ending diversity pruning\n"); fflush(stdout);

  // last, add random entries back in.
  if (origin.cliqueBeamUniformSampleAmount != 0) {
    ceCliqueUniformSamplePrunedCliquePortion(origin,origNumCliqueValuesUsed);
  }


  // lastly, resize.
  // To reallocate or not to reallocate, that is the question.  here,
  // we just reallocate for now.
  // 
  // TODO: reallocate only if change is > some percentage (say 5%),
  // and export to command line.
  // e.g., if ((origNumCliqueValuesUsed - numCliqueValuesUsed) > 0.05*origNumCliqueValuesUsed)
  // TODO: if -probE option or gmtkDecode is running, no need to re-allocate here since this
  //       will soon be deleted anyway.
  // TODO: if we do a backward pass, we might not want to resize, and if an entry ends up
  //       becoming probable, re-insert the entry as a valid clique entry.

  // resize only if we've shrunk by more than about 1.6%
  // relative. Assume strength reduction will take care of fast
  // divide.


  //   if ((origNumCliqueValuesUsed - numCliqueValuesUsed) > origNumCliqueValuesUsed/64) {
  //     cliqueValues.resizeAndCopy(numCliqueValuesUsed);
  //   }
  // printf("DEBUG: old cond = %d, new cond = %d\n",((origNumCliqueValuesUsed - numCliqueValuesUsed) > origNumCliqueValuesUsed/64),((cliqueValues.size() - numCliqueValuesUsed) > cliqueValues.size()/64));
  if ((cliqueValues.size() - numCliqueValuesUsed) > cliqueValues.size()/64) {
    cliqueValues.resizeAndCopy(numCliqueValuesUsed);
  }


#ifdef USE_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL
  // if this code is enabled, it means that we cannot call pruning more than once.
  if (origin.packer.packedLen() > IMC_NWWOH) {
    // finally, insert surviving entries into global shared pool.
    insertLocalCliqueValuesIntoSharedPool(origin);
    // and free up the local buffer.
    origin.temporaryCliqueValuePool.resize(CLIQUE_VALUE_HOLDER_STARTING_SIZE*origin.packer.packedLen());
  }
#endif

}


#ifdef USE_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL
/*-
 *-----------------------------------------------------------------------
 * MaxCliqueTable::insertLocalCliqueValuesIntoSharedPool();
 *
 *    This routine takes all the clique entries and values, and inserts
 *    the values that live currently in the local clique value pool into
 *    the globally shared clique value pool.
 *
 * Preconditions:
 *   The clique must currently be set up so that the clique entries indicate 
 *   an index into the local clique value pool. Note that the clique value pointers
 *   are such that they give the actual integer index into the exact location
 *   of the current local clique value pool (i.e., it is the entry location
 *   multiplied in by the number of words per clique value).
 *
 *   Also, this routine must not be called unless (origin.packer.packedLen() > IMC_NWWOH) is true.
 *
 * Postconditions:
 *    The clique entries are now set so that the values are pointers into
 *    the potentially shared value location in the globally shared pool.
 *
 * Side Effects:
 *    changes the clique entry value pointers.
 *
 * Results:
 *     nothing
 *
 *-----------------------------------------------------------------------
 */
void 
MaxCliqueTable::insertLocalCliqueValuesIntoSharedPool(MaxClique& origin)
{

  // clique value length in bytes
  const unsigned cvlb = origin.packer.packedLenInBytes();

  unsigned k;
  for (k=0;k<numCliqueValuesUsed;k++) {
    // First, grab pointer to storge where next clique value would
    // be stored if it ends up being used.
    unsigned *pcv = origin.valueHolder.curCliqueValuePtr();
    // next, grab pointer to storge where the local clique value currently lives.
    unsigned *lcv = &origin.temporaryCliqueValuePool.ptr[cliqueValues.ptr[k].ival];
    // copy from the local location to the permanent storage pool.
    // TODO: have a packer optimized copy routine here.
    ::memcpy((void*)pcv,(void*)lcv,cvlb);
    // Look it up in the hash table.
    bool foundp;
    unsigned *key;
    key = origin.cliqueValueHashSet.insert(pcv,foundp);
    if (!foundp) {
      // if it was not found, need to claim this storage that we
      // just used.
      origin.valueHolder.allocateCurCliqueValue();
    } 
#ifdef TRACK_NUM_CLIQUE_VALS_SHARED	
    else
      numCliqueValuesShared++;
#endif
    // Save the pointer to whatever the hash table decided to use.
    cliqueValues.ptr[k].ptr = key;
  }

}
#endif


/*-
 *-----------------------------------------------------------------------
 * MaxCliqueTable::ceCliqueStatePrune(unsigned k)
 *
 *    Collect Evidence, Clique Prune: This routine will organize a 
 *    reviously instantiated clique table so that the top k entries
 *    are at the beginning of the table (in arbitrary order) and the
 *    rest of the entries start at the k+1'st position. It does
 *    this using a fast O(N) algorithm.
 *    It returns k, so that pruning can be implemented by calling
 *    this routine and taking/using only the first k entries on return.
 *    It does not do any change of array allocation, so that can
 *    be done by the caller, and only does array organization
 *    within the bounds of the array so given by the arguments.
 *
 *    We can also quickly find the k top entries for Viterbi decoding
 *    by calling this, and then sorting the top k entries.
 *
 * Preconditions:
 *   1) the value of the max clique 'maxCEValue' must have been
 *      computed already.
 *
 *   2) clique table must be created, meaning that either:
 *
 *        MaxCliqueTable::ceIterateAssignedNodesRecurse()
 *
 * Postconditions:
 *
 *    Clique table has been reorganized so that the k top entries are
 *    at the beginning of the table. Memory for it has NOT been
 *    re-allocated and no entries have been removed.
 *
 * Side Effects:
 *    reorganizes the clique table. Note that this routine does
 *    nothing if k >= the size of the table. 
 *
 * Results:
 *     min(k,length of the given array)
 *
 *-----------------------------------------------------------------------
 */

unsigned
MaxCliqueTable::ceCliqueStatePrune(const unsigned k,
				   CliqueValue* curCliqueVals,
				   const unsigned curNumCliqueValuesUsed)
{
  // k can't be larger than the number of clique entries.
  if (k == 0 || k >= curNumCliqueValuesUsed)
    return curNumCliqueValuesUsed;

  // inclusive range to process.
  unsigned lower = 0;
  unsigned upper = curNumCliqueValuesUsed-1;


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
    if (curCliqueVals[pl1].p < curCliqueVals[pl2].p) {
      if (curCliqueVals[pl2].p < curCliqueVals[pl3].p) {
	pl = pl2;
      } else {
	pl = pl3;
      }
    } else {
      if (curCliqueVals[pl1].p < curCliqueVals[pl3].p) {
	pl = pl1;
      } else {
	pl = pl3;
      }
    }

    // Uncomment to give a valid but fixed deterministic pivot just for testing.
    // pl = (lower+upper)/2;

    logpr pivot = curCliqueVals[pl].p;
    // printf("pivot location = %d, pivot = %f\n",pl,pivot);

    // swap pivot into first position, which is a valid position for
    // pivot, since pivot >= pivot.
    swap(curCliqueVals[pl],curCliqueVals[lower]);

    unsigned l = lower+1;
    unsigned u = upper;
    while (l < u) {
      if (curCliqueVals[l].p >= pivot) {
	l++;
      } else {
	swap(curCliqueVals[l],curCliqueVals[u]);
	u--;
      }
    }

    // 1) Now all entries with index < l are >= pivot,
    //   which means at this point we know we have
    //   the top l entries at indices [0,l]
    // 2) All entries with index > u are < pivot
    // 3) we have l == u
    // 4) The entry at index u==l is unknown however.


    if (curCliqueVals[l].p >= pivot) {
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
    swap(curCliqueVals[lower],curCliqueVals[l]);
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
	if (curCliqueVals[i].p == pivot) {
	  swap(curCliqueVals[i],curCliqueVals[ll]);
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

  return k;

}


/*
 * structure used only for diversity pruning
 */
struct DistClust {
  // the distance of this point to its cluster rep (its nearest center)
  unsigned distance;
  // the cluster number of this point's cluster (i.e., its nearest center).
  unsigned cluster;
};


/*
 * TODO: this next macro is probably generally useful, so ultimately
 * put it somewhere else.
 */

#ifdef USE_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL
#define EXTRACT_KEY_FROM_CVN_INTO_KEY_P(__cvn,__key_p) \
    if (origin.packer.packedLen() <= IMC_NWWOH) { \
      (__key_p) = &(cliqueValues.ptr[(__cvn)].val[0]); \
    } else { \
      (__key_p) = \
	&origin.temporaryCliqueValuePool.ptr[cliqueValues.ptr[(__cvn)].ival]; \
    }
#else
#define EXTRACT_KEY_FROM_CVN_INTO_KEY_P(__cvn,__key_p) \
    if (origin.packer.packedLen() <= IMC_NWWOH) { \
      (__key_p) = &(cliqueValues.ptr[(__cvn)].val[0]); \
    } else { \
      (__key_p) = (unsigned*)cliqueValues.ptr[(__cvn)].ptr; \
    }
#endif





void 
MaxCliqueTable::ceCliqueDiversityPrune(MaxClique& origin,const unsigned numClusters)
{
  /*
   * there are two forms of clustering algorithm. 
   * 
   * 1) The first one clusters only by diversity (i.e., it chooses for
   * the next center, the point that maximizes the min distance to any
   * cluster. That is, each point has a cluster center that is
   * assigned based on how close it is to that center (each point is
   * assigned to the nearest cluster center), and the next center is
   * the point that is maximally distant from its center. 
   *
   * 2) the next cluster algorithm is quite similar except that when
   * there are multiple points that are maximally distant from their
   * centers, rather than breaking ties arbitrarily, we choose the
   * point that has has the highest current score.
   *
   * To get the first algorithm, do not #define the macro DIVERSITY_PRUNE_SCORE_BASED
   * To get the second algorithm, #define the macro DIVERSITY_PRUNE_SCORE_BASED
   */

#define DIVERSITY_PRUNE_SCORE_BASED


  // k can't be larger than the number of clique entries.
  if ((origin.cliqueBeamClusterBeam == (-LZERO)  && 
       (MaxClique::cliqueBeamClusterPruningMaxStates 
	== NO_PRUNING_CLIQUEBEAMCLUSTERPRUNINGMAXSTATES))
      || numClusters <= 1 
      || numClusters >= numCliqueValuesUsed) {
    return;
    // TODO: note that setting the number of clusters to 1 might lead the user
    // to think that this is similar to non-clustered pruning, but this is
    // actually not the case as the number of clusters to 1 option, based on
    // the code here, turns off diversity pruning.
  }

  infoMsg(IM::Med+9,"Diversity/cluster pruning with state space = %d\n",numCliqueValuesUsed);

  /*
    - ideas to speed up: 
    -- factor malloc outside of this routine.
    -- use smaller data types for arrays
    -- use a faster distance function
    -- use a log k data structure ala Feder&Greene (but then constants might
    be worse, and this would be only for relatively large k)
    -- compile with -fno-exceptions (to possibly make routine calls faster)
    -- 

  */

  unsigned long* centers = new unsigned long[numClusters];

  // printf("begin allocate\n");fflush(stdout);
  DistClust* distClusts = new DistClust[numCliqueValuesUsed];
  // printf("end allocate\n");fflush(stdout);


  // find the entry with the max score. The first one is arbitrary.
  centers[0] = 0;
  logpr maxScore = cliqueValues.ptr[0].p;
  // We could also pick a random point as the first cluster rather
  // than the first one, but this makes each run of the inference different.
  // centers[0] = rnd.uniformOpen(numCliqueValuesUsed);

#ifdef DIVERSITY_PRUNE_SCORE_BASED
  // we need to start with the one that has the highest score.  I.e.,
  // the first center is one having the highest score.
  for (unsigned cvn=1;cvn<numCliqueValuesUsed;cvn++) {
    // TODO: unroll this.
    if (cliqueValues.ptr[cvn].p > maxScore) {
      maxScore = cliqueValues.ptr[cvn].p;
      centers[0] = cvn;       
    }
  }
#else
  // nothing to do here.
#endif

  // get the clique entry corresponding to centers[0] which at
  // this point is the clique entry with the max score.
  unsigned* center_key_p;
  EXTRACT_KEY_FROM_CVN_INTO_KEY_P(centers[0],center_key_p);

  unsigned maxDist = 0;
  unsigned maxIndx = 0;
  maxScore.set_to_zero();

  for (unsigned cvn=0;cvn<numCliqueValuesUsed;cvn++) {
    unsigned* current_key_p;
    EXTRACT_KEY_FROM_CVN_INTO_KEY_P(cvn,current_key_p);
    distClusts[cvn].cluster = 0;
    distClusts[cvn].distance = 
      origin.packer.use_distance(center_key_p,current_key_p);
    // find max while we're at it.
#ifdef DIVERSITY_PRUNE_SCORE_BASED
    if (distClusts[cvn].distance > maxDist
	|| (((distClusts[cvn].distance == maxDist && cliqueValues.ptr[cvn].p > maxScore)))) {
      maxIndx = cvn;
      maxDist = distClusts[cvn].distance;
      maxScore = cliqueValues.ptr[cvn].p;
    }
#else
    if (distClusts[cvn].distance > maxDist) {
      maxIndx = cvn;
      maxDist = distClusts[cvn].distance;
      maxScore = cliqueValues.ptr[cvn].p;
    }
#endif
  }

  unsigned k = 1;
  while (k < numClusters) {

    // assign new center using previously computed max index.
    centers[k] = maxIndx;

    // now we need to update all distances. Each min distance is such
    // that it is either the same (i.e., it's cluster center has not
    // changed) or it has decreased (it has a new cluster center). So
    // this means that the distance can never increase (the distances
    // are monotonically decreasing). Therefore, we only need to check
    // against the newly assigned cluster, not all the rest of them.

    EXTRACT_KEY_FROM_CVN_INTO_KEY_P(centers[k],center_key_p);

    // now compute the new max distance and cluster membership.
    maxDist = 0;
    maxIndx = 0;
    maxScore.set_to_zero();
    for (unsigned cvn=0;cvn<numCliqueValuesUsed;cvn++) {

      unsigned* current_key_p;
      EXTRACT_KEY_FROM_CVN_INTO_KEY_P(cvn,current_key_p);

      unsigned cur_dist = 
	origin.packer.use_distance(center_key_p,current_key_p);

      // now check if the distance from the current point to the new
      // center is lower than its current distance, and if it is then
      // this new point gets a new cluster. We break ties by keeping
      // the point assigned to the old cluster, but it might be better
      // to break ties randomly, or to assign it to the cluster based
      // on the cluster center's score.
      if (cur_dist < distClusts[cvn].distance) {
	distClusts[cvn].distance = cur_dist;
	distClusts[cvn].cluster = k;
      }

      // we need unfortunately to re-compute the new max from scratch
      // since the previous max could have been overwritten. One idea:
      //  1) if each cluster keeps track of its own max, and if a cluster
      //     isn't touched by the finding of a new cluster (i.e., the
      //     cluster stays intact), then the points in that cluster
      //     need not be max searched again, only the max point

#ifdef DIVERSITY_PRUNE_SCORE_BASED
      if (distClusts[cvn].distance > maxDist
	  || (((distClusts[cvn].distance == maxDist && cliqueValues.ptr[cvn].p > maxScore)))) {
	maxIndx = cvn;
	maxDist = distClusts[cvn].distance;
	maxScore = cliqueValues.ptr[cvn].p;
      }
#else
      if (distClusts[cvn].distance > maxDist) {
	maxIndx = cvn;
	maxDist = distClusts[cvn].distance;
	maxScore = cliqueValues.ptr[cvn].p;
      }
#endif
    }
    k++;
  }

  // so now we have all the cluster centers and the assignment of all
  // points to their clusters based on the nearest center.

  // Initialize cluster sizes.
  unsigned*  cluster_sizes = new unsigned[numClusters];
  for (k=0;k<numClusters;k++) {
    cluster_sizes[k] = 0;
  }


  const unsigned origNumCliqueValuesUsed = numCliqueValuesUsed;

  // printf("Trying beam pruning\n");

  // Next do normal beam pruning, this next step costs O(n).
  if (MaxClique::cliqueBeamClusterBeam != (-LZERO)) {

    // First, calculate the max score value within each cluster.
    // Note, default values of logpr are set to zero.
    // We also here calculate the cluster sizes.
    logpr* intra_cluster_max_values = new logpr[numClusters];
    for (unsigned cvn=0;cvn<numCliqueValuesUsed;cvn++) {
      const unsigned clust = distClusts[cvn].cluster;
      cluster_sizes[clust]++;
      if (cliqueValues.ptr[cvn].p > intra_cluster_max_values[clust]) {
	intra_cluster_max_values[clust] = cliqueValues.ptr[cvn].p;
      }
    }


    // now we have the clustering do the pruning. 
    // printf("now we have the clustering do the pruning. \n"); fflush(stdout);
    if (IM::messageGlb(IM::Med+9)) {
      printf("Clique cluster beam pruning: Orig cluster-%d sizes: ",numClusters);
      for (k=0;k<numClusters;k++) {
	printf("%u ",cluster_sizes[k]);
	// turn max values into the needed threshold
      }
      printf("\n");
    }


    // Next, we need to calculate the per-cluster pruning threshold:
    for (k=0;k<numClusters;k++) {
      // turn max values into the needed threshold
      intra_cluster_max_values[k].valref() = 
	intra_cluster_max_values[k].valref() - MaxClique::cliqueBeamClusterBeam;
    }

    // Next, we do the actual pruning, and we do this without reordering the
    // clique table entries.
    for (unsigned cvn=0;cvn<numCliqueValuesUsed;) {
      const unsigned clust = distClusts[cvn].cluster;
      if (cliqueValues.ptr[cvn].p < intra_cluster_max_values[clust]) {
	// then we do a prune of this entry.

	assert ( clust < numClusters );

	cluster_sizes[clust]--;

	// swap with last entry, and decrease numCliqueValuesUsed by
	// one.  We swap rather than assign just in case some other
	// method wants to add them back in (but this is slower). If
	// we are assured that this clustering method is the last one
	// called, we can turn the swap into an assignment.
	swap(cliqueValues.ptr[cvn],cliqueValues.ptr[numCliqueValuesUsed-1]);

	// WARNING!! distance values for current entry (cvn) will be invalid
	// after this next step since we only assign the cluster id
	// member in a DistClust object rather than the entire
	// structure. We do this for the sake of speed, and since the
	// distance value is no longer needed at this point and 
	// presumably below.
	distClusts[cvn].cluster = 
	  distClusts[numCliqueValuesUsed-1].cluster;

	numCliqueValuesUsed--;
      } else {
	// the entry 'cvn' survived, so no pruning.
	cvn++;
      }
    }
    // optionally print information done after beam pruning.
    infoMsg(IM::Med,"Clique cluster beam pruning: Original clique state space = %d, new clique state space = %d, %2.2f%% reduction\n",
	    origNumCliqueValuesUsed,
	    numCliqueValuesUsed,
	    100*(1.0 - (double)numCliqueValuesUsed/(double)(origNumCliqueValuesUsed>0?origNumCliqueValuesUsed:1)) );
    
    if (IM::messageGlb(IM::Med+9)) {
      printf("Clique cluster beam pruning: Post cluster-%d sizes: ",numClusters);
      for (k=0;k<numClusters;k++) {
	printf("%u ",cluster_sizes[k]);
      }
      printf("\n");
    }
    delete [] intra_cluster_max_values;
  } else {
    // still need to calculate cluster sizes
    for (unsigned cvn=0;cvn<numCliqueValuesUsed;cvn++) {
      const unsigned clust = distClusts[cvn].cluster;
      cluster_sizes[clust]++;
    }
  }
  
  // printf("************* trying div state pruning, a = %ul, b = %ul\n",MaxClique::cliqueBeamClusterPruningMaxStates,NO_PRUNING_CLIQUEBEAMCLUSTERPRUNINGMAXSTATES);

  if (MaxClique::cliqueBeamClusterPruningMaxStates < NO_PRUNING_CLIQUEBEAMCLUSTERPRUNINGMAXSTATES) {

    // Now we do k-beam pruning. The algorithm is to:
    //  1) 'sort' entire list by cluster number, O(n)
    //  2) run k pruning within each set, amortized O(n) cost.
    //  3) reorganize and put pruned guys at the bottom.
    // total pruning cost is O(n log n) (which is probably ok).

    // - need to modify k-prune to take an arguemnt to pointer to clique
    //  valus + length and return new length (and have that array is
    //  modified to have top k at the top).

    // we can use the STL to create a container class that sorts
    // one array using keys from the other.
    // Also note, if for given s and k, sk >= current size, then just return.


    // First, we need to sort the clique values using the clusters,
    // rather than the values, as the key. This will organize the table
    // so that cluster 0 will come first, cluster 1 will come next, and
    // so on. This can be done in O(n) time rather than doing a sort.
    // Once done, the entries will be in arbitrary order in each cluster.

    // now calculate the cluster sizes
    unsigned*  cluster_starts = new unsigned[numClusters+1];
    unsigned*  cluster_endp = new unsigned[numClusters];
    cluster_endp[0] = cluster_starts[0] = 0;
    for (k=1;k<numClusters;k++) {
      cluster_endp[k] = 
	cluster_starts[k] = cluster_starts[k-1] + cluster_sizes[k-1];
    }
    // include this one for convenience to the code below.
    cluster_starts[numClusters] = numCliqueValuesUsed;

    unsigned cvn=0;
    // pos_cvn_k is the cluster that the position 'cvn' *should* be.
    // Since we are sweeping increasing in cvn, cvn_pos_k starts
    // at zero and increments up to 'numClusters' and switches
    // at teh boundaries given by the 'cluster_starts' array.
    unsigned pos_cvn_k = 0; // this corresponds to cvn == 0.

    // we use a separate index since the loop never needs to occur
    // more than numCliqueValuesUsed times.
    for (unsigned i=0;i<numCliqueValuesUsed;i++) {

      // what is the cluster of current entry at position cvn
      const unsigned cur_cvn_k = distClusts[cvn].cluster;

      if (cur_cvn_k == pos_cvn_k) {
	// then the guy at postion cvn is already in the right partition
	// of the array.
	if (cvn == cluster_endp[cur_cvn_k]) {
	  // then we haven't accounted for this cluster member
	  // yet and we must do so now.
	  cluster_endp[cur_cvn_k]++;
	} else {
	  assert ( cvn < cluster_endp[cur_cvn_k] );
	  // otherwise, this cluster member must be here because
	  // we swaped it here before.
	}

	// move to next entry
	cvn++;

	// deduce and update what the cluster should be of the next cvn.
	if (cvn >= cluster_starts[pos_cvn_k+1])
	  pos_cvn_k++;
      } else {
	// swap the curent one to where it belongs and continue
	swap(cliqueValues.ptr[cvn],cliqueValues.ptr[cluster_endp[cur_cvn_k]]);

	// this next step does not swap the distances only the cluster
	// member items, so leaves the distances totally invalid. This
	// is ok, presumably since the distances are no longer being
	// used and may be a bit faster since there is less swapping
	// going on.
	swap(distClusts[cvn].cluster,
	     distClusts[cluster_endp[cur_cvn_k]].cluster);
	// we've got a new end placement
	cluster_endp[cur_cvn_k]++;
      }
    }

    // Next, we do k-beam pruning within each cluster. We call the k-pruning
    // routine which will leave the top k clique entries at the beginning
    // of each cluster and the rest will be after the top k. We place
    // the new sizes of the clusters in the 'cluster_endp' array -- note
    // that the sizes are at most k but could be < k in which case
    // no pruning is done of course. The entire amortized process is O(n)
    // since each sub-step is O(size of cluster k) and 
    // O(n) = \sum_k O(size of cluster k)
    unsigned newNumCliqueValuesUsed = 0;
    for (k=0;k<numClusters;k++) {
      // what is returned is the new cluster size, which is
      // the min of the original size or the within-cluster state space.
      cluster_endp[k] = 
	ceCliqueStatePrune(MaxClique::cliqueBeamClusterPruningMaxStates,
			   cliqueValues.ptr + cluster_starts[k],
			   cluster_sizes[k]);
      newNumCliqueValuesUsed += cluster_endp[k];
    }

    assert ( newNumCliqueValuesUsed <= numCliqueValuesUsed );

    // Last, we need to re-organize the entire clique table so that
    // all of the to-be pruned entries are at the end. We do the
    // sweep in place, leading to another O(n) algorithm.
    
    // For better cache behavior, we have two pointers scan the array,
    // one that is ahead of the other. The earlier pointer points to
    // the next element that needs to be pruned and the later pointer
    // points to the next element that should not be pruned.  we
    // always maintain that earlier_p < later_p as we scan the
    // array. Note that later_p could touch (and prefetch into cache)
    // some memory items that are later used by earlier_p (this
    // wouldn't be the case if we used an alternative approach where
    // we had two pointers, one initalized at the start and one at the
    // end of the array, and where the pointers appraoch each other
    // (respectively incrementing/decrementing). 

    if (newNumCliqueValuesUsed < numCliqueValuesUsed) {
      // We first find the first "hole", i.e., location where
      // something should be moved.
      for (k=0;k<numClusters;k++) {
	if (cluster_endp[k] < cluster_sizes[k])
	  break;
      }
      assert (k < numClusters);

      // next, do the reorganization.
      unsigned early_loc = cluster_starts[k] + cluster_endp[k];
      k++;
      unsigned later_loc = cluster_starts[k];
      unsigned later_loc_end;
      if (k < numClusters) 
	later_loc_end = cluster_starts[k] + cluster_endp[k];
      else 
	later_loc_end = numCliqueValuesUsed;
      while (early_loc < newNumCliqueValuesUsed) {
	// swap earlier and later pointers
	swap(cliqueValues.ptr[early_loc],
	     cliqueValues.ptr[later_loc]);

	// update earlier pointer
	early_loc++;

	// update later pointer
	later_loc++;
	if (later_loc == later_loc_end) {
	  k++;
	  if (k < numClusters) 
	    later_loc_end = cluster_starts[k] + cluster_endp[k];
	  else 
	    later_loc_end = numCliqueValuesUsed;
	}
      }
    }


    numCliqueValuesUsed = newNumCliqueValuesUsed;


    infoMsg(IM::Med,"Clique cluster state pruning: Original clique state space = %d, new clique state space = %d, %2.2f%% reduction\n",
	    origNumCliqueValuesUsed,
	    numCliqueValuesUsed,
	    100*(1.0 - (double)numCliqueValuesUsed/(double)(origNumCliqueValuesUsed>0?origNumCliqueValuesUsed:1)) );

    if (IM::messageGlb(IM::Med+9)) {
      printf("Clique cluster state pruning: Post cluster-%d sizes: ",numClusters);
      for (k=0;k<numClusters;k++) {
	printf("%u ",cluster_endp[k]);
      }
      printf("\n");
    }

    // free memory for k-beam cluster pruning
    delete [] cluster_starts;
    delete [] cluster_endp;

  }


  // TODO: could still do an r (percentage prune) pruning within each
  // cluster, or a mass-based prune.


  // free memory for general stuff.
  delete [] centers;
  delete [] distClusts;
  delete [] cluster_sizes;

}


/*-
 *-----------------------------------------------------------------------
 * MaxCliqueTable::ceCliqueMassPrune(double removeFraction, unsigned minSize)
 *
 *    Collect Evidence, Clique Mass Prune: This routine will prune away
 *    part of a previously instantiated clique so that it has only
 *    'fraction' of the clique mass left. It does this by sorting
 *    and choosing the top entries such that the mass is retained. This
 *    is based on an idea by Andrew McCallum (2005).
 *
 * Preconditions:
 *   1) clique table must be created, meaning that either:
 *
 *        MaxCliqueTable::ceIterateAssignedNodesRecurse()
 *
 * Postconditions:
 *    Clique table has been pruned, and memory for it has NOT been re-allocated to
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

unsigned
MaxCliqueTable::ceCliqueMassPrune(const double removeFraction,
				  const double exponentiate,
				  const double furtherBeam,
				  const unsigned minSize,
				  CliqueValue* curCliqueVals,
				  const unsigned curNumCliqueValuesUsed)
{

  if (removeFraction <= 0.0 || curNumCliqueValuesUsed <= minSize)
    return curNumCliqueValuesUsed;

  // sort all current clique values descending based on clique mass value.
  // Note, sorting won't be correct if the exponent is negative.
  if (exponentiate < 0) {
    error("ERROR: trying to do exponentiated mass clique pruning with a negative exponent (%e). Exponent must be non-negative for sensible pruning.\n",exponentiate);
  }
  sort(curCliqueVals,curCliqueVals + curNumCliqueValuesUsed,CliqueValueDescendingProbCompare());

  // logpr loc_maxCEValue = curCliqueVals[0].p;
  // printf("mass pruning: maxVal %.18e, minVal %.18e\n",
  // loc_maxCEValue.val(),curCliqueVals[curNumCliqueValuesUsed-1].p.val());

  logpr origSum = sumExponentiatedProbabilities(exponentiate); // /loc_maxCEValue;

  if (origSum.zero())
    return curNumCliqueValuesUsed;

  // logpr desiredSum = origSum* (1.0 - removeFraction);
  logpr desiredSum = (origSum - (origSum*removeFraction));
  logpr actualSum;

#if 0
  printf("DEBUG: mass pruning: origSum %.18e, desiredSum %.18e, diff %.18e\n",
   	 origSum.val(),
   	 desiredSum.val(),
   	 (origSum - desiredSum).val());
#endif
  
  unsigned k;
  for (k=0;k<curNumCliqueValuesUsed;k++) {

    actualSum += curCliqueVals[k].p.pow(exponentiate); // /loc_maxCEValue;

    // printf("k=%d: origSum = %.16e, desiredSum = %.16e, actualSum = %.16e\n",k,origSum.valref(),desiredSum.valref(),actualSum.valref());

    // could use either ">" or ">=" here.
    //   - use ">=" if you want more aggressive pruning (will probably
    //       want to use a bigger minSize in this case).
    //   - use ">" if you want less agressive pruning.
    // This can actually have a big effect since due to the dynamic
    // range of the elements in the clique, the rest of the clique
    // when added to the current sum might not cause any increment at
    // all (all it needs to be is different by less than the min
    // difference, see logp.h). If we were to use ">" here, then
    // we would continue to add all the rest of the clique while
    // not changing actualSum, so no pruning would be done. We
    // therefore use ">=".
    if (actualSum >= desiredSum) {
      // increment k so that it becomes a count rather than an index.
      k++;
      break;
    }
  }
  
  if (furtherBeam != 0.0 && k < curNumCliqueValuesUsed ) {
    logpr curMax = curCliqueVals[k].p;
    logpr threshold = curMax/logpr((void*)0,furtherBeam);
    while (++k < curNumCliqueValuesUsed) {
      if (curCliqueVals[k].p  < threshold)
	break;
    }
  }

  //   // optional code to calculate and print residual (stuff that is pruned away)
  //   unsigned r=k;
  //   logpr residualSum;
  //   for (r=k+1;r<curNumCliqueValuesUsed;r++) {
  //     residualSum += curCliqueVals[r].p; // /loc_maxCEValue;
  //   }
  //   printf("mass pruning: removeFraction %.18e, maxCV = %.18e, actualSum = %.18e, residualSum = %0.18e, actual+residual = %0.18e\n",
  // 	 removeFraction,loc_maxCEValue.val(),actualSum.val(),residualSum.val(),
  // 	 (actualSum+residualSum).val());
  
  
  if (k<minSize)
    k=min(minSize,curNumCliqueValuesUsed);

  const unsigned newStateSpace = min(k,curNumCliqueValuesUsed);
  infoMsg(IM::Med,"Clique mass-beam pruning: Original state space = %d (exp-mass=%e), new state space = %d, desired exp-mass = %e, actual exp-mass = %e\n",
	  curNumCliqueValuesUsed,
	  origSum.val(),
	  newStateSpace,
	  desiredSum.val(),
	  actualSum.val());
  return newStateSpace;

}


/*-
 *-----------------------------------------------------------------------
 * MaxCliqueTable::ceCliqueUniformSamplePrunedCliquePortion()
 *
 *    After pruning has occured, this routine will uniformly sample a
 *    certain portion of the entries that have been pruned away, and
 *    place them back in the clique. It does this uniformly, i.e.,
 *    irrespectively of the mass/score on each clique entry. The
 *    reason for this is that if pruning by score (which is what all
 *    of the pruning methods so far do) leads us to get a zero
 *    probability, then we add a bit of noise irrespective of the
 *    probabilities back in to the inference so that hopefully one of
 *    the entries will allow us to get a non-zero score at the very
 *    end of inference.
 *  
 *    TODO: select entries rather than uniformly at random, but based also on
 *          some form of "diversity", so that for the number of re-samples, we
 *          get as different a set of entries as possible.
 *    
 *
 * Preconditions:
 *   1) clique table must be created, meaning that either:
 *
 *        MaxCliqueTable::ceIterateAssignedNodesRecurse()
 *
 *   We assume that the clique has been pruned in the past, but that
 *   the memory for the clique still exists in the clique table (cliqueValues), meaning
 *   that no memory reallocation has yet been done based on prunning. We also
 *   assume that the argument 'origNumCliqueValuesUsed' contains the size of the
 *   clique before any pruning has been done, and that 'numCliqueValuesUsed' gives
 *   the number of current clique values used, and that all entries between entry
 *   numCliqueValuesUsed and origNumCliqueValuesUsed-1 inclusive are previously pruned
 *   cliques entries.
 *   
 *
 * Postconditions:
 *    Clique table has been "un-pruned", and memory for it has NOT been re-allocated to
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
MaxCliqueTable::ceCliqueUniformSamplePrunedCliquePortion(MaxClique& origin,
							 const unsigned origNumCliqueValuesUsed)
{

  const unsigned numCliqueValuesUsedBeforeSampling = numCliqueValuesUsed;

  if (origin.cliqueBeamUniformSampleAmount == 0.0)
    return;
  else if (origin.cliqueBeamUniformSampleAmount == 1.0) {
    numCliqueValuesUsed = origNumCliqueValuesUsed;
  } else {

    unsigned numEntriesPruned = origNumCliqueValuesUsed - numCliqueValuesUsed;
    if (numEntriesPruned != 0) {
      unsigned numToSample;
      if (origin.cliqueBeamUniformSampleAmount < 1.0) {
	numToSample = (unsigned)(origin.cliqueBeamUniformSampleAmount*(double)numEntriesPruned);
      } else { // > 1.0
	numToSample = (unsigned)origin.cliqueBeamUniformSampleAmount;
      }

      numToSample = min(numToSample,numEntriesPruned);
      if (numToSample == 0) {
	; // do nothing
      } else if (numToSample == numEntriesPruned) {
	numCliqueValuesUsed = origNumCliqueValuesUsed;
      } else {
	while (numToSample > 0) {

	  unsigned entry = rnd.uniform(--numEntriesPruned);

	  // swap the entry to the end of the current clique.
	  swap(cliqueValues[numCliqueValuesUsed],cliqueValues[numCliqueValuesUsed + entry]);
    
	  numToSample --;
	  numCliqueValuesUsed++;
	}
      }
    }
  }

  infoMsg(IM::Med,"Clique uniform sampling: Upped state space from %d to %d, before pruning state space was %d\n",
	  numCliqueValuesUsedBeforeSampling,numCliqueValuesUsed,origNumCliqueValuesUsed);

}




/*-
 *-----------------------------------------------------------------------
 * MaxCliqueTable::sumProbabilities()
 *
 *    Simply sum up the probabilities of all elements in this clique
 *    and return the results.
 *
 * Preconditions:
 *      Clique data structures must be created.
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      None
 *
 * Results:
 *     sum of probabilities of all elements in the clique.
 *
 *-----------------------------------------------------------------------
 */
logpr
MaxCliqueTable::
sumProbabilities()
{
  logpr p;
  if (numCliqueValuesUsed > 0) {
    // We directly assign first one rather than adding to initialized
    // zero so that logpr's log(0) floating point value is preserved
    // (see logp.h).
    p = cliqueValues.ptr[0].p;
    for (unsigned i=1;i<numCliqueValuesUsed;i++)
      p += cliqueValues.ptr[i].p;
  }
  return p;
}



/*-
 *-----------------------------------------------------------------------
 * MaxCliqueTable::sumExponentiatedProbabilities()
 *
 *    Simply sum up the exponentiated probabilities of all elements in this clique
 *    and return the results.
 *
 * Preconditions:
 *      Clique data structures must be created.
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      None
 *
 * Results:
 *     sum of probabilities of all elements in the clique.
 *
 *-----------------------------------------------------------------------
 */
logpr
MaxCliqueTable::
sumExponentiatedProbabilities(double exponent)
{
  logpr p;
  if (numCliqueValuesUsed > 0) {
    // We directly assign first one rather than adding to initialized
    // zero so that logpr's log(0) floating point value is preserved.
    p = cliqueValues.ptr[0].p.pow(exponent);
    for (unsigned i=1;i<numCliqueValuesUsed;i++)
      p += cliqueValues.ptr[i].p.pow(exponent);
  }
  return p;
}



/*-
 *-----------------------------------------------------------------------
 * MaxCliqueTable::cliqueEntropy()
 *
 *    Compute the clique entropy (base 2) of the forcibly normalized clique.
 *
 * Preconditions:
 *      Clique data structures must be created.
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      None
 *
 * Results:
 *     sum of probabilities of all elements in the clique.
 *
 *-----------------------------------------------------------------------
 */
double
MaxCliqueTable::
cliqueEntropy()
{
  logpr sum = sumProbabilities();
  double H = 0.0;
  if (numCliqueValuesUsed > 0) {
    logpr tmp = cliqueValues.ptr[0].p/sum;
    H = tmp.unlog() * tmp.val();
    for (unsigned i=1;i<numCliqueValuesUsed;i++) {
      logpr tmp = cliqueValues.ptr[i].p/sum;
      H += tmp.unlog() * tmp.val();
    }
  }
  // convert to entropy and log base 2.
  return - H / logpr::internal_log(2.0);
}




/*-
 *-----------------------------------------------------------------------
 * MaxCliqueTable::cliqueDiversity()
 *
 *    Compute the clique diversity of the clique.
 *
 * Preconditions:
 *      Clique data structures must be created.
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      None
 *
 * Results:
 *     sum of probabilities of all elements in the clique.
 *
 *-----------------------------------------------------------------------
 */
double
MaxCliqueTable::
cliqueDiversity(MaxClique& origin)
{
  double D = 0.0;
  if (numCliqueValuesUsed > 0) {
    for (unsigned i=1;i<numCliqueValuesUsed;i++) {
      for (unsigned j=i;j<numCliqueValuesUsed;j++) {
	unsigned* i_key_p; 
	unsigned* j_key_p; 
#ifdef USE_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL
	i_key_p =
	  &origin.temporaryCliqueValuePool.ptr[cliqueValues.ptr[i].ival];
	j_key_p =
	  &origin.temporaryCliqueValuePool.ptr[cliqueValues.ptr[j].ival];
#else
	i_key_p = (unsigned*)cliqueValues.ptr[i].ptr;
	j_key_p = (unsigned*)cliqueValues.ptr[j].ptr;
#endif	
	double cur_dist = 
	  origin.packer.use_distance(i_key_p,j_key_p);
	D += cur_dist;
      }
    }
  }
  D = D /(  ((double)numCliqueValuesUsed)* ((double)(numCliqueValuesUsed+1))/2.0 );
  return D;
}




/*-
 *-----------------------------------------------------------------------
 * MaxCliqueTable::reportMemoryUsageTo()
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
MaxCliqueTable::
reportMemoryUsageTo(MaxClique& origin,FILE *f)
{
  // Memory: IC=Inference Clique, CT=Clique Table
  fprintf(f,"*MEM:IC CT(used=%lu,all=%lu=%luMB), ",
	  (unsigned long)numCliqueValuesUsed,
	  (unsigned long)cliqueValues.size(),
	  (unsigned long)((1+(sizeof(CliqueValue)*(unsigned long)cliqueValues.size())/(1024ul*1024ul))));
  origin.reportMemoryUsageTo(f);
  fprintf(f,"\n");
}




/*-
 *-----------------------------------------------------------------------
 * MaxCliqueTable::maxProbability()
 *
 *    This routine just computes the max clique probabiltiy and
 *    optionally sets the clique to the clique value to the one
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
MaxCliqueTable::
maxProbability(MaxCliqueTable::SharedLocalStructure& sharedStructure,
	       bool setCliqueToMaxValue)
{
  // syntactic convenience variables.
  MaxClique& origin = *(sharedStructure.origin);

  // check for empty clique and if so, return zero.
  if (numCliqueValuesUsed == 0)
    return logpr();
  
  if (origin.hashableNodes.size() == 0) {
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

    if (setCliqueToMaxValue) {
      // store the current table entry for the max clique.
      back_max_cvn = max_cvn;

      const bool imc_nwwoh_p = (origin.packer.packedLen() <= IMC_NWWOH); 
      if (imc_nwwoh_p) {
	origin.packer.unpack((unsigned*)&(cliqueValues.ptr[max_cvn].val[0]),
			     (unsigned**)sharedStructure.discreteValuePtrs.ptr);
      } else {
	origin.packer.unpack((unsigned*)cliqueValues.ptr[max_cvn].ptr,
			     (unsigned**)sharedStructure.discreteValuePtrs.ptr);
      }
      for (unsigned j=0;j<sharedStructure.fDeterminableNodes.size();j++) {
	RV* rv = sharedStructure.fDeterminableNodes.ptr[j];
	RV2DRV(rv)->assignDeterministicChild();
      }

    }

    return max_cvn_score;
  }

}
logpr
MaxCliqueTable::
maxProb()
{
  // check for empty clique and if so, return zero.
  if (numCliqueValuesUsed == 0)
    return logpr();
  logpr mx = cliqueValues.ptr[0].p;
  // find the max score clique entry
  for (unsigned cvn=1;cvn<numCliqueValuesUsed;cvn++) {
    if (cliqueValues.ptr[cvn].p > mx) {
      mx = cliqueValues.ptr[cvn].p;
    }
  }
  return mx;
}


/*-
 *-----------------------------------------------------------------------
 * MaxCliqueTable::printCliqueEntries()
 *
 *    Simply prints out all elements in the clique, giving both the RV
 *    values and the (not necessarily a probability) score for the
 *    corresponding clique entry. If the normalize option is given,
 *    then the scores can rightfully be interpreted as probabilities,
 *    and if CE/DE has been called, then this will produce true marginal
 *    probabilities over the variables in the clique.
 *
 * Preconditions:
 *      Clique data structures must be created.
 *
 * Postconditions:
 *      none
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
MaxCliqueTable::
printCliqueEntries(MaxCliqueTable::SharedLocalStructure& sharedStructure,
		   FILE *f,const char*str, 
		   const bool normalize,
		   const bool justPrintEntropy)
{
  MaxClique& origin = *(sharedStructure.origin);

  fprintf(f,"--------\n");
  if (str != NULL)
    fprintf(f,"%s ",str);
  
  // TODO: also have a cliqueDiversity routine that can compute the
  // variability/diversity in a clique in terms of its random
  // variables.

  fprintf(f,"Printing Clique with %d variables, %d entries, H=%e\n",
	  sharedStructure.fNodes.size(),numCliqueValuesUsed,cliqueEntropy());
  if (justPrintEntropy)
    return;
  logpr sum;
  if (normalize)
    sum = sumProbabilities();
  const bool imc_nwwoh_p = (origin.packer.packedLen() <= IMC_NWWOH);
  const bool clique_has_hidden_vars = (origin.hashableNodes.size() > 0);
  for (unsigned cvn=0;cvn<numCliqueValuesUsed;cvn++) {

    if (clique_has_hidden_vars) {
      if (imc_nwwoh_p) {
	origin.packer.unpack((unsigned*)&(cliqueValues.ptr[cvn].val[0]),
			     (unsigned**)sharedStructure.discreteValuePtrs.ptr);
      } else {

#ifdef USE_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL
	origin.packer.unpack((unsigned*)&origin.temporaryCliqueValuePool.ptr[cliqueValues.ptr[cvn].ival],
			     (unsigned**)sharedStructure.discreteValuePtrs.ptr);
#else
	origin.packer.unpack((unsigned*)cliqueValues.ptr[cvn].ptr,
			     (unsigned**)sharedStructure.discreteValuePtrs.ptr);
#endif
      }
      for (unsigned j=0;j<sharedStructure.fDeterminableNodes.size();j++) {
	RV* rv = sharedStructure.fDeterminableNodes.ptr[j];
	RV2DRV(rv)->assignDeterministicChild();
      }
    }
    if (normalize) {
      // then print the exponentiated probability
      fprintf(f,"%d: %.8e ",cvn,(cliqueValues.ptr[cvn].p/sum).unlog());
    } else {
      // print the log value directly
      fprintf(f,"%d: %f ",cvn,cliqueValues.ptr[cvn].p.valref());
    }
    printRVSetAndValues(f,sharedStructure.fNodes);
  }
}



/*-
 *-----------------------------------------------------------------------
 * MaxCliqueTable::emIncrement()
 *
 *      Do the work of EM increment. That is, for each clique element,
 *      load the values into the variables and then for each assigned prob node
 *      call that nodes emIncrement() routine thus realizing EM training.
 *      The routine has the option of using either a local normalization
 *      (useful when doing double precision 64-bit floating point) or a local
 *      normalization (perhaps better for signle precision, especially on
 *      longer utterances).
 *
 * Preconditions:
 *      Clique data structures must be created.
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      Potentially all assigned probability nodes accumulators are changed.
 *
 * Results:
 *     none
 *
 *-----------------------------------------------------------------------
 */
void
MaxCliqueTable::
emIncrement(MaxCliqueTable::SharedLocalStructure& sharedStructure,
	    const logpr probE,
	    const bool localCliqueNormalization,
	    const double emTrainingBeam)
{
  MaxClique& origin = *(sharedStructure.origin);

  // recompute here each time, shouldn't take too long
  // TODO: re-compute this once for each inference clique.
  sArray< RV*> fAssignedProbNodes;
  unsigned numAssignedProbNodes = 0;
  for (unsigned nodeNumber = 0; 
       nodeNumber < sharedStructure.fSortedAssignedNodes.size(); 
       nodeNumber ++ ) {
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
  for (unsigned nodeNumber = 0; nodeNumber < 
	 sharedStructure.fSortedAssignedNodes.size(); nodeNumber ++ ) {
    if (origin.dispositionSortedAssignedNodes[nodeNumber] == MaxClique::AN_CPT_ITERATION_COMPUTE_AND_APPLY_PROB ||
	origin.dispositionSortedAssignedNodes[nodeNumber] == MaxClique::AN_COMPUTE_AND_APPLY_PROB) {
      RV* rv = sharedStructure.fSortedAssignedNodes[nodeNumber];
      fAssignedProbNodes[numAssignedProbNodes++] = rv;
    }
  }

  logpr locProbE((void*)0);
  if (localCliqueNormalization) {
    // this case is probably better/safer when using 32-bit single precision IEEE floating point for logpr.
    locProbE = sumProbabilities();
  } else {
    // probably ok to do this when using double precision IEEE floating point.
    locProbE = probE;
  }

  const bool imc_nwwoh_p = (origin.packer.packedLen() <= IMC_NWWOH);
  const bool clique_has_hidden_vars = (origin.hashableNodes.size() > 0);
  
  logpr beamThreshold((void*)0);
  if (emTrainingBeam != (-LZERO)) {
    // then we do clique table pruning right here rather
    // than a separate call to ceCliqueBeamPrune().
    // break into the logp to avoid unnecessary zero checking.
    beamThreshold.valref() = - emTrainingBeam;
  } else {
    // set beam threshold to a value that will never cause real
    // pruning (i.e., we always prune posteriors that are almost
    // zero).
    beamThreshold.set_to_almost_zero();
  }
  // create unnormalized beam threshold.
  beamThreshold *= locProbE;

  // now go through updating each thing
  for (unsigned cvn=0;cvn<numCliqueValuesUsed;cvn++) {

    // EM pruning here based on unnormalized posterior. Don't bother
    // with things that are below threshold.
    if (cliqueValues.ptr[cvn].p <= beamThreshold)
      continue;

    // if still here, then create the posterior to update the
    // parameters.
    logpr posterior = cliqueValues.ptr[cvn].p/locProbE;

    // printf("EM training, cvn=%d, log(posterior) = %f\n",cvn,posterior.valref());

    if (clique_has_hidden_vars) {
      if (imc_nwwoh_p) {
	origin.packer.unpack((unsigned*)&(cliqueValues.ptr[cvn].val[0]),
			     (unsigned**)sharedStructure.discreteValuePtrs.ptr);
      } else {
	origin.packer.unpack((unsigned*)cliqueValues.ptr[cvn].ptr,
			     (unsigned**)sharedStructure.discreteValuePtrs.ptr);
      }
      for (unsigned j=0;j<sharedStructure.fDeterminableNodes.size();j++) {
	RV* rv = sharedStructure.fDeterminableNodes.ptr[j];
	RV2DRV(rv)->assignDeterministicChild();
      }

    }
    // Increment all assigned probability nodes.  
    // 
    // TODO: integerate out all but cont. variables parents so we
    // don't multiply increment those varialbes for the same parent
    // values (waisting time). Alternatively, this can be done by
    // creating a sub-clique hanging off of this clique which contains
    // only the observed variables and its parents.
    for (unsigned nodeNumber = 0; nodeNumber < fAssignedProbNodes.size(); nodeNumber ++ ) {
      RV* rv = fAssignedProbNodes[nodeNumber];
      rv->emIncrement(posterior);
    }
  }

}



/*-
 *-----------------------------------------------------------------------
 * MaxCliqueTable::deReceiveFromIncommingSeparator()
 *
 *      We are in the DE phase, and the separator that during the CE phase
 *      we sent a message out to is now ready with a backwards message back
 *      to this clique. We iterate through all clique entries and lookup
 *      the corresponding "incomming separator" entry, multiplying in its score
 *      into the corresponding clique entry.
 *
 *      Note that some of the incomming separator values might be zero (for
 *      reason based on higher up in the JT). Therefore, at the end of this
 *      routine we do another DE pruning pass to remove clique entries
 *      that have become zero. Since one separator entry might correspond to
 *      many clique entries, a zero separator entry might cause a lot of backwards
 *      clique pruning. 
 *
 *      Note: A condition might occur where we've got a clique entry without
 *      a corresponding incomming separator entry. On first thought, this
 *      shouldn't occur, but in actuallity it can occur for one of two reasons:
 *
 *     1) if sbeam is turned on, it's possible that the separator entry corresponding to this
 *        clique entry was pruned away. This can happen since separator pruning is later
 *        then the clique which created (projected down into) the separator. The solution
 *        here is to increase or entirely turn off sbeam pruning.
 *     
 *     2) During the island algorithm and using k-pruning (fixed k clique size), there
 *        were ties in the clique probability and since we use a random median to find the top
 *        k entries, we didn't get the same top k entries the 2nd (or 3rd, 4th, etc. depending 
 *        on the island algorithm's 'base' and 'lst' parameters) time that we created the 
 *        clique and then pruned it down to k entries. The solution we employ is to
 *        seed the random number generator in k-pruning to something that is clique dependent.
 *        This is now being done in MaxCliqueTable::ceCliqueStatePrune() above. Search for string 
 *        K-BEAM-SEED elsewhere in this file for further information. Therefore, this 
 *        case should not occur.
 *
 *        The solution (Marked by the keyword "SEPCLIQUEZERO" below) is to just
 *        prune this clique entry away, for case 1 since case 2 will not
 *        happen. 
 *
 *        Another solution would be to remove the clique-dependent seed for case
 *        2, but pruning would then get ugly. This is because if the clique was
 *        all ties, this would mean that we prune all but the intersection of
 *        the survivors of both k-prunings, and this might be small or empty.
 *        In the best of cases, this would probably be much smaller than k, so
 *        we instead just do the clique-dependent seed which is much easier.
 *
 * Preconditions:
 *      Clique data structures and separator must be created.
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      Potentially all assigned probability nodes accumulators are changed.
 *
 * Results:
 *     none
 *
 *-----------------------------------------------------------------------
 */
void 
MaxCliqueTable::
deReceiveFromIncommingSeparator(MaxCliqueTable::SharedLocalStructure& sharedStructure,
				ConditionalSeparatorTable* separatorTableArray,
				ConditionalSeparatorTable::SharedLocalStructure* sepSharedStructureArray)
{
  MaxClique& origin = *(sharedStructure.origin);
  deReceiveFromIncommingSeparator(sharedStructure,
				  separatorTableArray[origin.ceSendSeparator],
				  sepSharedStructureArray[origin.ceSendSeparator]
				  );
}
void 
MaxCliqueTable::
deReceiveFromIncommingSeparator(MaxCliqueTable::SharedLocalStructure& sharedStructure,
				ConditionalSeparatorTable& sep,
				ConditionalSeparatorTable::SharedLocalStructure& sepSharedStructure)
{
  if (JunctionTree::viterbiScore) {
    return deReceiveFromIncommingSeparatorViterbi(sharedStructure,
						  sep,
						  sepSharedStructure);
  }

  // syntactic convenience variables.
  MaxClique& origin = *(sharedStructure.origin);
  SeparatorClique& sepOrigin = 
    *(sepSharedStructure.origin);  

  // keep a local variable copy of this around to avoid potential dereferencing.
  ConditionalSeparatorTable::AISeparatorValue * const
    sepSeparatorValuesPtr = sep.separatorValues->ptr; 

  if (origin.hashableNodes.size() == 0) {
    // do the observed clique case up front right here so we don't
    // need to keep checking below.
    ConditionalSeparatorTable::AISeparatorValue& sv
      = sepSeparatorValuesPtr[0];
    cliqueValues.ptr[0].p *= sv.remValues.ptr[0].bp();
    return;
  }

  // allocate some temporary storage for packed separator values.
  // 128 words is *much* bigger than any possible packed clique value
  // will take on, but it is easy/fast to allocate on the stack right now.
  // unsigned packedVal[128];
  // but just in case, we assert.
  assert ((sepOrigin.hAccumulatedIntersection.size() == 0)
	  ||
	  (sepOrigin.accPacker.packedLen() < 128)
	  );
  assert ((sepOrigin.hRemainder.size() == 0) 
	  ||
	  (sepOrigin.remPacker.packedLen() < 128 )
	  );
  // If this assertion fails (at some time in the future, probably in
  // the year 2150), then it is fine to increase 128 to something larger.

  // cache check here.
  const bool imc_nwwoh_p = (origin.packer.packedLen() <= IMC_NWWOH);
  for (unsigned cvn=0;cvn<numCliqueValuesUsed;cvn++) {{

      // TODO: optimize away this conditional check.
      if (imc_nwwoh_p) {
	origin.packer.unpack((unsigned*)&(cliqueValues.ptr[cvn].val[0]),
			     (unsigned**)sharedStructure.discreteValuePtrs.ptr);
      } else {
	origin.packer.unpack((unsigned*)cliqueValues.ptr[cvn].ptr,
			     (unsigned**)sharedStructure.discreteValuePtrs.ptr);
      }
      for (unsigned j=0;j<sharedStructure.fDeterminableNodes.size();j++) {
	RV* rv = sharedStructure.fDeterminableNodes.ptr[j];
	RV2DRV(rv)->assignDeterministicChild();
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
      if (sepOrigin.hAccumulatedIntersection.size() > 0) {
	// an accumulated intersection exists.

	sepOrigin.accPacker.pack((unsigned**)sepSharedStructure.accDiscreteValuePtrs.ptr,
				 &CliqueBuffer::packedVal[0]);
	unsigned* accIndexp =
	  sep.iAccHashMap->find(&CliqueBuffer::packedVal[0]);

	// case SEPCLIQUEZERO, see comments in routine heading
	if ( accIndexp == NULL ) {
	  // Then separator entry got pruned away. Force prune of clique
	  // entry as well.
	  cliqueValues.ptr[cvn].p.set_to_zero();
	  goto next_iteration;
	}

	accIndex = *accIndexp;

	// TODO: optimize this check out of loop.
	if (sepSharedStructure.remDiscreteValuePtrs.size() == 0) {
	  // 2) AI exists and REM doesnt exist
	  // Then this separator is entirely covered by one or 
	  // more other separators earlier in the order.

	  // go ahead and insert it here to the 1st entry (entry 0).

	  // handy reference for readability.
	  ConditionalSeparatorTable::AISeparatorValue& sv
	    = sepSeparatorValuesPtr[accIndex];

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

      if (sepSharedStructure.remDiscreteValuePtrs.size() > 0) {
	// if we're here, then we must have some remainder pointers.
	// Do the remainder exists in this separator.
	// either:
	//   1) AI exists and REM exist
	//     or
	//   3) AI does not exist (accIndex == 0), but REM exists
	// 

	// keep handy reference for readability.
	ConditionalSeparatorTable::AISeparatorValue& sv
	  = sepSeparatorValuesPtr[accIndex];

	sepOrigin.remPacker.pack((unsigned**)sepSharedStructure.remDiscreteValuePtrs.ptr,
				 &CliqueBuffer::packedVal[0]);

	unsigned* remIndexp =
	  sv.iRemHashMap.find(&CliqueBuffer::packedVal[0]);

	// case SEPCLIQUEZERO, see comments in routine heading
	if ( remIndexp == NULL ) {
	  // Then separator entry got pruned away. Force prune of clique
	  // entry as well.
	  cliqueValues.ptr[cvn].p.set_to_zero();
	  goto next_iteration;
	  // warning("ERROR: During distribute evidence, found clique entry without corresponding incomming separator entry. Try increasing -sbeam beam, or just use -cbeam without the -sbeam option.\nClique contains:");
	  // printRVSetAndValues(stdout,fNodes);

	}

	// We've finally got the sep entry. Multiply it it into the
	// current clique value.
	cliqueValues.ptr[cvn].p *= sv.remValues.ptr[*remIndexp].bp();
      } else {
	// Either separator is all observed, or the separator
	// is completely contained in the accumulated intersection.
	// In either case, we multiply in its one value.

	ConditionalSeparatorTable::AISeparatorValue& sv
	  = sepSeparatorValuesPtr[accIndex];

	if (sv.numRemValuesUsed == 1) {
	  // We've finally got the sep entry. Multiply it it into the
	  // current clique value.
	  cliqueValues.ptr[cvn].p *= sv.remValues.ptr[0].bp();
	} else {
	  // case SEPCLIQUEZERO, see comments in routine heading
	  // Then separator entry got pruned away. Force prune of clique
	  // entry as well.
	  cliqueValues.ptr[cvn].p.set_to_zero();
	  goto next_iteration;
	}
      }

    }
  next_iteration:
    ;    
  }

  // Backwards pruning: Fixed backwards/distribute evidence beam
  // pruning here. Since we know log(prob(E)), we can do fairly
  // accurate pruning now. This will be useful particulalry when many
  // of the bp() values are zero (i.e., in this case, we have zero
  // compression). We just prune out the zeros for now.  
  // TODO: integrate this pruning into the above loop instead.  
  // TODO: export backwards beam width to command line & integrate
  //    with -ebeam
  {
    const unsigned origNumCliqueValuesUsed = numCliqueValuesUsed;
    for (unsigned cvn=0;cvn<numCliqueValuesUsed;) {
      if (cliqueValues.ptr[cvn].p.essentially_zero()) {

	// copy the last entry (that might be good) to the current
	// position (which has just been pruned). Do an assignment
	// since the pruned entries are never going to be needed again
	// (they're deallocated below).
	cliqueValues.ptr[cvn] = cliqueValues.ptr[--numCliqueValuesUsed];

	// alternatively, we could swap with last entry, and decrease
	// numCliqueValuesUsed by one (if for some reason we want to
	// possibly use these zero entries someday).
	// swap(cliqueValues.ptr[cvn],cliqueValues.ptr[numCliqueValuesUsed-1]);
	// numCliqueValuesUsed--;

      } else {
	cvn++;
      }
    }
    infoMsg(IM::High-1,"DE Clique Receive: (old,new) clique state space = (%d,%d).\n",
	    origNumCliqueValuesUsed,numCliqueValuesUsed);
    // TODO: resize only if size difference is large.
    if (numCliqueValuesUsed < origNumCliqueValuesUsed)
      cliqueValues.resizeAndCopy(numCliqueValuesUsed);
  }

}




/*-
 *-----------------------------------------------------------------------
 * MaxCliqueTable::deReceiveFromIncommingSeparatorViterbi()
 *
 *      We are in the DE phase, and the separator that during the CE
 *      phase we sent a message out to is now ready with a backwards
 *      message back to this clique. 
 * 
 *      This routine here is the Viterbi version of
 *      deReceiveFromIncommingSeparator() above.  In this case, we
 *      simply look up the separators back pointer index and use it to
 *      choose the current clique entry. Moreover, we unpack that
 *      clique entry leaving the RVs assigned to what is given in the
 *      clique entry.
 *
 *      An alternative strategy would be to assume that the RVs of the
 *      separator currently point to (i.e., hold the value of) the
 *      viterbi entry containing the back pointer index. We would then
 *      look up the separator entry corresponding to the currently set
 *      RV values via hashing, and choose the current clique entry
 *      accordingly. This, however, interacts poorly with the island
 *      algorithm, and would not generalize to n-best, so we do not do
 *      this here.
 *
 *      Yet another strategy would be to look up the separator
 *      backpointer index and use it to choose a clique entry and then
 *      store the index of that clique entry in the clique (rather
 *      than using the random variable values) since this would
 *      greatly simplify things like the island algorithm and shifting
 *      (but would require a bit more memory).
 *
 *
 * Preconditions:
 *      Clique data structures and separator must be created.
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      Potentially all assigned probability nodes accumulators are changed.
 *
 * Results:
 *     none
 *
 *-----------------------------------------------------------------------
 */
void 
MaxCliqueTable::
deReceiveFromIncommingSeparatorViterbi(MaxCliqueTable::SharedLocalStructure& sharedStructure,
				       ConditionalSeparatorTable& sep,
				       ConditionalSeparatorTable::SharedLocalStructure& sepSharedStructure)
{
  MaxClique& origin = *(sharedStructure.origin);

  if (origin.hashableNodes.size() == 0) {
    // All clique values already set to their max (and only) settings,
    // so no need to do anything.
    return;
  }

  // keep a local variable copy of this around to avoid potential dereferencing.
  ConditionalSeparatorTable::AISeparatorValue * const
    sepSeparatorValuesPtr = sep.separatorValues->ptr; 

  // cache check here.
  const bool imc_nwwoh_p = (origin.packer.packedLen() <= IMC_NWWOH);
  // grab backpointer from seps forward pointer entry directly.
  ConditionalSeparatorTable::AISeparatorValue& sv
    = sepSeparatorValuesPtr[sep.forwPointer.viterbiAccIndex];
  unsigned cvn = sv.remValues.ptr[sep.forwPointer.viterbiRemIndex].backPointer;

  // store the current table entry for the max clique.
  back_max_cvn = cvn;

  // unpack clique value 'cvn' into corresponding random variables and expand
  // any deterministic values.
  if (imc_nwwoh_p) {
    origin.packer.unpack((unsigned*)&(cliqueValues.ptr[cvn].val[0]),
			 (unsigned**)sharedStructure.discreteValuePtrs.ptr);
  } else {
    origin.packer.unpack((unsigned*)cliqueValues.ptr[cvn].ptr,
			 (unsigned**)sharedStructure.discreteValuePtrs.ptr);
  }
  for (unsigned j=0;j<sharedStructure.fDeterminableNodes.size();j++) {
    RV* rv = sharedStructure.fDeterminableNodes.ptr[j];
    RV2DRV(rv)->assignDeterministicChild();
  }

  // printf("*** max RV of clique values set from parent separator ***\n");
  // printRVSetAndValues(stdout,sharedStructure.fNodes,true);
  // fflush(stdout);
}



/*-
 *-----------------------------------------------------------------------
 * MaxCliqueTable::deScatterToOutgoingSeparators()
 *
 *      We are in the DE phase, and we've received the separator
 *      message fror the separator we sent a message out during CE.
 *      Now it is our turn to send (project down to) to the separators 
 *      that we received incomming messages during the CE phase.
 *
 *      Specifically, we "scatter" out to the now outgoing separators which are the same as the
 *      "incomming" separators in the collect evidence stage, so we use that
 *      array here directly here. We do this by iterating through all clique entries
 *      and sending it back out to all CE-incomming separators (unless the separator
 *      is a VE separator in which case it is just a constaint that
 *      we have already acounted for).
 *
 * Preconditions:
 *      Clique data structures and separator must be created.
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      Potentially all assigned probability nodes accumulators are changed.
 *
 * Results:
 *     none
 *
 *-----------------------------------------------------------------------
 */
void 
MaxCliqueTable::
deScatterToOutgoingSeparators(MaxCliqueTable::SharedLocalStructure& sharedStructure,
			      ConditionalSeparatorTable* separatorTableArray,
			      ConditionalSeparatorTable::SharedLocalStructure* sepSharedStructureArray)
{

  if (JunctionTree::viterbiScore) {
    // While we might think that there is nothing to do in this case
    // since if the RVs associated with the current clique have been
    // set to the appropriate clique table entry, the associated
    // separator RVs have also been set. But in some cases (namely the
    // island algorithm) it is possible for the separator to get
    // changed by an additional forward pass in which case we would be
    // using the wrong value.
    // 
    // Specifically, to get island decoding working, we need to have
    // separator keep track of which entry is current max rather than
    // having it assume that its values are set from the 
    // clique that, on CE, gathers in that separator and, on DE,
    // scatters out.
    // 
    // Consider the following example:
    //  
    //    The vertical bars mark the partition boundaries 
    //    of partitions P0, P1, ...
    //
    //        |           |           |   island   |  ...              
    //  P0    |     P1    |     P2    |     P3     |  ...              
    //        |           |           |            |  ...
    //    C0 -- s0 -- C1 -- s1 -- C2 -- s2 -- C3 --
    //             1>    2>    3>    <4          <b  
    //                         <5  
    // 
    //  P3 is an island in an island algorithm. We have gone through
    //  and created and deleted partitions P0, P1, P2, and P3, but we
    //  store and save P3 as an island.  We then move to the right of
    //  P3, and come back, but we then need to reconstruct partitions
    //  P0 - P3. In other words:
    // 
    //  1> and 3> are ceGatherIntoRoot msgs
    //  2> is a ceSendForwardsCrossPartitions msg
    //  <4 is a deReceveFromIncommingSep msg
    //      and <5 is a deScatterOutofRoot, but has not yet happened.
    //  <b is a deReceveFromIncommingSep msg from before (thus
    //     marked with a 'b').
    // 
    // The potential problem is that <4 might assume that it's 's2' is currently
    // set to appropriate values (max of C3 done by <b), but message 3>
    // changes clique C2 and thus changes s2 so that it no longer
    // necessarily holds max value, since s2 holds whever is 
    // currently set by 3>.
    // 
    // One solution (that adds memory) is to add two variables in
    // InferenceSep, i.e., two indices to get at the Viterbi values
    // for current separator.
    //
    //      unsigned viterbiAccIndex;
    //      unsigned viterbiRemIndex;
    //
    // to keep track of which entries in seps (such as s2) are max.
    // Then <4 would restore those values for each separator.  and
    // deScatterOutofRoot would compute them (since deScatterOutofRoot
    // is called at the time C3 is at its true max value). Note that
    // this works between partitions since the separator between
    // two partitions is contained in the right partition.
    //
    // Note also that this can probably be easily extended to the
    // N-best list case by having multiple index pairs.
    deScatterToOutgoingSeparatorsViterbi(sharedStructure,separatorTableArray,sepSharedStructureArray);
    return; 
  }
  MaxClique& origin = *(sharedStructure.origin);

  if (origin.ceReceiveSeparators.size() == 0)
    return;


  // Note. All separator .bp values have already been initialized to
  // zero when the structure containing them 'a RemainderValue' was
  // constructed. All memory reallocations will have preserved these
  // initializations, so there is no need to scan through initializing
  // bp to zero here.

  if (origin.hashableNodes.size() == 0) {
    // printf("<<<<<==== Observed clique in backwards pass\n");

    // Do the observed clique case up front right here so we don't
    // need to keep checking below. Here, the clique is observed which
    // means that all connecting separators are also observed. We just
    // sweep through all separators updating the values.
    for (unsigned sepNumber=0;sepNumber<origin.ceReceiveSeparators.size();sepNumber++) {
      ConditionalSeparatorTable& sep = 
	separatorTableArray[origin.ceReceiveSeparators[sepNumber]];
      ConditionalSeparatorTable::AISeparatorValue& sv
	= sep.separatorValues->ptr[0];
      // can use assignment rather than += here since there is only one value.
      sv.remValues.ptr[0].bp() = cliqueValues.ptr[0].p;      
      sv.numRemValuesUsed = 1;
    }
  } else {

    // Allocate some temporary storage for packed separator values.
    // 128 words is *much* bigger than any possible packed clique value
    // will take on, but it is easy/fast to allocate on the stack right now.
    // unsigned packedVal[128];

    infoMsg(IM::High-1,"DE Clique state space = %d.\n",numCliqueValuesUsed);

    // cache check here.
    const bool imc_nwwoh_p = (origin.packer.packedLen() <= IMC_NWWOH); 
    for (unsigned cvn=0;cvn<numCliqueValuesUsed;cvn++) {

      // TODO: optimize away this conditional check.
      if (imc_nwwoh_p) {
	origin.packer.unpack((unsigned*)&(cliqueValues.ptr[cvn].val[0]),
			     (unsigned**)sharedStructure.discreteValuePtrs.ptr);
      } else {
	origin.packer.unpack((unsigned*)cliqueValues.ptr[cvn].ptr,
			     (unsigned**)sharedStructure.discreteValuePtrs.ptr);
      }
      for (unsigned j=0;j<sharedStructure.fDeterminableNodes.size();j++) {
	RV* rv = sharedStructure.fDeterminableNodes.ptr[j];
	RV2DRV(rv)->assignDeterministicChild();
      }

      // now we iterate through all the separators.
      for (unsigned sepNumber=0;sepNumber<origin.ceReceiveSeparators.size();sepNumber++) {

	// get a handy reference to the current separator
	ConditionalSeparatorTable& sep = 
	  separatorTableArray[origin.ceReceiveSeparators[sepNumber]];
	ConditionalSeparatorTable::SharedLocalStructure& sepSharedStructure = 
	  sepSharedStructureArray[origin.ceReceiveSeparators[sepNumber]];
	SeparatorClique& sepOrigin = 
	  *(sepSharedStructure.origin);

	// don't distribute to VE separators or to one that is being skipped.
	if (sep.veSeparator() || sepOrigin.skipMe)
	  continue;

	// keep a local variable copy of this around to avoid potential dereferencing.
	ConditionalSeparatorTable::AISeparatorValue * const
	  sepSeparatorValuesPtr = sep.separatorValues->ptr; 

	// If these assertions fail (at some time in the future, probably in
	// the year 2150), then it is fine to increase 128 to something larger.
	// In fact, 128 is so large, lets not even do the assert.
	// assert ( sepOrigin.accPacker.packedLen() < 128 );
	// assert ( sepOrigin.remPacker.packedLen() < 128 );


	/*
	 * There are 3 cases.
	 * 1) AI exists and REM exist
	 * 2) AI exists and REM doesnt exist
	 * 3) AI does not exist, but REM exists
	 * AI not exist and REM not exist can't occur.
	 */

	unsigned accIndex;
	// TODO: optimize this check away out of loop.
	if (sepOrigin.hAccumulatedIntersection.size() > 0) {
	  // an accumulated intersection exists.

	  sepOrigin.accPacker.pack((unsigned**)sepSharedStructure.accDiscreteValuePtrs.ptr,
				   &CliqueBuffer::packedVal[0]);
	  unsigned* accIndexp =
	    sep.iAccHashMap->find(&CliqueBuffer::packedVal[0]);

	  // we should always find something or else something is wrong.
	  assert ( accIndexp != NULL ); 
	  accIndex = *accIndexp;

	  // TODO: optimize this check out of loop.
	  if (sepSharedStructure.remDiscreteValuePtrs.size() == 0) {
	    // 2) AI exists and REM doesnt exist
	    // Then this separator is entirely covered by one or 
	    // more other separators earlier in the order.

	    // go ahead and insert it here to the 1st entry (entry 0).

	    // handy reference for readability.
	    ConditionalSeparatorTable::AISeparatorValue& sv
	      = sepSeparatorValuesPtr[accIndex];

	    // Add in this clique value's probability.  Note that bp was
	    // initialized during forward pass.
	    sv.remValues.ptr[0].bp() += cliqueValues.ptr[cvn].p;
	    // done, move on to next separator.
	    continue; 
	  } // else, we continue on below.
	} else {
	  // no accumulated intersection exists, everything
	  // is in the remainder.
	  accIndex = 0;
	}

	if (sepSharedStructure.remDiscreteValuePtrs.size() > 0) {
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
	  ConditionalSeparatorTable::AISeparatorValue& sv
	    = sepSeparatorValuesPtr[accIndex];
	
	  sepOrigin.remPacker.pack((unsigned**)sepSharedStructure.remDiscreteValuePtrs.ptr,
				   &CliqueBuffer::packedVal[0]);

	  unsigned* remIndexp =
	    sv.iRemHashMap.find(&CliqueBuffer::packedVal[0]);

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
	  ConditionalSeparatorTable::AISeparatorValue& sv
	    = sepSeparatorValuesPtr[accIndex];

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
    ConditionalSeparatorTable& sep = 
      separatorTableArray[origin.ceReceiveSeparators[sepNumber]];
    SeparatorClique& sepOrigin = 
      *(sepSharedStructureArray[origin.ceReceiveSeparators[sepNumber]].origin);


    // don't distribute to VE separators or to one that is being skipped.
    if (sep.veSeparator() || sepOrigin.skipMe)
      continue;

    // keep a local variable copy of this around to avoid potential dereferencing.
    ConditionalSeparatorTable::AISeparatorValue * const
      sepSeparatorValuesPtr = sep.separatorValues->ptr; 

    for (unsigned aiNo=0;aiNo < sep.numSeparatorValuesUsed; aiNo ++) {
      ConditionalSeparatorTable::AISeparatorValue* aisep = &(sepSeparatorValuesPtr[aiNo]);
      for (unsigned remNo=0; remNo < aisep->numRemValuesUsed; remNo++) {
	ConditionalSeparatorTable::RemainderValue* sep_entry = &(aisep->remValues.ptr[remNo]);
	// We remove p from bp since bp will already have a factor of
	// p in it. We do this by dividing it out.
	// -
	//
	// We must make sure that if CE stage is entirely zero (i.e.,
	// zero decoding), we do not run DE stage, as in that case it
	// might be the case that sep_entry->p == 0.
	//
	// -
	// In the normal case (CE != 0), we could do direct value
	// reference subtraction in log domain (corresponding to
	// divison in original domain) to ensure that compiler creates
	// no temporaries. In other words, this operation could be either:
	//
	//        sep_entry->bp = sep_entry->bp / sep_entry->p;
	// or 
	//        sep_entry->bp.valref() = sep_entry->bp.valref() - sep_entry->p.valref(); 
	// 
	// We do slower version for now until we are certain this is
	// debugged: We assume here that (!sep_entry->p.zero()) is
	// true since we pruned all zero p's above. If we didn't
	// prune, then 'sep_entry->p == zero' would imply that
	// 'sep_entry->bp == zero', and we would need to do a
	// check. Note that this pruning always occurs, regardless of
	// beam.
	sep_entry->bp() /= sep_entry->p;
      }
    }
  }
}



/*-
 *-----------------------------------------------------------------------
 * MaxCliqueTable::deScatterToOutgoingSeparatorsViterbi()
 *
 *      Viterbi version of deScatterToOutgoingSeparators(). We assume
 *      here that all RV values in the separator are currently set to
 *      appropriate max value the current clique. This is done by
 *      actually assuming that the parent clique of the separators is
 *      currently assigned to the max value.
 *
 *      What we do here is find the index (two indices actually) in
 *      each CE-incomming separator corresponding to the current max
 *      clique assignment, and store this index in the separator
 *      itself. These indices can thus at a later time be used to look
 *      up the separator entries corresponding to the currently (at
 *      the time of this routine call) assigned clique RV values, but
 *      regardless of what the RV values happen to be assigned to at
 *      that later time (and they might be different, say, during the
 *      island algorithm).
 *
 * Preconditions:
 *      Clique data structures and separator must be created. RVs within clique 
 *      are assumed to be set to their maximum (or DE Viterbi max) values.
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      changes all CE-incomming separators.
 *
 * Results:
 *     none
 *
 *-----------------------------------------------------------------------
 */
void 
MaxCliqueTable::
deScatterToOutgoingSeparatorsViterbi(MaxCliqueTable::SharedLocalStructure& sharedStructure,
				     ConditionalSeparatorTable* separatorTableArray,
				     ConditionalSeparatorTable::SharedLocalStructure* sepSharedStructureArray)
{

  // allocate some temporary storage for packed separator values.
  // 128 words is *much* bigger than any possible packed clique value
  // will take on, but it is easy/fast to allocate on the stack right now.
  // unsigned packedVal[128];

  MaxClique& origin = *(sharedStructure.origin);

  // restore parent clique to the cvn stored in the clique table.

  // unpack clique value 'back_max_cvn' into corresponding random variables and expand
  // any deterministic values.
  const bool imc_nwwoh_p = (origin.packer.packedLen() <= IMC_NWWOH); 
  if (imc_nwwoh_p) {
    origin.packer.unpack((unsigned*)&(cliqueValues.ptr[back_max_cvn].val[0]),
			 (unsigned**)sharedStructure.discreteValuePtrs.ptr);
  } else {
    origin.packer.unpack((unsigned*)cliqueValues.ptr[back_max_cvn].ptr,
			 (unsigned**)sharedStructure.discreteValuePtrs.ptr);
  }
  for (unsigned j=0;j<sharedStructure.fDeterminableNodes.size();j++) {
    RV* rv = sharedStructure.fDeterminableNodes.ptr[j];
    RV2DRV(rv)->assignDeterministicChild();
  }

  // now we iterate through all the separators.
  for (unsigned sepNumber=0;sepNumber<origin.ceReceiveSeparators.size();sepNumber++) {
    // get a handy reference to the current separator
    ConditionalSeparatorTable& sep = 
      separatorTableArray[origin.ceReceiveSeparators[sepNumber]];
    ConditionalSeparatorTable::SharedLocalStructure& sepSharedStructure = 
      sepSharedStructureArray[origin.ceReceiveSeparators[sepNumber]];
    SeparatorClique& sepOrigin = 
      *(sepSharedStructure.origin);

    // don't distribute to VE separators or to one that is being skipped.
    if (sepOrigin.veSeparator || sepOrigin.skipMe)
      continue;

    // keep a local variable copy of this around to avoid potential dereferencing.
    ConditionalSeparatorTable::AISeparatorValue * const
      sepSeparatorValuesPtr = sep.separatorValues->ptr; 
    
    /*
     * There are 3 cases.
     * 1) AI exists and REM exist
     * 2) AI exists and REM doesnt exist
     * 3) AI does not exist, but REM exists
     * AI not exist and REM not exist can't occur.
     */

    unsigned accIndex;
    // TODO: optimize this check away out of loop.
    if (sepOrigin.hAccumulatedIntersection.size() > 0) {
      // an accumulated intersection exists.

      sepOrigin.accPacker.pack((unsigned**)sepSharedStructure.accDiscreteValuePtrs.ptr,
			       &CliqueBuffer::packedVal[0]);
      unsigned* accIndexp =
	sep.iAccHashMap->find(&CliqueBuffer::packedVal[0]);

      // we should always find something or else something is wrong.
      assert ( accIndexp != NULL ); 
      accIndex = *accIndexp;

    } else {
      // no accumulated intersection exists, everything
      // is in the remainder.
      accIndex = 0;
    }

    if (sepSharedStructure.remDiscreteValuePtrs.size() == 0) {
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
      ConditionalSeparatorTable::AISeparatorValue& sv
	= sepSeparatorValuesPtr[accIndex];
	
      sepOrigin.remPacker.pack((unsigned**)sepSharedStructure.remDiscreteValuePtrs.ptr,
			       &CliqueBuffer::packedVal[0]);

      unsigned* remIndexp =
	sv.iRemHashMap.find(&CliqueBuffer::packedVal[0]);

      if (remIndexp == NULL ) {
	// print out the rvs.
	fprintf(stderr,"ERROR: can't find separator rvs values from parent clique in fwrd hash table. Separator rv values follow.\n");
	printRVSetAndValues(stderr,sepSharedStructure.fNodes,true);
	fprintf(stderr,"Clique random variables follow:\n");
	printRVSetAndValues(stderr,sharedStructure.fNodes,true);
	assert ( remIndexp != NULL );
      }
	
      // We've finally got the sep entry.  Store the sep entry's id.
      sep.forwPointer.viterbiAccIndex = accIndex;
      sep.forwPointer.viterbiRemIndex = *remIndexp;
    }
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
  :  veSeparator(false),
     separatorValueSpaceManager(1,     // starting size
				SEPARATOR_VALUE_SPACE_MANAGER_GROWTH_RATE,   // growth rate
				1,     // growth addition
				SEPARATOR_VALUE_SPACE_MANAGER_DECAY_RATE),   // decay rate 
     remainderValueSpaceManager(1,     // starting size
				REMAINDER_VALUE_SPACE_MANAGER_GROWTH_RATE,   // growth rate
				1,     // growth addition
				REMAINDER_VALUE_SPACE_MANAGER_DECAY_RATE)    // decay rate
     
{
  nodes.clear();

  // not a ve sep clique.
  veSepClique = NULL;

  skipMe = false;

  // create a set of nodes that is the intersection of the two
  set_intersection(c1.nodes.begin(),c1.nodes.end(),
		   c2.nodes.begin(),c2.nodes.end(),
		   inserter(nodes,nodes.end()));
  // assert (nodes.size() > 0);
  
}

SeparatorClique::~SeparatorClique()
{
  if (veSeparator) {
    delete veSepClique;
  }
}


#if 0
// This constructor is not used for now.
SeparatorClique::SeparatorClique(SeparatorClique& from_sep,
				 vector <RV*>& newRvs,
				 map < RVInfo::rvParent, unsigned >& ppf,
				 const unsigned int frameDelta)
  :  separatorValueSpaceManager(1,     // starting size
				2.0,   // growth rate
				1,     // growth addition
				0.90), // decay rate 
     remainderValueSpaceManager(1,     // starting size
				2.0,   // growth rate
				1,     // growth addition
				0.90)  // decay rate 
{
  set<RV*>::iterator it;

  skipMe = false;

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
#endif








/*-
 *-----------------------------------------------------------------------
 * SeparatorClique::clearSeparatorValueCache()
 *   
 *   Clear out hash table and other memories for use between segments.
 *
 * Preconditions:
 *   Data structures must have been set up.
 *
 * Postconditions:
 *   All shared hash tables and value holders for this clique are cleared if 'force'
 *   is true.
 *
 * Side Effects:
 *     internal data structures will change.
 *
 * Results:
 *     none
 *
 *
 *-----------------------------------------------------------------------
 */
void
SeparatorClique::clearSeparatorValueCache(bool force) 
{
  if (force && accPacker.packedLen() > ISC_NWWOH_AI) {
    accValueHolder.prepare();
    accSepValHashSet.clear(AI_SEP_VALUE_HOLDER_STARTING_SIZE);
  }
  if (force && remPacker.packedLen() > ISC_NWWOH_RM) { 
    remValueHolder.prepare();
    remSepValHashSet.clear(REM_SEP_VALUE_HOLDER_STARTING_SIZE);
  }
  // shrink space asked for by clique values. 
  separatorValueSpaceManager.decay();
  remainderValueSpaceManager.decay(); 
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
    // VE separators have a different notion of hidden (see below).

    if ((!veSeparator && rv->hidden())
	||
	(veSeparator && !rv->discreteObservedImmediate()))
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

    if ((!veSeparator && rv->hidden())
	||
	(veSeparator && !rv->discreteObservedImmediate()))
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

  ////////////////////////////////////////////////////////////////////////
  // VE Separator support
  ////////////////////////////////////////////////////////////////////////

  // create any VE separators
  if (veSeparator) {

    ConditionalSeparatorTable::SharedLocalStructure veSepCliqueSharedStructure;

    // Need also to create the shared ConditionalSeparatorTable that
    // will be used for all instances of this SeparatorClique, this is
    // done below.

    // sanity check.
    assert ( veSepInfo.child != NULL && veSepInfo.parents.size() == veSepInfo.child->allParents.size() );

    // sort parents by increasing cardinality, for better branch
    // prediction.
    // printf("before sort");printRVSetAndCards(stdout,veSepInfo.parents);
    sort(veSepInfo.parents.begin(),veSepInfo.parents.end(),ParentCardinalityCompare());
    // printf("after sort");printRVSetAndCards(stdout,veSepInfo.parents);

    if (veSepInfo.grandChild == NULL) {
      ///////////////////////////
      // the PC case.

      ObsDiscRV* odc = (ObsDiscRV*)veSepInfo.child;

      // set any observed variables to their observed values
      // the child is guaranteed to be immediate const observed
      odc->setToObservedValue();
      // set any observed parents to their observed values.
      for (unsigned i=0;i<odc->allParents.size();i++) {
	// only discrete immediate observed RV parents are not stored,
	// the rest might vary per frame so we need to compute
	// the table for them with all possible values [0,card-1].
	if (odc->allParents[i]->discreteObservedImmediate()) {
	  ObsDiscRV*odrv = (ObsDiscRV*)odc->allParents[i];
	  odrv->setToObservedValue();
	}
      }

      // create a vector of just the hidden parents
      vector<RV*> hiddenParents;
      for (unsigned i = 0; i < odc->allParents.size() ; i++ ) {
	RV* rv = odc->allParents[i];
	// code will need to change when using continuous hidden
	// variables. Only consider a variable hidden if it is not a
	// discrete immediate observed, since any such varible might
	// still have varying values over time.
	if (!rv->discreteObservedImmediate())
	  hiddenParents.push_back(rv);
      }

      assert ( hiddenParents.size() == hAccumulatedIntersection.size() + hRemainder.size() );

      // Generate the table of parent values satisfying the observed
      // child, packing them into an array.
      PackCliqueValue parentPacker(hiddenParents);
      
      // create an array of direct pointers to the discrete hidden RV
      // values within the hidden nodes just computed above.
      sArray < DiscRVType*> hiddenNodeValPtrs;
      hiddenNodeValPtrs.resize(hiddenParents.size());
      for (unsigned i=0;i<hiddenParents.size();i++) {
	DiscRV* drv = RV2DRV(hiddenParents[i]);
	hiddenNodeValPtrs[i] = &(drv->val);
      }

      // create an array where we place the packed parent values that
      // satisfy the child.
      sArray < unsigned > packedParentVals(parentPacker.packedLen());

      unsigned num = 0;

      if (generatingVESeparatorTables == true) {

	if (message(Low)) {
	  float logProdCard = 
	    log10((double)RV2DRV(odc->allParents[0])->cardinality);
	  for (unsigned i=1;i<odc->allParents.size();i++) {
	    if (!odc->allParents[i]->discreteObservedImmediate())
	      logProdCard += log10((double)RV2DRV(odc->allParents[i])->cardinality);
	  }
	  printf("VE separator PC computation: Iterating 10^%f possible parent vals.\n",logProdCard);
	}
	odc->computeParentsSatisfyingChild(0,
					   veSepInfo.parents,
					   hiddenParents,
					   parentPacker,
					   hiddenNodeValPtrs,
					   odc,
					   packedParentVals,
					   num);
	if (num == 0) {
	  error("ERROR: found VE separator with 0 entries. No way for any parents to satisfy %s(%d)=const with non-zero probability. Makes all observations have probability 0.\n",
		odc->name().c_str(),odc->frame());
	}

	if (veSeparatorFile != NULL) {
	  // then we need to save this information to a file for next time.
	  unsigned tmp;
	  bool writeError = false;

	  // write a bunch of ID information.
	  // write a zero for case PC
	  tmp = 0;
	  if (!fwrite(&tmp,sizeof(tmp),1,veSeparatorFile))
	    writeError = true;

	  // write the cardinality of the child.
	  tmp = odc->cardinality;
	  if (!fwrite(&tmp,sizeof(tmp),1,veSeparatorFile))
	    writeError = true;

	  // write the number and cardinalities of the parents
	  tmp = odc->allParents.size();
	  if (!fwrite(&tmp,sizeof(tmp),1,veSeparatorFile))
	    writeError = true;

	  for (unsigned i=0; i < odc->allParents.size(); i++) {
	    tmp = RV2DRV(odc->allParents[i])->cardinality;
	    if (!fwrite(&tmp,sizeof(tmp),1,veSeparatorFile))
	      writeError = true;
	  }

	  // write the number of elements 
	  if (!fwrite(&num,sizeof(num),1,veSeparatorFile))
	    writeError = true;

	  // write the size of the elements
	  tmp = parentPacker.packedLen();
	  if (!fwrite(&tmp,sizeof(tmp),1,veSeparatorFile))
	    writeError = true;

	  
	  // finally write the set of parent values
	  if (!fwrite(packedParentVals.ptr,parentPacker.packedLen()*sizeof(unsigned),num,veSeparatorFile))
	    writeError = true;
	  
	  if (writeError)
	    error("ERROR: writing to PC VE separator file (%s)\n",veSeparatorFileName);


	}
      } else {
	// we must have a file to read from here.
	assert (veSeparatorFile != NULL);
	unsigned tmp;
	unsigned corrupt = 0;

	// read in and check the ID information for this separator information.
	if (fread(&tmp,sizeof(tmp),1,veSeparatorFile) != 1)
	  corrupt = 1;
	if (!corrupt && tmp != 0)
	  corrupt = 2;
	if (!corrupt) {
	  if (fread(&tmp,sizeof(tmp),1,veSeparatorFile) != 1)
	    corrupt = 3;
	}
	if (!corrupt && tmp != odc->cardinality)
	  corrupt = 4;
	if (!corrupt) {
	  if (fread(&tmp,sizeof(tmp),1,veSeparatorFile) != 1)
	    corrupt = 5;
	}
	if (!corrupt && tmp != odc->allParents.size())
	  corrupt = 6;
	if (!corrupt) {
	  for (unsigned i=0; i < odc->allParents.size(); i++) {
	    if (fread(&tmp,sizeof(tmp),1,veSeparatorFile) != 1) {
	      corrupt = 7;
	      break;
	    }
	    if (!corrupt && tmp != RV2DRV(odc->allParents[i])->cardinality) {
	      corrupt = 8;
	      break;
	    }
	  }
	}
	if (!corrupt) {
	  if (fread(&num,sizeof(num),1,veSeparatorFile) != 1)
	    corrupt = 9; 
	  // we can't check num.
	}
	if (!corrupt) {
	  if (fread(&tmp,sizeof(tmp),1,veSeparatorFile) != 1)
	    corrupt = 10;
	}
	if (!corrupt && tmp != parentPacker.packedLen())
	  corrupt = 11;
	packedParentVals.resize(num*parentPacker.packedLen());
	if (!corrupt) {
	  if (fread(packedParentVals.ptr,parentPacker.packedLen()*sizeof(unsigned),num,veSeparatorFile) != num)
	    corrupt = 12;
	}

	if (corrupt)
	  error("ERROR: corrupt/wrong PC VE separator file (%s) with respect to current structure/triangulation/command options. Reason %d.\n",veSeparatorFileName,corrupt);
      }

      infoMsg(Low,"VE separator PC generation: %d parent vals satisfying this case.\n",num);

      // Now create an ConditionalSeparatorTable and insert all
      // of the packed parent values we generated or got above.

      veSepClique = new ConditionalSeparatorTable(*this,
						  veSepCliqueSharedStructure);

#ifdef CHECK_VE_SEP_UNIQUE
      set < vector < unsigned > > tmpset;
#endif      

      for (unsigned i=0;i<num;i++) {
	parentPacker.unpack(&(packedParentVals.ptr[i*parentPacker.packedLen()]),
			    hiddenNodeValPtrs.ptr);
#ifdef CHECK_VE_SEP_UNIQUE
	{
	  vector < unsigned > foo(hiddenNodeValPtrs.size());
	  for (unsigned i=0;i<foo.size();i++) {
	    foo[i] = *(hiddenNodeValPtrs[i]);
	  }
	  if (tmpset.find(foo) != tmpset.end())
	    printf("Entry alreay in set\n");
	  else
	    tmpset.insert(foo);
	}
#endif
	if (message(Max+5)) {
	  printf("Inserting into PC VE separator:"); printRVSetAndValues(stdout,hiddenParents);
	}
	// insert current RV values into the separator. 
	veSepClique->insert(*this,
			    veSepCliqueSharedStructure);
      }

    } else {
      //////////////////////////////////////////////
      // the PCG case, grand-child is not NULL.

      // We know grandChild is observed.

      // set any observed parents to their observed values
      // Child and/or an parents of child might be observed. Need to check
      // those cases.
      // observed discrete grand child
      ObsDiscRV* odgc = (ObsDiscRV*)veSepInfo.grandChild;
      // dicrete child
      DiscRV* dc = (DiscRV*)veSepInfo.child;

      // the grandchild is guaranteed to be immediate const observed
      odgc->setToObservedValue();

      // child is going to be hidden, otherwise not PCG case.
      // if (dc->discreteObservedImmediate()) {
      //   ObsDiscRV* odc = (ObsDiscRV*)dc;
      //   odc->setToObservedValue();
      // }

      // set any observed parents to their observed values.
      for (unsigned i=0;i<dc->allParents.size();i++) {
	// only discrete immediate observed RV parents are not stored,
	// the rest might vary per frame so we need to compute
	// the table for them with all possible values [0,card-1].
	if (dc->allParents[i]->discreteObservedImmediate()) {
	  ObsDiscRV*odrv = (ObsDiscRV*)dc->allParents[i];
	  odrv->setToObservedValue();
	}
      }

      // create a vector of just the hidden (parents and child)
      vector<RV*> hiddenParents;
      for (unsigned i = 0; i < dc->allParents.size() ; i++ ) {
	RV* rv = dc->allParents[i];
	// code will need to change when using continuous hidden
	// variables.  only consider a variable hidden if it is not a
	// discrete immediate observed, since any such varible might
	// still have varying values over time.
	if (!rv->discreteObservedImmediate())
	  hiddenParents.push_back(rv);
      }
      // guaranteed not to be observed
      assert ( dc->hidden() );
      hiddenParents.push_back(dc);

      assert ( hiddenParents.size() == hAccumulatedIntersection.size() + hRemainder.size() );

      // Generate the table of parent values and child satisfying the
      // observed grandchild, packing them into an array.
      PackCliqueValue parentPacker(hiddenParents);
      
      // create an array of direct pointers to the discrete hidden
      // RV values within the hidden nodes just computed above.n 
      sArray < DiscRVType*> hiddenNodeValPtrs;
      hiddenNodeValPtrs.resize(hiddenParents.size());
      for (unsigned i=0;i<hiddenParents.size();i++) {
	DiscRV* drv = RV2DRV(hiddenParents[i]);
	hiddenNodeValPtrs[i] = &(drv->val);
      }

      // create an array where we place the packed parent values that
      // satisfy the child.
      sArray < unsigned > packedParentVals(parentPacker.packedLen());

      unsigned num = 0;

      if (generatingVESeparatorTables == true) {
	if (message(Low)) {
	  float logProdCard = 
	    log10((double)RV2DRV(dc->allParents[0])->cardinality);
	  for (unsigned i=1;i<dc->allParents.size();i++) {
	    if (!dc->allParents[i]->discreteObservedImmediate())
	      logProdCard += log10((double)RV2DRV(dc->allParents[i])->cardinality);
	  }
	  printf("VE separator PCG computation: Iterating 10^%f possible parent vals.\n",logProdCard);
	}

	dc->computeParentsChildSatisfyingGrandChild(0,
						    veSepInfo.parents,
						    hiddenParents,
						    parentPacker,
						    hiddenNodeValPtrs,
						    dc,
						    odgc,
						    packedParentVals,
						    num);

	if (num == 0) {
	  error("ERROR: found VE separator with 0 entries. No way for any parents to satisfy %s(%d)=const with non-zero probability. Makes all observations have probability 0.\n",
		odgc->name().c_str(),odgc->frame());
	}

	if (veSeparatorFile != NULL) {
	  // then we need to save this information to a file for next time.
	  unsigned tmp;
	  bool writeError = false;

	  // write a bunch of ID information.
	  // write a zero for case PC
	  tmp = 1;
	  if (!fwrite(&tmp,sizeof(tmp),1,veSeparatorFile))
	    writeError = true;

	  // write the cardinality of the grandchild.
	  tmp = odgc->cardinality;
	  if (!fwrite(&tmp,sizeof(tmp),1,veSeparatorFile))
	    writeError = true;

	  // write the cardinality of the child.
	  tmp = dc->cardinality;
	  if (!fwrite(&tmp,sizeof(tmp),1,veSeparatorFile))
	    writeError = true;

	  // write the number and cardinalities of the parents
	  tmp = dc->allParents.size();
	  if (!fwrite(&tmp,sizeof(tmp),1,veSeparatorFile))
	    writeError = true;

	  for (unsigned i=0; i < dc->allParents.size(); i++) {
	    tmp = RV2DRV(dc->allParents[i])->cardinality;
	    if (!fwrite(&tmp,sizeof(tmp),1,veSeparatorFile))
	      writeError = true;
	  }

	  // write the number of elements 
	  if (!fwrite(&num,sizeof(num),1,veSeparatorFile))
	    writeError = true;

	  // write the size of the elements
	  tmp = parentPacker.packedLen();
	  if (!fwrite(&tmp,sizeof(tmp),1,veSeparatorFile))
	    writeError = true;

	  
	  // finally write the set of parent values
	  if (!fwrite(packedParentVals.ptr,parentPacker.packedLen()*sizeof(unsigned),num,veSeparatorFile))
	    writeError = true;
	  
	  if (writeError)
	    error("ERROR: writing to PCG VE separator file (%s)\n",veSeparatorFileName);


	}
      } else {
	// we must have a file to read from here.
	assert (veSeparatorFile != NULL);
	unsigned tmp;
	unsigned corrupt = 0;

	// read in and check the ID information for this separator information.
	if (fread(&tmp,sizeof(tmp),1,veSeparatorFile) != 1)
	  corrupt = 1;
	if (!corrupt && tmp != 1)
	  corrupt = 2;
	if (!corrupt) {
	  if (fread(&tmp,sizeof(tmp),1,veSeparatorFile) != 1)
	    corrupt = 3;
	}
	if (!corrupt && tmp != odgc->cardinality)
	  corrupt = 4;
	if (!corrupt) {
	  if (fread(&tmp,sizeof(tmp),1,veSeparatorFile) != 1)
	    corrupt = 5;
	}
	if (!corrupt && tmp != dc->cardinality)
	  corrupt = 6;
	if (!corrupt) {
	  if (fread(&tmp,sizeof(tmp),1,veSeparatorFile) != 1)
	    corrupt = 7;
	}
	if (!corrupt && tmp != dc->allParents.size())
	  corrupt = 8;
	if (!corrupt) {
	  for (unsigned i=0; i < dc->allParents.size(); i++) {
	    if (fread(&tmp,sizeof(tmp),1,veSeparatorFile) != 1) {
	      corrupt = 9;
	      break;
	    }
	    if (!corrupt && tmp != RV2DRV(dc->allParents[i])->cardinality) {
	      corrupt = 10;
	      break;
	    }
	  }
	}
	if (!corrupt) {
	  if (fread(&num,sizeof(num),1,veSeparatorFile) != 1)
	    corrupt = 11; 
	  // we can't check num.
	}
	if (!corrupt) {
	  if (fread(&tmp,sizeof(tmp),1,veSeparatorFile) != 1)
	    corrupt = 12;
	}
	if (!corrupt && tmp != parentPacker.packedLen())
	  corrupt = 13;
	packedParentVals.resize(num*parentPacker.packedLen());
	if (!corrupt) {
	  if (fread(packedParentVals.ptr,parentPacker.packedLen()*sizeof(unsigned),num,veSeparatorFile) != num)
	    corrupt = 14;
	}

	if (corrupt)
	  error("ERROR: corrupt/wrong PCG VE separator file (%s) with respect to current structure/triangulation/command options. Reason %d.\n",veSeparatorFileName,corrupt);

      }

      infoMsg(Low,"VE separator PCG generation: %d (parent,child) combinations satisfying this case.\n",num);

      // now create an ConditionalSeparatorTable and insert all
      // of the packed parent values we generated or got above.

      veSepClique = new ConditionalSeparatorTable(*this,veSepCliqueSharedStructure);
      for (unsigned i=0;i<num;i++) {
	parentPacker.unpack(&(packedParentVals.ptr[i*parentPacker.packedLen()]),
			    hiddenNodeValPtrs.ptr);
	if (message(Max+5)) {
	  printf("Inserting into PCG VE separator:"); printRVSetAndValues(stdout,hiddenParents);
	}
	// insert current RV values into the separator. 
	// TODO: make this return bool to indicate if found or not.
	veSepClique->insert(*this,
			    veSepCliqueSharedStructure);
      }
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
  fprintf(f,"%sSeparator information: %d acc packed bits (%d words, %d splits), %d rem packed bits (%d words, %d splits)\n",
	  (veSeparator?"VE ":""),
	  accPacker.packedLenBits(),accPacker.packedLen(),accPacker.numSplits(),
	  remPacker.packedLenBits(),remPacker.packedLen(),remPacker.numSplits());

  fprintf(f,"%ld Nodes: ",(unsigned long)nodes.size()); printRVSetAndCards(f,nodes);
  fprintf(f,"%ld Acc Inter: ",(unsigned long)accumulatedIntersection.size()); printRVSet(f,accumulatedIntersection);  
  fprintf(f,"%ld Hid Acc Inter: ",(unsigned long)hAccumulatedIntersection.size()); printRVSet(f,hAccumulatedIntersection);  
  fprintf(f,"%ld remainder: ",(unsigned long)remainder.size()); printRVSet(f,remainder);  
  fprintf(f,"%ld hRemainder: ",(unsigned long)hRemainder.size()); printRVSet(f,hRemainder);

  if (veSeparator) {
    fprintf(f,"%ld VE sep parents: ",(unsigned long)veSepInfo.parents.size()); printRVSet(f,veSepInfo.parents);
    fprintf(f,"VE sep child: "); veSepInfo.child->printNameFrame(f);
    if (veSepInfo.grandChild != NULL) {
      fprintf(f,"VE sep grand child: "); veSepInfo.grandChild->printNameFrame(f);
    }
  }

}



/*-
 *-----------------------------------------------------------------------
 * SeparatorClique::reportMemoryUsageTo()
 *
 *    Report current memory usage of this inference sep clique *and*
 *    the origin sep clique to the file in units of MBs.
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
SeparatorClique::
reportMemoryUsageTo(FILE *f)
{
  // Memory: OSC=Origin Separator Clique, AI = accumulated intersection, RM = remainder
  if (accPacker.packedLen() > ISC_NWWOH_AI) {
    fprintf(f,"OSC(AIVH=%luMB,AIHS=%luMB,",
	    1+accValueHolder.bytesRequested()/(1024*1024),
	    1+accSepValHashSet.bytesRequested()/(1024*1024));
  } else {
    fprintf(f,"OSC(AIVH=0MB,AIHS=0MB,");
  }
  if (remPacker.packedLen() > ISC_NWWOH_RM) { 
    fprintf(f,"RMVH=%luMB,RMHS=%luMB)",
	    1+remValueHolder.bytesRequested()/(1024*1024),
	    1+remSepValHashSet.bytesRequested()/(1024*1024));
  } else {
    fprintf(f,"RMVH=0MB,RMHS=0MB)");
  }
}



////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        ConditionalSeparatorTable support
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////



/*-
 *-----------------------------------------------------------------------
 * ConditionalSeparatorTable::ConditionalSeparatorTable()
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
ConditionalSeparatorTable::
ConditionalSeparatorTable(SeparatorClique& origin)
  : separatorValues(NULL),iAccHashMap(NULL)
{
  init(origin);
}

void ConditionalSeparatorTable::init(SeparatorClique& origin) 
{

  if (origin.veSeparator) {
    // For VE separators, our origin contains the separator tables
    // already pre-generated and constant accross all instances of
    // this separator. Therefore, we special case here and just copy
    // over the pointer values tha are stored in the VE separator that
    // is contained in the origin.
    //
    // We need to make sure not to delete these values however when
    // this object is deleted since normally the separator tables will
    // be reused many times.
    separatorValues = origin.veSepClique->separatorValues;
    iAccHashMap =  origin.veSepClique->iAccHashMap;
    numSeparatorValuesUsed = origin.veSepClique->numSeparatorValuesUsed;
    setToVeSeparatorId();
  } else {

    clearInferenceMemory();

    separatorValues = new cArray< AISeparatorValue >;
    // allocate at one value for now.
    if (origin.hAccumulatedIntersection.size() == 0) {
      // in this case, we'll only need one and never more.
      separatorValues->resize(1);
      // there will always be one used value here.
      numSeparatorValuesUsed = 1;
      if (origin.hRemainder.size() > 0) {
	// So we have no accumulated intersection, and a remainder which
	// means we are in a good position to predict the size of the
	// (necessarily single) remainder vectors from the previous
	// times we used it. Therefore, we do just that, but only in
	// this case.
	separatorValues->ptr[0].remValues.resize(origin.remainderValueSpaceManager.currentSize());
	new (&separatorValues->ptr[0].iRemHashMap)VHashMapUnsignedUnsignedKeyUpdatable
	  (origin.remPacker.packedLen(),origin.remainderValueSpaceManager.currentSize());
	// new (&separatorValues->ptr[0].iRemHashMap)VHashMapUnsignedUnsignedKeyUpdatable(origin.remPacker.packedLen(),REM_HASH_MAP_STARTING_SIZE);
      } else {
	// The separator consists of all observed nodes. 
	// Search in file for key string
	// "ALLOCATE_REMVALUES_ALL_OBSERVED'
	// to find where the nec. single entry is allocated.
	// @@@ try allocating it here
	separatorValues->ptr[0].remValues.resize(1);
	separatorValues->ptr[0].numRemValuesUsed = 0;	
      }
    } else {
      // start with something a bit larger
      // TODO: optimize this.
      const unsigned starting_size = origin.separatorValueSpaceManager.currentSize(); // 3,2000;
      separatorValues->resize(starting_size);
      if (origin.hRemainder.size() > 0) {
	for (unsigned i=0;i<starting_size;i++) {
	  // need to re-construct individual hash tables.
	  new (&separatorValues->ptr[i].iRemHashMap)VHashMapUnsignedUnsignedKeyUpdatable
	    (origin.remPacker.packedLen(),REM_HASH_MAP_STARTING_SIZE);
	  // TODO: while we potentially could preallocate default size
	  // of separatorValues->ptr[i].remValues.resize(default); here,
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
      iAccHashMap = new VHashMapUnsignedUnsignedKeyUpdatable
	(origin.accPacker.packedLen(),starting_size); // 2
      numSeparatorValuesUsed = 0;
    }
  }
}

ConditionalSeparatorTable::SharedLocalStructure::
SharedLocalStructure(SeparatorClique& _origin,
		     vector <RV*>& newRvs,
		     map < RVInfo::rvParent, unsigned >& ppf,
		     const unsigned int frameDelta)
{

  origin = &_origin;

  set<RV*>::iterator it;

  // clone over nodes RVs.
  fNodes.resize(origin->nodes.size());
  unsigned i=0;
  for (it = origin->nodes.begin();
       it != origin->nodes.end();
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
  fAccumulatedIntersection.resize(origin->accumulatedIntersection.size());
  for (it = origin->accumulatedIntersection.begin();
       it != origin->accumulatedIntersection.end();
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
    // BUG!! (i needs to be incremented, check other cases too, but then this variable is never used).
    fAccumulatedIntersection[i] = nrv;
  }

  i=0;
  fRemainder.resize(origin->remainder.size());
  for (it = origin->remainder.begin();
       it != origin->remainder.end();
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
  accDiscreteValuePtrs.resize(origin->hAccumulatedIntersection.size());
  for (i=0;i<accDiscreteValuePtrs.size();i++) {
    // get the hidden rv for this location
    RV* rv = origin->hAccumulatedIntersection[i];;
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
  remDiscreteValuePtrs.resize(origin->hRemainder.size());
  for (i=0;i<remDiscreteValuePtrs.size();i++) {
    // get the hidden rv for this location
    RV* rv = origin->hRemainder[i];
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
}



// Next, version of the above constructor to create specifically for a
// VE separator that only lives and stays in its origin separator, and
// used in SeparatorClique::prepareForUnrolling(), where the set of
// RVs are the same as the origin separator, and there is no frame
// delta. This separator is then used to create copies of the
// pointers in the table data structures when we need
// a ve separator that is iterated against a set of random variables.
// 
// A inference vec separator is used just for its values at every
// possible time frame.  See comments above.
//
// this routine assumes that the sepSharedStructure argument is empty
// and that member variables can be removed/resized.

ConditionalSeparatorTable
::ConditionalSeparatorTable(SeparatorClique& origin,
			    ConditionalSeparatorTable::SharedLocalStructure& sepSharedStructure)
  : separatorValues(NULL),iAccHashMap(NULL)
{

  set<RV*>::iterator it;

  sepSharedStructure.origin = &origin;

  // this table is a ve separator
  assert ( origin.veSeparator );
  setToVeSeparatorId();

  // clone over nodes RVs.
  sepSharedStructure.fNodes.resize(origin.nodes.size());
  unsigned i=0;
  for (it = origin.nodes.begin();
       it != origin.nodes.end();
       it++) {
    RV* rv = (*it);
    sepSharedStructure.fNodes[i++] = rv;
  }

  i=0;
  sepSharedStructure.fAccumulatedIntersection.resize(origin.accumulatedIntersection.size());
  for (it = origin.accumulatedIntersection.begin();
       it != origin.accumulatedIntersection.end();
       it++) {
    RV* rv = (*it);
    sepSharedStructure.fAccumulatedIntersection[i] = rv;
  }

  i=0;
  sepSharedStructure.fRemainder.resize(origin.remainder.size());
  for (it = origin.remainder.begin();
       it != origin.remainder.end();
       it++) {
    RV* rv = (*it);
    sepSharedStructure.fRemainder[i++] = rv;
  }

  // Separator accumulated intersection values only store/hash values
  // of hidden (thus necessarily discrete) variables since they are
  // the only thing that change.
  sepSharedStructure.accDiscreteValuePtrs.resize(origin.hAccumulatedIntersection.size());
  for (i=0;i<sepSharedStructure.accDiscreteValuePtrs.size();i++) {
    // get the hidden rv for this location
    RV* rv = origin.hAccumulatedIntersection[i];;
    DiscRV* drv = (DiscRV*)rv;
    // grab a pointer directly to its value for easy access later.
    sepSharedStructure.accDiscreteValuePtrs[i] = &(drv->val);
  }

  // Separator remainder values only store/hash values of hidden (thus
  // necessarily discrete) variables since they are the only thing
  // that change.
  sepSharedStructure.remDiscreteValuePtrs.resize(origin.hRemainder.size());
  for (i=0;i<sepSharedStructure.remDiscreteValuePtrs.size();i++) {
    // get the hidden rv for this location
    RV* rv = origin.hRemainder[i];
    // hidden nodes are always discrete (in this version).
    DiscRV* drv = (DiscRV*)rv;
    // grab a pointer directly to its value for easy access later.
    sepSharedStructure.remDiscreteValuePtrs[i] = &(drv->val);
  }

  separatorValues = new cArray< AISeparatorValue >;
  // allocate at one value for now.
  if (origin.hAccumulatedIntersection.size() == 0) {
    // in this case, we'll only need one and never more.
    separatorValues->resize(1);
    // there will always be one used value here.
    numSeparatorValuesUsed = 1;
    if (origin.hRemainder.size() > 0) {
      // So we have no accumulated intersection, and a remainder which
      // means we are in a good position to predict the size of the
      // (necessarily single) remainder vectors from the previous
      // times we used it. Therefore, we do just that, but only in
      // this case.
      separatorValues->ptr[0].remValues.resize(origin.remainderValueSpaceManager.currentSize());
      new (&separatorValues->ptr[0].iRemHashMap)VHashMapUnsignedUnsignedKeyUpdatable
	(origin.remPacker.packedLen(),origin.remainderValueSpaceManager.currentSize());
      // new (&separatorValues->ptr[0].iRemHashMap)VHashMapUnsignedUnsignedKeyUpdatable(origin.remPacker.packedLen(),REM_HASH_MAP_STARTING_SIZE);
    } else {
      // The separator consists of all observed nodes. 
      // Search in file for key string
      // "ALLOCATE_REMVALUES_ALL_OBSERVED'
      // to find where the nec. single entry is allocated.
      // @@@ try allocating it here
      separatorValues->ptr[0].remValues.resize(1);
      separatorValues->ptr[0].numRemValuesUsed = 0;

    }
  } else {
    // start with something a bit larger
    // TODO: optimize this.
    const unsigned starting_size = origin.separatorValueSpaceManager.currentSize(); // 3,2000;
    separatorValues->resize(starting_size);
    if (origin.hRemainder.size() > 0) {
      for (unsigned i=0;i<starting_size;i++) {
	// need to re-construct individual hash tables.
	new (&separatorValues->ptr[i].iRemHashMap)VHashMapUnsignedUnsignedKeyUpdatable
	  (origin.remPacker.packedLen(),REM_HASH_MAP_STARTING_SIZE);
	// TODO: while we potentially could preallocate default size
	// of separatorValues->ptr[i].remValues.resize(default); here,
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
    iAccHashMap = new VHashMapUnsignedUnsignedKeyUpdatable
      (origin.accPacker.packedLen(),starting_size); // 2
    numSeparatorValuesUsed = 0;
  }

}


/*-
 *-----------------------------------------------------------------------
 * ConditionalSeparatorTable::insert()
 *
 *    Insert whatever the current RV values are set to into the current
 *    inference separator (based in their RV values).
 *
 * Preconditions:
 *   1) separator tables must be created, meaning that the same
 *      conditions must holed as just before ceSendToOutgoingSeparator() is
 *      called.
 *
 * Postconditions:
 *    Separator table has entry added, according to current accumulated intersection
 *    and remainder, and current RV values.
 *
 * Side Effects:
 *    potential memory allocations and hash table adjustments.
 *
 * Results:
 *     nothing
 *
 *-----------------------------------------------------------------------
 */
void ConditionalSeparatorTable::insert(SeparatorClique& origin,
				       ConditionalSeparatorTable::SharedLocalStructure& sepSharedStructure)
{

  // keep a local variable copy of this around to avoid potential
  // dereferencing.  This one cannot be const since it might change
  // during a resize, in which case we need to reassign this variable.
  ConditionalSeparatorTable::AISeparatorValue * 
    sepSeparatorValuesPtr = separatorValues->ptr; 

  // there must be someplace to insert in this case.
  assert (!(origin.hAccumulatedIntersection.size() == 0 && origin.hRemainder.size() == 0));

  // precompute some constants.
  const bool isc_nwwoh_ai_p = (origin.accPacker.packedLen() <= ISC_NWWOH_AI);
  const bool isc_nwwoh_rm_p = (origin.remPacker.packedLen() <= ISC_NWWOH_RM);
  const bool sep_origin_hAccumulatedIntersection_exists_p =
    (origin.hAccumulatedIntersection.size() > 0);
  const bool sep_remDiscreteValuePtrs_exists_p = 
    (sepSharedStructure.remDiscreteValuePtrs.size() > 0);
  
  unsigned accIndex;
  if (sep_origin_hAccumulatedIntersection_exists_p) { 
    // an accumulated intersection exists.

    // make sure there is at least one available accumulated intersection entry
    assert ( numSeparatorValuesUsed <= separatorValues->size());
    if (numSeparatorValuesUsed >= separatorValues->size()) {
      
      infoMsg(Max+5,"ac-rsz,");

      const unsigned old_size = separatorValues->size();
      // TODO: optimize this size re-allocation.
      if (numSeparatorValuesUsed >= origin.separatorValueSpaceManager.currentSize()) 
	origin.separatorValueSpaceManager.advanceToNextSize();
      separatorValues->resizeAndCopy(origin.separatorValueSpaceManager.currentSize()); 
      sepSeparatorValuesPtr = separatorValues->ptr;
      if (isc_nwwoh_ai_p) {
	// Then the above resize just invalided all our pointers to
	// keys (which in this case are compressed RV values for the
	// acc inter), but it did not invalidate the hash items (which
	// in this case are the array indices in the accumulated
	// intersection corresponding to a given compressed acc intr
	// RV values). We thus go through and correct the key pointers
	// within the hash table.  TODO: think of a better way to do
	// this that also looses no efficiency.
	for (unsigned i=0;i<iAccHashMap->tableSize();i++) {
	  if (!iAccHashMap->tableEmpty(i)) {
	    iAccHashMap->tableKey(i)
	      = &(sepSeparatorValuesPtr[iAccHashMap->tableItem(i)].val[0]);
	  }
	}
      }
      const unsigned new_size = separatorValues->size();
      // if (remDiscreteValuePtrs.size() > 0) {
      if (sep_remDiscreteValuePtrs_exists_p) {
	for (unsigned i=old_size;i<new_size;i++) {
	  // re-construct hash tables only for new entries.
	  new (&sepSeparatorValuesPtr[i].iRemHashMap)
	    VHashMapUnsignedUnsignedKeyUpdatable
	    (origin.remPacker.packedLen(),REM_HASH_MAP_STARTING_SIZE);
	  // TODO: potentially preallocate default size of  
	  // separatorValues->ptr[i].remValues.resize(default);
	  // TODO: potentially create zero size here, and only
	  //       grow bigger when we start adding things.
	}
      }
    }
      
    unsigned *accKey;
    // TODO: optimize this check out of loop.
    if (isc_nwwoh_ai_p) {
      accKey = &(sepSeparatorValuesPtr[numSeparatorValuesUsed].val[0]);
      origin.accPacker.pack((unsigned**)sepSharedStructure.accDiscreteValuePtrs.ptr,
			    accKey);
    } else {
      accKey = origin.accValueHolder.curCliqueValuePtr();
      origin.accPacker.pack((unsigned**)sepSharedStructure.accDiscreteValuePtrs.ptr,
			    accKey);
      // check if this value combination already lives in
      // origin's value holder hash table and if so, use that.
      bool foundp;
      accKey = origin.accSepValHashSet.insert(accKey,foundp);
      if (!foundp) {
	// only allocate a new value if it was inserted.
	origin.accValueHolder.allocateCurCliqueValue();
      }
      // store the pointer in case we use it.
      sepSeparatorValuesPtr[numSeparatorValuesUsed].ptr = accKey;
    }

      
    bool foundp;
    unsigned* accIndexp =
      iAccHashMap->insert(accKey,
			  numSeparatorValuesUsed,
			  foundp);

    if (!foundp) {
      //  add the values we just used. 
      numSeparatorValuesUsed++;
    }
    accIndex = *accIndexp;

    infoMsg(Max+5,"inst:ai=%d,ky=%X,",accIndex,*accKey);

    // TODO: optimize this check out of loop.
    // if (remDiscreteValuePtrs.size() == 0) {
    if (!sep_remDiscreteValuePtrs_exists_p) {
      // 2) AI exists and REM doesnt exist
      // Then this separator is entirely covered by one or 
      // more other separators earlier in the order.

      // go ahead and insert it here to the 1st entry (entry 0).

      // handy reference for readability.
      ConditionalSeparatorTable::AISeparatorValue& sv
	= sepSeparatorValuesPtr[accIndex];

      // Accumulate the clique's
      // probability into this separator's probability.
      if (sv.remValues.size() < 1) {
	// This must be first time for this entry.
	// Search for tag 'ALLOCATE_REMVALUES_OPTION' in this file for where else this could
	// be done.
	sv.remValues.resize(1);
	sv.numRemValuesUsed = 1;	  
	// initialize and assign.
	sv.remValues.ptr[0].p = 1.0; // probability being assigned is unity
	// if (JunctionTree::viterbiScore)
	// sv.remValues.ptr[0].backPointer = cvn;
      } else {
	// already there so must have hit before.
	// for now, die with an assertion if this case occurs (can't insert twice)
	assert ( 0 );
      }
      return;
    }

  } else {
    accIndex = 0;
    infoMsg(Max+5,"inst:ai=%d,",accIndex);
  }




  // If we're here, then we are guaranteed must have some remainder
  // pointers, i.e., we could do:
  //    assert (remDiscreteValuePtrs.size() > 0);
  // So, the remainder exists in this separator.
  // either:
  //   1) AI exists and REM exist
  //     or
  //   3) AI does not exist (accIndex == 0), but REM exists
  // 

  // keep handy reference for readability.
  ConditionalSeparatorTable::AISeparatorValue& sv
    = sepSeparatorValuesPtr[accIndex];
    
  // make sure there is at least one available entry
  assert (sv.numRemValuesUsed <= sv.remValues.size());
  if (sv.numRemValuesUsed >= sv.remValues.size()) {

    infoMsg(Max+5,"rm-rsz,u=%d,f=%d,",sv.numRemValuesUsed,sv.remValues.size());

    // TODO: optimize this growth rate.
    // start small but grow fast.
    // sv.remValues.resizeAndCopy(1+sv.remValues.size()*2); // *3
    sv.remValues.resizeAndCopy(origin.remainderValueSpaceManager.nextSizeFrom(sv.remValues.size()));
    origin.remainderValueSpaceManager.setCurrentAllocationSizeIfLarger(sv.remValues.size());

    infoMsg(Max+5,"=%d,",sv.remValues.size());

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
    origin.remPacker.pack((unsigned**)sepSharedStructure.remDiscreteValuePtrs.ptr,
			  remKey);
  } else {
    // grab pointer to next packed clique value to be used.
    remKey = origin.remValueHolder.curCliqueValuePtr();
    origin.remPacker.pack((unsigned**)sepSharedStructure.remDiscreteValuePtrs.ptr,
			  remKey);
    // check if this value combination already lives in
    // origin's value holder hash table and if so, use that.
    bool foundp;
    remKey = origin.remSepValHashSet.insert(remKey,foundp);
    if (!foundp) {
      // only allocate a new value if it was inserted.
      origin.remValueHolder.allocateCurCliqueValue();
    }
    // store the pointer in case we use it.
    sv.remValues.ptr[sv.numRemValuesUsed].ptr = remKey;
  }

  bool foundp;
  unsigned* remIndexp =
    sv.iRemHashMap.insert(remKey,
			  sv.numRemValuesUsed,
			  foundp);

  infoMsg(Max+5,"inst:ky=%X,rii=%d\n",*remKey,*remIndexp);

  if (!foundp) {
    // add the values we just used. 
    sv.numRemValuesUsed++;
  } else {
    // already found
    infoMsg(Max+5,"already found\n");
    fflush(stdout);
    fflush(stderr);
    assert ( 0 );
  }

  // We've finally got the entry, assign the probability as unity.
  sv.remValues.ptr[*remIndexp].p = 1.0;


}


set <RV*> 
ConditionalSeparatorTable::SharedLocalStructure::returnRVsAsSet()
{
  set<RV*> rc;
  for (unsigned i=0;i<fNodes.size();i++) {
    rc.insert(fNodes[i]);
  }
  return rc;
}



/*-
 *-----------------------------------------------------------------------
 * ConditionalSeparatorTable::ceSeparatorPrune()
 *
 *    Collect Evidence, Separator Prune: This routine will prune away
 *    part of a previously instantiated separator based on the current
 *    separator beam width.
 *
 * Preconditions:
 *   1) separator table must be created, meaning that either:
 *
 *        MaxCliqueTable::ceSendToOutgoingSeparator()
 *      must have been called sending a message (projection downto) this separator.
 *      Separator must not be a VE separator.
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
ConditionalSeparatorTable::ceSeparatorPrune(SeparatorClique& origin)
{

  // keep a local variable copy of this around to avoid potential dereferencing.
  AISeparatorValue * const
    separatorValuesPtr = separatorValues->ptr; 

  if (origin.separatorBeam != (-LZERO)) {
    // only do this if separator beam pruning is not turned off.

    // we shouldn't have to check this since we should never be pruning
    // a VE separator.
    // if (veSeparator())
    // return;

    // compute max and current state space.
    logpr maxCEsepValue;
    unsigned originalTotalStateSpace = 0;
    for (unsigned asv=0;asv<numSeparatorValuesUsed;asv++) {
      originalTotalStateSpace += separatorValuesPtr[asv].numRemValuesUsed;
      for (unsigned rsv=0;rsv<separatorValuesPtr[asv].numRemValuesUsed;rsv++) {
	if (separatorValuesPtr[asv].remValues.ptr[rsv].p > maxCEsepValue)
	  maxCEsepValue = separatorValuesPtr[asv].remValues.ptr[rsv].p;
      }
    }

    // Check for all observed case, or a case where there is only one
    // entry, in which case we never do anything.
    if (originalTotalStateSpace == 1)
      return;

    // check if we have a zero separator, and if we do, print message.
    // Do this before separator pruning.
    if (originalTotalStateSpace == 0) {
      infoMsg(IM::Mod,"WARNING: ZERO SEPARATOR: separator with no entries. Final probability will be zero.\n");
      // nothing to prune, so we return.
      return;
    }

    // create an ininitialized variable
    logpr beamThreshold((void*)0);
    // break into the logp to avoid unnecessary zero checking.
    beamThreshold.valref() = maxCEsepValue.valref() - origin.separatorBeam;

    // pointers to the ht keys for the two entries.
    unsigned** ht_prune_key_p=NULL;
    unsigned** ht_swap_key_p=NULL;

    // go through and shrink guys less than maximum.
    unsigned newTotalStateSpace = 0;  
    for (unsigned asv=0;asv<numSeparatorValuesUsed;asv++) {
      const unsigned origNumRemValuesUsed = separatorValuesPtr[asv].numRemValuesUsed;
      for (unsigned rsv=0;rsv<separatorValuesPtr[asv].numRemValuesUsed;) {
	if (separatorValuesPtr[asv].remValues.ptr[rsv].p < beamThreshold) {

	  if (separatorValuesPtr[asv].numRemValuesUsed > 1) {


	    // We prune away entry for rsv, by swapping it in last
	    // position. Here, however, it is not as easy as with a clique
	    // separator as we have also to deal with the hash
	    // table. Specifically, we need to swap index entries in hash
	    // table as well. Note that we can not remove the hash table
	    // entry for the one that got pruned away without re-hashing
	    // the entire hash table. The reason is that if the entry that
	    // got removed was a collision for another entry that is in
	    // the table, then removing the collision will make the other
	    // entry inaccessible. Therefore, for now, the hash table
	    // does not shrink while the table does.
	    // TODO: test if it is better to just prune here and just
	    // rehash everything.

	    // the index of the entry being swapped with the
	    // one that is being pruned.
	    const unsigned swap_index = separatorValuesPtr[asv].numRemValuesUsed-1;

	    // First, get pointers to hash table index values for the two
	    // entries corresponding to the one we are prunning
	    // and the one ware swapping it with.
	    unsigned* prune_index_p;
	    unsigned* swap_index_p;

	    // the keys for the two entries.
	    unsigned* prune_key_p;
	    unsigned* swap_key_p;


	    if (origin.remPacker.packedLen() <= ISC_NWWOH_RM) {
	      prune_key_p = &(separatorValuesPtr[asv].remValues.ptr[rsv].val[0]);
	      swap_key_p = &(separatorValuesPtr[asv].remValues.ptr[swap_index].val[0]);
	    } else {
	      prune_key_p = separatorValuesPtr[asv].remValues.ptr[rsv].ptr;
	      swap_key_p = separatorValuesPtr[asv].remValues.ptr[swap_index].ptr;
	    }
	
	    prune_index_p =  separatorValuesPtr[asv].iRemHashMap.find(prune_key_p,ht_prune_key_p);
	    // it must exist
	    assert ( prune_index_p != NULL );
	    swap_index_p =  separatorValuesPtr[asv].iRemHashMap.find(swap_key_p,ht_swap_key_p);
	    // it must exist
	    assert ( swap_index_p != NULL );

	    // swap the entries in the separator remainder (rv_val,
	    // prob) table. We can't do this any earlier than here
	    // since the hash finding above uses pointers to these
	    // entries.
	    swap(separatorValuesPtr[asv].remValues.ptr[rsv],
		 separatorValuesPtr[asv].remValues.ptr[swap_index]);

	    // and swap the hash table item index values (i.e., the ht items are integer indices
	    // into the separator remainder (rv_val, prob) table. We want the hash table entry for the
	    // item that we are not pruning to now point to the separator remainder entry that is not
	    // being pruned away (rather than its old entry, which is no longer being used). This
	    // does not change the position in the hash table of this entry, rather it only changes
	    // what the hash table entry is pointing back to in the separator remainder table.
	    swap((*prune_index_p),(*swap_index_p));

	    // and swap the hash table keys if they are pointers to
	    // the arrays which just got swapped. In other words, when
	    // the separator remainder value is small enough to fit in
	    // one machine word, then the hash table key values will
	    // point directly into the slot in the separator remainder
	    // (rv_val, prob) table rather than pointing to some
	    // globally shared value pool. In this case, since the
	    // separator remainder table (containing rv values
	    // directly rather than pointers) has changed, we need to
	    // adjust the hash table so that its pointer to key is
	    // appropriate.
	    if (origin.remPacker.packedLen() <= ISC_NWWOH_RM) {
	      // printf("foobarbaz");
	      swap((*ht_prune_key_p),(*ht_swap_key_p));
	    }

	    // TODO: we don't need to swap above, rather we just need
	    // to update the slot that is not getting pruned away.

	    // decrease values
	    separatorValuesPtr[asv].numRemValuesUsed--;

	    // TODO: the above needs to be looked at for 64bit (64-bit) case.

	  } else {
	    // then separatorValuesPtr[asv].numRemValuesUsed == 1. This
	    // will happen under two conditiosn.
	    //  1) There is no remainder, meaning (sepOrigin.hRemainder.size() == 0), 
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
	    separatorValuesPtr[asv].numRemValuesUsed = 0;

	  }

	} else {
	  rsv++;
	}
      }
      newTotalStateSpace += separatorValuesPtr[asv].numRemValuesUsed;
      if (separatorValuesPtr[asv].numRemValuesUsed < origNumRemValuesUsed) {
	if (separatorValuesPtr[asv].numRemValuesUsed == 0 ) {
	  // should/could remove accumulator entry here as well. 
	}
	// - re-allocate memory & adjust hash table.
	// - separatorValuesPtr[asv].remValues
	// - possibly re-hash hash tables if necessary.
      }
    }

    infoMsg(IM::Med,"Separator beam pruning, Max cv = %f, thres = %f. Original sep state space = %d, new sep state space = %d\n",
	    maxCEsepValue.valref(),
	    beamThreshold.valref(),
	    originalTotalStateSpace,newTotalStateSpace);

  }

#if 0

  // reallocate memory so that nothing is wasted.
  for (unsigned asv=0;asv<numSeparatorValuesUsed;asv++) {
    // shrink down if we can gain more than about 1.6% of size.
    if ((separatorValuesPtr[asv].remValues.size() - separatorValuesPtr[asv].numRemValuesUsed) > separatorValuesPtr[asv].remValues.size()/64) {
      separatorValuesPtr[asv].remValues.resizeAndCopy(separatorValuesPtr[asv].numRemValuesUsed);
    }
  }

  // TODO: get the rest of this code working, will need in some cases
  // to re-hash the above hash tables. Not sure this is worth it.

#endif 

}



/*-
 *-----------------------------------------------------------------------
 * ConditionalSeparatorTable::reportMemoryUsageTo()
 *
 *    Report current memory usage of this inference sep clique *and*
 *    the origin sep clique to the file in units of MBs.
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
ConditionalSeparatorTable::
reportMemoryUsageTo(SeparatorClique& origin, FILE *f)
{
  // Memory: IC=Inference Separator  Clique, AI=accumulated intersection
  fprintf(f,"*MEM:ISC AI(used=%lu,all=%lu=%luMB),AIH(%luMB)",
	  (unsigned long)numSeparatorValuesUsed,
	  (unsigned long)separatorValues->size(),
	  (unsigned long)((1+(sizeof(AISeparatorValue)*(unsigned long)separatorValues->size())/(1024ul*1024ul))),
	  (iAccHashMap != NULL)? (unsigned long)(1 + iAccHashMap->bytesRequested()/(1024ul*1024ul)) : 0 );

  // sum up stats from the remainders
  unsigned long remUsed = 0;
  unsigned long remAllocated = 0;
  unsigned long remHashAllocated = 0;

  for (unsigned long i =0; i< numSeparatorValuesUsed; i++) {
    remUsed += separatorValues->ptr[i].numRemValuesUsed;
    remAllocated += separatorValues->ptr[i].remValues.size();
    remHashAllocated += separatorValues->ptr[i].iRemHashMap.bytesRequested();
  }
  // note: remUsed is also equal to the state space size of the separator.

  // REM = remainder
  fprintf(f,"REM(used=%lu,all=%lu=%luMB),RH(%luMB),",
	  (unsigned long)remUsed,
	  (unsigned long)remAllocated,
	  (unsigned long)((1+(sizeof(RemainderValue)*(unsigned long)remAllocated)/(1024ul*1024ul))),
	  (unsigned long)((1+remHashAllocated/(1024ul*1024ul))));
  origin.reportMemoryUsageTo(f);
  fprintf(f,"\n");

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

  // newSize *MUST* be a multiple of 'cliqueValueSize' or else.
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


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        FactorClique support
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Constructors
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


FactorClique::FactorClique(FactorInfo& _factorInfo,
			   vector <RV*>& unrolled_rvs,
			   map < RVInfo::rvParent, unsigned > ppf,
			   const unsigned offset)
{
  nodes = getRVVec(unrolled_rvs,ppf,
		   _factorInfo.variables,
		   _factorInfo.frame + offset);

  orderedNodes = getRVOVec(unrolled_rvs,ppf,
			   _factorInfo.variables,
			   _factorInfo.frame + offset);

  factorInfo = &_factorInfo;

  assert( orderedNodes.size() == nodes.size() );

}




////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        InferenceFactorClique support
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Constructors
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////





/*-
 *-----------------------------------------------------------------------
 * InferenceFactorClique::InferenceFactorClique()
 *    Clone constructor with frame delta to create a clone but under an unrolling.
 *    I.e., this isn't really a normal constructor, this is a contructor that
 *    sets up a clone of the FactorClique given by the argument from_factor. The
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
InferenceFactorClique::InferenceFactorClique(FactorClique& from_factor,
					     vector <RV*>& newRvs,
					     map < RVInfo::rvParent, unsigned >& ppf,
					     const unsigned int frameDelta)
  : origin(from_factor)
{

  // clone over nodes RVs.
  fOrderedNodes.resize(origin.orderedNodes.size());

  for (unsigned i=0;i<origin.orderedNodes.size();i++) {
    RV* rv = origin.orderedNodes[i];
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    

    if ( ppf.find(rvp) == ppf.end() ) {
      coredump("ERROR: can't find rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	       rv->name().c_str(),rv->frame(),frameDelta,
	       rvp.first.c_str(),rvp.second);
    }

    RV* nrv = newRvs[ppf[rvp]];
    fOrderedNodes[i] = nrv;
  }


}


/////////////////////////////////////////////	
/// END OF FILE
/////////////////////////////////////////////
