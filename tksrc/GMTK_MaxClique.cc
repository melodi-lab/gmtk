/*-
 * GMTK_MaxClique.cc
 *     maxClique support
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
#include "GMTK_RandomVariable.h"
#include "GMTK_DiscreteRandomVariable.h"
#include "GMTK_ContinuousRandomVariable.h"
#include "GMTK_GM.h"
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
 * separator beam width, for separator-based beam pruning.  Default value is
 * very large (1.0/0.0 = est. of infty) meaning that we do no beam
 * pruning.
 *
 */
double
SeparatorClique::separatorBeam=(-LZERO);


// TODO: put this function somewhere more generally available.
static void
printRVSetAndValues(FILE*f,sArray<RandomVariable*>& locset) 
{
  bool first = true;
  for (unsigned i=0;i<locset.size();i++) {
    RandomVariable* rv = locset[i];
    if (!first)
      fprintf(f,", ");
    fprintf(f,"%s(%d)=",rv->name().c_str(),rv->frame());
    if (!rv->discrete) {
      fprintf(f,"C");
    } else {
      DiscreteRandomVariable* drv = (DiscreteRandomVariable*)rv;
      fprintf(f,"%d",drv->val);
    }
    first = false;
  }
  fprintf(f,"\n");
}
static void
printRVSet(FILE*f,sArray<RandomVariable*>& locset) 
{
  bool first = true;
  for (unsigned i=0;i<locset.size();i++) {
    RandomVariable* rv = locset[i];
    if (!first)
      fprintf(f,", ");
    fprintf(f,"%s(%d)",rv->name().c_str(),rv->frame());
    first = false;
  }
  fprintf(f,"\n");
}
static void
printRVSet(FILE*f,set<RandomVariable*>& locset)
{
  bool first = true;
  set<RandomVariable*>::iterator it;
  for (it = locset.begin();
       it != locset.end();it++) {
    RandomVariable* rv = (*it);
    if (!first)
      fprintf(f,",");
    fprintf(f,"%s(%d)",rv->name().c_str(),rv->frame());
    first = false;
  }
  fprintf(f,"\n");
}


static void
printRVSet(FILE*f,vector<RandomVariable*>& locvec)
{
  bool first = true;
  for (unsigned i=0;i<locvec.size();i++) {
    RandomVariable* rv = locvec[i];
    if (!first)
      fprintf(f,",");
    fprintf(f,"%s(%d)",rv->name().c_str(),rv->frame());
    first = false;
  }
  fprintf(f,"\n");
}



// TODO: put this in misc support
static void
psp(FILE*f,const int numSpaceChars)
{
  int tmp = numSpaceChars;
  while (tmp--)
    fprintf(f," ");
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
			   vector <RandomVariable*>& newRvs,
			   map < RVInfo::rvParent, unsigned >& ppf,
			   const unsigned int frameDelta)
{

  set<RandomVariable*>::iterator it;
  

  // clone over nodes RVs.
  for (it = from_clique.nodes.begin();
       it != from_clique.nodes.end();
       it++) {
    RandomVariable* rv = (*it);
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    

    if ( ppf.find(rvp) == ppf.end() ) {
      coredump("ERROR: can't find rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
    }

    RandomVariable* nrv = newRvs[ppf[rvp]];
    nodes.insert(nrv);
  }

  // and clone over assigned nodes and sorted assigned nodes
  sortedAssignedNodes.reserve(from_clique.sortedAssignedNodes.size());
  for (unsigned i=0;i<from_clique.sortedAssignedNodes.size();i++) {
    RandomVariable* rv = from_clique.sortedAssignedNodes[i];
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    

    if ( ppf.find(rvp) == ppf.end() ) {
      coredump("ERROR: can't find assigned rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
    }

    RandomVariable* nrv = newRvs[ppf[rvp]];
    assignedNodes.insert(nrv);
    sortedAssignedNodes.push_back(nrv);
  }

  // do unassignedIteratedNodes
  for (it = from_clique.unassignedIteratedNodes.begin();
       it != from_clique.unassignedIteratedNodes.end();
       it++) {
    RandomVariable* rv = (*it);
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    

    if ( ppf.find(rvp) == ppf.end() ) {
      coredump("ERROR: can't find assigned rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
    }
    RandomVariable* nrv = newRvs[ppf[rvp]];
    unassignedIteratedNodes.insert(nrv);
  }

  // do cumulativeAssignedNodes
  for (it = from_clique.cumulativeAssignedNodes.begin();
       it != from_clique.cumulativeAssignedNodes.end();
       it++) {
    RandomVariable* rv = (*it);
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    

    if ( ppf.find(rvp) == ppf.end() ) {
      coredump("ERROR: can't find assigned rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
    }
    RandomVariable* nrv = newRvs[ppf[rvp]];
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
MaxClique::makeComplete(const set<RandomVariable*> &rvs)
{
  // just go through each rv and union its neighbors
  // set with all of rvs.

  for (set<RandomVariable*>::iterator i=rvs.begin();
       i != rvs.end(); i++) {
    set<RandomVariable*> res;
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
computeWeight(const set<RandomVariable*>& nodes,
	      const RandomVariable* node,
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
    if (node->discrete && node->hidden) {
      DiscreteRandomVariable *const drv = (DiscreteRandomVariable*)node;
      // weight changes only if node is not deterministic (Lauritzen CG inference).
      if (useDeterminism && drv->sparse()) {
	// then there is a possibility that this node does not affect
	// the state space, as long as all of this nodes parents are
	// in the clique.  The variable 'truly_sparse' indicates that.
	bool truly_sparse = true;
	for (unsigned i=0;i<drv->allPossibleParents.size();i++) {
	  if (nodes.find(drv->allPossibleParents[i]) == nodes.end()) {
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
    } else if (!node->discrete) {
      // node is continuous observed.
      ContinuousRandomVariable *crv = (ContinuousRandomVariable*)node;
      tmp_weight += crv->dimensionality()*continuousObservationPerFeaturePenalty;
    }
  }
  // Next, get weight of all 'nodes'
  for (set<RandomVariable*>::iterator j=nodes.begin();
       j != nodes.end();
       j++) {
    RandomVariable *const rv = (*j);
    // First get cardinality of 'node', but if
    // it is continuous or observed, it does not change the weight.
    // TODO: The assumption here (for now) is that all continuous variables
    // are observed. This will change in a future version (Lauritzen CG inference).
    if (rv->discrete && rv->hidden) {
      DiscreteRandomVariable *const drv = (DiscreteRandomVariable*)rv;
      if (useDeterminism && drv->sparse()) {
	// then there is a possibility that this node does not affect
	// the state space, as long as all of this nodes parents are
	// in the clique.  The variable 'truly_sparse' indicates that.
	bool truly_sparse = true;
	for (unsigned i=0;i<drv->allPossibleParents.size();i++) {
	  if ((nodes.find(drv->allPossibleParents[i]) == nodes.end())
	      &&
	      (drv->allPossibleParents[i] != node)) {
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
    } else if (!rv->discrete) {
      // node is continuous observed.
      ContinuousRandomVariable *crv = (ContinuousRandomVariable*)rv;
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
computeWeightWithExclusion(const set<RandomVariable*>& nodes,
			   const set<RandomVariable*>& unassignedIteratedNodes,
			   const set<RandomVariable*>& unionSepNodes,
			   const bool useDeterminism)
{
  // compute weight in log base 10 so that
  //   1) we don't overflow
  //   2) base 10 is an easy to understand magnitude rating of state space.

  float tmp_weight = 0;
  // Next, get weight of all 'nodes'
  for (set<RandomVariable*>::iterator j=nodes.begin();
       j != nodes.end();
       j++) {
    RandomVariable *const rv = (*j);
    // First get cardinality of 'node', but if
    // it is continuous or observed, it does not change the weight.
    // TODO: The assumption here (for now) is that all continuous variables
    // are observed. This will change in a future version (Lauritzen CG inference).
    if (rv->discrete && rv->hidden) {
      DiscreteRandomVariable *const drv = (DiscreteRandomVariable*)rv;
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
	  for (unsigned i=0;i<drv->allPossibleParents.size();i++) {
	    if (nodes.find(drv->allPossibleParents[i]) == nodes.end()) {
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
    } else if (!rv->discrete) {
      // node is continuous observed.
      ContinuousRandomVariable *crv = (ContinuousRandomVariable*)rv;
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
    assigned_nodes: The set of nodes that are assigned (i.e., used
                    to compute a probability) to this clique.
    unassigned_iterated_nodes: nodes in this clique that   
                     are not assigned and are not incomming
                     separator nodes (so they need to be iterated in
                     full, based in their cardinality).
    separator_nodes: Nodes that are part of an incomming separator
                     during the collect evidence stage of inference. 
                     The cost of these nodes (if sparse) depends on two things.

                     If the separator nodes are also assigned here, we pay
                     full cardinality (since we're doing separator driven
                     clique instantiation).
                  
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
MaxClique::
computeWeightInJunctionTree(const set<RandomVariable*>& nodes,
			    const set<RandomVariable*>& assignedNodes,
			    const set<RandomVariable*>& cumulativeAssignedNodes,
			    const set<RandomVariable*>& unassignedIteratedNodes,
			    const set<RandomVariable*>& cumulativeUnassignedIteratedNodes,
			    const set<RandomVariable*>& separatorNodes,
			    const set<RandomVariable*>& unassignedInPartition,
			    const bool upperBound,
			    const bool useDeterminism)
{
  // compute weight in log base 10 so that
  //   1) we don't overflow
  //   2) base 10 is an easy to understand magnitude rating of state space.

  float tmp_weight = 0;
  // Next, get weight of all 'nodes'
  for (set<RandomVariable*>::iterator j=nodes.begin();
       j != nodes.end();
       j++) {
    RandomVariable *const rv = (*j);
    // First get cardinality of 'node', but if
    // it is continuous or observed, it does not change the weight.
    // TODO: The assumption here (for now) is that all continuous variables
    // are observed. This will change in a future version (Lauritzen CG inference).
    if (rv->discrete && rv->hidden) {
      DiscreteRandomVariable *const drv = (DiscreteRandomVariable*)rv;
      // see comments above for description and rational of this algorithm
      if (!useDeterminism || !drv->sparse()) {
	// charge full amount since not sparse.
	tmp_weight += log10((double)drv->cardinality);
      } else if (separatorNodes.find(rv) != separatorNodes.end()) {
	// separator node case.
	if (unassignedInPartition.find(rv) != unassignedInPartition.end()) {
	  // separator, unassigned in the current partition:

	  if (upperBound)
	    tmp_weight += log10((double)drv->cardinality);
	  else {
	    // Then unasigned in this partition. We assume that the node
	    // has been assigned in another (say previous) partition and
	    // we don't charge full amount. Note that this could cause
	    // the estimate to be LOWER than the true weight.
	    tmp_weight += log10((double)drv->useCardinality());
	    // cludge/hack:
	    // (min(card,unavailble_parents_prod_card) + use_card)/2 (but in log domain).
	    // tmp_weight += (min(rv->log10ProductCardOfParentsNotContainedInSet(separatorNodes),
	    // log10((double)drv->cardinality)) + log10((double)drv->useCardinality()))/2.0;
	    // tmp_weight += log10((double)drv->useCardinality());
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
	    // Charge low amount since it has been assigned in some
	    // previous clique, and at least one of the separators will
	    // kill off the zero prob entries. This could cause the
	    // estimate to be LOWER than the true weight.
	    tmp_weight += log10((double)drv->useCardinality());
	  }
	}
      } else {
	// non separator node case.
	if (unassignedInPartition.find(rv) != unassignedInPartition.end()) {
	  // non separator, unassigned in partition.

	  if (upperBound) 
	    tmp_weight += log10((double)drv->cardinality);
	  else {
	    // Then unasigned in this partition. We assume that the
	    // node has been assigned in another (say previous) partition
	    // and we don't charge full amount. Note that
	    // this could cause the estimate to be LOWER than
	    // the true weight.
	    tmp_weight += log10((double)drv->useCardinality());
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
    } else if (!rv->discrete) {
      // node is continuous observed.
      ContinuousRandomVariable *crv = (ContinuousRandomVariable*)rv;
      tmp_weight += crv->dimensionality()*continuousObservationPerFeaturePenalty;
    } 
  }
  return tmp_weight;
}


#if 0
float
MaxClique::
computeWeightInJunctionTree(const set<RandomVariable*>& nodes,
			    const set<RandomVariable*>& assignedNodes,
			    const set<RandomVariable*>& unassignedIteratedNodes,
			    const set<RandomVariable*>& separatorNodes,
			    const set<RandomVariable*>& cumulativeAssignedNodes,
			    const bool useDeterminism)
{
  // compute weight in log base 10 so that
  //   1) we don't overflow
  //   2) base 10 is an easy to understand magnitude rating of state space.

  float tmp_weight = 0;
  // Next, get weight of all 'nodes'
  for (set<RandomVariable*>::iterator j=nodes.begin();
       j != nodes.end();
       j++) {
    RandomVariable *const rv = (*j);
    // First get cardinality of 'node', but if
    // it is continuous or observed, it does not change the weight.
    // TODO: The assumption here (for now) is that all continuous variables
    // are observed. This will change in a future version (Lauritzen CG inference).
    if (rv->discrete && rv->hidden) {
      DiscreteRandomVariable *const drv = (DiscreteRandomVariable*)rv;
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
	      // @@@@ This is the problem case here. we need to estimate the
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
  set<RandomVariable*>::iterator it;
  for (it = nodes.begin();
       it != nodes.end();
       it++) {
    RandomVariable* rv = (*it);
    if (rv->hidden)
      hiddenNodes.push_back(rv);
  }

  // setup and re-construct packer
  assert (packer.unPackedLen() == 0); // make sure it is empty.
  new (&packer) PackCliqueValue(hiddenNodes);

  // ensure that we have something to store.
  assert (packer.packedLen() > 0);

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

  // @@@ set to small size now to test out re-allocation schemes.
  allocationUnitChunkSize = 1;
  // allocationUnitChunkSize = 10000;

  if (packer.packedLen() > IMC_NWWOH) {
    // setup value hodler
    new (&valueHolder) CliqueValueHolder(packer.packedLen(),
					 allocationUnitChunkSize,
					 1.25);
    // set up common clique hash tables 
    // TODO: add appropriate default staring hash sizes.
    // new (&cliqueValueHashSet) vhash_set< unsigned > (packer.packedLen(),2);
    new (&cliqueValueHashSet) vhash_set< unsigned > (packer.packedLen(),2); // 10000
  } else {
    // then no need to do a hash table at all, just store the packed
    // values in a local integer.
  }
}


/*-
 *-----------------------------------------------------------------------
 * MaxClique::computeAssignedNodesToIterate()
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
    RandomVariable*rv = sortedAssignedNodes[i];
    if (!rv->hidden) {
      // observed
      if (assignedProbNodes.find(rv) != assignedProbNodes.end()) {
	dispositionSortedAssignedNodes[i] = AN_CPT_ITERATION_COMPUTE_AND_APPLY_PROB;
      } else {
	if (rv->discrete) {
	  // discrete observed variable. Prune here if sparse.
	  DiscreteRandomVariable* drv = (DiscreteRandomVariable*)rv;	  
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
      DiscreteRandomVariable* drv = (DiscreteRandomVariable*)rv;
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
  


#if 0


  set<RandomVariable*> res;
  set_intersection(assignedNodes.begin(),assignedNodes.end(),
		   cumulativeUnassignedIteratedNodes.begin(),cumulativeUnassignedIteratedNodes.end(),
		   inserter(res,res.end()));
  if (res.size() == 0)
    return;

  // too bad, we've got assigned nodes in this clique that
  // we do not iterate, so we must add a check in iteration.

  iterateSortedAssignedNodesP.resize(sortedAssignedNodes.size());
  for (unsigned i=0;i<sortedAssignedNodes.size();i++) {
    if ((cumulativeUnassignedIteratedNodes.find(sortedAssignedNodes[i]) != cumulativeUnassignedIteratedNodes.end())) {
      // found, no need to iterate (which is bad since this means that
      // we'll be iterating this assigned node via the separator).
      iterateSortedAssignedNodesP[i] = false;
    }  else {
      // not found, need to iterate (which is good)
      iterateSortedAssignedNodesP[i] = true;
    }
  }
#endif
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
MaxClique::printAllJTInfo(FILE*f,const unsigned indent,const set<RandomVariable*>& unassignedInPartition)
{

  // TODO: also print out nubmer of bits for acc and rem.

  psp(f,indent*2);
  fprintf(f,"Clique information: %d packed bits, %d unsigned words, weight = %f, jt_weight = %f\n",
	  packer.packedLenBits(),packer.packedLen(),weight(),weightInJunctionTree(unassignedInPartition));



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
				       vector <RandomVariable*>& newRvs,
				       map < RVInfo::rvParent, unsigned >& ppf,
				       const unsigned int frameDelta)
  : origin(from_clique)
{

  set<RandomVariable*>::iterator it;

  // clone over nodes RVs.
  fNodes.resize(from_clique.nodes.size());
  unsigned i=0;
  for (it = from_clique.nodes.begin();
       it != from_clique.nodes.end();
       it++) {
    RandomVariable* rv = (*it);
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    

    if ( ppf.find(rvp) == ppf.end() ) {
      coredump("ERROR: can't find rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
    }

    RandomVariable* nrv = newRvs[ppf[rvp]];
    fNodes[i++] = nrv;
  }

  // and clone over assigned nodes and sorted assigned nodes
  fSortedAssignedNodes.resize(from_clique.sortedAssignedNodes.size());
  for (i=0;i<from_clique.sortedAssignedNodes.size();i++) {
    RandomVariable* rv = from_clique.sortedAssignedNodes[i];
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    

    if ( ppf.find(rvp) == ppf.end() ) {
      coredump("ERROR: can't find assigned rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
    }

    RandomVariable* nrv = newRvs[ppf[rvp]];
    fSortedAssignedNodes[i] = nrv;
  }

  // do unassignedIteratedNodes
  i=0;
  fUnassignedIteratedNodes.resize(from_clique.unassignedIteratedNodes.size());
  for (it = from_clique.unassignedIteratedNodes.begin();
       it != from_clique.unassignedIteratedNodes.end();
       it++) {
    RandomVariable* rv = (*it);
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    

    if ( ppf.find(rvp) == ppf.end() ) {
      coredump("ERROR: can't find assigned rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
    }
    RandomVariable* nrv = newRvs[ppf[rvp]];
    fUnassignedIteratedNodes[i++] = nrv;
  }


  // do unassignedNodes
  if (MaxClique::ceSeparatorDrivenInference == false) {
    i=0;
    fUnassignedNodes.resize(from_clique.unassignedNodes.size());
    for (it = from_clique.unassignedNodes.begin();
	 it != from_clique.unassignedNodes.end();
	 it++) {
      RandomVariable* rv = (*it);
      RVInfo::rvParent rvp;
      rvp.first = rv->name();
      rvp.second = rv->frame()+frameDelta;    

      if ( ppf.find(rvp) == ppf.end() ) {
	coredump("ERROR: can't find assigned rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
		 rv->name().c_str(),rv->frame(),frameDelta,
		 rvp.first.c_str(),rvp.second);
      }
      RandomVariable* nrv = newRvs[ppf[rvp]];
      fUnassignedNodes[i++] = nrv;
    }
  }


  // Clique values only store/hash values of hidden (thus necessarily
  // discrete) variables since they are the only thing that change.
  discreteValuePtrs.resize(from_clique.hiddenNodes.size());
  for (i=0;i<discreteValuePtrs.size();i++) {
    // get the unrolled rv for this hidden node
    RandomVariable* rv = from_clique.hiddenNodes[i];
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    
    if ( ppf.find(rvp) == ppf.end() ) {
      warning("ERROR: can't find assigned rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
      assert ( ppf.find(rvp) != ppf.end() );
    }
    RandomVariable* nrv = newRvs[ppf[rvp]];

    // hidden nodes are always discrete (in this version).
    DiscreteRandomVariable* drv = 
      (DiscreteRandomVariable*)nrv;
    // grab a pointer directly to its value for easy access later.
    discreteValuePtrs[i] = &(drv->val);
  }

  numCliqueValuesUsed = 0;
  maxCEValue.set_to_zero();

  // TODO: optimize this.
  cliqueValues.resize(3); // 10000

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

  // if we're still here, we do separator driven inference.
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
 *   All the cliques incomming separators *must* have been created and are ready
 *   to be used to produce the clique table.
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
    // TODO: check this, as this condition might fail if
    // we've completely pruned away parent clique.
    assert ( sep.separatorValues.size() == 1);
    assert ( sep.origin.hRemainder.size() > 0 );
    sepValueNumber = 0;
  }

  // Iterate through remainder of separator Q: should we do some
  // pruning here as well?
  if (sep.origin.hRemainder.size() == 0) {
    // only one remainder entry (in position 0) and also no need to
    // unpack since all has been covered by accumulated intersection
    // set above in a previous separator. Just continue on with single value.

    if (message(Giga)) {
      fprintf(stdout,"Separator iteration nounpack, sepNumber =%d, part sepNo = %d,p = %f, nodes:",
	      sepNumber,origin.ceReceiveSeparators[sepNumber],p.val());
      printRVSetAndValues(stdout,sep.fNodes);
    }

    // we should have one entry only.
    assert ( sep.separatorValues.ptr[sepValueNumber].remValues.size() == 1 );

    // Continue down with new probability value.
    // Search for tag 'ALLOCATE_REMVALUES_OPTION' in this file for
    // more info why remValues.ptr[0] exists.
    // TODO: separator prune here.
    ceIterateSeparators(part,sepNumber+1,
			p*
			sep.separatorValues.ptr[sepValueNumber].remValues.ptr[0].p);

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
	  fprintf(stdout,"Separator iteration %d, sepNumber =%d, part sepNo = %d,p = %f, nodes:",
		  i,
		  sepNumber,origin.ceReceiveSeparators[sepNumber],p.val());
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
	  fprintf(stdout,"pSeparator iteration %d, sepNumber =%d, part sepNo = %d,p = %f, nodes:",
		  i,
		  sepNumber,origin.ceReceiveSeparators[sepNumber],p.val());
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
 * MaxClique::ceIterateAssignedNodes()
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
 *   same for the unassigned nodes (if any).
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
InferenceMaxClique::ceIterateAssignedNodes(JT_InferencePartition& part,
					   const unsigned nodeNumber,
					   const logpr p)
{

  if (nodeNumber == fSortedAssignedNodes.size()) {
    // time to store clique value and total probability, p is
    // current clique probability.

    // printf("ceIterateAssignedNodes: nodeNumber = %d, p = %f,",nodeNumber,p.val());
    // printRVSet(fNodes);

    // keep track of the max clique probability right here.
    if (p > maxCEValue)
      maxCEValue = p;

    if (numCliqueValuesUsed >= cliqueValues.size()) {
      // TODO: optimize this.
      cliqueValues.resizeAndCopy(cliqueValues.size()*2);
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

    return;
  }
  RandomVariable* rv = fSortedAssignedNodes[nodeNumber];
  // do the loop right here

  infoMsg(Giga,"Starting assigned iteration of rv %s(%d), nodeNumber=%d, p = %f\n",
	  rv->name().c_str(),rv->frame(),nodeNumber,p.val());

  switch (origin.dispositionSortedAssignedNodes[nodeNumber]) {
  case MaxClique::AN_CPT_ITERATION_COMPUTE_AND_APPLY_PROB: 
    {
      rv->clampFirstValue();
      do {
	// At each step, we compute probability
	logpr cur_p = rv->probGivenParents();
	if (message(Giga)) {
	  if (!rv->discrete) {
	    infoMsg(Giga,"  Assigned iteration and prob application of rv %s(%d)=C, nodeNumber =%d, p = %f\n",
		    rv->name().c_str(),rv->frame(),nodeNumber,p.val());
	  } else {
	    // DiscreteRandomVariable* drv = (DiscreteRandomVariable*)rv;
	    infoMsg(Giga,"  Assigned CPT iteration and prob application of rv %s(%d)=%d, nodeNumber =%d, p = %f\n",
		    rv->name().c_str(),rv->frame(),rv->val,nodeNumber,p.val());
	  }
	}
	// if at any step, we get zero, then back out.
	if (!cur_p.essentially_zero()) {
	  // Continue, updating probability by cur_p, contributing this
	  // probability to the clique potential.
	  ceIterateAssignedNodes(part,nodeNumber+1,p*cur_p);
	}
      } while (rv->clampNextValue());
    }
    break;

  case MaxClique::AN_CPT_ITERATION_COMPUTE_PROB_REMOVE_ZEROS:
    {
      rv->clampFirstValue();
      do {
	// At each step, we compute probability
	logpr cur_p = rv->probGivenParents();
	if (message(Giga)) {
	  if (!rv->discrete) {
	    infoMsg(Giga,"  Assigned iteration of rv %s(%d)=C, nodeNumber =%d, p = %f\n",
		    rv->name().c_str(),rv->frame(),nodeNumber,p.val());
	  } else {
	    // DiscreteRandomVariable* drv = (DiscreteRandomVariable*)rv;
	    infoMsg(Giga,"  Assigned CPT iteration and zero removal of rv %s(%d)=%d, nodeNumber =%d, p = %f\n",
		    rv->name().c_str(),rv->frame(),rv->val,nodeNumber,p.val());
	  }
	}
	// if at any step, we get zero, then back out.
	if (!cur_p.essentially_zero()) {
	  // Continue, do not update probability!!
	  ceIterateAssignedNodes(part,nodeNumber+1,p);
	}
      } while (rv->clampNextValue());
    }
    break;


  case MaxClique::AN_CARD_ITERATION:
    {
      DiscreteRandomVariable* drv = (DiscreteRandomVariable*)rv;
      // do the loop right here
      drv->val = 0;
      do {
	infoMsg(Giga,"  Assigned card iteration of rv %s(%d)=%d, nodeNumber = %d, p = %f\n",
		rv->name().c_str(),rv->frame(),rv->val,nodeNumber,p.val());
	// Continue, do not update probability!!
	ceIterateAssignedNodes(part,nodeNumber+1,p);
      } while (++drv->val < drv->cardinality);
    }
    break;

  case MaxClique::AN_COMPUTE_AND_APPLY_PROB:
    {
      // TODO: Make more efficient version of this, based on the type of
      // RV.
      logpr cur_p = rv->probGivenParentsWSetup();
      // if at any step, we get zero, then back out.
      if (!cur_p.essentially_zero()) {
	// Continue, updating probability by cur_p.
	ceIterateAssignedNodes(part,nodeNumber+1,p*cur_p);
      }
    }
    break;

  case MaxClique::AN_CONTINUE:
    ceIterateAssignedNodes(part,nodeNumber+1,p);
    break;

  case MaxClique::AN_CONTINUE_COMPUTE_PROB_REMOVE_ZEROS:
    {
      // TODO: Make more efficient version of this, based on the type of
      // RV.
      logpr cur_p = rv->probGivenParentsWSetup();
      if (!cur_p.essentially_zero()) {
	// Continue, do not update probability!!
	ceIterateAssignedNodes(part,nodeNumber+1,p);
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
  RandomVariable* rv = fUnassignedIteratedNodes[nodeNumber];
  infoMsg(Giga,"Starting Unassigned iteration of rv %s(%d), nodeNumber = %d, p = %f\n",
	  rv->name().c_str(),rv->frame(),nodeNumber,p.val());

  if (rv->hidden) {
    // only discrete RVs can be hidden for now.
    DiscreteRandomVariable* drv = (DiscreteRandomVariable*)rv;
    // do the loop right here
    drv->val = 0;
    do {
      infoMsg(Giga,"  Unassigned iteration of rv %s(%d)=%d, nodeNumber = %d, p = %f\n",
	      rv->name().c_str(),rv->frame(),rv->val,nodeNumber,p.val());
      // continue on, effectively multiplying p by unity.
      ceIterateUnassignedIteratedNodes(part,nodeNumber+1,p);
    } while (++drv->val < drv->cardinality);
  } else {
    // observed, either discrete or continuous
    if (rv->discrete) {
      DiscreteRandomVariable* drv = (DiscreteRandomVariable*)rv;
      // TODO: for observed variables, do this once at the begining
      // before any looping here.
      drv->setToObservedValue();
      infoMsg(Giga,"  Unassigned iteration of rv %s(%d)=%d, nodeNumber = %d, p = %f\n",
	      rv->name().c_str(),rv->frame(),rv->val,nodeNumber,p.val());
    } else {
      // nothing to do since we get continuous observed
      // value indirectly
      infoMsg(Giga,"  Unassigned iteration of rv %s(%d)=C, nodeNumber = %d, p = %f\n",
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

  // create an ininitialized variable using dummy argument
  logpr beamThreshold((void*)0);

  if (origin.cliqueBeam != (-LZERO)) {
    // then we do clique table pruning right here.
    // break into the logp to avoid unnecessary zero checking.
    beamThreshold.valref() = maxCEValue.valref() - origin.cliqueBeam;
  } else {
    // set beam threshold to a value that will never cause pruning.
    beamThreshold.set_to_zero();
  }

  const unsigned origNumCliqueValuesUsed = numCliqueValuesUsed;  
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
    if (origin.packer.packedLen() <= IMC_NWWOH) {
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
    if (sep.origin.hAccumulatedIntersection.size() > 0) {
      // an accumulated intersection exists.

      // make sure there is at least one available entry
      assert ( sep.numSeparatorValuesUsed <= sep.separatorValues.size());
      if (sep.numSeparatorValuesUsed >= sep.separatorValues.size()) {
	const unsigned old_size = sep.separatorValues.size();
	// TODO: optimize this size re-allocation.
	sep.separatorValues.resizeAndCopy(sep.separatorValues.size()*2);
	if (sep.origin.accPacker.packedLen() <= ISC_NWWOH_AI) {
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
	if (sep.remDiscreteValuePtrs.size() > 0) {
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
      if (sep.origin.accPacker.packedLen() <= ISC_NWWOH_AI) {
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
      if (sep.remDiscreteValuePtrs.size() == 0) {
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
	} else {
	  // already there so must have hit before.
	  // we thus accumulate.
	  sv.remValues.ptr[0].p += cliqueValues.ptr[cvn].p;
	}

	goto next_iteration;
      }

    } else {
      accIndex = 0;
    }

    // if we're here, then we must have some remainder
    // pointers.
    // TODO: remove assertion when debugged.
    assert (sep.remDiscreteValuePtrs.size() > 0);

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
    
    // make sure there is at least one available entry
    assert (sv.numRemValuesUsed <= sv.remValues.size());
    if (sv.numRemValuesUsed >= sv.remValues.size()) {
      // TODO: optimize this growth rate.
      // start small but grow fast.
      sv.remValues.resizeAndCopy(1+sv.remValues.size()*2); // *3
      if (sep.origin.remPacker.packedLen() <= ISC_NWWOH_RM) {
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
    if (sep.origin.remPacker.packedLen() <= ISC_NWWOH_RM) {
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
    sv.remValues.ptr[*remIndexp].p += cliqueValues.ptr[cvn].p;

  }
  next_iteration:
  cvn++;
  }

  
  // To reallocate or not to reallocate, that is the question.  here,
  // we just reallocate for now.
  // 
  // TODO: reallocate only if change is > some percentage (say 5%),
  // and export to command line.
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
 *    same beam width).
 *
 * Preconditions:
 *   1) the value of the max clique 'maxCEValue' must have been
 *      computed already.
 *
 *   2) clique table must be created, meaning that either:
 *
 *        InferenceMaxClique::ceIterateAssignedNodes()
 *   or
 *        InferenceMaxClique::ceIterateAssignedNodesCliqueDriven
 *   must have just been called.
 *
 * Postconditions:
 *    Clique table has been pruned, and memory for it has been re-allocated to
 *    fit the smaller size.
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



//////////////
// Clique driven version of gather from incomming separators
/////////////

void
InferenceMaxClique::ceGatherFromIncommingSeparatorsCliqueDriven(JT_InferencePartition& part)
{
  assert (MaxClique::ceSeparatorDrivenInference == false); 
  logpr p = 1.0;
  if (origin.unassignedNodes.size() == 0) {
    ceIterateAssignedNodesCliqueDriven(part,0,p);
  } else {
    ceIterateUnassignedNodesCliqueDriven(part,0,p);
  }
}



void
InferenceMaxClique::ceIterateAssignedNodesCliqueDriven(JT_InferencePartition& part,
						       const unsigned nodeNumber,
						       logpr p)
{
  if (nodeNumber == fSortedAssignedNodes.size()) {
    // time to store clique value and total probability, p is
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

      // if we're here, then we must have some remainder
      // pointers.
      // TODO: remove assertion when debugged.
      assert (sep.remDiscreteValuePtrs.size() > 0);

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
    }


    // keep track of the max clique probability right here.
    if (p > maxCEValue)
      maxCEValue = p;

    if (numCliqueValuesUsed >= cliqueValues.size()) {
      // TODO: optimize this.
      cliqueValues.resizeAndCopy(cliqueValues.size()*2);
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

    return;
  }
  RandomVariable* rv = fSortedAssignedNodes[nodeNumber];
  // do the loop right here

  infoMsg(Giga,"Starting assigned iteration of rv %s(%d), nodeNumber=%d, p = %f\n",
	  rv->name().c_str(),rv->frame(),nodeNumber,p.val());

  switch (origin.dispositionSortedAssignedNodes[nodeNumber]) {
  case MaxClique::AN_NOTSEP_PROB_SPARSEDENSE:
  case MaxClique::AN_SEP_PROB_SPARSEDENSE:
    {
      rv->clampFirstValue();
      do {
	// At each step, we compute probability
	logpr cur_p = rv->probGivenParents();
	if (message(Giga)) {
	  if (!rv->discrete) {
	    infoMsg(Giga,"  Assigned iteration and prob application of rv %s(%d)=C, nodeNumber =%d, p = %f\n",
		    rv->name().c_str(),rv->frame(),nodeNumber,p.val());
	  } else {
	    // DiscreteRandomVariable* drv = (DiscreteRandomVariable*)rv;
	    infoMsg(Giga,"  Assigned CPT iteration and prob application of rv %s(%d)=%d, nodeNumber =%d, p = %f\n",
		    rv->name().c_str(),rv->frame(),rv->val,nodeNumber,p.val());
	  }
	}
	// if at any step, we get zero, then back out.
	if (!cur_p.essentially_zero()) {
	  // Continue, updating probability by cur_p, contributing this
	  // probability to the clique potential.
	  ceIterateAssignedNodesCliqueDriven(part,nodeNumber+1,p*cur_p);
	}
      } while (rv->clampNextValue());
    }
    break;


  case MaxClique::AN_NOTSEP_NOTPROB_DENSE:
    {
      DiscreteRandomVariable* drv = (DiscreteRandomVariable*)rv;
      // do the loop right here
      drv->val = 0;
      do {
	infoMsg(Giga,"  Assigned card iteration of rv %s(%d)=%d, nodeNumber = %d, p = %f\n",
		rv->name().c_str(),rv->frame(),rv->val,nodeNumber,p.val());
	// Continue, do not update probability!!
	ceIterateAssignedNodesCliqueDriven(part,nodeNumber+1,p);
      } while (++drv->val < drv->cardinality);
    }
    break;


  default: // all other cases, we do CPT iteration removing zeros without updating probabilities.
    {
      rv->clampFirstValue();
      do {
	// At each step, we compute probability
	logpr cur_p = rv->probGivenParents();
	if (message(Giga)) {
	  if (!rv->discrete) {
	    infoMsg(Giga,"  Assigned iteration of rv %s(%d)=C, nodeNumber =%d, p = %f\n",
		    rv->name().c_str(),rv->frame(),nodeNumber,p.val());
	  } else {
	    // DiscreteRandomVariable* drv = (DiscreteRandomVariable*)rv;
	    infoMsg(Giga,"  Assigned CPT iteration and zero removal of rv %s(%d)=%d, nodeNumber =%d, p = %f\n",
		    rv->name().c_str(),rv->frame(),rv->val,nodeNumber,p.val());
	  }
	}
	// if at any step, we get zero, then back out.
	if (!cur_p.essentially_zero()) {
	  // Continue, do not update probability!!
	  ceIterateAssignedNodesCliqueDriven(part,nodeNumber+1,p);
	}
      } while (rv->clampNextValue());
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
  RandomVariable* rv = fUnassignedNodes[nodeNumber];
  infoMsg(Giga,"Starting Unassigned iteration of rv %s(%d), nodeNumber = %d, p = %f\n",
	  rv->name().c_str(),rv->frame(),nodeNumber,p.val());

  if (rv->hidden) {
    // only discrete RVs can be hidden for now.
    DiscreteRandomVariable* drv = (DiscreteRandomVariable*)rv;
    // do the loop right here
    drv->val = 0;
    do {
      infoMsg(Giga,"  Unassigned iteration of rv %s(%d)=%d, nodeNumber = %d, p = %f\n",
	      rv->name().c_str(),rv->frame(),rv->val,nodeNumber,p.val());
      // continue on, effectively multiplying p by unity.
      ceIterateUnassignedNodesCliqueDriven(part,nodeNumber+1,p);
    } while (++drv->val < drv->cardinality);
  } else {
    // observed, either discrete or continuous
    if (rv->discrete) {
      DiscreteRandomVariable* drv = (DiscreteRandomVariable*)rv;
      // TODO: for observed variables, do this once at the begining
      // before any looping here.
      drv->setToObservedValue();
      infoMsg(Giga,"  Unassigned iteration of rv %s(%d)=%d, nodeNumber = %d, p = %f\n",
	      rv->name().c_str(),rv->frame(),rv->val,nodeNumber,p.val());
    } else {
      // nothing to do since we get continuous observed
      // value indirectly
      infoMsg(Giga,"  Unassigned iteration of rv %s(%d)=C, nodeNumber = %d, p = %f\n",
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
  for (unsigned i=0;i<numCliqueValuesUsed;i++) {
    p += cliqueValues.ptr[i].p;
  }
  return p;
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

  for (unsigned cvn=0;cvn<numCliqueValuesUsed;cvn++) {{
    // TODO: beam pruning
    // Prune away based on forward computed beam
    //if (cliqueValues.ptr[cvn].p < beamThreshold) {
    //     continue;
    //}


    // TODO: optimize away this conditional check. (and/or use const
    // local variable to indicate it wont change)
    if (origin.packer.packedLen() <= IMC_NWWOH) {
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
	cliqueValues.ptr[cvn].p *= sv.remValues.ptr[0].bp;
	// done
	goto next_iteration;
      }
    } else {
      // no accumulated intersection exists, everything
      // is in the remainder.
      accIndex = 0;
    }

    // if we're here, then we must have some remainder
    // pointers.
    // TODO: remove assertion when debugged.
    assert (sep.remDiscreteValuePtrs.size() > 0);

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

    // We've finally got the sep entry. Multiply it it into the
    // current clique value.
    cliqueValues.ptr[cvn].p *= sv.remValues.ptr[*remIndexp].bp;

  }
  next_iteration:
  ;    
  }
}


// Scatter out to the outgoing separators which are the same as the
// receive separators in the collect evidence stage, so we use that
// array here.
void 
InferenceMaxClique::
deScatterToOutgoingSeparators(JT_InferencePartition& part)
{
  if (origin.ceReceiveSeparators.size() == 0)
    return;

  // Note. All separator .bp values have already been initialized to
  // zero when the structure containing them 'a RemainderValue' was
  // constructed. All memory reallocations will have preserved these
  // initializations, so there is no need to scan through initializing
  // bp to zero here.

  
  // allocate some temporary storage for packed separator values.
  // 128 words is *much* bigger than any possible packed clique value
  // will take on, but it is easy/fast to allocate on the stack right now.
  unsigned packedVal[128];
 
  for (unsigned cvn=0;cvn<numCliqueValuesUsed;cvn++) {
    // TODO: beam pruning
    // Prune away based on forward computed beam
    //if (cliqueValues[cvn].p < beamThreshold) {
    //     continue;
    //}


    // TODO: optimize away this conditional check. (and/or use const
    // local variable to indicate it wont change)
    if (origin.packer.packedLen() <= IMC_NWWOH) {
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
	  sv.remValues.ptr[0].bp += cliqueValues.ptr[cvn].p;
	  // done, move on to next separator.
	  continue; 
	}
      } else {
	// no accumulated intersection exists, everything
	// is in the remainder.
	accIndex = 0;
      }

      // if we're here, then we must have some remainder
      // pointers.
      // TODO: remove assertion when debugged.
      assert (sep.remDiscreteValuePtrs.size() > 0);

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
      sv.remValues.ptr[*remIndexp].bp += cliqueValues.ptr[cvn].p;
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
	// Do direct value reference subtraction in log domain
	// (corresponding to divison in original domain) to ensure
	// that compiler creates no temporaries. In other words, this
	// operation could be "rv->bp = rv->bp / rv->p;"
	//    rv->bp.valref() = rv->bp.valref() - rv->p.valref();
	// Do slower version for now until debugged:
	// We assume here that (!rv->p.zero()) is true since
	// we pruned all zero p's above. If we didn't prune,
	// then rv->p == zero would imply that rv->bp == zero,
	// and we would need to do a check. Note that this
	// pruning always occurs, regardless of beam.
	rv->bp /= rv->p;
      }
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
{
  nodes.clear();

  // create a set of nodes that is the intersection of the two
  set_intersection(c1.nodes.begin(),c1.nodes.end(),
		   c2.nodes.begin(),c2.nodes.end(),
		   inserter(nodes,nodes.end()));
  assert (nodes.size() > 0);
  
}



SeparatorClique::SeparatorClique(SeparatorClique& from_sep,
				 vector <RandomVariable*>& newRvs,
				 map < RVInfo::rvParent, unsigned >& ppf,
				 const unsigned int frameDelta)
{
  set<RandomVariable*>::iterator it;

  // clone over nodes RVs and accumulated intersection.
  for (it = from_sep.nodes.begin();
       it != from_sep.nodes.end();
       it++) {
    RandomVariable* rv = (*it);
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    

    if ( ppf.find(rvp) == ppf.end() ) {
      coredump("INTERNAL ERROR: can't find rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
    }

    RandomVariable* nrv = newRvs[ppf[rvp]];
    nodes.insert(nrv);
  }

  // copy over accumulated intersection
  for (it=from_sep.accumulatedIntersection.begin();
       it != from_sep.accumulatedIntersection.end();
       it++) {
    RandomVariable* rv = (*it);
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    


    if ( ppf.find(rvp) == ppf.end() ) {
      coredump("INTERNAL ERROR: can't find rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
    }

    RandomVariable* nrv = newRvs[ppf[rvp]];
    accumulatedIntersection.insert(nrv);
  }

  // and 'remainder'
  for (it=from_sep.remainder.begin();
       it != from_sep.remainder.end();
       it++) {
    RandomVariable* rv = (*it);
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    

    if ( ppf.find(rvp) == ppf.end() ) {
      coredump("INTERNAL ERROR: can't find rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
    }

    RandomVariable* nrv = newRvs[ppf[rvp]];
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
  set<RandomVariable*>::iterator it;
  for (it = accumulatedIntersection.begin();
       it != accumulatedIntersection.end();
       it++) {
    RandomVariable* rv = (*it);
    if (rv->hidden)
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
					      2, // 5000
					      1.25);
      // TODO: optimize starting size.
      new (&accSepValHashSet) vhash_set< unsigned > (accPacker.packedLen(),2); // 5000
    }
  }

  // create a vector form of the variables.
  for (it = remainder.begin();
       it != remainder.end();
       it++) {
    RandomVariable* rv = (*it);
    if (rv->hidden)
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
					      2, // 5000
					      1.25);
      // TODO: optimize starting size
      new (&remSepValHashSet) vhash_set< unsigned > (remPacker.packedLen(),2); // 5000
    }
  }

  // make sure we have at least one hidden variable in separator.
  if (hAccumulatedIntersection.size() == 0 && hRemainder.size() == 0) {
    // We should never have a separator clique consisting of only and
    // entirely observed variables. If we do, it means that one part
    // of the JT is entirely cut off from the rest since cliques are
    // independent of each other given sep sets.  TODO: put an error
    // message in JT creation code if this case occurs.
    warning("ERROR: separator clique in junction tree consists entirely of observed values. Invalid graph.");
    warning("Separator clique has %d nodes:",nodes.size()); 
    printRVSet(stderr,nodes);
    error("...EXITING...");
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
						   vector <RandomVariable*>& newRvs,
						   map < RVInfo::rvParent, unsigned >& ppf,
						   const unsigned int frameDelta)
  : origin(from_clique)
{

  set<RandomVariable*>::iterator it;

  // clone over nodes RVs.
  fNodes.resize(origin.nodes.size());
  unsigned i=0;
  for (it = origin.nodes.begin();
       it != origin.nodes.end();
       it++) {
    RandomVariable* rv = (*it);
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    


    if ( ppf.find(rvp) == ppf.end() ) {
      coredump("ERROR: can't find rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
    }

    RandomVariable* nrv = newRvs[ppf[rvp]];
    fNodes[i++] = nrv;
  }

  i=0;
  fAccumulatedIntersection.resize(origin.accumulatedIntersection.size());
  for (it = origin.accumulatedIntersection.begin();
       it != origin.accumulatedIntersection.end();
       it++) {
    RandomVariable* rv = (*it);
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    


    if ( ppf.find(rvp) == ppf.end() ) {
      coredump("ERROR: can't find assigned rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
    }

    RandomVariable* nrv = newRvs[ppf[rvp]];
    fAccumulatedIntersection[i] = nrv;
  }

  i=0;
  fRemainder.resize(origin.remainder.size());
  for (it = origin.remainder.begin();
       it != origin.remainder.end();
       it++) {
    RandomVariable* rv = (*it);
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    


    if ( ppf.find(rvp) == ppf.end() ) {
      coredump("ERROR: can't find assigned rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
    }
    RandomVariable* nrv = newRvs[ppf[rvp]];
    fRemainder[i++] = nrv;
  }

  // Separator accumulated intersection values only store/hash values
  // of hidden (thus necessarily discrete) variables since they are
  // the only thing that change.
  accDiscreteValuePtrs.resize(origin.hAccumulatedIntersection.size());
  for (i=0;i<accDiscreteValuePtrs.size();i++) {
    // get the hidden rv for this location
    RandomVariable* rv = origin.hAccumulatedIntersection[i];;
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    

    if ( ppf.find(rvp) == ppf.end() ) {
      coredump("ERROR: can't find assigned rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
    }
    RandomVariable* nrv = newRvs[ppf[rvp]];

    // hidden nodes are always discrete (in this version).
    DiscreteRandomVariable* drv = 
      (DiscreteRandomVariable*)nrv;

    // grab a pointer directly to its value for easy access later.
    accDiscreteValuePtrs[i] = &(drv->val);
  }

  // Separator remainder values only store/hash values of hidden (thus
  // necessarily discrete) variables since they are the only thing
  // that change.
  remDiscreteValuePtrs.resize(origin.hRemainder.size());
  for (i=0;i<remDiscreteValuePtrs.size();i++) {
    // get the hidden rv for this location
    RandomVariable* rv = origin.hRemainder[i];
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    

    if ( ppf.find(rvp) == ppf.end() ) {
      coredump("ERROR: can't find assigned rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
    }
    RandomVariable* nrv = newRvs[ppf[rvp]];

    // hidden nodes are always discrete (in this version).
    DiscreteRandomVariable* drv = 
      (DiscreteRandomVariable*)nrv;

    // grab a pointer directly to its value for easy access later.
    remDiscreteValuePtrs[i] = &(drv->val);
  }

  // allocate at one value for now.
  if (origin.hAccumulatedIntersection.size() == 0) {
    // in this case, we'll only need one and never more.
    separatorValues.resize(1);
    // there will always be one used value here.
    numSeparatorValuesUsed = 1;
    new (&separatorValues.ptr[0].iRemHashMap)VHashMapUnsignedUnsignedKeyUpdatable
      (origin.remPacker.packedLen(),2);
  } else {
    // start with something a bit larger
    // TODO: optimize this.
    const unsigned starting_size = 3; // 2000;
    separatorValues.resize(starting_size);
    if (origin.hRemainder.size() > 0) {
      for (unsigned i=0;i<starting_size;i++) {
	// need to re-construct individual hash tables.
	new (&separatorValues.ptr[i].iRemHashMap)VHashMapUnsignedUnsignedKeyUpdatable
	  (origin.remPacker.packedLen(),2);
	// TODO: potentially preallocate default size of  separatorValues.ptr[i].remValues.resize(default);
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

  // compute max
  logpr maxCEsepValue;
  unsigned originalTotalStateSpace = 0;
  for (unsigned asv=0;asv<numSeparatorValuesUsed;asv++) {
    originalTotalStateSpace += separatorValues.ptr[asv].numRemValuesUsed;
    for (unsigned rsv=0;rsv<separatorValues.ptr[asv].numRemValuesUsed;rsv++) {
      if (separatorValues.ptr[asv].remValues.ptr[rsv].p > maxCEsepValue)
	maxCEsepValue = separatorValues.ptr[asv].remValues.ptr[rsv].p;
    }
  }

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
  : cliqueValueSize(_cliqueValueSize),growthFactor(_growthFactor),
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
	

