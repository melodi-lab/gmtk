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
 *   In some of the structures the choice of when
 *   to use a hash table vs. an inlined value is done
 *   by either using an unsigned* or a unsigned. Create
 *   options where we do 'unsigned val[LEN]' and where
 *   only if packed value is > LEN words do we resort to
 *   the hash table.
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

    // TODO: ultimately turn this just into an assert
    if ( ppf.find(rvp) == ppf.end() ) {
      warning("ERROR: can't find rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
      assert ( ppf.find(rvp) != ppf.end() );
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

    // TODO: ultimately turn this just into an assert
    if ( ppf.find(rvp) == ppf.end() ) {
      warning("ERROR: can't find assigned rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
      assert ( ppf.find(rvp) != ppf.end() );
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

    // TODO: ultimately turn this just into an assert
    if ( ppf.find(rvp) == ppf.end() ) {
      warning("ERROR: can't find assigned rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
      assert ( ppf.find(rvp) != ppf.end() );
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

    // TODO: ultimately turn this just into an assert
    if ( ppf.find(rvp) == ppf.end() ) {
      warning("ERROR: can't find assigned rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
      assert ( ppf.find(rvp) != ppf.end() );
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
 *     none
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
      if (useDeterminism && drv->deterministic()) {
	// then there is a possibility that this node
	// does not affect the state space, as long
	// as all of this nodes parents are in the clique.
	bool truly_deterministic = true;
	for (unsigned i=0;i<drv->allPossibleParents.size();i++) {
	  if (nodes.find(drv->allPossibleParents[i]) == nodes.end()) {
	    // found a parent that is not in 'node' set so the
	    // node would not truly be deterministic here.
	    truly_deterministic = false;
	    break;
	  }
	}
	if (!truly_deterministic)
	  tmp_weight += log10((double)drv->useCardinality());	
      } else
	tmp_weight += log10((double)drv->useCardinality());
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
      if (useDeterminism && drv->deterministic()) {
	// then there is a possibility that this node
	// does not affect the state space, as long
	// as all of this nodes parents are in the clique.
	bool truly_deterministic = true;
	for (unsigned i=0;i<drv->allPossibleParents.size();i++) {
	  if ((nodes.find(drv->allPossibleParents[i]) == nodes.end())
	      &&
	      (drv->allPossibleParents[i] != node)) {
	    // found a parent that is not in 'node' set so the
	    // node would not truly be deterministic here.
	    truly_deterministic = false;
	    break;
	  }
	}
	if (!truly_deterministic)
	  tmp_weight += log10((double)drv->useCardinality());
      } else
	tmp_weight += log10((double)drv->useCardinality());
    }
  }
  return tmp_weight;
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

  // TODO: optimize initial size and growth factor.  Compute an
  // estimate of the state space of this clique for a starting
  // allocation of the value holder.  Take 1/4 of the estimated weight
  // for now.
  allocationUnitChunkSize =
    (unsigned)(::pow(10,weight() - ::log10(4.0)));
  // lower bound.
  if (allocationUnitChunkSize < 16)
    allocationUnitChunkSize = 16;

  if (packer.packedLen() > 1) {
    // setup value hodler
    new (&valueHolder) CliqueValueHolder(hiddenNodes.size(),
					 allocationUnitChunkSize,
					 1.25);
    // set up common clique hash tables
    new (&hashTable) vhash_set< unsigned > (hiddenNodes.size());
  } else {
    // then no need to do a hash table at all, just store the packed
    // values in a local integer.
  }
}



////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        InferenceMaxClique support
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

    // TODO: ultimately turn this just into an assert
    if ( ppf.find(rvp) == ppf.end() ) {
      warning("ERROR: can't find rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
      assert ( ppf.find(rvp) != ppf.end() );
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

    // TODO: ultimately turn this just into an assert
    if ( ppf.find(rvp) == ppf.end() ) {
      warning("ERROR: can't find assigned rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
      assert ( ppf.find(rvp) != ppf.end() );
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

    // TODO: ultimately turn this just into an assert
    if ( ppf.find(rvp) == ppf.end() ) {
      warning("ERROR: can't find assigned rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
      assert ( ppf.find(rvp) != ppf.end() );
    }
    RandomVariable* nrv = newRvs[ppf[rvp]];
    fUnassignedIteratedNodes[i++] = nrv;
  }

  discreteValuePtrs.resize(from_clique.hiddenNodes.size());
  for (i=0;i<discreteValuePtrs.size();i++) {
    // hidden nodes are always discrete (in this version).
    DiscreteRandomVariable* drv = 
      (DiscreteRandomVariable*)from_clique.hiddenNodes[i];
    // grab a pointer directly to its value for easy access later.
    discreteValuePtrs[i] = &(drv->val);
  }

  numCliqueValuesUsed = 0;
  maxCEValue.set_to_zero();

  // TODO: optimize this.
  cliqueValues.resize(1000);

}


/*-
 *-----------------------------------------------------------------------
 * MaxClique::CollectEvidenceFromSeparators()
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
void
InferenceMaxClique::collectEvidenceFromSeparators(JT_InferencePartition& part)
{
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
  unsigned sepValueNumber;
  if (sep.fAccumulatedIntersection.size() > 0) {
    // look up existing intersected values to see if we have a match
    // and only proceed if we do.
    
    // find index if it exists.
    unsigned tmp;
    unsigned *key;
    unsigned *data;
    // TODO: optimize this check out of loop
    if (sep.origin.accPacker.packedLen() == 1) {
      key = &tmp;
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

    data = sep.iAccHashTable.find(key);
    if (data == NULL) {
      // Then not found in this separator, so it must have zero
      // probability. We continue with the next value of the previous
      // separator.
      return;
    } else {
      // need to further iterate.
      sepValueNumber = *data;
    }

  } else {
    // TODO: check this, as this condition might fail if
    // we've completely pruned away parent clique.
    assert ( sep.separatorValues.size() == 1);
    sepValueNumber = 0;
  }

  // iterate through remainder of separator
  for (unsigned i=0;i< sep.separatorValues[sepValueNumber].numRemValuesUsed; i++) {
    
    // TODO: optimize away this check here by duplicating code.
    if (sep.origin.remPacker.packedLen() == 1) {
      sep.origin.remPacker.unpack(
				  (unsigned*)&(sep.separatorValues[sepValueNumber].remValues[i].val),
				  (unsigned**)sep.remDiscreteValuePtrs.ptr);
    } else {
      sep.origin.remPacker.unpack(
				  (unsigned*)sep.separatorValues[sepValueNumber].remValues[i].ptr,
				  (unsigned**)sep.remDiscreteValuePtrs.ptr);
    }
    // continue down with new probability value.
    // Q: should we do some pruning here as well?
    ceIterateSeparators(part,sepNumber+1,
			p*
			sep.separatorValues[sepValueNumber].remValues[i].p);
  }

}



void
InferenceMaxClique::ceIterateAssignedNodes(JT_InferencePartition& part,
					   const unsigned nodeNumber,
					   const logpr p)
{
  if (nodeNumber == fSortedAssignedNodes.size()) {
    // time to store clique value and total probability, p is
    // current clique probability.

    // keep track of the max clique probability right here.
    if (p > maxCEValue)
      maxCEValue = p;

    if (numCliqueValuesUsed >= cliqueValues.size()) {
      // TODO: optimize this.
      cliqueValues.resizeAndCopy(cliqueValues.size()*2);
    }

    // TODO: figure out if it is possible to get around doint this
    // check (or to help branch prediction predict it, since it will
    // be different for differnet cliques). Answer: We can do this
    // check once at beginning of iteration of assigned nodes, and
    // have two versions of this code.
    if (origin.packer.packedLen() == 1) {
      // pack the clique values directly into place
      origin.packer.pack(
			 (unsigned**)discreteValuePtrs.ptr,
			 (unsigned*)&(cliqueValues[numCliqueValuesUsed].val));
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
      key = origin.hashTable.insert(pcv,foundp);
      if (!foundp) {
	// if it was not found, need to claim this storage that we
	// just used.
	origin.valueHolder.allocateCurCliqueValue();
      }
      // Save the pointer to whatever the hash table decided to use.
      cliqueValues[numCliqueValuesUsed].ptr = key;
    }
    // save the probability
    cliqueValues[numCliqueValuesUsed].p = p;
    numCliqueValuesUsed++;

    return;
  }
  RandomVariable* rv = fSortedAssignedNodes[nodeNumber];
  // do the loop right here
  rv->clampFirstValue();
  do {
    // At each step, we compute probability
    logpr cur_p = rv->probGivenParents();
    // if at any step, we get zero, then back out.
    if (!p.essentially_zero()) {
      // Q: should be do more severe pruning here as well as beam?
      // 
      // Continue, updating probability by cur_p.
      ceIterateAssignedNodes(part,nodeNumber+1,p*cur_p);
    }
  } while (rv->clampNextValue());
}


void
InferenceMaxClique::ceIterateUnassignedIteratedNodes(JT_InferencePartition& part,
						     const unsigned nodeNumber,
						     const logpr p)
{
  if (nodeNumber == fUnassignedIteratedNodes.size()) {
    ceIterateAssignedNodes(part,0,p);
  }
  RandomVariable* rv = fUnassignedIteratedNodes[nodeNumber];
  if (rv->hidden) {
    // only discrete RVs can be hidden for now.
    DiscreteRandomVariable* drv = (DiscreteRandomVariable*)rv;
    // do the loop right here
    drv->val = 0;
    do {
      ceIterateUnassignedIteratedNodes(part,nodeNumber+1,p);
    } while (++drv->val < drv->cardinality);
  } else {
    // observed, either discrete or continuous
    if (rv->discrete) {
      DiscreteRandomVariable* drv = (DiscreteRandomVariable*)rv;
      // TODO: for observed variables, do this once at the begining
      // before any looping here.
      drv->setToObservedValue();
    } else {
      // nothing to do since we get continuous observed
      // value indirectly
    }
    ceIterateUnassignedIteratedNodes(part,nodeNumber+1,p);
  }
}

/*
 * we have now a fully instantiated clique. Iterate
 * through the values that are above beam and instantiate
 * the outgoing separator with those values.
 *
 */

void 
InferenceMaxClique::
ceCollectToSeparator(JT_InferencePartition& part)
{
  ceCollectToSeparator(part,
		       part.separatorCliques[origin.ceSendSeparator]);
}
void 
InferenceMaxClique::
ceCollectToSeparator(JT_InferencePartition& part,
		     InferenceSeparatorClique& sep)
{
  
  for (unsigned i=0;i<numCliqueValuesUsed;i++) {{
    // TODO: beam pruning
    //if (cliqueValues[i].p < beamThreshold) {
    //     continue;
    //}

    // TODO: optimiae away this conditional check.
    if (origin.packer.packedLen() == 1) {
      origin.packer.unpack((unsigned*)&(cliqueValues[i].val),
			   (unsigned**)discreteValuePtrs.ptr);
    } else {
      origin.packer.unpack((unsigned*)cliqueValues[i].ptr,
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
    if (sep.origin.accumulatedIntersection.size() > 0) {
      // an accumulated intersection exists.

      // make sure there is at least one available entry
      if (sep.numSeparatorValuesUsed >= sep.separatorValues.size()) {
	const unsigned old_size = sep.separatorValues.size();
	sep.separatorValues.resizeAndCopy(sep.separatorValues.size()*2);
	const unsigned new_size = sep.separatorValues.size();
	// need to re-construct hash tables for new entries.
	for (unsigned i=old_size;i<new_size;i++) {
	  new (&sep.separatorValues[i].iRemHashTable)vhash_map< unsigned, unsigned >
	    (sep.origin.remPacker.packedLen());
	}
      }
      
      unsigned *accKey;
      // TODO: optimize this check out of loop.
      if (sep.origin.accPacker.packedLen() == 1) {
	accKey = &(sep.separatorValues[sep.numSeparatorValuesUsed].val);
	sep.origin.accPacker.pack((unsigned**)sep.accDiscreteValuePtrs.ptr,
				  accKey);
      } else {
	accKey = sep.origin.accValueHolder.curCliqueValuePtr();
	sep.origin.accPacker.pack((unsigned**)sep.accDiscreteValuePtrs.ptr,
				  accKey);
	// check if this value combination already lives in
	// origin's value holder hash table and if so, use that.
	bool foundp;
	accKey = sep.origin.accHashTable.insert(accKey,foundp);
	if (!foundp) {
	  // only allocate a new value if it was inserted.
	  sep.origin.accValueHolder.allocateCurCliqueValue();
	}
	// store the pointer in case we use it.
	sep.separatorValues[sep.numSeparatorValuesUsed].ptr = accKey;
      }
      
      bool foundp;
      unsigned* accIndexp =
	sep.iAccHashTable.insert(accKey,
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

	// go ahead and insert it here to the 1st entry.

	// keep handy reference for readability.
	InferenceSeparatorClique::AISeparatorValue& sv
	  = sep.separatorValues[accIndex];

	// Accumulate the clique's
	// probability into this separator's probability.
	if (sv.remValues.size() < 1) {
	  // this must be first time for this entry.
	  sv.remValues.resize(1);
	  sv.numRemValuesUsed = 1;	  
	  sv.remValues[0].p = cliqueValues[i].p;
	} else {
	  sv.remValues[0].p += cliqueValues[i].p;
	}

	goto next_iteration;
      }

    } else {
      accIndex = 0;
    }
    // Do the remainder exists in this separator.
    // 
    // either:
    //   1) AI exists and REM exist
    //     or
    //   3) AI does not exist, but REM exists
    // 

    // keep handy reference for readability.
    InferenceSeparatorClique::AISeparatorValue& sv
      = sep.separatorValues[accIndex];
    
    // make sure there is at least one available entry
    if (sv.numRemValuesUsed >= sv.remValues.size()) {
      sv.remValues.resizeAndCopy(sv.remValues.size()*2);
    }


    unsigned *remKey;
    // pack relevant variable values
    // TODO: optimize away this check.
    if (sep.origin.remPacker.packedLen() == 1) {
      // grab pointer to next location to be used in this case.
      remKey = &(sv.remValues[sv.numRemValuesUsed].val);
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
      remKey = sep.origin.remHashTable.insert(remKey,foundp);
      if (!foundp) {
	// only allocate a new value if it was inserted.
	sep.origin.remValueHolder.allocateCurCliqueValue();
      }
      // store the pointer in case we use it.
      sv.remValues[sv.numRemValuesUsed].ptr = remKey;
    }

    bool foundp;
    unsigned* remIndexp =
      sv.iRemHashTable.insert(
			      remKey,
			      sv.numRemValuesUsed,
			      foundp);
    if (!foundp) {
      // add the values we just used. 
      sv.numRemValuesUsed++;
    }

    // We've finally got the entry, so accumulate the clique's
    // probability into this separator's probability.
    sv.remValues[*remIndexp].p += cliqueValues[i].p;

  }
  next_iteration:
  ;
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

    // TODO: ultimately turn this just into an assert
    if ( ppf.find(rvp) == ppf.end() ) {
      warning("INTERNAL ERROR: can't find rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
      assert ( ppf.find(rvp) != ppf.end() );
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

    // TODO: ultimately turn this just into an assert
    if ( ppf.find(rvp) == ppf.end() ) {
      warning("INTERNAL ERROR: can't find rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
      assert ( ppf.find(rvp) != ppf.end() );
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

    // TODO: ultimately turn this just into an assert
    if ( ppf.find(rvp) == ppf.end() ) {
      warning("INTERNAL ERROR: can't find rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
      assert ( ppf.find(rvp) != ppf.end() );
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
  set<RandomVariable*>::iterator it;
  for (it = accumulatedIntersection.begin();
       it != accumulatedIntersection.end();
       it++) {
    RandomVariable* rv = (*it);
    if (rv->hidden)
      hAccumulatedIntersection.push_back(rv);
  }

  new (&accPacker) PackCliqueValue(hAccumulatedIntersection);

  if (accPacker.packedLen() > 1) {
    // only setup hash table if the packed accumulated insersection
    // set is larger than one machine word (unsigned).
    new (&accValueHolder) CliqueValueHolder(hAccumulatedIntersection.size(),
					    // TODO: optimize this 1000 value.
					    1000,
					    1.25);
    new (&accHashTable) vhash_set< unsigned > (hAccumulatedIntersection.size());
  }

  for (it = remainder.begin();
       it != remainder.end();
       it++) {
    RandomVariable* rv = (*it);
    if (rv->hidden)
      hRemainder.push_back(rv);
  }

  new (&remPacker) PackCliqueValue(hRemainder);

  if (remPacker.packedLen() > 1) { 
    // Only setup hash table if the packed remainder set is larger
    // than one machine word (unsigned).
    new (&remValueHolder) CliqueValueHolder(hRemainder.size(),
					    // TODO: optimize this
					    1000,
					    1.25);
    new (&remHashTable) vhash_set< unsigned > (hRemainder.size());
  }

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

    // TODO: ultimately turn this just into an assert
    if ( ppf.find(rvp) == ppf.end() ) {
      warning("ERROR: can't find rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
      assert ( ppf.find(rvp) != ppf.end() );
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

    // TODO: ultimately turn this just into an assert
    if ( ppf.find(rvp) == ppf.end() ) {
      warning("ERROR: can't find assigned rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
      assert ( ppf.find(rvp) != ppf.end() );
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

    // TODO: ultimately turn this just into an assert
    if ( ppf.find(rvp) == ppf.end() ) {
      warning("ERROR: can't find assigned rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	    rv->name().c_str(),rv->frame(),frameDelta,
	    rvp.first.c_str(),rvp.second);
      assert ( ppf.find(rvp) != ppf.end() );
    }
    RandomVariable* nrv = newRvs[ppf[rvp]];
    fRemainder[i++] = nrv;
  }

  accDiscreteValuePtrs.resize(origin.hAccumulatedIntersection.size());
  for (i=0;i<accDiscreteValuePtrs.size();i++) {
    // hidden nodes are always discrete (in this version).
    DiscreteRandomVariable* drv = 
      (DiscreteRandomVariable*)origin.hAccumulatedIntersection[i];
    // grab a pointer directly to its value for easy access later.
    accDiscreteValuePtrs[i] = &(drv->val);
  }

  remDiscreteValuePtrs.resize(origin.hRemainder.size());
  for (i=0;i<remDiscreteValuePtrs.size();i++) {
    // hidden nodes are always discrete (in this version).
    DiscreteRandomVariable* drv = 
      (DiscreteRandomVariable*)origin.hRemainder[i];
    // grab a pointer directly to its value for easy access later.
    remDiscreteValuePtrs[i] = &(drv->val);
  }

  // allocate at one value for now.
  if (origin.accumulatedIntersection.size() == 0) {
    // in this case, we'll only need one and never more.
    separatorValues.resize(1);
    new (&separatorValues[0].iRemHashTable)vhash_map< unsigned, unsigned >
      (origin.remPacker.packedLen());
  } else {
    // start with something a bit larger
    // TODO: optimize this.
    const unsigned starting_size = 1000;
    separatorValues.resize(starting_size);
    // need to re-construct individual hash tables.
    for (unsigned i=0;i<starting_size;i++) {
      new (&separatorValues[i].iRemHashTable)vhash_map< unsigned, unsigned >
	(origin.remPacker.packedLen());
    }
    // need to re-construct the hash table.
    new (&iAccHashTable) vhash_map< unsigned, unsigned >
      (origin.accPacker.packedLen());
  }
  numSeparatorValuesUsed = 0;



}


/*-
 *-----------------------------------------------------------------------
 * MaxClique::CollectEvidenceToSeparator()
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
  values.resize(1);

  unsigned newSize = cliqueValueSize*allocationUnitChunkSize;
  values[values.size()-1].resize(newSize);

  curAllocationPosition = values[values.size()-1].ptr;
  curAllocationEnd = values[values.size()-1].ptr + newSize;
}


void
CliqueValueHolder::clear()
{
  for (unsigned i=0;i<values.size();i++) {
    values[i].clear();
  }
  values.clear();
}

void
CliqueValueHolder::allocateCurCliqueValue()
{

  curAllocationPosition += cliqueValueSize;

  // first to a cheap and fast allocation of new clique value storage
  if (curAllocationPosition != curAllocationEnd) {
    return;
  } 

  // if here, we need to allocate another chunk add a new chunk so we
  // don't need to re-copy all the existing ones already.
  values.resizeAndCopy(values.size()+1);

  unsigned newSize = unsigned(cliqueValueSize*allocationUnitChunkSize*
			      ::pow(growthFactor,values.size()-1));

  values[values.size()-1].resize(newSize);

  curAllocationPosition = values[values.size()-1].ptr;
  curAllocationEnd = values[values.size()-1].ptr + newSize;

}
