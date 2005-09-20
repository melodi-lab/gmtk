/*-
 * GMTK_LatticeNodeCPT.cc
 *
 *  Written by Gang Ji <gang@ee.washington.edu>
 * 
 *  $Header$
 * 
 * Copyright (c) 2001, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about
 * the suitability of this software for any purpose.  It is provided
 * "as is" without express or implied warranty.
 *
 */


#include "GMTK_LatticeNodeCPT.h"
#include "GMTK_CPT.h"
#include "GMTK_DiscRV.h"


/*-
 *-----------------------------------------------------------------------
 * LatticeNodeCPT::LatticeNodeCPT()
 *     default constructor
 *
 * Results:
 *     no results
 *-----------------------------------------------------------------------
 */
LatticeNodeCPT::LatticeNodeCPT() : CPT(di_LatticeNodeCPT) , _latticeAdt(NULL) {
  // some values are fixed
  // the first parent must be previous node and second one must be
  // word transition
  _numParents = 2;
  cardinalities.resize(2);
}


/*-
 *-----------------------------------------------------------------------
 * LatticeNodeCPT::~LatticeNodeCPT()
 *     default destructor
 *
 * Results:
 *     no results
 *-----------------------------------------------------------------------
 */
LatticeNodeCPT::~LatticeNodeCPT() {
}


/*-
 *-----------------------------------------------------------------------
 * LatticeNodeCPT::becomeAwareOfParentValuesAndIterBegin(vector< RV* >& parents, iterator &it, DiscRV* drv, logpr& p)
 *     assign parent value and begin an iterator
 *
 * Results:
 *     no results
 * 
 * Notes:
 *     need more implementation on different scores
 *-----------------------------------------------------------------------
 */
void LatticeNodeCPT::becomeAwareOfParentValuesAndIterBegin(vector< RV* >& parents, iterator &it, DiscRV* drv, logpr& p) {
  // For simplity reason, here cardinality of lattice nodes can be
  // bigger than number of real nodes in some particular lattice.
  // This is because in iterable lattices, different lattices can
  // have different number of nodes.  But in master file, we can
  // just specify the max of those.
  if ( RV2DRV(parents[0])->val > _latticeAdt->_end ) {
    it.internalStatePtr = NULL;
    p.set_to_zero();
    return;
  }

  // case on the current frame index
  if ( drv->frame() < _latticeAdt->_latticeNodes[RV2DRV(parents[0])->val].startFrame ) {
    if ( RV2DRV(parents[1])->val ) {
      // word transition not allowed
      it.internalStatePtr = NULL;
      p.set_to_zero();
    } else {
      // no word transition, just copy values from previous frame
      it.internalStatePtr = NULL;
      drv->val = RV2DRV(parents[0])->val;
      p.set_to_one();
    }
  } else if ( drv->frame() > _latticeAdt->_latticeNodes[RV2DRV(parents[0])->val].endFrame ) {
    // we force this cannot happen by setting prob to zero
    it.internalStatePtr = NULL;
    p.set_to_zero();
  } else {
    if ( RV2DRV(parents[1])->val ) {
      // iterate next lattice nodes
      // find the out going edge
      LatticeADT::LatticeNode &node = _latticeAdt->_latticeNodes[RV2DRV(parents[0])->val];

      // if the lattice node is the end, return prob zero
      if ( node.edges.totalNumberEntries() == 0 ) {
	it.internalStatePtr = NULL;
	p.set_to_zero();
	return;
      }

      // check out the out-going edges based on parent value
      shash_map_iter<unsigned, LatticeADT::LatticeEdge>::iterator *pit = new shash_map_iter<unsigned, LatticeADT::LatticeEdge>::iterator();

      // now the current node can have next transition
      // find the correct iterators
      node.edges.begin(*pit);

      // set up the internal state for iterator
      it.internalStatePtr = (void*)pit;
      it.internalState = RV2DRV(parents[0])->val;
      it.drv = drv;

      drv->val = pit->key();
      p = (**pit).prob_score;
    } else {
      // no word transition, copy value
      it.internalStatePtr = NULL;
      drv->val = RV2DRV(parents[0])->val;
      p.set_to_one();
    }
  }
}


/*-
 *-----------------------------------------------------------------------
 * LatticeNodeCPT::probGivenParents(vector< RV* >& parents, DiscRV* drv)
 *     calculate the probability with known parent/child values
 *
 * Results:
 *     log probability of child given the parents
 * 
 * Notes:
 *     need more implementation on different scores
 *-----------------------------------------------------------------------
 */
logpr LatticeNodeCPT::probGivenParents(vector< RV* >& parents, DiscRV* drv) {
	LatticeADT::LatticeEdge* outEdge = _latticeAdt->_latticeNodes[RV2DRV(parents[0])->val].edges.find(drv->val);

	if ( outEdge == NULL )
		return logpr(0.0);
	else
		return outEdge->prob_score;	// TODO: implement other scores
}


/*-
 *-----------------------------------------------------------------------
 * LatticeNodeCPT::next(iterator &it, logpr& p)
 *     proceed to next iterator
 *
 * Results:
 *     true if the next exists
 * 
 * Notes:
 *     need more implementation on different scores
 *-----------------------------------------------------------------------
 */
bool LatticeNodeCPT::next(iterator &it, logpr& p) {
  shash_map_iter<unsigned, LatticeADT::LatticeEdge>::iterator* pit = (shash_map_iter<unsigned, LatticeADT::LatticeEdge>::iterator*) it.internalStatePtr;

  // check whether pit is null
  if ( pit == NULL ) {
    p.set_to_zero();
    return false;
  }

  // find the next available in the tree
  if ( pit->next() ) {
    // set up the values
    it.drv->val = pit->key();
    p = (**pit).prob_score;		// TODO: implement other scores

    return true;
  } else {
    // we didn't find anything satisfying the frame constrain
    delete pit;
    it.internalStatePtr = NULL;
    p.set_to_zero();
    return false;
  }
}


/*-
 *-----------------------------------------------------------------------
 * LatticeNodeCPT::setLatticeADT(const LatticeADT &latticeAdt)
 *     set the lattice ADT to be used for this CPT
 *
 * Results:
 *     none.
 *-----------------------------------------------------------------------
 */
void LatticeNodeCPT::setLatticeADT(const LatticeADT &latticeAdt) {
  _latticeAdt = &latticeAdt;
  // first parent is noade
  _card = cardinalities[0] = _latticeAdt->_nodeCardinality;
  // second parent is word transition
  cardinalities[1] = 2;
}
