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
	_numParents = 1;
	cardinalities.resize(1);
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
	// check out the out-going edges based on parent value
	shash_map_iter<unsigned, LatticeADT::LatticeEdge>::iterator *pit = new shash_map_iter<unsigned, LatticeADT::LatticeEdge>::iterator();

	LatticeADT::LatticeNode &node = _latticeAdt->_latticeNodes[RV2DRV(parents[0])->val];

	// if the lattice node is the end, return prob zero
	if ( node.edges.totalNumberEntries() == 0 ) {
		p.set_to_zero();
		return;
	}

	// now the current node can have next transition
	// find the correct iterators
	node.edges.begin(*pit);

	// check the frame number and parent frame index to see whether
	// transition is allowed at this frame
	if ( drv->frame() >= _latticeAdt->_latticeNodes[pit->key()].startFrame
			&& drv->frame() <= _latticeAdt->_latticeNodes[pit->key()].endFrame ) {
		// easy case
		// set up the iternal state for iterators
		it.internalStatePtr = (void*)pit;
		it.drv = drv;

		// set up the values
		drv->val = pit->key();
		p = (**pit).prob_score;		// TODO: implement other scores

		return;
	} else {
		// need to find the next available in the tree
		while ( pit->next() ) {
			if ( drv->frame() >= _latticeAdt->_latticeNodes[pit->key()].startFrame && drv->frame() <= _latticeAdt->_latticeNodes[pit->key()].endFrame ) {
				// this is what we need
				// set up the iternal state for iterators
				it.internalStatePtr = (void*)pit;
				it.drv = drv;

				// set up the values
				drv->val = pit->key();
				p = (**pit).prob_score;		// TODO: implement other scores

				return;
			}
		}

		// we didn't find anything satisfying the frame constrain
		delete pit;
		it.internalStatePtr = NULL;
		p.set_to_zero();
		return;
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
	// althout this shouldn't happen
	// 1. if the end user didn't check prob value in begin
	//    the system will not know that next shouldn't be called
	if ( pit == NULL ) {
		// theoretically, this should throw an error
		// but the way inference engine is designed now, this case
		// might still happen
		p.set_to_zero();
		return false;
		//error("Error: LatticeNodeCPT::next shouldn't be called when begin returns a zero probability");
	}

	// need to find the next available in the tree
	while ( pit->next() ) {
		if ( it.drv->frame() >= _latticeAdt->_latticeNodes[pit->key()].startFrame && it.drv->frame() <= _latticeAdt->_latticeNodes[pit->key()].endFrame ) {
			// set up the values
			it.drv->val = pit->key();
			p = (**pit).prob_score;		// TODO: implement other scores

			return true;
		}
	}

	// we didn't find anything satisfying the frame constrain
	delete pit;
	it.internalStatePtr = NULL;
	p.set_to_zero();
	return false;
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
	_card = cardinalities[0] = _latticeAdt->_nodeCardinality;
}
