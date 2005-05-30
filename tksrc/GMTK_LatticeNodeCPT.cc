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

LatticeNodeCPT::LatticeNodeCPT() : CPT(di_LatticeNodeCPT) , _latticeAdt(NULL) {
	// some values are fixed
	_numParents = 1;
	cardinalities.resize(1);
}


LatticeNodeCPT::~LatticeNodeCPT() {
}


void LatticeNodeCPT::becomeAwareOfParentValuesAndIterBegin(vector< RV* >& parents, iterator &it, DiscRV* drv, logpr& p) {
	// check out the out-going edges based on parent value
	shash_map_iter<unsigned, LatticeADT::LatticeEdge>::iterator *pit = new shash_map_iter<unsigned, LatticeADT::LatticeEdge>::iterator();
	assert(_latticeAdt->_latticeNodes[RV2DRV(parents[0])->val].edges.totalNumberEntries() != 0);
	_latticeAdt->_latticeNodes[RV2DRV(parents[0])->val].edges.begin(*pit);

	// set up the iternal state for iterators
	it.internalStatePtr = (void*)pit;
	it.drv = drv;

	// set up the values
	drv->val = pit->key();
	p = (**pit).ac_score;		// TODO: instead of acoustic score, implement other scores
}


logpr LatticeNodeCPT::probGivenParents(vector< RV* >& parents, DiscRV* drv) {
	LatticeADT::LatticeEdge* outEdge = _latticeAdt->_latticeNodes[RV2DRV(parents[0])->val].edges.find(drv->val);


	if ( outEdge == NULL )
		return logpr(0.0);
	else
		return outEdge->ac_score;	// TODO: implement other scores
}


bool LatticeNodeCPT::next(iterator &it, logpr& p) {
	shash_map_iter<unsigned, LatticeADT::LatticeEdge>::iterator* pit = (shash_map_iter<unsigned, LatticeADT::LatticeEdge>::iterator*) it.internalStatePtr;
	if ( pit->next() ) {
		it.drv->val = pit->key();
		p = (**pit).ac_score;		// TODO: instead of acoustic score, implement other scores
		return true;
	} else {
		delete pit;
		return false;
	}
}


void LatticeNodeCPT::setLatticeADT(const LatticeADT &latticeAdt) {
	_latticeAdt = &latticeAdt;
	_card = cardinalities[0] = _latticeAdt->_nodeCardinality;
}


