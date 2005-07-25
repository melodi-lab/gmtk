/*-
 * GMTK_LatticeEdgeCPT.cc
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


#include "GMTK_LatticeEdgeCPT.h"
#include "GMTK_CPT.h"
#include "GMTK_DiscRV.h"

/*-
 *-----------------------------------------------------------------------
 * Function LatticeEdgeCPT::LatticeEdgeCPT
 *      Default constructor.
 *
 * Results:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
LatticeEdgeCPT::LatticeEdgeCPT() : CPT(di_LatticeEdgeCPT), _latticeAdt(NULL) {
	// some values are fixed
	_numParents = 2;
	cardinalities.resize(2);
}


/*-
 *-----------------------------------------------------------------------
 * Function LatticeEdgeCPT::~LatticeEdgeCPT
 *      Default destructor.
 *
 * Results:
 *      None.
 *-----------------------------------------------------------------------
 */
LatticeEdgeCPT::~LatticeEdgeCPT() {
}


/*-
 *-----------------------------------------------------------------------
 * Function LatticeEdgeCPT::becomeAwareOfParentValuesAndIterBegin
 *      Begin an iterator with the parents values known.
 *
 * Results:
 *      None.
 *-----------------------------------------------------------------------
 */
void LatticeEdgeCPT::becomeAwareOfParentValuesAndIterBegin(vector< RV* >& parents, iterator &it, DiscRV* drv, logpr& p) {
	LatticeADT::LatticeEdge* outEdge = _latticeAdt->_latticeNodes[RV2DRV(parents[0])->val].edges.find(RV2DRV(parents[1])->val);
	if ( outEdge == NULL ) {
		// this should report an error because I have set prob in
		// latticeNodeCPT to zero.
		//error("Error: lattice edge CPT is not synchronized with node cpt");
		// I changed it to be zero prob.
		it.drv = drv;
		drv->val = 0;	// some junk number
		p.set_to_zero();
		return;
	}

	it.drv = drv;
	drv->val = outEdge->emissionId;
	p.set_to_one();
}


/*-
 *-----------------------------------------------------------------------
 * Function LatticeEdgeCPT::probGivenParents
 *      Get probability with parents and child value known.
 *
 * Results:
 *      Probability.
 *
 *-----------------------------------------------------------------------
 */
logpr LatticeEdgeCPT::probGivenParents(vector< RV* >& parents, DiscRV* drv) {
	LatticeADT::LatticeEdge* outEdge = _latticeAdt->_latticeNodes[RV2DRV(parents[0])->val].edges.find(RV2DRV(parents[1])->val);
	if ( outEdge == NULL || outEdge->emissionId != drv->val ) {
		return logpr(0.0);
	} else {
		return logpr(1.0);
	}
}


/*-
 *-----------------------------------------------------------------------
 * Function LatticeEdgeCPT::next
 *      Advance an iterator.
 *
 * Results:
 *      Always return false since this is a deterministic mapping.
 *-----------------------------------------------------------------------
 */
bool LatticeEdgeCPT::next(iterator &it, logpr& p) {
	p.set_to_zero();
	return false;
}


/*-
 *-----------------------------------------------------------------------
 * Function LatticeEdgeCPT::setLatticeADT
 *      Set the lattice adt.
 *
 * Results:
 *      None.
 *-----------------------------------------------------------------------
 */
void LatticeEdgeCPT::setLatticeADT(const LatticeADT &latticeAdt) {
        _latticeAdt = &latticeAdt;
        cardinalities[0] = cardinalities[1] = _latticeAdt->_nodeCardinality;
	_card = _latticeAdt->_wordCardinality;
}

