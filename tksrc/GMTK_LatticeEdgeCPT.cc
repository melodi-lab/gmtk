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
	  assert ( RV2DRV(parents[0])->val == RV2DRV(parents[1])->val );
	  // If this occurs, it means that there is no outgoing edge
	  // in the lattice starting from the node with value parent0->val.
	  // We are guaranteed that the corresponding LatticeNodeCPT, when
	  // this occurs, will do two things:
	  //    1) it will give this event zero probability
	  //    2) the value given to the child in LatticeNodeCPT will
	  //       be the same node value as the parent, meaning that
	  //       parent0->val == parent1->val.
	  // Note that if the graph is evaluated topologically, then
	  // the LatticeNodeCPT will be evaluated first (before this)
	  // and inference should prune away before we ever get to this
	  // case. For some triangluations, however, the LatticeEdgeCPT might
	  // be evaluated first, so we need to cover this case here.

	  // ultimately, since the LatticeNodeCPT gives this zero probabilty,
	  // we know this case will get trimmed away, but for now we just
	  // give some junk values.

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
 * Function LatticeEdgeCPT::assignDeterministicChild
 *
 *      Assign the value of the child random variable to be the
 *      determinstic function of the current parent assignments. We
 *      can do this since this CPT is "determinstic", meaning that
 *      given a set of parent assignments, there is only one child
 *      assignment that exists that has non zero (and unity)
 *      probability.
 *
 * Results:
 *      none
 *-----------------------------------------------------------------------
 */
void LatticeEdgeCPT::assignDeterministicChild( vector < RV* >& parents, DiscRV* drv )
{
  LatticeADT::LatticeEdge* outEdge = _latticeAdt->_latticeNodes[RV2DRV(parents[0])->val].edges.find(RV2DRV(parents[1])->val);
  if ( outEdge == NULL ) {
    // For documentation on this case, see the routine:
    // LatticeEdgeCPT::becomeAwareOfParentValuesAndIterBegin()
    drv->val = 0; // some junk number
  } else
    drv->val = outEdge->emissionId;
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

