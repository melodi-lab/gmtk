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
  // parents are current lattice node and previous lattice node.
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
 * Function LatticeEdgeCPT::probGivenParents
 *   Get probability with parents and child value known.
 *   Here, parents[0] is the previous lattice node, parents[1] is the current
 *   lattice node, and we map from a pair of lattice nodes to a lattice edge list
 *   which is one or more lattice edges.
 *
 * Results:
 *      Probability.
 *
 *-----------------------------------------------------------------------
 */
logpr LatticeEdgeCPT::probGivenParents(vector< RV* >& parents, DiscRV* drv) {
  LatticeADT::LatticeEdgeList* outEdges
    = _latticeAdt->_latticeNodes[RV2DRV(parents[0])->val].edges.find(RV2DRV(parents[1])->val);
  if ( outEdges == NULL) {
    // then this is an impossible pair of nodes, so give it zero
    // probability.
    return logpr(0.0);
  }
  // we've got an edge, but how many?  We've got to search to find
  // one with a matching drv->val.
  for (unsigned edge_ctr=0;edge_ctr < outEdges->num_edges; edge_ctr ++ ) {
    LatticeADT::LatticeEdge &edge = outEdges->edge_array[edge_ctr];
    if (edge.emissionId == drv->val) {
      return edge.gmtk_score;
    }
  }

  // still here? return 0, since not found.
  return logpr(0.0);
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
void LatticeEdgeCPT::becomeAwareOfParentValuesAndIterBegin(vector< RV* >& parents, 
							   iterator &it, 
							   DiscRV* drv, 
							   logpr& p) 
{
  LatticeADT::LatticeEdgeList* outEdge 
    = _latticeAdt->_latticeNodes[RV2DRV(parents[0])->val].edges.find(RV2DRV(parents[1])->val);
  if ( outEdge == NULL ) {
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
    // Set to zero (0), a dummy number that has no meaning.
    drv->val = 0;	
    p.set_to_zero();
    return;
  }

  // get first edge of list, we are guaranteed that there is at least
  // one edge in the array.
  assert ( outEdge->edge_array.size() > 0 );
  it.drv = drv;
  it.uInternalState = 0;
  it.internalStatePtr = outEdge;
  drv->val = outEdge->edge_array[it.uInternalState].emissionId;
  p = outEdge->edge_array[it.uInternalState].gmtk_score;
}



/*-
 *-----------------------------------------------------------------------
 * Function LatticeEdgeCPT::next
 *      Advance an iterator.
 *
 * Results:
 *      True if there is a next one.
 *-----------------------------------------------------------------------
 */
bool LatticeEdgeCPT::next(iterator &it, logpr& p)
{
  LatticeADT::LatticeEdgeList* outEdge =
    (LatticeADT::LatticeEdgeList*) it.internalStatePtr;
  if (it.uInternalState+1 < outEdge->edge_array.size()) {
    it.uInternalState++;
    it.drv->val = outEdge->edge_array[it.uInternalState].emissionId;
    p = outEdge->edge_array[it.uInternalState].gmtk_score;
    return true;
  } else {
    p.set_to_zero();
    return false;
  }
}


#if 0
// this is no longer determinisic 

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
    // this shouldn't happen under the current inference since zero
    // score clique eventsa are never inserted into the clique (so
    // we'll never need to reconstruct a clique entry that would have
    // zero probability). If the outEdge is null, this would have caused
    // the clique entry to have zero probability.

    // Just for consistency, include the following code.
    // For documentation on this case, see the routine:
    // LatticeEdgeCPT::becomeAwareOfParentValuesAndIterBegin()
    drv->val = 0; // some junk number
  } else
    drv->val = outEdge->emissionId;
}
#endif


/*-
 *-----------------------------------------------------------------------
 * Function LatticeEdgeCPT::setLatticeADT
 *      Set the lattice adt.
 *
 * Results:
 *      None.
 *-----------------------------------------------------------------------
 */
void LatticeEdgeCPT::setLatticeADT(const LatticeADT &latticeAdt) 
{
  _latticeAdt = &latticeAdt;

  // In addition to lattice ADT, also set the cardinalties of the
  // parents.

  // Typically, cardinalities comes from structure file or master
  // file. But in this case, we support iterable lattice CPTs which
  // will have different number of lattice nodes for each one.
  cardinalities[0] = cardinalities[1] = _latticeAdt->_nodeCardinality;
  _card = _latticeAdt->_wordCardinality;
}

