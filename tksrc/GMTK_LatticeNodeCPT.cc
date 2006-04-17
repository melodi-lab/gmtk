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
  // 
  // We assume that there are either 2 or 3 parents,
  // in both cases, 
  //    drv is the current lattice node.
  // In the case of 2 parents, (where time is normally obtained from 'drv', the current lattice node)
  //      parent[0] is the previous lattice node
  //      parent[1] is the  "word transition" (or variable that is acting like such a construct)
  // In the case of 3 parents, 
  //      parent[0] is the previous lattice node
  //      parent[1] is the the "word transition" (or variable that is acting like such a construct)
  //      parent[2] is the "time observation", namely it is a variable that is presumably observed that
  //                keeps track of the time frame to use (rather than using the time frame of 'drv').
  // 
  // For simplity reason, here cardinality of lattice nodes can be
  // bigger than number of real nodes in some particular lattice.
  // This is because in iterable lattices, different lattices can
  // have different number of nodes.  But in master file, we can
  // just specify the max of those.

  // initialize it to something always at least valid.
  drv->val = 0;

  if ( RV2DRV(parents[0])->val > _latticeAdt->_end ) {
    it.internalStatePtr = NULL;
    p.set_to_zero();
    return;
  }

  if ( _latticeAdt->useTimeParent() )
  {
    // use time parent to check time

    // first, compute lat_time, the time of the previous lattice node
    // rounded to the closest frame. parent[0] contains the value of
    // the previous lattice node, and we need to do a lookup in the
    // lattce to find the actual previous lattice node to get its
    // time.
    unsigned lat_time = (unsigned)round(_latticeAdt->_frameRate * _latticeAdt->_latticeNodes[RV2DRV(parents[0])->val].time);

    // case on the current time
    if ( RV2DRV(parents[2])->val < lat_time ) {
      // then the current time is less than the previous lattice node.
      if ( RV2DRV(parents[1])->val ) {
	// a (word) transition is being hypothesized, but we do not allow it. Give it a
	// zero probability.
	it.internalStatePtr = NULL;
	p.set_to_zero();
      } else {
	// a (word) transition is not being hypothesized, so
	// we copy the previous lattice node's value (from previous frame) to drv.
	it.internalStatePtr = NULL;
	drv->val = RV2DRV(parents[0])->val;
	p.set_to_one();
      }
    } else if ( RV2DRV(parents[2])->val > lat_time ) {
      // Then, the current time (i.e., parent[2]'s time value) is
      // already ahead (i.e., after, later) of when a transition for
      // the previous lattice node value may occur. We also want to
      // force this not to happen by setting prob to zero. While
      // you might think this might be valid, only allow jumping
      // from the prevous lattice node when the time is exactly
      // equal to the previous lattice nodes.

      it.internalStatePtr = NULL;
      p.set_to_zero();
    } else {
      // (RV2DRV(parents[2])->val == lat_time), which means that the
      // current time is right at the point that we allow the previous
      // lattice node to jump to the next set of possible lattice
      // nodes.


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
	p = (**pit).gmtk_score;
      } else {
	// no word transition, copy value
	it.internalStatePtr = NULL;
	drv->val = RV2DRV(parents[0])->val;
	p.set_to_one();
      }
    }
  } else {
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
        p = (**pit).gmtk_score;
      } else {
        // no word transition, copy value
        it.internalStatePtr = NULL;
        drv->val = RV2DRV(parents[0])->val;
        p.set_to_one();
      }
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
		return outEdge->gmtk_score;
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
    p = (**pit).gmtk_score;

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

  // check whether time is also used as parent
  if ( _latticeAdt->useTimeParent() ) {
    _numParents = 3;
    cardinalities.resize(3);
    cardinalities[2] = _latticeAdt->_timeCardinality;
  } else {
    _numParents = 2;
    cardinalities.resize(2);
  }

  // first parent is noade
  _card = cardinalities[0] = _latticeAdt->_nodeCardinality;
  // second parent is word transition
  cardinalities[1] = 2;
}
