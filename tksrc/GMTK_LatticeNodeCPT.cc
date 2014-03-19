/*-
 * GMTK_LatticeNodeCPT.cc
 *
 *  Written by Gang Ji <gang@ee.washington.edu> and Jeff Bilmes <bilmes@ee.washington.edu>
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
 * LatticeNodeCPT::probGivenParents(vector< RV* >& parents, DiscRV* drv)
 *     calculate the probability with known parent/child values
 *     parents[0] has previous node, and drv has current node value.
 *
 * Results:
 *     log probability of child given the parents
 * 
 * Notes:
 *     need more implementation on different scores
 *-----------------------------------------------------------------------
 */
logpr LatticeNodeCPT::probGivenParents(vector< RV* >& parents, DiscRV* drv)
{
  logpr rc;
  if ( _latticeAdt->useTimeParent() ) {
    // see below for documentation.

    unsigned lat_time = 
      (unsigned)round(_latticeAdt->_frameRate * 
		      _latticeAdt->_latticeNodes[RV2DRV(parents[0])->val].time
		      );

    // case on the current time.
    // Same as no time parent case, we also allow some relaxation
    // of time.
    if (_latticeAdt->_timeMarks && RV2DRV(parents[2])->val < lat_time - _latticeAdt->_frameRelax ) {
      // then the current time is less than the previous lattice node.
      if ( RV2DRV(parents[1])->val ) {
	// a (word) transition is being hypothesized, but we do not allow it. Give it a
	// zero probability.
	rc.set_to_zero();
      } else {
	// a (word) transition is not being hypothesized, so
	// we copy the previous lattice node's value (from previous frame) to drv.
	if (drv->val == RV2DRV(parents[0])->val)
	  rc.set_to_one();
	else 
	  rc.set_to_zero();
      }
    } else if (_latticeAdt->_timeMarks && RV2DRV(parents[2])->val > lat_time + _latticeAdt->_frameRelax ) {
      // Then, the current time (i.e., parent[2]'s time value) is
      // already ahead (i.e., after, later) of when a transition for
      // the previous lattice node value may occur. We also want to
      // force this not to happen by setting prob to zero. While
      // you might think this might be valid, only allow jumping
      // from the prevous lattice node when the time is exactly
      // equal to the previous lattice nodes.
      rc.set_to_zero();
    } else {
      // In range for one reason or another.

      // (RV2DRV(parents[2])->val == lat_time), which means that the
      // current time is right at (or in the range of) the point that
      // we allow the previous lattice node to jump to the next set of
      // possible lattice nodes.

      // Next, we check that the 'transition' variable (word
      // transition in the 2006 paper) is set, and it is stored in
      // parent[1]. Note we assume that parent[1] is a binary
      // variable, and use its value as a boolean int.
      if ( RV2DRV(parents[1])->val ) {
	// Then a transition is being asked for.

	// iterate next lattice nodes find the out going edge
	LatticeADT::LatticeNode &node = _latticeAdt->_latticeNodes[RV2DRV(parents[0])->val];

	// if the lattice node is the end, return prob zero
	if ( node.edges.totalNumberEntries() == 0 ) {
	  rc.set_to_zero();
	} else {
	  LatticeADT::LatticeEdgeList* outEdges
	    = _latticeAdt->_latticeNodes[RV2DRV(parents[0])->val].edges.find(drv->val);
	  if ( outEdges == NULL )
	    rc.set_to_zero();
	  else {
	    if (LatticeADT::_latticeNodeUseMaxScore) {
	      rc = outEdges->max_gmtk_score;
	    } else {
	      rc.set_to_one();
	    }
	  }
	}
      } else {
	// No word transition is being asked for, so we copy the value
	// value from the parent lattice node to the child node in
	// accordance with the event that no transition occurs.

	if (drv->val == RV2DRV(parents[0])->val)
	  rc.set_to_one();
	else 
	  rc.set_to_zero();
      }
    }
  } else {
    // In this case, we are not using the time parent, and we assume
    // that the time value comes from 'drv'.

    // Note: we don't adjust for frame relax here as that's already
    // been done in LatticeADT::resetFrameIndices(unsigned numFrames).

    // case on the current frame index
    if ( _latticeAdt->_timeMarks && drv->frame() < _latticeAdt->_latticeNodes[RV2DRV(parents[0])->val].startFrame ) {
      if ( RV2DRV(parents[1])->val ) {
	// a (word) transition is being hypothesized, but we do not allow it. Give it a
	// zero probability.
	rc.set_to_zero();
      } else {
	// a (word) transition is not being hypothesized, so
	// we copy the previous lattice node's value (from previous frame) to drv.
	if (drv->val == RV2DRV(parents[0])->val)
	  rc.set_to_one();
	else
	  rc.set_to_zero();
      }
    } else if ( _latticeAdt->_timeMarks && drv->frame() > _latticeAdt->_latticeNodes[RV2DRV(parents[0])->val].endFrame ) {
      // we force this cannot happen by setting prob to zero
      rc.set_to_zero();
    } else {
      // in range.

      if ( RV2DRV(parents[1])->val ) {
	// find the out going edge
	LatticeADT::LatticeNode &node = _latticeAdt->_latticeNodes[RV2DRV(parents[0])->val];
	// if the lattice node is the end, return prob zero
	if ( node.edges.totalNumberEntries() == 0 ) {
	  rc.set_to_zero();
	} else {
	  LatticeADT::LatticeEdgeList* outEdges
	    = _latticeAdt->_latticeNodes[RV2DRV(parents[0])->val].edges.find(drv->val);
	  if ( outEdges == NULL )
	    rc.set_to_zero();
	  else {
	    if (LatticeADT::_latticeNodeUseMaxScore) {
	      rc = outEdges->max_gmtk_score;
	    } else {
	      rc.set_to_one();
	    }
	  }
	}
      } else {
	// No word transition is being asked for, so we score a copy
	// of the value value from the parent lattice node to the
	// child node in accordance with the event that no
	// transition occurs.
	if (drv->val == RV2DRV(parents[0])->val)
	  rc.set_to_one();
	else 
	  rc.set_to_zero();
      }
    }
  }
  return rc;
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
void LatticeNodeCPT::becomeAwareOfParentValuesAndIterBegin(vector< RV* >& parents, 
							   iterator &it, 
							   DiscRV* drv, 
							   logpr& p) 
{
  // 
  // Note that there can be either 2 or 3 parents:
  // 
  // In both cases: 
  //    drv is the current lattice node.
  // 
  // In the case of 2 parents, (where time is normally obtained from 'drv', the current lattice node):
  //    parent[0] is the previous lattice node
  //    parent[1] is the  "word transition" (or variable that is acting like such a construct)
  // 
  // In the case of 3 parents, 
  //    parent[0] is the previous lattice node
  //    parent[1] is the the "word transition" (or variable that is acting like such a construct)
  //    parent[2] is the "time observation", namely it is a variable that is presumably observed that
  //              keeps track of the time frame to use (rather than using the time frame of 'drv').
  // 
  // For simplicity, the cardinality of lattice nodes can be bigger
  // than number of real nodes in some particular lattice.  This is
  // because in iterable lattices, different lattices can have
  // different number of nodes.  But in master file, we can just
  // specify the max of those.

  // initialize it to something always at least valid.
  drv->val = 0;

  if ( _latticeAdt->useTimeParent() )
  {
    // use time parent to check time. In this case,
    // _ignoreLatticeNodeTimeMarks has no effect.


    // Since we're using a time parent, the number of frames in the
    // observation file does not indicate the number of true frames of
    // the segment, so we are not able to compute a frame rate.  We do
    // the node time relaxation inline in the below.  See
    // LatticeADT::resetFrameIndices(unsigned numFrames) for the case
    // of when this adjustment is done without a time parent.


    // First, compute lat_time, the time of the previous lattice node
    // rounded to the closest frame. parent[0] contains the value of
    // the previous lattice node, and we need to do a lookup in the
    // lattce to find the actual previous lattice node to get its
    // time.
    unsigned lat_time = 
      (unsigned)round(_latticeAdt->_frameRate * 
		      _latticeAdt->_latticeNodes[RV2DRV(parents[0])->val].time
		      );



    // case on the current time.
    // Same as no time parent case, we also allow some relaxation
    // of time.
    if (_latticeAdt->_timeMarks && RV2DRV(parents[2])->val < lat_time - _latticeAdt->_frameRelax ) {
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
    } 



    else if (_latticeAdt->_timeMarks && RV2DRV(parents[2])->val > lat_time + 1000){ //_latticeAdt->_frameRelax ) {
      // Then, the current time (i.e., parent[2]'s time value) is
      // already ahead (i.e., after, later) of when a transition for
      // the previous lattice node value may occur. We also want to
      // force this not to happen by setting prob to zero. While
      // you might think this might be valid, only allow jumping
      // from the prevous lattice node when the time is exactly
      // equal to the previous lattice nodes.

      it.internalStatePtr = NULL;
      p.set_to_zero();

    } 


    else {
      // In range for one reason or another.

      // (RV2DRV(parents[2])->val == lat_time), which means that the
      // current time is right at (or in the range of) the point that
      // we allow the previous lattice node to jump to the next set of
      // possible lattice nodes.

      // Next, we check that the 'transition' variable (word
      // transition in the 2006 paper) is set, and it is stored in
      // parent[1]. Note we assume that parent[1] is a binary
      // variable, and use its value as a boolean int.
      if ( RV2DRV(parents[1])->val ) {
	// Then a transition is being asked for.

	// iterate next lattice nodes find the out going edge
	LatticeADT::LatticeNode &node = _latticeAdt->_latticeNodes[RV2DRV(parents[0])->val];

	// if the lattice node is the end, return prob zero
	if ( node.edges.totalNumberEntries() == 0 ) {
	  it.internalStatePtr = NULL;
	  p.set_to_zero();
	  return;
	}

	// check out the out-going edges based on parent value
	shash_map_iter<unsigned, LatticeADT::LatticeEdgeList>::iterator *pit 
	  = new shash_map_iter<unsigned, LatticeADT::LatticeEdgeList>::iterator();


    unsigned delta = RV2DRV(parents[1])->val;

    unsigned node_id = RV2DRV(parents[0])->val;

/*
    unsigned target_id = node_id + delta;

    if(target_id > 12) {
        it.internalStatePtr = NULL;
        p.set_to_zero();
        return;
    }

*/


    //printf("node id: %u, delta: %u, cpt size: %u, max_delta: %u\n", node_id, delta, node.larger_cpts.size(), node.max_delta);

    if(delta >= node.max_delta || node.larger_cpts[delta] == NULL) {
        it.internalStatePtr = NULL;
        p.set_to_zero();
        return;
    }




    //LatticeADT::LatticeNode &node_before = _latticeAdt->_latticeNodes[target_id - 1];

    //if(node.larger_cpts[delta].totalNumberEntries() <= 0) {
    

    /*
    bool new_cpt = false;
    if(node.larger_cpts[delta] == NULL) {

        new_cpt = true;
        node.larger_cpt = node.larger_cpts[delta];

        printf("one larger cpt size: %u\n", node.larger_cpt->totalNumberEntries());
    }
    */

    
    //node.larger_cpt = &node.larger_cpts[delta];
    node.larger_cpt = node.larger_cpts[delta];

    //node.larger_cpt = &node_before.edges;
    


    node.larger_cpt->begin(*pit);
        



	// now the current node can have next transition
	// find the correct iterators
	//node.edges.begin(*pit);

	// set up the internal state for iterator
	it.internalStatePtr = (void*)pit;
	it.internalState = RV2DRV(parents[0])->val;
	it.drv = drv;

	drv->val = pit->key();

	if (LatticeADT::_latticeNodeUseMaxScore) {
	  p = (**pit).max_gmtk_score;
	} else {
	  // then the various lattice edges handle the scoring
	  p.set_to_one();
	}
      } else {

	// No word transition is being asked for, so we copy the value
	// value from the parent lattice node to the child node in
	// accordance with the event that no transition occurs.

	it.internalStatePtr = NULL;
	drv->val = RV2DRV(parents[0])->val;
	p.set_to_one();

      }
    }
  } else {
    // In this case, we are not using the time parent, and we assume
    // that the time value comes from 'drv'.

    // Note: we don't adjust for frame relax here as that's already
    // been done in LatticeADT::resetFrameIndices(unsigned numFrames).

    // case on the current frame index
    if ( _latticeAdt->_timeMarks && drv->frame() < _latticeAdt->_latticeNodes[RV2DRV(parents[0])->val].startFrame ) {
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
    } else if ( _latticeAdt->_timeMarks && drv->frame() > _latticeAdt->_latticeNodes[RV2DRV(parents[0])->val].endFrame ) {
      // we force this cannot happen by setting prob to zero
      it.internalStatePtr = NULL;
      p.set_to_zero();
    } else {
      // in range.

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
        shash_map_iter<unsigned, LatticeADT::LatticeEdgeList>::iterator *pit 
	  = new shash_map_iter<unsigned, LatticeADT::LatticeEdgeList>::iterator();



        //TODO: create a list of hypothesized nodes, corresponding to DELTA (jump amount)
        //Let parents[1] be the value of DELTA
        //Note that for now time frames don't work for this case, and edge weights are not considered either

        unsigned delta = RV2DRV(parents[1])->val;

        if(delta >= node.max_delta || node.larger_cpts[delta] == NULL) {
            it.internalStatePtr = NULL;
            p.set_to_zero();
            return;
        }

        node.larger_cpt = node.larger_cpts[delta];        
        node.larger_cpt->begin(*pit);



        // now the current node can have next transition
        // find the correct iterators

        //node.edges.begin(*pit);

        // set up the internal state for iterator
        it.internalStatePtr = (void*)pit;
        it.internalState = RV2DRV(parents[0])->val;
        it.drv = drv;

        drv->val = pit->key();

	if (LatticeADT::_latticeNodeUseMaxScore) {
	  p = (**pit).max_gmtk_score;
	} else {
	  // then the various lattice edges handle the scoring
	  p.set_to_one();
	}

      } else {

	// No word transition is being asked for, so we copy the value
	// value from the parent lattice node to the child node in
	// accordance with the event that no transition occurs.

        it.internalStatePtr = NULL;
        drv->val = RV2DRV(parents[0])->val;
        p.set_to_one();
      }
    }
  }
}


/*

void LatticeNodeCPT::addMap2Map(shash_map_iter<unsigned, LatticeADT::LatticeEdgeList> & map1, shash_map_iter<unsigned, LatticeADT::LatticeEdgeList> & map2) {
    shash_map_iter<unsigned, LatticeADT::LatticeEdgeList>::iterator it;
    map1.begin(it);
    
    do {
        LatticeADT::LatticeEdgeList outEdge = (LatticeADT::LatticeEdgeList) (*it);
        unsigned child_id = it.key();

        for(unsigned k=0; k<outEdge.edge_array.size(); k++) {

            addNodeEdge2Map(map2, child_id, outEdge.edge_array[k]);
        }
        
    } while(it.next());
}



void LatticeNodeCPT::dfsOnNode(LatticeADT::LatticeNode & cur_node, unsigned delta, shash_map_iter<unsigned, LatticeADT::LatticeEdgeList> & result) {

    if(cur_node.edges.totalNumberEntries() <= 0) return;

    shash_map_iter<unsigned, LatticeADT::LatticeEdgeList>::iterator it;

    //printf("dfs num entries: %d\n", cur_node.edges.totalNumberEntries());
    cur_node.edges.begin(it);
    
    //TODO: return_structure to include delta information?

     do {
        LatticeADT::LatticeEdgeList outEdge = (LatticeADT::LatticeEdgeList) (*it);
        unsigned child_id = it.key();

        //printf("\tloop on child_id %u: ", child_id);

        //For each edge to child_node
        for(unsigned k=0; k<outEdge.edge_array.size(); k++) {
            unsigned length = outEdge.edge_array[k].length;

            LatticeADT::LatticeNode &child_node = _latticeAdt->_latticeNodes[child_id];
    
            //we can explore more nodes
            if(delta > length) {
                dfsOnNode(child_node, delta-length, result);
                //queryLargerCPT(child_node, delta-length, result);
            }
            //we cannot go further beyond the child node
            else {
                addNodeEdge2Map(result, child_id, outEdge.edge_array[k]);
            }
        }
    } while(it.next());


}



void LatticeNodeCPT::addNodeEdge2Map(shash_map_iter<unsigned, LatticeADT::LatticeEdgeList> & added, unsigned child_id, const LatticeADT::LatticeEdge & edge) {

    LatticeADT::LatticeEdgeList* edge_listp = added.find(child_id);

    //printf(" adding child id %u\n", child_id);
    
    if (edge_listp == NULL) {
        //printf("insert new\n");

        // need to insert new edge list
        LatticeADT::LatticeEdgeList edge_list;
        // start with only one edge (conserve memory)
        edge_list.edge_array.resize(1);
        // place edge (using copy).
        edge_list.edge_array[0] = edge;
        // record edge.
        edge_list.num_edges  = 1;
        // and lastly, insert it in the hash table
        added.insert(child_id, edge_list);

    } else {
        // edge list already there, add to the end.
        assert ( edge_listp->num_edges > 0 );

        for(unsigned i=0; i<edge_listp->num_edges; i++) {
            if(&edge_listp->edge_array[i] == &edge) {
                printf("same");
                return;
            }
        }

        // make sure there is room for at least one more, use a
        // conservative growth rate
        const float growth_rate = 1.25;
        edge_listp->edge_array.growByFIfNeededAndCopy(growth_rate,
						    edge_listp->num_edges+1);

        // place edge (using copy).
        edge_listp->edge_array[edge_listp->num_edges] = edge;
        edge_listp->num_edges++;
    }
}



void LatticeNodeCPT::setGMTKScoresOnLargerCPT(LatticeADT::LatticeNode & cur_node) {

    if(cur_node.larger_cpt->totalNumberEntries() > 0) {
        shash_map_iter<unsigned, LatticeADT::LatticeEdgeList>::iterator it;
        cur_node.larger_cpt->begin(it);
        
        do {
            LatticeADT::LatticeEdgeList	&edge_list = (*it);
	        edge_list.max_gmtk_score.set_to_zero();

            for (unsigned edge_ctr=0;edge_ctr < edge_list.num_edges; edge_ctr ++ ) {
	              LatticeADT::LatticeEdge &edge = edge_list.edge_array[edge_ctr];

	            if (edge.gmtk_score > edge_list.max_gmtk_score) 
	                edge_list.max_gmtk_score = edge.gmtk_score;
            }

            if (LatticeADT::_latticeNodeUseMaxScore) {
                // need to renormize the edge scores.
                for (unsigned edge_ctr=0;edge_ctr < edge_list.num_edges; edge_ctr ++ ) {
                    LatticeADT::LatticeEdge &edge = edge_list.edge_array[edge_ctr];
                    edge.gmtk_score = edge.gmtk_score / edge_list.max_gmtk_score;
                }
	        } else {
	            // need to reset the max scores to 1 so that they are
	            // effectively not used.
	            edge_list.max_gmtk_score.set_to_one();
            }

        } while(it.next());
    }
}
*/

void LatticeNodeCPT::printLargerCPT(shash_map_iter<unsigned, LatticeADT::LatticeEdgeList> & cpt, LatticeADT::LatticeNode & cur_node, unsigned cur_node_id) {

    if(cpt.totalNumberEntries() <= 0) return;

    shash_map_iter<unsigned, LatticeADT::LatticeEdgeList>::iterator it;

    cpt.begin(it);
    
    printf("Printing Larger CPT for node: %u\n", cur_node_id);
    //printf("num entries: %d\n", cur_node.edges.totalNumberEntries());

    do {
        LatticeADT::LatticeEdgeList outEdge = (LatticeADT::LatticeEdgeList) *it;
        unsigned child_id = it.key();

        printf("\t%u\n", child_id);
    } while(it.next());

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
   shash_map_iter<unsigned, LatticeADT::LatticeEdgeList>::iterator* pit 
    = (shash_map_iter<unsigned, LatticeADT::LatticeEdgeList>::iterator*) it.internalStatePtr;

  // check whether pit is null
  if ( pit == NULL ) {
    p.set_to_zero();
    return false;
  }

  // No need to check the time constrains in the lattice
  // transitions.  The reason is once the time passed
  // checking in BeginIterate, it already satisfies the constrains.

  // find the next available in the tree
  if ( pit->next() ) {

    // set up the values
    it.drv->val = pit->key();

    if (LatticeADT::_latticeNodeUseMaxScore) {
      p = (**pit).max_gmtk_score;
    } else {
      p.set_to_one();
    }

    return true;
  } else {
    // we're done with all the next nodes for this current node
    // so we free things up.
    delete pit;
    // and return that we're done.
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

  // in addtiont to setting lattice ADT, also need to set
  // up the number of parents and cardinalities of the parents.

  // typically, cardinalities comes from structure file or master
  // file.  But in this case, we hope to support iterable lattice
  // CPTs which will have different number of lattice nodes for
  // each one.

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

    //To enable multiple jumps, we can have arbitary positive value for transition
    //TODO: interface with user and check errors
    cardinalities[1] = _latticeAdt->_transitionCardinality;

    printf("Transition Cardinality: %u\n", cardinalities[1]);

}
