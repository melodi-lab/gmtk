/*-
 * GMTK_JunctionTree.cc
 *     Junction Tree
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
#include <new>

#include "general.h"
#include "error.h"
#include "debug.h"
#include "rand.h"

#include "GMTK_FileParser.h"
#include "GMTK_RV.h"
#include "GMTK_DiscRV.h"
#include "GMTK_GMTemplate.h"
#include "GMTK_JunctionTree.h"
#include "GMTK_GMParms.h"
#include "GMTK_PackCliqueValue.h"

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Static variables used by classes
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

VCID("$Header$")


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Partition structure support
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////



PartitionStructures::PartitionStructures(JT_Partition& from_part,
					 vector <RV*>& newRvs,
					 map < RVInfo::rvParent, unsigned >& ppf,
					 const unsigned int frameDelta,
					 const bool has_li_separator)
  : origin(from_part)
{

  // first allocate space with empty (and unusable) entries
  maxCliquesSharedStructure.resize(origin.cliques.size());
  separatorCliquesSharedStructure.resize(origin.separators.size());

  // then actually re-construct the objects in the array appropriately.
  for (unsigned i=0;i<maxCliquesSharedStructure.size();i++) {
    new (&maxCliquesSharedStructure[i]) 
      MaxCliqueTable::SharedLocalStructure(origin.cliques[i],
					   newRvs,ppf,frameDelta);
  }
  for (unsigned i=0;i<separatorCliquesSharedStructure.size();i++) {
    new (&separatorCliquesSharedStructure[i]) 
      ConditionalSeparatorTable::SharedLocalStructure(origin.separators[i],
						      newRvs,ppf,frameDelta);
  }


#if 0
  for (unsigned i=0;i<factorCliques.size();i++) {
    new (&factorCliques[i]) InferenceFactorClique(origin.factorCliques[i],
						  newRvs,ppf,frameDelta);
  }
#endif

  // grab the union of all RVs 
  allrvs.clear();
  for (unsigned i=0;i<maxCliquesSharedStructure.size();i++) {
    set<RV*> cur_set = maxCliquesSharedStructure[i].returnRVsAsSet();
    copy(cur_set.begin(),cur_set.end(),
	 inserter(allrvs,allrvs.end()));
  }
  // Since the final LI separator is guaranteed to be in the last
  // posistion of the separator array, we can just subtract that out (if it
  // exists, but a passed in argument tells us that).

  if (has_li_separator) {
    set <RV*> lirvs = 
      separatorCliquesSharedStructure[separatorCliquesSharedStructure.size()-1].returnRVsAsSet();
    set <RV*> res;
    set_difference(allrvs.begin(),allrvs.end(),
		   lirvs.begin(),lirvs.end(),
		   inserter(res,res.end()));
    allrvs = res;
  }

  if (JunctionTree::viterbiScore == true) {

    // set up a few members that are needed for computing and storing
    // values of this partition.

    // Now, store a copy of all the random variables in this
    // partition, minus the left interface separator. This needs to be
    // stored in an array (i.e., it needs to be ordered) so that we can
    // pack the variales in a partition all at once, rather than for
    // each clique.
    // 
  
    // now need to create the packer for the hidden variables.
    hidRVVector.clear();
    for (set <RV*>::iterator it = allrvs.begin();
	 it != allrvs.end(); it++) {
      RV* rv = (*it);
      // TODO: could save even more space by storing only the values
      // that can't be derived from the other values!!
      if (rv->hidden())
	hidRVVector.push_back(rv);
    }
    // printf("%d %d\n",allrvs.size(),hidRVVector.size());
    if (hidRVVector.size() > 0) {
      new (&packer) PackCliqueValue(hidRVVector);
  
      hrvValuePtrs.resize(hidRVVector.size());
      for (unsigned i=0; i<hidRVVector.size(); i++) {
	// hidden nodes are always discrete (in this version).
	DiscRV* drv = 
	  (DiscRV*)hidRVVector[i];
	hrvValuePtrs[i] = &(drv->val);
      }
    } else {
      // we are going to need to check the zero size condition
      // when we use the packer since it needs strict > 0
    }
  }

}


/*
 * return RVS and any of their observed parents as a set.
 *
 */
set <RV*>
PartitionStructures::returnRVsAndTheirObservedParentsAsSet()
{
  // TODO: do this once and cache it.
  // grab the union of all RVs in part1 and part2.
  set<RV*> rc;
  for (unsigned i=0;i<maxCliquesSharedStructure.size();i++) {
    set<RV*> cur_set = maxCliquesSharedStructure[i].returnRVsAndTheirObservedParentsAsSet();
    copy(cur_set.begin(),cur_set.end(),
	 inserter(rc,rc.end()));
  }
  return rc;
}


#if 0  

/*-
 *-----------------------------------------------------------------------
 * PartitionStructures::adjustFramesBy()
 *
 *    Go through each clique in the partition and update the frame numbers of
 *    each clique by the adjustment (which could be positive or negative) by
 *    the amoun given as the argument.
 *
 *       
 *
 * Preconditions:
 *    The cliques, separators, and facor cliques must be an
 *    valid instantiated table. All data structures must be set up.
 *
 * Postconditions:
 *   Frames of all cliques, separators, and factor cliques are changed.
 *
 * Side Effects:
 *   This will update the internal random variables of the underlying
 *   variables held by the cliques.
 *
 * Results:
 *   None
 *
 *-----------------------------------------------------------------------
 */
void PartitionStructures::adjustFramesBy(PartitionStructures& part1,
					   PartitionStructures& part2,
					   const int frameDelta)
{


  // grab the union of all RVs in part1 and part2.
  set<RV*> all_rvs;
  for (unsigned i=0;i<part1.maxCliques.size();i++) {


    maxCliques[i].adjustFrameBy(frameDelta);
  }

	    set_union(empty.begin(),empty.end(),
		      sep_j.nodes.begin(),sep_j.nodes.end(),
		      inserter(sep_union_set,sep_union_set.end()));  


union_1_2_to_3(const set<RV*>& A,
	       const set<RV*>& B,
	       set<RV*>& C,
	       bool do_not_clear = false)  


  // Note that since separators and factor cliques are contained in
  // maxcliques, we need only to adjust the cliques. Also, we assume
  // here that each of a maxClique, separator, and factorclique is
  // such that the frameDelta changes only the set of random variables
  // associated with it.

  // should just grab the un


}

#endif					   



/////////////////////////////////////////////	
/// END OF FILE
/////////////////////////////////////////////
