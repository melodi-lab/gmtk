/*-
 * GMTK_GMTemplate.cc
 *     manipulations of a GM template
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2001, < fill in later >
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

#include <iterator>
#include <map>
#include <set>

#include "general.h"
#include "error.h"
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
#include "GMTK_MixGaussians.h"
#include "GMTK_ObservationMatrix.h"

VCID("$Header$");

////////////////////////////////////////////////////////////////////
//        General create, read, destroy routines 
////////////////////////////////////////////////////////////////////



/*-
 *-----------------------------------------------------------------------
 * GMTemplate::dropEdgeDirections()
 *      Create an undirected graph in 'rvs' by essentially
 *      dropping the edge directions (i.e., the graph is the same
 *      except it is not directed)
 *
 * Preconditions:
 *      rvs must contain a valid graph
 *
 * Postconditions:
 *      rvs now represent an undirected graph
 *
 * Side Effects:
 *     none
 *
 * Results:
 *     returns true if everything works, return false if something
 *     goes wrong.
 *
 *-----------------------------------------------------------------------
 */
bool
GMTemplate::dropEdgeDirections()
{
  for (unsigned i=0;i<rvs.size();i++) {
    rvs[i]->createNeighborsFromParentsChildren();
  }
  return true;
}


/*-
 *-----------------------------------------------------------------------
 * GMTemplate::moralize()
 *      Moralize the graph that currently resides in member variable rvs.
 *
 * Preconditions:
 *      - dropEdgeDirections() must be called first.
 *      - rvs must contain a valid graph
 *
 * Postconditions:
 *      graph is moralized.
 *
 * Side Effects:
 *     random variables in rvs are changed.
 *
 * Results:
 *     returns true if everything works, return false if something
 *     goes wrong.
 *
 *-----------------------------------------------------------------------
 */
bool
GMTemplate::moralize()
{
  for (unsigned i=0;i<rvs.size();i++) {
    rvs[i]->moralize();    
  }
  return true;
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
GMTemplate::makeComplete(set<RandomVariable*> &rvs)
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
    (*i)->neighbors = res;
  }

  return;
}



/*-
 *-----------------------------------------------------------------------
 * GMTemplate::findBestLeftInterface()
 *   NOTE: If this routine gets updated, so should findBestRightInterface()
 *   Given a once unrolled graph, P,C1,C2,C3,E
 *   find the best left interface within C2 starting at the
 *   "standard" or initial left interface between C1 and C2.
 *   Note that the routine only uses C1,C2, and C3 and
 *   P and E are not needed.
 *   
 *
 * Preconditions:
 *     Graph must be valid (i.e., unroller should pass graph w/o
 *                                problem).
 *     Graph must be already moralized.
 *     and *must* unrolled exactly two times.
 *     This means graph must be in the form P,C1,C2,C3,E
 *     where P = prologue, 
 *           C1 = first chunk
 *           C2 = 2nd chunk
 *           C3 = 3nd chunk
 *           E =  epilogue
 *
 *     Note that P or E (but not both) could be empty (this
 *     is why we do this procedure on the unrolled-by-2 graph
 *     rather than on the unrolled-by- 0 or 1 graph).
 *
 * Postconditions:
 *     Return values have best left-interface found.
 *
 * Side Effects:
 *      None
 *
 *
 * Results:
 *    Put best left interface from C2 into C_l
 *    and place any additional variables from C2 to the left of C_l
 *    into 'left_C_l';
 *
 *
 *-----------------------------------------------------------------------
 */
void
GMTemplate::findBestLeftInterface(const set<RandomVariable*> &C1,
				  const set<RandomVariable*> &C2,
				  const set<RandomVariable*> &C3,
				  set<RandomVariable*> &left_C_l,
				  set<RandomVariable*> &C_l)
{

  // left interface
  C_l.clear();

  // First, construct the basic left interface (i.e., left
  // interface of C2).

  // go through through set C1, and pick out all neighbors
  // of variables in set C1 that live in C2.
  set<RandomVariable*>::iterator c1_iter;
  for (c1_iter = C1.begin(); c1_iter != C1.end(); c1_iter ++) {
    // go through all neighbors of nodes in C1
    RandomVariable *cur_rv = (*c1_iter);
    // neighbor iterator
    set<RandomVariable*>::iterator n_iter;
    for (n_iter = cur_rv->neighbors.begin();
	 n_iter != cur_rv->neighbors.end();
	 n_iter ++) {
      if (C2.find((*n_iter)) != C2.end()) {
	// found a neighbor of cur_rv in C2, so
	// it must be in C_l
	C_l.insert((*n_iter));
      }
    }
  }
  left_C_l.clear();
  
  printf("Size of basic left interface C_l = %d\n",C_l.size());
  printf("Size of left_C_l = %d\n",left_C_l.size());
  {
    printf("Left interface nodes include:");
    set<RandomVariable*>::iterator i;    
    for (i=C_l.begin();i!=C_l.end();i++) {
      printf(" %s(%d)",
	     (*i)->name().c_str(),
	     (*i)->frame());
	     
    }
    printf("\n");
  }


  // best ones found so far
  set<RandomVariable*> best_left_C_l = left_C_l;
  set<RandomVariable*> best_C_l = C_l;
  set< set<RandomVariable*> > setset;

  findBestLeftInterface(left_C_l,
			C_l,
			C2,
			C3,
			setset,
			best_left_C_l,
			best_C_l);

  
  printf("Size of best left interface = %d\n",best_C_l.size());
  printf("Size of best_left_C_l = %d\n",best_left_C_l.size());
  {
    printf("Best left interface nodes include:");
    set<RandomVariable*>::iterator i;    
    for (i=best_C_l.begin();i!=best_C_l.end();i++) {
      printf(" %s(%d)",
	     (*i)->name().c_str(),
	     (*i)->frame());
    }
    printf("\n");

  }

  left_C_l = best_left_C_l;
  C_l = best_C_l;

}


/*-
 *-----------------------------------------------------------------------
 * GMTemplate::findBestLeftInterface()
 *    recursive helper function for the first call findBestLeftInterface()
 *    See that routine for documentation.
 *
 * Preconditions:
 *
 * Postconditions:
 *
 * Side Effects:
 *
 * Results:
 *
 *-----------------------------------------------------------------------
 */
void
GMTemplate::
findBestLeftInterface(set<RandomVariable*> &left_C_l,
		      set<RandomVariable*> &C_l,
		      const set<RandomVariable*> &C2,
		      const set<RandomVariable*> &C3,
		      set< set<RandomVariable*> >& setset,
		      set<RandomVariable*> &best_left_C_l,
		      set<RandomVariable*> &best_C_l)
{
  set<RandomVariable*>::iterator v;  // vertex
  // consider all v in the current C_l as candidates
  // to be moved left.
  for (v = C_l.begin(); v != C_l.end(); v ++) {
    // TODO: rather than "for all nodes in C_l", we could
    // do a random subset of nodes in C_l to speed this up if
    // it takes too long. But note that this is only run once
    // per graph so it will be beneficial to do this since
    // its cost is ammortized over the many runs of the graph

    // if v has neighbors in C3
    // next;
    // do a set intersection.
    set<RandomVariable*> res;
    set_intersection((*v)->neighbors.begin(),
		     (*v)->neighbors.end(),
		     C3.begin(),C3.end(),
		     inserter(res, res.end()));
    if (res.size() != 0)
      continue;

    // take v from C_l and place it in left_C_l    
    set<RandomVariable*> next_left_C_l = left_C_l;
    next_left_C_l.insert((*v));

    // and add all neighbors of v that are 
    // in C2\next_left_C_l to next_C_l
    set<RandomVariable*> next_C_l;
    set<RandomVariable*> tmp;
    set_intersection((*v)->neighbors.begin(),
		     (*v)->neighbors.end(),
		     C2.begin(),C2.end(),
		     inserter(tmp,tmp.end()));
    res.clear();    
    set_difference(tmp.begin(),tmp.end(),
		   next_left_C_l.begin(),next_left_C_l.end(),
		   inserter(res,res.end()));

    set_union(res.begin(),res.end(),
	      C_l.begin(),C_l.end(),
	      inserter(next_C_l,next_C_l.end()));
    next_C_l.erase((*v));

    if (setset.find(next_C_l) != setset.end())
      continue; // check if memoized, if so, no need to go further.
    
    // check size of candiate interface, and keep a copy of
    // it if it is less then the one we have seen so far.
    if (next_C_l.size() < best_C_l.size()) {
      best_left_C_l = next_left_C_l;
      best_C_l = next_C_l;
    } 

    // memoize
    setset.insert(next_C_l);

    findBestLeftInterface(next_left_C_l,
			  next_C_l,
			  C2,C3,setset,
			  best_left_C_l,best_C_l);

  }
}




/*-
 *-----------------------------------------------------------------------
 * GMTemplate::findBestRightInterface()
 *   just like findBestLeftInterface() but looks for the right interface.
 *   NOTE: If this routine gets updated, so should findBestLeftInterface()
 *
 * Preconditions:
 *   Original graph *must* be unrolled 2 times for this to work.
 *
 * Postconditions:
 *
 * Side Effects:
 *
 *
 * Results:
 *    Put best right interface from C2 into C_r
 *    and place any additional variables from C2 to the right of C_r
 *    into 'right_C_r';
 *
 *
 *-----------------------------------------------------------------------
 */
void
GMTemplate::findBestRightInterface(const set<RandomVariable*> &C1,
				   const set<RandomVariable*> &C2,
				   const set<RandomVariable*> &C3,
				   set<RandomVariable*> &right_C_r,
				   set<RandomVariable*> &C_r)
{


  // right interface
  C_r.clear();

  // First, construct the basic right interface (i.e., right
  // interface of C2).

  // go through through set C3, and pick out all neighbors
  // of variables in set C3 that live in C2.
  set<RandomVariable*>::iterator c3_iter;
  for (c3_iter = C3.begin(); c3_iter != C3.end(); c3_iter ++) {
    // go through all neighbors of nodes in C1
    RandomVariable *cur_rv = (*c3_iter);
    // neighbor iterator
    set<RandomVariable*>::iterator n_iter;
    for (n_iter = cur_rv->neighbors.begin();
	 n_iter != cur_rv->neighbors.end();
	 n_iter ++) {
      if (C2.find((*n_iter)) != C2.end()) {
	// found a neighbor of cur_rv in C2, so
	// it must be in C_r
	C_r.insert((*n_iter));
      }
    }
  }
  right_C_r.clear();
  
  printf("Size of basic right interface C_r = %d\n",C_r.size());
  printf("Size of right_C_r = %d\n",right_C_r.size());
  {
    printf("Right interface nodes include:");
    set<RandomVariable*>::iterator i;    
    for (i=C_r.begin();i!=C_r.end();i++) {
      printf(" %s(%d)",
	     (*i)->name().c_str(),
	     (*i)->frame());
	     
    }
    printf("\n");
  }


  // best ones found so far
  set<RandomVariable*> best_right_C_r = right_C_r;
  set<RandomVariable*> best_C_r = C_r;
  set< set<RandomVariable*> > setset;

  findBestRightInterface(right_C_r,
			 C_r,
			 C2,
			 C1,
			 setset,
			 best_right_C_r,
			 best_C_r);

  
  printf("Size of best right interface = %d\n",best_C_r.size());
  printf("Size of best_right_C_r = %d\n",best_right_C_r.size());
  {
    printf("Best right interface nodes include:");
    set<RandomVariable*>::iterator i;    
    for (i=best_C_r.begin();i!=best_C_r.end();i++) {
      printf(" %s(%d)",
	     (*i)->name().c_str(),
	     (*i)->frame());
    }
    printf("\n");

  }    

  right_C_r = best_right_C_r;
  C_r = best_C_r;

}




/*-
 *-----------------------------------------------------------------------
 * GMTemplate::findBestRightInterface()
 *    recursive helper function for the first call findBestRightInterface()
 *    See that routine for documentation.
 *
 * Preconditions:
 *
 * Postconditions:
 *
 * Side Effects:
 *
 * Results:
 *
 *-----------------------------------------------------------------------
 */
void
GMTemplate::
findBestRightInterface(set<RandomVariable*> &right_C_r,
		       set<RandomVariable*> &C_r,
		       const set<RandomVariable*> &C2,
		       const set<RandomVariable*> &C1,
		       set< set<RandomVariable*> >& setset,
		       set<RandomVariable*> &best_right_C_r,
		       set<RandomVariable*> &best_C_r)
{
  set<RandomVariable*>::iterator v;  // vertex
  // consider all v in the current C_r as candidates
  // to be moved right.
  for (v = C_r.begin(); v != C_r.end(); v ++) {

    // if v has neighbors in C1
    // next;
    // do a set intersection.
    set<RandomVariable*> res;
    set_intersection((*v)->neighbors.begin(),
		     (*v)->neighbors.end(),
		     C1.begin(),C1.end(),
		     inserter(res, res.end()));
    if (res.size() != 0)
      continue;

    // take v from C_r and place it in right_C_r
    set<RandomVariable*> next_right_C_r = right_C_r;
    next_right_C_r.insert((*v));

    // and add all neighbors of v that are 
    // in C2\next_right_C_r to next_C_r
    set<RandomVariable*> next_C_r;
    set<RandomVariable*> tmp;
    set_intersection((*v)->neighbors.begin(),
		     (*v)->neighbors.end(),
		     C2.begin(),C2.end(),
		     inserter(tmp,tmp.end()));
    res.clear();    
    set_difference(tmp.begin(),tmp.end(),
		   next_right_C_r.begin(),next_right_C_r.end(),
		   inserter(res,res.end()));

    set_union(res.begin(),res.end(),
	      C_r.begin(),C_r.end(),
	      inserter(next_C_r,next_C_r.end()));
    next_C_r.erase((*v));

    if (setset.find(next_C_r) != setset.end())
      continue; // check if memoized, and if so, no need to go further.
    
    // check size of candiate interface, and keep a copy of
    // it if it is less then the one we have seen so far.
    if (next_C_r.size() < best_C_r.size()) {
      best_right_C_r = next_right_C_r;
      best_C_r = next_C_r;
    } 

    // memoize
    setset.insert(next_C_r);

    findBestLeftInterface(next_right_C_r,
			  next_C_r,
			  C2,C1,setset,
			  best_right_C_r,best_C_r);

  }
}


/*-
 *-----------------------------------------------------------------------
 * GMTemplate::triangulatePartitions()
 *   Create the three partitions and triangulate them.
 *
 * Preconditions:
 *
 * Postconditions:
 *
 * Side Effects:
 *
 * Results:
 *
 *
 *-----------------------------------------------------------------------
 */
void
GMTemplate::
triangulatePartitions()
{

  vector <RandomVariable*> unroll2_rvs;
  fp.unroll(2,unroll2_rvs);
  // drop all the edge directions
  for (unsigned i=0;i<unroll2_rvs.size();i++) {
    unroll2_rvs[i]->createNeighborsFromParentsChildren();
  }
  // moralize graph
  for (unsigned i=0;i<unroll2_rvs.size();i++) {
    unroll2_rvs[i]->moralize();    
  }
  // create sets P, C1, C2, C3, and E, from graph unrolled 2 times
  // prologue
  set<RandomVariable*> P_u2;
  // 1st chunk
  set<RandomVariable*> C1_u2;
  // 2nd chunk
  set<RandomVariable*> C2_u2;
  // 3rd chunk
  set<RandomVariable*> C3_u2;
  // epilogue
  set<RandomVariable*> E_u2;
  int start_index_of_C2_u2 = -1;
  for (unsigned i=0;i<unroll2_rvs.size();i++) {
    if (unroll2_rvs[i]->frame() < firstChunkFrame)
      P_u2.insert(unroll2_rvs[i]);
    else if (unroll2_rvs[i]->frame() <= lastChunkFrame)
      C1_u2.insert(unroll2_rvs[i]);
    else if (unroll2_rvs[i]->frame() <= lastChunkFrame+chunkNumFrames) {
      C2_u2.insert(unroll2_rvs[i]);
      if (start_index_of_C2_u2 == -1)
	start_index_of_C2_u2 = i;
    } else if (unroll2_rvs[i]->frame() <= lastChunkFrame+2*chunkNumFrames)
      C3_u2.insert(unroll2_rvs[i]);
    else 
      E_u2.insert(unroll2_rvs[i]);
  }

  assert (C1_u2.size() == C2_u2.size());
  assert (C2_u2.size() == C3_u2.size());
  printf("Size of (P,C1,C2,C3,E) = (%d,%d,%d,%d,%d)\n",
	 P_u2.size(),C1_u2.size(),C2_u2.size(),C3_u2.size(),E_u2.size());


  // this next routine gives us the best left interface that
  // exists from within the chunk C2_u2 and places
  // it in C_l_u2, and everything to the 'left' of C_l_u2
  // that still lies within C2_u2 is placed in left_C_l_u2
  set<RandomVariable*> left_C_l_u2C2;
  set<RandomVariable*> C_l_u2C2;
  findBestLeftInterface(C1_u2,C2_u2,C3_u2,left_C_l_u2C2,C_l_u2C2);


  // create sets P', C1', C2', and E', from graph unrolled 1 time
  vector <RandomVariable*> unroll1_rvs;
  fp.unroll(1,unroll1_rvs);
  for (unsigned i=0;i<unroll1_rvs.size();i++) {
    unroll1_rvs[i]->createNeighborsFromParentsChildren();
  }
  for (unsigned i=0;i<unroll1_rvs.size();i++) {
    unroll1_rvs[i]->moralize();    
  }
  set<RandomVariable*> P_u1;
  set<RandomVariable*> C1_u1;
  set<RandomVariable*> C2_u1;
  set<RandomVariable*> E_u1;
  int start_index_of_C1_u1 = -1;
  int start_index_of_C2_u1 = -1;
  for (unsigned i=0;i<unroll1_rvs.size();i++) {
    if (unroll1_rvs[i]->frame() < firstChunkFrame)
      P_u1.insert(unroll1_rvs[i]);
    else if (unroll1_rvs[i]->frame() <= lastChunkFrame) {
      C1_u1.insert(unroll1_rvs[i]);
      if (start_index_of_C1_u1 == -1)
	start_index_of_C1_u1 = i;
    } else if (unroll1_rvs[i]->frame() <= lastChunkFrame+chunkNumFrames) {
      C2_u1.insert(unroll1_rvs[i]);
      if (start_index_of_C2_u1 == -1)
	start_index_of_C2_u1 = i;
    } else 
      E_u1.insert(unroll1_rvs[i]);
  }

  assert (C1_u1.size() == C2_u1.size());
  assert (C1_u1.size() == C2_u2.size());

  // creat mapping from C2_u2 (for which we now
  // have an appropriate interface vars) to nodes C1_u1 and C2_u1
  // which we are going to use to set up the 3-way partition
  // to ultimately triangulate.
  map < RandomVariable*, RandomVariable* > C2_u2_to_C1_u1;
  map < RandomVariable*, RandomVariable* > C2_u2_to_C2_u1;
  for (unsigned i=0;i<C2_u2.size();i++) {
    C2_u2_to_C1_u1[unroll2_rvs[i+start_index_of_C2_u2]]
      = unroll1_rvs[i+start_index_of_C1_u1];
    C2_u2_to_C2_u1[unroll2_rvs[i+start_index_of_C2_u2]]
      = unroll1_rvs[i+start_index_of_C2_u1];

  }

  // now we need to make a bunch of sets to be unioned
  // together to get the partitions.
  set<RandomVariable*> C_l_u1C1;
  set<RandomVariable*> C_l_u1C2;
  for (set<RandomVariable*>::iterator i = C_l_u2C2.begin();
       i!= C_l_u2C2.end(); i++) {

    C_l_u1C1.insert(C2_u2_to_C1_u1[(*i)]);
    C_l_u1C2.insert(C2_u2_to_C2_u1[(*i)]);					  
  }

  // These (C_l_u1C1 and C_l_u1C2) are the interfaces which are forced
  // to be complete (i.e., part of a maxclique).
  // TODO: these next steps will ruin the clique_size = 2 property
  //     of the 'snake' structure. The todo is to get this working 
  //     with that (and similar) structures.
  /*
   *     idea: make complete components in C_l *only* if they
   *           are connected via nodes/edges within either (P + left_C_l)
   *           or preceeding C.
   *           What we will then have is a collection of cliques for 
   *           the interface(s). In this case, we glue together
   *           the corresponding sets of cliques. Right interface
   *           algorithm should be similar (and use E rather than P).
   *           But might both left and right interface need to be used in the
   *           same repeated chunk in this case to get clique_size=2???
   *          
   */
  makeComplete(C_l_u1C1);
  makeComplete(C_l_u1C2);

  set<RandomVariable*> left_C_l_u1C1;
  set<RandomVariable*> left_C_l_u1C2;
  for (set<RandomVariable*>::iterator i = left_C_l_u2C2.begin();
       i != left_C_l_u2C2.end(); i++) {

    left_C_l_u1C1.insert(C2_u2_to_C1_u1[(*i)]);
    left_C_l_u1C2.insert(C2_u2_to_C2_u1[(*i)]);					  
  }
  
  // finally, create the sets A, B, and C
  // which are to be triangulated separately.
  //  A = P' + C1'(left_C_l) + C1'(C_l)
  //  B = C1'\C1'(left_C_l) + C2'(left_C_l) + C2'(C_l)
  //  C = C2'\C2'(left_C_l) + E

  set<RandomVariable*> A = P_u1;
  set<RandomVariable*> B;
  set<RandomVariable*> C = E_u1;

  // Finish A
  set_union(left_C_l_u1C1.begin(),left_C_l_u1C1.end(),
	    C_l_u1C1.begin(),C_l_u1C1.end(),
	    inserter(A,A.end()));

  // B
  set_union(left_C_l_u1C2.begin(),left_C_l_u1C2.end(),
	    C_l_u1C2.begin(),C_l_u1C2.end(),
	    set_difference(C1_u1.begin(),C1_u1.end(),
			   left_C_l_u1C1.begin(),left_C_l_u1C1.end(),
			   inserter(B,B.end())));
  // finish C
  set_difference(C2_u1.begin(),C2_u1.end(),
		 left_C_l_u1C2.begin(),left_C_l_u1C2.end(),
		 inserter(C,C.end()));
		 

  printf("---\nSet A is:\n");
  for (set<RandomVariable*>::iterator i=A.begin();
       i != A.end(); i++) {
    RandomVariable* rv = (*i);
    printf("%s(%d)\n",rv->name().c_str(),rv->frame());
  }

  printf("---\nSet B is:\n");
  for (set<RandomVariable*>::iterator i=B.begin();
       i != B.end(); i++) {
    RandomVariable* rv = (*i);
    printf("%s(%d)\n",rv->name().c_str(),rv->frame());
  }


  printf("---\nSet C is:\n");
  for (set<RandomVariable*>::iterator i=C.begin();
       i != C.end(); i++) {
    RandomVariable* rv = (*i);
    printf("%s(%d)\n",rv->name().c_str(),rv->frame());
  }

  set<RandomVariable*> Ac;
  set<RandomVariable*> Bc;
  set<RandomVariable*> Cc;

  clone(A,Ac);
  clone(B,Bc);
  clone(C,Cc);

  set<MaxClique> Acliques;
  set<MaxClique> Bcliques;
  set<MaxClique> Ccliques;

  triangulate(Ac,Acliques,MIN_WEIGHT);
  triangulate(Bc,Bcliques,MIN_WEIGHT);
  triangulate(Cc,Ccliques,MIN_WEIGHT);

}



/*-
 *-----------------------------------------------------------------------
 * GMTemplate::clone()
 *   Clone a set of random variables from 'in' to 'out'. The cloned
 *   variables have the same set of parents and children
 *   but the 'neighbors' member is entirely new and points
 *   only to the corresponding neighbors in the new 'out' set.
 *   Therefore, parents,children should not be used
 *   among cloned set.
 *   TODO: this needs to be changed, and completed.  
 *
 * Preconditions:
 *   'in' is a set of random variables to be cloned.
 *
 * Postconditions:
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
GMTemplate::
clone(set<RandomVariable*>& in, 
      set<RandomVariable*>& out
      )
{
  
  map < RandomVariable*, RandomVariable* > in_to_out;
  for (set<RandomVariable*>::iterator i=in.begin();
       i != in.end(); i++) {
    RandomVariable*rv = (*i)->cloneWithoutParents();
    out.insert(rv);
    in_to_out[(*i)] = rv;
  }

  // just do neighbors for now, don't bother with parents, children,
  // and so on.
  // Set the neighbors of out to be the correctly associated
  // variables, but do not include neighbors that are
  // not in the current set (i.e., dissociate with any
  // other possible portion of the network)

  for (set<RandomVariable*>::iterator i=in.begin();
       i != in.end(); i++) {
    RandomVariable*rv = in_to_out[(*i)];
    set<RandomVariable*> tmp;
    for (set<RandomVariable*>::iterator j = rv->neighbors.begin();
	 j != rv->neighbors.end(); j++) {    
      if (in.find((*j)) != in.end()) {
	// then it is included in this set.
	tmp.insert(in_to_out[(*j)]);
      }
    }
    rv->neighbors = tmp;
  }
}



/*-
 *-----------------------------------------------------------------------
 * GMTemplate::triangulate()
 *   Triangulate a set of nodes using one of three very
 *   simple elimination heuristics (min weight, min size, or min fill).
 *   For a good description of these heuristics, see
 *      D. Rose, 1970.
 *    
 * Preconditions:
 *
 * Postconditions:
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
GMTemplate::
triangulate(
	    // input nodes to triangulate
	    set<RandomVariable*> nodes,
	    // nodes ordered according to resulting triangulation
	    set<RandomVariable*>& orderedNodes,  
	    // resulting max cliques
	    set<MaxClique>& cliques,             
	    // triangulation heuristic
	    const TriangulateHeuristic th        
	    )
{
  
  // Create a priority queue data structure using a multimap, so that
  // those nodes with the lowest 'X' where 'X' is the heuristic
  // can be accessed immediately in the front, and where
  // whenever a node gets eliminated, only its neighbors weight
  // can be recalculated easily.

  set<RandomVariable*> unorderedNodes = nodes;
  orderedNodes.clear();

  multimap<double,RandomVariable*> unorderedNodes;
  // first go through and create initial data structure

  for (set<RandomVariable*>::iterator i = nodes.begin();
       i != nodes.end();
       i++) {
    double weight;
    
    // different Heuristics
    if (th == MIN_SIZE) {
      weight = (*i)->neighbors.size();
    } else if (th == MIN_WEIGHT) {
      weight = 1;
      for (set<RandomVariable*>::iterator j=(*i)->neighbors.begin();
	   j != (*i)->neighbors.end();
	   j++) {
	RandomVariable *const rv = (*j);
	if (rv->discrete && rv->hidden) {
	  DiscreteRandomVariable *const drv = rv;
	  if (!drv->deterministic()) {
	    weight *= drv->cardinality;
	    // TODO: change to line below when
	    // reading in parameters.
	    // weight *= drv->averageCardinality();
	  }
	}
      }
    } else { // (th == MIN_FILLIN)


    }




  }

  while (unorderedNodes.size() > 0) {
    // go through and update multi-map
    

  }





}








////////////////////////////////////////////////////////////////////
//        Test Driver
////////////////////////////////////////////////////////////////////

#ifdef MAIN


#endif
