/*
 * GMTK_BoundaryTriangulate.cc
 *   The GMTK Boundary Algorithm and Graph Triangulation Support Routines
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu> & Chris Bartels <bartels@ee.washington.edu>
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



#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern long int lrand48(void);

#include <algorithm>
#include <iterator>
#include <list>
#include <map>
#include <set>

#include "debug.h"
#include "error.h"
#include "general.h"
#include "rand.h"

#include "GMTK_BoundaryTriangulate.h"
#include "GMTK_ContinuousRandomVariable.h"
#include "GMTK_DiscreteRandomVariable.h"
#include "GMTK_FileParser.h"
#include "GMTK_MDCPT.h"
#include "GMTK_Mixture.h"
#include "GMTK_GM.h"
#include "GMTK_GMParms.h"
#include "GMTK_MSCPT.h"
#include "GMTK_MTCPT.h"
#include "GMTK_ObservationMatrix.h"
#include "GMTK_RandomVariable.h"
#include "GMTK_JunctionTree.h"
#include "GMTK_GraphicalModel.h"

VCID("$Header$");

#ifndef MAXFLOAT
#define MAXFLOAT (3.40282347e+38F)
#endif


/*******************************************************************************
 *******************************************************************************
 *******************************************************************************
 **
 **           General Support Routines
 **
 *******************************************************************************
 *******************************************************************************
 *******************************************************************************
 */

// TODO: move this somewhere else, generally accessible.
static void
printRVSet(FILE*f,set<RandomVariable*>& locset)
{
  bool first = true;
  set<RandomVariable*>::iterator it;
  for (it = locset.begin();
       it != locset.end();it++) {
    RandomVariable* rv = (*it);
    if (!first)
      fprintf(f,",");
    fprintf(f,"%s(%d)",rv->name().c_str(),rv->frame());
    first = false;
  }
  fprintf(f,"\n");
}
static void
printRVSet(FILE*f,vector<RandomVariable*>& locvec)
{
  bool first = true;
  for (unsigned i=0;i<locvec.size();i++) {
    RandomVariable* rv = locvec[i];
    if (!first)
      fprintf(f,",");
    fprintf(f,"%s(%d)",rv->name().c_str(),rv->frame());
    first = false;
  }
  fprintf(f,"\n");
}



/*-
 *-----------------------------------------------------------------------
 * Partition::saveCurrentNeighbors()
 *   Save a copy of the current neighbors of the partition.
 *
 * Preconditions:
 *   The corresponding partition must be instantiated.
 *
 * Postconditions:
 *   a copy of the current neighborset will be made
 *
 * Side Effects:
 *   Any previous neighbor set will be lost.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
BoundaryTriangulate::
saveCurrentNeighbors(const set<RandomVariable*> nodes,vector<nghbrPairType>& orgnl_nghbrs)
{
  set<RandomVariable*>::iterator crrnt_node; 
  set<RandomVariable*>::iterator end_node; 
  orgnl_nghbrs.clear();
  for ( crrnt_node=nodes.begin(), end_node=nodes.end();
        crrnt_node != end_node;
        ++crrnt_node ) {
    orgnl_nghbrs.push_back(
      nghbrPairType((*crrnt_node), (*crrnt_node)->neighbors) );
  }
}


/*-
 *-----------------------------------------------------------------------
 * Partition::saveCurrentNeighbors()
 *   Save a copy of the current neighbors of the partition.
 *
 * Preconditions:
 *   The corresponding partition must be instantiated.
 *
 * Postconditions:
 *   a copy of the current neighborset will be made
 *
 * Side Effects:
 *   Any previous neighbor set will be lost.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
BoundaryTriangulate::
saveCurrentNeighbors(
  vector<triangulateNode>& nodes,
  vector<triangulateNghbrPairType>& orgnl_nghbrs
  )
{
  vector<triangulateNode>::iterator crrnt_node; 
  vector<triangulateNode>::iterator end_node; 

  orgnl_nghbrs.clear();
  for ( crrnt_node=nodes.begin(), end_node=nodes.end();
        crrnt_node != end_node;
        ++crrnt_node ) {
    orgnl_nghbrs.push_back(
      triangulateNghbrPairType(&(*crrnt_node), (*crrnt_node).neighbors) );
  }
}


/*-
 *-----------------------------------------------------------------------
 * Partition::restoreNeighbors()
 *   Restore graph's neighbors to state on previous save.
 *   This can be used as an "unTriangulate" routine.
 *
 * Preconditions:
 *   The corresponding partition must be instantiated. 
 *   For this routine to have any effect, previous neighbors must have 
 *     already been saved.
 *
 * Postconditions:
 *   Old neighbors restored.
 *
 * Side Effects:
 *   Any previous neighbor set will be lost.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
BoundaryTriangulate::
restoreNeighbors(vector<nghbrPairType>& orgnl_nghbrs)
{
  vector<pair<RandomVariable*, set<RandomVariable*> > >::iterator crrnt_node;
  vector<pair<RandomVariable*, set<RandomVariable*> > >::iterator end_node;
  for( end_node=orgnl_nghbrs.end(), crrnt_node=orgnl_nghbrs.begin();
       crrnt_node!=end_node;
       ++crrnt_node )
  {
    (*crrnt_node).first->neighbors = (*crrnt_node).second;
  }
}


/*-
 *-----------------------------------------------------------------------
 * Partition::restoreNeighbors()
 *   Restore graph's neighbors to state on previous save.
 *   This can be used as an "unTriangulate" routine.
 *
 * Preconditions:
 *   The corresponding partition must be instantiated. 
 *   For this routine to have any effect, previous neighbors must have 
 *     already been saved.
 *
 * Postconditions:
 *   Old neighbors restored.
 *
 * Side Effects:
 *   Any previous neighbor set will be lost.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
BoundaryTriangulate::
restoreNeighbors(vector<triangulateNghbrPairType>& orgnl_nghbrs)
{
  vector<pair<triangulateNode*, vector<triangulateNode*> > >::iterator 
    crrnt_node;
  vector<pair<triangulateNode*, vector<triangulateNode*> > >::iterator  
    end_node;

  for( end_node=orgnl_nghbrs.end(), crrnt_node=orgnl_nghbrs.begin();
       crrnt_node!=end_node;
       ++crrnt_node )
  {
    (*crrnt_node).first->neighbors = (*crrnt_node).second;
  }
}


/*-
 *-----------------------------------------------------------------------
 * BoundaryTriangulate::deleteNodes()
 *   Given a set of random variables, delete the nodes pointed to by
 *   the set.
 *
 * Preconditions:
 *   'nodes' is a valid set of node pointers
 *
 * Postconditions:
 *   all RV*'s in 'nodes' have been deleted. The set should thereafter
 *   immediately be deleted or filled with new nodes
 *
 * Side Effects:
 *     Will delete all variables pointed to by the RV*'s within the set.
 *
 * Results:
 *     none
 *
 *
 *-----------------------------------------------------------------------
 */
void
BoundaryTriangulate::
deleteNodes(const set<RandomVariable*>& nodes)
{
  for (set<RandomVariable*>::iterator i = nodes.begin();
       i != nodes.end(); i++) 
    delete (*i);
}



/*-
 *-----------------------------------------------------------------------
 * BoundaryTriangulate::parseTriHeuristicString()
 *      parses a triheuristic string returning a triHeuristic structure
 *
 * Preconditions:
 *      - String should contain valid tri heuristic, see code
 *        for definition
 *
 * Postconditions:
 *      - structure has been filled in.
 *
 * Side Effects:
 *      none
 *
 * Results:
 *     filled in vector via argument
 *
 *-----------------------------------------------------------------------
 */
void
BoundaryTriangulate::parseTriHeuristicString(const string& tri_heur_str,
					     TriangulateHeuristics& tri_heur)
{
  tri_heur.init();
  if (tri_heur_str.size() == 0) {
    return;
  } else {
    // first see if tri_heur_str starts with '<num>-'
    // to get number of trials.
    char *endp;
    long numTrials = -1;
    const char *startp = tri_heur_str.c_str();
    numTrials = strtol(startp,&endp,10);
    if (endp != startp)
      // then there is an iteration number
      tri_heur.numberTrials = numTrials;
    else
      tri_heur.numberTrials = 1;
    if (*endp == '-')
      endp++;
    if (!*endp || tri_heur.numberTrials <= 0)
      error("ERROR: bad triangulation heuristic string given in '%s'\n",startp);

    if (strncmp(endp, "anneal", strlen(endp)) == 0) {
      tri_heur.style = TS_ANNEALING;
    } else if (strncmp(endp, "exhaustive", strlen(endp)) == 0) {
      tri_heur.style = TS_EXHAUSTIVE;
    } else if (strncmp(endp, "pre-edge-all", strlen(endp)) == 0) {
      tri_heur.style = TS_PRE_EDGE_ALL;
    } else if (strncmp(endp, "pre-edge-lo", strlen(endp)) == 0) {
      tri_heur.style = TS_PRE_EDGE_LO;
    } else if (strncmp(endp, "pre-edge-random", strlen(endp)) == 0) {
      tri_heur.style = TS_PRE_EDGE_RANDOM;
    } else if (strncmp(endp, "heuristics", strlen(endp)) == 0) {
      tri_heur.style = TS_ELIMINATION_HEURISTICS;
    } else if (strncmp(endp, "frontier", strlen(endp)) == 0) {
      tri_heur.style = TS_FRONTIER;
    } else if (strcmp(endp, "MCS") == 0) { 
      // MCS must be specified with full "MCS" string, so use strcmp
      // rather than strncmp here.
      tri_heur.style = TS_MCS;
    } else if (strncmp(endp, "completed", strlen(endp)) == 0) {
      tri_heur.style = TS_COMPLETED;
    } else {
      tri_heur.style = TS_BASIC;
      tri_heur.heuristic_vector.clear();
      tri_heur.basic_method_string = string(endp);
      while (*endp) {
	switch (*endp) {
	case 'S':
	  tri_heur.heuristic_vector.push_back(TH_MIN_SIZE);
	  break;
	case 'T':
	  tri_heur.heuristic_vector.push_back(TH_MIN_TIMEFRAME);
	  break;
	case 'X':
	  tri_heur.heuristic_vector.push_back(TH_MAX_TIMEFRAME);
	  break;
	case 'F':
	  tri_heur.heuristic_vector.push_back(TH_MIN_FILLIN);
	  break;
	case 'W':
	  tri_heur.heuristic_vector.push_back(TH_MIN_WEIGHT);
	  break;
	case 'E':
	  tri_heur.heuristic_vector.push_back(TH_MIN_ENTROPY);
	  break;
	case 'P':
	  tri_heur.heuristic_vector.push_back(TH_MIN_POSITION_IN_FILE);
	  break;
	case 'H':
	  tri_heur.heuristic_vector.push_back(TH_MIN_HINT);
	  break;
	case 'N':
	  tri_heur.heuristic_vector.push_back(TH_MIN_WEIGHT_NO_D);
	  break;
	case 'R':
	  tri_heur.heuristic_vector.push_back(TH_RANDOM);
	  break;
	default:
	  error("ERROR: Unknown triangulation heuristic given '%c' in string '%s'\n",
		*endp,tri_heur_str.c_str());
	  break;
	}
	endp++;
      }
    }
  }
  return;
}


/*-
 *-----------------------------------------------------------------------
 * BoundaryTriangulate::createVectorBoundaryHeuristic()
 *      create a vector of prioritized boundary heuristics based
 *      on a string that is passed in.
 *
 * Preconditions:
 *      - String should contain set of valid heuristcs, see code
 *        for definition
 *
 * Postconditions:
 *      - vector has been filled in.
 *
 * Side Effects:
 *      none
 *
 * Results:
 *     filled in vector via argument
 *
 *-----------------------------------------------------------------------
 */
void
BoundaryTriangulate::createVectorBoundaryHeuristic(const string& bnd_heur_str,
					  vector<BoundaryHeuristic>& bnd_heur_v)
{
  if (bnd_heur_str.size() == 0) {
    // default case.
    // first by weight
    bnd_heur_v.push_back(IH_MIN_WEIGHT); 
    // then by fill in if weight in tie
    bnd_heur_v.push_back(IH_MIN_FILLIN); 
  } else {
    for (unsigned i=0;i<bnd_heur_str.size();i++) {
      switch (bnd_heur_str[i]) {
      case 'S':
	bnd_heur_v.push_back(IH_MIN_SIZE);
	break;
      case 'F':
	bnd_heur_v.push_back(IH_MIN_FILLIN);
	break;
      case 'W':
	bnd_heur_v.push_back(IH_MIN_WEIGHT);
	break;
      case 'E':
	bnd_heur_v.push_back(IH_MIN_ENTROPY);
	break;
      case 'C':
	bnd_heur_v.push_back(IH_MIN_MAX_C_CLIQUE);
	break;
      case 'M':
	bnd_heur_v.push_back(IH_MIN_MAX_CLIQUE);
	break;
      case 'A':
	bnd_heur_v.push_back(IH_MIN_STATE_SPACE);
	break;
      case 'Q':
	bnd_heur_v.push_back(IH_MIN_C_STATE_SPACE);
	break;
      case 'N':
	bnd_heur_v.push_back(IH_MIN_WEIGHT_NO_D);
	break;

      case 'p':
	bnd_heur_v.push_back(IH_MIN_MIN_POSITION_IN_FILE);
	break;
      case 'P':
	bnd_heur_v.push_back(IH_MIN_MAX_POSITION_IN_FILE);
	break;

      case 't':
	bnd_heur_v.push_back(IH_MIN_MIN_TIMEFRAME);
	break;
      case 'T':
	bnd_heur_v.push_back(IH_MIN_MAX_TIMEFRAME);
	break;

      case 'h':
	bnd_heur_v.push_back(IH_MIN_MIN_HINT);
	break;
      case 'H':
	bnd_heur_v.push_back(IH_MIN_MAX_HINT);
	break;

      default:
	error("ERROR: Unknown triangulation heuristic given '%c' in string '%s'\n",
	      bnd_heur_str[i],bnd_heur_str.c_str());
	break;
      }
    }
  }
  return;
}



/*-
 *-----------------------------------------------------------------------
 * BoundaryTriangulate::computeFillIn()
 *   Computes the number of edges that would need to be added
 *   among 'nodes' to make 'nodes' complete.
 *
 * Preconditions:
 *   Set of nodes must be valid meaning that it has valid neighbors,
 *   parents, and children member variables.
 *
 * Postconditions:
 *   computed weight is provided.
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
int
BoundaryTriangulate::
computeFillIn(const set<RandomVariable*>& nodes) 
{

  int fill_in = 0;
  for (set<RandomVariable*>::iterator j=nodes.begin();
       j != nodes.end();
       j++) {

    // TODO: figure out if there is a way to just to compute
    // the size of the set intersection rather than to
    // actually produce the set intersection and then use its size.
    set<RandomVariable*> tmp;
    set_intersection(nodes.begin(),nodes.end(),
		     (*j)->neighbors.begin(),(*j)->neighbors.end(),
		     inserter(tmp,tmp.end()));

    // Nodes i and j should share the same neighbors except for
    // node i has j as a neighbor and node j has i as a
    // neighbor.  Therefore, if fill in is zero, then 
    //    tmp.size() == (*i)->neighbors.size() - 1
    // (we subtract 1 since we don't include the neighbor
    //  in the fill in, it is included in neighbors 
    //  but not tmp).
    // In otherwords, we have:

    fill_in += (nodes.size() - 1 - tmp.size());
  }
  // counted each edge twice, so fix that (although not 
  // strictly necessary since we could just compute with 2*fill_in,
  // we do this, however, since the user will be happier. 
  fill_in /= 2;
  return fill_in;

}



/*-
 *-----------------------------------------------------------------------
 * BoundaryTriangulate::graphWeight()
 *   Calculate the weight of a graph from a vector of cliques.  
 *
 * Preconditions:
 *   none
 * 
 * Postconditions:
 *   none
 *
 * Side Effects:
 *   none
 *
 * Results:
 *   The log base 10 weight of the graph
 *-----------------------------------------------------------------------
 */
double 
BoundaryTriangulate::
graphWeight(vector<MaxClique>& cliques)
{
  vector<MaxClique>::iterator crrnt_clique; 
  vector<MaxClique>::iterator end_clique; 
  double crrnt_weight;
  double weight = -1;

  for( end_clique=cliques.end(), crrnt_clique=cliques.begin(); 
       crrnt_clique != end_clique;
       ++crrnt_clique )
  {
    crrnt_weight = MaxClique::computeWeight(crrnt_clique->nodes);

    if (weight < 0) {
      weight = crrnt_weight; 
    }
    else {
      weight = weight + log10(1+pow(10,crrnt_weight-weight));
    } 
  } 

  return(weight);
}



/*-
 *-----------------------------------------------------------------------
 * BoundaryTriangulate::graphWeight()
 *   Calculate the weight of a graph from a vector of cliques.  This
 *   version also allows for an estimate of the JT weight (i.e., how
 *   much this set of cliques would cost when set up in a junction
 *   tree and used during the collect evidence stage in inference).
 *
 * Preconditions:
 *   none
 * 
 * Postconditions:
 *   none
 *
 * Side Effects:
 *   Changes some of the JT members of maxcliques.
 *
 * Results:
 *   The log base 10 weight of the graph
 *-----------------------------------------------------------------------
 */
double 
BoundaryTriangulate::
graphWeight(vector<MaxClique>& cliques,		     
	    const bool useJTWeight,
	    const set<RandomVariable*>& interfaceNodes)
{
  if (!useJTWeight)
    return graphWeight(cliques);
  else
    return JunctionTree::junctionTreeWeight(cliques,interfaceNodes);
}



/*******************************************************************************
 *******************************************************************************
 *******************************************************************************
 **
 **         Main Partition Routines
 **
 *******************************************************************************
 *******************************************************************************
 *******************************************************************************
 */


/*-
 *-----------------------------------------------------------------------
 * BoundaryTriangulate::findPartitions()
 *  Create the three partitions (P,C,E) of the template using
 *  the heuristics supplied.
 *     fh = a string with triangulation heuristics to use (in order)
 *     flr = force the use of either left or right, rather than use min.
 *     findBestBoundary = T/F if to use the exponential boundary finding alg.
 *
 * Preconditions:
 *   Object must be instantiated and have the use of the information
 *   in a valid FileParser object (which stores the parsed structure file)
 *   Arguments must indicate valid heuristics to use.
 *
 * Postconditions:
 *   Arguments Pc, Cc, and Ec now contain partitions in a separate
 *   graph from the FileParser (i.e., file parser information is
 *   not disturbed)
 *
 * Side Effects:
 *   none
 *
 * Results:
 *   Pc, Cc, and Ec
 *
 *
 *-----------------------------------------------------------------------
 */
void
BoundaryTriangulate::
findPartitions(// boundary quality heuristic
	       const string& bnd_heur_str,  
	       // force use of only left or right interface
	       const string& flr, 
	       // triangualtion heuristic to use
	       const string& tri_heur_str, 
	       // should we run the exponential find best interface
	       const bool findBestBoundary,
	       // the resulting template
	       GMTemplate& gm_template)
{

  //
  // M = number of chunks in which boundary algorithm is allowed to exist
  // i.e., a boundary is searched for within M repeated chunks, where M >= 1.
  // Note that M and S put constraints on the number of time frames
  // of the observations (i.e., the utterance length in time frames).
  // Specificaly, if N = number of frames of utterance, 
  //   then we must have N = length(P) + length(E) + (M+k*S)*length(C)
  //   where k = 1,2,3,4, ... is some integer >= 1.
  // Therefore, making M and/or S larger reduces the number valid possible utterance lengths.

  // First create a network (called u2, but could be called "find
  // boundary") that is used to find the boundary in. This network is
  // the template unrolled M+1 times. M = Max number simultaneous
  // chunks in which interface boundary may exist. By unrolling M+1
  // times, we create a network with M+2 chunks. The first and last
  // one are used to ensure boundary algorithm works (since we are not
  // guaranteed that P or E will exist), and the middle M chunks are
  // used to find the boundary.

  vector <RandomVariable*> unroll2_rvs;
  fp.unroll(M+1,unroll2_rvs);
  // drop all the edge directions
  for (unsigned i=0;i<unroll2_rvs.size();i++) {
    unroll2_rvs[i]->createNeighborsFromParentsChildren();
  }
  // moralize graph
  for (unsigned i=0;i<unroll2_rvs.size();i++) {
    unroll2_rvs[i]->moralize();
  }
  // create sets P, C1, C2, C3, and E, from graph unrolled M+1 times
  // prologue
  set<RandomVariable*> P_u2;
  // 1st chunk, 1 chunk long
  set<RandomVariable*> C1_u2;
  // 2nd chunk, M chunks long
  set<RandomVariable*> C2_u2;
  // 3rd chunk, 1 chunk long
  set<RandomVariable*> C3_u2;
  // epilogue
  set<RandomVariable*> E_u2;
  int start_index_of_C1_u2 = -1;
  int start_index_of_C2_u2 = -1;
  int start_index_of_C3_u2 = -1;
  int start_index_of_E_u2 = -1;
  for (unsigned i=0;i<unroll2_rvs.size();i++) {
    if (unroll2_rvs[i]->frame() < fp.firstChunkFrame())
      // prologue
      P_u2.insert(unroll2_rvs[i]);
    else if (unroll2_rvs[i]->frame() <= fp.lastChunkFrame()) {
      // 1st chunk, 1 chunk long
      C1_u2.insert(unroll2_rvs[i]);
      if (start_index_of_C1_u2 == -1)
	start_index_of_C1_u2 = i;
    } else if (unroll2_rvs[i]->frame() <= fp.lastChunkFrame()+M*fp.numFramesInC()) {
      // 2nd chunk, M chunks long
      C2_u2.insert(unroll2_rvs[i]);
      if (start_index_of_C2_u2 == -1)
	start_index_of_C2_u2 = i;
    } else if (unroll2_rvs[i]->frame() <= fp.lastChunkFrame()+(M+1)*fp.numFramesInC()) {
      // 3rd chunk, 1 chunk long
      C3_u2.insert(unroll2_rvs[i]);
      if (start_index_of_C3_u2 == -1)
	start_index_of_C3_u2 = i;
    } else {
      // epilogue
      E_u2.insert(unroll2_rvs[i]);
      if (start_index_of_E_u2 == -1)
	start_index_of_E_u2 = i;
    }
  }

  assert (M*C1_u2.size() == C2_u2.size());
  assert (C2_u2.size() == M*C3_u2.size());
  infoMsg(High,"Size of (P,C1,C2,C3,E) = (%d,%d,%d,%d,%d)\n",
	  P_u2.size(),C1_u2.size(),C2_u2.size(),C3_u2.size(),E_u2.size());


  ///////////////////////////////////////////////////////////////////
  // Next, create the network from which the new P, C, and E will
  // be formed. I.e., this new P, C, and E will be the graph partitions
  // that are ultimately triangulated, where the new C will
  // be unrolled as appropriate to get a full network. We
  // call the network from which P, C, and E are created u1 (but
  // could call it "partition")

  // When we have an S and an M parameter, , number of chunks needed
  // is M + S, so we should unroll (M+S-1) times.  We start by
  // creating sets P', C1', C2', and E', from the graph unrolled
  // (S+M-1) time(s), where each hyper-chunk C1' and C2' are M frames long
  // Note that 
  //    1) if S==M, Cextra will be empty, and C1' and C2' will not overlap.
  //    2) if M > S, then C1' and C2' will overlap certain variables.
  //    3) if S > M, C1' and C2' will not overlap, and there will be chunks 
  //       between C1' and C2', and they'll be placed in Cextra

  vector <RandomVariable*> unroll1_rvs;
  fp.unroll(M+S-1,unroll1_rvs);
  for (unsigned i=0;i<unroll1_rvs.size();i++) {
    unroll1_rvs[i]->createNeighborsFromParentsChildren();
  }
  for (unsigned i=0;i<unroll1_rvs.size();i++) {
    unroll1_rvs[i]->moralize();
  }
  set<RandomVariable*> P_u1;
  set<RandomVariable*> C1_u1;
  set<RandomVariable*> C2_u1;
  set<RandomVariable*> Cextra_u1;
  set<RandomVariable*> E_u1;
  int start_index_of_C1_u1 = -1;
  int start_index_of_C2_u1 = -1;
  int start_index_of_E_u1 = -1;
  for (unsigned i=0;i<unroll1_rvs.size();i++) {
    // Note that there are some casts to (int) in the below below
    // since it might be the case that M = 0, and if unsigned
    // comparisions are used, the condition could fail inappropriately
    if (unroll1_rvs[i]->frame() < fp.firstChunkFrame())
      P_u1.insert(unroll1_rvs[i]);
    if ((unroll1_rvs[i]->frame() >= fp.firstChunkFrame()) &&
	((int)unroll1_rvs[i]->frame() <= (int)fp.lastChunkFrame() + ((int)M-1)*(int)fp.numFramesInC())) {
      C1_u1.insert(unroll1_rvs[i]);
      if (start_index_of_C1_u1 == -1)
	start_index_of_C1_u1 = i;
    }
    if (((int)unroll1_rvs[i]->frame() > (int)fp.lastChunkFrame() + ((int)M-1)*(int)fp.numFramesInC()) &&
	(unroll1_rvs[i]->frame() < fp.firstChunkFrame() + S*fp.numFramesInC())) {
      // this should only happen when S > M
      assert ( S > M );
      Cextra_u1.insert(unroll1_rvs[i]);
    }
    if ((unroll1_rvs[i]->frame() >= fp.firstChunkFrame() + S*fp.numFramesInC()) &&
	((int)unroll1_rvs[i]->frame() <= (int)fp.lastChunkFrame() + ((int)M+(int)S-1)*(int)fp.numFramesInC())) {
      C2_u1.insert(unroll1_rvs[i]);
      if (start_index_of_C2_u1 == -1)
	start_index_of_C2_u1 = i;
    }
    if ((int)unroll1_rvs[i]->frame() > (int)fp.lastChunkFrame() + ((int)M+(int)S-1)*(int)fp.numFramesInC()) {
      E_u1.insert(unroll1_rvs[i]);
      if (start_index_of_E_u1 == -1)
	start_index_of_E_u1 = i;
    }
  }

  // sanity checks
  assert (C1_u1.size() == C2_u1.size());
  assert (C1_u1.size() == C2_u2.size());
  assert ((S<=M) || (Cextra_u1.size() == C1_u2.size()*(S-M)));
  infoMsg(High,"Size of (P,C1,C2,E) = (%d,%d,%d,%d)\n",
	  P_u1.size(),C1_u1.size(),C2_u1.size(),E_u1.size());




  // create mapping from C2_u2 (for which we now
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

  // allocate space for results
  // left of the left interface in u2C2
  set<RandomVariable*> left_C_l_u2C2;
  // the left interface in u2C2
  set<RandomVariable*> C_l_u2C2;
  // right of the right interface in u2C2
  set<RandomVariable*> right_C_r_u2C2;
  // the right interface in u2C2
  set<RandomVariable*> C_r_u2C2;

  vector<BoundaryHeuristic> bnd_heur_v;
  createVectorBoundaryHeuristic(bnd_heur_str,bnd_heur_v);

  TriangulateHeuristics tri_heur;
  parseTriHeuristicString(tri_heur_str,tri_heur);

  vector<float> best_L_score;
  vector<float> best_R_score;

  if (flr.size() == 0 || (flr.size() > 0 && toupper(flr[0]) != 'R')) {
    // find best left interface
    findingLeftInterface = true;
    infoMsg(Tiny,"---\nFinding BEST LEFT interface\n");
    

    set<RandomVariable*> C2_1_u2; // first chunk in C2_u2
    if (M == 0) {
      // TODO: remove this, as we require M >= 1
      // find interface in basic template
      assert ( C2_u2.size() == 0 );
      findBestInterface(C1_u2,C3_u2,C2_1_u2,C3_u2,
			left_C_l_u2C2,C_l_u2C2,best_L_score,
			bnd_heur_v,
			false,
			// find best boundary args
			tri_heur,
			P_u1,
			C1_u1,
			Cextra_u1,
			C2_u1,
			E_u1,
			C2_u2_to_C1_u1,
			C2_u2_to_C2_u1
			);
    } else { 
      if (M == 1) {
	// make it empty signaling that we don't bother to check it in this case.
	C2_1_u2.clear();
      } else {
	for (set<RandomVariable*>::iterator i=C2_u2.begin();
	     i != C2_u2.end();i++) {
	  if ((*i)->frame() > fp.lastChunkFrame() && (*i)->frame() <= fp.lastChunkFrame()+fp.numFramesInC())
	    C2_1_u2.insert((*i));
	}
      }
      findBestInterface(C1_u2,C2_u2,C2_1_u2,C3_u2,
			left_C_l_u2C2,C_l_u2C2,best_L_score,
			bnd_heur_v,
			findBestBoundary,
			// find best boundary args
			tri_heur,
			P_u1,
			C1_u1,
			Cextra_u1,
			C2_u1,
			E_u1,
			C2_u2_to_C1_u1,
			C2_u2_to_C2_u1
			);
    }
  }
  if (flr.size() == 0 || (flr.size() > 0 && toupper(flr[0]) != 'L')) {
    // find best right interface
    findingLeftInterface = false;
    infoMsg(Tiny,"---\nFinding BEST RIGHT interface\n");

    set<RandomVariable*> C2_l_u2; // last chunk of C2_u2
    if (M == 0) {
      // TODO: remove this, as we require M >= 1
      // find interface in basic template
      assert ( C2_u2.size() == 0 );
      findBestInterface(C3_u2,C1_u2,C2_l_u2,C1_u2,
			right_C_r_u2C2,C_r_u2C2,best_R_score,
			bnd_heur_v,
			false,
			// find best boundary args
			tri_heur,
			E_u1,
			C2_u1,
			Cextra_u1,
			C1_u1,
			P_u1,
			C2_u2_to_C2_u1,
			C2_u2_to_C1_u1
			);
    } else {
      if (M == 1) {
	// make it empty signaling that we don't bother to check it in this case.
	C2_l_u2.clear();
      } else {
	for (set<RandomVariable*>::iterator i=C2_u2.begin();
	     i != C2_u2.end();i++) {
	  if ((*i)->frame() > (fp.lastChunkFrame()+(M-1)*fp.numFramesInC()) 
	      && (*i)->frame() <= (fp.lastChunkFrame()+M*fp.numFramesInC()))
	    C2_l_u2.insert((*i));
	}
      }
      findBestInterface(C3_u2,C2_u2,C2_l_u2,C1_u2,
			right_C_r_u2C2,C_r_u2C2,best_R_score,
			bnd_heur_v,
			findBestBoundary,
			// find best boundary args
			tri_heur,
			E_u1,
			C2_u1,
			Cextra_u1,
			C1_u1,
			P_u1,
			C2_u2_to_C2_u1,
			C2_u2_to_C1_u1
			);
    }
  }

  // Now find the partitions (i.e., left or right) corresponding
  // the interface which had minimum size, prefering the left
  // interface if there is a tie.
  if ((flr.size() == 0 && best_L_score <= best_R_score)
      || 
      (flr.size() > 0 && toupper(flr[0]) == 'L')) {
    // this next routine gives us the best left interface that
    // exists from within the chunk C2_u2 and places
    // it in C_l_u2, and everything to the 'left' of C_l_u2
    // that still lies within C2_u2 is placed in left_C_l_u2
    findingLeftInterface = true;
    infoMsg(Tiny,"---\nUsing left interface to define partitions\n");

    findInterfacePartitions(P_u1,
			    C1_u1,
			    Cextra_u1,
			    C2_u1,
			    E_u1,
			    C2_u2_to_C1_u1,
			    C2_u2_to_C2_u1,
			    left_C_l_u2C2,
			    C_l_u2C2,
			    gm_template);

    // Write information about how boundary was created to a string,
    // but make sure there is no white space in the string.
    char buff[2048];
    sprintf(buff,"Left_interface:Run_Bdry_Alg(%c),Bnd_Heurs(%s),TravFrac(%f),Tri_Heur(%s)",
	    (findBestBoundary? 'T' : 'F'),
	    bnd_heur_str.c_str(),
	    boundaryTraverseFraction,
	    tri_heur_str.c_str());
    gm_template.boundaryMethod = buff;
    gm_template.leftInterface = true; 
  } else {
    // find right interface partitions
    findingLeftInterface = false;
    infoMsg(Nano,"---\nUsing right interface to define partitions\n");


    // create a mapping from C1_u2 to C1_u1 and from C1_u2 to P_u1
    // (where a correspondence exists)

    findInterfacePartitions(E_u1,
			    C2_u1,
			    Cextra_u1,
			    C1_u1,
			    P_u1,
			    C2_u2_to_C2_u1,
			    C2_u2_to_C1_u1,
			    right_C_r_u2C2,
			    C_r_u2C2,
			    gm_template);

    // Write information about how boundary was created to a string,
    // but make sure there is no white space in the string.
    char buff[2048];
    sprintf(buff,"Right_interface:Run_Bdry_Alg(%c),Bnd_Heurs(%s),TravFrac(%f),Tri_Heur(%s)",
	    (findBestBoundary? 'T' : 'F'),
	    bnd_heur_str.c_str(),
	    boundaryTraverseFraction,
	    tri_heur_str.c_str());
    gm_template.boundaryMethod = buff;
    gm_template.leftInterface = false; 

  }
}


/*-
 *-----------------------------------------------------------------------
 * BoundaryTriangulate::findInterfacePartitions()
 *   Create the three partitions, either left or right depending
 *   on the order of the arguments given.
 *
 * For the left interface, we create new P,C, and E variable sets where
 *  where P = modified prologue
 *  where C = modified chunk to repeat
 *  where E = modified epilogue to repeat
 *
 * Preconditions:
 *   The variable findingLeftInterface must be set to
 *   appropriate value before calling this routine.
 *
 * Postconditions:
 *
 * Side Effects:
 *   Makes the interfaces in C1_u1 and C2_u1 complete.
 *
 * Results:
 *
 *
 *-----------------------------------------------------------------------
 */
void
BoundaryTriangulate::
findInterfacePartitions(
 // input variables
 const set<RandomVariable*>& P_u1,
 const set<RandomVariable*>& C1_u1,
 const set<RandomVariable*>& Cextra_u1, // non-empty only when S > M
 const set<RandomVariable*>& C2_u1,
 const set<RandomVariable*>& E_u1,
 // these next 2 should be const, but there is no "op[] const"
 map < RandomVariable*, RandomVariable* >& C2_u2_to_C1_u1,
 map < RandomVariable*, RandomVariable* >& C2_u2_to_C2_u1,
 const set<RandomVariable*>& left_C_l_u2C2,
 const set<RandomVariable*>& C_l_u2C2,
 // output variables
 GMTemplate& gm_template)

{
  // now we need to make a bunch of sets to be unioned
  // together to get the partitions.
  set<RandomVariable*> C_l_u1C1;
  set<RandomVariable*> C_l_u1C2;

  set<RandomVariable*> left_C_l_u1C1;
  set<RandomVariable*> left_C_l_u1C2;
  if (M > 0) {
    for (set<RandomVariable*>::iterator i = C_l_u2C2.begin();
	 i!= C_l_u2C2.end(); i++) {

      assert (C2_u2_to_C1_u1.find((*i)) != C2_u2_to_C1_u1.end());
      assert (C2_u2_to_C2_u1.find((*i)) != C2_u2_to_C2_u1.end());

      C_l_u1C1.insert(C2_u2_to_C1_u1[(*i)]);
      C_l_u1C2.insert(C2_u2_to_C2_u1[(*i)]);
    }

    for (set<RandomVariable*>::iterator i = left_C_l_u2C2.begin();
	 i != left_C_l_u2C2.end(); i++) {
      
      assert (C2_u2_to_C1_u1.find((*i)) != C2_u2_to_C1_u1.end());
      assert (C2_u2_to_C2_u1.find((*i)) != C2_u2_to_C2_u1.end());
    
      left_C_l_u1C1.insert(C2_u2_to_C1_u1[(*i)]);
      left_C_l_u1C2.insert(C2_u2_to_C2_u1[(*i)]);
    }
  } else {
    // M == 0
    error("findInterfacePartitions: M==0 case not implemented, use unroll and triangulate\n");
    // this case is needed presumably because we want to do something
    // where the original GMTK template is unrolled 0 times. Rather
    // than doing this here, we will just use the unroll and triangulate
    // method and then do inference on the resulting cliques.
  }

  
  // Finally, create the modified sets P, C, and E
  //   where P = modified prologue
  //   where C = modified chunk to repeat
  //   where E = modified epilogue to repeat.
  // which are to be triangulated separately. 
  // For the *left* interface, the modified sets are defined
  // as follows:
  //    P = P_u1 + C1_u1(left_C_l_u2) + C1_u1(C_l_u2)
  //    E = C2_u1\C2_u1(left_C_l_u1) + E_u1
  //    C = ((C1_u1 + C2_u1) \ (P + E)) + C2_u1(C_l_u2) + C1_u1(C_l_u2) + Cextra_u1
  // and the symmetric definitions apply for the right interface.
  // We use the left interface definitions in this code and
  // assume the caller calles with inverted arguments
  // to get the right interface behavior.

  // Finish P
  set<RandomVariable*> P = P_u1;
  set_union(left_C_l_u1C1.begin(),left_C_l_u1C1.end(),
	    C_l_u1C1.begin(),C_l_u1C1.end(),
	    inserter(P,P.end()));

  // Finish E
  set<RandomVariable*> E = E_u1;
  set_difference(C2_u1.begin(),C2_u1.end(),
		 left_C_l_u1C2.begin(),left_C_l_u1C2.end(),
		 inserter(E,E.end()));

  // Finish C
  set<RandomVariable*> C;
  set<RandomVariable*> tmp1;
  set<RandomVariable*> tmp2;
  set<RandomVariable*> tmp3;
  set_union(C1_u1.begin(),C1_u1.end(),
	    C2_u1.begin(),C2_u1.end(),
	    inserter(tmp1,tmp1.end()));
  set_union(P.begin(),P.end(),
	    E.begin(),E.end(),
	    inserter(tmp2,tmp2.end()));
  set_union(C_l_u1C2.begin(),C_l_u1C2.end(),
	    C_l_u1C1.begin(),C_l_u1C1.end(),
	    inserter(tmp3,tmp3.end()));
  set_union(Cextra_u1.begin(),Cextra_u1.end(),
	    tmp3.begin(),tmp3.end(),
	    inserter(C,C.end()));
  set_difference(tmp1.begin(),tmp1.end(),
		 tmp2.begin(),tmp2.end(),
		 inserter(C,C.end()));

  // create the template with these partition definitions.
  if (findingLeftInterface)
    gm_template.createPartitions(P,C,E,C_l_u1C1,C_l_u1C2);
  else // we are finding right interface, need to swap arguments. 
    gm_template.createPartitions(E,C,P,C_l_u1C2,C_l_u1C1);

}


/*******************************************************************************
 *******************************************************************************
 *******************************************************************************
 **
 **         Main Triangulation Routines
 **
 *******************************************************************************
 *******************************************************************************
 *******************************************************************************
 */


/*-
 *-----------------------------------------------------------------------
 * BoundaryTriangulate::triangulate()
 *   Triangulate an entire GMTemplate, using the same method for each partition.
 *    
 * Preconditions:
 *   Graphs in GMTemplate must be a valid undirected model. This means if the graph
 *   was originally a directed model, it must have been properly
 *   moralized and their 'neighbors' structure is valid. It is also
 *   assumed that the parents of each r.v. are valid but that they
 *   only poiint to variables within the set 'nodes' (i.e., parents
 *   must not point out of this set).
 *
 *   Also, it is assumed that graphs are not yet triangulated. If they
 *   are, the result will still be triangulated, but it might add more
 *   edgeds to the already triangulated graph, depending on the
 *   triangulation heuristic given.
 *
 * Postconditions:
 *   Resulting graphs are  now triangulated, and cliques are stored (not in RIP order)
 *
 * Side Effects:
 *   Will change neighbors members of variables.
 *
 * Results:
 *     none
 *
 *
 *-----------------------------------------------------------------------
 */
void
BoundaryTriangulate::
triangulate(const string& tri_heur_str,
	    bool jtWeight,
	    GMTemplate& gm_template,
	    bool doP,
	    bool doC,
	    bool doE)
{
  TriangulateHeuristics tri_heur;
  vector<nghbrPairType> orgnl_P_nghbrs;
  vector<nghbrPairType> orgnl_C_nghbrs;
  vector<nghbrPairType> orgnl_E_nghbrs;
  string best_P_method_str;
  string best_C_method_str;
  string best_E_method_str;
  double best_P_weight = HUGE_VAL;
  double best_C_weight = HUGE_VAL;
  double best_E_weight = HUGE_VAL;
  const set <RandomVariable*> emptySet;

  parseTriHeuristicString(tri_heur_str,tri_heur);

  if (doP)
    saveCurrentNeighbors(gm_template.P,orgnl_P_nghbrs);
  if (doC)
    saveCurrentNeighbors(gm_template.C,orgnl_C_nghbrs);
  if (doE)
    saveCurrentNeighbors(gm_template.E,orgnl_E_nghbrs);

  if (doP) {
    infoMsg(IM::Tiny, "---\nTriangulating P:\n");
    triangulate(gm_template.P.nodes,jtWeight,gm_template.PCInterface_in_P,tri_heur,orgnl_P_nghbrs,gm_template.P.cliques,gm_template.P.triMethod,best_P_weight);
  }
  if (doC) {
    infoMsg(IM::Tiny, "---\nTriangulating C:\n");
    triangulate(gm_template.C.nodes,jtWeight,gm_template.CEInterface_in_C,tri_heur,orgnl_C_nghbrs,gm_template.C.cliques,gm_template.C.triMethod,best_C_weight);
  }
  if (doE) {
    infoMsg(IM::Tiny, "---\nTriangulating E:\n");
    triangulate(gm_template.E.nodes,jtWeight,emptySet,tri_heur,orgnl_E_nghbrs,gm_template.E.cliques,gm_template.E.triMethod,best_E_weight);
  }

  ////////////////////////////////////////////////////////////////////////
  // Return with the best triangulations found, which is
  // be stored within the template at this point.
  ////////////////////////////////////////////////////////////////////////
  if (doP)
    restoreNeighbors(orgnl_P_nghbrs);
  if (doC)
    restoreNeighbors(orgnl_C_nghbrs);
  if (doE)
    restoreNeighbors(orgnl_E_nghbrs);
  gm_template.triangulatePartitionsByCliqueCompletion();

}



/*-
 *-----------------------------------------------------------------------
 * BoundaryTriangulate::triangulate()
 *   The actual triangulation that does the work of triangulation.
 *  
 *    
 * Preconditions:
 *   Graph must be a valid undirected model. This means if the graph
 *   was originally a directed model, it must have been properly
 *   moralized and their 'neighbors' structure is valid. It is also
 *   assumed that the parents of each r.v. are valid but that they
 *   only poiint to variables within the set 'nodes' (i.e., parents
 *   must not point out of this set).
 *
 * Postconditions:
 *   Resulting graph is now triangulated, and cliques are returned (not in RIP order)
 *
 * Side Effects:
 *   Will change neighbors members of variables, even if the current
 *   triangulation does not beat best_weight.
 *
 * Results:
 *     none
 *
 *
 *-----------------------------------------------------------------------
 */
void
BoundaryTriangulate::
triangulate(// input: nodes to be triangulated
	    const set<RandomVariable*>& nodes,
	    // use JT weight rather than sum of weight
	    const bool jtWeight,
	    // nodes that a JT root must contain (ok to be empty).
	    const set<RandomVariable*>& nodesRootMustContain,
	    // triangulation heuristic method
	    const TriangulateHeuristics& tri_heur,
	    // original neighbor structures
	    vector<nghbrPairType>& orgnl_nghbrs,
	    // output: resulting max cliques
	    vector<MaxClique>& best_cliques,
	    // output: string giving resulting method used
	    string& best_meth_str,
	    // weight to best
	    double& best_weight)

{
  infoMsg(Huge,"\nBEGINNING TRIANGULATION --- \n");

  vector<MaxClique>       cliques;
  vector<RandomVariable*> order;
  double                  weight;
  string                  meth_str;

  // compute the real best weight for a set of current
  // cliques, if the weight has not already been computed.
  if (best_weight == HUGE_VAL && best_cliques.size() > 0) {
    best_weight = graphWeight(best_cliques,jtWeight,nodesRootMustContain);
  }


  for (unsigned trial = 0;trial<tri_heur.numberTrials;trial++) {
    string annealing_str;
    char   buff[BUFSIZ];

    sprintf(buff,"%d",trial);
    restoreNeighbors(orgnl_nghbrs);
    if (tri_heur.style == TS_ANNEALING) {
      triangulateSimulatedAnnealing(nodes, jtWeight, nodesRootMustContain,
				    cliques, 
				    order, 
                                    annealing_str );

      meth_str = string(buff) + "-" + "annealing-" + annealing_str;
    } else if (tri_heur.style == TS_EXHAUSTIVE) {
      triangulateExhaustiveSearch( nodes, jtWeight, nodesRootMustContain,
				   orgnl_nghbrs, cliques ); 
      meth_str = string(buff) + "-" + "exhaustive";
    } else if (tri_heur.style == TS_PRE_EDGE_ALL) {
      preEdgeAdditionElimination( nodes, jtWeight, nodesRootMustContain, 
        ALL_EDGES, cliques, meth_str );
    } else if (tri_heur.style == TS_PRE_EDGE_LO) {
      preEdgeAdditionElimination( nodes, jtWeight, nodesRootMustContain, 
        LOCALLY_OPTIMAL_EDGES, cliques, meth_str );
    } else if (tri_heur.style == TS_PRE_EDGE_RANDOM) {
      preEdgeAdditionElimination( nodes, jtWeight, nodesRootMustContain, 
        RANDOM_EDGES, cliques, meth_str );
    } else if (tri_heur.style == TS_ELIMINATION_HEURISTICS) {
      vector<nghbrPairType>   orgnl_nghbrs;
      saveCurrentNeighbors(nodes, orgnl_nghbrs);

      tryEliminationHeuristics( nodes, jtWeight, nodesRootMustContain, 
        orgnl_nghbrs, cliques, meth_str );
    } else if (tri_heur.style == TS_MCS) {
      triangulateMaximumCardinalitySearch(nodes, 
					  cliques, order );
      meth_str = string(buff) + "-" + "MCS";
    } else if (tri_heur.style == TS_FRONTIER) {
      triangulateFrontier(nodes, 
			  cliques);
      meth_str = string(buff) + "-" + "FRONTIER";
    } else if (tri_heur.style == TS_COMPLETED) {
      triangulateCompletePartition( nodes, cliques );
      meth_str = string(buff) + "-" + "completed";
    } else if (tri_heur.style == TS_BASIC) {
      basicTriangulate( nodes, tri_heur.heuristic_vector, order, cliques);
      meth_str = string(buff) + "-" + tri_heur.basic_method_string;
    } else {
      // shouldn't happen
      assert (0);
    }

    weight = graphWeight(cliques,jtWeight,nodesRootMustContain);
    if (weight < best_weight) {

      /////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////
      // @@@ There is a bug due to the stubbed out = operator in the MaxClique 
      //   class the following line which clears the cliques appears to work 
      //   around the bug. @@@
      // TODO: ultimately take this out, but keep in for now until
      // we do code restructuring (1/20/2004).
      best_cliques.clear();
      /////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////



      best_cliques = cliques;
      best_weight  = weight;
      best_meth_str = meth_str;

      infoMsg(IM::Tiny, "***New Best: %10s %-10f\n", meth_str.c_str(), weight);

      cliques.clear();
      order.clear();
    }
    if (timer && timer->Expired()) {
      infoMsg(IM::Tiny, "Exiting triangulation algorithm because time has expired\n");
      break;
    }
  }

  infoMsg(Huge,"\nENDING TRIANGULATION --- \n");
}


/*-
 *-----------------------------------------------------------------------
 * BoundaryTriangulate::basicTriangulate()
 *   The actual basic triangulation that does the work of triangulation.
 *  
 *   This routine will triangulate a set of nodes using any
 *   combination of a number if different (but simple) triangulation
 *   heuristics such as (min weight, min size, min fill, etc.).  For a
 *   good description of these heuristics, see D. Rose et. al, 1970,
 *   1976. The routine also allows for other heuristics to be used
 *   such as eliminate the earlier nodes (temporal order) first, or
 *   eliminate the nodes in order that they appear in the structure
 *   file (sometimes this simple constrained triangulation will work
 *   better than the "intelligent" heuristics, such as for certain
 *   lattice structures). The routine allows heuristics to be
 *   prioritized and combined, so that if there is a tie with the
 *   first heuristic, the second will be used, and if there is still a
 *   tie, the third one will be used, etc.
 *    
 * Preconditions:
 *   Graph must be a valid undirected model and untriangulated. This means if the graph
 *   was originally a directed model, it must have been properly
 *   moralized and their 'neighbors' structure is valid. It is also
 *   assumed that the parents of each r.v. are valid but that they
 *   only poiint to variables within the set 'nodes' (i.e., parents
 *   must not point out of this set).
 *
 * Postconditions:
 *   Resulting graph is now triangulated.
 *
 * Side Effects:
 *   Might (and probably will unless graph is already triangulated and
 *   you get lucky by having found the perfect elimination order)
 *   change neighbors members of variables.
 *
 * Results:
 *     none
 *
 *
 *-----------------------------------------------------------------------
 */
void
BoundaryTriangulate::
basicTriangulate(// input: nodes to triangulate
		 const set<RandomVariable*>& nodes,
		 // input: triangulation heuristic
		 const vector<BasicTriangulateHeuristic>& th_v,
		 // output: nodes ordered according to resulting elimination
		 vector<RandomVariable*>& orderedNodes,  
		 // output: resulting max cliques
		 vector<MaxClique>& cliques,
		 // input: find the cliques as well
		 const bool findCliques
		 )
{
  const unsigned num_nodes = nodes.size();

  infoMsg(Huge,"\nBEGINNING BASIC TRIANGULATION --- \n");

  // need to use some form of heuristic
  assert ( th_v.size() != 0 );

  cliques.clear();

  // orderedNodes == already eliminated nodes. 
  orderedNodes.clear();

  // Also keep ordered (eliminated) nodes as a set for easy
  // intersection, with other node sets.
  set<RandomVariable*> orderedNodesSet;

  
  // Approach: essentially, create a priority queue data structure
  // using a multimap. In this case, those nodes with the lowest 'X'
  // where 'X' is the combined prioritized weight heuristics can be
  // accessed immediately in the front of the map. Also, whenever a
  // node gets eliminated, *ONLY* its neighbor's weights are
  // recalculated. With this data structure it is possible and
  // efficient to do so. This is because a multimap (used to simulate
  // a priority queue) has the ability to efficiently remove stuff
  // from the middle.
  multimap< vector<float> ,RandomVariable*> unorderedNodes;

  // We also need to keep a map to be able to remove elements of
  // 'unorderedNodes' when a node is eliminated. I.e., we need to be
  // able to map back from a RV* directly to its entry in the priority
  // queue so that when a node is eliminated, its neighbors can be
  // removed from the queue (since their weight is now invalid) and
  // then (only) their weight can be recalculated anew.
  map<RandomVariable*, multimap< vector<float>,RandomVariable*>::iterator > 
    rv2unNodesMap;

  // Also, create a set of nodes which are the ones whose weight
  // needs to be updated in the priority queue 'unorderedNodes'.
  // We begin by updating the weights of all nodes.
  set<RandomVariable*> nodesToUpdate = nodes;

  do {

    for (set<RandomVariable*>::iterator i = nodesToUpdate.begin();
	 i != nodesToUpdate.end();
	 i++) {

      infoMsg(Huge,"TR: computing weight of node %s(%d)\n",
	      (*i)->name().c_str(),(*i)->frame());

      // Create a vector with (weight,fillin,timeframe, etc.)
      // and choose in increasing lexigraphic order in that
      // tuple.
      vector<float> weight;

      // Create activeNeighbors, which contains only those neighbors
      // that are active and are not yet eliminated.
      set<RandomVariable*> activeNeighbors;
      set_difference((*i)->neighbors.begin(),(*i)->neighbors.end(),
		     orderedNodesSet.begin(),orderedNodesSet.end(),
		     inserter(activeNeighbors,activeNeighbors.end()));

      for (unsigned thi=0;thi<th_v.size();thi++) {

	const BasicTriangulateHeuristic th = th_v[thi];

	if (th == TH_MIN_WEIGHT || th == TH_MIN_WEIGHT_NO_D) {
	  float tmp_weight = MaxClique::computeWeight(activeNeighbors,(*i),
					   (th == TH_MIN_WEIGHT));
	  weight.push_back(tmp_weight);
	  infoMsg(Huge,"  node has weight = %f\n",tmp_weight);
	} else if (th == TH_MIN_FILLIN) {
	  int fill_in = computeFillIn(activeNeighbors);
	  weight.push_back((float)fill_in);
	  infoMsg(Huge,"  node has fill_in = %d\n",fill_in);
	} else if (th == TH_MIN_TIMEFRAME) {
	  weight.push_back((*i)->frame());
	  infoMsg(Huge,"  node has time frame = %d\n",(*i)->frame());
	} else if (th == TH_MAX_TIMEFRAME) {
	  weight.push_back( - ((*i)->frame()));
	  infoMsg(Huge,"  node has (neg) time frame = -%d\n",(*i)->frame());
	} else if (th == TH_MIN_SIZE) {
	  weight.push_back((float)activeNeighbors.size());
	  infoMsg(Huge,"  node has active neighbor size = %d\n",
		  activeNeighbors.size());
	} else if (th == TH_MIN_POSITION_IN_FILE) {
	  weight.push_back((float)(*i)->rv_info.variablePositionInStrFile);
	  infoMsg(Huge,"  node has position in file = %d\n",
		  (*i)->rv_info.variablePositionInStrFile);
	} else if (th == TH_MIN_HINT) {
	  weight.push_back((float)(*i)->rv_info.eliminationOrderHint);
	  infoMsg(Huge,"  node has elimination order hint = %f\n",
		  (*i)->rv_info.eliminationOrderHint);
	} else if (th == TH_RANDOM) {
	  float tmp = rnd.drand48();
	  weight.push_back(tmp);
	  infoMsg(Huge,"  node has random value = %f\n",tmp);
	} else
	  warning("Warning: unimplemented triangulation heuristic (ignored)\n");
      }

      pair< vector<float>,RandomVariable*> p(weight,(*i));
      rv2unNodesMap[(*i)] = (unorderedNodes.insert(p));
    }

    if (message(Huge)) {
      // go through and print sorted guys.
      printf("Order of nodes to be eliminated\n");
      for (multimap< vector<float> ,RandomVariable*>::iterator m
	     = unorderedNodes.begin();
	   m != unorderedNodes.end(); m++) {
	printf("Node %s(%d) with weights:",
	       (*m).second->name().c_str(),
	       (*m).second->frame());
	for (unsigned l = 0; l < (*m).first.size(); l++ )
	  printf(" %f,",(*m).first[l]);
	printf("\n");
      }
    }


    // Go through the updated multi-map, grab an iterator pair to
    // iterate between all nodes that have the lowest weight. This
    // utilizes the fact that the multimap stores values in ascending
    // order based on key (in this case the weight), and so the first
    // one (i.e., mm.begin() ) should have the lowest weight.
    pair< 
      multimap< vector<float>,RandomVariable*>::iterator, 
      multimap< vector<float>,RandomVariable*>::iterator 
      > ip = unorderedNodes.equal_range( (*(unorderedNodes.begin())).first );

    const int d = distance(ip.first,ip.second);
    if (d == 1) {
      // then there is only one node and we eliminate that
      // so do nothing
    } else {
      // choose a random one, don't re-seed rng as that
      // should be done one time for the program via command line arguments.
      RAND rnd(false);
      int val = rnd.uniform(d-1);
      // TODO: there must be a better way to do this next step!
      while (val--)
	ip.first++;
    }

    // ip.first now points to the pair containing the random variable that
    // we eliminate.
    RandomVariable *rv = (*(ip.first)).second;

    if (message(Huge)) {
      printf("\nEliminating node %s(%d) with weights:",
	     rv->name().c_str(),rv->frame());
      for (unsigned l = 0; l < (*(ip.first)).first.size(); l++ )
	printf(" %f,",(*(ip.first)).first[l]);
      printf("\n");
    }


    // connect all neighbors of r.v. excluding nodes in 'orderedNodesSet'.
    rv->connectNeighbors(orderedNodesSet);

    // find the cliques if they are asked for.
    if (findCliques) {
      // check here if this node + its neighbors is a subset of
      // previous maxcliques. If it is not a subset of any previous
      // maxclique, then this node and its neighbors is a new
      // maxclique.
      set<RandomVariable*> candidateMaxClique;
      set_difference(rv->neighbors.begin(),rv->neighbors.end(),
		     orderedNodesSet.begin(),orderedNodesSet.end(),
		     inserter(candidateMaxClique,candidateMaxClique.end()));
      candidateMaxClique.insert(rv);
      bool is_max_clique = true;
      for (unsigned i=0;i < cliques.size(); i++) {
	if (includes(cliques[i].nodes.begin(),cliques[i].nodes.end(),
		     candidateMaxClique.begin(),candidateMaxClique.end())) {
	  // then found a 'proven' maxclique that includes our current
	  // candidate, so the candidate cannot be a maxclique
	  is_max_clique = false;
	  break;
	}
      }
      if (is_max_clique) {
	if (message(Huge)) {
	  // print out clique information.
	  printf("Found a max clique of size %d while eliminating node %s(%d):",
		 candidateMaxClique.size(),
		 rv->name().c_str(),rv->frame());
	  for (set<RandomVariable*>::iterator j=candidateMaxClique.begin();
	       j != candidateMaxClique.end(); j++) {
	    RandomVariable* rv = (*j);
	    printf(" %s(%d)",rv->name().c_str(), rv->frame());
	  }
	  printf("\n");
	}
 	cliques.push_back(MaxClique(candidateMaxClique));
      }
    }
      
    // insert node into ordered list
    orderedNodes.push_back(rv);

    // insert node into ordered set
    orderedNodesSet.insert(rv);

    // erase node from priority queue
    unorderedNodes.erase(ip.first);

    // only update not-yet-eliminated nodes that could possibly have
    // been effected by the current node 'rv' being eliminated. I.e.,
    // create the set of active neighbors of rv.
    nodesToUpdate.clear();
    set_difference(rv->neighbors.begin(),rv->neighbors.end(),
		   orderedNodesSet.begin(),orderedNodesSet.end(),
		   inserter(nodesToUpdate,nodesToUpdate.end()));

    // erase active neighbors of nodes since they will need to be
    // recomputed above.
    for (set<RandomVariable*>::iterator n = nodesToUpdate.begin();
	 n != nodesToUpdate.end();
	 n++) {
      unorderedNodes.erase(rv2unNodesMap[(*n)]);
    }

    // continue until all nodes are eliminated.
  } while (orderedNodesSet.size() < num_nodes);

  infoMsg(Huge,"\nENDING BASIC TRIANGULATION --- \n");
}


/*-
 *-----------------------------------------------------------------------
 * BoundaryTriangulate::
 *  fillAccordingToCliques 
 *  
 * Preconditions:
 *  If one wants the graph to match the given cliques exactly, the 
 *  RandomVariables should not have any neighbor which are present in the 
 *  input cliques as this procedure does not remove any edges.  Typically 
 *  this procedure is preceeded by resotring the original graph neighbors 
 *  or removing all neighbors, and the vector of cliques corresponds to a
 *  triangulation.
 *
 * Postconditions:
 *   Graph is triangulated according to the best weight found.
 *   
 * Side Effects:
 *
 * Results:
 *     none
 *
 *
 *-----------------------------------------------------------------------
 */
void
BoundaryTriangulate::
fillAccordingToCliques(
  const vector<MaxClique>& cliques
  ) 
{
  //////////////////////////////////////////////////////////////////////////
  // Convert to triangulate the RandomVariable structures
  //////////////////////////////////////////////////////////////////////////
  vector<MaxClique>::const_iterator crrnt_clique;
  vector<MaxClique>::const_iterator end_clique;

  for ( crrnt_clique = cliques.begin(), 
        end_clique   = cliques.end();
        crrnt_clique != end_clique;
        ++crrnt_clique ) {
    MaxClique::makeComplete((*crrnt_clique).nodes);
  }
}


/*-
 *-----------------------------------------------------------------------
 * BoundaryTriangulate::preEdgeAdditionElimination()
 *   Triangulates a graph by adding extra edges followed by the 
 *   elimination heuristics. 
 *  
 * Preconditions:
 *   Each variable in the set of nodes must have valid parent and 
 *   neighbor members and the parents/neighbors must only point to other 
 *   nodes in the set. 
 *
 * Postconditions:
 *   Graph is triangulated according to the best weight found.
 *   The triangulation is stored in 'cliques'.
 *   The method that found the best triangulation is written to best_method. 
 *
 * Side Effects:
 *   none
 *
 * Results:
 *   none
 *-----------------------------------------------------------------------
 */
void
BoundaryTriangulate::
preEdgeAdditionElimination(
  const set<RandomVariable*>& nodes,           // input: nodes to triangulate
  const bool                  jtWeight,        // If true, use jtWeight 
  const set<RandomVariable*>& nodesRootMustContain,
  const extraEdgeHeuristic    edge_heuristic, 
  vector<MaxClique>&          cliques,         // output: resulting max cliques
  string&                     best_method
  )
{
  vector<triangulateNode> triangulate_nodes;
  vector<nghbrPairType>   orgnl_nghbrs;
  string tri_method;

  infoMsg(IM::Tiny, "--- adding ancestral edges ---\n");

  //////////////////////////////////////////////////////////////////////
  // Initialize
  //////////////////////////////////////////////////////////////////////
  fillTriangulateNodeStructures( nodes, triangulate_nodes );

  //////////////////////////////////////////////////////////////////////
  // Add the extra edges accoring to 'edge_heuristic' 
  //////////////////////////////////////////////////////////////////////
  addExtraEdgesToGraph( triangulate_nodes, edge_heuristic );

  //////////////////////////////////////////////////////////////////////
  // Now eliminate using the heuristics 
  //////////////////////////////////////////////////////////////////////
  saveCurrentNeighbors(nodes, orgnl_nghbrs);
  tryEliminationHeuristics( nodes, jtWeight, nodesRootMustContain, orgnl_nghbrs,
    cliques, tri_method );

  //////////////////////////////////////////////////////////////////////
  // Triangulate according to best triangulation 
  //////////////////////////////////////////////////////////////////////
  restoreNeighbors( orgnl_nghbrs );
  fillAccordingToCliques( cliques );

  //////////////////////////////////////////////////////////////////////
  // Record the best method 
  //////////////////////////////////////////////////////////////////////
  best_method = "pre-edge-";
  switch (edge_heuristic) {
    case ALL_EDGES:
      best_method += "all";
      break; 
    case LOCALLY_OPTIMAL_EDGES:
      best_method += "lo";
      break; 
    case RANDOM_EDGES:
      best_method += "random";
      break; 
    default:
      assert(0);
      break;
  }
  best_method += "-" + tri_method;  
}


/*-
 *-----------------------------------------------------------------------
 * BoundaryTriangulate::fillParentChildLists
 *   Fill containers in the triangulateNode structure that list each 
 *   node's parents and nonChildren.  nonChildren are the union of the 
 *   parents and undirected neighbors.  Note that these lists will be 
 *   different from the ones contained the associated RandomVariable 
 *   because the containers created here don't include nodes from outside 
 *   of the partition. 
 * 
 * Preconditions:
 *   The nodes given as input are initialized 
 *
 * Postconditions:
 *   The parents and nonChildren containers of each node are filled in.
 * 
 * Side Effects:
 *   none
 *
 * Results:
 *   none
 *-----------------------------------------------------------------------
 */
void
BoundaryTriangulate::
fillParentChildLists(
  vector<triangulateNode>& nodes
  )
{
  vector<triangulateNode>::iterator crrnt_node;
  vector<triangulateNode>::iterator end_node;
  triangulateNeighborType::iterator crrnt_nghbr;
  triangulateNeighborType::iterator end_nghbr;
  vector<RandomVariable*>::iterator found_node;

  for( crrnt_node = nodes.begin(),
       end_node   = nodes.end();
       crrnt_node != end_node;
       ++crrnt_node ) {

    for( crrnt_nghbr = (*crrnt_node).neighbors.begin(), 
         end_nghbr   = (*crrnt_node).neighbors.end(); 
         crrnt_nghbr != end_nghbr;
         ++crrnt_nghbr ) {

      found_node = find( 
        (*crrnt_node).randomVariable->allPossibleParents.begin(), 
        (*crrnt_node).randomVariable->allPossibleParents.end(),
        (*crrnt_nghbr)->randomVariable );
 
      //////////////////////////////////////////////////////////////
      // If a parent, add to both parents and non-children
      //////////////////////////////////////////////////////////////
      if (found_node != 
          (*crrnt_node).randomVariable->allPossibleParents.end()) {
        (*crrnt_node).parents.push_back( *crrnt_nghbr );
        (*crrnt_node).nonChildren.push_back( *crrnt_nghbr );
      }
      //////////////////////////////////////////////////////////////
      // If not a parent, check if it is a non-child
      //////////////////////////////////////////////////////////////
      else {
        found_node = find( 
          (*crrnt_node).randomVariable->allPossibleChildren.begin(), 
          (*crrnt_node).randomVariable->allPossibleChildren.end(),
          (*crrnt_nghbr)->randomVariable ); 

        if (found_node == 
            (*crrnt_node).randomVariable->allPossibleChildren.end()) {
          (*crrnt_node).nonChildren.push_back( *crrnt_nghbr );
        }
      }
    }
  } 

}


/*-
 *-----------------------------------------------------------------------
 * BoundaryTriangulate::addExtraEdgesToGraph()
 *   Add extra edges to the graph according to a particular heuristic
 *  
 * Preconditions:
 *   Each variable in the set of nodes must have valid parent and 
 *   neighbor members and the parents/neighbors must only point to other 
 *   nodes in the set. 
 *
 * Postconditions:
 *   The nodes may have additional neighbors.
 * 
 * Side Effects:
 *   The parents and nonChildren members of the triangulateNode 
 *   structures are filled in. 
 *
 * Results:
 *   none
 *-----------------------------------------------------------------------
 */
void
BoundaryTriangulate::
addExtraEdgesToGraph(
  vector<triangulateNode>& nodes,
  const extraEdgeHeuristic edge_heuristic 
  )
{
  vector<triangulateNode>::iterator crrnt_node;
  vector<triangulateNode>::iterator end_node;
  vector<triangulateNode>::iterator crrnt_mark;
  vector<triangulateNode>::iterator end_mark;
  triangulateNeighborType::iterator crrnt_nghbr;
  triangulateNeighborType::iterator end_nghbr;
  vector<edge> extra_edges;

  fillParentChildLists( nodes );

  //////////////////////////////////////////////////////////////////////////
  // Iterate through all nodes
  //////////////////////////////////////////////////////////////////////////
  for( crrnt_node = nodes.begin(),
       end_node   = nodes.end();
       crrnt_node != end_node;
       ++crrnt_node ) {

    //////////////////////////////////////////////////////////////////////////
    // Mark all nodes as not visited 
    //////////////////////////////////////////////////////////////////////////
    for( crrnt_mark = nodes.begin(),
         end_mark   = nodes.end();
         crrnt_mark != end_mark;
         ++crrnt_mark ) {
      (*crrnt_mark).marked = false;
    }

    (*crrnt_node).marked = true;

    //////////////////////////////////////////////////////////////////////////
    // Find all nodes which are parents or undirected neighbors, and are 
    //   deterministic or sparse 
    //////////////////////////////////////////////////////////////////////////
    for( crrnt_nghbr = (*crrnt_node).nonChildren.begin(), 
         end_nghbr   = (*crrnt_node).nonChildren.end(); 
         crrnt_nghbr != end_nghbr;
         ++crrnt_nghbr ) {

      (*crrnt_nghbr)->marked = true;

      if  ((*crrnt_nghbr)->parents.size() > 0)  { 
        if ( ((*crrnt_nghbr)->randomVariable->deterministic()) || 
             (((*crrnt_nghbr)->randomVariable->discrete) && 
              (((DiscreteRandomVariable*)
               ((*crrnt_nghbr)->randomVariable))->sparse())) ) {

          addEdgesToNode((*crrnt_nghbr)->parents, *crrnt_nghbr, &(*crrnt_node), 
            edge_heuristic, extra_edges);
        }
      }
    }
  } 

  addEdges(extra_edges);
}


/*-
 *-----------------------------------------------------------------------
 * BoundaryTriangulate::addEdgesToNode
 *   Choose if an edge should be added from a set of nodes to another 
 *   node.  This is intended to be a set of parents of child connecting  
 *   to a grandchild for use in addExtraEdgesToGraph.  If the edge is 
 *   added, this procedure calls itself recursively. 
 *  
 * Preconditions:
 *   Each variable in the set of nodes must have valid parent and 
 *   neighbor members and the parents/neighbors must only point to other 
 *   nodes in the set. 
 * 
 * Postconditions:
 *   The choice of edges is stored in 'extra_edges'.  These will be from 
 *   'grandchild' to any other node in the graph.  The neighbor sets of 
 *   the nodes are not modified.  
 *
 * Side Effects:
 *   none
 *
 * Results:
 *   none
 *-----------------------------------------------------------------------
 */
void
BoundaryTriangulate::
addEdgesToNode(
  triangulateNeighborType& parent_set,  
  triangulateNode* const   child, 
  triangulateNode* const   grandchild,
  const extraEdgeHeuristic edge_heuristic, 
  vector<edge>&            extra_edges
  )
{
  triangulateNeighborType::iterator crrnt_prnt;
  triangulateNeighborType::iterator end_prnt;
  float parent_weight, child_weight; 
  float no_edge_weight, with_edge_weight;  
  set<RandomVariable*> parents_child, child_grandchild, all;  
  RAND rndm_nmbr(0);
  bool add_edge;

  //////////////////////////////////////////////////////////////////////////
  // Choose if the edges should be added according to edge_heuristic
  //////////////////////////////////////////////////////////////////////////
  switch (edge_heuristic) {

    case ALL_EDGES:
      add_edge = true;
      break;

    case RANDOM_EDGES:
      add_edge = rndm_nmbr.uniform( 1 ); 
      break;

    case LOCALLY_OPTIMAL_EDGES:
      for( crrnt_prnt = parent_set.begin(), end_prnt = parent_set.end();
           crrnt_prnt != end_prnt;
         ++crrnt_prnt ) {
        parents_child.insert( (*crrnt_prnt)->randomVariable );
      }
      parents_child.insert( child->randomVariable );

      child_grandchild.insert( child->randomVariable );
      child_grandchild.insert( grandchild->randomVariable );

      all = parents_child;
      all.insert( grandchild->randomVariable );

      parent_weight = MaxClique::computeWeight(parents_child);  
      child_weight  = MaxClique::computeWeight(child_grandchild);

      no_edge_weight = parent_weight + 
        log10(1+pow(10,child_weight-parent_weight));
      with_edge_weight = MaxClique::computeWeight(all); 
      if (with_edge_weight < no_edge_weight) { 
        add_edge = true;
      }
      else {
        add_edge = false;
      }
      break;
  
    default:
      assert(0);
      break;
  }

  //////////////////////////////////////////////////////////////////////////
  // Add the edges if the criteria was met 
  //////////////////////////////////////////////////////////////////////////
  if (add_edge) {

    //////////////////////////////////////////////////////////////////////////
    // Add the edge from each node in parent_set to grandchild 
    //////////////////////////////////////////////////////////////////////////
    for( crrnt_prnt = parent_set.begin(), end_prnt = parent_set.end();
         crrnt_prnt != end_prnt;
         ++crrnt_prnt ) {

      extra_edges.push_back( edge(*crrnt_prnt, grandchild) ); 
      (*crrnt_prnt)->marked = true;

      ////////////////////////////////////////////////////////////////////////
      // The current parent just gained an undirected neighbor, so recurse if 
      // the parent is deterministic or sparse. 
      ////////////////////////////////////////////////////////////////////////
      if  ((*crrnt_prnt)->parents.size() > 0)  { 
        if ( ((*crrnt_prnt)->randomVariable->deterministic()) || 
             (((*crrnt_prnt)->randomVariable->discrete) && 
              (((DiscreteRandomVariable*)
               ((*crrnt_prnt)->randomVariable))->sparse())) ) {
          addEdgesToNode( (*crrnt_prnt)->parents, *crrnt_prnt, grandchild, 
            edge_heuristic, extra_edges );        
        }
      }
    }
  }

}


/*-
 *-----------------------------------------------------------------------
 * BoundaryTriangulate::addEdges
 *   Add neighbors to RandomVariables according to a vector edge 
 *   structures.  Node that the edge class contains triangulateNode's, 
 *   but this proceedure adds the neighbors to the associated 
 *   RandomVariables. 
 * 
 * Preconditions:
 *   none
 *
 * Postconditions:
 *   The neighbors of the RandomVariable associated with each 
 *   triangulateNode in each edge will change. 
 *   
 * Side Effects:
 *   none
 *
 * Results:
 *   none
 *-----------------------------------------------------------------------
 */
void
BoundaryTriangulate::
addEdges(
  const vector<edge>& extra_edges
  )
{
  vector<edge>::const_iterator crrnt_edge; 
  vector<edge>::const_iterator end_edge; 

  for ( crrnt_edge = extra_edges.begin(),
        end_edge   = extra_edges.end();
        crrnt_edge != end_edge;
        ++crrnt_edge ) {

    (*crrnt_edge).first()->randomVariable->neighbors.insert(
      (*crrnt_edge).second()->randomVariable);
    (*crrnt_edge).second()->randomVariable->neighbors.insert(
      (*crrnt_edge).first()->randomVariable);
  }
}


/*-
 *-----------------------------------------------------------------------
 * BoundaryTriangulate::triangulateSimulatedAnnealing()
 *   Simulated Annealing is stochastic decent search through elimination 
 *   orderings.
 *  
 *   Simulated  annealing randomly permutes two nodes in the elimination 
 *   ordering, this move is accepted if the weight is improved and 
 *   accepted with some probability proportional to a "temperature" 
 *   parameter.  This function controls the change of the temperature 
 *   parameter, and the random permutations are handled by annealChain.   
 *
 *   The temperature begins high enough that most moves are accepted.  
 *   The temperature is gradually lowered so that fewer and fewer 
 *   ascent moves are accepted.  The temperature is lowered by the 
 *   schedule given by Aarts and Laarhoven:
 *  
 * crrnt_tmprtr=crrnt_tmprtr*(1+(crrnt_tmprtr*log(1+distance))/(3*std_dev))^-1
 *
 *   Annealing finishes when no moves are accepted in a chain or when the 
 *   stop criterion given by Varanelli is fulfilled: 
 *
 *     ((mean - best_this_weight) / std_dev) < stop_ratio
 *
 *   References:
 *   E.H.L Aarts and P.J.M. van Laarhoven, "A New Polynomial-Time Cooling 
 *   Schedule," Proc.  IEEE ICCAD-85, Santa Clara, CA, 206-208, 1985.
 *
 *   U. Kjærulff.  "Optimal decomposition of probabilistic networks by 
 *   simulated annealing." Statistics and Computing, 2(7-17), 1992.
 *
 *   J. Varanelli, "On the Acceleration of Simulated Annealing," PhD thesis,
 *   University of Virginia, Department of Computer Science
 *
 * Preconditions:
 *   Each variable in the set of nodes must have valid parent and 
 *   neighbor members and the parents/neighbors must only point to other 
 *   nodes in the set. 
 * 
 * Postconditions:
 *   The graph is triangulated to the lowest weight triangulation found
 *   at any point during the search.  The corresponding order is stored 
 *   in best_order.
 *
 * Side Effects:
 *   Neighbor members of each random variable can be changed.
 *
 * Results:
 *
 *-----------------------------------------------------------------------
 */
void
BoundaryTriangulate::
triangulateSimulatedAnnealing(
  const set<RandomVariable*>& nodes,
  const bool jtWeight,
  const set<RandomVariable*>& nodesRootMustContain,
  vector<MaxClique>&          best_cliques,
  vector<RandomVariable*>&    best_order,
  string&                     parameter_string 
  )
{
  ///////////////////////////////////////////////////////////////////////
  // Simulated annealing parameters:
  //
  // distance - controls amount of change in annealing temperature, a 
  //   smaller value gives a slower schedule 
  //     
  // stop_ratio - controls when annealing stops, a smaller value allows 
  //   annealing to run longer 
  //     
  // chain_length - number of permutations to try before lowering the 
  //   temperature.  A larger number brings annealing closer to a steady
  //   state distribution at each temperature, a smaller number anneals 
  //   faster.
  ///////////////////////////////////////////////////////////////////////
  const double   distance = 0.01;
  const double   stop_ratio = 0.0001;
  const unsigned chain_length = 1000;

  ///////////////////////////////////////////////////////////////////////
  // Record parameters in a comment string
  ///////////////////////////////////////////////////////////////////////
  enum {
    string_length = 16
  };
  char parameter[string_length];    

  sprintf( parameter, "%f", distance );
  assert( strlen(parameter)<string_length ); 
  parameter_string = "Distance:" + string(parameter); 
  sprintf( parameter, "%f", stop_ratio);
  assert( strlen(parameter)<string_length ); 
  parameter_string = parameter_string + "  Stop Ratio:" + string(parameter); 
  sprintf( parameter, "%d", chain_length);
  assert( strlen(parameter)<string_length ); 
  parameter_string = parameter_string + "  Chain Length:" + string(parameter); 
  
  ///////////////////////////////////////////////////////////////////////
  // Local Variables
  ///////////////////////////////////////////////////////////////////////
  vector<triangulateNode>          triangulate_nodes;
  vector<triangulateNode*>         crrnt_order;
  vector<triangulateNode*>         triangulate_best_order;
  list<vector<triangulateNode*> >  cliques;
  vector<MaxClique>                rv_cliques;
  vector<triangulateNode*>         dummy_order;
  vector<triangulateNghbrPairType> orgnl_nghbrs;

  back_insert_iterator<vector<triangulateNode*> > bi_crrnt(crrnt_order);
  back_insert_iterator<vector<triangulateNode*> > 
    bi_best(triangulate_best_order);
  vector<triangulateNode>::iterator  crrnt_node;  
  vector<triangulateNode>::iterator  end_node;  
  vector<triangulateNode*>::iterator crrnt_np;  
  vector<triangulateNode*>::iterator end_np;  

  double     crrnt_tmprtr;        // Current annealing temperature

  double     best_graph_weight;   // Best overal graph weight
  double     best_this_weight;    // Best graph weight on most recent trial

  double     weight_sum = 0;      // Sum of weights (for variance calculation)
  double     weight_sqr_sum = 0;  // Sum of weights^2 (for variance calculation)
  double     mean, variance, std_dev;

  double     ratio = 0;           // Stop ratio for current run 

  unsigned   i;
  unsigned   moves_accepted;

  bool       exit;

  ///////////////////////////////////////////////////////////////////
  // Initialize data structures 
  ///////////////////////////////////////////////////////////////////
  infoMsg(IM::Tiny, "Annealing:\n");

  fillTriangulateNodeStructures( nodes, triangulate_nodes );

  saveCurrentNeighbors( triangulate_nodes, orgnl_nghbrs );  

  for( crrnt_node = triangulate_nodes.begin(), 
       end_node   = triangulate_nodes.end(); 
       crrnt_node != end_node; 
       ++crrnt_node ) { 

    crrnt_order.push_back( &(*crrnt_node) );    
  }

  ///////////////////////////////////////////////////////////////////
  // Begin with random elimination order  
  ///////////////////////////////////////////////////////////////////
  random_shuffle (crrnt_order.begin(), crrnt_order.end());
  copy (crrnt_order.begin(), crrnt_order.end(), bi_best);

  fillInComputation( crrnt_order );
  maximumCardinalitySearch( triangulate_nodes, cliques, dummy_order, false );

  listVectorCliquetoVectorSetClique( cliques, best_cliques );
  best_graph_weight = graphWeight(best_cliques,jtWeight,nodesRootMustContain);
  weight_sum = best_graph_weight;
  weight_sqr_sum = best_graph_weight*best_graph_weight;

  infoMsg(IM::Tiny, "  Starting weight:  %f\n", best_graph_weight); 

  ///////////////////////////////////////////////////////////////////
  // Run a markov chain at high temperature to get statistics 
  // for starting temperature. 
  ///////////////////////////////////////////////////////////////////

  moves_accepted = annealChain(
    triangulate_nodes, 
    jtWeight,
    nodesRootMustContain,
    crrnt_order,
    triangulate_best_order,
    best_graph_weight,
    best_this_weight,
    HUGE_VAL,
    chain_length,
    weight_sum,         
    weight_sqr_sum,
    orgnl_nghbrs);         

  ++moves_accepted; 
  mean = weight_sum/moves_accepted; 
  variance = (weight_sqr_sum-(weight_sum*weight_sum)/(double)moves_accepted)/
              (double)moves_accepted;
  std_dev = sqrt(variance);

  crrnt_tmprtr = 1.0+std_dev;

  infoMsg(IM::Tiny, "  Initial Temperature: %f\n", crrnt_tmprtr);
  infoMsg(IM::Tiny, "  Distance:   %f\n", distance );
  infoMsg(IM::Tiny, "  Stop Ratio: %f\n", stop_ratio);
 
  ////////////////////////////////////////////////////////////
  // Begin with best order 
  ////////////////////////////////////////////////////////////
  crrnt_order.clear();
  copy( triangulate_best_order.begin(), triangulate_best_order.end(), 
    bi_crrnt );

  ////////////////////////////////////////////////////////////
  // Loop until stop condition found 
  ////////////////////////////////////////////////////////////
  exit = false;
  i = 0;
  do {
    ////////////////////////////////////////////////////////////////
    // Run a markov chain at current temperature 
    ////////////////////////////////////////////////////////////////
    weight_sum = 0;
    weight_sqr_sum = 0;

    moves_accepted = annealChain(
		       triangulate_nodes, 
		       jtWeight,
		       nodesRootMustContain,
                       crrnt_order,
                       triangulate_best_order,
                       best_graph_weight,
                       best_this_weight,
                       crrnt_tmprtr,
                       chain_length,
                       weight_sum,         
                       weight_sqr_sum,
		       orgnl_nghbrs);         

    ////////////////////////////////////////////////////////////////
    // Lower the temperature and calculate current stop ratio 
    ////////////////////////////////////////////////////////////////
    if (moves_accepted > 0) {

      ++moves_accepted; 
      mean = weight_sum/moves_accepted; 
      variance = ( weight_sqr_sum - (weight_sum*weight_sum)/
                   (double)moves_accepted ) / (double)moves_accepted;

      std_dev = sqrt(variance);

      crrnt_tmprtr = crrnt_tmprtr /
        (1 + (crrnt_tmprtr*log(1+distance))/(3*std_dev));  

      ratio = (mean - best_this_weight) / std_dev; 
      ++i;
    }
    ////////////////////////////////////////////////////////////////
    // Exit prematurely if no moves in the chain were accepted 
    ////////////////////////////////////////////////////////////////
    else { 
      exit = true; 
    }

    ////////////////////////////////////////////////////////////////
    // Exit prematurely if the timer has expired 
    ////////////////////////////////////////////////////////////////
    if (timer && timer->Expired()) {
      exit = true;
      infoMsg(IM::Tiny, "Exiting Annealing Because Time Has Expired: %d\n", i); 
    }  

  } while ((! exit) && (ratio > stop_ratio));

  ////////////////////////////////////////////////////////////////
  // Exit with the best triangulated graph found 
  ////////////////////////////////////////////////////////////////
  infoMsg(IM::Tiny, "  Exiting on iteration: %d\n", i); 

  best_order.clear();
  for( crrnt_np = triangulate_best_order.begin(), 
       end_np   = triangulate_best_order.end(); 
       crrnt_np != end_np; 
       ++crrnt_np ) {

    best_order.push_back( (*crrnt_np)->randomVariable ); 
  }

  best_cliques.clear();
  triangulateElimination( nodes, best_order, best_cliques);
  infoMsg(IM::Tiny, "  Annealing Best Weight: %f\n", best_graph_weight ); 
}


/*-
 *-----------------------------------------------------------------------
 * annealChain 
 *   Use annealing for a Markov Chain of specified length at a constant
 *   temperature.  Helper function to triangulateSimulatedAnnealing.
 *   The chain is generated by purmutation two of the nodes in the 
 *   elimination order, the move is accepted if it results in a lower 
 *   weight and the move is accepted with some probability if it results
 *   in a higher weight.  
 *
 * Preconditions:
 *   Each variable in the set of nodes must have valid parent and 
 *   neighbor members and the parents/neighbors must only point to other 
 *   nodes in the set. 
 *
 * Postconditions:
 *   crrnt_order       - the last order of the last accepted perumutation
 *   best_order        - the order giving the best weight found 
 *   best_graph_weight - the best graph weight found overall
 *   best_this_weight  - the best graph weight found in this chain
 *   weight_sum        - sum of all the weights of accepted triangulations 
 *   weight_sqr_sum    - sum of all the weight^2 of accepted triangulations 
 *
 * Side Effects:
 *   The partition is triangulated to the last triangulation in the 
 *   current order.  Neighbor members of each random variable can be 
 *   changed.
 *
 * Results:
 *   The number of permutations that were accepted 
 *
 *-----------------------------------------------------------------------
 */
unsigned
BoundaryTriangulate::
annealChain(
  vector<triangulateNode>&          nodes,
  const bool jtWeight,
  const set<RandomVariable*>& nodesRootMustContain,
  vector<triangulateNode*>&         crrnt_order,
  vector<triangulateNode*>&         best_order,
  double&                           best_graph_weight,
  double&                           best_this_weight,
  double                            temperature,
  unsigned                          iterations,
  double&                           weight_sum,         
  double&                           weight_sqr_sum,         
  vector<triangulateNghbrPairType>& orgnl_nghbrs
  )
{  
  vector<MaxClique>               rv_cliques;
  list<vector<triangulateNode*> > list_cliques;
  vector<triangulateNode*>        dummy_order;

  back_insert_iterator<vector<triangulateNode*> > bi_best(best_order);

  triangulateNode* tmp_node; 
  RAND             rndm_nmbr(0);
  double           crrnt_graph_weight;
  double           prvs_graph_weight;
  double           tmprtr_penalty;
  unsigned         first_index  = 0;    
  unsigned         second_index = 0;    
  unsigned         moves_accepted;
  unsigned         i;
  bool             accepted;
  
  crrnt_graph_weight = HUGE_VAL;
  prvs_graph_weight  = HUGE_VAL;
  best_this_weight   = HUGE_VAL;
  moves_accepted     = 0;

  ///////////////////////////////////////////////////////////////////////
  // Permute nodes 'iterations' times  
  ///////////////////////////////////////////////////////////////////////
  for(i=0; i<iterations; ++i) {

    ////////////////////////////////////////////////////////////////
    // Randomly swap two nodes 
    ////////////////////////////////////////////////////////////////
    if (crrnt_order.size() >= 2) {
      do { 
        first_index = rndm_nmbr.uniform( crrnt_order.size() - 1 );
        assert( first_index < crrnt_order.size() );
        second_index = rndm_nmbr.uniform( crrnt_order.size() - 1 );
        assert( second_index < crrnt_order.size() );
      } while (first_index == second_index);

      tmp_node = crrnt_order[first_index];
      crrnt_order[first_index]  = crrnt_order[second_index];
      crrnt_order[second_index] = tmp_node; 
    } 

    ////////////////////////////////////////////////////////////////
    // Calculate new graph weight 
    ////////////////////////////////////////////////////////////////
    restoreNeighbors(orgnl_nghbrs);
    fillInComputation( crrnt_order );
    maximumCardinalitySearch( nodes, list_cliques, dummy_order, false );
    listVectorCliquetoVectorSetClique( list_cliques, rv_cliques );
    crrnt_graph_weight = graphWeight(rv_cliques,jtWeight,nodesRootMustContain);

    ////////////////////////////////////////////////////////////////
    // Check if it is the best ordering so far 
    ////////////////////////////////////////////////////////////////
    if (crrnt_graph_weight < best_graph_weight) {
      best_graph_weight = crrnt_graph_weight;

      best_order.clear();
      copy(crrnt_order.begin(), crrnt_order.end(), bi_best);
    }

    ////////////////////////////////////////////////////////////////
    // If new weight is better always accept it.  If weight is 
    //  worse, accept it with a probability calculated from the 
    //  temperature. 
    ////////////////////////////////////////////////////////////////
    accepted = true;
    if ((crrnt_graph_weight-prvs_graph_weight) > 1e-8) {
      tmprtr_penalty = temperature * log(rndm_nmbr.drand48());

      if ( crrnt_graph_weight > (prvs_graph_weight-tmprtr_penalty)) {
        tmp_node = crrnt_order[first_index];
        crrnt_order[first_index]  = crrnt_order[second_index];
        crrnt_order[second_index] = tmp_node;
        accepted = false;
      }
    } 

    if (accepted == true) {
      prvs_graph_weight = crrnt_graph_weight;
      weight_sum += prvs_graph_weight;
      weight_sqr_sum += prvs_graph_weight*prvs_graph_weight;
      ++moves_accepted;
    }
  
    if (prvs_graph_weight < best_this_weight) {
      best_this_weight = prvs_graph_weight; 
    }
  }
 
  return(moves_accepted);
}


/*-
 *-----------------------------------------------------------------------
 * triangulateNode::triangulateNode (default constructor) 
 *
 * Preconditions:
 *   none
 *
 * Postconditions:
 *   none
 *
 * Side Effects:
 *   none
 *
 * Results:
 *
 *-----------------------------------------------------------------------
 */
BoundaryTriangulate::
triangulateNode::
triangulateNode(
  void 
  )
  : randomVariable( NULL ), 
    nodeList( NULL ),
    cardinality( 0 ),
    position( 0 ),
    eliminated( false ),
    marked( false ),
    previousNode( NULL ),
    nextNode( NULL )
{
}


/*-
 *-----------------------------------------------------------------------
 * triangulateNode::triangulateNode (constructor) 
 *
 * Preconditions:
 *   none
 *
 * Postconditions:
 *   none
 *
 * Side Effects:
 *   none
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
BoundaryTriangulate::
triangulateNode::
triangulateNode(
  RandomVariable* random_variable 
  )
  : randomVariable( random_variable ), 
    nodeList( NULL ),
    cardinality( 0 ),
    position( 0 ),
    eliminated( false ),
    marked( false ),
    previousNode( NULL ),
    nextNode( NULL )
{
}


/*-
 *-----------------------------------------------------------------------
 * triangulateNodeList::triangulateNodeList (constructor)
 *
 * Preconditions:
 *   none
 *
 * Postconditions:
 *   none
 *
 * Side Effects:
 *   none
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
BoundaryTriangulate::
triangulateNodeList::
triangulateNodeList()
  : last( NULL ),
    list_length( 0 )
{
}


/*-
 *-----------------------------------------------------------------------
 * triangulateNodeList::push_back
 *
 * Preconditions:
 *   none
 *
 * Postconditions:
 *   none
 *
 * Side Effects:
 *   none
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
BoundaryTriangulate::
triangulateNodeList::
push_back(
  triangulateNode* node 
  )
{
  node->previousNode = last; 
  node->nextNode     = NULL;
  node->nodeList     = this;

  if (last != NULL) { 
    last->nextNode = node;
  }

  last = node;
  list_length++;
}


/*-
 *-----------------------------------------------------------------------
 * triangulateNodeList::pop_back
 *
 * Preconditions:
 *   none
 *
 * Postconditions:
 *   The last node is removed from the list 
 *
 * Side Effects:
 *   none
 *
 * Results:
 *   A pointer to the (former) last node is returned  
 *-----------------------------------------------------------------------
 */
BoundaryTriangulate::triangulateNode*
BoundaryTriangulate::
triangulateNodeList::
pop_back()
{
  triangulateNode* deleted; 

  if (last->previousNode != NULL) {
    last->previousNode->nextNode = NULL; 
  }

  deleted = last;
  last    = last->previousNode;

  deleted->nextNode     = NULL;
  deleted->previousNode = NULL;
  deleted->nodeList     = NULL;
  
  list_length--;

  return(deleted);
}


/*-
 *-----------------------------------------------------------------------
 * triangulateNodeList::erase
 *
 * Preconditions:
 *   none
 *
 * Postconditions:
 *   The given node is removed from the list.  The memory for the node 
 *   is not-deallocated.
 *
 * Side Effects:
 *   none
 *
 * Results:
 *   none
 *-----------------------------------------------------------------------
 */
void
BoundaryTriangulate::
triangulateNodeList::
erase(
  triangulateNode* node 
  )
{
  assert(node->nodeList == this); 

  if (node == last) {
    last = last->previousNode;
  }

  if (node->previousNode != NULL) {
    node->previousNode->nextNode = node->nextNode;
  }

  if (node->nextNode != NULL) {
    node->nextNode->previousNode = node->previousNode;
  }

  node->nextNode     = NULL;
  node->previousNode = NULL;
  node->nodeList     = NULL;
  
  list_length--;
}


/*-
 *-----------------------------------------------------------------------
 * triangulateNode::operator[]
 *
 * Preconditions:
 *   none
 *
 * Postconditions:
 *   none
 *
 * Side Effects:
 *   none
 *
 * Results:
 *   A pointer to the i'th node in the list.  
 *
 * Complexity:
 *   This is a O(N) operation. 
 *-----------------------------------------------------------------------
 */
BoundaryTriangulate::triangulateNode* 
BoundaryTriangulate::
triangulateNodeList::
operator[] (
  unsigned i 
  )
{
  triangulateNode* node; 
  unsigned count; 

  assert(i < list_length);

  node = last;
  for( count=list_length-1; count>i; count-- ) {
    assert(node != NULL);
    node = node->previousNode;
  } 

  return(node);
}


/*-
 *-----------------------------------------------------------------------
 * fillTriangulateNodeStructures
 *
 * Preconditions:
 *   none 
 *
 * Postconditions:
 *   The graph structure given by a set of RandomVariable*'s is copied 
 *   into a set of triangulateNode's.  
 *
 * Side Effects:
 *   none 
 *
 * Results:
 *   none
 *
 * Complexity:
 *   O(N+E*log(N)) where E is the number of edges and N is the number 
 *   of nodes.
 *-----------------------------------------------------------------------
 */
void
BoundaryTriangulate::
fillTriangulateNodeStructures( 
  const set<RandomVariable*>& orgnl_nodes,
  vector<triangulateNode>&    new_nodes 
  )
{
  map<RandomVariable*, triangulateNode*> rv_to_mcs; 
  set<RandomVariable*>::iterator         crrnt_node;
  set<RandomVariable*>::iterator         end_node; 
  vector<triangulateNode>::iterator      crrnt_triangulate; 
  vector<triangulateNode>::iterator      end_triangulate; 
  triangulateNode                        new_node;

  ////////////////////////////////////////////////////////////////////////
  // Create a triangulateNode instance for each random variable 
  ////////////////////////////////////////////////////////////////////////
  new_nodes.clear();

  for (crrnt_node = orgnl_nodes.begin(),  
       end_node   = orgnl_nodes.end();  
       crrnt_node != end_node;
       ++crrnt_node) { 

    new_node.randomVariable = *crrnt_node;
    new_nodes.push_back( new_node );
  }

  ////////////////////////////////////////////////////////////////////////
  // Create map from the randomVariable* to its triangulateNode 
  ////////////////////////////////////////////////////////////////////////

  for (crrnt_triangulate = new_nodes.begin(),  
       end_triangulate   = new_nodes.end();  
       crrnt_triangulate != end_triangulate;
       ++crrnt_triangulate) {
 
    rv_to_mcs[(*crrnt_triangulate).randomVariable] = &(*crrnt_triangulate);
  }

  ////////////////////////////////////////////////////////////////////////
  // Create neighbor sets composed of triangulateNode's which match the 
  // original sets of RandomVariable*'s    
  ////////////////////////////////////////////////////////////////////////
  set<RandomVariable*>::iterator crrnt_nghbr;
  set<RandomVariable*>::iterator end_nghbr; 

  for (crrnt_node = orgnl_nodes.begin(),  
       end_node   = orgnl_nodes.end();  
       crrnt_node != end_node;
       ++crrnt_node) {

    for (crrnt_nghbr = (*crrnt_node)->neighbors.begin(),  
         end_nghbr   = (*crrnt_node)->neighbors.end();  
         crrnt_nghbr != end_nghbr;
         ++crrnt_nghbr) {

      rv_to_mcs[*crrnt_node]->neighbors.push_back( rv_to_mcs[*crrnt_nghbr] ); 
    }  
  }
}


/*-
 *-----------------------------------------------------------------------
 * maximumCardinalitySearch
 *   
 *   Calculated a perfect ordering on the graph, if it exists.  If the 
 *   order is perfect, the maximal cliques of the graph are determined.
 *   The procedure can be used as a first step in an O(N+E) chordality 
 *   test.  The procedure can also be used as a heuristic search for an 
 *   optimal elimination order.  If the given graph is not triangulated 
 *   the cliques correspond to the maximal cliques that will occur if 
 *   elimination is run on the graph using the given order.  
 *
 *   The function calculated the order from back to front.  At each step 
 *   it chooses the node with the largest number of previously ordered 
 *   neighbors. 
 *
 *   If randomize_order = true, when there is a tie for the node with the 
 *   largest number of previously ordered neighbors the tie is broken 
 *   randomly.  In this case, the function runs in O(N^2). 
 *
 *   If randomize_order = false, ties are broken in a deterministic
 *   manner and the function runs in O(N+E).  
 *
 *   Maximum cardinality search was originally described in:
 *     R. E. Tarjan and M. Yannakakis.  "Simple linear time algorithm to 
 *     test chordality of graphs, test acyclicity of hypergraphs, and 
 *     selectively reduce acyclic hypergraphs."  SIAM J. Comput., 
 *     13:566--579, 1984.  
 *
 * Preconditions:
 *   Each variable in the set of nodes must have valid parent and 
 *   neighbor members and the parents/neighbors must only point to other 
 *   nodes in the set. 
 *
 * Postconditions:
 *   order is cleared and replaced with an order based on the give nodes. 
 *  
 *   cliques is cleared and replaced with a new list based on the given 
 *   nodes.   
 *
 * Side Effects:
 *   none 
 *
 * Results:
 *   none 
 *
 *-----------------------------------------------------------------------
 */
void 
BoundaryTriangulate::
maximumCardinalitySearch( 
  vector<triangulateNode>&         nodes, 
  list<vector<triangulateNode*> >& cliques,
  vector<triangulateNode*>&        order,
  bool                             randomize_order 
  )
{
  vector<triangulateNodeList> card_set;

  vector<triangulateNode*> empty_clique; 
  list<vector<triangulateNode*> >::iterator prvs_clique;
  list<vector<triangulateNode*> >::iterator crrnt_clique;
  list<vector<triangulateNode*> >::iterator tmp_clique;

  vector<triangulateNode>::iterator  crrnt_node; 
  vector<triangulateNode>::iterator  end_node; 
  vector<triangulateNode*>::iterator crrnt_nghbr; 
  vector<triangulateNode*>::iterator end_nghbr; 

  unsigned         max_cardinality;    // largest cardinality found so far
  int              index;              // index of selected variable 
  triangulateNode* elmnt_node;         // node selected for elimination
  unsigned         nmbr_nodes_ordered; // count of number of nodes ordered
  RAND             rndm_nmbr(0);       // random number class

  ////////////////////////////////////////////////////////////////////
  // Put all nodes in cardinality set 0 
  ////////////////////////////////////////////////////////////////////
  card_set.resize( nodes.size() );

  for (crrnt_node = nodes.begin(), 
       end_node = nodes.end(); 
       crrnt_node != nodes.end();
       ++crrnt_node ) {

    (*crrnt_node).cardinality = 0;
    (*crrnt_node).eliminated  = false;
    card_set[0].push_back( &(*crrnt_node) ); 
  }

  max_cardinality = 0;
  
  ////////////////////////////////////////////////////////////////////
  // Set up the current clique and previous clique.
  ////////////////////////////////////////////////////////////////////
  cliques.clear();
  empty_clique.reserve( nodes.size() );
 
  cliques.push_back(empty_clique);
  prvs_clique = cliques.end();
  --prvs_clique;
  cliques.push_back(empty_clique);
  crrnt_clique = cliques.end();
  --crrnt_clique;

  ////////////////////////////////////////////////////////////////////
  // Iterate through all of the nodes 
  ////////////////////////////////////////////////////////////////////
  order.clear();

  for ( nmbr_nodes_ordered = 0;
        nmbr_nodes_ordered < nodes.size();
        nmbr_nodes_ordered++ ) {

    ////////////////////////////////////////////////////////////////////
    // Find unordered node with the largest number of ordered neighbors
    ////////////////////////////////////////////////////////////////////
    if (! randomize_order) {

      elmnt_node = card_set[max_cardinality].pop_back();
    }
    else { 
      ////////////////////////////////////////////////////////////////
      // Choose randomly among all nodes tied for maximum cardinality
      ////////////////////////////////////////////////////////////////
      index = rndm_nmbr.uniform( card_set[max_cardinality].size() - 1 );

      ////////////////////////////////////////////////////////////////
      // Get the pointer referred to by the index, be warned that 
      //  there is an O(N) operation hidden in the second [] 
      ////////////////////////////////////////////////////////////////
      elmnt_node = card_set[max_cardinality][index];  

      card_set[max_cardinality].erase(elmnt_node);
    }

    assert(elmnt_node != NULL);
    order.push_back(elmnt_node); 
    elmnt_node->eliminated = true;

    //////////////////////////////////////////////////////////////////
    // a) Move the node's ordered neighbors into a higher cardinality 
    //    set 
    // b) Remember the node's unordered neighbors as a potential max
    //    clique 
    //////////////////////////////////////////////////////////////////
    (*crrnt_clique).push_back( elmnt_node );
   
    for (crrnt_nghbr = elmnt_node->neighbors.begin(), 
         end_nghbr   = elmnt_node->neighbors.end(); 
         crrnt_nghbr != end_nghbr;
         ++crrnt_nghbr) {

      if ((*crrnt_nghbr)->eliminated == false) {

        card_set[(*crrnt_nghbr)->cardinality].erase( *crrnt_nghbr );
        (*crrnt_nghbr)->cardinality++; 
        card_set[(*crrnt_nghbr)->cardinality].push_back( *crrnt_nghbr );
      }
      else {

        (*crrnt_clique).push_back( *crrnt_nghbr );   
      }
    }

    ////////////////////////////////////////////////////////////////////
    // If the previous node had more or the same number of unordered 
    //   neighbors prvs_clique is a maximal clique.  Set prvs_clique as
    //   the crrnt_clique and start a new current clique.
    ////////////////////////////////////////////////////////////////////
    if ((*prvs_clique).size() >= (*crrnt_clique).size()) {

      prvs_clique = crrnt_clique;
      cliques.push_back(empty_clique);
      crrnt_clique = cliques.end();
      --crrnt_clique; 
    }
    ////////////////////////////////////////////////////////////////////
    // Maximumal clique is not detected, so swap prvs_clique and  
    // crrnt_clique and clear the new crrnt_clique.  This sets the 
    // previous clique to be the current clique and clears next 
    // iteration's current clique without deallocating or allocating any
    // memory. 
    ////////////////////////////////////////////////////////////////////
    else {
      cliques.splice(prvs_clique, cliques, crrnt_clique );

      tmp_clique   = prvs_clique;
      prvs_clique  = crrnt_clique;
      crrnt_clique = tmp_clique;
      (*crrnt_clique).clear();
    }

    //////////////////////////////////////////////////////////////////
    // Set the maximum cardinality  
    //////////////////////////////////////////////////////////////////
    ++max_cardinality;
    while( card_set[max_cardinality].size() == 0 ) {
      --max_cardinality;
    }

  }

  //////////////////////////////////////////////////////////////////
  // Remove the current clique (which will just be the last node 
  // and not a maximal clique)
  //////////////////////////////////////////////////////////////////
  cliques.erase( crrnt_clique ); 

  //////////////////////////////////////////////////////////////////
  // Put the order in forward elimination order 
  //////////////////////////////////////////////////////////////////
  reverse( order.begin(), order.end() );
}  


/*-
 *-----------------------------------------------------------------------
 * fillInComputation
 *  
 *   Triangulates a graph according to a given elimination order in 
 *   O(N+E') time, where N is the number of nodes and E' is the number 
 *   of original graph edges plus the number of fill in edges. 
 *
 *   This algorithm is given in:
 *     R. E. Tarjan and M. Yannakakis.  "Simple linear time algorithm to 
 *     test chordality of graphs, test acyclicity of hypergraphs, and 
 *     selectively reduce acyclic hypergraphs."  SIAM J. Comput., 
 *     13:566--579, 1984.  
 *
 * Preconditions:
 *   Each variable in the set of nodes must have valid parent and 
 *   neighbor members and the parents/neighbors must only point to other 
 *   nodes in the set. 
 *
 * Postconditions:
 *   The graph is triangulated according to the elimination order 
 *   ordered_nodes. 
 *
 * Side Effects:
 *   Neighbor members of each random variable can be changed.
 *
 * Results:
 *   none 
 *
 *-----------------------------------------------------------------------
 */
void
BoundaryTriangulate::
fillInComputation( 
  vector<triangulateNode*>& ordered_nodes 
  )
{
  const unsigned invalid_index = ~0;

  typedef pair<triangulateNode*, triangulateNode*> edge; 
  
  vector<triangulateNode*> follower; 
  vector<unsigned>         index;
  vector<edge>             edges;  
 
  vector<triangulateNode*>::iterator crrnt_node;
  vector<triangulateNode*>::iterator end_node;
  vector<triangulateNode*>::iterator crrnt_nghbr;
  vector<triangulateNode*>::iterator end_nghbr;
  unsigned i;
 
  triangulateNode* new_nghbr;  

  ////////////////////////////////////////////////////////////////////
  // Set up follower and index arrays 
  ////////////////////////////////////////////////////////////////////
  follower.resize(ordered_nodes.size());
  index.resize(ordered_nodes.size(), invalid_index);

  ////////////////////////////////////////////////////////////////////
  // Record each node's position in the order into the triangulateNode 
  // structure so that positions can be determined from the pointers 
  ////////////////////////////////////////////////////////////////////
  for (unsigned j = 0; j<ordered_nodes.size(); j++) {

    ordered_nodes[j]->position = j; 
  }

  ////////////////////////////////////////////////////////////////////
  // Iterate through all nodes 
  ////////////////////////////////////////////////////////////////////
  for (crrnt_node = ordered_nodes.begin(), 
       end_node   = ordered_nodes.end(),
       i = 0;
       crrnt_node != end_node;
       ++i, ++crrnt_node ) {

    follower[i] = *crrnt_node; 
    index[i]    = i;

    ////////////////////////////////////////////////////////////////////
    // Find edges {crrnt_node, new_nghbr} with (crrnt_node<new_nghbr) 
    // such that there is a vertex crrnt_nghbr with with 
    // {crrnt_nghbr, new_nghbr} is a graph edge and 
    // follower^i(crrnt_nghbr)=crrnt_node  
    ////////////////////////////////////////////////////////////////////
    for (crrnt_nghbr = (*crrnt_node)->neighbors.begin(), 
         end_nghbr   = (*crrnt_node)->neighbors.end();
         crrnt_nghbr != end_nghbr; 
         ++crrnt_nghbr ) {

      if ((*crrnt_nghbr)->position < i) {

        new_nghbr = *crrnt_nghbr; 
        assert( index[new_nghbr->position] != invalid_index );

        while ( index[new_nghbr->position] < i ) { 

          index[new_nghbr->position] = i;
          edges.push_back( edge(*crrnt_node, new_nghbr) );
          new_nghbr = follower[new_nghbr->position]; 
        }

        if (follower[new_nghbr->position] == new_nghbr) {

          follower[new_nghbr->position] = *crrnt_node;
        } 
      } 
    } 
  }   

  ////////////////////////////////////////////////////////////////////
  // Remove the previous edge sets 
  ////////////////////////////////////////////////////////////////////
  vector<edge>::iterator crrnt_edge;  
  vector<edge>::iterator end_edge;  
 
  for (crrnt_node = ordered_nodes.begin(), 
       end_node   = ordered_nodes.end();
       crrnt_node != end_node;
       ++crrnt_node ) {

    (*crrnt_node)->neighbors.clear();  
  }
 
  ////////////////////////////////////////////////////////////////////
  // Add the new edge sets to the graph 
  ////////////////////////////////////////////////////////////////////
  for (crrnt_edge = edges.begin(), 
       end_edge   = edges.end();
       crrnt_edge != end_edge;
       ++crrnt_edge ) {

    (*crrnt_edge).first->neighbors.push_back( (*crrnt_edge).second ); 
    (*crrnt_edge).second->neighbors.push_back( (*crrnt_edge).first ); 
  }
}
 

/*-
 *-----------------------------------------------------------------------
 * testZeroFillIn
 *  
 *  Determines if an elimination order is zero fill in.  This algorithm
 *  runs in O(N+E) time.  Combined with maximum cardinality search this
 *  function provides the second step of a chordality test. 
 *
 *  The algorithm was given in:
 *    R. E. Tarjan and M. Yannakakis.  "Simple linear time algorithm to 
 *    test chordality of graphs, test acyclicity of hypergraphs, and 
 *    selectively reduce acyclic hypergraphs."  SIAM J. Comput., 
 *    13:566--579, 1984.  
 *
 * Preconditions:
 *   Each variable in the set of nodes must have valid parent and 
 *   neighbor members and the parents/neighbors must only point to other 
 *   nodes in the set. 
 *
 * Postconditions:
 *   none
 *
 * Side Effects:
 *   none 
 *
 * Results:
 *   true if the graph is triangulated   
 *
 *-----------------------------------------------------------------------
 */
bool
BoundaryTriangulate::
testZeroFillIn( 
  vector<triangulateNode*>& ordered_nodes 
  )
{
  const unsigned invalid_index = ~0;

  vector<triangulateNode*> follower; 
  vector<unsigned> index;
  
  vector<triangulateNode*>::iterator crrnt_node;
  vector<triangulateNode*>::iterator end_node;
  vector<triangulateNode*>::iterator crrnt_nghbr;
  vector<triangulateNode*>::iterator end_nghbr;
  unsigned i;
  bool chordal = true;

  ////////////////////////////////////////////////////////////////////
  // Set up follower and index arrays 
  ////////////////////////////////////////////////////////////////////
  follower.resize(ordered_nodes.size());
  index.resize(ordered_nodes.size(), invalid_index);

  ////////////////////////////////////////////////////////////////////
  // Record each node's position in the order into the triangulateNode 
  // structure so that positions can be determined from the pointers 
  ////////////////////////////////////////////////////////////////////
  for (unsigned j = 0; j<ordered_nodes.size(); j++) {

    ordered_nodes[j]->position = j; 
  }

  ////////////////////////////////////////////////////////////////////
  // Iterate through all nodes 
  ////////////////////////////////////////////////////////////////////
  for (crrnt_node = ordered_nodes.begin(), 
       end_node   = ordered_nodes.end(),
       i = 0;
       (crrnt_node != end_node) && (chordal == true); 
       ++i, ++crrnt_node ) {

    assert( (*crrnt_node)->position == i ); 
    follower[i] = *crrnt_node; 
    index[i]    = i;

    ////////////////////////////////////////////////////////////////////
    // Calculate the followers and indexes of ordered neighbors 
    ////////////////////////////////////////////////////////////////////
    for (crrnt_nghbr = (*crrnt_node)->neighbors.begin(), 
         end_nghbr   = (*crrnt_node)->neighbors.end();
         crrnt_nghbr != end_nghbr; 
         ++crrnt_nghbr ) {

      assert( (*crrnt_nghbr)->position != invalid_index ); 

      if ((*crrnt_nghbr)->position < i) {

        index[(*crrnt_nghbr)->position] = i;

        if ( follower[(*crrnt_nghbr)->position] == *crrnt_nghbr ) {
          follower[(*crrnt_nghbr)->position] = *crrnt_node;
        }           
      }
    }

    ////////////////////////////////////////////////////////////////////
    // Graph is not triangulated if any ordered neighbors did not have
    // their indexes set. 
    ////////////////////////////////////////////////////////////////////
    for (crrnt_nghbr = (*crrnt_node)->neighbors.begin(), 
         end_nghbr   = (*crrnt_node)->neighbors.end();
         crrnt_nghbr != end_nghbr; 
         ++crrnt_nghbr ) {

      if ((*crrnt_nghbr)->position < i) { 

        if ( index[(follower[(*crrnt_nghbr)->position])->position] < i ) {
          chordal = false;
        } 
      }
    }
  } 

  return(chordal);
}


/*-
 *-----------------------------------------------------------------------
 * triangulateMaximumCardinalitySearch
 *
 *   The graph is triangulated with an elimination order as determined by 
 *   maximum cardinality search,  The maximal cliques are stored in RIP 
 *   order.  If the input graph is already triangulated it will not be
 *   changed by this function. 
 *
 * Preconditions:
 *   Each variable in the set of nodes must have valid parent and 
 *   neighbor members and the parents/neighbors must only point to other 
 *   nodes in the set. 
 *
 * Postconditions:
 *   The graph is triangulated, the maximal cliques of the triangulated
 *   graph are stored in the parameter 'cliques', and the elimination order 
 *   used is stored in the parameter 'order'.
 *
 * Side Effects:
 *   none 
 *
 * Results:
 *   none 
 *
 * Complexity:
 *   This procedure runs O((N^2)log(N)), the cost of converting the 
 *   list of vector cliques to the vector of set cliques. 
 * 
 * (Also see triangulateMCSIfNotTriangulated and maximumCardinalitySearch) 
 *-----------------------------------------------------------------------
 */
void 
BoundaryTriangulate::
triangulateMaximumCardinalitySearch( 
  const set<RandomVariable*>& nodes,
  vector<MaxClique>&          cliques,
  vector<RandomVariable*>&    order
  )
{
  vector<triangulateNode>           triangulate_nodes; 
  list<vector<triangulateNode*> >   list_cliques;
  vector<triangulateNode*>          triangulate_order; 
  vector<triangulateNode*>          dummy_order; 
  vector<triangulateNode>::iterator crrnt_node; 
  vector<triangulateNode>::iterator end_node; 

  fillTriangulateNodeStructures( nodes, triangulate_nodes );

  ////////////////////////////////////////////////////////////////////////// 
  // Find a triangulation using maximum cardinality search 
  ////////////////////////////////////////////////////////////////////////// 
  maximumCardinalitySearch( triangulate_nodes, list_cliques, triangulate_order, 
    true);
  fillInComputation( triangulate_order );

  ////////////////////////////////////////////////////////////////////////// 
  // Now use maximum cardinality search to get the cliques of the 
  // triangulated graph 
  ////////////////////////////////////////////////////////////////////////// 
  maximumCardinalitySearch( triangulate_nodes, list_cliques, dummy_order, 
    false);

  ////////////////////////////////////////////////////////////////////////// 
  // Convert to MaxCliques and triangulate the RandomVariable structures 
  ////////////////////////////////////////////////////////////////////////// 
  listVectorCliquetoVectorSetClique( list_cliques, cliques );
  fillAccordingToCliques( cliques );

  ////////////////////////////////////////////////////////////////////////// 
  // Store the order 
  ////////////////////////////////////////////////////////////////////////// 
  for ( crrnt_node = triangulate_nodes.begin(), 
        end_node   = triangulate_nodes.end();  
        crrnt_node != end_node ;
        ++crrnt_node) {
    order.push_back( (*crrnt_node).randomVariable );
  }
}


/*-
 *-----------------------------------------------------------------------
 * chordalityTest 
 *   
 * Runs in O(N+E*log(N)), the cost of converting from random variables to 
 * triangulateNodes
 *
 * Preconditions:
 *   Each variable in the set of nodes must have valid parent and 
 *   neighbor members and the parents/neighbors must only point to other 
 *   nodes in the set. 
 *
 * Postconditions:
 *   none 
 *
 * Side Effects:
 *   none 
 *
 * Results:
 *   true if graph is chordal, false if not 
 *
 * Complexity:
 *   O(N+E*log(N)), the cost of converting the RandomVariables to 
 *   triangulateNodes.  
 *-----------------------------------------------------------------------
 */
bool
BoundaryTriangulate::
chordalityTest( 
  const set<RandomVariable*>& nodes
  )
{
  vector<triangulateNode>         triangulate_nodes; 
  list<vector<triangulateNode*> > list_cliques;
  vector<triangulateNode*>        order; 
  bool                            chordal;

  fillTriangulateNodeStructures( nodes, triangulate_nodes );
  maximumCardinalitySearch( triangulate_nodes, list_cliques, order, false);
  chordal = testZeroFillIn( order );

  return(chordal);
}


/*-
 *-----------------------------------------------------------------------
 * getCliques 
 *   
 * Preconditions:
 *   Each variable in the set of nodes must have valid parent and 
 *   neighbor members and the parents/neighbors must only point to other 
 *   nodes in the set. 
 *
 * Postconditions:
 *   If the graph is chordal, the maximal cliques are stored in the 
 *   cliques variable in RIP order.  If the graph is not chordal, the 
 *   cliques are not modified. 
 *
 * Side Effects:
 *   none 
 *
 * Results:
 *   true if graph is chordal, false if not 
 *
 * Complexity:
 *   This procedure runs O((N^2)log(N)), the cost of converting the 
 *   list of vector cliques to the vector of set cliques. 
 *-----------------------------------------------------------------------
 */
bool
BoundaryTriangulate::
getCliques( 
  const set<RandomVariable*>& nodes,
  vector<MaxClique>&          cliques
  )
{
  vector<triangulateNode>         triangulate_nodes; 
  list<vector<triangulateNode*> > list_cliques;
  vector<triangulateNode*>        order; 
  bool                            chordal;

  fillTriangulateNodeStructures( nodes, triangulate_nodes );
  maximumCardinalitySearch( triangulate_nodes, list_cliques, order, false);
  chordal = testZeroFillIn( order );

  if (chordal) {
    listVectorCliquetoVectorSetClique( list_cliques, cliques );
  }

  return(chordal);
}


/*-
 *-----------------------------------------------------------------------
 * triangulateMCSIfNotTriangulated
 *
 *   The graph is triangulated with an elimination order as determined by 
 *   maximum cardinality search,  The maximal cliques are stored in RIP 
 *   order.  If the input graph is already triangulated it will not be
 *   changed by this function. 
 *
 *   This procedure performs almost the same function as 
 *   triangulateMaximumCardinalitySearch.  This version has a return 
 *   value telling the user if the original graph was triangulated, does 
 *   not return the elimination order, will be a bit faster if the 
 *   input graph is already triangualted.  This version will always 
 *   return the same triangulation, whereas 
 *   triangulateMaximumCardinalitySearch will randomize if possible. 
 *
 * Preconditions:
 *   Each variable in the set of nodes must have valid parent and 
 *   neighbor members and the parents/neighbors must only point to other 
 *   nodes in the set. 
 *
 * Postconditions:
 *   The graph is triangulated, and the maximal cliques of the triangulated
 *   graph are stored in the parameter 'cliques'.
 *
 * Side Effects:
 *   none 
 *
 * Results:
 *   true if input graph is chordal, false if not 
 *
 * Complexity:
 *   This procedure runs O((N^2)log(N)), the cost of converting the 
 *   list of vector cliques to the vector of set cliques. 
 *-----------------------------------------------------------------------
 */
bool
BoundaryTriangulate::
triangulateMCSIfNotTriangulated( 
  const set<RandomVariable*>& nodes,
  vector<MaxClique>&          cliques
  )
{
  vector<triangulateNode>         triangulate_nodes; 
  list<vector<triangulateNode*> > list_cliques;
  vector<triangulateNode*>        order; 
  bool                            chordal;

  fillTriangulateNodeStructures( nodes, triangulate_nodes );
  maximumCardinalitySearch( triangulate_nodes, list_cliques, order, false);
  chordal = testZeroFillIn( order );

  if (!chordal) {
    ////////////////////////////////////////////////////////////////////////// 
    // Triangulate the graph according to the order found using MCS 
    ////////////////////////////////////////////////////////////////////////// 
    fillInComputation( order );

    ////////////////////////////////////////////////////////////////////////// 
    // Use maximum cardinality search to get the cliques of the 
    // triangulated graph 
    ////////////////////////////////////////////////////////////////////////// 
    maximumCardinalitySearch( triangulate_nodes, list_cliques, order, 
      false);
  }

  listVectorCliquetoVectorSetClique( list_cliques, cliques );

  ////////////////////////////////////////////////////////////////////////// 
  // Convert to MaxCliques and triangulate the RandomVariable structures 
  ////////////////////////////////////////////////////////////////////////// 
  listVectorCliquetoVectorSetClique( list_cliques, cliques );
  fillAccordingToCliques( cliques );

  return(chordal);
}


/*-
 *-----------------------------------------------------------------------
 * triangulateCompletePartition 
 *   adds eges between every node in the graph
 *
 * Preconditions:
 *   Each variable in the set of nodes must have valid parent and 
 *   neighbor members and the parents/neighbors must only point to other 
 *   nodes in the set. 
 *
 * Postconditions:
 *   Triangulates by adding edge between every node in the graph.  
 *
 * Side Effects:
 *   There are edges between every node in the graph.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
BoundaryTriangulate::
triangulateCompletePartition( 
  const set<RandomVariable*>& nodes,
  vector<MaxClique>&          cliques
  )
{
  ////////////////////////////////////////////////////////////////////
  // Add an edge between every combination of nodes 
  ////////////////////////////////////////////////////////////////////
  MaxClique::makeComplete(nodes);

  cliques.clear();
  cliques.push_back(MaxClique(nodes));

  return;
}



/*-
 *-----------------------------------------------------------------------
 * triangulateFrontier
 *
 *   Triangulate using the Frontier algorithm, a procedure that works
 *   entirely on the directed graph. The Frontier algorithm is defined
 *   in:
 *
 *       "A forward-backward algorithm for inference in Bayesian
 *       networks and an empirical comparison with HMMs", G. Zweig,
 *       Master's Thesis, Dept. of Computer Science, U.C. Berkeley,
 *       May 9th, 1996.
 *
 *   We don't use the Frontier alg. for inference (as it was
 *   originally defined, it was an algorithm to specify the order of
 *   variables to marginalize a big joint probability
 *   distribution). Here, rather, the essentials of the Frontier
 *   algorithm are extracted just to perform a graph triangulation for
 *   us that can be evaluated and compared with the other graph
 *   triangulation heuristics implemented in GMTK.
 *  
 *   Note, that Frontier requires a topological ordering of the nodes.
 *   In this implementation, we first compute a random topological
 *   ordering. This allows this routine to be called many times, each
 *   time it will produce a different triangulation.
 *
 *   Note that the graph doesn't need to be moralized here, as
 *   Frontier does that for us when it selects cliques, but it won't
 *   change things or hurt if it already is.  Also, note that Frontier
 *   sometimes (but rarely) will not triangulate the graph with
 *   respect to the cumpulsory interface completion edges that have at
 *   this point been added to the graph. If Frontier misses those
 *   edges, then we do a quick MCS pass to fix this up. In practice,
 *   however, this does not happen very often, and also it depends on
 *   the current boundary. If we spent time finding a good boundary
 *   then the cases when Frontier does this are likely to be quite
 *   poor triangulations.
 *   
 *
 * Preconditions:
 *   Each variable in the set of nodes must have valid parent and 
 *   neighbor members and the parents/neighbors must only point to other 
 *   nodes in the set. 
 *
 * Postconditions:
 *   Triangulates by running the Frontier algorithm.
 *
 * Side Effects:
 *   There are edges between every node in each clique in the graph.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
BoundaryTriangulate::
triangulateFrontier(const set<RandomVariable*>& nodes,
		    vector<MaxClique>&          cliques
		    )
{
  cliques.clear();
  vector <RandomVariable*> sortedNodes;

  GraphicalModel::topologicalSortRandom(nodes,nodes,sortedNodes);
  if (message(High)) {
    infoMsg(High,"Frontier: Sorted Nodes:");
    printRVSet(stdout,sortedNodes);
  }

  set <RandomVariable*> frontier;
  set <RandomVariable*> cumulativeFrontier;
  for (unsigned i=0;i<sortedNodes.size();i++) {
    RandomVariable*nrv = sortedNodes[i];
    if (message(High)) {
      infoMsg(High,"Frontier: Current nodes");
      printRVSet(stdout,frontier);
    }
    
    set <RandomVariable*>::iterator it;
    set <RandomVariable*>::iterator it_end = frontier.end();
    set <RandomVariable*> toRemove;
    for (it=frontier.begin();it!=it_end;it++) {
      RandomVariable*rv = (*it);
      bool childrenInFrontier = true;
      // Check if all of rv's children are in the
      // cumulative frontier.
      infoMsg(High,"Frontier: Node %s(%d) has %d children\n",
	      rv->name().c_str(),rv->frame(),rv->allPossibleChildren.size());

      for (unsigned c=0;c<rv->allPossibleChildren.size();c++) {
	RandomVariable* child = rv->allPossibleChildren[c];
	// only consider nodes within this partition set.
	if (nodes.find(child) == nodes.end()) {
	  infoMsg(High,"Frontier: Child %d %s(%d) not in nodes\n",
		  c,child->name().c_str(),child->frame());
	  continue;
	}
 	if (cumulativeFrontier.find(child) == cumulativeFrontier.end()) {
	  infoMsg(High,"Frontier: Child %d %s(%d) not in frontier\n",
		  c,child->name().c_str(),child->frame());
	  childrenInFrontier = false;
	  break;
	}
      }
      if (childrenInFrontier)
	toRemove.insert(rv);
    }
    if (toRemove.size() > 0) {
      // we have a clique
      if (message(High)) {
	infoMsg(High,"Frontier: Clique:");
	printRVSet(stdout,frontier);
      }
      MaxClique::makeComplete(frontier);
      set <RandomVariable*> res;
      if (message(High)) {
	infoMsg(High,"Frontier: Removing:");
	printRVSet(stdout,toRemove);
      }

      set_difference(frontier.begin(),frontier.end(),
		     toRemove.begin(),toRemove.end(),
		     inserter(res,res.end()));
      frontier = res;
    }
    infoMsg(High,"Frontier: adding node %s(%d)\n",
	    nrv->name().c_str(),nrv->frame());
    // insert new members
    frontier.insert(nrv);
    cumulativeFrontier.insert(nrv);
  }
  // last one is always a clique
  if (message(High)) {
    infoMsg(High,"Frontier: Clique:");
    printRVSet(stdout,frontier);
  }
  MaxClique::makeComplete(frontier);

  // Sometimes, but rarely, frontier might not triangulate the graph
  // since it doesn't know about the extra compulsory edges for the
  // completing the left and right interfaces of the current
  // partition.
  // 
  // In order to make sure, we run a MCS pass adding edges in case
  // Frontier didn't enclose right interface with a clique and it
  // isn't triangulated.  Note that if Frontier indeed triangulated
  // the graph (such that the extra forced completion of the
  // interfaces are included in cliques or the result is
  // triangulated), then this next step is guaranteed not to change
  // the graph at all.
  if (!triangulateMCSIfNotTriangulated(nodes,cliques))
    infoMsg(High,"Frontier: MCS used to fix up frontier triangulation\n");

  return;
}




/*-
 *-----------------------------------------------------------------------
 * listVectorCliquetoVectorSetClique
 *
 * Preconditions:
 *   none
 * 
 * Postconditions:
 *   The vector of MaxCliques matches the list of vectors of 
 *   triangulateNode pointers.  This is useful for converting a graph 
 *   defined by triangulateNode's into a graph defined by RandomVariable's
 *
 * Side Effects:
 *   none
 *
 * Results:
 *   none
 *
 * Complexity:
 *   O((N^2)log(N))
 *-----------------------------------------------------------------------
 */
void  
BoundaryTriangulate::
listVectorCliquetoVectorSetClique(
  const list<vector<triangulateNode*> >& lv_cliques,
  vector<MaxClique>&                     vs_cliques
  )
{
  vector<MaxClique>::iterator crrnt_vs_clique;

  list<vector<triangulateNode*> >::const_iterator crrnt_lv_clique; 
  list<vector<triangulateNode*> >::const_iterator end_lv_clique;

  vector<triangulateNode*>::const_iterator crrnt_node;
  vector<triangulateNode*>::const_iterator end_node;
 
  set<RandomVariable*> empty_RV_set;
  MaxClique empty_MaxClique( empty_RV_set );

  ////////////////////////////////////////////////////////////////////
  // Iterate through original data structure and construct the new 
  ////////////////////////////////////////////////////////////////////
  vs_cliques.clear();

  for( crrnt_lv_clique = lv_cliques.begin(),  
       end_lv_clique   = lv_cliques.end();  
       crrnt_lv_clique != end_lv_clique;   
       ++crrnt_lv_clique ) {

    vs_cliques.push_back(empty_MaxClique);
    crrnt_vs_clique = vs_cliques.end();
    --crrnt_vs_clique; 
     
    for( crrnt_node = (*crrnt_lv_clique).begin(),  
         end_node   = (*crrnt_lv_clique).end();  
         crrnt_node != end_node;   
         ++crrnt_node ) {

      (*crrnt_vs_clique).nodes.insert( (*crrnt_node)->randomVariable );   
    }
  }

}


/*-
 *-----------------------------------------------------------------------
 * triangulateExhaustiveSearch 
 *   Tries every possible combination of edges to determine find a 
 *   minimum weight triangulation.
 * 
 * Preconditions:
 *   Each variable in the set of nodes must have valid parent and 
 *   neighbor members and the parents/neighbors must only point to other 
 *   nodes in the set. 
 *
 * Postconditions:
 *   The minimum weight triangulation is stored in best_cliques.  Graph 
 *   is triangulated.
 *
 * Side Effects:
 *   Neighbor members of each random variable can be changed.  
 *
 * Results:
 *      none
 *-----------------------------------------------------------------------
 */
void
BoundaryTriangulate::
triangulateExhaustiveSearch( 
  const set<RandomVariable*>&  nodes,
  const bool                   jtWeight,
  const set<RandomVariable*>&  nodesRootMustContain,
  const vector<nghbrPairType>& orgnl_nghbrs,
  vector<MaxClique>&           best_cliques   
  )
{
  // define local constants used in this routine.
  const unsigned NONE = 0;
  const unsigned FILL_IN = 1;
  const unsigned EDGE = 2; 
  const unsigned UNUSED = 3;

  vector<vector<unsigned> > adjacency;  
  vector<MaxClique>         vector_cliques;

  unsigned nmbr_empty, nmbr_nodes;
  unsigned i, j;
  unsigned row, col;

  vector<triangulateNode>            triangulate_nodes; 

  list<vector<triangulateNode*> > cliques;
  list<vector<triangulateNode*> > best_list_cliques;

  list<vector<triangulateNode*> >::iterator crrnt_clique;
  list<vector<triangulateNode*> >::iterator end_clique;

  vector<triangulateNode*>           order; 
  vector<triangulateNode*>::iterator crrnt_node;
  vector<triangulateNode*>::iterator end_node;
  vector<triangulateNode*>::iterator nghbr; 

  ////////////////////////////////////////////////////////////////////
  // Reserve space in the containers 
  ////////////////////////////////////////////////////////////////////
  vector_cliques.reserve( nodes.size() );
  triangulate_nodes.reserve( nodes.size() ); 
  order.reserve( nodes.size() ); 

  ////////////////////////////////////////////////////////////////////
  // Create triagulateNode object from the RandomVariable set 
  ////////////////////////////////////////////////////////////////////
  fillTriangulateNodeStructures( nodes, triangulate_nodes );

  ////////////////////////////////////////////////////////////////////
  // Build the adjacency matrix 
  ////////////////////////////////////////////////////////////////////

  nmbr_empty = 0;
  nmbr_nodes = nodes.size();
  adjacency.resize( nmbr_nodes );

  for( i=0; i<nmbr_nodes; i++) {

    adjacency[i].resize( nmbr_nodes, UNUSED );

    for( j=i+1; j<nmbr_nodes; j++) {
  
      nghbr = find( triangulate_nodes[i].neighbors.begin(), 
        triangulate_nodes[i].neighbors.end(), &triangulate_nodes[j] );   
      if (nghbr == triangulate_nodes[i].neighbors.end()) {
         adjacency[i][j] = NONE;
         nmbr_empty++; 
      }
      else {
        adjacency[i][j] = EDGE;
      } 
    }
  } 
  
  ////////////////////////////////////////////////////////////////////
  // Display node key and initial matrix 
  ////////////////////////////////////////////////////////////////////
  infoMsg(IM::Tiny, "----------\n"); 
        
  for( i=0; i<nmbr_nodes; i++) {
    infoMsg(IM::Tiny, "[%d] %s(%d)\n", i, 
      triangulate_nodes[i].randomVariable->name().c_str(), 
      triangulate_nodes[i].randomVariable->timeIndex );
  } 
  infoMsg(IM::Tiny, "\n");

  infoMsg(IM::Tiny, "   ");
  for( i=0; i<nmbr_nodes; i++) {
    infoMsg(IM::Tiny, "[%2d]", i);
  }
  infoMsg(IM::Tiny, "\n");

  for( i=0; i<nmbr_nodes; i++) {
    infoMsg(IM::Tiny, "[%2d]", i );
    for( j=0; j<nmbr_nodes; j++) {
      infoMsg(IM::Tiny," %2d ", adjacency[i][j]);
    }
    infoMsg(IM::Tiny, "\n");
  }  

  ////////////////////////////////////////////////////////////////////
  // Begin Search 
  ////////////////////////////////////////////////////////////////////
  double best_weight = HUGE_VAL; 
  double weight; 
  bool done = false;
  bool chordal;
  double nmbr_cmbntns;
  unsigned crrnt_trial;

  if (nmbr_empty < 1024) {
    nmbr_cmbntns = pow(2.0, (double)nmbr_empty);
  }
  else {
    nmbr_cmbntns = -1; 
  }
 
  crrnt_trial  = 0;

  if (nmbr_empty < 1024) {
    infoMsg(IM::Tiny, "%d nodes, %d missing edges, %0e combinations\n", 
      nmbr_nodes , nmbr_empty, nmbr_cmbntns );
  } 
  else {
    infoMsg(IM::Tiny, "%d nodes, %d missing edges, 2^%d combinations\n", 
      nmbr_nodes , nmbr_empty, nmbr_empty );
  } 

  while (!done) {

    crrnt_trial++;
    if ((crrnt_trial % 0x10000) == 0) {

      if (nmbr_empty < 1024) { 
        infoMsg(IM::Tiny, "%e of %0e\n", (double)crrnt_trial, nmbr_cmbntns );
      }
      else {
        infoMsg(IM::Tiny, "%e of 2^%d\n", (double)crrnt_trial, nmbr_empty);
      }
    }
 
    //////////////////////////////////////////////////////////////////
    // Test current configuration  
    //////////////////////////////////////////////////////////////////
    maximumCardinalitySearch( triangulate_nodes, cliques, order, false);
    chordal = testZeroFillIn( order );

    if (chordal) {

      listVectorCliquetoVectorSetClique( cliques, vector_cliques );
      weight = graphWeight(vector_cliques,jtWeight,nodesRootMustContain);

      if (weight < best_weight) {
        infoMsg(IM::Tiny, "----- New Best: %f -----\n", weight); 

        infoMsg(IM::Tiny, "    ");
        for( i=0; i<nmbr_nodes; i++) {
          infoMsg(IM::Tiny, "[%2d]", i );
        } 
        infoMsg(IM::Tiny, "\n"); 

        for( i=0; i<nmbr_nodes; i++) {
          infoMsg(IM::Tiny, "[%2d]", i );
          for( j=0; j<nmbr_nodes; j++) {
            infoMsg(IM::Tiny, " %2d ", adjacency[i][j]);
          } 
          infoMsg(IM::Tiny, "\n"); 
        } 

        for( crrnt_clique = cliques.begin(), 
             end_clique   = cliques.end();
             crrnt_clique != end_clique; 
             ++crrnt_clique ) {

          infoMsg(IM::Tiny, "--- Clique ---\n"); 
          for( crrnt_node = (*crrnt_clique).begin(), 
	       end_node   = (*crrnt_clique).end(); 
               crrnt_node != end_node;
	       crrnt_node++ ) { 
            infoMsg(IM::Tiny, "%s(%d)\n", 
              (*crrnt_node)->randomVariable->name().c_str(), 
              (*crrnt_node)->randomVariable->timeIndex); 
          }
        } 
        infoMsg(IM::Tiny, "--------------------------\n");

        best_list_cliques = cliques;
        best_weight = weight;
      } 
    } 

    //////////////////////////////////////////////////////////////////
    // Set up next edge configuration 
    //////////////////////////////////////////////////////////////////
    row = 0;    
    col = 1;

    //////////////////////////////////////////////////////////////
    // Find next empty edge 
    //////////////////////////////////////////////////////////////
    while ( ((adjacency[row][col] == EDGE) || 
             (adjacency[row][col] == FILL_IN)) &&
            (row < nmbr_nodes) ) {
 
      col++;
      if (col >= nmbr_nodes) {
        row++;
        col = row+1;
      }
    }

    //////////////////////////////////////////////////////////////
    // Set edge, erase prior edges 
    //////////////////////////////////////////////////////////////
    if (col >= nmbr_nodes) { 
      done = true;
    } 
    else {

      adjacency[row][col] = FILL_IN;
      triangulate_nodes[row].neighbors.push_back(&(triangulate_nodes[col]));
      triangulate_nodes[col].neighbors.push_back(&(triangulate_nodes[row]));

      i = 0;
      j = 1;

      while ((i != row) || (j != col)) {
        if (adjacency[i][j] == FILL_IN) {
          adjacency[i][j] = NONE;

          crrnt_node = find( triangulate_nodes[j].neighbors.begin(),
            triangulate_nodes[j].neighbors.end(), &triangulate_nodes[i] );
          triangulate_nodes[j].neighbors.erase(crrnt_node);

          crrnt_node = find( triangulate_nodes[i].neighbors.begin(),
            triangulate_nodes[i].neighbors.end(), &triangulate_nodes[j] );
          triangulate_nodes[i].neighbors.erase(crrnt_node);
        }

        j++;
        if (j >= nmbr_nodes) {
          i++;
          j = i+1;
        }
      }
    } 

    //////////////////////////////////////////////////////////////
    // Check if timer has expired 
    //////////////////////////////////////////////////////////////
    if (timer && timer->Expired()) { 
      infoMsg(IM::Tiny, 
        "Time expired before completion of exhaustive search\n");
      done = true;
    }

  }

  //////////////////////////////////////////////////////////////////////////
  // Convert to MaxClique and triangulate the RandomVariable structures
  //////////////////////////////////////////////////////////////////////////
  listVectorCliquetoVectorSetClique( best_list_cliques, best_cliques );
  fillAccordingToCliques( best_cliques );

  infoMsg(IM::Tiny, "---->Tested:  %d\n", crrnt_trial); 
}


/*-
 *-----------------------------------------------------------------------
 * BoundaryTriangulate::triangulateElimination()
 *   Triangulate a set of nodes using an elimination order
 *    
 * Preconditions:
 *   Graph must be a valid undirected model. This means if the graph
 *   was originally a directed model, it must have been properly
 *   moralized and their 'neighbors' structure is valid. It is also
 *   assumed that the parents of each r.v. are valid but that they
 *   only point to variables within the set 'nodes' (i.e., parents
 *   must not point out of this set).
 *
 * Postconditions:
 *   Resulting graph is now triangulated.
 *
 * Side Effects:
 *   Changes rv's neighbors variables.
 *
 * Results:
 *     none
 *
 *-----------------------------------------------------------------------
 */
void
BoundaryTriangulate::
triangulateElimination(// input: nodes to be triangulated
                       const set<RandomVariable*> nodes,
                       // elimination order 
                       vector<RandomVariable*> orderedNodes,  
                       // output: resulting max cliques
                       vector<MaxClique>& cliques
                       )
{
  // Keep ordered eliminated nodes as a set for easy
  // intersection, with other node sets.
  set<RandomVariable*> orderedNodesSet;

  // Triangulate and make cliques
  for (unsigned i=0;i<orderedNodes.size();i++) {
    RandomVariable* rv = orderedNodes[i];

    
    // connect all neighbors of r.v. excluding nodes in 'orderedNodesSet'.
    rv->connectNeighbors(orderedNodesSet); 
    

    // check here if this node + its neighbors is a subset
    // of previous maxcliques. If it is not a subset of any previous
    // maxclique, then this node and its neighbors is a 
    // new maxclique.
    set<RandomVariable*> candidateMaxClique;
    set_difference(rv->neighbors.begin(),rv->neighbors.end(),
                   orderedNodesSet.begin(),orderedNodesSet.end(),
                   inserter(candidateMaxClique,candidateMaxClique.end()));
    candidateMaxClique.insert(rv);
    bool is_max_clique = true;
    for (unsigned i=0;i < cliques.size(); i++) {
      if (includes(cliques[i].nodes.begin(),cliques[i].nodes.end(),
                   candidateMaxClique.begin(),candidateMaxClique.end())) {
        // then found a 'proven' maxclique that includes our current
        // candidate, so the candidate cannot be a maxclique
        is_max_clique = false;
        break;
      }
    }
    if (is_max_clique) {
      cliques.push_back(MaxClique(candidateMaxClique));
    }

    // insert node into ordered set
    orderedNodesSet.insert(rv);
  }
}


/*-
 *-----------------------------------------------------------------------
 * BoundaryTriangulate::unrollAndTriangulate(th,int)
 *  Just unroll the graph a given number of times, triangulate it
 *  using the supplied heuristics, and report the results.
 *  This routine corresponds to unconstrained triangulation.
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
BoundaryTriangulate::
unrollAndTriangulate(// triangulate heuristics
		     const string& tri_heur_str,  
		     // number of times it should be unrolled
		     const unsigned numTimes)
{
  TriangulateHeuristics tri_heur;
  parseTriHeuristicString(tri_heur_str,tri_heur);
  const set <RandomVariable*> emptySet;

  if (numTimes >= 0) {
    vector <RandomVariable*> rvs;
    set <RandomVariable*> rvsSet;
    fp.unroll(numTimes,rvs);
    for (unsigned i=0;i<rvs.size();i++) {
      rvs[i]->createNeighborsFromParentsChildren();
    }
    for (unsigned i=0;i<rvs.size();i++) {
      rvs[i]->moralize();
      rvsSet.insert(rvs[i]);
    }
    vector<MaxClique> cliques;
    vector<nghbrPairType> orgnl_nghbrs;
    saveCurrentNeighbors(rvsSet,orgnl_nghbrs);
    string best_meth_str;
    double best_weight = HUGE_VAL;
    triangulate(rvsSet,false,emptySet,tri_heur,orgnl_nghbrs,cliques,best_meth_str,best_weight);
    unsigned maxSize = 0;
    float maxSizeCliqueWeight=0;
    float maxWeight = -1.0;
    float totalWeight = -1.0; // starting flag
    unsigned maxWeightCliqueSize=0;

    // TODO: just print out the resulting information for now. Ultimately
    // return the cliques and do inference with them.

    printf("Cliques from graph unrolled %d times\n",numTimes);
    for (unsigned i=0;i<cliques.size();i++) {

      float curWeight = MaxClique::computeWeight(cliques[i].nodes);
      if (curWeight > maxWeight) {
	maxWeight = curWeight;
	maxWeightCliqueSize = cliques[i].nodes.size();
      }
      if (totalWeight == -1.0)
	totalWeight = curWeight;
      else
	totalWeight = totalWeight + log10(1+pow(10,curWeight-totalWeight));
      if (cliques[i].nodes.size() > maxSize) {
	maxSize = cliques[i].nodes.size();
	maxSizeCliqueWeight = curWeight;
      }
      printf("%d : %d  %f\n",i,
	     cliques[i].nodes.size(),curWeight);
      for (set<RandomVariable*>::iterator j=cliques[i].nodes.begin();
	   j != cliques[i].nodes.end(); j++) {
	RandomVariable* rv = (*j);
	printf("   %s(%d)\n",rv->name().c_str(),rv->frame());
      }
    }
    printf("When unrolling %d times, max size clique = %d (with a weight of %f) and max weight of a clique = %f (with a size of %d). Total state space = %f\n",
	   numTimes,
	   maxSize,
	   maxSizeCliqueWeight,
	   maxWeight,
	   maxWeightCliqueSize,
	   totalWeight);
    printf("\n");
  }
}


/*-
 *-----------------------------------------------------------------------
 * BoundaryTriangulate::ensurePartitionsAreChordal(gm_template)
 *   ensure that the partitions in the given template are chordal, and die 
 *   with an error if not.
 *
 * Preconditions:
 *   Template should be fully instantiated
 *
 * Postconditions:
 *   Partitions in template are guaranteed to be chordal, if program has 
 *   not died.
 *
 * Side Effects:
 *   Re-writes the maxcliques in the current gm_template (i.e.,
 *   they'll still be maxcliques, but now they come from and are
 *   ordered by MCS).
 *
 * Results:
 *   none
 *-----------------------------------------------------------------------
 */
void
BoundaryTriangulate::
ensurePartitionsAreChordal(GMTemplate& gm_template)
{
  vector<triangulateNode>         triangulate_nodes; 
  list<vector<triangulateNode*> > list_cliques;
  vector<triangulateNode*>        order; 
  bool p_chordal;
  bool c_chordal;
  bool e_chordal;
 
  fillTriangulateNodeStructures( gm_template.P.nodes, triangulate_nodes );
  maximumCardinalitySearch( triangulate_nodes, list_cliques, order, false);
  p_chordal = testZeroFillIn( order );

  triangulate_nodes.clear();
  fillTriangulateNodeStructures( gm_template.C.nodes, triangulate_nodes );
  maximumCardinalitySearch( triangulate_nodes, list_cliques, order, false);
  c_chordal = testZeroFillIn( order );

  triangulate_nodes.clear();
  fillTriangulateNodeStructures( gm_template.E.nodes, triangulate_nodes );
  maximumCardinalitySearch( triangulate_nodes, list_cliques, order, false);
  e_chordal = testZeroFillIn( order );

  if (!p_chordal || !c_chordal || !e_chordal) {
    error("ERROR: Program exiting since the following partitions are not chordal:%s%s%s",
	  (p_chordal?"":" P"),
	  (c_chordal?"":" C"),
	  (e_chordal?"":" E"));
  }
}


/*******************************************************************************
 *******************************************************************************
 *******************************************************************************
 **
 **        Support for Anytime Triangulation.
 **
 *******************************************************************************
 *******************************************************************************
 *******************************************************************************
 */


/*-
 *-----------------------------------------------------------------------
 * BoundaryTriangulate::anyTimeTriangulate()
 *  Triangulate the P,C, and E partitions using an anytime procedure
 *
 * Preconditions:
 *   Each partition must corresond to a valid and separate GM. Each
 *   variable in each GM must have a valid parent and neighbor members
 *   and the parents/neighbors must only point to other members of a
 *   given partition.
 *   Current class timer must be a valid pointer to a timer.
 * 
 * Postconditions:
 *   Each of the partitions have been triangulated to the graph giving
 *   the best weight that was found. 
 *
 * Side Effects:
 *   Neighbor members of each random variable can be changed.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
BoundaryTriangulate::
anyTimeTriangulate(GMTemplate& gm_template,
		   const bool jtWeight,
		   const bool doP,
		   const bool doC,
		   const bool doE)
{
  assert (timer != NULL);

  nghbrPairType nghbr_pair;
  vector<nghbrPairType> orgnl_P_nghbrs;
  vector<nghbrPairType> orgnl_C_nghbrs;
  vector<nghbrPairType> orgnl_E_nghbrs;
  string best_P_method_str;
  string best_C_method_str;
  string best_E_method_str;
  double best_P_weight = HUGE_VAL;
  double best_C_weight = HUGE_VAL;
  double best_E_weight = HUGE_VAL;
  const set <RandomVariable*> emptySet;

  ////////////////////////////////////////////////////////////////////////
  // Save the untriangulated graphs so that they can be quickly restored
  // for multiple iterations
  ////////////////////////////////////////////////////////////////////////
  if (doP) saveCurrentNeighbors(gm_template.P,orgnl_P_nghbrs);
  if (doC) saveCurrentNeighbors(gm_template.C,orgnl_C_nghbrs);
  if (doE) saveCurrentNeighbors(gm_template.E,orgnl_E_nghbrs);

  ////////////////////////////////////////////////////////////////////////
  // Triangulate P and E partitions using a basic heuristic so that 
  // something valid is in place in case these never complete before timer 
  // expires.
  ////////////////////////////////////////////////////////////////////////

  if (doP) {
    infoMsg(IM::Tiny, "---\nPreliminary Triangulation of P:\n");
    triangulate( gm_template.P.nodes, jtWeight, 
		 gm_template.PCInterface_in_P,
		 "FWH", orgnl_P_nghbrs, gm_template.P.cliques,
		 gm_template.P.triMethod, best_P_weight );
  }

  if (doE) {
    infoMsg(IM::Tiny, "---\nPreliminary Triangulation of E:\n");
    triangulate( gm_template.E.nodes, jtWeight, 
		 emptySet,
		 "FWH", orgnl_E_nghbrs, gm_template.E.cliques, 
		 gm_template.E.triMethod, best_E_weight ); 
  }

  ////////////////////////////////////////////////////////////////////////
  // Triangulate using a variety of heuristic searches 
  ////////////////////////////////////////////////////////////////////////

  // Like above, always do at least C, so we return something valid.
  if (doC) {
    infoMsg(IM::Tiny, "---\nTriangulating C using Heuristics:\n");
    best_C_weight = tryEliminationHeuristics( gm_template.C.nodes, jtWeight, 
      gm_template.CEInterface_in_C, orgnl_C_nghbrs, gm_template.C.cliques, 
      gm_template.C.triMethod );
    best_C_weight = tryNonEliminationHeuristics( gm_template.C.nodes, jtWeight, 
      gm_template.CEInterface_in_C, orgnl_C_nghbrs, gm_template.C.cliques, 
      gm_template.C.triMethod );
  }
 
  if (doP && !timer->Expired()) { 
    infoMsg(IM::Tiny, "---\nTriangulating P using Heuristics:\n");
    best_P_weight = tryEliminationHeuristics( gm_template.P.nodes, 
      jtWeight, gm_template.PCInterface_in_P, orgnl_P_nghbrs, 
      gm_template.P.cliques, gm_template.P.triMethod );
    best_P_weight = tryNonEliminationHeuristics( gm_template.P.nodes, 
      jtWeight, gm_template.PCInterface_in_P, orgnl_P_nghbrs, 
      gm_template.P.cliques, gm_template.P.triMethod );
  }
  
  if (doE && !timer->Expired()) { 
    infoMsg(IM::Tiny, "---\nTriangulating E using Heuristics:\n");
    best_E_weight = tryEliminationHeuristics( gm_template.E.nodes, jtWeight, 
      emptySet, orgnl_E_nghbrs, gm_template.E.cliques, gm_template.E.triMethod);
    best_E_weight = tryNonEliminationHeuristics( gm_template.E.nodes, jtWeight, 
      emptySet, orgnl_E_nghbrs, gm_template.E.cliques, gm_template.E.triMethod);
  }
  
  infoMsg(IM::Tiny, "Time Remaining: %d\n", (int)timer->SecondsLeft() ); 

  ////////////////////////////////////////////////////////////////////////
  // Triangulate using simulated annealing 
  ////////////////////////////////////////////////////////////////////////
  if (doC && timer->SecondsLeft() > 10) {
    // Do C first since it is more important.
    infoMsg(IM::Tiny, "---\nTriangulating C using Simulated Annealing:\n");

    triangulate( gm_template.C.nodes, jtWeight, gm_template.CEInterface_in_C,
		 "anneal", orgnl_C_nghbrs, gm_template.C.cliques, 
		 gm_template.C.triMethod, best_C_weight ); 

    infoMsg(IM::Tiny, "Time Remaining: %d\n", (int)timer->SecondsLeft() ); 
  }

  if (doP && timer->SecondsLeft() > 10) {
    infoMsg(IM::Tiny, "---\nTriangulating P using Simulated Annealing:\n");

    triangulate( gm_template.P.nodes, jtWeight, gm_template.PCInterface_in_P,
		 "anneal", orgnl_P_nghbrs, gm_template.P.cliques, 
		 gm_template.P.triMethod, best_P_weight ); 

    infoMsg(IM::Tiny, "Time Remaining: %d\n", 
      (int)timer->SecondsLeft()); 
  }

  if (doE && timer->SecondsLeft() > 10) {
    infoMsg(IM::Tiny, "---\nTriangulating E using Simulated Annealing:\n");

    triangulate( gm_template.E.nodes, jtWeight, emptySet,
		 "anneal", orgnl_E_nghbrs, gm_template.E.cliques, 
		 gm_template.E.triMethod, best_E_weight ); 

    infoMsg(IM::Tiny, "Time Remaining: %d\n", (int)timer->SecondsLeft() ); 
  }

  ////////////////////////////////////////////////////////////////////////
  // Triangulate using exhaustive search
  ////////////////////////////////////////////////////////////////////////

  if (doC && timer->SecondsLeft() > 10) {
    infoMsg(IM::Tiny, "Triangulating C using Exhaustive Search:\n");
    
    triangulate( gm_template.C.nodes, jtWeight, gm_template.CEInterface_in_C,
		 "exhaustive", orgnl_C_nghbrs, gm_template.C.cliques, 
		 gm_template.C.triMethod, best_C_weight ); 

    infoMsg(IM::Tiny, "Time Remaining: %d\n", (int)timer->SecondsLeft() ); 
  }

  if (doP && timer->SecondsLeft() > 10) {
    infoMsg(IM::Tiny, "Triangulating P using Exhaustive Search:\n");

    triangulate( gm_template.P.nodes, jtWeight, gm_template.PCInterface_in_P,
		 "exhaustive", orgnl_P_nghbrs, gm_template.P.cliques, 
		 gm_template.P.triMethod, best_P_weight ); 

    infoMsg(IM::Tiny, "Time Remaining: %d\n", (int)timer->SecondsLeft() ); 
  }

  if (doE && timer->SecondsLeft() > 10) {
    infoMsg(IM::Tiny, "Triangulating E using Exhaustive Search:\n");

    triangulate( gm_template.E.nodes, jtWeight, emptySet,
		 "exhaustive", orgnl_E_nghbrs, gm_template.E.cliques, 
		 gm_template.E.triMethod, best_E_weight ); 

    infoMsg(IM::Tiny, "Time Remaining: %d\n", (int)timer->SecondsLeft() ); 
  }

  ////////////////////////////////////////////////////////////////////////
  // Return with the best triangulations found, which is
  // be stored within the template at this point.
  ////////////////////////////////////////////////////////////////////////
  if (doP) restoreNeighbors(orgnl_P_nghbrs);
  if (doC) restoreNeighbors(orgnl_C_nghbrs);
  if (doE) restoreNeighbors(orgnl_E_nghbrs);
  gm_template.triangulatePartitionsByCliqueCompletion();

}


/*-
 *-----------------------------------------------------------------------
 * BoundaryTriangulate::tryEliminationHeuristics
 *  Triangulate a partition using elimination with multiple iterations of 
 *  a variety of heuristics.  Use overal graph weight to judge ultimate 
 *  triangulation quality.
 *
 *  There are two ways to judget the quality of the triangulation. 
 *  If jtWeight = false, use just sum of weights in each resulting clique.
 *  If jtWeight = true, approximate JT quality by forming JT and computing 
 *  JT weight.
 *
 * Preconditions:
 *   Each variable in the set of nodes must have valid parent and 
 *   neighbor members and the parents/neighbors must only point to other 
 *   members in the set. 
 * 
 * Postconditions:
 *   The triangulation with the best weight is stored in 'cliques' 
 *
 * Side Effects:
 *   The partition is triangulated to some triangulation, not 
 *   necessarily corresponding to the best order.  Neighbor members of 
 *   each random variable can be changed.
 *
 * Results:
 *   The best weight found 
 *-----------------------------------------------------------------------
 */
double 
BoundaryTriangulate::
tryEliminationHeuristics(
  const set<RandomVariable*>& nodes,
  const bool                  jtWeight,
  const set<RandomVariable*>& nrmc,        // nrmc = nodes root must contain
  vector<nghbrPairType>&      orgnl_nghbrs,
  vector<MaxClique>&          cliques,
  string&                     tri_method
  )
{
  double best_weight = HUGE_VAL;

  ///////////////////////////////////////////////////////////////////////////// 
  // Try Weight, Fill, Size
  triangulate( nodes, jtWeight, nrmc, "20-WFS",  
	       orgnl_nghbrs, cliques, tri_method, best_weight );
  triangulate( nodes, jtWeight, nrmc, "20-HWFS", 
	       orgnl_nghbrs, cliques, tri_method, best_weight );

  ///////////////////////////////////////////////////////////////////////////// 
  // Try Fill, Weight, Size 
  triangulate( nodes, jtWeight, nrmc, "20-FWS",  
	       orgnl_nghbrs, cliques, tri_method, best_weight );
  triangulate( nodes, jtWeight, nrmc, "20-HFWS", 
	       orgnl_nghbrs, cliques, tri_method, best_weight ); 

  ///////////////////////////////////////////////////////////////////////////// 
  // Try Weight
  triangulate( nodes, jtWeight, nrmc, "30-W",  
	       orgnl_nghbrs, cliques, tri_method, best_weight );
  triangulate( nodes, jtWeight, nrmc, "20-HW", 
	       orgnl_nghbrs, cliques, tri_method, best_weight );
  triangulate( nodes, jtWeight, nrmc, "20-TW", 
	       orgnl_nghbrs, cliques, tri_method, best_weight );

  ///////////////////////////////////////////////////////////////////////////// 
  // Try Fill 
  triangulate( nodes, jtWeight, nrmc, "30-F",  
	       orgnl_nghbrs, cliques, tri_method, best_weight );
  triangulate( nodes, jtWeight, nrmc, "20-HF", 
	       orgnl_nghbrs, cliques, tri_method, best_weight );
  triangulate( nodes, jtWeight, nrmc, "20-TF", 
	       orgnl_nghbrs, cliques, tri_method, best_weight );

  ///////////////////////////////////////////////////////////////////////////// 
  // Try Size 
  triangulate( nodes, jtWeight, nrmc, "30-S",  
	       orgnl_nghbrs, cliques, tri_method, best_weight ); 
  triangulate( nodes, jtWeight, nrmc, "20-HS", 
	       orgnl_nghbrs, cliques, tri_method, best_weight ); 
  triangulate( nodes, jtWeight, nrmc, "20-ST", 
	       orgnl_nghbrs, cliques, tri_method, best_weight ); 

  ///////////////////////////////////////////////////////////////////////////// 
  // Try Time
  triangulate( nodes, jtWeight, nrmc, "100-T",  
	       orgnl_nghbrs, cliques, tri_method, best_weight ); 
  triangulate( nodes, jtWeight, nrmc, "30-TS", 
	       orgnl_nghbrs, cliques, tri_method, best_weight ); 
  triangulate( nodes, jtWeight, nrmc, "30-TW", 
	       orgnl_nghbrs, cliques, tri_method, best_weight ); 
  triangulate( nodes, jtWeight, nrmc, "30-TF", 
	       orgnl_nghbrs, cliques, tri_method, best_weight ); 


  ///////////////////////////////////////////////////////////////////////////// 
  // Try Position
  triangulate( nodes, jtWeight, nrmc, "P",  
	       orgnl_nghbrs, cliques, tri_method, best_weight ); 
  triangulate( nodes, jtWeight, nrmc, "30-WP", 
	       orgnl_nghbrs, cliques, tri_method, best_weight ); 
  triangulate( nodes, jtWeight, nrmc, "30-FP", 
	       orgnl_nghbrs, cliques, tri_method, best_weight ); 
  triangulate( nodes, jtWeight, nrmc, "30-SP", 
	       orgnl_nghbrs, cliques, tri_method, best_weight ); 

  ///////////////////////////////////////////////////////////////////////////// 
  // Try Maximum Cardinality Search
  triangulate( nodes, jtWeight, nrmc, "50-MCS", 
	       orgnl_nghbrs, cliques, tri_method, best_weight ); 

  ///////////////////////////////////////////////////////////////////////////// 
  // Try purely random
  triangulate( nodes, jtWeight, nrmc, "500-R", 
	       orgnl_nghbrs, cliques, tri_method, best_weight );

  ///////////////////////////////////////////////////////////////////////////// 
  // Try hints only 
  triangulate( nodes, jtWeight, nrmc, "30-H", 
	       orgnl_nghbrs, cliques, tri_method, best_weight );

  return(best_weight);
}


/*-
 *-----------------------------------------------------------------------
 * BoundaryTriangulate::tryNonEliminationHeuristics
 *  Triangulate a partition using elimination with multiple iterations of 
 *  a variety of heuristics.  All of the heuristics in this proceedure 
 *  have the potential to give triangulations not obtainable through any
 *  elimination order.  Use overal graph weight to judge ultimate 
 *  triangulation quality.
 *
 *  There are two ways to judget the quality of the triangulation. 
 *  If jtWeight = false, use just sum of weights in each resulting clique.
 *  If jtWeight = true, approximate JT quality by forming JT and computing 
 *    JT weight.
 *
 * Preconditions:
 *   Each variable in the set of nodes must have valid parent and 
 *   neighbor members and the parents/neighbors must only point to other 
 *   members in the set. 
 * 
 * Postconditions:
 *   The triangulation with the best weight is stored in 'cliques' 
 *
 * Side Effects:
 *   The partition is triangulated to some triangulation, not 
 *   necessarily corresponding to the best order.  Neighbor members of 
 *   each random variable can be changed.
 *
 * Results:
 *   The best weight found 
 *-----------------------------------------------------------------------
 */
double 
BoundaryTriangulate::
tryNonEliminationHeuristics(
  const set<RandomVariable*>& nodes,
  const bool                  jtWeight,
  const set<RandomVariable*>& nrmc,         // nrmc = nodes root must contain
  vector<nghbrPairType>&      orgnl_nghbrs,
  vector<MaxClique>&          cliques,
  string&                     tri_method
  )
{
  double best_weight = HUGE_VAL;

  ///////////////////////////////////////////////////////////////////////////// 
  // Try adding all ancestral edges, followed by elimination heuristics 
  triangulate( nodes, jtWeight, nrmc, "pre-edge-all", 
	       orgnl_nghbrs, cliques, tri_method, best_weight );

  ///////////////////////////////////////////////////////////////////////////// 
  // Try adding locally optimal ancestral edges, followed by elimination 
  // heuristics 
  triangulate( nodes, jtWeight, nrmc, "pre-edge-lo", 
	       orgnl_nghbrs, cliques, tri_method, best_weight );

  ///////////////////////////////////////////////////////////////////////////// 
  // Try adding random optimal ancestral edges, followed by elimination 
  // heuristics 
  triangulate( nodes, jtWeight, nrmc, "10-pre-edge-random", 
	       orgnl_nghbrs, cliques, tri_method, best_weight );

  ///////////////////////////////////////////////////////////////////////////// 
  // Try frontier algorithm
  triangulate( nodes, jtWeight, nrmc, "500-frontier", 
	       orgnl_nghbrs, cliques, tri_method, best_weight );

  ///////////////////////////////////////////////////////////////////////////// 
  // Try simply completing the partition (this can work well if many 
  // deterministic nodes exist in the partition)
  triangulate( nodes, jtWeight, nrmc, "completed", orgnl_nghbrs, cliques, 
    tri_method, best_weight );

  return(best_weight);
}


/*******************************************************************************
 *******************************************************************************
 *******************************************************************************
 **
 **         Boundary Algorithm Routines
 **
 *******************************************************************************
 *******************************************************************************
 *******************************************************************************
 */



/*-
 *-----------------------------------------------------------------------
 * BoundaryTriangulate::interfaceScore()
 *      Compute the 'score' of a set of variables that are to be a
 *      candidate interface. In the best of cases, it is a simple easy
 *      score based on the vars size, fill-in, or weight.  In the
 *      "worst" of cases, the score is based on an entire
 *      triangulation of the variables that would result from this
 *      interface.
 *      
 *
 * Preconditions:
 *      - set of variables must be instantiated and their
 *        neighbors member been filled in.
 *      - score heuristic vector instantiated.
 *
 * Postconditions:
 *      - score in vector is returned.
 *
 * Side Effects:
 *      none
 *
 * Results:
 *     - score in vector returned.
 *
 *-----------------------------------------------------------------------
 */
void
BoundaryTriangulate::interfaceScore(
 // the interface heuristic used to score the interface
 const vector<BoundaryHeuristic>& bnd_heur_v,
 // the interface itself
 const set<RandomVariable*>& C_l,
 // --------------------------------------------------------------
 // The next 8 input arguments are used only with the optimal
 // interface algorithm when the IH_MIN_MAX_C_CLIQUE or
 // IH_MIN_MAX_CLIQUE heuristics are used:
 // triangulation heuristic
 // Variables to the left (or right) of the interface
 const set<RandomVariable*>& left_C_l,
 // the triangulation heuristic 
 const TriangulateHeuristics& tri_heur,
 // The network unrolled 1 time
 const set<RandomVariable*>& P_u1,
 const set<RandomVariable*>& C1_u1,
 const set<RandomVariable*>& Cextra_u1,
 const set<RandomVariable*>& C2_u1,
 const set<RandomVariable*>& E_u1,
 // Mappings from C2 in the twice unrolled network to C1 and C2
 // in the once unrolled network.
 // (these next 2 should be const, but there is no "op[] const")
 map < RandomVariable*, RandomVariable* >& C2_u2_to_C1_u1,
 map < RandomVariable*, RandomVariable* >& C2_u2_to_C2_u1,
 // --------------------------------------------------------------
 // output score
 vector<float>& score)
{

  if (message(Moderate)) {
    set<RandomVariable*>::iterator i;    
    printf("  --- Cur Interface:");
    for (i=C_l.begin();i!=C_l.end();i++) {
      printf(" %s(%d)",
	     (*i)->name().c_str(),
	     (*i)->frame());
      
    }
    if (message(VerbosityLevels(Moderate + (Increment>>1)))) {
      // also print out left_C_l
      printf("  Left of C_l:");
      for (i=left_C_l.begin();i!=left_C_l.end();i++) {
	printf(" %s(%d)",
	       (*i)->name().c_str(),
	       (*i)->frame());
      }
    }
    printf("\n");
  }
  score.clear();
  for (unsigned fhi=0;fhi<bnd_heur_v.size();fhi++) {
    const BoundaryHeuristic fh = bnd_heur_v[fhi];
    if (fh == IH_MIN_WEIGHT || fh == IH_MIN_WEIGHT_NO_D) {
      float tmp_weight = MaxClique::computeWeight(C_l,NULL,
				       (fh == IH_MIN_WEIGHT));
      score.push_back(tmp_weight);
      infoMsg(Low,"  Interface Score: set has weight = %f\n",tmp_weight);
    } else if (fh == IH_MIN_FILLIN) {
      int fill_in = computeFillIn(C_l);
      score.push_back((float)fill_in);
      infoMsg(Low,"  Interface Score: set has fill_in = %d\n",fill_in);
    } else if (fh == IH_MIN_SIZE) {
      score.push_back((float)C_l.size());
      infoMsg(Low,"  Interface Score: set has size = %d\n",
	      C_l.size());
    } else if (fh == IH_MIN_MIN_POSITION_IN_FILE) {
      set<RandomVariable*>::iterator i;    
      unsigned val = ~0x0;
      for (i=C_l.begin();i!=C_l.end();i++) {
	if (unsigned((*i)->rv_info.variablePositionInStrFile) < val)
	  val = (*i)->rv_info.variablePositionInStrFile;
      }
      score.push_back((float)val);
      infoMsg(Low,"  Interface Score: set has min pos = %d\n",val);
    } else if (fh == IH_MIN_MAX_POSITION_IN_FILE) {
      set<RandomVariable*>::iterator i;    
      unsigned val = 0;
      for (i=C_l.begin();i!=C_l.end();i++) {
	if (unsigned((*i)->rv_info.variablePositionInStrFile) > val)
	  val = (*i)->rv_info.variablePositionInStrFile;
      }
      score.push_back((float)val);
      infoMsg(Low,"  Interface Score: set has max pos = %d\n",val);
    } else if (fh == IH_MIN_MIN_TIMEFRAME) {
      set<RandomVariable*>::iterator i;    
      unsigned val = ~0x0;
      for (i=C_l.begin();i!=C_l.end();i++) {
	if ((*i)->frame() < val)
	  val = (*i)->frame();
      }
      score.push_back((float)val);
      infoMsg(Low,"  Interface Score: set has min timeframe = %d\n",val);
    } else if (fh == IH_MIN_MAX_TIMEFRAME) {
      set<RandomVariable*>::iterator i;    
      unsigned val = 0;
      for (i=C_l.begin();i!=C_l.end();i++) {
	if ((*i)->frame() > val)
	  val = (*i)->frame();
      }
      score.push_back((float)val);
      infoMsg(Low,"  Interface Score: set has max timeframe = %d\n",val);
    } else if (fh == IH_MIN_MIN_HINT) {
      set<RandomVariable*>::iterator i;    
      float val = MAXFLOAT;
      for (i=C_l.begin();i!=C_l.end();i++) {
	if ((*i)->rv_info.eliminationOrderHint < val)
	  val = (*i)->rv_info.eliminationOrderHint;
      }
      score.push_back(val);
      infoMsg(Low,"  Interface Score: set has min hint = %f\n",val);
    } else if (fh == IH_MIN_MAX_HINT) {
      set<RandomVariable*>::iterator i;    
      float val = -MAXFLOAT;
      for (i=C_l.begin();i!=C_l.end();i++) {
	if ((*i)->rv_info.eliminationOrderHint > val)
	  val = (*i)->rv_info.eliminationOrderHint;
      }
      score.push_back(val);
      infoMsg(Low,"  Interface Score: set has max hint = %f\n",val);
    } else if (fh == IH_MIN_MAX_CLIQUE || fh == IH_MIN_MAX_C_CLIQUE ||
	       fh == IH_MIN_STATE_SPACE || fh == IH_MIN_C_STATE_SPACE) {
      // This is the expensive one, need to form a set of partitions,
      // given the current interface, triangulate that partition set,
      // and then compute the score of the worst clique, and fill the
      // score variable above with this worst scoring clique (i.e., we
      // find the interface that has the best worst-case performance.

      // TODO: some of these values could be cached to speed this
      // up a bit.

      GMTemplate gm_template(fp,M,S);
      findInterfacePartitions(P_u1,
			      C1_u1,
			      Cextra_u1,
			      C2_u1,
			      E_u1,
			      C2_u2_to_C1_u1,
			      C2_u2_to_C2_u1,
			      left_C_l,
			      C_l,
			      gm_template);

      if (fh == IH_MIN_MAX_CLIQUE || IH_MIN_STATE_SPACE) {
	// then we do the entire graph, all three partitions
	vector<nghbrPairType> orgnl_P_nghbrs;
	vector<nghbrPairType> orgnl_C_nghbrs;
	vector<nghbrPairType> orgnl_E_nghbrs;
	string best_P_method_str;
	string best_C_method_str;
	string best_E_method_str;
	double best_P_weight = HUGE_VAL;
	double best_C_weight = HUGE_VAL;
	double best_E_weight = HUGE_VAL;

	saveCurrentNeighbors(gm_template.P,orgnl_P_nghbrs);
	saveCurrentNeighbors(gm_template.C,orgnl_C_nghbrs);
	saveCurrentNeighbors(gm_template.E,orgnl_E_nghbrs);

	// TODO: use version that includes jtWeight
	triangulate(gm_template.P.nodes,tri_heur,orgnl_P_nghbrs,gm_template.P.cliques,
		    gm_template.P.triMethod,best_P_weight);
	triangulate(gm_template.C.nodes,tri_heur,orgnl_C_nghbrs,gm_template.C.cliques,
		    gm_template.C.triMethod,best_C_weight);
	triangulate(gm_template.E.nodes,tri_heur,orgnl_E_nghbrs,gm_template.E.cliques,
		    gm_template.E.triMethod,best_E_weight);

	if (fh == IH_MIN_STATE_SPACE) {
	  // just sum up the weights
	  double weight = best_P_weight;
	  weight += log10(1+pow(10,best_C_weight-weight));
	  weight += log10(1+pow(10,best_E_weight-weight));
	  score.push_back(weight);
	  infoMsg(Low,"  Interface Score: resulting graph with this interface has total weight = %f\n",
		  weight);
	} else {
	  // need to find the clique with max weight.
	  float maxWeight = -1.0;
	  for (unsigned i=0;i<gm_template.P.cliques.size();i++) {
	    float curWeight = MaxClique::computeWeight(gm_template.P.cliques[i].nodes);
	    if (curWeight > maxWeight) maxWeight = curWeight;
	  }
	  for (unsigned i=0;i<gm_template.C.cliques.size();i++) {
	    float curWeight = MaxClique::computeWeight(gm_template.C.cliques[i].nodes);
	    if (curWeight > maxWeight) maxWeight = curWeight;
	  }
	  for (unsigned i=0;i<gm_template.E.cliques.size();i++) {
	    float curWeight = MaxClique::computeWeight(gm_template.E.cliques[i].nodes);
	    if (curWeight > maxWeight) maxWeight = curWeight;
	  }
	  score.push_back(maxWeight);
	  infoMsg(Low,"  Interface Score: with this interface, resulting graph has largest clique weight = %f\n",
		  maxWeight);
	}
      } else {
	// then we do almost same, but just on C
	vector<nghbrPairType> orgnl_C_nghbrs;
	string best_C_method_str;
	double best_C_weight;
	saveCurrentNeighbors(gm_template.C,orgnl_C_nghbrs);
	// TODO: use version that includes jtWeight
	triangulate(gm_template.C.nodes,tri_heur,orgnl_C_nghbrs,gm_template.C.cliques,best_C_method_str,best_C_weight);

	if (fh == IH_MIN_C_STATE_SPACE) {
	  // just sum up the weights
	  score.push_back(best_C_weight);
	  infoMsg(Low,"  Interface Score: resulting graph with this interface has total C-weight = %f\n",
		  best_C_weight);
	} else {
	  // need to find the clique with max weight.
	  float maxWeight = -1.0;
	  for (unsigned i=0;i<gm_template.C.cliques.size();i++) {
	    float curWeight = MaxClique::computeWeight(gm_template.C.cliques[i].nodes);
	    if (curWeight > maxWeight) maxWeight = curWeight;
	  }
	  score.push_back(maxWeight);
	  infoMsg(Low,"  Interface Score: with this interface, resulting C partition has largest clique weight = %f\n",
		  maxWeight);
	}
      }

    } else
      warning("Warning: invalid variable set score given. Ignored\n");
  }

  return;
}




/*-
 *-----------------------------------------------------------------------
 * BoundaryTriangulate::findBestInterface()
 *
 *  An exponential-time routine to find best left (or right)
 *  interfaces and corresponding partitions. Since the operations for
 *  finding the best left and right interfaces are symmetric, the case
 *  of if we are searching for the best left interface or right
 *  interface is determined entirely based on the arguments that are
 *  passed into these routines (i.e., just take the mirror image of
 *  the arguments for the right interface).  Note that for simplicity,
 *  the names have been defined in terms of the left interface.
 *
 *  The routine is given portions of an unrolled graph, P,C1,C2,C3,E,
 *  where C1 is one chunk, C2 is M chunks, and C3 is one chunk. The
 *  routine finds the best left (or right) interface entirely within
 *  C2 by searching for the boundary with the best interface.  The
 *  Boundary may span M chunks. The algorithm starts at the boundary
 *  corresponding to the "standard" or initial left interface between
 *  C1 and C2. Nodes in the left interface are shifted left across the
 *  boundary (thereby advancing the boundary). When a node is shifted
 *  left, then a new boundary forms. Note that the routine only uses
 *  the partitions C1,C2, and C3 as P and E are not needed (and can
 *  even be empty).
 *
 *
 * For left interface:
 *   Given a the unrolled graph, P,C1,C2,C3,E, find the best left
 *   interface within C2 starting at the boundary corresponding to the
 *   "standard" or initial left interface between C1 and C2.  Note
 *   that the routine only uses C1,C2, and C3 and P and E are not
 *   needed.
 *
 * For right interface:
 *   Given an unrolled graph, P,C1,C2,C3,E, find the best right
 *   interface within C2 starting at the "standard" or initial right
 *   interface between C3 and C2.  Note that the routine only uses
 *   C1,C2, and C3 and P and E are not needed. In order to get
 *   the behavior for right interface, the routine is called by
 *   swapping the arguments for C1 and C3 (see the caller).
 *
 *
 * Preconditions:
 *     Graph must be valid (i.e., unroller should pass graph w/o
 *                                problem).
 *     Graph must be already moralized.
 *     and *must* unrolled exactly M+1 times.
 *     This means graph must be in the form P,C1,C2,C3,E
 *     where P = prologue, 
 *           C1 = first chunk
 *           C2 = 2nd section, M chunks long
 *           C3 = last chunk
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
 *    Put best left (resp. right) interface from C2 into C_l
 *    and place any additional variables from C2 to the left (resp.
 *    right) of C_l into 'left_C_l';
 *
 *
 *-----------------------------------------------------------------------
 */
void
BoundaryTriangulate::findBestInterface(
 // first chunk of twice unrolled graph
 const set<RandomVariable*> &C1,
 // second chunk of twice unrolled graph
 const set<RandomVariable*> &C2,
 // first chunk of C2, empty when M=1
 const set<RandomVariable*> &C2_1,
 // third chunk of twice unrolled graph
 const set<RandomVariable*> &C3,
 // nodes to the "left" of the left interface within C2
 set<RandomVariable*> &left_C_l,
 // the starting left interface
 set<RandomVariable*> &C_l,
 // the resulting score of the best interface
 vector<float>& best_score,
 // what should be used to judge the quality of the interface
 const vector<BoundaryHeuristic>& bnd_heur_v,
 // true if we should use the exponential time optimal interface algorithm
 const bool recurse,
 // --------------------------------------------------------------
 // The next 7 input arguments are used only with the optimal
 // interface algorithm when the IH_MIN_MAX_C_CLIQUE or
 // IH_MIN_MAX_CLIQUE heuristics are used:
 // triangulation heuristic
 const TriangulateHeuristics& tri_heur,
 // The network unrolled 1 time (TODO: change to M+S-1)
 const set<RandomVariable*>& P_u1,
 const set<RandomVariable*>& C1_u1,
 const set<RandomVariable*>& Cextra_u1,
 const set<RandomVariable*>& C2_u1,
 const set<RandomVariable*>& E_u1,
 // Mappings from C2 in the twice unrolled network to C1 and C2
 // in the once unrolled network.
 // (these next 2 should be const, but there is no "op[] const" in STL)
 map < RandomVariable*, RandomVariable* >& C2_u2_to_C1_u1,
 map < RandomVariable*, RandomVariable* >& C2_u2_to_C2_u1
 // end of input arguments (finally)
 )
{

  // left interface
  C_l.clear();

  // First, construct the basic left interface (i.e., left
  // interface of C2).

  // go through through set C1, and pick out all neighbors
  // of variables in set C1 that live in C2, and these neighbors
  // become the initial left interface C_l
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
  // Since the interface is currently the first one, 
  // the set left of the left interface within C2 is now empty.
  left_C_l.clear();

  // Note that the partition "boundary" is the border/line that cuts
  // the edges connecting nodes C_l and nodes left_C_l.

  // how good is this interface.
  interfaceScore(bnd_heur_v,C_l,left_C_l,
		 tri_heur,
		 P_u1,C1_u1,Cextra_u1,C2_u1,E_u1,
		 C2_u2_to_C1_u1,C2_u2_to_C2_u1,
		 best_score);

  if (message(Tiny)) {
    printf("  Size of basic interface = %d\n",C_l.size());
    printf("  Score of basic interface =");
    for (unsigned i=0;i<best_score.size();i++)
      printf(" %f ",best_score[i]);
    printf("\n");
    // Note that this next print always prints the string "left_C_l"
    // but if we are running the right interface version, it should
    // print the string "right_C_l" instead. The algorithm however
    // does not know if it is running left or right interface.
    infoMsg(Med,"  Size of left_C_l = %d\n",left_C_l.size());
    {
      printf("  Interface nodes include:");
      set<RandomVariable*>::iterator i;    
      for (i=C_l.begin();i!=C_l.end();i++) {
	printf(" %s(%d)",
	       (*i)->name().c_str(),
	       (*i)->frame());
	     
      }
      printf("\n");
    }
  }

  // start exponential recursion to find the truly best interface.
  if (recurse) {
    // best ones found so far
    set<RandomVariable*> best_left_C_l = left_C_l;
    set<RandomVariable*> best_C_l = C_l;
    set< set<RandomVariable*> > setset;
    // call recursive routine.
    boundaryRecursionDepth = 0;
    findBestInterface(left_C_l,
		      C_l,
		      C2,
		      C2_1,
		      C3,
		      setset,
		      best_left_C_l,
		      best_C_l,
		      best_score,bnd_heur_v,
		      // 
		      tri_heur,
		      P_u1,C1_u1,Cextra_u1,C2_u1,E_u1,
		      C2_u2_to_C1_u1,C2_u2_to_C2_u1);
    if (message(Tiny)) {
      printf("  Size of best interface = %d\n",best_C_l.size());
      printf("  Score of best interface =");
      for (unsigned i=0;i<best_score.size();i++)
	printf(" %f ",best_score[i]);
      printf("\n");
      infoMsg(Med,"  Size of best_left_C_l = %d\n",best_left_C_l.size());
      {
	printf("  Best interface nodes include:");
	set<RandomVariable*>::iterator i;    
	for (i=best_C_l.begin();i!=best_C_l.end();i++) {
	  printf(" %s(%d)",
		 (*i)->name().c_str(),
		 (*i)->frame());
	}
	printf("\n");

      }
    }
    left_C_l = best_left_C_l;
    C_l = best_C_l;
  }
}
/*-
 *-----------------------------------------------------------------------
 * BoundaryTriangulate::findBestInterface() 
 *    Recursive helper function for the * first version of
 *    findBestInterface(). See that routine above for * documentation.
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
BoundaryTriangulate::
findBestInterface(
  const set<RandomVariable*> &left_C_l,
  const set<RandomVariable*> &C_l,
  const set<RandomVariable*> &C2,
  const set<RandomVariable*> &C2_1,
  const set<RandomVariable*> &C3,
  set< set<RandomVariable*> >& setset,
  set<RandomVariable*> &best_left_C_l,
  set<RandomVariable*> &best_C_l,
  vector<float>& best_score,
  const vector<BoundaryHeuristic>& bnd_heur_v,
  // --------------------------------------------------------------
  // The next 7 input arguments are used only with the optimal
  // interface algorithm when the IH_MIN_MAX_C_CLIQUE or
  // IH_MIN_MAX_CLIQUE heuristics are used:
  const TriangulateHeuristics& tri_heur,
  const set<RandomVariable*>& P_u1,
  const set<RandomVariable*>& C1_u1,
  const set<RandomVariable*>& Cextra_u1,
  const set<RandomVariable*>& C2_u1,
  const set<RandomVariable*>& E_u1,
  // these next 2 should be const, but there is no "op[] const"
  map < RandomVariable*, RandomVariable* >& C2_u2_to_C1_u1,
  map < RandomVariable*, RandomVariable* >& C2_u2_to_C2_u1
)
{
  if (timer && timer->Expired()) {
    infoMsg(IM::Low, "Exiting boundary algorithm because time has expired\n");
    return;
  }
  boundaryRecursionDepth++;
  set<RandomVariable*>::iterator v;  // vertex
  // consider all v in the current C_l as candidates to be moved
  // left. Essentially, we advance the "partition boundary" to the right
  // by taking nodes that are right of the boundary and moving those
  // nodes (one by one) to the left of the boundary. This is done
  // recursively, and memoization is employed to ensure that we don't
  // redundantly explore the same partition boundary.
  // In all cases, 
  //  1) the boundary is a set of edges (and is not explicitly
  //     represented in the code).
  //  2) C_l are the nodes adjacent and to and directly to the right of the current boundary 
  //  3) left_C_l are the nodes to the left of the current boundary (i.e., they
  //     are former interface nodes that have been shifted left accross the boundary).

  // All nodes in C1 must be to the left of the boundary (ensured by
  // start condition and by the nature of the algoritm).  All nodes to
  // the right of the boundary must still be in C2, or otherwise
  // we would try to produce a boundary with > M chunks. This could
  // lead to invalid networks since the left interface structure in E
  // might not be the same as the left interface structure in C (and this
  // is the reason for unrolling u2 by M+1, to get an extra chunk at the
  // beginning and end of of the chunks in which we do a boundary search).
  // This condition is ensured by "condition ***" below.

  // If C2 consists of multiple chunks (say C2_1, C2_2, etc.) then we
  // should never have that all of C2_1 lies completely to the left of
  // the boundary, since this would be redundant (i.e., a boundary
  // between C1 and C2_1 would be the same edges shifted in time as
  // the boundary between C2_1 and C2_2). This is ensured by
  // "condition ###" below where further comments on this point are
  // given.
  for (v = C_l.begin(); v != C_l.end(); v ++) {
    // TODO: 
    //      Rather than "for all nodes in C_l", we could do a random
    //      subset of nodes in C_l to speed this up if it takes too
    //      long for certain graphs.  Alternatively, a greedy strategy
    //      could be employed where rather than recursing for each v,
    //      we could first choose the best v and then recurse on that.
    //      But note that this is only run once per graph so it will
    //      be beneficial to do this since its cost would probably be
    //      ammortized over the many runs of inference with the graph.

    // continue with probabilty 'boundaryTraverseFraction'
    if (!rnd.coin(boundaryTraverseFraction))
      continue;

    set<RandomVariable*> res;

    // Condition ***: if v has neighbors in C3 (via set intersection),
    // then continue since if v was moved left, we would end up with a
    // boundary > M chunks. This would be invalid since there
    // would be a node left of C_l that connects directly to the right
    // of C_l. This might also lead to problems since we are not guaranteed
    // that the basic left interface in E is the same as the basic left interface
    // in C.
    set_intersection((*v)->neighbors.begin(),
		     (*v)->neighbors.end(),
		     C3.begin(),C3.end(),
		     inserter(res, res.end()));
    if (res.size() != 0)
      continue;

    // Then it is ok to remove v from C_l and move v to the left.
    // take v from C_l and place it in left_C_l
    set<RandomVariable*> next_left_C_l = left_C_l;
    next_left_C_l.insert((*v));

    // Only do this check if C2_1 is non-empty. If it is empty, we
    // assume the check is not needed (say because C2 is itself only
    // one chunk wide).
    if (C2_1.size() > 0) {
      // Condition ###: next check to make sure that we haven't used
      // up all the nodes in C2_1, because if we did we would be in
      // the (M - 1) case (i.e., the only other possible boundaries
      // would span M-1 chunks, skipping the first one). Any
      // boundaries that we return thus will always have a left
      // interface with at least one node in the first chunk of C2, which
      // we call C2_1. Note that it is possible to return a boundary
      // that spans fewer than M chunks (which means that the optimal
      // interface would have been found for a smaller M).

      // We ensure this condition by making sure that C2_1 is not a
      // proper subset of (left_C_l U {v}) = next_left_C_l.

      set<RandomVariable*> tmp;      
      set_difference(C2_1.begin(),C2_1.end(),
		     next_left_C_l.begin(),next_left_C_l.end(),
		     inserter(tmp,tmp.end()));
      if (tmp.size() == 0)
	continue;
    }

    if (message(VerbosityLevels(High))) {
      printf("  Recursion depth = %d\n",boundaryRecursionDepth);
      printf("  Moving node %s(%d) over boundary\n",
	     (*v)->name().c_str(),(*v)->frame());
    }

    // Next, create the new left interface, next_C_l by adding all
    // neighbors of v that are in C3 U C2\next_left_C_l to an
    // initially empty next_C_l.
    set<RandomVariable*> next_C_l;
    set<RandomVariable*> tmp;

    // TODO: remove this next check and fix comments above, since the
    // check is redundant now.
    // add neighbors of v that are in C3
    set_intersection((*v)->neighbors.begin(),
		     (*v)->neighbors.end(),
		     C3.begin(),C3.end(),
		     inserter(tmp,tmp.end()));
    assert ( tmp.size() == 0 );

    // add neighbors of v that are in C2
    set_intersection((*v)->neighbors.begin(),
		     (*v)->neighbors.end(),
		     C2.begin(),C2.end(),
		     inserter(tmp,tmp.end()));

    // remove nodes that are in next_left_C_l
    res.clear();    
    set_difference(tmp.begin(),tmp.end(),
		   next_left_C_l.begin(),next_left_C_l.end(),
		   inserter(res,res.end()));
    // add this to the previous left interface
    set_union(res.begin(),res.end(),
	      C_l.begin(),C_l.end(),
	      inserter(next_C_l,next_C_l.end()));
    // and lastly, remove v from the new left interface since that was
    // added just above.
    next_C_l.erase((*v));

    // check if we've seen this left interface already.
    if (!noBoundaryMemoize) {
      if (setset.find(next_C_l) != setset.end()) {
	if (message(VerbosityLevels(High))) {
	  printf("  Continuing: found existing interface:");
	  set<RandomVariable*>::iterator i;    
	  for (i=next_C_l.begin();i!=next_C_l.end();i++) {
	    printf(" %s(%d)",
		   (*i)->name().c_str(),
		   (*i)->frame());
	  }	
	  printf("\n");
	}
	continue; // check if memoized, if so, no need to go further.
      }

      // memoize
      setset.insert(next_C_l);
    }

    // score the new left interface
    vector<float> next_score;


    interfaceScore(bnd_heur_v,next_C_l,next_left_C_l,
		   tri_heur,
		   P_u1,C1_u1,Cextra_u1,C2_u1,E_u1,
		   C2_u2_to_C1_u1,C2_u2_to_C2_u1,
		   next_score);
    // check size of candiate interface, and keep a copy of it if it
    // is strictly less then the one we have seen so far.
    if (next_score < best_score) {
      best_left_C_l = next_left_C_l;
      best_C_l = next_C_l;
      best_score = next_score;
    } 

    // recurse
    findBestInterface(next_left_C_l,
		      next_C_l,
		      C2,C2_1,C3,setset,
		      best_left_C_l,best_C_l,best_score,
		      bnd_heur_v,
		      tri_heur,
		      P_u1,C1_u1,Cextra_u1,C2_u1,E_u1,
		      C2_u2_to_C1_u1,C2_u2_to_C2_u1);

  }

  boundaryRecursionDepth--;
}

