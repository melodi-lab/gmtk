/*-
 * GMTK_GMTemplate.cc
 *     manipulations of a GM template
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
////////////////////////////////////////////////////////////////////
//        Main Partition Routines
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

/*-
 *-----------------------------------------------------------------------
 * GMTemplate::findPartitions()
 *  Create the three partitions (P,C,E) of the template using
 *  the heuristics supplied.
 *     fh = a string with triangulation heuristics to use (in order)
 *     flr = force the use of either left or right, rather than use min.
 *     findBestFace = T/F if to use the exponential face finding alg.
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
GMTemplate::
findPartitions(const string& fh,  // face quality heuristic
	       const string& flr, // force left or right
	       const string& th, // triangualtion heuristic
	       const bool findBestFace,
	       set<RandomVariable*>& Pc,
	       set<RandomVariable*>& Cc,
	       set<RandomVariable*>& Ec)
{
  const int debug = 0;

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
  int start_index_of_C1_u2 = -1;
  int start_index_of_C2_u2 = -1;
  int start_index_of_C3_u2 = -1;
  for (unsigned i=0;i<unroll2_rvs.size();i++) {
    if (unroll2_rvs[i]->frame() < firstChunkFrame)
      P_u2.insert(unroll2_rvs[i]);
    else if (unroll2_rvs[i]->frame() <= lastChunkFrame) {
      C1_u2.insert(unroll2_rvs[i]);
      if (start_index_of_C1_u2 == -1)
	start_index_of_C1_u2 = i;
    } else if (unroll2_rvs[i]->frame() <= lastChunkFrame+chunkNumFrames) {
      C2_u2.insert(unroll2_rvs[i]);
      if (start_index_of_C2_u2 == -1)
	start_index_of_C2_u2 = i;
    } else if (unroll2_rvs[i]->frame() <= lastChunkFrame+2*chunkNumFrames) {
      C3_u2.insert(unroll2_rvs[i]);
      if (start_index_of_C3_u2 == -1)
	start_index_of_C3_u2 = i;
    } else 
      E_u2.insert(unroll2_rvs[i]);
  }

  assert (C1_u2.size() == C2_u2.size());
  assert (C2_u2.size() == C3_u2.size());
  if (debug > 0) 
    printf("Size of (P,C1,C2,C3,E) = (%d,%d,%d,%d,%d)\n",
	   P_u2.size(),C1_u2.size(),C2_u2.size(),C3_u2.size(),E_u2.size());


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
  set<RandomVariable*> left_C_l_u2C2;
  set<RandomVariable*> C_l_u2C2;
  set<RandomVariable*> right_C_r_u2C2;
  set<RandomVariable*> C_r_u2C2;

  vector<InterfaceHeuristic> fh_v;
  createVectorInterfaceHeuristic(fh,fh_v);

  vector<TriangulateHeuristic> th_v;
  createVectorTriHeuristic(th,th_v);

  if (flr.size() == 0 || (flr.size() > 0 && toupper(flr[0]) != 'R')) {
    // find best left interface
    findBestInterface(C1_u2,C2_u2,C3_u2,
		      left_C_l_u2C2,C_l_u2C2,fh_v,
		      findBestFace,
		      // find best face args
		      th_v,
		      P_u1,
		      C1_u1,
		      C2_u1,
		      E_u1,
		      C2_u2_to_C1_u1,
		      C2_u2_to_C2_u1
		      );
  }
  if (flr.size() == 0 || (flr.size() > 0 && toupper(flr[0]) != 'L')) {
    // find best right interface
    findBestInterface(C3_u2,C2_u2,C1_u2,
		      right_C_r_u2C2,C_r_u2C2,fh_v,
		      findBestFace,
		      // find best face args
		      th_v,
		      E_u1,
		      C2_u1,
		      C1_u1,
		      P_u1,
		      C2_u2_to_C2_u1,
		      C2_u2_to_C1_u1
		      );
  }

  if ((flr.size() == 0 && C_l_u2C2.size() <= C_r_u2C2.size())
      || 
      (flr.size() > 0 && toupper(flr[0]) == 'L')) {
    // this next routine gives us the best left interface that
    // exists from within the chunk C2_u2 and places
    // it in C_l_u2, and everything to the 'left' of C_l_u2
    // that still lies within C2_u2 is placed in left_C_l_u2
    findInterfacePartitions(P_u1,
			    C1_u1,
			    C2_u1,
			    E_u1,
			    C2_u2_to_C1_u1,
			    C2_u2_to_C2_u1,
			    left_C_l_u2C2,
			    C_l_u2C2,
			    Pc,
			    Cc,
			    Ec);
  } else {
    // find right interface partitions
    findInterfacePartitions(E_u1,
			    C2_u1,
			    C1_u1,
			    P_u1,
			    C2_u2_to_C2_u1,
			    C2_u2_to_C1_u1,
			    right_C_r_u2C2,
			    C_r_u2C2,
			    Ec,
			    Cc,
			    Pc);
  }
}




/*-
 *-----------------------------------------------------------------------
 * GMTemplate::findPartitions()
 *  Create the three partitions (P,C,E) of the template using
 *  the given information stored in the input file.
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
GMTemplate::
findPartitions(iDataStreamFile& is,
	       set<RandomVariable*>& Pc,
	       set<RandomVariable*>& Cc,
	       set<RandomVariable*>& Ec)
{
  vector <RandomVariable*> unroll1_rvs;
  map < RVInfo::rvParent, unsigned > positions;
  fp.unroll(1,unroll1_rvs,positions);

  // need to moralize.
  for (unsigned i=0;i<unroll1_rvs.size();i++) {
    unroll1_rvs[i]->createNeighborsFromParentsChildren();
  }
  for (unsigned i=0;i<unroll1_rvs.size();i++) {
    unroll1_rvs[i]->moralize();    
  }

  unsigned setSize;
  set<RandomVariable*> P;
  set<RandomVariable*> C;
  set<RandomVariable*> E;

  is.read(setSize,"set size");
  for (unsigned i=0;i<setSize;i++) {
    RVInfo::rvParent par;
    is.read(par.first,"parent name");
    is.read(par.second,"parent position");

    map < RVInfo::rvParent, unsigned >::iterator loc;
    loc = positions.find(par);
    if (loc == positions.end())
      error("ERROR: P partition information in file %s is invalid for given graph structure\n",
	    is.fileName());
    P.insert(unroll1_rvs[(*loc).second]);
  }

  is.read(setSize,"set size");
  for (unsigned i=0;i<setSize;i++) {
    RVInfo::rvParent par;
    is.read(par.first,"parent name");
    is.read(par.second,"parent position");

    map < RVInfo::rvParent, unsigned >::iterator loc;
    loc = positions.find(par);
    if (loc == positions.end())
      error("ERROR: C partition information in file %s is invalid for given graph structure\n",
	    is.fileName());
    C.insert(unroll1_rvs[(*loc).second]);
  }


  is.read(setSize,"set size");
  for (unsigned i=0;i<setSize;i++) {
    RVInfo::rvParent par;
    is.read(par.first,"parent name");
    is.read(par.second,"parent position");

    map < RVInfo::rvParent, unsigned >::iterator loc;
    loc = positions.find(par);
    if (loc == positions.end())
      error("ERROR: E partition information in file %s is invalid for given graph structure\n",
	    is.fileName());
    E.insert(unroll1_rvs[(*loc).second]);
  }

  clone(P,Pc);
  clone(C,Cc);
  clone(E,Ec);

}



/*-
 *-----------------------------------------------------------------------
 * GMTemplate::storePartitions()
 *  Store the given argument partitions into a file for later
 *  retreival. This routine writes out the information both
 *  in a more human readable format (as comments preceeded by
 *  a coment character) and in machine readable form to be read
 *  in again (lines that do not begin with comment characters).
 *
 * Preconditions:
 *   Each partition must corresond to a valid and separte GM. Each
 *   variable in each GM must have a valid parent and neighbor members
 *   and the parents/neighbors must only point to other members of a
 *   given partition.
 * 
 * Postconditions:
 *   Each of the partitions have been printed
 *
 * Side Effects:
 *   none, other than changing the file pointer of os
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
GMTemplate::
storePartitions(oDataStreamFile& os,
		const set<RandomVariable*>& Pc,
		const set<RandomVariable*>& Cc,
		const set<RandomVariable*>& Ec)
{

  string buffer;
  char buff[2048];

  // First write it out in human readable form as a comment.
  os.nl();
  os.writeComment("---\n");
  os.writeComment("--- P partition information: variables and their neighbors\n");
  buffer.clear();
  for (set<RandomVariable*>::iterator i=Pc.begin();
       i != Pc.end(); i++) {
    RandomVariable* rv = (*i);
    sprintf(buff,"%s(%d) :",rv->name().c_str(),rv->frame());
    buffer = buff;
    for (set<RandomVariable*>::iterator j=rv->neighbors.begin();
	 j != rv->neighbors.end(); j++) {
      sprintf(buff," %s(%d),",
		      (*j)->name().c_str(),(*j)->frame());
      buffer += buff;
    }
    os.writeComment("%s\n",buffer.c_str());
  }
  // Then write it out in machine readable form not as a comment
  os.writeComment("--- P partition definition\n");
  os.write(Pc.size());
  for (set<RandomVariable*>::iterator i = Pc.begin();
       i != Pc.end(); i++) {
    RandomVariable* rv = (*i);
    os.write(rv->name().c_str(),"rv name");
    os.write(rv->frame(),"rv frame");
  }
  os.nl();

  // First write it out in human readable form as a comment.
  os.nl();
  os.writeComment("---\n");
  os.writeComment("--- C partition information: variables and their neighbors\n");
  buffer.clear();
  for (set<RandomVariable*>::iterator i=Cc.begin();
       i != Cc.end(); i++) {
    RandomVariable* rv = (*i);
    sprintf(buff,"%s(%d) :",rv->name().c_str(),rv->frame());
    buffer = buff;
    for (set<RandomVariable*>::iterator j=rv->neighbors.begin();
	 j != rv->neighbors.end(); j++) {
      sprintf(buff," %s(%d),",
		      (*j)->name().c_str(),(*j)->frame());
      buffer += buff;
    }
    os.writeComment("%s\n",buffer.c_str());
  }
  // Then write it out in machine readable form not as a comment
  os.writeComment("--- C partition definition\n");
  os.write(Cc.size());
  for (set<RandomVariable*>::iterator i = Cc.begin();
       i != Cc.end(); i++) {
    RandomVariable* rv = (*i);
    os.write(rv->name().c_str(),"rv name");
    os.write(rv->frame(),"rv frame");
  }
  os.nl();


  // First write it out in human readable form as a comment.
  os.nl();
  os.writeComment("---\n");
  os.writeComment("--- E partition information: variables and their neighbors\n");
  buffer.clear();
  for (set<RandomVariable*>::iterator i=Ec.begin();
       i != Ec.end(); i++) {
    RandomVariable* rv = (*i);
    sprintf(buff,"%s(%d) :",rv->name().c_str(),rv->frame());
    buffer = buff;
    for (set<RandomVariable*>::iterator j=rv->neighbors.begin();
	 j != rv->neighbors.end(); j++) {
      sprintf(buff," %s(%d),",
		      (*j)->name().c_str(),(*j)->frame());
      buffer += buff;
    }
    os.writeComment("%s\n",buffer.c_str());
  }
  // Then write it out in machine readable form not as a comment
  os.writeComment("--- E partition definition\n");
  os.write(Ec.size());
  for (set<RandomVariable*>::iterator i = Ec.begin();
       i != Ec.end(); i++) {
    RandomVariable* rv = (*i);
    os.write(rv->name().c_str(),"rv name");
    os.write(rv->frame(),"rv frame");
  }
  os.nl();

}


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Main Triangulation Routines
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

/*-
 *-----------------------------------------------------------------------
 * GMTemplate::triangulatePartitions()
 *  Given the argument partitions, triangulate each
 *  one using the given heuristics.
 *
 * Preconditions:
 *   Each partition must corresond to a valid and separte GM. Each
 *   variable in each GM must have a valid parent and neighbor members
 *   and the parents/neighbors must only point to other members of a
 *   given partition.
 * 
 * Postconditions:
 *   Each of the partitions have been triangulated.
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
GMTemplate::
triangulatePartitions(const string& th,
		      set<RandomVariable*>& P,
		      set<RandomVariable*>& C,
		      set<RandomVariable*>& E,
		      vector<MaxClique>& Pcliques,
		      vector<MaxClique>& Ccliques,
		      vector<MaxClique>& Ecliques,
		      vector<RandomVariable*>& Pordered,
		      vector<RandomVariable*>& Cordered,
		      vector<RandomVariable*>& Eordered)
{
  vector<TriangulateHeuristic> th_v;
  createVectorTriHeuristic(th,th_v);
  triangulatePartitions(th_v,P,C,E,Pcliques,Ccliques,Ecliques,
			Pordered,Cordered,Eordered);
}
void
GMTemplate::
triangulatePartitions(const vector<TriangulateHeuristic> th_v,
		      set<RandomVariable*>& P,
		      set<RandomVariable*>& C,
		      set<RandomVariable*>& E,
		      vector<MaxClique>& Pcliques,
		      vector<MaxClique>& Ccliques,
		      vector<MaxClique>& Ecliques,
		      vector<RandomVariable*>& Pordered,
		      vector<RandomVariable*>& Cordered,
		      vector<RandomVariable*>& Eordered)
{

  basicTriangulate(P,th_v,
		   Pordered,Pcliques);
  basicTriangulate(C,th_v,
		   Cordered,Ccliques);
  basicTriangulate(E,th_v,
		   Eordered,Ecliques);
}



/*-
 *-----------------------------------------------------------------------
 * GMTemplate::triangulatePartitions()
 *  Given the argument partitions and the information in the files, triangulate each
 *  one using the given heuristics.
 *
 * Preconditions:
 *   Each partition must corresond to a valid and separte GM. Each
 *   variable in each GM must have a valid parent and neighbor members
 *   and the parents/neighbors must only point to other members of a
 *   given partition.
 * 
 * Postconditions:
 *   Each of the partitions have been triangulated.
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
GMTemplate::
triangulatePartitions(iDataStreamFile& is,
		      set<RandomVariable*>& P,
		      set<RandomVariable*>& C,
		      set<RandomVariable*>& E,
		      vector<MaxClique>& Pcliques,
		      vector<MaxClique>& Ccliques,
		      vector<MaxClique>& Ecliques,
		      vector<RandomVariable*>& Pordered,
		      vector<RandomVariable*>& Cordered,
		      vector<RandomVariable*>& Eordered)
{
  basicTriangulate(is,P,Pordered,Pcliques);
  basicTriangulate(is,C,Cordered,Ccliques);
  basicTriangulate(is,E,Eordered,Ecliques);
}





/*-
 *-----------------------------------------------------------------------
 * GMTemplate::storePartitionTriangulation()
 *   Given a set of partitions and a triangulation (elimination
 *   ordering) given by the arguments, store that
 *   information in file pointed to by 'os'
 *
 *   There are two versions of this routine, the first one also writes
 *   out cliques that result from using this elimination order
 *   (using the standard algorithm where a complete set is 
 *    taken to be a maxclique if it is not a subset of a previously
 *    identified clique). It writes out the cliques as comments,
 *    meaning that lines are preceeded by a comment character.
 *   The second version of the routine does not write any out
 *   clique information (as it is not given cliques to write out).
 *
 * Preconditions:
 * 
 * Postconditions:
 *
 * Side Effects:
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
GMTemplate::
storePartitionTriangulation(oDataStreamFile& os,
			    const set<RandomVariable*>& P,
			    const set<RandomVariable*>& C,
			    const set<RandomVariable*>& E,
			    const vector<MaxClique>& Pcliques,
			    const vector<MaxClique>& Ccliques,
			    const vector<MaxClique>& Ecliques,
			    const vector<RandomVariable*>& Pordered,
			    const vector<RandomVariable*>& Cordered,
			    const vector<RandomVariable*>& Eordered)
{

  // First write out the cliques for the user
  os.nl();
  os.writeComment("---- P Partitions Cliques and their weights\n");
  float maxWeight = -1.0;
  for (unsigned i=0;i<Pcliques.size();i++) {
    float curWeight = computeWeight(Pcliques[i].nodes);
    if (curWeight > maxWeight) maxWeight = curWeight;
    os.writeComment("%d : %d  %f\n",
		    i,
		    Pcliques[i].nodes.size(),curWeight);
    for (set<RandomVariable*>::iterator j=Pcliques[i].nodes.begin();
	 j != Pcliques[i].nodes.end(); j++) {
      RandomVariable* rv = (*j);
      os.writeComment("   %s(%d)\n",rv->name().c_str(),rv->frame());
    }
  }
  os.writeComment("Maximum clique state space = 1e%f\n",maxWeight);
  // Then write out the elimination order that achieves these
  // cliques.
  os.writeComment("P Parition Elimination Order\n");
  os.write(Pordered.size());
  for (unsigned int i = 0; i < Pordered.size(); i++) {
    RandomVariable* rv = Pordered[i];
    os.write(rv->name().c_str(),"rv name");
    os.write(rv->frame(),"rv frame");
  }
  os.nl();



  // First write out the cliques for the user
  os.nl();
  os.writeComment("---- C Partitions Cliques and their weights\n");
  maxWeight = -1.0;
  for (unsigned i=0;i<Ccliques.size();i++) {
    float curWeight = computeWeight(Ccliques[i].nodes);
    if (curWeight > maxWeight) maxWeight = curWeight;
    os.writeComment("%d : %d  %f\n",
		    i,
		    Ccliques[i].nodes.size(),curWeight);
    for (set<RandomVariable*>::iterator j=Ccliques[i].nodes.begin();
	 j != Ccliques[i].nodes.end(); j++) {
      RandomVariable* rv = (*j);
      os.writeComment("   %s(%d)\n",rv->name().c_str(),rv->frame());
    }
  }
  os.writeComment("Maximum clique state space = 1e%f\n",maxWeight);
  // Then write out the elimination order that achieves these
  // cliques.
  os.writeComment("C Parition Elimination Order\n");
  os.write(Cordered.size());
  for (unsigned int i = 0; i < Cordered.size(); i++) {
    RandomVariable* rv = Cordered[i];
    os.write(rv->name().c_str(),"rv name");
    os.write(rv->frame(),"rv frame");
  }
  os.nl();



  // First write out the cliques for the user
  os.nl();
  os.writeComment("---- E Partitions Cliques and their weights\n");
  maxWeight = -1.0;
  for (unsigned i=0;i<Ecliques.size();i++) {
    float curWeight = computeWeight(Ecliques[i].nodes);
    if (curWeight > maxWeight) maxWeight = curWeight;
    os.writeComment("%d : %d  %f\n",
		    i,
		    Ecliques[i].nodes.size(),curWeight);
    for (set<RandomVariable*>::iterator j=Ecliques[i].nodes.begin();
	 j != Ecliques[i].nodes.end(); j++) {
      RandomVariable* rv = (*j);
      os.writeComment("   %s(%d)\n",rv->name().c_str(),rv->frame());
    }
  }
  os.writeComment("Maximum clique state space = 1e%f\n",maxWeight);
  // Then write out the elimination order that achieves these
  // cliques.
  os.writeComment("E Parition Elimination Order\n");
  os.write(Eordered.size());
  for (unsigned int i = 0; i < Eordered.size(); i++) {
    RandomVariable* rv = Eordered[i];
    os.write(rv->name().c_str(),"rv name");
    os.write(rv->frame(),"rv frame");
  }
  os.nl();


#if 0
  os.writeComment("C Parition Elimination Order\n");
  os.write(Cordered.size());
  for (unsigned int i = 0; i < Cordered.size(); i++) {
    RandomVariable* rv = Cordered[i];
    os.write(rv->name().c_str(),"rv name");
    os.write(rv->frame(),"rv frame");
  }
  os.nl();


  os.writeComment("E Parition Elimination Order\n");
  os.write(Eordered.size());
  for (unsigned int i = 0; i < Eordered.size(); i++) {
    RandomVariable* rv = Eordered[i];
    os.write(rv->name().c_str(),"rv name");
    os.write(rv->frame(),"rv frame");
  }
  os.nl();
#endif

}
void
GMTemplate::
storePartitionTriangulation(oDataStreamFile& os,
			    const set<RandomVariable*>& P,
			    const set<RandomVariable*>& C,
			    const set<RandomVariable*>& E,
			    const vector<RandomVariable*>& Pordered,
			    const vector<RandomVariable*>& Cordered,
			    const vector<RandomVariable*>& Eordered)
{
  // Write out the elimination orders
  os.nl();
  os.writeComment("P Parition Elimination Order\n");
  os.write(Pordered.size());
  for (unsigned int i = 0; i < Pordered.size(); i++) {
    RandomVariable* rv = Pordered[i];
    os.write(rv->name().c_str(),"rv name");
    os.write(rv->frame(),"rv frame");
  }
  os.nl();


  os.writeComment("C Parition Elimination Order\n");
  os.write(Cordered.size());
  for (unsigned int i = 0; i < Cordered.size(); i++) {
    RandomVariable* rv = Cordered[i];
    os.write(rv->name().c_str(),"rv name");
    os.write(rv->frame(),"rv frame");
  }
  os.nl();


  os.writeComment("E Parition Elimination Order\n");
  os.write(Eordered.size());
  for (unsigned int i = 0; i < Eordered.size(); i++) {
    RandomVariable* rv = Eordered[i];
    os.write(rv->name().c_str(),"rv name");
    os.write(rv->frame(),"rv frame");
  }
  os.nl();
}





/*-
 *-----------------------------------------------------------------------
 * GMTemplate::unrollAndTriangulate(th,int)
 *  Just unroll the graph a given number of times and triangulate it
 *  using the heuristics, reporting the results.
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
unrollAndTriangulate(const string& th,  // triangulate heuristics
		     const unsigned numTimes)
{
  vector<TriangulateHeuristic> th_v;
  createVectorTriHeuristic(th,th_v);

  if (numTimes > 0) {
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
    vector<RandomVariable*> ordered;
    basicTriangulate(rvsSet,th_v,
		     ordered,cliques);
    // TODO: just print out for now. Ultimately
    // return the cliques and do inference with them.
    printf("Cliques from graph unrolled %d times\n",numTimes);
    for (unsigned i=0;i<cliques.size();i++) {
      printf("%d : %d  %f\n",i,
	     cliques[i].nodes.size(),
	     computeWeight(cliques[i].nodes));
      for (set<RandomVariable*>::iterator j=cliques[i].nodes.begin();
	   j != cliques[i].nodes.end(); j++) {
	RandomVariable* rv = (*j);
	printf("   %s(%d)\n",rv->name().c_str(),rv->frame());
      }
    }
    printf("\n");
  }
}


/*-
 *-----------------------------------------------------------------------
 * GMTemplate::basicTriangulate()
 *   Triangulate a set of nodes using any combination of a number
 *   if different (but simple) triangulation heuristics such as
 *   (min weight, min size, min fill, etc.).
 *   For a good description of these heuristics, see
 *      D. Rose et. al, 1970, 1976
 *    
 * Preconditions:
 *   Graph must be a valid undirected model. This means
 *   if the graph was originally a directed model, it must
 *   have been properly moralized and their 'neighbors' structure
 *   is valid. It is also assumed that
 *   the parents of each r.v. are valid but that they
 *   only poiint to variables within the set 'nodes' (i.e.,
 *   parents must not point out of this set).
 *
 * Postconditions:
 *   Resulting graph is now triangulated.
 *
 * Side Effects:
 *   Might (and probably will unless graph is already triangulated
 *   and you get lucky by having the found elimination order is perfect) 
 *   change neighbors members of variables.
 *
 * Results:
 *     none
 *
 *
 *-----------------------------------------------------------------------
 */
void
GMTemplate::
basicTriangulate(// input: nodes to triangulate
		 set<RandomVariable*> nodes,
		 // input: triangulation heuristic
		 const vector<TriangulateHeuristic>& th_v,
		 // output: nodes ordered according to resulting elimination
		 vector<RandomVariable*>& orderedNodes,  
		 // output: resulting max cliques
		 vector<MaxClique>& cliques
		 )
{

  const int debug = 0;
  const unsigned num_nodes = nodes.size();

  if (debug > 0)
    printf("\nBEGINNING TRIANGULATION --- \n");

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
  // where 'X' is the combined weight heuristic can be accessed
  // immediately in the front of the map. Also, whenever a node gets
  // eliminated, *ONLY* its neighbor's weights are recalculated. With
  // this data structure it is possible and efficient to do so. This
  // is because a multimap (used to simulate a priority queue) has the
  // ability to remove stuff from the middle.
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

      if (debug > 0)
	printf("TR: computing weight of node %s(%d)\n",
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

	const TriangulateHeuristic th = th_v[thi];

	if (th == TH_MIN_WEIGHT) {
	  float tmp_weight = computeWeight(activeNeighbors,(*i));
	  weight.push_back(tmp_weight);
	  if (debug > 0)
	    printf("  node has weight = %f\n",tmp_weight);
	} else if (th == TH_MIN_FILLIN) {
	  int fill_in = computeFillIn(activeNeighbors);
	  weight.push_back((float)fill_in);
	  if (debug > 0)
	    printf("  node has fill_in = %d\n",fill_in);
	} else if (th == TH_MIN_TIMEFRAME) {
	  weight.push_back((*i)->frame());
	  if (debug > 0)
	    printf("  node has time frame = %d\n",(*i)->frame());
	} else if (th == TH_MIN_SIZE) {
	  weight.push_back((float)activeNeighbors.size());
	  if (debug > 0)
	    printf("  node has active neighbor size = %d\n",
		   activeNeighbors.size());
	} else if (th == TH_MIN_POSITION_IN_FILE) {
	  weight.push_back((float)(*i)->rv_info.variablePositionInStrFile);
	  if (debug > 0)
	    printf("  node has position in file = %d\n",
		   (*i)->rv_info.variablePositionInStrFile);
	} else if (th == TH_MIN_HINT) {
	  weight.push_back((float)(*i)->rv_info.eliminationOrderHint);
	  if (debug > 0)
	    printf("  node has elimination order hint = %f\n",
		   (*i)->rv_info.eliminationOrderHint);
	} else
	  warning("Warning: unimplemented triangulation heuristic (ignored)\n");
      }

      pair< vector<float>,RandomVariable*> p(weight,(*i));
      rv2unNodesMap[(*i)] = (unorderedNodes.insert(p));
    }

    if (debug > 0) {
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


    // go through the updated multi-map, 
    // grab an iterator pair to iterate between
    // all nodes that have the lowest weight. This utilizes
    // the fact that the multimap stores values in ascending
    // order based on key (in this case the weight), and so
    // the first one (i.e., mm.begin() ) should have the lowest weight.
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

    if (debug > 0) {
      printf("\nEliminating node %s(%d) with weights:",
	     rv->name().c_str(),rv->frame());
      for (unsigned l = 0; l < (*(ip.first)).first.size(); l++ )
	printf(" %f,",(*(ip.first)).first[l]);
      printf("\n");
    }


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

    // insert node into ordered list
    orderedNodes.push_back(rv);

    // insert node into ordered set
    orderedNodesSet.insert(rv);

    // erase node from priority queue
    unorderedNodes.erase(ip.first);

    // only update not-yet-eliminated nodes that could possibly
    // have been effected by the current node 'rv' being
    // eliminated. I.e., create the set of active neighbors
    // of rv.
    nodesToUpdate.clear();
    set_difference(rv->neighbors.begin(),rv->neighbors.end(),
		   orderedNodesSet.begin(),orderedNodesSet.end(),
		   inserter(nodesToUpdate,nodesToUpdate.end()));

    // erase active neighbors of nodes since they
    // will need to be recomputed above.
    for (set<RandomVariable*>::iterator n = nodesToUpdate.begin();
	 n != nodesToUpdate.end();
	 n++) {
      unorderedNodes.erase(rv2unNodesMap[(*n)]);
    }

    // continue until all nodes are eliminated.
  } while (orderedNodesSet.size() < num_nodes);

}






/*-
 *-----------------------------------------------------------------------
 * GMTemplate::basicTriangulate()
 *   Triangulate a set of nodes using the elimination order information 
 *   given in file at the current position.
 *    
 * Preconditions:
 *   Graph must be a valid undirected model. This means
 *   if the graph was originally a directed model, it must
 *   have been properly moralized and their 'neighbors' structure
 *   is valid. It is also assumed that
 *   the parents of each r.v. are valid but that they
 *   only point to variables within the set 'nodes' (i.e.,
 *   parents must not point out of this set).
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
 *
 *-----------------------------------------------------------------------
 */
void
GMTemplate::
basicTriangulate(iDataStreamFile& is,
		 set<RandomVariable*> nodes,
		 // output: nodes ordered according to resulting elimination
		 vector<RandomVariable*>& orderedNodes,  
		 // output: resulting max cliques
		 vector<MaxClique>& cliques
		 )
{
  cliques.clear();
  // orderedNodes == already eliminated nodes. 
  orderedNodes.clear();
  
  // create a map for easy access to set of nodes
  map < RVInfo::rvParent, RandomVariable* > namePos2Var;
  for (set<RandomVariable*>::iterator i=nodes.begin();
       i != nodes.end(); i++) {
    RandomVariable* rv = (*i);
    RVInfo::rvParent par;
    par.first = rv->name();
    par.second = rv->frame();
    namePos2Var[par] = rv;
  }

  // now read in the information.
  unsigned setSize;
  is.read(setSize,"set size");
  orderedNodes.resize(setSize);
  for (unsigned i=0;i<setSize;i++) {
    RVInfo::rvParent par;
    is.read(par.first,"parent name");
    is.read(par.second,"parent position");

    map < RVInfo::rvParent, RandomVariable* >::iterator loc;
    loc = namePos2Var.find(par);
    if (loc == namePos2Var.end())
	error("ERROR: elimination order list in file %s is not valid for given structure file.\n",
	      is.fileName());

    RandomVariable* rv = (*loc).second;
    orderedNodes[i] = rv;
  }

  // Also keep ordered (eliminated) nodes as a set for easy
  // intersection, with other node sets.
  set<RandomVariable*> orderedNodesSet;

  // finally triangulate and make cliques
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



////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        General Support Routines for Triangulation              //
////////////////////////////////////////////////////////////////////
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
GMTemplate::dropEdgeDirections(vector <RandomVariable*>& rvs)
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
GMTemplate::moralize(vector <RandomVariable*>& rvs)
{
  for (unsigned i=0;i<rvs.size();i++) {
    rvs[i]->moralize();    
  }
  return true;
}



/*-
 *-----------------------------------------------------------------------
 * GMTemplate::computeWeight()
 *   Computes the log base 10 weight of a set of nodes (i.e.,
 *   the union of 'node' and 'nodes', ignores 'node' if 'node == NULL').
 *   
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
float
GMTemplate::
computeWeight(const set<RandomVariable*>& nodes,
	      const RandomVariable* node) 
	      
{
  // compute weight in log base 10 so that
  //   1) we don't overflow
  //   2) base 10 is an easy to understand magnitude rating of state space.

  float tmp_weight = 0;
  // First get cardinality of 'node', but if
  // it is continuous or observed, it does not change the weight.
  // TODO: The assumption here (for now) is that all continuous variables
  // are observed. This will change in a future version.
  if (node != NULL) {
    if (node->discrete && node->hidden) {
      DiscreteRandomVariable *const drv = (DiscreteRandomVariable*)node;
      // weight changes only if node is not deterministic (Lauritzen CG inference).
      if (drv->deterministic()) {
	// then there is a possibility that this node
	// does not affect the state space, as long
	// as all of this nodes parents are in the clique.
	bool truly_deterministic = true;
	for (unsigned i=0;i<drv->allPossibleParents.size();i++) {
	  if (nodes.find(drv->allPossibleParents[i]) == nodes.end()) {
	    // found a parent that is not in 'node' set so the
	    // node would not truly be deterministic here.
	    truly_deterministic = false;
	    break;
	  }
	}
	if (!truly_deterministic)
	  tmp_weight += log10((double)drv->useCardinality());	
      } else
	tmp_weight += log10((double)drv->useCardinality());
    }
  }
  // Next, get weight of all 'nodes'
  for (set<RandomVariable*>::iterator j=nodes.begin();
       j != nodes.end();
       j++) {
    RandomVariable *const rv = (*j);
    // First get cardinality of 'node', but if
    // it is continuous or observed, it does not change the weight.
    // TODO: The assumption here (for now) is that all continuous variables
    // are observed. This will change in a future version (Lauritzen CG inference).
    if (rv->discrete && rv->hidden) {
      DiscreteRandomVariable *const drv = (DiscreteRandomVariable*)rv;
      if (drv->deterministic()) {
	// then there is a possibility that this node
	// does not affect the state space, as long
	// as all of this nodes parents are in the clique.
	bool truly_deterministic = true;
	for (unsigned i=0;i<drv->allPossibleParents.size();i++) {
	  if (nodes.find(drv->allPossibleParents[i]) == nodes.end()) {
	    // found a parent that is not in 'node' set so the
	    // node would not truly be deterministic here.
	    truly_deterministic = false;
	    break;
	  }
	}
	if (!truly_deterministic)
	  tmp_weight += log10((double)drv->useCardinality());	
      } else
	tmp_weight += log10((double)drv->useCardinality());
    }
  }
  return tmp_weight;
}


/*-
 *-----------------------------------------------------------------------
 * GMTemplate::computeFillIn()
 *   Computes the number of edges that would need to be added
 *   among 'nodes' to make 'nodes' complete.
 *   
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
GMTemplate::
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

    /*
      printf("    Neighbor %s(%d) needs %d for fill in\n",
      (*j)->name().c_str(),(*j)->frame(),
      (nodes.size() - 1 - tmp.size()));
      printf("    Neighbors of %s(%d) are:",
      (*j)->name().c_str(),(*j)->frame());
      for (set<RandomVariable*>::iterator k=(*j)->neighbors.begin();
      k != (*j)->neighbors.end();
      k++) {
      printf("%s(%d) ",
      (*k)->name().c_str(),(*k)->frame());
      }
      printf("\n");
    */

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
 * GMTemplate::variableSetScore()
 *      Compute the 'score' of a set of variables, where the 
 *      score is based on one or more of the set-based triangulation
 *      heuristics.
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
GMTemplate::variableSetScore(const vector<InterfaceHeuristic>& fh_v,
			     const set<RandomVariable*>& varSet,
			     vector<float>& score)
{
  const int debug=0;

  score.clear();
  for (unsigned fhi=0;fhi<fh_v.size();fhi++) {
    const InterfaceHeuristic fh = fh_v[fhi];
    if (fh == IH_MIN_WEIGHT) {
      float tmp_weight = computeWeight(varSet);
      score.push_back(tmp_weight);
      if (debug > 0)
	printf("  set has weight = %f\n",tmp_weight);
    } else if (fh == IH_MIN_FILLIN) {
      int fill_in = computeFillIn(varSet);
      score.push_back((float)fill_in);
      if (debug > 0)
	printf("  set has fill_in = %d\n",fill_in);
    } else if (fh == IH_MIN_SIZE) {
      score.push_back((float)varSet.size());
      if (debug > 0)
	printf("  set has size = %d\n",
	       varSet.size());
    } else
      warning("Warning: invalid variable set score given. Ignored\n");
  }

  return;
}


/*-
 *-----------------------------------------------------------------------
 * GMTemplate::interfaceScore()
 *      Compute the 'score' of a set of variables that are to be
 *      a candidate interface. In the best of cases, it is a
 *      simple easy score based on the vars size, fill-in, or weight.
 *      In the "worst" of cases, the score is based on an entire
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
GMTemplate::interfaceScore(
 const vector<InterfaceHeuristic>& fh_v,
 const set<RandomVariable*>& C_l,
 // more input variables for use when MIN_CLIQUE is active
 const set<RandomVariable*>& left_C_l,
 const vector<TriangulateHeuristic>& th_v,
 const set<RandomVariable*>& P_u1,
 const set<RandomVariable*>& C1_u1,
 const set<RandomVariable*>& C2_u1,
 const set<RandomVariable*>& E_u1,
 // these next 2 should be const, but there is no "op[] const"
 map < RandomVariable*, RandomVariable* >& C2_u2_to_C1_u1,
 map < RandomVariable*, RandomVariable* >& C2_u2_to_C2_u1,
 // output score
 vector<float>& score)
{
  const int debug=1;
  score.clear();
  for (unsigned fhi=0;fhi<fh_v.size();fhi++) {
    const InterfaceHeuristic fh = fh_v[fhi];
    if (fh == IH_MIN_WEIGHT) {
      float tmp_weight = computeWeight(C_l);
      score.push_back(tmp_weight);
      if (debug > 0)
	printf("  set has weight = %f\n",tmp_weight);
    } else if (fh == IH_MIN_FILLIN) {
      int fill_in = computeFillIn(C_l);
      score.push_back((float)fill_in);
      if (debug > 0)
	printf("  set has fill_in = %d\n",fill_in);
    } else if (fh == IH_MIN_SIZE) {
      score.push_back((float)C_l.size());
      if (debug > 0)
	printf("  set has size = %d\n",
	       C_l.size());
    } else if (fh == IH_MIN_MAX_CLIQUE) {
      // This is the expensive one, need to form a set of partitions,
      // given the current interface, triangulate that partition set, and then
      // compute the score of the worst best clique, and fill the score
      // variable above with this worst scoring clique (i.e., 
      // we find the interface that has the best worst-case performance.

      set<RandomVariable*> Pc;
      set<RandomVariable*> Cc;
      set<RandomVariable*> Ec;
      vector<MaxClique> Pcliques;
      vector<MaxClique> Ccliques;
      vector<MaxClique> Ecliques;
      vector<RandomVariable*> Pordered;
      vector<RandomVariable*> Cordered;
      vector<RandomVariable*> Eordered;

      findInterfacePartitions(P_u1,
			      C1_u1,
			      C2_u1,
			      E_u1,
			      C2_u2_to_C1_u1,
			      C2_u2_to_C2_u1,
			      left_C_l,
			      C_l,
			      Pc,
			      Cc,
			      Ec);

      triangulatePartitions(th_v,
			    Pc,Cc,Ec,
			    Pcliques,Ccliques,Ecliques,
			    Pordered,Cordered,Eordered);

      // Now got cliques compute worst score using
      // the weight of a clique as the score mechanism.
      float maxWeight = -1.0;
      for (unsigned i=0;i<Pcliques.size();i++) {
	float curWeight = computeWeight(Pcliques[i].nodes);
	printf("   --- P curWeight = %f\n",curWeight);
	if (curWeight > maxWeight) maxWeight = curWeight;
      }
      for (unsigned i=0;i<Ccliques.size();i++) {
	float curWeight = computeWeight(Ccliques[i].nodes);
	printf("   --- C curWeight = %f\n",curWeight);
	if (curWeight > maxWeight) maxWeight = curWeight;
      }
      for (unsigned i=0;i<Ecliques.size();i++) {
	float curWeight = computeWeight(Ecliques[i].nodes);
	printf("   --- E curWeight = %f\n",curWeight);
	if (curWeight > maxWeight) maxWeight = curWeight;
      }
      score.push_back(maxWeight);

      if (debug > 0)
	printf("  set has max clique weight = %f\n",maxWeight);

      deleteNodes(Pc);
      deleteNodes(Cc);
      deleteNodes(Ec);      

    }  else
      warning("Warning: invalid variable set score given. Ignored\n");
  }

  return;
}



/*-
 *-----------------------------------------------------------------------
 * GMTemplate::clone()
 *   Clone a set of random variables from 'in' to 'out'. The cloned
 *   variables have parents, children, and neighbors consisting
 *   of only the members that are in the corresponding 'in' set (i.e.,
 *   if any parents, children, neighbors, in 'in' pointed to variables
 *   outside of 'in', then the corresponding variables in 'out' do
 *   not contain those parents,children,neighbors. 
 *
 * Preconditions:
 *   'in' is a set of random variables to be cloned.
 *
 * Postconditions:
 *   'out' is a clone of 'in' but parents,children,neighbors only
 *   point to members within the set 'out'
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
clone(const set<RandomVariable*>& in, 
      set<RandomVariable*>& out)
{
  
  map < RandomVariable*, RandomVariable* > in_to_out;
  for (set<RandomVariable*>::iterator i=in.begin();
       i != in.end(); i++) {

    // sanity check
    assert ( (*i)->neighbors.find((*i)) == (*i)->neighbors.end() );

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

    RandomVariable*rv = (*i);

    // first set up new neighbors for in_to_out[rv]
    set<RandomVariable*> tmp;
    for (set<RandomVariable*>::iterator j = rv->neighbors.begin();
	 j != rv->neighbors.end(); j++) {    
      if (in_to_out.find((*j)) != in_to_out.end()) {
	// then it is included in this set.
	tmp.insert(in_to_out[(*j)]);
      }
    }
    in_to_out[rv]->neighbors = tmp;
    // assertion to make sure that no node has itself as neighbor.
    assert( in_to_out[rv]->neighbors.find(in_to_out[rv])
	    == in_to_out[rv]->neighbors.end() );

    // next, set new sparents for in_to_out[rv]
    vector<RandomVariable *> sParents;
    for (unsigned l=0;l<rv->switchingParents.size();l++) {
      if (in_to_out.find(rv->switchingParents[l]) != in_to_out.end())
	sParents.push_back(in_to_out[rv->switchingParents[l]]);
    }

    // next, set conditional parents
    vector< vector < RandomVariable* > > cParentsList;
    cParentsList.resize(rv->conditionalParentsList.size());
    for (unsigned l=0;l<rv->conditionalParentsList.size();l++) {
      for (unsigned m=0;m<rv->conditionalParentsList[l].size();m++) {
	if (in_to_out.find(rv->conditionalParentsList[l][m]) 
	    != in_to_out.end())
	  cParentsList[l].push_back(	
              in_to_out[rv->conditionalParentsList[l][m]]
	      );
      }
    }
    in_to_out[rv]->setParents(sParents,cParentsList);
  }
}



/*-
 *-----------------------------------------------------------------------
 * GMTemplate::deleteNodes()
 *   Given a set of random variables, delete the nodes pointed
 *   to by the set. 
 *
 * Preconditions:
 *   'nodes' is a valid set of node pointers
 *
 * Postconditions:
 *   all RV*'s in 'nodes' have been deleted. The set should thereafter
 *   immediately be deleted or filled with new nodes
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
deleteNodes(const set<RandomVariable*>& nodes)
{
  for (set<RandomVariable*>::iterator i = nodes.begin();
       i != nodes.end(); i++) 
    delete (*i);
}




/*-
 *-----------------------------------------------------------------------
 * GMTemplate::createVectorTriHeuristic()
 *      create a vector of triangluation heuristics based
 *      on a string that is passed in.
 *
 * Preconditions:
 *      - String should contain set of heuristcs, see code
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
GMTemplate::createVectorTriHeuristic(const string& th,
				     vector<TriangulateHeuristic>& th_v)
{
  if (th.size() == 0) {
    // default case.
    // first by weight
    th_v.push_back(TH_MIN_WEIGHT); 
    // then by fill in if weight in tie
    th_v.push_back(TH_MIN_FILLIN); 
    // and lastly by time frame (earliest first)
    // (but note that this is not valid to judge a face)
    th_v.push_back(TH_MIN_TIMEFRAME); 
  } else {
    for (unsigned i=0;i<th.size();i++) {
      switch (th[i]) {
      case 'S':
	th_v.push_back(TH_MIN_SIZE);
	break;
      case 'T':
	th_v.push_back(TH_MIN_TIMEFRAME);
	break;
      case 'F':
	th_v.push_back(TH_MIN_FILLIN);
	break;
      case 'W':
	th_v.push_back(TH_MIN_WEIGHT);
	break;
      case 'E':
	th_v.push_back(TH_MIN_ENTROPY);
	break;
      case 'P':
	th_v.push_back(TH_MIN_POSITION_IN_FILE);
	break;
      case 'H':
	th_v.push_back(TH_MIN_HINT);
	break;
      default:
	error("ERROR: Unknown triangulation heuristic given '%c' in string '%s'\n",
	      th[i],th.c_str());
	break;
      }
    }
  }
  return;
}



/*-
 *-----------------------------------------------------------------------
 * GMTemplate::createVectorInterfaceHeuristic()
 *      create a vector of interface heuristics based
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
GMTemplate::createVectorInterfaceHeuristic(const string& fh,
					   vector<InterfaceHeuristic>& fh_v)
{
  if (fh.size() == 0) {
    // default case.
    // first by weight
    fh_v.push_back(IH_MIN_WEIGHT); 
    // then by fill in if weight in tie
    fh_v.push_back(IH_MIN_FILLIN); 
  } else {
    for (unsigned i=0;i<fh.size();i++) {
      switch (fh[i]) {
      case 'S':
	fh_v.push_back(IH_MIN_SIZE);
	break;
      case 'F':
	fh_v.push_back(IH_MIN_FILLIN);
	break;
      case 'W':
	fh_v.push_back(IH_MIN_WEIGHT);
	break;
      case 'E':
	fh_v.push_back(IH_MIN_ENTROPY);
	break;
      case 'C':
	fh_v.push_back(IH_MIN_MAX_CLIQUE);
	break;
      default:
	error("ERROR: Unknown triangulation heuristic given '%c' in string '%s'\n",
	      fh[i],fh.c_str());
	break;
      }
    }
  }
  return;
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
    // make sure self is not its own neighbor
    res.erase((*i));
    (*i)->neighbors = res;
  }

  return;
}


////////////////////////////////////////////////////////////////////
//        Support Routines for P,C,E Interface Computation
////////////////////////////////////////////////////////////////////



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
GMTemplate::findBestLeftInterface(
           const set<RandomVariable*> &C1,
	   const set<RandomVariable*> &C2,
	   const set<RandomVariable*> &C3,



	   set<RandomVariable*> &left_C_l,
	   set<RandomVariable*> &C_l,
	   const vector<InterfaceHeuristic>& fh_v,
	   const bool recurse)
{
  const int debug = 1;

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

  vector<float> best_score;
  variableSetScore(fh_v,C_l,best_score);

  if (debug > 0) {
    printf("Size of basic left interface C_l = %d\n",C_l.size());
    printf("Score of basic left interface C_l =");
    for (unsigned i=0;i<best_score.size();i++)
      printf(" %f ",best_score[i]);
    printf("\n");
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
  }

  // start recursion to find the truly best interface.
  if (recurse) {
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
			  best_C_l,
			  best_score,fh_v);
    if (debug > 0) {
      printf("Size of best left interface = %d\n",best_C_l.size());
      printf("Score of best left interface =");
      for (unsigned i=0;i<best_score.size();i++)
	printf(" %f ",best_score[i]);
      printf("\n");
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
    }
    left_C_l = best_left_C_l;
    C_l = best_C_l;
  }


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
findBestLeftInterface(const set<RandomVariable*> &left_C_l,
		      const set<RandomVariable*> &C_l,
		      const set<RandomVariable*> &C2,
		      const set<RandomVariable*> &C3,
		      set< set<RandomVariable*> >& setset,
		      set<RandomVariable*> &best_left_C_l,
		      set<RandomVariable*> &best_C_l,
		      vector<float>& best_score,
		      const vector<InterfaceHeuristic>& fh_v)
{
  set<RandomVariable*>::iterator v;  // vertex
  // consider all v in the current C_l as candidates
  // to be moved left.
  for (v = C_l.begin(); v != C_l.end(); v ++) {
    // TODO: rather than "for all nodes in C_l", we could
    // do a random subset of nodes in C_l to speed this up if
    // it takes too long. But note that this is only run once
    // per graph so it will be beneficial to do this since
    // its cost might be ammortized over the many runs of the graph

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

    // memoize
    setset.insert(next_C_l);


    vector<float> next_score;
    variableSetScore(fh_v,next_C_l,next_score);
    // check size of candiate interface, and keep a copy of
    // it if it is less then the one we have seen so far.
    if (next_score < best_score) {
      best_left_C_l = next_left_C_l;
      best_C_l = next_C_l;
      best_score = next_score;
    } 

    findBestLeftInterface(next_left_C_l,
			  next_C_l,
			  C2,C3,setset,
			  best_left_C_l,best_C_l,best_score,
			  fh_v);

  }
}


/*-
 *-----------------------------------------------------------------------
 * GMTemplate::findLeftInterfacePartitions()
 *   Create the three partitions and triangulate them. This
 *   routine is essentialy the left-interface specific portion
 *   of routine triangulatePartitions(). See that routine
 *   for more documentation.
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
findLeftInterfacePartitions(
 // input variables
 const set<RandomVariable*>& P_u1,
 const set<RandomVariable*>& C1_u1,
 const set<RandomVariable*>& C2_u1,
 const set<RandomVariable*>& E_u1,
 // these next 2 should be const, but there is no "op[] const"
 map < RandomVariable*, RandomVariable* >& C2_u2_to_C1_u1,
 map < RandomVariable*, RandomVariable* >& C2_u2_to_C2_u1,
 const set<RandomVariable*>& left_C_l_u2C2,
 const set<RandomVariable*>& C_l_u2C2,
 // output variables
 set<RandomVariable*>& Pc,
 set<RandomVariable*>& Cc,
 set<RandomVariable*>& Ec
)
{

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
   *      alternatively (and easier): for snake, just use the unconstrained
   *           triangulation method (which works perfectly for snake).
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
  
  // finally, create the modified sets P, C, and E
  // where P = modified prologue
  // where C = modified chunk to repeat
  // where E = modified epilogue to repeat
  // which are to be triangulated separately.
  //  P = P' + C1'(left_C_l) + C1'(C_l)
  //  C = C1'\C1'(left_C_l) + C2'(left_C_l) + C2'(C_l)
  //  E = C2'\C2'(left_C_l) + E'

  set<RandomVariable*> P = P_u1;
  set<RandomVariable*> C;
  set<RandomVariable*> E = E_u1;

  // Finish P
  set_union(left_C_l_u1C1.begin(),left_C_l_u1C1.end(),
	    C_l_u1C1.begin(),C_l_u1C1.end(),
	    inserter(P,P.end()));

  // C
  set_union(left_C_l_u1C2.begin(),left_C_l_u1C2.end(),
	    C_l_u1C2.begin(),C_l_u1C2.end(),
	    set_difference(C1_u1.begin(),C1_u1.end(),
			   left_C_l_u1C1.begin(),left_C_l_u1C1.end(),
			   inserter(C,C.end())));
  // finish E
  set_difference(C2_u1.begin(),C2_u1.end(),
		 left_C_l_u1C2.begin(),left_C_l_u1C2.end(),
		 inserter(E,E.end()));
		 
#if 0  
  printf("---\nSet P is:\n");
  for (set<RandomVariable*>::iterator i=P.begin();
       i != P.end(); i++) {
    RandomVariable* rv = (*i);
    printf("%s(%d) :",rv->name().c_str(),rv->frame());
    for (set<RandomVariable*>::iterator j=rv->neighbors.begin();
	 j != rv->neighbors.end(); j++) {
      printf(" %s(%d),",
	     (*j)->name().c_str(),(*j)->frame());

    }
    printf("\n");
  }

  printf("---\nSet C is:\n");
  for (set<RandomVariable*>::iterator i=C.begin();
       i != C.end(); i++) {
    RandomVariable* rv = (*i);
    printf("%s(%d) :",rv->name().c_str(),rv->frame());
    for (set<RandomVariable*>::iterator j=rv->neighbors.begin();
	 j != rv->neighbors.end(); j++) {
      printf(" %s(%d),",
	     (*j)->name().c_str(),(*j)->frame());

    }
    printf("\n");
  }


  printf("---\nSet E is:\n");
  for (set<RandomVariable*>::iterator i=E.begin();
       i != E.end(); i++) {
    RandomVariable* rv = (*i);
    printf("%s(%d) :",rv->name().c_str(),rv->frame());
    for (set<RandomVariable*>::iterator j=rv->neighbors.begin();
	 j != rv->neighbors.end(); j++) {
      printf(" %s(%d),",
	     (*j)->name().c_str(),(*j)->frame());

    }
    printf("\n");
  }
#endif

  clone(P,Pc);
  clone(C,Cc);
  clone(E,Ec);

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
				   set<RandomVariable*> &C_r,
				   const vector<InterfaceHeuristic>& fh_v,
				   const bool recurse)
{
  const int debug = 1;

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

  vector<float> best_score;
  variableSetScore(fh_v,C_r,best_score);

  if (debug > 0) {
    printf("Size of basic right interface C_r = %d\n",C_r.size());
    printf("Score of basic right interface C_r =");
    for (unsigned i=0;i<best_score.size();i++)
      printf(" %f ",best_score[i]);
    printf("\n");
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
  }


  // start recursion to find the truly best interface.
  if (recurse) {
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
			   best_C_r,
			   best_score,fh_v);

    if (debug > 0) {
      printf("Size of best right interface = %d\n",best_C_r.size());
      printf("Score of best right interface =");
      for (unsigned i=0;i<best_score.size();i++)
	printf(" %f ",best_score[i]);
      printf("\n");
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
    }

    right_C_r = best_right_C_r;
    C_r = best_C_r;
  }

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
findBestRightInterface(const set<RandomVariable*> &right_C_r,
		       const set<RandomVariable*> &C_r,
		       const set<RandomVariable*> &C2,
		       const set<RandomVariable*> &C1,
		       set< set<RandomVariable*> >& setset,
		       set<RandomVariable*> &best_right_C_r,
		       set<RandomVariable*> &best_C_r,
		       vector<float>& best_score,
		       const vector<InterfaceHeuristic>& fh_v)
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

    // memoize
    setset.insert(next_C_r);

    vector<float> next_score;
    variableSetScore(fh_v,next_C_r,next_score);
    // check size of candiate interface, and keep a copy of
    // it if it is less then the one we have seen so far.
    if (next_score < best_score) {
      best_right_C_r = next_right_C_r;
      best_C_r = next_C_r;
      best_score = next_score;
    } 

    findBestLeftInterface(next_right_C_r,
			  next_C_r,
			  C2,C1,setset,
			  best_right_C_r,best_C_r,best_score,
			  fh_v);

  }
}



/*-
 *-----------------------------------------------------------------------
 * GMTemplate::findRightInterfacePartitions()
 *   Create the three partitions and triangulate them. This
 *   routine is essentialy the left-interface specific portion
 *   of routine triangulatePartitions(). See that routine
 *   for more documentation.
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
findRightInterfacePartitions(
 // input variables
 const set<RandomVariable*>& P_u1,
 const set<RandomVariable*>& C1_u1,
 const set<RandomVariable*>& C2_u1,
 const set<RandomVariable*>& E_u1,
 map < RandomVariable*, RandomVariable* >& C2_u2_to_C1_u1,
 map < RandomVariable*, RandomVariable* >& C2_u2_to_C2_u1,
 const set<RandomVariable*>& right_C_r_u2C2,
 const set<RandomVariable*>& C_r_u2C2,
 // output variables
 set<RandomVariable*>& Pc,
 set<RandomVariable*>& Cc,
 set<RandomVariable*>& Ec
)
{

  // now we need to make a bunch of sets to be unioned
  // together to get the partitions.
  set<RandomVariable*> C_r_u1C1;
  set<RandomVariable*> C_r_u1C2;
  for (set<RandomVariable*>::iterator i = C_r_u2C2.begin();
       i!= C_r_u2C2.end(); i++) {

    C_r_u1C1.insert(C2_u2_to_C1_u1[(*i)]);
    C_r_u1C2.insert(C2_u2_to_C2_u1[(*i)]);					  
  }

  // These (C_r_u1C1 and C_r_u1C2) are the interfaces which are forced
  // to be complete (i.e., part of a maxclique).
  // TODO: these next steps will ruin the clique_size = 2 property
  //     of the 'snake' structure. The todo is to get this working 
  //     with that (and similar) structures.
  /*
   *     idea: make complete components in C_r *only* if they
   *           are connected via nodes/edges within either (P + left_C_r)
   *           or preceeding C.
   *           What we will then have is a collection of cliques for 
   *           the interface(s). In this case, we glue together
   *           the corresponding sets of cliques. Right interface
   *           algorithm should be similar (and use E rather than P).
   *           But might both left and right interface need to be used in the
   *           same repeated chunk in this case to get clique_size=2???
   *          
   */
  makeComplete(C_r_u1C1);
  makeComplete(C_r_u1C2);

  set<RandomVariable*> right_C_r_u1C1;
  set<RandomVariable*> right_C_r_u1C2;
  for (set<RandomVariable*>::iterator i = right_C_r_u2C2.begin();
       i != right_C_r_u2C2.end(); i++) {

    right_C_r_u1C1.insert(C2_u2_to_C1_u1[(*i)]);
    right_C_r_u1C2.insert(C2_u2_to_C2_u1[(*i)]);					  
  }
  
  // finally, create the modified sets P, C, and E
  // where P = modified prologue
  // where C = modified chunk to repeat
  // where E = modified epilogue to repeat
  // which are to be triangulated separately.
  //  P = P' + C1'\C1'(right_C_r)
  //  C = C1'(C_r) + C1'(right_C_r) + C2'\C2'(right_C_r)
  //  E = C2'(C_r) + C2'(right_C_l) + E'


  set<RandomVariable*> P = P_u1;
  set<RandomVariable*> C;
  set<RandomVariable*> E = E_u1;

  // Finish E
  set_union(right_C_r_u1C2.begin(),right_C_r_u1C2.end(),
	    C_r_u1C2.begin(),C_r_u1C2.end(),
	    inserter(E,E.end()));

  // C
  set_union(right_C_r_u1C1.begin(),right_C_r_u1C1.end(),
	    C_r_u1C1.begin(),C_r_u1C1.end(),
	    set_difference(C2_u1.begin(),C2_u1.end(),
			   right_C_r_u1C2.begin(),right_C_r_u1C2.end(),
			   inserter(C,C.end())));
  // finish P
  set_difference(C1_u1.begin(),C1_u1.end(),
		 right_C_r_u1C1.begin(),right_C_r_u1C1.end(),
		 inserter(P,P.end()));
		 

#if 0
  printf("---\nSet P is:\n");
  for (set<RandomVariable*>::iterator i=P.begin();
       i != P.end(); i++) {
    RandomVariable* rv = (*i);
    printf("%s(%d) :",rv->name().c_str(),rv->frame());
    for (set<RandomVariable*>::iterator j=rv->neighbors.begin();
	 j != rv->neighbors.end(); j++) {
      printf(" %s(%d),",
	     (*j)->name().c_str(),(*j)->frame());

    }
    printf("\n");
  }

  printf("---\nSet C is:\n");
  for (set<RandomVariable*>::iterator i=C.begin();
       i != C.end(); i++) {
    RandomVariable* rv = (*i);
    printf("%s(%d) :",rv->name().c_str(),rv->frame());
    for (set<RandomVariable*>::iterator j=rv->neighbors.begin();
	 j != rv->neighbors.end(); j++) {
      printf(" %s(%d),",
	     (*j)->name().c_str(),(*j)->frame());

    }
    printf("\n");
  }


  printf("---\nSet E is:\n");
  for (set<RandomVariable*>::iterator i=E.begin();
       i != E.end(); i++) {
    RandomVariable* rv = (*i);
    printf("%s(%d) :",rv->name().c_str(),rv->frame());
    for (set<RandomVariable*>::iterator j=rv->neighbors.begin();
	 j != rv->neighbors.end(); j++) {
      printf(" %s(%d),",
	     (*j)->name().c_str(),(*j)->frame());

    }
    printf("\n");
  }
#endif

  clone(P,Pc);
  clone(C,Cc);
  clone(E,Ec);

}



/*-
 *-----------------------------------------------------------------------
 * GMTemplate::findBestInterface()
 *
 *  An exponential-time routine to find best left/right interfaces
 *  and corresponding partitions. Since the operations
 *  for finding the best left and right interfaces are symmetric,
 *  the determiniation of if we are searching for the best
 *  left interface or right interface is determined entirely
 *  based on the arguments that are passed into these routines. 
 *  Note that for simplicity, the names have been defined
 *  in terms of the left interface, but that is only for ease of
 *  understanding.
 *
 * For left interface:
 *   Given a twice unrolled graph, P,C1,C2,C3,E
 *   find the best left interface within C2 starting at the
 *   "standard" or initial left interface between C1 and C2.
 *   Note that the routine only uses C1,C2, and C3 and
 *   P and E are not needed.
 *
 * For right interface:
 *   Given a twice unrolled graph, P,C1,C2,C3,E
 *   find the best right interface within C2 starting at the
 *   "standard" or initial right interface between C3 and C2.
 *   Note that the routine only uses C1,C2, and C3 and
 *   P and E are not needed.
 *
 *
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
 *    Put best left (resp. right) interface from C2 into C_l
 *    and place any additional variables from C2 to the left (resp.
 *    right) of C_l into 'left_C_l';
 *
 *
 *-----------------------------------------------------------------------
 */
void
GMTemplate::findBestInterface(
 const set<RandomVariable*> &C1,
 const set<RandomVariable*> &C2,
 const set<RandomVariable*> &C3,
 set<RandomVariable*> &left_C_l,
 set<RandomVariable*> &C_l,
 const vector<InterfaceHeuristic>& fh_v,
 const bool recurse,
 // more input variables
 const vector<TriangulateHeuristic>& th_v,
 const set<RandomVariable*>& P_u1,
 const set<RandomVariable*>& C1_u1,
 const set<RandomVariable*>& C2_u1,
 const set<RandomVariable*>& E_u1,
 // these next 2 should be const, but there is no "op[] const"
 map < RandomVariable*, RandomVariable* >& C2_u2_to_C1_u1,
 map < RandomVariable*, RandomVariable* >& C2_u2_to_C2_u1
 )
{
  const int debug = 1;

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

  vector<float> best_score;
  interfaceScore(fh_v,C_l,left_C_l,
		 th_v,
		 P_u1,C1_u1,C2_u1,E_u1,
		 C2_u2_to_C1_u1,C2_u2_to_C2_u1,
		 best_score);

  if (debug > 0) {
    printf("Size of basic left interface C_l = %d\n",C_l.size());
    printf("Score of basic left interface C_l =");
    for (unsigned i=0;i<best_score.size();i++)
      printf(" %f ",best_score[i]);
    printf("\n");
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
  }

  // start recursion to find the truly best interface.
  if (recurse) {
    // best ones found so far
    set<RandomVariable*> best_left_C_l = left_C_l;
    set<RandomVariable*> best_C_l = C_l;
    set< set<RandomVariable*> > setset;
    findBestInterface(left_C_l,
		      C_l,
		      C2,
		      C3,
		      setset,
		      best_left_C_l,
		      best_C_l,
		      best_score,fh_v,
		      th_v,
		      P_u1,C1_u1,C2_u1,E_u1,
		      C2_u2_to_C1_u1,C2_u2_to_C2_u1);
    if (debug > 0) {
      printf("Size of best left interface = %d\n",best_C_l.size());
      printf("Score of best left interface =");
      for (unsigned i=0;i<best_score.size();i++)
	printf(" %f ",best_score[i]);
      printf("\n");
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
    }
    left_C_l = best_left_C_l;
    C_l = best_C_l;
  }
}

/*-
 *-----------------------------------------------------------------------
 * GMTemplate::findBestInterface()
 *    recursive helper function for the first call findBestInterface()
 *    See that routine above for documentation.
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
findBestInterface(
  const set<RandomVariable*> &left_C_l,
  const set<RandomVariable*> &C_l,
  const set<RandomVariable*> &C2,
  const set<RandomVariable*> &C3,
  set< set<RandomVariable*> >& setset,
  set<RandomVariable*> &best_left_C_l,
  set<RandomVariable*> &best_C_l,
  vector<float>& best_score,
  const vector<InterfaceHeuristic>& fh_v,
  // more input variables
  const vector<TriangulateHeuristic>& th_v,
  const set<RandomVariable*>& P_u1,
  const set<RandomVariable*>& C1_u1,
  const set<RandomVariable*>& C2_u1,
  const set<RandomVariable*>& E_u1,
  // these next 2 should be const, but there is no "op[] const"
  map < RandomVariable*, RandomVariable* >& C2_u2_to_C1_u1,
  map < RandomVariable*, RandomVariable* >& C2_u2_to_C2_u1
)
{
  set<RandomVariable*>::iterator v;  // vertex
  // consider all v in the current C_l as candidates
  // to be moved left.
  for (v = C_l.begin(); v != C_l.end(); v ++) {
    // TODO: rather than "for all nodes in C_l", we could
    // do a random subset of nodes in C_l to speed this up if
    // it takes too long. But note that this is only run once
    // per graph so it will be beneficial to do this since
    // its cost might be ammortized over the many runs of the graph

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

    // memoize
    setset.insert(next_C_l);

    vector<float> next_score;
    interfaceScore(fh_v,next_C_l,next_left_C_l,
		   th_v,
		   P_u1,C1_u1,C2_u1,E_u1,
		   C2_u2_to_C1_u1,C2_u2_to_C2_u1,
		   next_score);

    // check size of candiate interface, and keep a copy of
    // it if it is less then the one we have seen so far.
    if (next_score < best_score) {
      best_left_C_l = next_left_C_l;
      best_C_l = next_C_l;
      best_score = next_score;
    } 

    findBestInterface(next_left_C_l,
		      next_C_l,
		      C2,C3,setset,
		      best_left_C_l,best_C_l,best_score,
		      fh_v,
		      th_v,
		      P_u1,C1_u1,C2_u1,E_u1,
		      C2_u2_to_C1_u1,C2_u2_to_C2_u1);

  }
}

/*-
 *-----------------------------------------------------------------------
 * GMTemplate::findInterfacePartitions()
 *   Create the three partitions, either left or right depending
 *   on the order of the arguments given.
 *
 * For the left interface, we create new P,C, and E variable sets where
 *  where P = modified prologue
 *  where C = modified chunk to repeat
 *  where E = modified epilogue to repeat
 *  which are to be triangulated separately.
 *   P = P' + C1'(left_C_l) + C1'(C_l)
 *   C = C1'\C1'(left_C_l) + C2'(left_C_l) + C2'(C_l)
 *   E = C2'\C2'(left_C_l) + E'
 *
 * For the right interface,  we create new P,C, and E variable sets where
 *   where P = modified prologue
 *   where C = modified chunk to repeat
 *   where E = modified epilogue to repeat
 *   which are to be triangulated separately.
 *    P = P' + C1'\C1'(right_C_r)
 *    C = C1'(C_r) + C1'(right_C_r) + C2'\C2'(right_C_r)
 *    E = C2'(C_r) + C2'(right_C_l) + E'
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
findInterfacePartitions(
 // input variables
 const set<RandomVariable*>& P_u1,
 const set<RandomVariable*>& C1_u1,
 const set<RandomVariable*>& C2_u1,
 const set<RandomVariable*>& E_u1,
 // these next 2 should be const, but there is no "op[] const"
 map < RandomVariable*, RandomVariable* >& C2_u2_to_C1_u1,
 map < RandomVariable*, RandomVariable* >& C2_u2_to_C2_u1,
 const set<RandomVariable*>& left_C_l_u2C2,
 const set<RandomVariable*>& C_l_u2C2,
 // output variables
 set<RandomVariable*>& Pc,
 set<RandomVariable*>& Cc,
 set<RandomVariable*>& Ec
)
{

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
   *      alternatively (and easier): for snake, just use the unconstrained
   *           triangulation method (which works perfectly for snake).
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
  
  // finally, create the modified sets P, C, and E
  // where P = modified prologue
  // where C = modified chunk to repeat
  // where E = modified epilogue to repeat
  // which are to be triangulated separately.
  //  P = P' + C1'(left_C_l) + C1'(C_l)
  //  C = C1'\C1'(left_C_l) + C2'(left_C_l) + C2'(C_l)
  //  E = C2'\C2'(left_C_l) + E'

  set<RandomVariable*> P = P_u1;
  set<RandomVariable*> C;
  set<RandomVariable*> E = E_u1;

  // Finish P
  set_union(left_C_l_u1C1.begin(),left_C_l_u1C1.end(),
	    C_l_u1C1.begin(),C_l_u1C1.end(),
	    inserter(P,P.end()));

  // C
  set_union(left_C_l_u1C2.begin(),left_C_l_u1C2.end(),
	    C_l_u1C2.begin(),C_l_u1C2.end(),
	    set_difference(C1_u1.begin(),C1_u1.end(),
			   left_C_l_u1C1.begin(),left_C_l_u1C1.end(),
			   inserter(C,C.end())));
  // finish E
  set_difference(C2_u1.begin(),C2_u1.end(),
		 left_C_l_u1C2.begin(),left_C_l_u1C2.end(),
		 inserter(E,E.end()));

  clone(P,Pc);
  clone(C,Cc);
  clone(E,Ec);

}



////////////////////////////////////////////////////////////////////
//        Test Driver
////////////////////////////////////////////////////////////////////




#ifdef MAIN


#endif
