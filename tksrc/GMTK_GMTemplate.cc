/*-
 * GMTK_GMTemplate.cc
 *     Basic GM Template and Basic Triangulation Routines
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
#include "debug.h"
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
//        Static variables used in this file
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

static const string P_partition_name("P_PARTITION");
static const string C_partition_name("C_PARTITION");
static const string E_partition_name("E_PARTITION");
static const string PC_interface_name("PC_PARTITION");
static const string CE_interface_name("CE_PARTITION");

static const string P_elimination_order("P_ELIMINATION_ORDER");
static const string C_elimination_order("C_ELIMINATION_ORDER");
static const string E_elimination_order("E_ELIMINATION_ORDER");


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
findPartitions(// face quality heuristic
	       const string& fh,  
	       // force use of only left or right interface
	       const string& flr, 
	       // triangualtion heuristic to use
	       const string& th, 
	       // should we run the exponential find best interface
	       const bool findBestFace,
	       // number of chunks in which to find interface boundary (M>=1)
	       const unsigned M,
	       // the resulting new prologue
	       set<RandomVariable*>& Pc,
	       // the resulting new chunk
	       set<RandomVariable*>& Cc,
	       // the resulting new epilogue
	       set<RandomVariable*>& Ec,
	       // interface between P and C
	       set<RandomVariable*>& PCInterface,
	       // interface between C and E
	       set<RandomVariable*>& CEInterface)
{

  // M = number of chunks in which interface/pipeline algorithm can operate.
  // i.e., a boundary is searched for in M repeated chunks, where M >= 1.
  // Note that M puts a constraints on the number of time frames
  // of the observations (i.e., the utterance length in time frames).
  // Specificaly, if N = number of frames of utterance,
  //  then we must have N = length(P) + length(E) + k*M*length(C)
  //  where k = 2,3,4, ... is some integer >= 2.
  // Therefore, making M larger reduces the number valid possible utterance lengths.
  assert ( M >= 1 );

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
  for (unsigned i=0;i<unroll2_rvs.size();i++) {
    if (unroll2_rvs[i]->frame() < firstChunkFrame)
      P_u2.insert(unroll2_rvs[i]);
    else if (unroll2_rvs[i]->frame() <= lastChunkFrame) {
      C1_u2.insert(unroll2_rvs[i]);
      if (start_index_of_C1_u2 == -1)
	start_index_of_C1_u2 = i;
    } else if (unroll2_rvs[i]->frame() <= lastChunkFrame+M*chunkNumFrames) {
      C2_u2.insert(unroll2_rvs[i]);
      if (start_index_of_C2_u2 == -1)
	start_index_of_C2_u2 = i;
    } else if (unroll2_rvs[i]->frame() <= lastChunkFrame+(M+1)*chunkNumFrames) {
      C3_u2.insert(unroll2_rvs[i]);
      if (start_index_of_C3_u2 == -1)
	start_index_of_C3_u2 = i;
    } else 
      E_u2.insert(unroll2_rvs[i]);
  }

  assert (M*C1_u2.size() == C2_u2.size());
  assert (C2_u2.size() == M*C3_u2.size());
  infoMsg(Low,"Size of (P,C1,C2,C3,E) = (%d,%d,%d,%d,%d)\n",
	  P_u2.size(),C1_u2.size(),C2_u2.size(),C3_u2.size(),E_u2.size());


  // create sets P', C1', C2', and E', from graph unrolled 2*M-1 time(s)
  // each hyper-chunk C1' and C2' are M frames long.
  vector <RandomVariable*> unroll1_rvs;
  fp.unroll(2*M-1,unroll1_rvs);
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
    else if (unroll1_rvs[i]->frame() <= lastChunkFrame + (M-1)*chunkNumFrames) {
      C1_u1.insert(unroll1_rvs[i]);
      if (start_index_of_C1_u1 == -1)
	start_index_of_C1_u1 = i;
    } else if (unroll1_rvs[i]->frame() <= lastChunkFrame+(2*M-1)*chunkNumFrames) {
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

  vector<float> best_L_score;
  vector<float> best_R_score;

  if (flr.size() == 0 || (flr.size() > 0 && toupper(flr[0]) != 'R')) {
    // find best left interface
    infoMsg(Tiny,"---\nFinding BEST LEFT interface\n");

    set<RandomVariable*> C2_1_u2; // first chunk in C2_u2
    if (M == 1) {
      // make it empty signaling that we don't bother to check it in this case.
      C2_1_u2.clear();
    } else {
      for (set<RandomVariable*>::iterator i=C2_u2.begin();
	   i != C2_u2.end();i++) {
	if ((*i)->frame() > lastChunkFrame && (*i)->frame() <= lastChunkFrame+chunkNumFrames)
	  C2_1_u2.insert((*i));
      }
    }
    findBestInterface(C1_u2,C2_u2,C2_1_u2,C3_u2,
		      left_C_l_u2C2,C_l_u2C2,best_L_score,
		      fh_v,
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
    infoMsg(Tiny,"---\nFinding BEST RIGHT interface\n");

    set<RandomVariable*> C2_l_u2; // last chunk of C2_u2
    if (M == 1) {
      // make it empty signaling that we don't bother to check it in this case.
      C2_l_u2.clear();
    } else {
      for (set<RandomVariable*>::iterator i=C2_u2.begin();
	   i != C2_u2.end();i++) {
	if ((*i)->frame() > (lastChunkFrame+(M-1)*chunkNumFrames) && (*i)->frame() <= (lastChunkFrame+M*chunkNumFrames))
	  C2_l_u2.insert((*i));
      }
    }
    findBestInterface(C3_u2,C2_u2,C2_l_u2,C1_u2,
		      right_C_r_u2C2,C_r_u2C2,best_R_score,
		      fh_v,
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
    infoMsg(Nano,"---\nUsing left interface to define partitions\n");
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
			    Ec,
			    PCInterface,
			    CEInterface);
  } else {
    // find right interface partitions
    infoMsg(Nano,"---\nUsing right interface to define partitions\n");
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
			    Pc,
			    CEInterface,
			    PCInterface);
  }
}



/*-
 *-----------------------------------------------------------------------
 * GMTemplate::findPartitions()
 *  Create the three partitions (P,C,E) of the template using the
 *  given information stored in the input file. See
 *  findPartitions() routine above for argument definitions.
 *
 * Preconditions:
 *   Object must be instantiated and have the use of the information
 *   in a valid FileParser object (which stores the parsed structure
 *   file) Arguments must indicate valid heuristics to use.
 *
 * Postconditions:
 *   Arguments Pc, Cc, and Ec now contain partitions in a separate
 *   graph from the FileParser (i.e., file parser information is not
 *   disturbed)
 *
 * Side Effects:
 *   none
 *
 * Results:
 *   Pc, Cc, and Ec as arguments.
 *
 *
 *-----------------------------------------------------------------------
 */
void
GMTemplate::
findPartitions(iDataStreamFile& is,
	       unsigned &M,
	       set<RandomVariable*>& Pc,
	       set<RandomVariable*>& Cc,
	       set<RandomVariable*>& Ec,
	       set<RandomVariable*>& PCIc,
	       set<RandomVariable*>& CEIc)
{


  is.read(M,"M value");
  if (M == 0)
    error("ERROR: M (number of chunks in which to find interface boundary) must be >= 1\n");

  vector <RandomVariable*> unroll1_rvs;
  map < RVInfo::rvParent, unsigned > positions;
  fp.unroll(2*M-1,unroll1_rvs,positions);

  // need to moralize.
  for (unsigned i=0;i<unroll1_rvs.size();i++) {
    unroll1_rvs[i]->createNeighborsFromParentsChildren();
  }
  for (unsigned i=0;i<unroll1_rvs.size();i++) {
    unroll1_rvs[i]->moralize();    
  }

  unsigned setSize;
  string str_tmp;
  set<RandomVariable*> P;
  set<RandomVariable*> C;
  set<RandomVariable*> E;

  is.read(str_tmp,"name");
  if (str_tmp != P_partition_name)
    error("ERROR: P partition information in file %s is invalid for given graph structure\n",is.fileName());
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

  is.read(str_tmp,"name");
  if (str_tmp != C_partition_name)
    error("ERROR: C partition information in file %s is invalid for given graph structure\n",is.fileName());
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

  is.read(str_tmp,"name");
  if (str_tmp != E_partition_name)
    error("ERROR: E partition information in file %s is invalid for given graph structure\n",is.fileName());
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
  
  // next, read in the interface definitions.

  set<RandomVariable*> PCInterface;
  set<RandomVariable*> CEInterface;

  // get PC interface
  is.read(str_tmp,"name");
  if (str_tmp != PC_interface_name)
    error("ERROR: PC interface information in file %s is invalid for given graph structure\n",is.fileName());
  is.read(setSize,"set size");
  for (unsigned i=0;i<setSize;i++) {
    RVInfo::rvParent par;
    is.read(par.first,"parent name");
    is.read(par.second,"parent position");

    map < RVInfo::rvParent, unsigned >::iterator loc;
    loc = positions.find(par);
    if (loc == positions.end())
      error("ERROR: PC interface information in file %s is invalid for given graph structure\n",
	    is.fileName());
    PCInterface.insert(unroll1_rvs[(*loc).second]);
  }

  // get CE interface
  is.read(str_tmp,"name");
  if (str_tmp != CE_interface_name)
    error("ERROR: CE interface information in file %s is invalid for given graph structure\n",is.fileName());
  is.read(setSize,"set size");
  for (unsigned i=0;i<setSize;i++) {
    RVInfo::rvParent par;
    is.read(par.first,"parent name");
    is.read(par.second,"parent position");

    map < RVInfo::rvParent, unsigned >::iterator loc;
    loc = positions.find(par);
    if (loc == positions.end())
      error("ERROR: CE interface information in file %s is invalid for given graph structure\n",
	    is.fileName());
    CEInterface.insert(unroll1_rvs[(*loc).second]);
  }

  // finally create a new variable set for each, make
  // the interfaces complete, and finish up.

  map < RandomVariable*, RandomVariable* > P_in_to_out;
  map < RandomVariable*, RandomVariable* > C_in_to_out;
  map < RandomVariable*, RandomVariable* > E_in_to_out;
  set < RandomVariable* > tmp;

  setUpClonedPartitionGraph(P,C,E,Pc,Cc,Ec,P_in_to_out,C_in_to_out,E_in_to_out);

  tmp.clear();
  for (set<RandomVariable*>::iterator i=PCInterface.begin();
       i != PCInterface.end();i++) {
    tmp.insert(P_in_to_out[(*i)]);
  }
  makeComplete(tmp);
  PCIc = tmp;

  tmp.clear();
  for (set<RandomVariable*>::iterator i=PCInterface.begin();
       i != PCInterface.end();i++) {
    tmp.insert(C_in_to_out[(*i)]);
  }
  makeComplete(tmp);
  tmp.clear();
  for (set<RandomVariable*>::iterator i=CEInterface.begin();
       i != CEInterface.end();i++) {
    tmp.insert(C_in_to_out[(*i)]);
  }
  makeComplete(tmp);  
  CEIc = tmp;

  tmp.clear();
  for (set<RandomVariable*>::iterator i=CEInterface.begin();
       i != CEInterface.end();i++) {
    tmp.insert(E_in_to_out[(*i)]);
  }
  makeComplete(tmp);  
}



/*-
 *-----------------------------------------------------------------------
 * GMTemplate::storePartitions()
 *
 *  Store the given argument partitions into a file for later
 *  retreival. This routine writes out the information both in a more
 *  human readable format (as comments preceeded by a coment character
 *  which includes other useful information) and in machine readable
 *  form to be read in again (lines that do not begin with comment
 *  characters).  The information written includes;
 *
 *    P partition
 *    C partition
 *    E partition
 *    PC interface
 *    CE interface
 *
 * Preconditions:
 *
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
storePartitions(oDataStreamFile& os,const GMInfo& info) 
{

  string buffer;
  char buff[2048];

  // number of chunks in which to find interface boundary
  os.nl();
  os.writeComment("---\n");
  os.writeComment("--- M, number of chunks in which to find interface boundary\n");
  os.write(info.M);
  os.nl();

  // next write it out in human readable form as a comment.
  os.nl();
  os.writeComment("---\n");
  os.writeComment("--- P partition information: variables and their neighbors\n");
  buffer.clear();
  for (set<RandomVariable*>::iterator i=info.P.begin();
       i != info.P.end(); i++) {
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
  os.write(P_partition_name);
  os.write(info.P.size());
  for (set<RandomVariable*>::iterator i = info.P.begin();
       i != info.P.end(); i++) {
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
  for (set<RandomVariable*>::iterator i=info.C.begin();
       i != info.C.end(); i++) {
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
  os.write(C_partition_name);
  os.write(info.C.size());
  for (set<RandomVariable*>::iterator i = info.C.begin();
       i != info.C.end(); i++) {
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
  for (set<RandomVariable*>::iterator i=info.E.begin();
       i != info.E.end(); i++) {
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
  os.write(E_partition_name);
  os.write(info.E.size());
  for (set<RandomVariable*>::iterator i = info.E.begin();
       i != info.E.end(); i++) {
    RandomVariable* rv = (*i);
    os.write(rv->name().c_str(),"rv name");
    os.write(rv->frame(),"rv frame");
  }
  os.nl();


  // First write it out in human readable form as a comment.
  os.nl();
  os.writeComment("---\n");
  os.writeComment("--- PC information : variables and their neighbors\n");
  buffer.clear();
  for (set<RandomVariable*>::iterator i=info.PCInterface.begin();
       i != info.PCInterface.end(); i++) {
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
  os.writeComment("--- PC interface definition\n");
  os.write(PC_interface_name);
  os.write(info.PCInterface.size());
  for (set<RandomVariable*>::iterator i = info.PCInterface.begin();
       i != info.PCInterface.end(); i++) {
    RandomVariable* rv = (*i);
    os.write(rv->name().c_str(),"rv name");
    os.write(rv->frame(),"rv frame");
  }
  os.nl();


  // First write it out in human readable form as a comment.
  os.nl();
  os.writeComment("---\n");
  os.writeComment("--- CE information : variables and their neighbors\n");
  buffer.clear();
  for (set<RandomVariable*>::iterator i=info.CEInterface.begin();
       i != info.CEInterface.end(); i++) {
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
  os.writeComment("--- CE interface definition\n");
  os.write(CE_interface_name);
  os.write(info.CEInterface.size());
  for (set<RandomVariable*>::iterator i = info.CEInterface.begin();
       i != info.CEInterface.end(); i++) {
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
triangulatePartitions(const vector<TriangulateHeuristic>& th_v,
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
 *  Given the argument partitions and the information in the files,
 *  triangulate each one using the given heuristics.
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
  string tmp;

  is.read(tmp,"name");
  if (tmp != P_elimination_order)
    error("ERROR: P elimination order interface information in file %s is invalid for given graph structure\n",is.fileName());
  basicTriangulate(is,P,Pordered,Pcliques);

  is.read(tmp,"name");
  if (tmp != C_elimination_order)
    error("ERROR: C elimination order interface information in file %s is invalid for given graph structure\n",is.fileName());
  basicTriangulate(is,C,Cordered,Ccliques);

  is.read(tmp,"name");
  if (tmp != E_elimination_order)
    error("ERROR: E elimination order interface information in file %s is invalid for given graph structure\n",is.fileName());
  basicTriangulate(is,E,Eordered,Ecliques);
}





/*-
 *-----------------------------------------------------------------------
 * GMTemplate::storePartitionTriangulation()
 *   Given a set of partitions and a triangulation (elimination
 *   ordering) given by the arguments, store that
 *   information in file pointed to by 'os'
v *
 *   There are two versions of this routine, the first one also writes
 *   out cliques that result from using this elimination order (using
 *   the standard algorithm where a complete set is taken to be a
 *   maxclique if it is not a subset of a previously identified
 *   clique). It writes out the cliques as comments, meaning that
 *   lines are preceeded by a comment character.  The second version
 *   of the routine does not write any out clique information (as it
 *   is not given cliques to write out).
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
			    const unsigned M,
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
  os.writeComment("---\n");
  os.writeComment("---- P Partitions Cliques and their weights\n");
  float maxWeight = -1.0;
  float totalWeight = -1.0; // starting flag
  for (unsigned i=0;i<Pcliques.size();i++) {
    float curWeight = computeWeight(Pcliques[i].nodes);
    if (curWeight > maxWeight) maxWeight = curWeight;
    if (totalWeight == -1.0)
      totalWeight = curWeight;
    else
      totalWeight = totalWeight + log10(1+pow(10,curWeight-totalWeight));
    os.writeComment("%d : %d  %f\n",
		    i,
		    Pcliques[i].nodes.size(),curWeight);
    for (set<RandomVariable*>::iterator j=Pcliques[i].nodes.begin();
	 j != Pcliques[i].nodes.end(); j++) {
      RandomVariable* rv = (*j);
      os.writeComment("   %s(%d)\n",rv->name().c_str(),rv->frame());
    }
  }
  os.writeComment("Maximum clique state space = 1e%f, total state space = 1e%f\n",maxWeight,totalWeight);
  // Then write out the elimination order that achieves these
  // cliques.
  os.writeComment("P Parition Elimination Order\n");
  os.write(P_elimination_order);
  os.write(Pordered.size());
  for (unsigned int i = 0; i < Pordered.size(); i++) {
    RandomVariable* rv = Pordered[i];
    os.write(rv->name().c_str(),"rv name");
    os.write(rv->frame(),"rv frame");
  }
  os.nl();



  // First write out the cliques for the user
  os.nl();
  os.writeComment("---\n");
  os.writeComment("---- C Partitions Cliques and their weights\n");
  maxWeight = -1.0;
  totalWeight = -1.0;
  for (unsigned i=0;i<Ccliques.size();i++) {
    float curWeight = computeWeight(Ccliques[i].nodes);
    if (curWeight > maxWeight) maxWeight = curWeight;
    if (totalWeight == -1.0)
      totalWeight = curWeight;
    else
      totalWeight = totalWeight + log10(1+pow(10,curWeight-totalWeight));
    os.writeComment("%d : %d  %f\n",
		    i,
		    Ccliques[i].nodes.size(),curWeight);
    for (set<RandomVariable*>::iterator j=Ccliques[i].nodes.begin();
	 j != Ccliques[i].nodes.end(); j++) {
      RandomVariable* rv = (*j);
      os.writeComment("   %s(%d)\n",rv->name().c_str(),rv->frame());
    }
  }
  os.writeComment("Max clique state space = 1e%f, total Cx%d state space = 1e%f, per-chunk total C state space = 1e%f\n",
		  maxWeight,
		  M,
		  totalWeight,
		  totalWeight - log10((double)M));
  // Then write out the elimination order that achieves these
  // cliques.
  os.writeComment("C Parition Elimination Order\n");
  os.write(C_elimination_order);
  os.write(Cordered.size());
  for (unsigned int i = 0; i < Cordered.size(); i++) {
    RandomVariable* rv = Cordered[i];
    os.write(rv->name().c_str(),"rv name");
    os.write(rv->frame(),"rv frame");
  }
  os.nl();



  // First write out the cliques for the user
  os.nl();
  os.writeComment("---\n");
  os.writeComment("---- E Partitions Cliques and their weights\n");
  maxWeight = -1.0;
  totalWeight = -1.0;
  for (unsigned i=0;i<Ecliques.size();i++) {
    float curWeight = computeWeight(Ecliques[i].nodes);
    if (curWeight > maxWeight) maxWeight = curWeight;
    if (totalWeight == -1.0)
      totalWeight = curWeight;
    else
      totalWeight = totalWeight + log10(1+pow(10,curWeight-totalWeight));
    os.writeComment("%d : %d  %f\n",
		    i,
		    Ecliques[i].nodes.size(),curWeight);
    for (set<RandomVariable*>::iterator j=Ecliques[i].nodes.begin();
	 j != Ecliques[i].nodes.end(); j++) {
      RandomVariable* rv = (*j);
      os.writeComment("   %s(%d)\n",rv->name().c_str(),rv->frame());
    }
  }
  os.writeComment("Maximum clique state space = 1e%f, total state space = 1e%f\n",maxWeight,totalWeight);
  // Then write out the elimination order that achieves these
  // cliques.
  os.writeComment("E Parition Elimination Order\n");
  os.write(E_elimination_order);
  os.write(Eordered.size());
  for (unsigned int i = 0; i < Eordered.size(); i++) {
    RandomVariable* rv = Eordered[i];
    os.write(rv->name().c_str(),"rv name");
    os.write(rv->frame(),"rv frame");
  }
  os.nl();
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
  os.writeComment("---\n");
  os.writeComment("P Parition Elimination Order\n");
  os.write(P_elimination_order);
  os.write(Pordered.size());
  for (unsigned int i = 0; i < Pordered.size(); i++) {
    RandomVariable* rv = Pordered[i];
    os.write(rv->name().c_str(),"rv name");
    os.write(rv->frame(),"rv frame");
  }
  os.nl();

  os.writeComment("---\n");
  os.writeComment("C Parition Elimination Order\n");
  os.write(C_elimination_order);
  os.write(Cordered.size());
  for (unsigned int i = 0; i < Cordered.size(); i++) {
    RandomVariable* rv = Cordered[i];
    os.write(rv->name().c_str(),"rv name");
    os.write(rv->frame(),"rv frame");
  }
  os.nl();

  os.writeComment("---\n");
  os.writeComment("E Parition Elimination Order\n");
  os.write(E_elimination_order);
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
GMTemplate::
unrollAndTriangulate(// triangulate heuristics
		     const string& th,  
		     // number of times it should be unrolled
		     const unsigned numTimes)
{
  vector<TriangulateHeuristic> th_v;
  createVectorTriHeuristic(th,th_v);

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
    vector<RandomVariable*> ordered;
    basicTriangulate(rvsSet,th_v,
		     ordered,cliques);
    // TODO: just print out for now. Ultimately
    // return the cliques and do inference with them.
    unsigned maxSize = 0;
    float maxSizeCliqueWeight;
    float maxWeight = -1.0;
    float totalWeight = -1.0; // starting flag
    unsigned maxWeightCliqueSize;
    printf("Cliques from graph unrolled %d times\n",numTimes);
    for (unsigned i=0;i<cliques.size();i++) {

      float curWeight = computeWeight(cliques[i].nodes);
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
 * GMTemplate::basicTriangulate()
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
 *   Graph must be a valid undirected model. This means if the graph
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
GMTemplate::
basicTriangulate(// input: nodes to triangulate
		 const set<RandomVariable*> nodes,
		 // input: triangulation heuristic
		 const vector<TriangulateHeuristic>& th_v,
		 // output: nodes ordered according to resulting elimination
		 vector<RandomVariable*>& orderedNodes,  
		 // output: resulting max cliques
		 vector<MaxClique>& cliques,
		 // input: find the cliques as well
		 const bool findCliques
		 )
{
  const unsigned num_nodes = nodes.size();



  infoMsg(Huge,"\nBEGINNING TRIANGULATION --- \n");

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

	const TriangulateHeuristic th = th_v[thi];

	if (th == TH_MIN_WEIGHT || th == TH_MIN_WEIGHT_NO_D) {
	  float tmp_weight = computeWeight(activeNeighbors,(*i),
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
}




/*-
 *-----------------------------------------------------------------------
 * GMTemplate::basicTriangulate()
 *   Triangulate a set of nodes using the elimination order
 *   information given in file at the current position.
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
 *
 *-----------------------------------------------------------------------
 */
void
GMTemplate::
basicTriangulate(// input data stream
		 iDataStreamFile& is,
		 // input: nodes to be triangulated
		 const set<RandomVariable*> nodes,
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

  // Triangulate using the given elimination order
  triangulateElimination(nodes, orderedNodes, cliques);

}



/*-
 *-----------------------------------------------------------------------
 * GMTemplate::triangulateElimination()
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
 *
 *-----------------------------------------------------------------------
 */
void
GMTemplate::
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
	      const RandomVariable* node,
	      const bool useDeterminism)
	      
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
      if (useDeterminism && drv->deterministic()) {
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
      if (useDeterminism && drv->deterministic()) {
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
 *      Compute the 'score' of a set of variables, where the score is
 *      based on one or more of the set-based triangulation
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
  score.clear();
  for (unsigned fhi=0;fhi<fh_v.size();fhi++) {
    const InterfaceHeuristic fh = fh_v[fhi];
    if (fh == IH_MIN_WEIGHT || fh == IH_MIN_WEIGHT_NO_D) {
      float tmp_weight = computeWeight(varSet,NULL,
				       (fh == IH_MIN_WEIGHT));
      score.push_back(tmp_weight);
      infoMsg(Huge,"  variableSetScore: set has weight = %f\n",tmp_weight);
    } else if (fh == IH_MIN_FILLIN) {
      int fill_in = computeFillIn(varSet);
      score.push_back((float)fill_in);
      infoMsg(Huge,"  variableSetScore: set has fill_in = %d\n",fill_in);
    } else if (fh == IH_MIN_SIZE) {
      score.push_back((float)varSet.size());
      infoMsg(Huge,"  variableSetScore: set has size = %d\n",
	      varSet.size());
    } else
      warning("Warning: invalid variable set score given. Ignored\n");
  }

  return;
}


/*-
 *-----------------------------------------------------------------------
 * GMTemplate::interfaceScore()
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
GMTemplate::interfaceScore(
 // the interface heuristic used to score the interface
 const vector<InterfaceHeuristic>& fh_v,
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
 const vector<TriangulateHeuristic>& th_v,
 // The network unrolled 1 time
 const set<RandomVariable*>& P_u1,
 const set<RandomVariable*>& C1_u1,
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
  score.clear();
  for (unsigned fhi=0;fhi<fh_v.size();fhi++) {
    const InterfaceHeuristic fh = fh_v[fhi];
    if (fh == IH_MIN_WEIGHT || fh == IH_MIN_WEIGHT_NO_D) {
      float tmp_weight = computeWeight(C_l,NULL,
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
    } else if (fh == IH_MIN_MAX_CLIQUE || fh == IH_MIN_MAX_C_CLIQUE ||
	       fh == IH_MIN_STATE_SPACE || fh == IH_MIN_C_STATE_SPACE) {
      // This is the expensive one, need to form a set of partitions,
      // given the current interface, triangulate that partition set,
      // and then compute the score of the worst clique, and fill the
      // score variable above with this worst scoring clique (i.e., we
      // find the interface that has the best worst-case performance.

      // TODO: some of these values could be cached to speed this
      // up a bit.

      set<RandomVariable*> Pc;
      set<RandomVariable*> Cc;
      set<RandomVariable*> Ec;
      set<RandomVariable*> PCInterface;
      set<RandomVariable*> CEInterface;
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
			      Ec,
			      PCInterface,
			      CEInterface);

      triangulatePartitions(th_v,
			    Pc,Cc,Ec,
			    Pcliques,Ccliques,Ecliques,
			    Pordered,Cordered,Eordered);

      // Now got cliques compute worst score using
      // the weight of a clique as the score mechanism.
      float maxWeight = -1.0;
      float totalWeight = -1.0; // starting flag
      for (unsigned i=0;i<Ccliques.size();i++) {
	float curWeight = computeWeight(Ccliques[i].nodes);
	if (curWeight > maxWeight) maxWeight = curWeight;
	if (totalWeight == -1.0)
	  totalWeight = curWeight;
	else
	  totalWeight = totalWeight + log10(1+pow(10,curWeight-totalWeight));
      }
      if (fh == IH_MIN_MAX_CLIQUE || fh == IH_MIN_STATE_SPACE) {
	for (unsigned i=0;i<Pcliques.size();i++) {
	  float curWeight = computeWeight(Pcliques[i].nodes);
	  if (curWeight > maxWeight) maxWeight = curWeight;
	  totalWeight = totalWeight + log10(1+pow(10,curWeight-totalWeight));
	}
	for (unsigned i=0;i<Ecliques.size();i++) {
	  float curWeight = computeWeight(Ecliques[i].nodes);
	  if (curWeight > maxWeight) maxWeight = curWeight;
	  totalWeight = totalWeight + log10(1+pow(10,curWeight-totalWeight));
	}
      }
      if (fh == IH_MIN_MAX_CLIQUE || fh == IH_MIN_MAX_C_CLIQUE) {
	score.push_back(maxWeight);
	infoMsg(Low,"  Interface Score: set has max %sclique weight = %f\n",
		(fh == IH_MIN_MAX_C_CLIQUE ?"C ":""),
		maxWeight);
      } else {
	score.push_back(totalWeight);
	infoMsg(Low,"  Interface Score: set has total %sweight = %f\n",
		(fh == IH_MIN_C_STATE_SPACE?"C ":""),
		totalWeight);
      }

      deleteNodes(Pc);
      deleteNodes(Cc);
      deleteNodes(Ec);      

    } else
      warning("Warning: invalid variable set score given. Ignored\n");
  }

  return;
}



/*-
 *-----------------------------------------------------------------------
 * GMTemplate::cloneWithtoutParents()
 *   Clone a set of random variables from 'in' to 'out'. The cloned
 *   variables have only empty parents, children, and neighbors structures.

 * Preconditions:
 *   'in' is a set of random variables to be cloned.
 *
 * Postconditions:
 *   'out' is a clone of 'in' but without parents,children,neighbors.
 *   
 *
 * Side Effects:
 *     none
 *
 * Results:
 *     returns the in_to_out mapping
 *
 *
 *-----------------------------------------------------------------------
 */
void
GMTemplate::
cloneWithoutParents(const set<RandomVariable*>& in, 
		    set<RandomVariable*>& out,
		    map < RandomVariable*, RandomVariable* >& in_to_out)
{
  in_to_out.clear();
  out.clear();

  for (set<RandomVariable*>::iterator i=in.begin();
       i != in.end(); i++) {

    // sanity check, to ensure a node is not its own neighbor
    assert ( (*i)->neighbors.find((*i)) == (*i)->neighbors.end() );

    RandomVariable*rv = (*i)->cloneWithoutParents();
    out.insert(rv);
    in_to_out[(*i)] = rv;

  }
}



/*-
 *-----------------------------------------------------------------------
 * GMTemplate::setPartitionChildrenNeighbors()
 *   This routine takes as input an original partition S from an unrolled
 *   graph, and its unfinished cloned version Sc (unfinished in that
 *   it has no parents, neighbors (and children) set up yet). It also
 *   takes a mapping from S to Sc, and the corresponding mappings from
 *   the other partitions whatever they are, called O1 and O1. 
 *   
 *   It then sets the neighbors structures for all of Sc from S. The neighbors
 *   are variables that are forced to be part of the partition S itself.
 *
 *   It then sets the parents variables for each variable in Sc. The parents
 *   (and also the children) might NOT be fully contained in Sc, so it needs
 *   to use the mappings for O1 and O2 to get the location of the correspondly
 *   cloned variables for those partitions.
 *
 * Preconditions:
 *   Sc is unfishined, i.e., variables have been cloned without parents
 *
 * Postconditions:
 *   Sc is finished, as described above.
 *   
 *
 * Side Effects:
 *   Changes parents, neighbors, and children of all variables in Sc
 *
 * Results:
 *     results returned via output variable Sc and its parents, neighbors, children
 *
 *
 *-----------------------------------------------------------------------
 */
void
GMTemplate::
setPartitionParentsChildrenNeighbors(const set<RandomVariable*>& S,
				     set<RandomVariable*>& Sc,
				     // next 3 should be const but ther eis no "op[] const"
				     map < RandomVariable*, RandomVariable* >& S_in_to_out,
				     map < RandomVariable*, RandomVariable* >& O1_in_to_out,
				     map < RandomVariable*, RandomVariable* >& O2_in_to_out)
{

  for (set<RandomVariable*>::iterator i=S.begin();i != S.end(); i++) {

    RandomVariable*rv = (*i);

    // first set up new neighbors for S_in_to_out[rv]
    // Note: Neighbors are defined to point ONLY TO OTHER VARIABLES
    // IN THE SET S_in_to_out (i.e., the members of the partition). Any
    // other parents or children are not contained in that set are not included.
    set<RandomVariable*> tmp;
    for (set<RandomVariable*>::iterator j = rv->neighbors.begin();
	 j != rv->neighbors.end(); j++) {    
      if (S_in_to_out.find((*j)) != S_in_to_out.end()) {
	// then it is included in this set.
	tmp.insert(S_in_to_out[(*j)]);
      }
    }
    S_in_to_out[rv]->neighbors = tmp;
    // assertion to make sure that no node has itself as neighbor.
    assert( S_in_to_out[rv]->neighbors.find(S_in_to_out[rv])
	    == S_in_to_out[rv]->neighbors.end() );

    // next, set new sparents for in_to_out[rv].
    // Note that the parents might be outside of the set S.
    vector<RandomVariable *> sParents;
    for (unsigned l=0;l<rv->switchingParents.size();l++) {
      // grab a copy for readability
      RandomVariable* const par = rv->switchingParents[l];
      if (S_in_to_out.find(par) != S_in_to_out.end())
	sParents.push_back(S_in_to_out[par]);
      else if (O1_in_to_out.find(par) != O1_in_to_out.end())
	sParents.push_back(O1_in_to_out[par]);
      else if (O2_in_to_out.find(par) != O2_in_to_out.end())
	sParents.push_back(O2_in_to_out[par]);
      else
	// this shouldn't happen since the parent should live
	// in one of the three partitions.
	assert ( 0 );
    }

    // next, set conditional parents
    vector< vector < RandomVariable* > > cParentsList;
    cParentsList.resize(rv->conditionalParentsList.size());
    for (unsigned l=0;l<rv->conditionalParentsList.size();l++) {
      for (unsigned m=0;m<rv->conditionalParentsList[l].size();m++) {
	// grab a copy for readability
	RandomVariable* const par = rv->conditionalParentsList[l][m];
	if (S_in_to_out.find(par) != S_in_to_out.end())
	  cParentsList[l].push_back(S_in_to_out[par]);
	else if (O1_in_to_out.find(par) != O1_in_to_out.end())
	  cParentsList[l].push_back(O1_in_to_out[par]);
	else if (O2_in_to_out.find(par) != O2_in_to_out.end())
	  cParentsList[l].push_back(O2_in_to_out[par]);
	else
	  // this shouldn't happen since the parent should live
	  // in one of the three partitions.
	  assert ( 0 );
      }
    }
    S_in_to_out[rv]->setParents(sParents,cParentsList);
  }
}


/*-
 *-----------------------------------------------------------------------
 * GMTemplate::setUpClonedPartitionGraph()
 *   Given a P,C, and E that has been cloned without parents, this 
 *   will set up the partitions as follows.
 *   only the members that are in the corresponding 'in' set (i.e., if
 *   any parents, children, neighbors, in 'in' pointed to variables
 *   outside of 'in', then the corresponding variables in 'out' do not
 *   contain those parents,children,neighbors.
 *
 *   There are two versions of this routine, one which returns
 *   the in to out variable map in case that might be useful.
 *
 * Preconditions:
 *   'in' is a set of random variables to be cloned.
 *
 * Postconditions:
 *   'out' is a clone of 'in' but without parents,children,neighbors.
 *   
 *
 * Side Effects:
 *     none
 *
 * Results:
 *     returns the in_to_out mapping
 *
 *
 *-----------------------------------------------------------------------
 */
void
GMTemplate::
setUpClonedPartitionGraph(const set<RandomVariable*>& P,
			  const set<RandomVariable*>& C,
			  const set<RandomVariable*>& E,
			  // cloned variables
			  set<RandomVariable*>& Pc,
			  set<RandomVariable*>& Cc,
			  set<RandomVariable*>& Ec,
			  // next 3 should be const but ther eis no "op[] const"
			  map < RandomVariable*, RandomVariable* >& P_in_to_out,
			  map < RandomVariable*, RandomVariable* >& C_in_to_out,
			  map < RandomVariable*, RandomVariable* >& E_in_to_out)
{

  // just do neighbors for now, don't bother with parents, children,
  // and so on.
  // Set the neighbors of out to be the correctly associated
  // variables, but do not include neighbors that are not in the
  // current set (i.e., dissociate with any other possible portion of
  // the network)

  cloneWithoutParents(P,Pc,P_in_to_out);
  cloneWithoutParents(C,Cc,C_in_to_out);
  cloneWithoutParents(E,Ec,E_in_to_out);

  setPartitionParentsChildrenNeighbors(P,Pc,P_in_to_out,C_in_to_out,E_in_to_out);
  setPartitionParentsChildrenNeighbors(C,Cc,C_in_to_out,P_in_to_out,E_in_to_out);
  setPartitionParentsChildrenNeighbors(E,Ec,E_in_to_out,P_in_to_out,C_in_to_out);

}





/*-
 *-----------------------------------------------------------------------
 * GMTemplate::deleteNodes()
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
      case 'N':
	th_v.push_back(TH_MIN_WEIGHT_NO_D);
	break;
      case 'R':
	th_v.push_back(TH_RANDOM);
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
	fh_v.push_back(IH_MIN_MAX_C_CLIQUE);
	break;
      case 'M':
	fh_v.push_back(IH_MIN_MAX_CLIQUE);
	break;
      case 'A':
	fh_v.push_back(IH_MIN_STATE_SPACE);
	break;
      case 'Q':
	fh_v.push_back(IH_MIN_C_STATE_SPACE);
	break;
      case 'N':
	fh_v.push_back(IH_MIN_WEIGHT_NO_D);
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
 * GMTemplate::findBestInterface()
 *
 *  An exponential-time routine to find best left (or right)
 *  interfaces and corresponding partitions. Since the operations for
 *  finding the best left and right interfaces are symmetric, the case
 *  of if we are searching for the best left interface or right
 *  interface is determined entirely based on the arguments that are
 *  passed into these routines (i.e., just take the mirror image of
 *  the arguments for the right interface).  Note that for simplicity,
 *  the names have been defined in terms of the left interface, but
 *  that is only for simplicity.
 *
 *  The routine is given portions of a twice unrolled graph,
 *  P,C1,C2,C3,E, and finds the best (left) interface within C2
 *  starting at the "standard" or initial left interface between C1
 *  and C2.  Note that the routine only uses the partitions C1,C2, and
 *  C3 as P and E are not needed.
 *
 *
 * For left interface:
 *   Given a twice unrolled graph, P,C1,C2,C3,E, find the best left
 *   interface within C2 starting at the "standard" or initial left
 *   interface between C1 and C2.  Note that the routine only uses
 *   C1,C2, and C3 and P and E are not needed.
 *
 * For right interface:
 *   Given a twice unrolled graph, P,C1,C2,C3,E, find the best right
 *   interface within C2 starting at the "standard" or initial right
 *   interface between C3 and C2.  Note that the routine only uses
 *   C1,C2, and C3 and P and E are not needed.
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
 const vector<InterfaceHeuristic>& fh_v,
 // true if we should use the exponential time optimal interface algorithm
 const bool recurse,
 // --------------------------------------------------------------
 // The next 7 input arguments are used only with the optimal
 // interface algorithm when the IH_MIN_MAX_C_CLIQUE or
 // IH_MIN_MAX_CLIQUE heuristics are used:
 // triangulation heuristic
 const vector<TriangulateHeuristic>& th_v,
 // The network unrolled 1 time
 const set<RandomVariable*>& P_u1,
 const set<RandomVariable*>& C1_u1,
 const set<RandomVariable*>& C2_u1,
 const set<RandomVariable*>& E_u1,
 // Mappings from C2 in the twice unrolled network to C1 and C2
 // in the once unrolled network.
 // (these next 2 should be const, but there is no "op[] const")
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

  // Note that the "partition boundary" is the border/line that cuts the edges
  // connecting nodes C_l and nodes left_C_l.

  interfaceScore(fh_v,C_l,left_C_l,
		 th_v,
		 P_u1,C1_u1,C2_u1,E_u1,
		 C2_u2_to_C1_u1,C2_u2_to_C2_u1,
		 best_score);

  if (message(Tiny)) {
    printf("Size of basic interface C_l = %d\n",C_l.size());
    printf("Score of basic interface C_l =");
    for (unsigned i=0;i<best_score.size();i++)
      printf(" %f ",best_score[i]);
    printf("\n");
    // printf("Size of remainder_C_l = %d\n",left_C_l.size());
    {
      printf("Interface nodes include:");
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
		      C2_1,
		      C3,
		      setset,
		      best_left_C_l,
		      best_C_l,
		      best_score,fh_v,
		      th_v,
		      P_u1,C1_u1,C2_u1,E_u1,
		      C2_u2_to_C1_u1,C2_u2_to_C2_u1);
    if (message(Tiny)) {
      printf("Size of best interface = %d\n",best_C_l.size());
      printf("Score of best interface =");
      for (unsigned i=0;i<best_score.size();i++)
	printf(" %f ",best_score[i]);
      printf("\n");
      // printf("Size of best_remainder_C_l = %d\n",best_left_C_l.size());
      {
	printf("Best interface nodes include:");
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
GMTemplate::
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
  // consider all v in the current C_l as candidates to be moved
  // left. Essentially, we advance the "partition boundary" to the right
  // by taking nodes that are right of the boundary and moving those
  // nodes (one by one) to the left of the boundary. This is done
  // recursively, and memoization is employed to ensure that we don't
  // redundantly explore the same partition boundary.
  // In all cases, 
  //  1) the boundary is a set of edges (and is not explicitly
  //     represented in the code.
  //  2) C_l are the nodes adjacent and to the right of the current boundary 
  //  3) left_C_l are the nodes adjacent and to the left of the current boundary.

  // All nodes in C1 must be to the left of the boundary (ensured by
  // start condition and by the nature of the algoritm).  All nodes in
  // C3 must be to the right of the boundary (ensured by "condition ***" below).

  // If C2 consists of multiple chunks (say C2_1, C2_2, etc.) then we
  // should never have that all of C2_1 lies completely to the left of
  // the boundary, since this would be redundant (i.e., a boundary
  // between C1 and C2_1 would be the same edges shifted in time as
  // the boundary between C2_1 and C2_2). This is ensured by
  // "condition ###" below.

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

    // Condition ***: if v has neighbors in C3 (via set intersection),
    // then continue since if v was moved left, we would end up with
    // an invalid interface (i.e., invalid since there would be a node
    // left of C_l that connects directly to the right of C_l).
    set<RandomVariable*> res;
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

    // only do this check if C2_1 is non-empty. If it is empty,
    // we assume the check is not needed.
    if (C2_1.size() > 0) {
      // Condition ###: next check to make sure that we haven't used up
      // all the nodes in C2_1. I.e., make sure that C2_1 is not a
      // proper subset of (left_C_l U {v}) = next_left_C_l.

      set<RandomVariable*> tmp;      
      set_difference(C2_1.begin(),C2_1.end(),
		     next_left_C_l.begin(),next_left_C_l.end(),
		     inserter(tmp,tmp.end()));
      if (tmp.size() == 0)
	continue;
    }

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
		      C2,C2_1,C3,setset,
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
 *   Makes the interfaces in C1_u1 and C2_u1 complete.
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
 set<RandomVariable*>& Ec,
 set<RandomVariable*>& PCInterface,
 set<RandomVariable*>& CEInterface)
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

  set<RandomVariable*> left_C_l_u1C1;
  set<RandomVariable*> left_C_l_u1C2;
  for (set<RandomVariable*>::iterator i = left_C_l_u2C2.begin();
       i != left_C_l_u2C2.end(); i++) {

    left_C_l_u1C1.insert(C2_u2_to_C1_u1[(*i)]);
    left_C_l_u1C2.insert(C2_u2_to_C2_u1[(*i)]);
  }
  
  // Finally, create the modified sets P, C, and E
  //   where P = modified prologue
  //   where C = modified chunk to repeat
  //   where E = modified epilogue to repeat.
  // which are to be triangulated separately. 
  // For the *left* interface, the modified sets are defined
  // as follows:
  //    P = P' + C1'(left_C_l) + C1'(C_l)
  //    C = C1'\C1'(left_C_l) + C2'(left_C_l) + C2'(C_l)
  //    E = C2'\C2'(left_C_l) + E'
  // and the symmetric definitions apply for the right interface.
  // We use the left interface definitions in this code and
  // assume the caller calles with inverted arguments
  // to get the right interface behavior.

  set<RandomVariable*> P = P_u1;
  set<RandomVariable*> C;
  set<RandomVariable*> E = E_u1;

  // Finish P
  set_union(left_C_l_u1C1.begin(),left_C_l_u1C1.end(),
	    C_l_u1C1.begin(),C_l_u1C1.end(),
	    inserter(P,P.end()));

  // Finish C
  set_union(left_C_l_u1C2.begin(),left_C_l_u1C2.end(),
	    C_l_u1C2.begin(),C_l_u1C2.end(),
	    set_difference(C1_u1.begin(),C1_u1.end(),
			   left_C_l_u1C1.begin(),left_C_l_u1C1.end(),
			   inserter(C,C.end())));
  // Finish E
  set_difference(C2_u1.begin(),C2_u1.end(),
		 left_C_l_u1C2.begin(),left_C_l_u1C2.end(),
		 inserter(E,E.end()));

  // These (C_l_u1C1 and C_l_u1C2) are the interfaces which are forced
  // to be complete (i.e., part of a maxclique). This might take
  // the form of of the code:
  //
  //   makeComplete(C_l_u1C1);
  //   makeComplete(C_l_u1C2);
  //   clone(P,Pc);
  //   clone(C,Cc);
  //   clone(E,Ec);
  //   return;
  //
  // but we want to make them complete only after cloning so that
  // this routine can be called again. Therefore, we do:

  map < RandomVariable*, RandomVariable* > P_in_to_out;
  map < RandomVariable*, RandomVariable* > C_in_to_out;
  map < RandomVariable*, RandomVariable* > E_in_to_out;
  set < RandomVariable* > tmp;

  setUpClonedPartitionGraph(P,C,E,Pc,Cc,Ec,P_in_to_out,C_in_to_out,E_in_to_out);

  // next make the left interface that lives in Pc complete
  tmp.clear();
  for (set<RandomVariable*>::iterator i = C_l_u1C1.begin();
       i != C_l_u1C1.end(); i++) {
    tmp.insert(C_in_to_out[(*i)]);
  }
  makeComplete(tmp);
  PCInterface = tmp;

  // next, make the left interface (part of P) that lives in Cc complete 
  tmp.clear();
  for (set<RandomVariable*>::iterator i = C_l_u1C1.begin();
       i != C_l_u1C1.end(); i++) {
    tmp.insert(C_in_to_out[(*i)]);
  }
  makeComplete(tmp);
  // next, make the left interface (part of C) that lives in Cc complete
  tmp.clear();
  for (set<RandomVariable*>::iterator i = C_l_u1C2.begin();
       i != C_l_u1C2.end(); i++) {
    tmp.insert(C_in_to_out[(*i)]);
  }
  makeComplete(tmp);
  CEInterface = tmp;

  // next, make the left interface that lives in Ec complete
  tmp.clear();
  for (set<RandomVariable*>::iterator i = C_l_u1C2.begin();
       i != C_l_u1C2.end(); i++) {
    tmp.insert(E_in_to_out[(*i)]);
  }
  makeComplete(tmp);

}



////////////////////////////////////////////////////////////////////
//        Test Driver
////////////////////////////////////////////////////////////////////




#ifdef MAIN


#endif
