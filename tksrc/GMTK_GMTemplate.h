/*
 * GMTK_GMTemplate.h
 *   Basic GM Template and Basic Triangulation Routines
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2003, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about the suitability of this software
 * for any purpose. It is provided "as is" without express or implied warranty.
 *
 */

#ifndef GMTK_GMTEMPLATE_H
#define GMTK_GMTEMPLATE_H

#include <vector>
#include <string>
#include <map>

#include <stdio.h>
#include <stdlib.h>

#include "GMTK_RandomVariable.h"
#include "GMTK_FileParser.h"
#include "GMTK_MaxClique.h"

class GraphicalModel;
class AnyTimeTriangulation;

struct GMInfo {
  // the modified prologue, chunk, and epilogue
  set<RandomVariable*> P;
  set<RandomVariable*> C;
  set<RandomVariable*> E;
  // Interface between P and C. The set could contain
  // variables within either P or C depending
  // on if left or right interface is created.
  set<RandomVariable*> PCInterface;
  // Interface between C and E. The set could contain
  // variables within either C or E depending
  // on if the left or right interface is created.
  set<RandomVariable*> CEInterface;
  // cliques of prologue, chunk, and epilogue respectively
  vector<MaxClique> Pcliques;
  vector<MaxClique> Ccliques;
  vector<MaxClique> Ecliques;
  // elimination order of prologue, chunk, and epilogue 
  // respectively, but *after* interfaces have already
  // been made complete (i.e., the elimination order
  // won't necessarily produce a triangulated graph
  // unless the interfaces are made complete).
  vector<RandomVariable*> Pordered;
  vector<RandomVariable*> Cordered;
  vector<RandomVariable*> Eordered;
};


class GMTemplate
{
  friend class FileParser;
  friend class GraphicalModel;
  friend class AnyTimeTriangulation;

  // the file parser for this model.
  FileParser& fp;

  // Keep a number of member variables here for convenience.
  // The number of frames in this template
  const unsigned numFrames;
  // Number of frames of prologue=P, chunk=C, epilogue=E
  const unsigned prologueNumFrames;
  const unsigned chunkNumFrames;
  const unsigned epilogueNumFrames;
  // First chunk frame
  const unsigned firstChunkFrame;
  // Last Chunk Frame
  const unsigned lastChunkFrame;

  
public:


  ////////////////////////////////////////////////////////////
  // constructors/destructors
  ////////////////////////////////////////////////////////////
  GMTemplate(FileParser& _fp) 
    : fp(_fp), 
      numFrames(_fp.numFrames()),
      prologueNumFrames(_fp.firstChunkFrame()),
      chunkNumFrames(_fp.lastChunkFrame() - _fp.firstChunkFrame() + 1),
      epilogueNumFrames(_fp.numFrames() - _fp.firstChunkFrame() - 1),
      firstChunkFrame(_fp.firstChunkFrame()),
      lastChunkFrame(_fp.lastChunkFrame()) { }
  ~GMTemplate() {}

  ////////////////////////////////////////////////////////////
  // data types and other stuff
  ////////////////////////////////////////////////////////////

  enum TriangulateHeuristic { /* S */ TH_MIN_SIZE = 1,        
			      /* T */ TH_MIN_TIMEFRAME = 2,
			      /* F */ TH_MIN_FILLIN = 3,
			      /* W */ TH_MIN_WEIGHT = 4,
			      // use average entropy in CPTs
			      /* E */ TH_MIN_ENTROPY = 5,
			      // use position in .str file
			      /* P */ TH_MIN_POSITION_IN_FILE = 6,
			      // use triangulation hint of variables
			      // in .str file given by user.
			      /* H */ TH_MIN_HINT = 7
  };

  enum InterfaceHeuristic { /* S */ IH_MIN_SIZE = 1,        
			    /* F */ IH_MIN_FILLIN = 2,
			    /* W */ IH_MIN_WEIGHT = 3,
			    // use average entropy in CPTs
			    /* E */ IH_MIN_ENTROPY = 4,
			    /* C */ IH_MIN_MAX_C_CLIQUE = 5,
			    /* M */ IH_MIN_MAX_CLIQUE = 6
  };


  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  // The Main Routines To Form optimal P,C,E partitions
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////


  // Given the template, compute the best partitions
  // using the heuristics that are provided. If 'findBestFace'
  // is true, this will run the exponential algorithm (which
  // isn't necessarily bad since for many structures
  // it is still efficient, but for some it might blow
  // up in amount of time to complete, so caller should
  // be forwarned)
  void findPartitions(const string& faceHeuristic,
		      const string& forceLeftRight,
		      const string& triHeuristic,
		      const bool findBestFace,
		      set<RandomVariable*>& P,
		      set<RandomVariable*>& C,
		      set<RandomVariable*>& E,
		      set<RandomVariable*>& PCI,
		      set<RandomVariable*>& CEI);
  void findPartitions(const string& faceHeuristic,
		      const string& forceLeftRight,
		      const string& triHeuristic,
		      const bool findBestFace,
		      GMInfo& info) {
    findPartitions(faceHeuristic,forceLeftRight,triHeuristic,
		   findBestFace,
		   info.P,info.C,info.E,
		   info.PCInterface,info.CEInterface);
  }



  // Find partitions using the information that
  // has been pre-computed and stored in file 'is'
  void findPartitions(iDataStreamFile& is,
		      set<RandomVariable*>& P,
		      set<RandomVariable*>& C,
		      set<RandomVariable*>& E,
		      set<RandomVariable*>& PCI,
		      set<RandomVariable*>& CEI);
  void findPartitions(iDataStreamFile& is,GMInfo& info) {
    findPartitions(is,
		   info.P,info.C,info.E,
		   info.PCInterface,info.CEInterface);
  }


  // Store partition information into file
  void storePartitions(oDataStreamFile& os,const GMInfo& info);
		       
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  // The Main Triangulation Routines
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////

  // Given a set of partitions, triangulate them using
  // the heuristics provided, returning
  // the resulting cliques and elimination order.
  void triangulatePartitions(const string& th,
			     set<RandomVariable*>& P,
			     set<RandomVariable*>& C,
			     set<RandomVariable*>& E,
			     vector<MaxClique>& Pcliques,
			     vector<MaxClique>& Ccliques,
			     vector<MaxClique>& Ecliques,
			     vector<RandomVariable*>& Pordered,
			     vector<RandomVariable*>& Cordered,
			     vector<RandomVariable*>& Eordered);
  void triangulatePartitions(const string& th,GMInfo& info) {
    triangulatePartitions(th,
			  info.P,info.C,info.E,
			  info.Pcliques,info.Ccliques,info.Ecliques,
			  info.Pordered,info.Cordered,info.Eordered);
  }
  void triangulatePartitions(const vector<TriangulateHeuristic>& th_v,
			     set<RandomVariable*>& P,
			     set<RandomVariable*>& C,
			     set<RandomVariable*>& E,
			     vector<MaxClique>& Pcliques,
			     vector<MaxClique>& Ccliques,
			     vector<MaxClique>& Ecliques,
			     vector<RandomVariable*>& Pordered,
			     vector<RandomVariable*>& Cordered,
			     vector<RandomVariable*>& Eordered);
  void triangulatePartitions(const vector<TriangulateHeuristic>& th_v,
			     GMInfo& info) {
    triangulatePartitions(th_v,
			  info.P,info.C,info.E,
			  info.Pcliques,info.Ccliques,info.Ecliques,
			  info.Pordered,info.Cordered,info.Eordered);
  }

  // Given a set of partitions, triangulate them using
  // the elimination order that is given in file 'is' at
  // current file position.
  void triangulatePartitions(iDataStreamFile& is,
			     set<RandomVariable*>& P,
			     set<RandomVariable*>& C,
			     set<RandomVariable*>& E,
			     vector<MaxClique>& Pcliques,
			     vector<MaxClique>& Ccliques,
			     vector<MaxClique>& Ecliques,
			     vector<RandomVariable*>& Pordered,
			     vector<RandomVariable*>& Cordered,
			     vector<RandomVariable*>& Eordered);
  void triangulatePartitions(iDataStreamFile& is,GMInfo& info) {
    triangulatePartitions(is,
			  info.P,info.C,info.E,
			  info.Pcliques,info.Ccliques,info.Ecliques,
			  info.Pordered,info.Cordered,info.Eordered);
  }


  // Given a set of partitions and a triangulation (elimination
  // ordering) given by the arguments, store that
  // information in file pointed to by 'os'
  void storePartitionTriangulation(oDataStreamFile& os,
				   const set<RandomVariable*>& P,
				   const set<RandomVariable*>& C,
				   const set<RandomVariable*>& E,
				   const vector<MaxClique>& Pcliques,
				   const vector<MaxClique>& Ccliques,
				   const vector<MaxClique>& Ecliques,
				   const vector<RandomVariable*>& Pordered,
				   const vector<RandomVariable*>& Cordered,
				   const vector<RandomVariable*>& Eordered);
  void storePartitionTriangulation(oDataStreamFile& os,GMInfo& info) {
    storePartitionTriangulation(os,
				info.P,info.C,info.E,
				info.Pcliques,info.Ccliques,info.Ecliques,
				info.Pordered,info.Cordered,info.Eordered);
  }

  // a version that doesn't write out clique information.
  void storePartitionTriangulation(oDataStreamFile& os,
				   const set<RandomVariable*>& P,
				   const set<RandomVariable*>& C,
				   const set<RandomVariable*>& E,
				   const vector<RandomVariable*>& Pordered,
				   const vector<RandomVariable*>& Cordered,
				   const vector<RandomVariable*>& Eordered);


  // Given the template, just unroll it flat-out a given number of
  // times and triangulate the result (possibly unconstrained
  // but depending on the heuristics given in 'th').
  void unrollAndTriangulate(const string& th,
			    const unsigned numTimes);


  // The main basic triangulation heuristic routine. Given a set of
  // nodes (that have valid 'neighbors' members, triangulate it using
  // the heuristic(s) given and return the resuting set of cliques in
  // no particular order (i.e., not in RIP order).  Note that more
  // than one heuristic can be used a a time in order to break
  // ties. As a last resort (if all the heuristics agree on the same
  // next set of nodes to eliminate) then a uniformly at random choice
  // is made. This routine is guaranteed to run fast an so can be used
  // in online mode. Also, this routine can be called multiple times
  // where the caller chooses the results with the best cliques --
  // this is because each call might produce a different clique set
  // via the internal randomness that can occur if a tie occurs.
  void basicTriangulate(const set<RandomVariable*> nodes,
			const vector<TriangulateHeuristic>& th_v,
			vector<RandomVariable*>& orderedNodes,
			vector<MaxClique>& cliques,
			const bool findCliques = true);


  // Basic triangulation, via elimination from a pre-existing order.
  // Given a set of nodes (that have valid 'neighbors' members,
  // triangulate it using the information given in the file, and
  // return the resuting set of cliques in no particular order (i.e.,
  // not in RIP order)
  void basicTriangulate(iDataStreamFile& is,
			const set<RandomVariable*> nodes,
			vector<RandomVariable*>& orderedNodes,
			vector<MaxClique>& cliques);


  //private:
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  // General Internal Support Routines
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////

  // creates an undirected graph by dropping edge directions.
  bool dropEdgeDirections(vector<RandomVariable*> &rvs);

  // moralize the variables in rvs.
  bool moralize(vector<RandomVariable*> &rvs);

  // compute the weight (log10(state space)) of a set
  // of variables if they were to be placed within one clique.
  float computeWeight(const set<RandomVariable*>& nodes,
		      const RandomVariable* node = NULL);

  // computes the fill in of a set of variables.
  int computeFillIn(const set<RandomVariable*>& nodes);

  // give the score for a set of variables using given heuristic
  void variableSetScore(const vector<InterfaceHeuristic>& th_v,
			const set<RandomVariable*>& varSet,
			vector<float>& score);

  // general routine to compute the score of a candidate interface.
  void interfaceScore(
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
	 vector<float>& score);

  // clone a set of variables
  void clone(const set<RandomVariable*>& in, set<RandomVariable*>& out); 
  void clone(const set<RandomVariable*>& in, set<RandomVariable*>& out,
	     map < RandomVariable*, RandomVariable* >& in_to_out);

  // delete a set of variables
  void deleteNodes(const set<RandomVariable*>& nodes);


  // given a string, create a vector of triangulation heuristics
  void createVectorTriHeuristic(const string& th,
				vector<TriangulateHeuristic>& th_v);
  void createVectorInterfaceHeuristic(const string& th,
				      vector<InterfaceHeuristic>& th_v);


  // make complete the set of random variables given.
  void makeComplete(set<RandomVariable*> &rvs);


  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  // Support Routines for P,C,E Interface 
  ////////////////////////////////////////////////////////////


  // exponential time routines to find best left/right interfaces
  // and corresponding partitions. Since the operations
  // for finding the best left and right interfaces are symmetric,
  // the determiniation of if we are searching for the best
  // left interface or right interface is determined entirely
  // based on the arguments that are passed into these routines. 
  // Note that for simplicity, the names have been defined
  // in terms of the left interface, but that is only for ease of
  // understanding.
  void findBestInterface(
	     const set<RandomVariable*> &C1,
	     const set<RandomVariable*> &C2,
	     const set<RandomVariable*> &C3,
	     set<RandomVariable*> &left_C_l,
	     set<RandomVariable*> &C_l,
	     const vector<InterfaceHeuristic>& fh_v,
	     const bool recurse,
	     const vector<TriangulateHeuristic>& th_v,
	     const set<RandomVariable*>& P_u1,
	     const set<RandomVariable*>& C1_u1,
	     const set<RandomVariable*>& C2_u1,
	     const set<RandomVariable*>& E_u1,
	     // these next 2 should be const, but there is no "op[] const"
	     map < RandomVariable*, RandomVariable* >& C2_u2_to_C1_u1,
	     map < RandomVariable*, RandomVariable* >& C2_u2_to_C2_u1);

  void findBestInterface(
             const set<RandomVariable*> &left_C_l,
	     const set<RandomVariable*> &C_l,
	     const set<RandomVariable*> &C2,
	     const set<RandomVariable*> &C3,
	     set< set<RandomVariable*> >& setset,
	     set<RandomVariable*> &best_left_C_l,
	     set<RandomVariable*> &best_C_l,
	     vector<float>& best_score,
	     const vector<InterfaceHeuristic>& fh_v,
	     const vector<TriangulateHeuristic>& th_v,
	     const set<RandomVariable*>& P_u1,
	     const set<RandomVariable*>& C1_u1,
	     const set<RandomVariable*>& C2_u1,
	     const set<RandomVariable*>& E_u1,
	     // these next 2 should be const, but there is no "op[] const"
	     map < RandomVariable*, RandomVariable* >& C2_u2_to_C1_u1,
	     map < RandomVariable*, RandomVariable* >& C2_u2_to_C2_u1);
  void findInterfacePartitions(
   // input params
   const set<RandomVariable*>& P_u1,
   const set<RandomVariable*>& C1_u1,
   const set<RandomVariable*>& C2_u1,
   const set<RandomVariable*>& E_u1,
   map < RandomVariable*, RandomVariable* >& C2_u2_to_C1_u1,
   map < RandomVariable*, RandomVariable* >& C2_u2_to_C2_u1,
   const set<RandomVariable*>& left_C_l_u2C2,
   const set<RandomVariable*>& C_l_u2C2,
   // output params
   set<RandomVariable*>& Pc,
   set<RandomVariable*>& Cc,
   set<RandomVariable*>& Ec,
   set<RandomVariable*>& PCI,
   set<RandomVariable*>& CEI
   );




};

#endif

