/*
 * GMTK_GMTemplate.h
 *   Basic GM Template and Basic Triangulation Routines.
 *   This includes code that is common to both triangulation and inference,
 *   so does not contain the more elaborate triangulation methods so that
 *   they are not appart of inference code.
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
 *
 * $Header$
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

#include "debug.h"

// class mention for forward references.
class GraphicalModel;
class BoundaryTriangulate;
class Partition;
class GMTemplate;

class Partition {

  friend class GMTemplate;
  friend class BoundaryTriangulate;

public:

  // variables comprising this partition.
  set<RandomVariable*> nodes;  

  // The cliques themselves, used to store the current triangulation
  // of each of the partitions.
  vector<MaxClique> cliques;
  
  // a string with information about the method used to form the cliques
  string triMethod;

  Partition() {}

  // Clone constructor from another Partition, but that uses a new set
  // of random variables, and adjusts the frame of each new set of
  // random variable with offset
  Partition(Partition& from_part,
	    vector <RandomVariable*>& newRvs,
	    map < RVInfo::rvParent, unsigned >& ppf,
	    const unsigned int frameDelta = 0);

  void clear() { nodes.clear(); cliques.clear(); triMethod.clear(); }

  void writeMaxCliques(oDataStreamFile& os);
  void readMaxCliques(iDataStreamFile& is);
  void triangulatePartitionsByCliqueCompletion();

};

#define GMTEMPLATE_UNINITIALIZED_MS  (~(unsigned)0)

class GMTemplate : public IM 
{

  friend class FileParser;
  friend class GraphicalModel;
  friend class Triangulate;
  friend class Partition;
  friend class BoundaryTriangulate;
  friend class JunctionTree;

  // the file parser for this model.
  FileParser& fp;

  // number of chunks in which to find interface boundary
  const unsigned M;

  // chunk skip, Number of chunks that should exist between boundaries
  const unsigned S;

  // a string with information about the boundary method
  string boundaryMethod;


private:
  // private support functions

  static const string P_partition_name;
  static const string C_partition_name;
  static const string E_partition_name;
  static const string PC_interface_name;
  static const string CE_interface_name;

  ////////////////////////////
  // clear up everything.
  void clear() {
    P.clear(); C.clear(); E.clear();
    PCInterface_in_P.clear(); 
    PCInterface_in_C.clear(); 
    CEInterface_in_C.clear();
    CEInterface_in_E.clear();
    boundaryMethod.clear();
  }


  // TODO: put this next routine in the RV .cc file.
  void makeComplete(set<RandomVariable*> &rvs);

  // support for P,C,and E partitions.
  void setUpClonedPartitionGraph(const set<RandomVariable*>& P,
				 const set<RandomVariable*>& C,
				 const set<RandomVariable*>& E,
				 // cloned variables
				 set<RandomVariable*>& Pc,
				 set<RandomVariable*>& Cc,
				 set<RandomVariable*>& Ec,
				 // next 3 should be const but ther eis no "op[] const"
				 map < RandomVariable*, RandomVariable* >& P_in_to_out,
				 map < RandomVariable*, RandomVariable* >& C_in_to_out,
				 map < RandomVariable*, RandomVariable* >& E_in_to_out);
  void setPartitionParentsChildrenNeighbors(const set<RandomVariable*>& S,
					    set<RandomVariable*>& Sc,
					    // next 3 should be const but ther eis no "op[] const"
					    map < RandomVariable*, RandomVariable* >& S_in_to_out,
					    map < RandomVariable*, RandomVariable* >& O1_in_to_out,
					    map < RandomVariable*, RandomVariable* >& O2_in_to_out);
  void cloneWithoutParents(const set<RandomVariable*>& in, 
			   set<RandomVariable*>& out,
			   map < RandomVariable*, RandomVariable* >& in_to_out);

  void writeMaxCliques(oDataStreamFile& os, const vector<MaxClique>& cliques);
  void readMaxCliques(iDataStreamFile& is, const set<RandomVariable*> nodes, vector<MaxClique>& cliques);

  void triangulatePartitionsByCliqueCompletion(vector<MaxClique>& cliques);

public:

  ////////////////////////////////////////////////////////////////////////////
  // The prologue, chunk, and epilogue partitions. Each partition
  // includes the interface to its neighboring partition (e.g., P
  // intersect C is the PC interface) but the variables in the
  // interface between P and C are not the same actuall variables
  // (they have the same name and frame number, but they are clones of
  // each other, so are different C++ RV objects). The reason for this
  // is that each partition can be triangulated separately without needing
  // to worry about what happens in the other partitions. Also, we are guaranteed
  // that when these are read in, the interfaces in each partition are complete.
  Partition P;
  Partition C;
  Partition E;

  // Interface between P and C, variables in P
  set<RandomVariable*> PCInterface_in_P;
  // Interface between P and C, variables in C
  set<RandomVariable*> PCInterface_in_C;
  // Interface between C and E, variables in C
  set<RandomVariable*> CEInterface_in_C;
  // Interface between C and E, variables in E
  set<RandomVariable*> CEInterface_in_E;

  // public interface

  // the extension of the file
  static const string fileExtension;

  ////////////////////////////////////////////////////////////
  // constructors/destructors
  ////////////////////////////////////////////////////////////
  GMTemplate(FileParser& arg_fp,
	     const unsigned arg_M,
	     const unsigned arg_S)
    : fp(arg_fp),M(arg_M),S(arg_S)
  {
    clear(); 
  }

  GMTemplate(FileParser& arg_fp)
    : fp(arg_fp),
      M(GMTEMPLATE_UNINITIALIZED_MS),
      S(GMTEMPLATE_UNINITIALIZED_MS)
  {
    clear();
  }

  GMTemplate(GMTemplate& t)
    : fp(t.fp)
  {
    error("not yet implemented");
  }
  
  GMTemplate& operator=(const GMTemplate&) 
  {
    error("not yet implemented");
    return *this;
  }

  ~GMTemplate() {}


  // returning M and S by their "proper" names.
  unsigned maxNumChunksInBoundary() { return M; }
  unsigned chunkSkip() { return S; }

  // Read partition information into file
  void createPartitions(const set<RandomVariable*>& P,
			const set<RandomVariable*>& C,
			const set<RandomVariable*>& E,
			const set<RandomVariable*>& PCInterface,
			const set<RandomVariable*>& CEInterface);

  // Write partition information into file
  void writePartitions(oDataStreamFile& os);

  // Read partition information into file
  void readPartitions(iDataStreamFile& is);

  // Write clique information into file
  void writeMaxCliques(oDataStreamFile& os);

  // Read clique information into file, and triangulate
  // the resulting paritions while reading the cliques.
  void readMaxCliques(iDataStreamFile& is);

  // given the cliques, triangulate the partitions.
  void triangulatePartitionsByCliqueCompletion();

  // routine to compute how much to unroll (for basic template P,C, E
  // and modified template P',C',E') and any adjustments to
  // the observation matrix needed (to define where start is).
  enum JustifyType { leftJustify, rightJustify, centerJustify };
  bool computeUnrollParamaters(const unsigned numFrames,
			       unsigned& basicTemplateUnrollAmount,
			       unsigned& modifiedTemplateUnrollAmount,
			       unsigned& numUsableFrames,
			       unsigned& frameStart,
			       const JustifyType justifyType=leftJustify);


};


#endif

