/*
 * GMTK_JunctionTree.h
 *   GMTK Junction Tree, for three partitions.
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

#ifndef GMTK_JUNCTIONTREE_H
#define GMTK_JUNCTIONTREE_H

#include <vector>
#include <string>
#include <map>

#include <stdio.h>
#include <stdlib.h>

#include "GMTK_RandomVariable.h"
#include "GMTK_FileParser.h"
#include "GMTK_GMTemplate.h"
#include "GMTK_MaxClique.h"


#include "debug.h"

// class mention for forward references.
class GraphicalModel;
class BoundaryTriangulate;
class Partition;
class GMTemplate;

class JunctionTree {

  friend class GMTemplate;
  friend class BoundaryTriangulate;

  // the base partitions from which all real partitions will be
  // created, where real partitions are the ones in which
  // inference will take place.
  vector <Partition> base_partitions;

  // the real partitions, where inference will take place
  // and which will be unrolled depending on the observation vector.
  vector <Partition> partitions;

  // identities of cliques in junction trees.
  // for P, P's right interface to C (a root)
  unsigned P_ri_to_C; 
  // for C
  //    C's left interface to P
  unsigned C_li_to_P;
  //    C's left interface to C
  unsigned C_li_to_C;
  //    C's right interface to C (a root)
  unsigned C_ri_to_C;
  //    C's right interface to E (a root)
  unsigned C_ri_to_E;
  // for E, E's left interface to C
  unsigned E_li_to_C;
  // root inside of E.
  unsigned E_root_clique;

  // Message passing orders for each partition.  Increasing index
  // order is 'collect evidence' phase from left to right in directio
  // of time, and decreasing order is 'distribute evidence' phase from
  // right to left in direction of time. Note that this assumes that
  // the overal root node in the JT is on the far right within E (which
  // might not be the best order).
  vector< pair<unsigned,unsigned> > P_to_C_message_order;
  vector< pair<unsigned,unsigned> > C_to_C_message_order;
  vector< pair<unsigned,unsigned> > C_to_E_message_order;
  vector< pair<unsigned,unsigned> > E_message_order;  


  // Helper routines that are private (only called by
  // other member functions of this class).

  void setUpMessagePassingOrderRecurse(Partition& part,
				       const unsigned root,
				       vector< pair<unsigned,unsigned> >&order);


  void assignRVToClique(Partition&part,
			const unsigned root,
			unsigned depth,
			RandomVariable* rv,
			set<RandomVariable*>& parSet,
			bool& assigned,
			multimap< vector<float>, unsigned >& scoreSet);


  // a version of unroll that starts with the gm_template and
  // fills up base_partitions.
  void base_unroll(unsigned k=1);

  void createDirectedGraphOfCliquesRecurse(Partition& part,
					   const unsigned root,
					   vector< bool >& visited);

  void getCumulativeAssignedNodes(Partition& part,
				  const unsigned root,
				  set<RandomVariable*> &res);

public:

  JunctionTree(GMTemplate& arg_gm_template)
    : fp(arg_gm_template.fp),
      gm_template(arg_gm_template)
  {
  }

  // the fixed file parser for this model, for unrolling, etc.
  FileParser& fp;
  // the fixed gm_template for this model
  GMTemplate& gm_template;

  struct Edge {
    unsigned clique1;
    unsigned clique2;
    unsigned weight;
  };

  // for sorting edges in decreasing weight order.
  struct EdgeCompare {  
    bool operator() (const Edge& a, 
		     const Edge& b) {
      return (a.weight) > (b.weight);
    }
  };

  // create a junction tree within a partition.
  void createPartitionJunctionTree(Partition& part);

  // routine to find the interface cliques of the partitions
  void computePartitionInterfaces();
  // routine to find the interface cliques of a partition
  void computePartitionInterface(Partition& part1,
				 unsigned int& part1_ric,
				 Partition& part2,
				 unsigned int& part2_lic);

  // create the three junction trees for the basic partitions.
  void createPartitionJunctionTrees() {
    createPartitionJunctionTree(gm_template.P);
    createPartitionJunctionTree(gm_template.C);
    createPartitionJunctionTree(gm_template.E);
  }

  // Set up internal structures for unrolled network k>=0 times,
  // where k is the number of times C' is repeated.
  void unroll(unsigned k);

  // For the three partitions, set up the different message passing orders
  // that are to be used.
  void setUpMessagePassingOrders();
  void setUpMessagePassingOrder(Partition& part,
				const unsigned root,
				vector< pair<unsigned,unsigned> >&order);


  // basic collect evidence phase on basic structures.
  void collectEvidence();
  void distributeEvidence();

  // root the JT
  void createDirectedGraphOfCliques();
  void createDirectedGraphOfCliques(Partition& part,
				    const unsigned root);



  // Assign probability giving random variables to cliques (i.e.,
  // these are assigned only to cliques such that the random variables
  // and *all* their parents live in the clique, plus some other
  // criterion in order to make message passing as efficient as
  // possible).
  void assignRVsToCliques();
  void assignRVsToCliques(Partition&part,
			  const unsigned rootClique);


  // actuall message routines.
  void collectMessage(MaxClique& from,MaxClique& to);
  void distributeMessage(MaxClique& from,MaxClique& to);

};


#endif

