/*
 * GMTK_JunctionTree.h
 *   GMTK Junction Tree. Exact inference support for GMTK.
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
class JunctionTree;

// TODO: perhaps create a subclass of maxClique at some point, rather than
// adding everything for exact inference to the base class.


// child class of partition that includes support for 
// doing exact inference.
class InferencePartition : public Partition {

  friend class JunctionTree;
public:
  
  // the separators for this partition, not including
  // the ones that go between interfaces.
  vector<SeparatorClique> separators;

  // create an empty one to be filled in later.
  InferencePartition() {}

  InferencePartition(Partition& from_part,
		     vector <RandomVariable*>& newRvs,
		     map < RVInfo::rvParent, unsigned >& ppf,
		     const unsigned int frameDelta = 0);

}; 


class JunctionTree {

  friend class GMTemplate;
  friend class BoundaryTriangulate;

  // the base partitions from which real unrolled things are cloned from.
  InferencePartition P1; 
  InferencePartition C1; 
  InferencePartition C2; 
  InferencePartition C3; 
  InferencePartition E1; 
  // extra one for unrolling 0 times.
  InferencePartition Cu0; 

  // between partitions, we need extra separator cliques that are
  // between the corresponding partitions interface cliques. These
  // are used to instantiate the separators. Note we need
  // four of them since we might be in case where we unroll 1x.
  SeparatorClique base_separator_P_C;
  SeparatorClique base_separator_C_C;
  SeparatorClique base_separator_C_E;

  // the real partitions, where inference will take place
  // and which will be unrolled depending on the observation vector.
  vector <InferencePartition> partitions;

  // after unrolling, these are the separators between the partitions.
  // partitionSeparators.size() = partitions.size() - 1;
  vector <SeparatorClique> partitionSeparators;


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

  // booleans telling if the interface cliques
  // of the two partitions are the same, meaning
  // we don't need both and can drop one (to save
  // a bit of computation).
  bool P_to_C_icliques_same;
  bool C_to_C_icliques_same;
  bool C_to_E_icliques_same;


  // Message passing orders for each partition.  Increasing index
  // order is 'collect evidence' phase from left to right in directio
  // of time, and decreasing order is 'distribute evidence' phase from
  // right to left in direction of time. Note that this assumes that
  // the overal root node in the JT is on the far right within E (which
  // might not be the best order).
  // NOTE: These are kept here rather than in the partitions,
  // since they are re-used for all cloned partitions.
  vector< pair<unsigned,unsigned> > P1_message_order;
  vector< pair<unsigned,unsigned> > C1_message_order;
  // C2' message order is the same as C1's message order
  vector< pair<unsigned,unsigned> > C3_message_order;
  vector< pair<unsigned,unsigned> > E1_message_order;  
  // Cu0's message order is the same as C3's message order


  // Helper routines that are private (only called by
  // other member functions of this class).

  void setUpMessagePassingOrderRecurse(InferencePartition& part,
				       const unsigned root,
				       vector< pair<unsigned,unsigned> >&order);


  void assignRVToClique(const char *const partName,
			InferencePartition&part,
			const unsigned root,
			unsigned depth,
			RandomVariable* rv,
			set<RandomVariable*>& parSet,
			bool& assigned,
			multimap< vector<float>, unsigned >& scoreSet);


  // a version of unroll that starts with the gm_template and
  // fills up base partitions.
  void base_unroll();

  void createDirectedGraphOfCliquesRecurse(InferencePartition& part,
					   const unsigned root,
					   vector< bool >& visited);

  void getCumulativeAssignedNodes(InferencePartition& part,
				  const unsigned root);

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


  // create the three junction trees for the basic partitions.
  void createPartitionJunctionTrees() {
    createPartitionJunctionTree(gm_template.P);
    createPartitionJunctionTree(gm_template.C);
    createPartitionJunctionTree(gm_template.E);
  }
  // create a junction tree within a partition.
  void createPartitionJunctionTree(Partition& part);

  // routine to find the interface cliques of the partitions
  void computePartitionInterfaces();
  // routine to find the interface cliques of a partition
  void computePartitionInterface(InferencePartition& part1,
				 unsigned int& part1_ric,
				 InferencePartition& part2,
				 unsigned int& part2_lic,
				 bool& icliques_same);


  // root the JT
  void createDirectedGraphOfCliques();
  void createDirectedGraphOfCliques(InferencePartition& part,
				    const unsigned root);



  // Assign probability giving random variables to cliques (i.e.,
  // these are assigned only to cliques such that the random variables
  // and *all* their parents live in the clique, plus some other
  // criterion in order to make message passing as efficient as
  // possible).
  void assignRVsToCliques();
  void assignRVsToCliques(const char *const partName,
			  InferencePartition&part,
			  const unsigned rootClique);


  // For the three partitions, set up the different message passing orders
  // that are to be used.
  void setUpMessagePassingOrders();
  void setUpMessagePassingOrder(InferencePartition& part,
				const unsigned root,
				vector< pair<unsigned,unsigned> >&order);

  // separator creation
  void createSeparators(InferencePartition& part,
			vector< pair<unsigned,unsigned> >&order);
  void createSeparators();


  // basic collect evidence phase on basic structures.
  void collectEvidence();
  void distributeEvidence();

  // actuall message routines.
  void collectMessage(MaxClique& from,MaxClique& to);
  void distributeMessage(MaxClique& from,MaxClique& to);

  // Set up internal structures for unrolled network k>=0 times,
  // where k is the number of times C' is repeated.
  void unroll(unsigned k);


};


#endif

