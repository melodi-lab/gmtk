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
class JT_Partition : public Partition {

  friend class JunctionTree;
public:
  
  // The separators for this partition.  If this is a P partition,
  // then all of the separators in this partition are between cliques
  // that live entirely within this partition.  If this is a C or an E
  // partition, then most of the sperators are between cliques that
  // live entirely within this partition.  The last separator in this
  // vector is the seperator between the right interface (RI) clique
  // of the adjacent partition on the left, and the left interface
  // (LI) clique of this partition. The reason for this is that, if
  // this is a P partition, there is no left interface separator, but
  // there is with a C or an E partition.
  vector<SeparatorClique> separators;

  // create an empty one to be filled in later.
  JT_Partition() {}

  // constructor
  JT_Partition(Partition& from_part,
		     vector <RandomVariable*>& newRvs,
		     map < RVInfo::rvParent, unsigned >& ppf,
		     const unsigned int frameDelta = 0);

}; 

// Still another version of a partition object. This one is used only
// in the last stage, during actuall inference and has no STL objects.
// When we do a final partition unroll, we have a list of these
// objects.  It is not inherited from JT_Partition because we do
// not want to store space for STL members.
class JT_InferencePartition {
  friend class JunctionTree;
public:

  // original partition that this has been cloned from.
  JT_Partition& origin;
  sArray< InferenceMaxClique > maxCliques;
  sArray< InferenceSeparatorClique > separatorCliques;

  // WARNING: constructor hack to create a very broken object with
  // non-functional reference objects (in order to create an array of
  // these objects and then initialize them later with appropriate
  // references). Do not use until after proper re-constructor.
  JT_InferencePartition() : origin(*((JT_Partition*)NULL)) {}
  // normal (or re-)constructor
  JT_InferencePartition(JT_Partition& _origin,
			vector <RandomVariable*>& newRvs,
			map < RVInfo::rvParent, unsigned >& ppf,
			const unsigned int frameDelta);
  // destructor
  ~JT_InferencePartition() {}

};


class JunctionTree {

  friend class GMTemplate;
  friend class BoundaryTriangulate;

  // the base partitions from which real unrolled things are cloned from.
  JT_Partition P1; 
  JT_Partition C1; 
  JT_Partition C2; 
  JT_Partition C3; 
  JT_Partition E1; 
  // extra one for unrolling 0 times.
  JT_Partition Cu0; 

  // between partitions, we need extra separator cliques that are
  // between the corresponding partitions interface cliques. These
  // are used to instantiate the separators. Note we need
  // four of them since we might be in case where we unroll 1x.
  // SeparatorClique base_separator_P_C;
  // SeparatorClique base_separator_C_C;
  // SeparatorClique base_separator_C_E;

  // the real partitions, where inference will take place
  // and which will be unrolled depending on the observation vector.
  sArray <JT_InferencePartition> jtIPartitions;

  // after unrolling, these are the separators between the partitions.
  // partitionSeparators.size() = partitions.size() - 1;
  // vector <SeparatorClique> partitionSeparators;


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
  // the overal root node in the JT is on the far right within E
  // (which might not be the best order).
  // NOTE: These are kept here rather than in the partitions,
  // since they are re-used for all cloned partitions.
  vector< pair<unsigned,unsigned> > P1_message_order;
  vector< unsigned > P1_leaf_cliques;
  vector< pair<unsigned,unsigned> > C1_message_order;
  vector< unsigned > C1_leaf_cliques;
  // C2's message order is the same as C1's message order since
  ///   they have same root (RI) clique.
  // C2's leaf cliques the same as C3's leaf cliques since
  //    they have same left interface clique.
  vector< pair<unsigned,unsigned> > C3_message_order;
  vector< unsigned > C3_leaf_cliques;
  vector< pair<unsigned,unsigned> > E1_message_order;  
  vector< unsigned > E1_leaf_cliques;
  // Cu0's message order is the same as C3's message order since
  //    they have same root (RI) clique.
  // Cu0's leaf cliques the same as C1's leaf cliques since
  //    they have same left interface clique.

  // Q: why is C3's message order not the same as C2's?
  // A: because C3 might have a different root than C2.


  // A version of unroll that starts with the gm_template and fills up
  // base partitions.
  void base_unroll();

  // Helper routines that are private (only called by other member
  // functions of this class). 
  void setUpMessagePassingOrderRecurse(JT_Partition& part,
				       const unsigned root,
				       vector< pair<unsigned,unsigned> >&order,
				       const unsigned excludeFromLeafCliques,
				       vector< unsigned>& leaf_cliques);
  void assignRVToClique(const char *const partName,
			JT_Partition&part,
			const unsigned root,
			unsigned depth,
			RandomVariable* rv,
			set<RandomVariable*>& parSet,
			bool& assigned,
			multimap< vector<float>, unsigned >& scoreSet);
  void createDirectedGraphOfCliquesRecurse(JT_Partition& part,
					   const unsigned root,
					   vector< bool >& visited);
  void getCumulativeAssignedNodes(JT_Partition& part,
				  const unsigned root);


  
public:

  // constructor
  JunctionTree(GMTemplate& arg_gm_template)
    : fp(arg_gm_template.fp),
      gm_template(arg_gm_template) {}

  // the fixed file parser for this model, for RV unrolling, etc.
  FileParser& fp;

  // The fixed gm_template for this model, contains the
  // pre-triangulated graph.
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
  void computePartitionInterface(JT_Partition& part1,
				 unsigned int& part1_ric,
				 JT_Partition& part2,
				 unsigned int& part2_lic,
				 bool& icliques_same);


  // root the JT
  void createDirectedGraphOfCliques();
  void createDirectedGraphOfCliques(JT_Partition& part,
				    const unsigned root);



  // Assign probability giving random variables to cliques (i.e.,
  // these are assigned only to cliques such that the random variables
  // and *all* their parents live in the clique, plus some other
  // criterion in order to make message passing as efficient as
  // possible).
  void assignRVsToCliques();
  void assignRVsToCliques(const char *const partName,
			  JT_Partition&part,
			  const unsigned rootClique);


  // For the three partitions, set up the different message passing
  // orders that are to be used. This basically just does a tree
  // traversal using the previously selected root.
  void setUpMessagePassingOrders();
  void setUpMessagePassingOrder(JT_Partition& part,
				const unsigned root,
				vector< pair<unsigned,unsigned> >&order,
				const unsigned excludeFromLeafCliques,
				vector< unsigned>& leaf_cliques);

  // Separator creation, meaning create the seperator objects
  // both within and between partitions. Given two neighboring
  // partitions L and R, the separator between the interface
  // cliques in L and R is contained in R.
  void createSeparators(JT_Partition& part,
			vector< pair<unsigned,unsigned> >&order);
  void createSeparators();


  // Separator iteration order and accumulated set intersection
  // creation for separator driven clique potential creation, and
  // also updates the seperators partial accumulator structure and
  // sets up cliques other variables.
  void computeSeparatorIterationOrder(MaxClique& clique,
				      JT_Partition& part);
  void computeSeparatorIterationOrders(JT_Partition& part);
  void computeSeparatorIterationOrders();

  // 
  // Print all information about the JT. Must
  // have had computeSeparatorIterationOrders() called
  // already.
  void printAllJTInfo(char* fileName);
  void printAllJTInfo(FILE* f,JT_Partition& part,const unsigned root);
  void printAllJTInfoCliques(FILE* f,JT_Partition& part,const unsigned root,const unsigned treeLevel);
  void printMessageOrder(FILE *f,vector< pair<unsigned,unsigned> >& message_order);


  // 
  // Do some last-minute data structure setup to prepare for
  // unrolling to work (such as preliminary and pre work for
  // leaving STL, etc.)
  void prepareForUnrolling(JT_Partition& part);
  void prepareForUnrolling();

  // Set up internal structures for unrolled network k>=0 times, where
  // k is the number of times C' is duplicated (so unroll by 0 means
  // the basic template-prime, unroll by one means two C's, etc.).
  // Unrolling only affects the non-STL data structures.
  // void unroll(unsigned k);

  // unroll for frames = numFrames 
  // Unrolling only affects the non-STL data structures.
  // Returns number of frames actually used, or 0 if invalid num frames.
  unsigned unroll(unsigned numFrames);

  // Perhaps make different unrolls for decoding, unroll for EM
  // training unroll for viterbi training, etc.
  // ...

  // basic collect evidence phase on basic structures.
  void collectEvidence();
  // void distributeEvidence() {}

  // compute P(E), probability of the evidence
  logpr probEvidence();

  // actuall message routines.
  // void collectMessage(MaxClique& from,MaxClique& to);
  // void distributeMessage(MaxClique& from,MaxClique& to);

};


#endif

