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

  void findInterfaceCliques(const set <RandomVariable*>& iNodes,
			    unsigned& iClique,
			    bool& iCliqueSameAsInterface);
public:


  // Interface nodes on the "left" of this partition. I.e., To find
  // the left interface clique, find a clique that is a superset of
  // these nodes. Empty if there is no such set (e.g., for a P
  // partition)
  set <RandomVariable*> liNodes;

  // Interface nodes on the "right" of this partition. I.e., to
  // compute the root clique of this partition, we find a clique that
  // is a superset of these nodes. Empty if there is no such set
  // (e.g., for an E partition).
  set <RandomVariable*> riNodes;

  // Nodes that are not assigned in this partition. If all nodes are
  // forward-time directed, we are guaranteed that they will be
  // assigned in the left adjacent partition. Similarly, if all nodes
  // are backward-time directed, the nodes are assigned in the right
  // adjacent partition. With a bi-directional graph, the nodes could
  // be assigned in either the left or right adjacent partition.
  set <RandomVariable*> unassignedInPartition;
  
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
	       const unsigned int frameDelta,
	       // the left and right interface variables for
	       // this JT partition Empty if doesn't exist
	       // (say for an P or E partition). These have
	       // their own frame deltas since they might be
	       // different.
	       const set <RandomVariable*>& from_liVars,
	       const unsigned int liFrameDelta,
	       const set <RandomVariable*>& from_riVars,
	       const unsigned int riFrameDelta,
	       // Information todo the mapping.
	       vector <RandomVariable*>& newRvs,
	       map < RVInfo::rvParent, unsigned >& ppf);

  JT_Partition(Partition& from_part,
	       const set <RandomVariable*>& from_liVars,
	       const set <RandomVariable*>& from_riVars);
  

  // returns the left and right interface clique. If not defined,
  // sets the variable to ~0x0.
  void findLInterfaceClique(unsigned& liClique,bool& liCliqueSameAsInterface);
  void findRInterfaceClique(unsigned& riClique,bool& riCliqueSameAsInterface);

  // return the index of the clique with max/min weight.
  unsigned cliqueWithMaxWeight();
  unsigned cliqueWithMinWeight();


  void clearDataMemory() {
    for (unsigned i=0;i<cliques.size();i++)
      cliques[i].clearDataMemory();
    for (unsigned i=0;i<separators.size();i++) 
      separators[i].clearDataMemory();
  }

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
  // ~JT_InferencePartition() {}
  
  // EM updating.
  void emIncrement(const logpr probE, 
		   const bool localCliqueNormalization = false);

};


class JunctionTree  {

  friend class GMTemplate;
  friend class BoundaryTriangulate;

  // The set of base partitions from which real unrolled things are cloned from.
  // When unrolling zero time, we get:
  //   u0: P1 Cu0 E1
  // When unrolling 1 or more times, the method depends on
  // if the template was created using either the left or right interface
  // method.
  // If template created using left interface method, we do:
  //  u0: P1 Cu0 E1
  //  u1: P1 Cu0 Co E1 
  //  u2: P1 Cu0 Co Co E1 
  //  u3: P1 Cu0 Co Co Co E1 
  //  u4: etc.
  // If template created using right interface method, we do:
  //  u0: P1 Cu0 E1
  //  u1: P1 Co Cu0 E1 
  //  u2: P1 Co Co Cu0 E1 
  //  u3: P1 Co Co Co Cu0 E1 
  //  u4: etc.
  
  JT_Partition P1; 
  JT_Partition Cu0;  // C when unrolling 0 times
  JT_Partition Co;   // C "other", depending on if right or left interface method is used.
  JT_Partition E1; 

  // Note, while we need extra separator cliques that are between the
  // corresponding partitions interface cliques, these separators will
  // live in the partition on the right of the separator. They will be
  // pointed to by the left interface clique in the partition on the
  // right.

  // The real partitions, where inference will take place and which
  // will be unrolled depending on the observation vector.
  sArray <JT_InferencePartition> jtIPartitions;

  ////////////////////////////////////////////////////////////////////////
  // Island algorithm support variables.
  // 
  // an array that is used by the Island algorithm for inference.
  struct PartitionInfo {
    // pointer to origin
    JT_Partition* JT;
    // offset to set random variables
    unsigned offset;
    // a pointer to the partition information
    JT_InferencePartition* p;
    // Message order for this partition.
    vector< pair<unsigned,unsigned> >* mo;
    // the name (type) of this partition, i.e., either
    // a P1, Cu0, Co, or an E1.
    char *nm;
    // clique number of right interface 
    unsigned ri;
    // clique number of left interface
    unsigned li;
  };
  // Partition Pointer Array of IPartitionInfo used by recursive Island algorithm.
  sArray <PartitionInfo> partPArray;
  // current set of unrolled random variables
  vector <RandomVariable*> cur_unrolled_rvs;
  // current mapping from 'name+frame' to integer index into unrolled_rvs.
  map < RVInfo::rvParent, unsigned > cur_ppf;
  // the evidence probability used during island algorithm.
  logpr cur_prob_evidence;
  ////////////////////////////////////////////////////////////////////////

  // Identities of cliques in junction trees: 
  // for P, 
  //    P's right  interface to C (a root in a JT partition)
  unsigned P_ri_to_C; 
  // for C
  //    C's left interface to P
  unsigned C_li_to_P;
  //    C's left interface to C
  unsigned C_li_to_C;
  //    C's right interface to C (a root in a JT partition)
  unsigned C_ri_to_C;
  //    C's right interface to E (a root in a JT partition)
  unsigned C_ri_to_E;
  // for E, E's left interface to C
  unsigned E_li_to_C;
  // root inside of E.
  unsigned E_root_clique;

  // Booleans telling if the interface cliques of the two partitions
  // are the same, meaning we don't need both and can drop one (to
  // save a bit of computation). These are currently computed but are
  // not yet used for anything.
  bool P_to_C_icliques_same;
  bool C_to_C_icliques_same;
  bool C_to_E_icliques_same;

  // Message passing orders for each partition.  Increasing index
  // order is 'collect evidence' phase from left to right in direction
  // of time, and decreasing order is 'distribute evidence' phase from
  // right to left in direction of time. Note that this assumes that
  // the overal root node in the JT is on the far right within E
  // (which might not be the best order).  NOTE: These are kept here
  // rather than in the partitions, since they are re-used for all
  // cloned partitions.
  vector< pair<unsigned,unsigned> > P1_message_order;
  vector< unsigned > P1_leaf_cliques;
  vector< pair<unsigned,unsigned> > Cu0_message_order;
  vector< unsigned > Cu0_leaf_cliques;
  vector< pair<unsigned,unsigned> > Co_message_order;
  vector< unsigned > Co_leaf_cliques;
  vector< pair<unsigned,unsigned> > E1_message_order;  
  vector< unsigned > E1_leaf_cliques;


  // A version of unroll that starts with the gm_template and fills up
  // base partitions.
  void base_unroll();

  // Helper routines that are private (only called by other member
  // functions of this class).
  static void setUpMessagePassingOrderRecurse(JT_Partition& part,
					      const unsigned root,
					      vector< pair<unsigned,unsigned> >&order,
					      const unsigned excludeFromLeafCliques,
					      vector< unsigned>& leaf_cliques);
  static void assignRVToClique(const char *const partName,
			       JT_Partition&part,
			       const unsigned root,
			       const unsigned depth,
			       RandomVariable* rv,
			       bool& alreadyAProbContributer,
			       set<RandomVariable*>& parSet,
			       multimap< vector<double>, unsigned >& scoreSet);
  static void createDirectedGraphOfCliquesRecurse(JT_Partition& part,
					   const unsigned root,
					   vector< bool >& visited);
  static void getCumulativeAssignedNodes(JT_Partition& part,
					 const unsigned root);
  static void getCumulativeUnassignedIteratedNodes(JT_Partition& part,const unsigned root);

  void ceGatherIntoRoot(JT_InferencePartition& part,
			 const unsigned root,
			 vector< pair<unsigned,unsigned> >& message_order,
			 const char*const part_type_name,
			 const unsigned part_num);

  void ceSendToNextPartition(JT_InferencePartition& previous_part,
			     const unsigned previous_part_root,
			     const char*const previous_part_type_name,
			     const unsigned previous_part_num,
			     JT_InferencePartition& next_part,
			     const unsigned next_part_leaf,
			     const char*const next_part_type_name,
			     const unsigned next_part_num);

  void deScatterOutofRoot(JT_InferencePartition& part,
			  const unsigned root,
			  vector< pair<unsigned,unsigned> >& message_order,
			  const char*const part_type_name,
			  const unsigned part_num);

  void deReceiveToPreviousPartition(JT_InferencePartition& next_part,
				    const unsigned next_part_leaf,
				    const char*const next_part_type_name,
				    const unsigned next_part_num,
				    JT_InferencePartition& previous_part,
				    const unsigned previous_part_root,
				    const char*const previous_part_type_name,
				    const unsigned previous_part_num);

  // Support routines for island algorithm inference.
  void ceGatherIntoRoot(const unsigned part);
  void createPartition(const unsigned part);
  void ceSendToNextPartition(const unsigned part,const unsigned nextPart);
  void deReceiveToPreviousPartition(const unsigned part,const unsigned prevPart);
  void deletePartition(const unsigned part);
  void deScatterOutofRoot(const unsigned part);
  logpr probEvidenceRoot(const unsigned part);
  logpr setRootToMaxCliqueValue(const unsigned part);
  void emIncrementIsland(const unsigned part,
			  const logpr probE, 
			  const bool localCliqueNormalization);
  void collectDistributeIslandRecurse(const unsigned start,
				      const unsigned end,
				      const unsigned base,
				      const unsigned linear_section_threshold,
				      const bool runEMalgorithm,
				      const bool runViterbiAlgorithm,
				      const bool localCliqueNormalization);
  void collectDistributeIslandBase(const unsigned start,
				   const unsigned end,
				   const bool runEMalgorithm,
				   const bool runViterbiAlgorithm,
				   const bool localCliqueNormalization);

public:

  // Set to true if the JT weight that we compute should be an upper
  // bound.  It is not guaranteed to be a tight upper bound, but is
  // guaranteed to at least be *an* upper bound.
  static bool jtWeightUpperBound;

  // When doing separator driven clique instantiation, if this
  // variable is true, we intersect the separators first before we
  // create the clique. If set to false, we take the product of the
  // cliques irrespective of clique intersection.
  static bool separatorIntersection;


  // When doing scoring (prob(evidence)), do we make compute the 'viterbi'
  // score (meaning the score of the most probable set of variables
  // or P(evidence,best_hidden)), or the full inference score,
  // namely \sum_hidden P(evidence,hidden)
  static bool viterbiScore;


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
    vector <double> weights;
  };

  // for sorting edges in decreasing weight order.
  struct EdgeCompare {  
    bool operator() (const Edge& a, 
		     const Edge& b) {
      return (a.weights) > (b.weights);
    }
  };


  // Call many of the routines below in the right order.
  void setUpDataStructures();

  // create the three junction trees for the basic partitions.
  void createPartitionJunctionTrees() {
    createPartitionJunctionTree(gm_template.P);
    createPartitionJunctionTree(gm_template.C);
    createPartitionJunctionTree(gm_template.E);
  }
  // create a junction tree within a partition.
  static void createPartitionJunctionTree(Partition& part);

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
  static void createDirectedGraphOfCliques(JT_Partition& part,
				    const unsigned root);



  // Assign probability giving random variables to cliques (i.e.,
  // these are assigned only to cliques such that the random variables
  // and *all* their parents live in the clique, plus some other
  // criterion in order to make message passing as efficient as
  // possible).
  void assignRVsToCliques();
  static void assignRVsToCliques(const char *const partName,
				 JT_Partition&part,
				 const unsigned rootClique);

  // determine and set the unassignedNodes in each clique
  // in each partition
  void computeUnassignedCliqueNodes();
  static void computeUnassignedCliqueNodes(JT_Partition& part);


  // For the three partitions, set up the different message passing
  // orders that are to be used. This basically just does a tree
  // traversal using the previously selected root.
  void setUpMessagePassingOrders();
  static void setUpMessagePassingOrder(JT_Partition& part,
				       const unsigned root,
				       vector< pair<unsigned,unsigned> >&order,
				       const unsigned excludeFromLeafCliques,
				       vector< unsigned>& leaf_cliques);

  // Separator creation, meaning create the seperator objects
  // both within and between partitions. Given two neighboring
  // partitions L and R, the separator between the interface
  // cliques in L and R is contained in R.
  static void createSeparators(JT_Partition& part,
			       vector< pair<unsigned,unsigned> >&order);
  void createSeparators();


  // Separator iteration order and accumulated set intersection
  // creation for separator driven clique potential creation, and
  // also updates the seperators partial accumulator structure and
  // sets up cliques other variables.
  static void computeSeparatorIterationOrder(MaxClique& clique,
					     JT_Partition& part);
  static void computeSeparatorIterationOrders(JT_Partition& part);
  void computeSeparatorIterationOrders();

  // Computes the preceding iterated unassigned nodes and therein the
  // set of assigned nodes in each clique that should/shouldn't be
  // iterated.
  void getCumulativeUnassignedIteratedNodes();


  // return an upper bound on the weight of the junction tree in the
  // given partition, where the JT weight is defined as the cost of
  // doing collect evidence on this JT.
  static double junctionTreeWeight(JT_Partition& part,
				   const unsigned rootClique);


  // Given a set of maxcliques for a partition, and an interface for
  // this (can be left right, or any set including empty, the only
  // condition is that it must be covered by at least one of the
  // cliques), compute the junction tree for this set and return the
  // estimated JT cost. This is a static routine so can be called from
  // anywhere.
  static double junctionTreeWeight(vector<MaxClique>& cliques,
				   const set<RandomVariable*>& interfaceNodes);
				   
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
  static void prepareForUnrolling(JT_Partition& part);
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

  // Set all random variables to their observed values either from an
  // observation matrix or to the values given in the file. unroll()
  // must be called first!!
  void setObservedRVs(vector <RandomVariable*>& rvs);


  // Perhaps make different unrolls for decoding, unroll for EM
  // training unroll for viterbi training, etc.
  // ...

  // basic collect evidence phase on basic structures.
  void collectEvidence();
  void distributeEvidence();
  // compute P(E), probability of the evidence, after collect evidence has been run.
  logpr probEvidence();
  // version that does unrolling, and const. memory.
  logpr probEvidence(const unsigned numFrames, unsigned& numUsableFrames);
  // version that does unrolling, and const. memory, & stops after timer interupt occurs.
  static bool probEvidenceTimeExpired;
  logpr probEvidenceTime(const unsigned numFrames, unsigned& numUsableFrames, unsigned &numPartitionsDone);

  // return the island's idea of the current prob of evidence
  logpr curProbEvidenceIsland() { return cur_prob_evidence; }


  // EM training increment, for use with collectEvidence and distributeEvidence.
  void emIncrement(const logpr probEvidence,
		   const bool localCliqueNormalization = false);

  // print P(E) to stdout using all cliques. After a ce,de stage,
  // all values should be the same.
  void printAllCliquesProbEvidence();

  // Set the root to the max clique value and return that value.
  logpr setRootToMaxCliqueValue();

  // Routine that calls collect/distribute evidence using the island
  // algorithm (i.e., log-space inference).
  void
  collectDistributeIsland(const unsigned numFrames,
			  unsigned& numUsableFrames,
			  const unsigned base,
			  const unsigned linear_section_threshold,
			  const bool runEMalgorithm = false,
			  const bool runViterbiAlgorithm = false,
			  const bool localCliqueNormalization = false);





  // actuall message routines.
  // void collectMessage(MaxClique& from,MaxClique& to);
  // void distributeMessage(MaxClique& from,MaxClique& to);

  // Returns a good approximation of the weight of (optionally) P, C,
  // and E in the variables pWeight, cWeight, and eWeight for
  // a given number of frame unrollings.
  double junctionTreeWeight(const bool includeP,
			    const bool includeC,
			    const bool includeE,
			    const unsigned numFrames,
			    double& pWeight,
			    double& cWeight,
			    double& eWeight);

  // Compute the 'junction tree weight' (roughly, the log10(cost of
  // doing inference)) for the set of cliques given in cliques. Note,
  // cliques *must* be a valid set of maxcliques of a junction tree --
  // if they are not, unexpected results are returned.
  double junctionTreeWeight(vector<MaxClique>& cliques);


  // used to clear out hash table memory between segments
  void clearDataMemory() {
    P1.clearDataMemory();
    Cu0.clearDataMemory();
    Co.clearDataMemory();
    E1.clearDataMemory();
  }

  // access to the current set of nodes.
  inline vector <RandomVariable*>& curNodes() { return cur_unrolled_rvs; }


};


#endif

