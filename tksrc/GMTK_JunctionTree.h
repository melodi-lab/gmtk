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

#include "bp_range.h"

#include "GMTK_RV.h"
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


  struct PriorityClique {
    unsigned clique;
    vector <double> weights;
    PriorityClique(unsigned c,vector <double> w) : clique(c), weights(w) {}
  };

  struct PriorityCliqueCompare {  
    // sort descending
    bool operator() (const PriorityClique& a, 
		     const PriorityClique& b) {
      return (a.weights) > (b.weights);
    }
  };


  void findInterfaceCliques(const set <RV*>& iNodes,
			    unsigned& iClique,
			    bool& iCliqueSameAsInterface,
			    const string priorityStr);
public:


  // Interface nodes on the "left" of this partition. I.e., To find
  // the left interface clique, find a clique that is a superset of
  // these nodes. Empty if there is no such set (e.g., for a P
  // partition)
  set <RV*> liNodes;

  // Interface nodes on the "right" of this partition. I.e., to
  // compute the root clique of this partition, we find a clique that
  // is a superset of these nodes. Empty if there is no such set
  // (e.g., for an E partition).
  set <RV*> riNodes;

  // Nodes that are not assigned in this partition. If all nodes are
  // forward-time directed, we are guaranteed that they will be
  // assigned in the left adjacent partition. Similarly, if all nodes
  // are backward-time directed, the nodes are assigned in the right
  // adjacent partition. With a bi-directional graph, the nodes could
  // be assigned in either the left or right adjacent partition.
  set <RV*> unassignedInPartition;
  
  // The separators for this partition.  If this is a P partition,
  // then all of the separators in this partition are between cliques
  // that live entirely within this partition.  If this is a C or an E
  // partition, then most of the sperators are between cliques that
  // live entirely within this partition.  The last separator in this
  // vector is guaranteed to be the seperator between the right
  // interface (RI) clique of the adjacent partition on the left, and
  // the left interface (LI) clique of this partition. The reason for
  // this is that, if this is a P partition, there is no left
  // interface separator, but there is with a C or an E partition.
  // Created in: JunctionTree::createSeparators(); This final
  // separator is called the LI separator.
  vector<SeparatorClique> separators;

  // The set of factor cliques (i.e., hard and soft constraints) that
  // live in this partition, corresponding to the 'factor' constructs
  // in the .str file.
  vector<FactorClique> factorCliques;

  void useLISeparator()  { separators[separators.size()-1].skipMe = false; }
  void skipLISeparator() { separators[separators.size()-1].skipMe = true; }

  // number of VE separators among the separators in this partition.
  unsigned numVEseps;

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
	       const set <RV*>& from_liVars,
	       const unsigned int liFrameDelta,
	       const set <RV*>& from_riVars,
	       const unsigned int riFrameDelta,
	       // Information todo the mapping.
	       vector <RV*>& newRvs,
	       map < RVInfo::rvParent, unsigned >& ppf);

  JT_Partition(Partition& from_part,
	       const set <RV*>& from_liVars,
	       const set <RV*>& from_riVars);
  

  // returns the left and right interface clique. If not defined,
  // sets the variable to ~0x0.
  void findLInterfaceClique(unsigned& liClique,bool& liCliqueSameAsInterface,
			    const string priorityStr);
  void findRInterfaceClique(unsigned& riClique,bool& riCliqueSameAsInterface,
			    const string priorityStr);

  // return the index of the clique with max/min weight.
  unsigned cliqueWithMaxWeight();
  unsigned cliqueWithMinWeight();


  void clearCliqueSepValueCache(bool force = false) {
    for (unsigned i=0;i<cliques.size();i++)
      cliques[i].clearCliqueValueCache(force);
    for (unsigned i=0;i<separators.size();i++) 
      separators[i].clearSeparatorValueCache(force);
  }

  void clearCliqueValueCache(bool force=false) {
    for (unsigned i=0;i<cliques.size();i++)
      cliques[i].clearCliqueValueCache(force);
  }

  void clearSeparatorValueCache(bool force=false) {
    for (unsigned i=0;i<separators.size();i++) 
      separators[i].clearSeparatorValueCache(force);
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
  sArray< InferenceFactorClique > factorCliques;

  // WARNING: constructor hack to create a very broken object with
  // non-functional reference objects (in order to create an array of
  // these objects and then initialize them later with appropriate
  // references). Do not use until after proper re-constructor.
  JT_InferencePartition() : origin(*((JT_Partition*)NULL)) {}
  // normal (or re-)constructor
  JT_InferencePartition(JT_Partition& _origin,
			vector <RV*>& newRvs,
			map < RVInfo::rvParent, unsigned >& ppf,
			const unsigned int frameDelta);
  // destructor
  // ~JT_InferencePartition() {}
  
  // EM updating.
  void emIncrement(const logpr probE, 
		   const bool localCliqueNormalization = false,
		   const double emTrainingBeam = -LZERO);

};


class JunctionTree {

  friend class GMTemplate;
  friend class BoundaryTriangulate;

  // The set of base partitions from which real unrolled things are cloned from.
  // When unrolling zero time, we get:
  //   u0: P' E'
  // When unrolling 1 or more times, the method depends on
  // if the template was created using either the left or right interface
  // method.
  // If template created using left interface method, we do:
  //  u0: P' E'
  //  u1: P' C' E' 
  //  u2: P' C' C' E' 
  //  u3: P' C' C' C' E'
  //  u4: etc.
  // Note that for left interface, an E1 contains M copies of 
  // the original chunk C, so u0 is the basic template.
  //
  // If template created using right interface method, we do:
  //  u0: P' E'
  //  u1: P' C' E' 
  //  u2: P' C' C' E' 
  //  u3: P' C' C' C' E'
  //  u4: etc.
  // in the right interface case, P' contains an original P and M
  // copies of C.  The next three variables hold partitions P', C',
  // and E', where, for the standard left interface and simple
  // boundary case, we have that P' = P, C' = C , and E' = [C
  // E]. These partitions use the *same* set of C++ random variables
  // (so the intersection of the random variables in P1 and Co will be
  // the interface).
  JT_Partition P1; 
  JT_Partition Co;   // C "other", depending on if right or left interface method is used.
  JT_Partition E1; 

  // The set of random variables corresponding to the union of the rvs
  // P1, Co, E1, corresponding to the template unrolled M+S-1
  // times. In other words, the number of C repetitions is M + S, so
  // we unroll one less than that, see
  // BoundaryTriangulate::findPartitions() for more information. These
  // are determined in base_unroll().
  vector <RV*> partition_unrolled_rvs; 
  // mapping from name(frame) to integer index into unrolled_rvs.
  map < RVInfo::rvParent, unsigned > partition_ppf;

  // The names of the above three partitions to use for printing,
  // debugging, etc.
  static char* P1_n;
  static char* Co_n;
  static char* E1_n;

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
    int offset;
    // a pointer to the partition information
    JT_InferencePartition* p;
    // Message order for this partition.
    vector< pair<unsigned,unsigned> >* mo;
    // the name (type) of this partition, i.e., either
    // a P1, Co, or an E1.
    char *nm;
    // clique number of right interface 
    unsigned ri;
    // clique number of left interface
    unsigned li;
  };
  // Partition Pointer Array of IPartitionInfo used by recursive Island algorithm.
  sArray <PartitionInfo> partPArray;
  // current set of unrolled random variables
  vector <RV*> cur_unrolled_rvs;
  // current mapping from 'name+frame' to integer index into unrolled_rvs.
  map < RVInfo::rvParent, unsigned > cur_ppf;
  // the evidence probability used during island algorithm.
  logpr cur_prob_evidence;
  // the EM training beam used for island training (TODO:, move this elsewhere, perhaps in clique)
  double curEMTrainingBeam;
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
  vector< pair<unsigned,unsigned> > Co_message_order;
  vector< unsigned > Co_leaf_cliques;
  vector< pair<unsigned,unsigned> > E1_message_order;  
  vector< unsigned > E1_leaf_cliques;

  // do a bit of setup for the upcomming inference round.
  void prepareForNextInferenceRound();

  // A version of unroll that starts with the gm_template and fills up
  // base partitions.
  void base_unroll();
  void insertFactorClique(FactorClique& factorClique,FactorInfo& factor);


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
			       RV* rv,
			       unsigned& numberOfTimesAssigned,
			       set<RV*>& parSet,
			       const bool allParentsObserved,
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
  void cePruneRootCliqueOfPartition(const unsigned part);
  void deReceiveToPreviousPartition(const unsigned part,const unsigned prevPart);
  void deletePartition(const unsigned part);
  void deScatterOutofRoot(const unsigned part);
  logpr probEvidenceRoot(const unsigned part);
  logpr setRootToMaxCliqueValue(const unsigned part);
  void emIncrementIsland(const unsigned part,
			 const logpr probE, 
			 const bool localCliqueNormalization);

  void printAllCliques(JT_InferencePartition& part,
		       const unsigned partNo,
		       const char *const nm,
		       BP_Range* rng,
		       FILE* f,
		       const bool normalize,
		       const bool justPrintEntropy = false);
  void printAllCliques(const unsigned part,FILE* f,
		       const bool normalize,
		       const bool justPrintEntropy = false);

  
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


  // Set to true if the JT should create extra separators for any
  // virtual evidence (VE) that might be usefully exploitable
  // computationally in the clique. I.e., we treat the parents
  // of an immediate observed child, where the child is a deterministic 
  // function of the parents, a a separator over the parents
  // to be intersected as normal with all the other separators.
  static unsigned useVESeparators;
  enum VESeparatorType { VESEP_PC = 0x1, VESEP_PCG = 0x2 }; 
  // booleans to indicate where ve-seps should be used.
  static unsigned veSeparatorWhere;
  enum VESeparatorWhere { VESEP_WHERE_P = 0x1, VESEP_WHERE_C = 0x2, VESEP_WHERE_E = 0x4 }; 

  
  // When doing inference if any kind, this variable determines
  // if we should clear the clique and separator value cache
  // between segments/utterances. It might be beneficial, for
  // example, to retain the value cache around between segments/utterances
  // if for example, there are many such values that are common. If
  // not, on the other hand, setting this to true will cause
  // an increase in memory use on each segment.
  static bool perSegmentClearCliqueValueCache;

  // Set to true if the JT weight that we compute should be an upper
  // bound.  It is not guaranteed to be a tight upper bound, but is
  // guaranteed to at least be *an* upper bound.
  static bool jtWeightUpperBound;

  // Set to true if the JT weight scoring mechansim should be more
  // conservative, meaning it should not underestimate the charge of a
  // node in a clique as much. See code for details.
  static bool jtWeightMoreConservative;

  // The priority string for selecting the next edge when constructing
  // the junction tree. Default is in .cc file, and see .cc file for
  // what options are supported.
  static char* junctionTreeMSTpriorityStr;

  // The priority string for selecting which clique of a partition
  // (from the set of valid ones) should be used as the partition
  // interface clique. See .cc file in routine findInterfaceClique()
  // for documentation.
  static char* interfaceCliquePriorityStr;

  // Set to > 0.0 if the JT weight that we compute should heavily
  // penalize any unassigned iterated nodes. Penalty = factor
  // that gets multiplied by number of unassigned iterated.
  static float jtWeightPenalizeUnassignedIterated;
  
  // scaling factors (must be > 0 and <= 1) corresponding
  // to how much the separator nodes' charge gets scaled. The more pruning
  // we do during inference, the lower this should go.
  static float jtWeightSparseNodeSepScale;
  static float jtWeightDenseNodeSepScale;


  // When doing scoring (prob(evidence)), do we make compute the 'viterbi'
  // score (meaning the score of the most probable set of variables
  // or P(evidence,best_hidden)), or the full inference score,
  // namely \sum_hidden P(evidence,hidden)
  static bool viterbiScore;

  // range of cliques within each partition to print out when doing
  // CE/DE inference. If these are NULL, then we print nothing.
  BP_Range* pPartCliquePrintRange;
  BP_Range* cPartCliquePrintRange;
  BP_Range* ePartCliquePrintRange;

  // constructor
  JunctionTree(GMTemplate& arg_gm_template)
    : curEMTrainingBeam(-LZERO),
      fp(arg_gm_template.fp),
      gm_template(arg_gm_template) 
  {
    pPartCliquePrintRange = cPartCliquePrintRange = ePartCliquePrintRange = NULL;
  }
  ~JunctionTree() {
    delete pPartCliquePrintRange;
    delete cPartCliquePrintRange;
    delete ePartCliquePrintRange;
  }

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
    // sort descending
    bool operator() (const Edge& a, 
		     const Edge& b) {
      return (a.weights) > (b.weights);
    }
  };


  // Call many of the routines below in the right order.
  void setUpDataStructures(const char* varPartitionAssignmentPrior,
			   const char *varCliqueAssignmentPrior);

  // create the three junction trees for the basic partitions.
  void createPartitionJunctionTrees(const string pStr = junctionTreeMSTpriorityStr) {
    createPartitionJunctionTree(gm_template.P,pStr);
    createPartitionJunctionTree(gm_template.C,pStr);
    createPartitionJunctionTree(gm_template.E,pStr);
  }
  // create a junction tree within a partition.
  static void createPartitionJunctionTree(Partition& part, 
					  const string pStr = junctionTreeMSTpriorityStr);

  // routine to find the interface cliques of the partitions
  void computePartitionInterfaces();
  // routine to create the factors in the appropriate partitions
  void createFactorCliques();
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
  void assignRVsToCliques(const char* varPartitionAssignmentPrior,
			  const char *varCliqueAssignmentPrior);
  static void assignRVsToCliques(const char *const partName,
				 JT_Partition&part,
				 const unsigned rootClique,
				 const char* varPartitionAssignmentPrior,
				 const char *varCliqueAssignmentPrior);



  void assignFactorsToCliques();
  void assignFactorsToCliques(JT_Partition& part);


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

  // create the virtual evidence separators
  static void createVESeparators(JT_Partition& part);


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

  // compute the assignment order for nodes in this
  // partition's cliques relative to each clique's incomming separators, and while
  // doing so, also set the dispositions for each of the resulting
  // nodes in each clique.
  void sortCliqueAssignedNodesAndComputeDispositions(const char *varCliqueAssignmentPrior);
  void sortCliqueAssignedNodesAndComputeDispositions(JT_Partition& part,
						     const char *varCliqueAssignmentPrior);


  // return an upper bound on the weight of the junction tree in the
  // given partition, where the JT weight is defined as the cost of
  // doing collect evidence on this JT.
  static double junctionTreeWeight(JT_Partition& part,
				   const unsigned rootClique,
				   set<RV*>* lp_nodes,
				   set<RV*>* rp_nodes);

  // Given a set of maxcliques for a partition, and an interface for
  // this (can be left right, or any set including empty, the only
  // condition is that it must be covered by at least one of the
  // cliques), compute the junction tree for this set and return the
  // estimated JT cost. This is a static routine so can be called from
  // anywhere.
  static double junctionTreeWeight(vector<MaxClique>& cliques,
				   const set<RV*>& interfaceNodes,
				   set<RV*>* lp_nodes,
				   set<RV*>* rp_nodes);
				   
  // 
  // Print all information about the JT. Must
  // have had computeSeparatorIterationOrders() called
  // already.
  void printAllJTInfo(char* fileName);
  void printAllJTInfo(FILE* f,JT_Partition& part,const unsigned root,
		      set <RV*>* lp_nodes,set <RV*>* rp_nodes);
  void printAllJTInfoCliques(FILE* f,JT_Partition& part,const unsigned root,const unsigned treeLevel,
			     set <RV*>* lp_nodes,set <RV*>* rp_nodes);
  void printMessageOrder(FILE *f,vector< pair<unsigned,unsigned> >& message_order);
  void printCurrentRVValues(FILE* f);
  void setCliquePrintRanges(char *p,char*c,char*e);
  void printAllCliques(FILE* f,const bool normalize,const bool justPrintEntropy);

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
  void setObservedRVs(vector <RV*>& rvs);


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
  logpr probEvidenceTime(const unsigned numFrames, unsigned& numUsableFrames, 
			 unsigned &numPartitionsDone, const bool noE = false);

  // return the island's idea of the current prob of evidence
  logpr curProbEvidenceIsland() { return cur_prob_evidence; }
  // set island's EM beam.
  void setCurEMTrainingBeam(const double b) { curEMTrainingBeam = b; }

  // EM training increment, for use with collectEvidence and distributeEvidence.
  void emIncrement(const logpr probEvidence,
		   const bool localCliqueNormalization = false,
		   const double emTrainingBeam = -LZERO);

  // print P(E) to stdout using all cliques. After a ce,de stage,
  // all values should be the same.
  void printProbEvidenceAccordingToAllCliques();

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
  void clearCliqueSepValueCache(bool force=false) {
    P1.clearCliqueSepValueCache(force);
    Co.clearCliqueSepValueCache(force);
    E1.clearCliqueSepValueCache(force);
  }

  // access to the current set of nodes.
  inline vector <RV*>& curNodes() { return cur_unrolled_rvs; }


};


#endif

