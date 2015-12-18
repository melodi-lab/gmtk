/*
 * GMTK_SectionScheduler.h Root class for the inference algorithms at
 *   the time series level.
 *
 * Written by Richard Rogers <rprogers@uw.edu>
 *
 * Copyright (C) 2015 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 * SectionScheduler subclasses handle inference over a series of sections. 
 *
 * The specific inference task subclasses of SectionScheduler are 
 * (X is observed, Q is hidden, T is segment length):
 *
 *  - ProbEvidenceTask       compute P(X_{0:T-1}), forward pass only, O(1) memory
 *  - ViterbiTask            compute argmax_{Q_{0:T-1}} P(Q_{0:T-1} | X_{0:T-1})
 *  - ForwardBackwardTask    compute P(Q_{0:T-1} | X_{0:T-1})
 *  - SmoothingTask          compute argmax_{Q_{t-\tau}} P(Q_{t-\tau} | X_{0:t})  
 *
 * The time-series algorithm subclasses of SectionScheduler are:
 *  - LinearSectionScheduler        Just iterate sequentially over sections
 *  - IslandSectionScheduler        Skip through the time series, leaving "islands" of partial results to save memory
 *  - ArchipelagosSectionScheduler  Parallel version of Island
 * 
 * The SectionScheduler subclasses multiply inherit the *Task classes representing the
 * inference tasks the time series level inference algorithms can perform. Not all 
 * SectionSchedulers support every inference task, e.g., IslandSectionScheduler
 * can't do online filtering/smoothing (OnlineSmoothingTask).
 *
 * The Sections are responsible for doing inference within themselves via whatever
 * SectionInferenceAlgorithm the Section subclass implements. The SectionScheduler
 * classes just pass messages between Sections via instances of InterfaceSeparator.
 *
 */


#ifndef GMTK_SECTIONSCHEDULER_H
#define GMTK_SECTIONSCHEDULER_H

#include <utility>
#include <vector>
#include <set>
#include <map>

#include "logp.h"

// TODO: what's the difference between range and bp_range?
#include "range.h"
#include "bp_range.h"
#include "mArray.h"

#include "fileParser.h"

#include "GMTK_ObservationSource.h"
#include "GMTK_FileParser.h"
#include "GMTK_RV.h"

#include "GMTK_PartitionStructures.h"
#include "GMTK_SectionTablesBase.h"

class SectionIterator;

class SectionScheduler {
  friend class SectionIterator;
  // MUST be friends with all SectionInferenceAlgorithm & SectionTablesBase subclasses
  friend class SectionInferenceAlgorithm;
  friend class SparseJoinInference;
  friend class PedagogicalInference;

  friend class SectionTablesBase;
  friend class SparseJoinSectionTables;
  friend class PedagogicalSectionTables;

 public:


  SectionScheduler(GMTemplate &gm_template, FileParser &fp, ObservationSource *obs_source) :
    p_clique_print_range(NULL), c_clique_print_range(NULL), e_clique_print_range(NULL),
    section_debug_range("all", 0, 0x7FFFFFFF), obs_source(obs_source), fp(fp), gm_template(gm_template)
  {
    assert(obs_source);
  }
    
  virtual ~SectionScheduler() {}


  // Initialize stuff at the model-level. See prepareForSegment() for segment-level initialization.
  // TODO: explain parameters
  virtual void setUpDataStructures(iDataStreamFile &tri_file,
				   char const *varSectionAssignmentPrior,
				   char const *varCliqueAssignmentPrior,
				   bool checkTriFileCards);


  // Prepare to do inference on segment of length T. Note that this
  // is the analog of the current JunctionTree::unroll() which sets up the
  // section_structure_array for at most  P' C' C' E' (4 sections, so
  // O(1) memory), and the section_table_array for between 0 and T
  // sections depending on ZeroTable, ShortTable, or LongTable. 
  // JT::unroll() also allocates O(T) memory (optionally memory mapped)
  // to store the the Viterbi values if they aren't being written to
  // a file.

  // Since we now have the ability to write Viterbi values to files, do
  // we want to drop the ability to store them in memory? This would
  // break things, as a Viterbi file would then be a required argument.
  // It would also be slower for applications that do fit in memory
  // because of the higher overhead of file I/O operations.

  // Maybe rename this to prepareForSegment() or something to avoid confusion
  // with O(T) space graph unrolling?
  // returns number of usable frames
  enum UnrollTableOptions { LongTable, ShortTable, ZeroTable, NoTouchTable };
  virtual unsigned unroll(unsigned numFrames,
			  const UnrollTableOptions tableOption = LongTable,
			  unsigned *totalNumberSections = NULL); 


  // Formerly JunctionTree::printAllJTInfo()
  virtual void printInferencePlanSummary(char const *fileName);

  // Formerly GMTemplate::reportScoreStats()
  virtual void reportScoreStats();


  // NOTE:  This assumes all inference algorithms will have clique-like things they
  //        need to print posteriors of

  // Set the range of selected clique #'s in P', C', E' for printing.
  // TODO: preconditions
  virtual void setCliquePrintRanges(char *p_range, char *c_range, char *e_range);

  // Print to f the order of the variables in each clique selected by setCliquePrintRanges().
  virtual void printCliqueOrders(FILE *f);

  // Returns the size (in # of floats) of the cliques selected by setCliquePrintRanges().
  virtual void getCliquePosteriorSize(unsigned &p_size, unsigned &c_size, unsigned &e_size);



  virtual void printAllCliques(FILE *f,const bool normalize, const bool unlog,
			       const bool justPrintEntropy,
			       ObservationFile *obs_file = NULL);


  // Set the range of section #s that get elevated verbosity
  // P' = 0, C' \in 1, ..., T=2, E' = T-1
  virtual void setSectionDebugRange(Range const &rng) {
    section_debug_range.SetLimits(rng.first(), rng.last()); 
    section_debug_range.SetDefStr(rng.GetDefStr()); 
  }


  // within-section junction tree finding should happen in a SectionInferenceAlgorithm ?

  // The priority string for selecting the next edge when constructing
  // the junction tree. Default is in .cc file, and see .cc file for
  // what options are supported.
  static const char* junctionTreeMSTpriorityStr;


  // generalize to factored interface - multiple "clique" interface

  // The priority string for selecting which clique of a partition
  // (from the set of valid ones) should be used as the partition
  // interface clique. See .cc file in routine findInterfaceClique()
  // for documentation.
  static const char* interfaceCliquePriorityStr;


  // move to SectionInferenceAlgorithm ?
  
  // Set to true if the JT should create extra separators for any
  // virtual evidence (VE) that might be usefully exploitable
  // computationally in the clique. I.e., we treat the parents
  // of an immediate observed child, where the child is a deterministic 
  // function of the parents, a a separator over the parents
  // to be intersected as normal with all the other separators.
  static unsigned useVESeparators;

  enum VESeparatorType { VESEP_PC = 0x1, VESEP_PCG = 0x2 }; 

  // booleans to indicate where ve-seps should be used.
  enum VESeparatorWhere { VESEP_WHERE_P = 0x1, VESEP_WHERE_C = 0x2, VESEP_WHERE_E = 0x4 }; 
  static unsigned veSeparatorWhere;





 
  // When doing inference of any kind, this variable determines
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

#if 0
  // already moved into SectionScheduler

  // The priority string for selecting the next edge when constructing
  // the junction tree. Default is in .cc file, and see .cc file for
  // what options are supported.
  static const char* junctionTreeMSTpriorityStr;

  // The priority string for selecting which clique of a partition
  // (from the set of valid ones) should be used as the partition
  // interface clique. See .cc file in routine findInterfaceClique()
  // for documentation.
  static const char* interfaceCliquePriorityStr;
#endif 

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

  // if true, use mmap() to allocate memory for C_partition_values,
  // otherwise use new
  static bool mmapViterbi;
  
  // if ture, do distribute evidence just within the modified section
  // to compute P(Q_t | X_{0:t})
  static bool sectionDoDist;


  // Viterbi printing triggers
  static char *pVitTrigger;
  static char *cVitTrigger;
  static char *eVitTrigger;

  static bool vitRunLength;

  // For O(1) memory inference, write Viterbi values to this file for
  // later printing by a separate program
  static bool  binaryViterbiSwap;
  static FILE *binaryViterbiFile;
  static char *binaryViterbiFilename;
  static gmtk_off_t binaryViterbiOffset;    // offset to start of current segment
  static gmtk_off_t nextViterbiOffset;      // offset to start of next segment

  // binary viterbi files should start with the cookie
#define GMTK_VITERBI_COOKIE        "GMTKVIT\n"
#define GMTK_VITERBI_COOKIE_NOLF   "GMTKVIT"
#define GMTK_VITERBI_COOKIE_LENGTH 8
  // cookie + BOM + k (for k-best) + # segements
#define GMTK_VITERBI_HEADER_SIZE (sizeof(unsigned) + \
                                  sizeof(unsigned) + \
                                  sizeof(unsigned) + \
                                  GMTK_VITERBI_COOKIE_LENGTH)


  // online filtering/smoothing needs to take some Viterbi code
  // paths but not others (particularly it should not allocate O(T)
  // memory for the Viterbi values, but it should call the Viterbi
  // versions of the MaxClique DE routines and setup the hidRVVector
  // in the PartitionStructures). JunctionTree::onlineViterbi is only
  // true in gmtkOnline so it can take the necessary code paths where
  // viterbiScore needs to be false to avoid the unwanted code paths.
  static bool onlineViterbi;

  // should printAllCliques() print scores or probabilities?
  static bool normalizePrintedCliques;





 protected:

  // The set of base sections from which real unrolled things are cloned from.
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
  // copies of C.  The next three variables hold sections P', C',
  // and E', where, for the standard left interface and simple
  // boundary case, we have that P' = P, C' = C , and E' = [C
  // E]. These sections use the *same* set of C++ random variables
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
  // are determined in create_base_partitions().
  vector <RV*> section_unrolled_rvs; 
  // mapping from name(frame) to integer index into unrolled_rvs.
  map < RVInfo::rvParent, unsigned > section_ppf;


  // The names of the above three sections to use for printing,
  // debugging, etc.
  static const char* P1_n;
  static const char* Co_n;
  static const char* E1_n;
  

  // generalize for factored interfaces - multiple "cliques"

  // Identities of cliques in junction trees: 
  // for P, 
  //    P's right  interface to C (a root in a JT section)
  unsigned P_ri_to_C; 
  //    The next one does not exist since we currently always do CE first.
  // unsigned P_li_clique; 
  // 
  // for C
  //    C's left interface to P
  unsigned C_li_to_P;
  //    C's left interface to C (same as C_li_to_P)
  unsigned C_li_to_C;
  //    C's right interface to C (a root in a JT section)
  unsigned C_ri_to_C;
  //    C's right interface to E (a root in a JT section) (same as C_ri_to_C)
  unsigned C_ri_to_E;
  // 
  // for E, 
  // E's left interface to C
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


  // within-section message ordering should be in a SectionInferenceAlgorithm

  // Message passing orders for each section.  Increasing index
  // order is 'collect evidence' phase from left to right in direction
  // of time, and decreasing order is 'distribute evidence' phase from
  // right to left in direction of time. Note that this assumes that
  // the overal root node in the JT is on the far right within E
  // (which might not be the best order).  NOTE: These are kept here
  // rather than in the sections, since they are re-used for all
  // cloned sections.
  vector< pair<unsigned,unsigned> > P1_message_order;
  vector< unsigned > P1_leaf_cliques;
  vector< pair<unsigned,unsigned> > Co_message_order;
  vector< unsigned > Co_leaf_cliques;
  vector< pair<unsigned,unsigned> > E1_message_order;  
  vector< unsigned > E1_leaf_cliques;


  // sets of random variables for frame shifting.

  // Trickiness: there are no uses of cur_{CC,CE}_rvs other than keeping
  //   them up-to-date as they're shifted. But the RVs in the sets are 
  //   from the SectionScheduler::section_structure_array, and their values
  //   are set by the adustFramesBy() calls as they're shifted.
  set <RV*> cur_CC_rvs;
  set <RV*> cur_CE_rvs;
  unsigned cur_cc_shift;
  unsigned cur_ce_shift;
  void shiftCCtoPosition(int pos, bool reset_observed=true);
  void shiftCCrelative(int delta) { shiftCCtoPosition(cur_cc_shift+delta); }
  void shiftCEtoPosition(int, bool reset_observed=true);
  void shiftOriginalVarstoPosition(vector<RV*> rvs, int pos, int &prevPos);
  void shiftCErelative(int delta) { shiftCEtoPosition(cur_cc_shift+delta); }
  void init_CC_CE_rvs(SectionIterator &it);
  void setCurrentInferenceShiftTo(SectionIterator &it, int pos);


  // Note, while we need extra separator cliques that are between the
  // corresponding sections interface cliques, these separators will
  // live in the section on the right of the separator. They will be
  // pointed to by the left interface clique in the section on the
  // right.


  // These will need adjusting per SectionInferenceAlgorithm...


  // The section structures that hold structures for a set of
  // RVs unrolled enough to cover any actual length DGM.
  sArray <PartitionStructures> section_structure_array;
  // The section tables that hold the actual clique/separator tables
  // for a DGM. This might be much longer than the
  // section_structure_array but is certainly no shorter.
  sArray <SectionTablesBase *> section_table_array;


  // range of cliques within each section to print out when doing
  // CE/DE inference. If these are NULL, then we print nothing.
  BP_Range* p_clique_print_range;
  BP_Range* c_clique_print_range;
  BP_Range* e_clique_print_range;

  Range section_debug_range;

  ObservationSource *obs_source; // & other common members?
  // inference task methods can dynamic cast to FileSource or StreamSource as needed?

  // the fixed file parser for this model, for RV unrolling, etc.
  FileParser   &fp;

  // The fixed gm_template for this model, contains the pre-triangulated graph.
  GMTemplate   &gm_template;


  // I think all below here belongs in SectionInferenceAlgorithm subclasses ?



  // 
  // Do some last-minute data structure setup to prepare for
  // unrolling to work (such as preliminary and pre work for
  // leaving STL, etc.)
  static void prepareForUnrolling(JT_Partition &section);
  void prepareForUnrolling();



  // A version of unroll that starts with the gm_template and fills up
  // base sections.
  void create_base_sections();
  void insertFactorClique(FactorClique& factorClique,FactorInfo& factor);


  // TODO: consider moving this to separate file?

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


 public:

  virtual void setUpJTDataStructures(const char* varSectionAssignmentPrior,
				     const char *varCliqueAssignmentPrior);

  // create the three junction trees for the basic sections.
  void createSectionJunctionTrees(const string pStr = junctionTreeMSTpriorityStr);

  // create a junction tree within a section.
  static void createSectionJunctionTree(Section& section, const string pStr = junctionTreeMSTpriorityStr);

  // routine to find the interface cliques of the sections
  void computeSectionInterfaces();

  // routine to create the factors in the appropriate sections
  void createFactorCliques();

  // routine to find the interface cliques of a section
  void computeSectionInterface(JT_Partition& section1,
			       unsigned int& section1_ric,
			       JT_Partition& section2,
			       unsigned int& section2_lic,
			       bool& icliques_same);


  // root the JT
  void createDirectedGraphOfCliques();
  static void createDirectedGraphOfCliquesRecurse(JT_Partition& section,
						  const unsigned root,
						  vector< bool >& visited);

  static void createDirectedGraphOfCliques(JT_Partition &section, const unsigned root);
  static void getCumulativeAssignedNodes(JT_Partition &section, const unsigned root);
  static void getCumulativeUnassignedIteratedNodes(JT_Partition &section,const unsigned root);


  static void setUpMessagePassingOrderRecurse(JT_Partition &section,
					      const unsigned root,
					      vector< pair<unsigned,unsigned> >&order,
					      const unsigned excludeFromLeafCliques,
					      vector< unsigned>& leaf_cliques);

  static void assignRVToClique(const char *const sectionName,
			       JT_Partition &section,
			       const unsigned root,
			       const unsigned depth,
			       RV* rv,
			       unsigned& numberOfTimesAssigned,
			       set<RV*>& parSet,
			       const bool allParentsObserved,
			       multimap< vector<double>, unsigned >& scoreSet);

  // Assign probability giving random variables to cliques (i.e.,
  // these are assigned only to cliques such that the random variables
  // and *all* their parents live in the clique, plus some other
  // criterion in order to make message passing as efficient as
  // possible).
  void assignRVsToCliques(const char* varSectionAssignmentPrior,
			  const char *varCliqueAssignmentPrior);

  static void assignRVsToCliques(const char *const sectionName,
				 JT_Partition&section,
				 const unsigned rootClique,
				 const char* varSectionAssignmentPrior,
				 const char *varCliqueAssignmentPrior);

  void assignFactorsToCliques();
  void assignFactorsToCliques(JT_Partition& section);

  // determine and set the unassignedNodes in each clique
  // in each partition
  void computeUnassignedCliqueNodes();
  static void computeUnassignedCliqueNodes(JT_Partition& part);

  // For the three sections, set up the different message passing
  // orders that are to be used. This basically just does a tree
  // traversal using the previously selected root.
  void setUpMessagePassingOrders();
  static void setUpMessagePassingOrder(JT_Partition& section,
				       const unsigned root,
				       vector< pair<unsigned,unsigned> >&order,
				       const unsigned excludeFromLeafCliques,
				       vector< unsigned>& leaf_cliques);

  // Separator creation, meaning create the seperator objects
  // both within and between sections. Given two neighboring
  // sections L and R, the separator between the interface
  // cliques in L and R is contained in R.
  static void createSeparators(JT_Partition& section, vector< pair<unsigned,unsigned> >&order);
  void createSeparators();

  // create the virtual evidence separators
  static void createVESeparators(JT_Partition& section);


  // Separator iteration order and accumulated set intersection
  // creation for separator driven clique potential creation, and
  // also updates the seperators partial accumulator structure and
  // sets up cliques other variables.
  static void computeSeparatorIterationOrder(MaxClique& clique, JT_Partition& section);
  static void computeSeparatorIterationOrders(JT_Partition& section);
  void computeSeparatorIterationOrders();

  // Computes the preceding iterated unassigned nodes and therein the
  // set of assigned nodes in each clique that should/shouldn't be
  // iterated.
  void getCumulativeUnassignedIteratedNodes();

  // compute the assignment order for nodes in this
  // section's cliques relative to each clique's incomming separators, and while
  // doing so, also set the dispositions for each of the resulting
  // nodes in each clique.
  void sortCliqueAssignedNodesAndComputeDispositions(const char *varCliqueAssignmentPrior);
  void sortCliqueAssignedNodesAndComputeDispositions(JT_Partition& section, const char *varCliqueAssignmentPrior);



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
				   

};

#endif
