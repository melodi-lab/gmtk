/*
 * GMTK_SectionScheduler.h
 *   Root class for the inference algorithms at the time series level.
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

#include <vector>
#include <map>

#include "logp.h"

// TODO: what's the difference between range and bp_range?
#include "range.h"
#include "bp_range.h"
#include "mArray.h"

#include "fileParser.h"

#include "GMTK_ObservationSource.h"
#include "GMTK_FileParser.h"

#include "GMTK_PtPsIterator.h"
#include "GMTK_PartitionStructures.h"
#include "GMTK_PartitionTables.h"
#include "GMTK_JT_Partition.h"
#include "GMTK_PtPsIterator.h"
#include "GMTK_SectionInferenceAlgorithm.h"

class SectionScheduler {

 public:

#if 0
  SectionScheduler(GMTemplate &gm_template, FileParser &fp, SectionInferenceAlgorithm *algorithm, ObservationSource *obs_source) :
    pCliquePrintRange(NULL), cCliquePrintRange(NULL), eCliquePrintRange(NULL),
    section_debug_range("all", 0, 0x7FFFFFFF), obs_source(obs_source), fp(fp), gm_template(gm_template), algorithm(algorithm)
  {}
#else
    JunctionTree *extremely_bogus_jt; // BOGUS
    
  SectionScheduler(GMTemplate &gm_template, FileParser &fp, SectionInferenceAlgorithm *algorithm, ObservationSource *obs_source) :
    extremely_bogus_jt(new JunctionTree(gm_template)), // TOTALLY BOGUS
    pCliquePrintRange(NULL), cCliquePrintRange(NULL), eCliquePrintRange(NULL),
    section_debug_range("all", 0, 0x7FFFFFFF), obs_source(obs_source), fp(fp), gm_template(gm_template),
      inference_it(*extremely_bogus_jt), algorithm(algorithm)
  {}
#endif
    
  virtual ~SectionScheduler() {}

  virtual JunctionTree *getJT() = 0; // BOGUS

  // Initialize stuff at the model-level. See prepareForSegment() for segment-level initialization.
  // TODO: explain parameters
  virtual void setUpDataStructures(iDataStreamFile &tri_file,
				   char const *varSectionAssignmentPrior,
				   char const *varCliqueAssignmentPrior,
				   bool checkTriFileCards) = 0;


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
  virtual void printInferencePlanSummary(char const *fileName) = 0;

  // Formerly GMTemplate::reportScoreStats()
  virtual void reportScoreStats() = 0;


  // NOTE:  This assumes all inference algorithms will have clique-like things they
  //        need to print posteriors of

  // Set the range of selected clique #'s in P', C', E' for printing.
  // TODO: preconditions
  virtual void setCliquePrintRanges(char *p_range, char *c_range, char *e_range) = 0;

  // Print to f the order of the variables in each clique selected by setCliquePrintRanges().
  virtual void printCliqueOrders(FILE *f) = 0;

  // Returns the size (in # of floats) of the cliques selected by setCliquePrintRanges().
  virtual void getCliquePosteriorSize(unsigned &p_size, unsigned &c_size, unsigned &e_size) = 0;

  // Set the range of section #s that get elevated verbosity
  // P' = 0, C' \in 1, ..., T=2, E' = T-1
  virtual void setSectionDebugRange(Range const &rng) {
    section_debug_range.SetLimits(rng.first(), rng.last()); 
    section_debug_range.SetDefStr(rng.GetDefStr()); 
  }

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
  
  // Note, while we need extra separator cliques that are between the
  // corresponding sections interface cliques, these separators will
  // live in the section on the right of the separator. They will be
  // pointed to by the left interface clique in the section on the
  // right.

  // The section structures that hold structures for a set of
  // RVs unrolled enough to cover any actual length DGM.
  sArray <PartitionStructures> section_structure_array;
  // The section tables that hold the actual clique/separator tables
  // for a DGM. This might be much longer than the
  // section_structure_array but is certainly no shorter.
  sArray <PartitionTables> section_table_array;


  // range of cliques within each section to print out when doing
  // CE/DE inference. If these are NULL, then we print nothing.
  BP_Range* pCliquePrintRange;
  BP_Range* cCliquePrintRange;
  BP_Range* eCliquePrintRange;

  Range section_debug_range;

  ObservationSource *obs_source; // & other common members?
  // inference task methods can dynamic cast to FileSource or StreamSource as needed?

  // the fixed file parser for this model, for RV unrolling, etc.
  FileParser   &fp;

  // The fixed gm_template for this model, contains the pre-triangulated graph.
  GMTemplate   &gm_template;

  PtPsIterator inference_it;
  SectionInferenceAlgorithm *algorithm;
};

#endif
