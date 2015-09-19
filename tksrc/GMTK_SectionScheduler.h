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

#include "logp.h"

// TODO: what's the difference between range and bp_range?
#include "range.h"
#include "bp_range.h"

#include "GMTK_SectionInferenceAlgorithm.h"

class SectionScheduler {

 public:

  virtual ~SectionScheduler() {}

  // Initialize stuff at the model-level. See prepareForSegment() for segment-level initialization.
  // TODO: explain parameters
  virtual void setUpDataStructures(char const *varSectionAssignmentPrior,
				   char const *varCliqueAssignmentPrior) = 0;


  // Prepare to do inference on segment of length T. Note that this
  // is the analog of the current JunctionTree::unroll() which sets up the
  // PartitionStructureArray for at most  P' C' C' E' (4 sections, so
  // O(1) memory), and the PartitionTableArray for between 0 and T
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
  virtual unsigned unroll(unsigned T) = 0; // returns number of usable frames


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
  virtual void setSectionDebugRange(Range &rng) {
    section_debug_range.SetLimits(rng.first(), rng.last()); 
    section_debug_range.SetDefStr(rng.GetDefStr()); 
  }

 protected:

  // range of cliques within each section to print out when doing
  // CE/DE inference. If these are NULL, then we print nothing.
  BP_Range* pCliquePrintRange;
  BP_Range* cCliquePrintRange;
  BP_Range* eCliquePrintRange;

  Range section_debug_range;

};

#endif
