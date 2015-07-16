/*
 * GMTK_SeriesInferenceAlgorithm.h
 *   Root class for the inference algorithms at the time series level.
 *
 * Written by Richard Rogers <rprogers@uw.edu>
 *
 * Copyright (C) 2015 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 * SeriesInferenceAlgorithm subclasses handle inference over a series of sections. 
 *
 * The specific inference task subclasses of SeriesInferenceAlgorithm are 
 * (X is observed, Q is hidden, T is segment length):
 *
 *  - ProbEvidenceTask       compute P(X_{0:T-1}), forward pass only, O(1) memory
 *  - ViterbiTask            compute argmax_{Q_{0:T-1}} P(Q_{0:T-1} | X_{0:T-1})
 *  - ForwardBackwardTask    compute P(Q_{0:T-1} | X_{0:T-1})
 *  - OnlineSmoothingTask    compute argmax_{Q_t} P(Q_t | X_{0:t+\tau})
 *
 * The time-series algorithm subclasses of SeriesInferenceAlgorithm are:
 *  - SequenceInference      Just iterate sequentially over sections
 *  - IslandInference        Skip through the time series, leaving "islands" of partial results to save memory
 *  - ArchipelagosInference  Parallel version of Island
 * 
 * The *Inference subclasses multiply inherit the *Task classes representing the
 * inference tasks the time series level inference algorithm can perform. Not all 
 * SeriesInferenceAlgorithms support every inference task, e.g., IslandInference 
 * can't do online filtering/smoothing.
 *
 * The Sections are responsible for doing inference within themselves via whatever
 * SectionInferenceAlgorithm the Section subclass implements. The SeriesInferenceAlgorithm
 * classes just pass messages between Sections via instances of InterfaceSeparator.
 *
 */


#ifndef GMTK_SERIESINFERENCEALGORITHM_H
#define GMTK_SERIESINFERENCEALGORITHM_H

#include "logp.h"

// TODO: what's the difference between range and bp_range?
#include "range.h"
#include "bp_range.h"

#include "GMTK_SectionInferenceAlgorithm.h"

class SeriesInferenceAlgorithm {

 public:

  SeriesInferenceAlgorithm() {}

  virtual ~SeriesInferenceAlgorithm() {}

  // Initialize stuff at the model-level. See prepareForSegment() for segment-level initialization.
  // TODO: explain parameters
  virtual void setUpDataStructures(char const *varPartitionAssignmentPrior,
				   char const *varCliqueAssignmentPrior) = 0;

  // Formerly JunctionTree::printAllJTInfo()
  virtual void printInferencePlanSummary(char const *fileName) = 0;

  // Formerly GMTemplate::reportScoreStats()
  virtual void reportScoreStats();


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

  // range of cliques within each partition to print out when doing
  // CE/DE inference. If these are NULL, then we print nothing.
  BP_Range* pPartCliquePrintRange;
  BP_Range* cPartCliquePrintRange;
  BP_Range* ePartCliquePrintRange;

  Range section_debug_range;

};

#endif
