/*
 * GMTK_LinearSectionScheduler.cc
 *
 * Written by Richard Rogers <rprogers@uw.edu>
 *
 * Copyright (C) 2015 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

#include "GMTK_LinearSectionScheduler.h"

  // Initialize stuff at the model-level. See prepareForSegment() for segment-level initialization.
  // TODO: explain parameters
void 
LinearSectionScheduler::setUpDataStructures(char const *varSectionAssignmentPrior,
				   char const *varCliqueAssignmentPrior)
{
}

  // Formerly JunctionTree::printAllJTInfo()
void 
LinearSectionScheduler::printInferencePlanSummary(char const *fileName) {
}

  // Formerly GMTemplate::reportScoreStats()
void 
LinearSectionScheduler::reportScoreStats() {
}


  // NOTE:  This assumes all inference algorithms will have clique-like things they
  //        need to print posteriors of

  // Set the range of selected clique #'s in P', C', E' for printing.
  // TODO: preconditions
void 
LinearSectionScheduler::setCliquePrintRanges(char *p_range, char *c_range, char *e_range) {
}

  // Print to f the order of the variables in each clique selected by setCliquePrintRanges().
void 
LinearSectionScheduler::printCliqueOrders(FILE *f) {
}

  // Returns the size (in # of floats) of the cliques selected by setCliquePrintRanges().
void 
LinearSectionScheduler::getCliquePosteriorSize(unsigned &p_size, unsigned &c_size, unsigned &e_size) {
}

