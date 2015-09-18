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


logpr
LinearSectionScheduler::probEvidence(unsigned *numUsableFrames,
				     unsigned *numSectionsDone,
				     const bool limitTime,
				     const bool noE, 
				     const bool cliquePosteriorNormalize,
				     const bool cliquePosteriorUnlog,
				     ObservationFile *posteriorFile)
{
  unsigned T; // # of sections
  unsigned nUsableFrames = algorithm->unroll(observation_file->numFrames() /*, ZeroTable, &T*/);
  if (numUsableFrames) *numUsableFrames = nUsableFrames;

  // do P'
  InterfaceSeparator msg = algorithm->computeForwardInterfaceSeparator(0);

  // do C'
  unsigned t;
  for (t=1; t < T-1; t+=1) {
    algorithm->receiveForwardInterfaceSeparator(t, msg);
    msg = algorithm->computeForwardInterfaceSeparator(t);
    //if (limitTime && probEvidenceTimeExpired) goto finished;
  }

  // do E'
  algorithm->receiveForwardInterfaceSeparator(t, msg);
  algorithm->computeForwardInterfaceSeparator(t);

 finished:
  
  if (numSectionsDone) *numSectionsDone = t;
  return algorithm->probEvidence(t);
}
  
logpr 
LinearSectionScheduler::forwardBackward(unsigned *numUsableFrames,
					const bool cliquePosteriorNormalize,
					const bool cliquePosteriorUnlog,
					ObservationFile *posteriorFile)
{
  logpr result;
  return result;
}

logpr 
LinearSectionScheduler::viterbi(unsigned *numUsableFrames,
				const bool cliquePosteriorNormalize,
				const bool cliquePosteriorUnlog,
				ObservationFile *posteriorFile)
{
  logpr result;
  return result;
}


logpr 
LinearSectionScheduler::smoothing(unsigned *numUsableFrames,
				  unsigned *numSectionsDone,
				  const bool noE,
				  FILE *f,
				  const bool printObserved,
				  regex_t *preg,
				  regex_t *creg,
				  regex_t *ereg,
				  ObservationFile *posteriorFile,
				  const bool cliquePosteriorNormalize,
				  const bool cliquePosteriorUnlog)
{
  logpr result;
  return result;
}
