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

#include "GMTK_FileSource.h"
#include "GMTK_StreamSource.h"



LinearSectionScheduler::~LinearSectionScheduler() {}


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

void 
LinearSectionScheduler::setSectionDebugRange(Range &rng) {
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
  FileSource *observation_source = dynamic_cast<FileSource *>(obs_source);
  assert(observation_source);

  unsigned T; // # of sections

  // MOVE UNROLL TO SectionScheduler
  unsigned nUsableFrames = unroll(observation_source->numFrames(), ZeroTable, &T);
  if (numUsableFrames) *numUsableFrames = nUsableFrames;

  SectionIterator inference_it(*this,T);
  //myjt->sparseJoinSegementInit(T);        // BOGUS

  init_CC_CE_rvs(inference_it);
  
  PartitionTables *cur_sect_tab = new PartitionTables(inference_it.cur_jt_section());
  
 // do P'
  SectionSeparator *msg = algorithm->computeForwardInterfaceSeparator(0, cur_sect_tab);

  // do C' C' ... E'
  unsigned t;
  for (t=1; t < T; t+=1) {
    delete cur_sect_tab;

    setCurrentInferenceShiftTo(inference_it, t);

    cur_sect_tab = new PartitionTables(inference_it.cur_jt_section());
    algorithm->receiveForwardInterfaceSeparator(inference_it, msg, cur_sect_tab);
    delete msg;
    msg = algorithm->computeForwardInterfaceSeparator(t, cur_sect_tab);
    //if (limitTime && probEvidenceTimeExpired) goto finished;
  }

  // do E'
  delete msg; // not if E' is the first section!

  //finished:
  
  if (numSectionsDone) *numSectionsDone = t;
  logpr probE = algorithm->probEvidence(t, cur_sect_tab);
  delete cur_sect_tab;
  return probE;
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
