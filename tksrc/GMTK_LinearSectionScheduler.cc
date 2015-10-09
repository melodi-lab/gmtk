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


logpr
LinearSectionScheduler::probEvidence(SectionInferenceAlgorithm *algorithm,
				     unsigned *numUsableFrames,
				     unsigned *numSectionsDone,
				     const bool limitTime,
				     const bool noE, 
				     const bool cliquePosteriorNormalize,
				     const bool cliquePosteriorUnlog,
				     ObservationFile *posteriorFile)
{
  assert(algorithm);
  FileSource *observation_source = dynamic_cast<FileSource *>(obs_source);
  assert(observation_source);

  unsigned T; // # of sections

  unsigned nUsableFrames = unroll(observation_source->numFrames(), ZeroTable, &T);
  if (numUsableFrames) *numUsableFrames = nUsableFrames;

  SectionIterator inference_it(*this,T);

  init_CC_CE_rvs(inference_it);
  
  PartitionTables *prev_sect_tab = NULL;
  PartitionTables *cur_sect_tab = new PartitionTables(inference_it.cur_jt_section());
  
 // do P'
  SectionSeparator *msg = algorithm->computeForwardInterfaceSeparator(inference_it, cur_sect_tab);

  // do C' C' ... E'
  unsigned t;
  for (t=1; t < T; t+=1) {

    setCurrentInferenceShiftTo(inference_it, t);

    delete prev_sect_tab;
    prev_sect_tab = cur_sect_tab; // msg points into prev_sect_tab now
    cur_sect_tab = new PartitionTables(inference_it.cur_jt_section());

    algorithm->receiveForwardInterfaceSeparator(inference_it, msg, cur_sect_tab);
    msg = algorithm->computeForwardInterfaceSeparator(inference_it, cur_sect_tab);
    // msg points into cur_sect_tab now

    //if (limitTime && probEvidenceTimeExpired) goto finished;
  }

  //finished:
  
  if (numSectionsDone) *numSectionsDone = t;
  logpr probE = algorithm->probEvidence(inference_it, cur_sect_tab);
  delete cur_sect_tab;
  return probE;
}
  
logpr 
LinearSectionScheduler::forwardBackward(SectionInferenceAlgorithm *algorithm,
					unsigned *numUsableFrames,
					const bool cliquePosteriorNormalize,
					const bool cliquePosteriorUnlog,
					ObservationFile *posteriorFile)
{
  assert(algorithm);
  logpr result;
  return result;
}

logpr 
LinearSectionScheduler::viterbi(SectionInferenceAlgorithm *algorithm,
				unsigned nBest,
				unsigned *numUsableFrames,
				const bool cliquePosteriorNormalize,
				const bool cliquePosteriorUnlog,
				ObservationFile *posteriorFile)
{
  assert(algorithm);
  logpr result;
  return result;
}


logpr 
LinearSectionScheduler::smoothing(SectionInferenceAlgorithm *algorithm,
				  unsigned nBest,
				  unsigned *numUsableFrames,
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
  assert(algorithm);
  logpr result;
  return result;
}
