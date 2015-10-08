/*
 * GMTK_LinearSectionScheduler.h
 *   Just go sequentially through the sections in time order, 
 *   applying the SectionInferenceAlgorithm
 *
 * Written by Richard Rogers <rprogers@uw.edu>
 *
 * Copyright (C) 2015 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */


#ifndef GMTK_LINEARSECTIONSCHEDULER_H
#define GMTK_LINEARSECTIONSCHEDULER_H

#include <assert.h>

#include "logp.h"

#include "GMTK_ObservationSource.h"

#include "GMTK_GMTemplate.h"
#include "GMTK_FileParser.h"

#include "GMTK_SectionScheduler.h"

#include "GMTK_ProbEvidenceTask.h"
#include "GMTK_ForwardBackwardTask.h"
#include "GMTK_ViterbiTask.h"
#include "GMTK_SmoothingTask.h"

#include "GMTK_SectionInferenceAlgorithm.h"


class LinearSectionScheduler : public SectionScheduler, 
                               public ProbEvidenceTask, 
                               public ForwardBackwardTask, 
                               public ViterbiTask, 
                               public SmoothingTask 
{

 public:

  LinearSectionScheduler(GMTemplate                &gm_template,
                         FileParser                &fp,
                         ObservationSource         *obs_source)
    : SectionScheduler(gm_template, fp, obs_source)
  {}

  virtual ~LinearSectionScheduler(); // putting dtor in .cpp seems to help with 'missing vtable' link errors...
  
  
  // Formerly JunctionTree::printAllJTInfo()
  void printInferencePlanSummary(char const *fileName);

  // Formerly GMTemplate::reportScoreStats()
  void reportScoreStats();


  // NOTE:  This assumes all inference algorithms will have clique-like things they
  //        need to print posteriors of

  // Set the range of selected clique #'s in P', C', E' for printing.
  // TODO: preconditions
  void setCliquePrintRanges(char *p_range, char *c_range, char *e_range);

  // Set range of sections that should produce extra debug info
  void setSectionDebugRange(Range &rng);

  // Print to f the order of the variables in each clique selected by setCliquePrintRanges().
  void printCliqueOrders(FILE *f);

  // Returns the size (in # of floats) of the cliques selected by setCliquePrintRanges().
  void getCliquePosteriorSize(unsigned &p_size, unsigned &c_size, unsigned &e_size);

  logpr probEvidence(SectionInferenceAlgorithm *algorithm,
		     unsigned *numUsableFrames = NULL,
		     unsigned *numSectionsDone = NULL,
		     const bool limitTime = false,
		     const bool noE = false, 
		     const bool cliquePosteriorNormalize = true,
		     const bool cliquePosteriorUnlog = true,
		     ObservationFile *posteriorFile = NULL);

  logpr forwardBackward(SectionInferenceAlgorithm *algorithm,
			unsigned *numUsableFrames = NULL,
			const bool cliquePosteriorNormalize = true,
			const bool cliquePosteriorUnlog = true,
			ObservationFile *posteriorFile = NULL);

  logpr viterbi(SectionInferenceAlgorithm *algorithm,
		unsigned nBest = 1,
		unsigned *numUsableFrames = NULL,
		const bool cliquePosteriorNormalize = true,
		const bool cliquePosteriorUnlog = true,
		ObservationFile *posteriorFile = NULL);

  logpr smoothing(SectionInferenceAlgorithm *algorithm,
		  unsigned nBest = 1,
		  unsigned *numUsableFrames = NULL,
		  unsigned *numSectionsDone=NULL,
		  const bool noE=false,
		  FILE *f=stdout,
		  const bool printObserved=false,
		  regex_t *preg=NULL,
		  regex_t *creg=NULL,
		  regex_t *ereg=NULL,
		  ObservationFile *posteriorFile = NULL,
		  const bool cliquePosteriorNormalize = true,
		  const bool cliquePosteriorUnlog = true);

 private:

};

#endif

