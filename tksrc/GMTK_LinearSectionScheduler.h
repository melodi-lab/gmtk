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

#include "GMTK_FileSource.h"
#include "GMTK_StreamSource.h"

#include "GMTK_GMTemplate.h"

#include "GMTK_SectionScheduler.h"

#include "GMTK_ProbEvidenceTask.h"
#include "GMTK_ForwardBackwardTask.h"
#include "GMTK_ViterbiTask.h"
#include "GMTK_SmoothingTask.h"

#include "GMTK_SectionInferenceAlgorithm.h"

class LinearSectionScheduler : public SectionScheduler, ProbEvidenceTask, ForwardBackwardTask, ViterbiTask, SmoothingTask {

 public:

  LinearSectionScheduler(GMTemplate                &gm_template,
			 SectionInferenceAlgorithm *algorithm,
			 FileSource                *observation_file = NULL,
			 StreamSource              *observation_stream = NULL)
    : gm_template(gm_template), algorithm(algorithm),
      observation_file(observation_file), observation_stream(observation_stream)
  {
    assert(algorithm);
    assert(observation_file || observation_stream);
  }

  virtual ~LinearSectionScheduler() {}


  // Initialize stuff at the model-level. See prepareForSegment() for segment-level initialization.
  // TODO: explain parameters
  void setUpDataStructures(char const *varSectionAssignmentPrior,
				   char const *varCliqueAssignmentPrior);

  // Formerly JunctionTree::printAllJTInfo()
  void printInferencePlanSummary(char const *fileName);

  // Formerly GMTemplate::reportScoreStats()
  void reportScoreStats();


  // NOTE:  This assumes all inference algorithms will have clique-like things they
  //        need to print posteriors of

  // Set the range of selected clique #'s in P', C', E' for printing.
  // TODO: preconditions
  void setCliquePrintRanges(char *p_range, char *c_range, char *e_range);

  // Print to f the order of the variables in each clique selected by setCliquePrintRanges().
  void printCliqueOrders(FILE *f);

  // Returns the size (in # of floats) of the cliques selected by setCliquePrintRanges().
  void getCliquePosteriorSize(unsigned &p_size, unsigned &c_size, unsigned &e_size);

  logpr probEvidence(unsigned *numUsableFrames = NULL,
		     unsigned *numSectionsDone = NULL,
		     const bool limitTime = false,
		     const bool noE = false, 
		     const bool cliquePosteriorNormalize = true,
		     const bool cliquePosteriorUnlog = true,
		     ObservationFile *posteriorFile = NULL);

  logpr forwardBackward(unsigned *numUsableFrames = NULL,
			const bool cliquePosteriorNormalize = true,
			const bool cliquePosteriorUnlog = true,
			ObservationFile *posteriorFile = NULL);

  logpr viterbi(unsigned *numUsableFrames = NULL,
		const bool cliquePosteriorNormalize = true,
		const bool cliquePosteriorUnlog = true,
		ObservationFile *posteriorFile = NULL);

  logpr smoothing(unsigned *numUsableFrames = NULL,
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

  GMTemplate   &gm_template;
  SectionInferenceAlgorithm *algorithm;
  FileSource   *observation_file;
  StreamSource *observation_stream;
  
};

#endif

