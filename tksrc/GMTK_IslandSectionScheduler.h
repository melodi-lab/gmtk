/*
 * GMTK_IslandSectionScheduler.h
 *   
 *
 * Written by Richard Rogers <rprogers@uw.edu>
 *
 * Copyright (C) 2015 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */


#ifndef GMTK_ISLANDSECTIONSCHEDULER_H
#define GMTK_ISLANDSECTIONSCHEDULER_H

#include "GMTK_SectionScheduler.h"

#include "GMTK_ProbEvidenceTask.h"
#include "GMTK_ForwardBackwardTask.h"
#include "GMTK_ViterbiTask.h"
#include "GMTK_SmoothingTask.h"

class IslandSectionScheduler : public SectionScheduler,
                               public ProbEvidenceTask, // just for completeness, LinearSectionScheduler is faster and takes less memory
                               public ForwardBackwardTask, 
                               public ViterbiTask 
{

 public:

  IslandSectionScheduler(GMTemplate                &gm_template, 
			 FileParser                &fp, 
			 ObservationSource         *obs_source)
    : SectionScheduler(gm_template, fp, obs_source)
  {}

  virtual ~IslandSectionScheduler() {}

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
};

#endif

