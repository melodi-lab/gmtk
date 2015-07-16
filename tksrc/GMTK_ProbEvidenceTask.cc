/*
 * GMTK_ProbEvidenceTask.cc
 *   compute P(X_{0:T-1}), forward pass only, O(1) memory
 *
 * Written by Richard Rogers <rprogers@uw.edu>
 *
 * Copyright (C) 2015 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

#include "GMTK_ProbEvidenceTask.h"

logpr 
ProbEvidenceTask::probEvidence(unsigned *numUsableFrames,
			       unsigned *numSectionsDone,
			       const bool limitTime,
			       const bool noE,
			       const bool cliquePosteriorNormalize,
			       const bool cliquePosteriorUnlog,
			       ObservationFile *posteriorFile)
{
  logpr result;
  return result;
}
