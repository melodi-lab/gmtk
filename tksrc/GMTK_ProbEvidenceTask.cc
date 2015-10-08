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

#include <assert.h>

#include "GMTK_ProbEvidenceTask.h"

// This method is here just to make the compiler happy.
// It should be overridden by the SectionScheduler subclass.

logpr 
ProbEvidenceTask::probEvidence(SectionInferenceAlgorithm *algorithm,
			       unsigned *numUsableFrames,
			       unsigned *numSectionsDone,
			       const bool limitTime,
			       const bool noE,
			       const bool cliquePosteriorNormalize,
			       const bool cliquePosteriorUnlog,
			       ObservationFile *posteriorFile)
{
  assert(false); // nothing should call this method
  // just here to make the vtable non-empty
  logpr result;
  return result;
}
