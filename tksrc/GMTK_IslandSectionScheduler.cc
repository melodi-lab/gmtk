/*
 * GMTK_IslandSectionScheduler.cc
 *   
 *
 * Written by Richard Rogers <rprogers@uw.edu>
 *
 * Copyright (C) 2015 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

#include "GMTK_IslandSectionScheduler.h"

IslandSectionScheduler::~IslandSectionScheduler() {}

  
logpr
IslandSectionScheduler::probEvidence(SectionInferenceAlgorithm *algorithm,
					   unsigned *numUsableFrames,
					   unsigned *numSectionsDone,
					   const bool limitTime,
					   const bool noE, 
					   const bool cliquePosteriorNormalize,
					   const bool cliquePosteriorUnlog,
					   ObservationFile *posteriorFile)
{
  assert(algorithm);
  logpr result;
  return result;
}

logpr 
IslandSectionScheduler::forwardBackward(SectionInferenceAlgorithm *algorithm,
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
IslandSectionScheduler::viterbi(SectionInferenceAlgorithm *algorithm,
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

