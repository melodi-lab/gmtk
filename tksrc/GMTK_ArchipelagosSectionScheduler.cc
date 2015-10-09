/*
 * GMTK_ArchipelagosSectionScheduler.cc
 *   
 *
 * Written by Richard Rogers <rprogers@uw.edu>
 *
 * Copyright (C) 2015 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

#include "GMTK_ArchipelagosSectionScheduler.h"

ArchipelagosSectionScheduler::~ArchipelagosSectionScheduler() {}

  
logpr
ArchipelagosSectionScheduler::probEvidence(SectionInferenceAlgorithm *algorithm,
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
ArchipelagosSectionScheduler::forwardBackward(SectionInferenceAlgorithm *algorithm,
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
ArchipelagosSectionScheduler::viterbi(SectionInferenceAlgorithm *algorithm,
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

