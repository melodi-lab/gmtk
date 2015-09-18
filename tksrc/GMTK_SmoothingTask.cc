/*
 * GMTK_SmoothingTask.cc
 *   compute argmax_{Q_{t-\tau}} P(Q_{t-\tau} | X_{0:t}) 
 *
 * Written by Richard Rogers <rprogers@uw.edu>
 *
 * Copyright (C) 2015 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

#include <assert.h>
#include <regex.h>

#include "GMTK_ObservationFile.h"

#include "GMTK_SmoothingTask.h"

// This method is here just to make the compiler happy.
// It should be overridden by the SectionScheduler subclass.

logpr 
SmoothingTask::smoothing(unsigned *numUsableFrames,
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
  assert(false); // nothing should call this method
  logpr result;
  return result;
}
