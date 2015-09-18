/*
 * GMTK_ForwardBackwardTask.cc
 *   compute P(Q_{0:T-1} | X_{0:T-1})
 *
 * Written by Richard Rogers <rprogers@uw.edu>
 *
 * Copyright (C) 2015 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

#include <assert.h>

#include "GMTK_ForwardBackwardTask.h"

// This method is here just to make the compiler happy.
// It should be overridden by the SectionScheduler subclass.

logpr 
ForwardBackwardTask::forwardBackward(unsigned *numUsableFrames,
				     const bool cliquePosteriorNormalize,
				     const bool cliquePosteriorUnlog,
				     ObservationFile *posteriorFile)
{
  assert(false); // nothing should call this method
  logpr result;
  return result;
}
