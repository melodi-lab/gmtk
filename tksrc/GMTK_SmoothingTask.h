/*
 * GMTK_SmoothingTask.h
 *    compute   argmax_{Q_{t-\tau}} P(Q_{t-\tau} | X_{0:t}) 
 *    for \tau = 0 (filtering) or \tau > 0 (smoothing)
 *
 * Written by Richard Rogers <rprogers@uw.edu>
 *
 * Copyright (C) 2015 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */


#ifndef GMTK_SMOOTHINGTASK_H
#define GMTK_SMOOTHINGTASK_H

#include "GMTK_StreamSource.h"

class SmoothingTask {

 public:

  virtual ~SmoothingTask() {}

  /*
   * For the currently active segment of the observation_source, compute 
   *  argmax_{Q_{t-\tau}} P(Q_{t-\tau} | X_{0:t}) 
   *
   * numUsableFrames returns the number of frames in the currently active
   *                 segment that were used for inference
   *
   * numSectionsDone returns the number of modified sections in the currently
   *                 active segment that were used for inference
   * 
   * if noE is true, do not perform inference on E' (return results from up
   *                 to the last C'
   *
   * Viterbi-style output is written to the file f
   * 
   * if printObserved is true, included the values of observed variables in posterior output
   * 
   * preg, creg, ereg are regular expressions for selecting the names of random variables
   *                  to include in the Viterbi-style output
   * 
   * if cliquePosteriorNormalize is true, normalize the posterior scores to sum to 1 (in log space)
   *
   * if cliquePosteriorUnlog is true, return the posterior probability (vs. log(probability))
   * 
   * if posteriorFile is non-NULL, write the clique posteriors to the posteriorFile
   *
   */
  virtual logpr smoothing(StreamSource &observation_source,
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
};

#endif
