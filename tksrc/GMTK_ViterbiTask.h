/*
 * GMTK_ViterbiTask.h
 *    compute argmax_{Q_{0:T-1}} P(Q_{0:T-1} | X_{0:T-1})
 *
 * Written by Richard Rogers <rprogers@uw.edu>
 *
 * Copyright (C) 2015 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */


#ifndef GMTK_VITERBITASK_H
#define GMTK_VITERBITASK_H

#include "GMTK_FileSource.h"

class ViterbiTask {

 public:

  virtual ~ViterbiTask() {}

  /*
   * For the currently active segment of the observation_source, compute 
   * argmax_{Q_{0:T-1}} P(Q_{0:T-1} | X_{0:T-1})
   *
   * observation_source is the source of the observed data
   * 
   * numUsableFrames returns the number of frames in the currently active
   *                 segment that were used for inference
   *
   * if cliquePosteriorNormalize is true, normalize the posterior scores to sum to 1 (in log space)
   *
   * if cliquePosteriorUnlog is true, return the posterior probability (vs. log(probability))
   * 
   * if posteriorFile is non-NULL, write the clique posteriors to the posteriorFile
   *
   */

  virtual logpr viterbi(FileSource &observation_source,
                        unsigned *numUsableFrames = NULL,
			const bool cliquePosteriorNormalize = true,
			const bool cliquePosteriorUnlog = true,
		       	ObservationFile *posteriorFile = NULL) = 0;
  
};

#endif
