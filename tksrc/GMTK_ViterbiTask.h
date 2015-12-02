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

#include "logp.h"
#include "mArray.h"
#include "sArray.h"

#include "GMTK_FileSource.h"

#include "GMTK_SectionIterator.h"
#include "GMTK_SectionInferenceAlgorithm.h"

class ViterbiTask {

 public:

  static ObservationFile *vitObsFile;
  static char *vitObsFileName;
  static char *vitObsListName;
  static const char *vitObsNameSeparator;
  static const char *vitObsFileFmt;
  static bool  vitObsFileSwap;
  vector<string> vitObsVariableNames;

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

  virtual logpr viterbi(SectionInferenceAlgorithm *algorithm,
			unsigned nBest=1,
			unsigned *numUsableFrames = NULL,
			const bool cliquePosteriorNormalize = true,
			const bool cliquePosteriorUnlog = true,
		       	ObservationFile *posteriorFile = NULL);

 protected:

  ////////////////////////////////////////////////////////////////////////
  // Support variables specific to Viterbi and N-best decoding
  // 
  sArray < unsigned > P_section_values;
  mArray < unsigned > C_section_values;
  sArray < unsigned > E_section_values;

  void recordSectiontionViterbiValue(SectionIterator &it);

  ////////////////////////////////////////////////////////////////////////

  
};

#endif
