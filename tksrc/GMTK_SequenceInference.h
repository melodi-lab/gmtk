/*
 * GMTK_SeriesnIference.h
 *   Just go sequentially through the sections applying the SectionInferenceAlgorithm
 *
 * Written by Richard Rogers <rprogers@uw.edu>
 *
 * Copyright (C) 2015 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */


#ifndef GMTK_SEQUENCEINFERENCE_H
#define GMTK_SEQUENCEINFERENCE_H

#include "GMTK_FileSource.h"

#include "GMTK_GMTemplate.h"

#include "GMTK_ProbEvidenceTask.h"
#include "GMTK_ForwardBackwardTask.h"
#include "GMTK_ViterbiTask.h"

class SequenceInference : public ProbEvidenceTask, ForwardBackwardTask, ViterbiTask {

 public:

  SequenceInferece(GMTemplate &gm_template, FileSource &observation_source) {}

  ~SequenceInference() {}

 protected:

};

#endif

