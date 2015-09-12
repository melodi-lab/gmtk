/*
 * GMTK_ArchipelagosSectionScheduler.h
 *   
 *
 * Written by Richard Rogers <rprogers@uw.edu>
 *
 * Copyright (C) 2015 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */


#ifndef GMTK_ARCHIPELAGOSSECTIONSCHEDULER_H
#define GMTK_ARCHIPELAGOSSECTIONSCHEDULER_H

#include "GMTK_FileSource.h"

#include "GMTK_GMTemplate.h"

#include "GMTK_ProbEvidenceTask.h"
#include "GMTK_ForwardBackwardTask.h"
#include "GMTK_ViterbiTask.h"
#include "GMTK_SmoothingTask.h"

class ArchipelagosSectionScheduler : public ProbEvidenceTask, ForwardBackwardTask, ViterbiTask {

 public:

  ArchipelagosSectionScheduler(GMTemplate &gm_template, FileSource &observation_source) {}

  ~ArchipelagosSectionScheduler() {}

};

#endif

