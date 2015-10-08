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

#include "GMTK_SectionScheduler.h"

#include "GMTK_ProbEvidenceTask.h"
#include "GMTK_ForwardBackwardTask.h"
#include "GMTK_ViterbiTask.h"
#include "GMTK_SmoothingTask.h"

class ArchipelagosSectionScheduler : public SectionScheduler,
                                     public ProbEvidenceTask, // just for completeness, LinearSectionScheduler is faster and takes less memory
                                     public ForwardBackwardTask, 
                                     public ViterbiTask 
{

 public:

  ArchipelagosSectionScheduler(GMTemplate                &gm_template, 
                               FileParser                &fp, 
                               ObservationSource         *obs_source)
    : SectionScheduler(gm_template, fp, obs_source)
  {}

  virtual ~ArchipelagosSectionScheduler() {}

};

#endif

