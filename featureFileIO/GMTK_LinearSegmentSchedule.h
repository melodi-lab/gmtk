/*
 * GMTK_LinearSegmentSchedule.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (C) 2015 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

#ifndef GMTK_LINEARSEGMENTSCHEDULE_H
#define GMTK_LINEARSEGMENTSCHEDULE_H

#include "GMTK_BaseSegmentSchedule.h"

// Just present the batches in order
class LinearSegmentSchedule : public BaseSegmentSchedule {
 public:

  LinearSegmentSchedule(FileSource *obs_source, char const *trrng_str, unsigned batch_size) 
    : BaseSegmentSchedule(obs_source, trrng_str, batch_size)
  { }

  // Returns a vector containing the segment numbers for the next batch
  void getBatch(vector<unsigned> &batch) {
    generateBatch(num_batches_dispensed * batch_size, batch);
  }

};

#endif // GMTK_LINEARSEGMENTSCHEDULE
