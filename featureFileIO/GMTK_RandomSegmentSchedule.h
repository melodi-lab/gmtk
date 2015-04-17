/*
 * GMTK_RandomSegmentSchedule.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (C) 2015 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

#ifndef GMTK_RANDOMSEGMENTSCHEDULE_H
#define GMTK_RANDOMSEGMENTSCHEDULE_H

#include "rand.h"
#include "GMTK_BaseSegmentSchedule.h"

// Random sample (with replacement) batch_size segments to
// form a batch. Note! this schedule does not produce batches
// of sequential segments - you get batch_size segments chosen
// completely at random in every batch.

class RandomSegmentSchedule : public BaseSegmentSchedule {
 private:
  unsigned num_segments;

 public:

  RandomSegmentSchedule(FileSource *obs_source, char const *trrng_str, unsigned batch_size) 
    : BaseSegmentSchedule(obs_source, trrng_str, batch_size), num_segments(trrng->length())
  {}

  // Returns a vector containing the segment numbers for the next batch
  void getBatch(vector<unsigned> &batch) {
    batch.resize(batch_size);
    for (unsigned i=0; i < batch_size; ++i) {
      batch[i] = trrng->index(rnd.uniformOpen(num_segments));
    }
  }

};

#endif // GMTK_RANDOMSEGMENTSCHEDULE
