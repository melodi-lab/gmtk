/*
 * GMTK_ShuffleSegmentSchedule.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (C) 2015 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

#ifndef GMTK_SHUFFLESEGMENTSCHEDULE_H
#define GMTK_SHUFFLESEGMENTSCHEDULE_H

#include "rand.h"
#include "GMTK_BaseSegmentSchedule.h"

// Knuth shuffle the batch order
class ShuffleSegmentSchedule : public BaseSegmentSchedule {
 private:
  vector<unsigned> batch_starting_segment; // first segment of each batch

 public:

  ShuffleSegmentSchedule(FileSource *obs_source, char const *trrng_str, unsigned batch_size) 
    : BaseSegmentSchedule(obs_source, trrng_str, batch_size)
  {
    batch_starting_segment.resize(num_training_units);
    Range::iterator s = trrng->begin();
    batch_starting_segment[0] = *s;
    for (unsigned i=1; i < num_training_units; ++i) {
      batch_starting_segment[i] = s.step_by(batch_size);
    }
    // Knuth shuffle
    for (unsigned i = num_training_units-1; i > 0; --i) {
      unsigned j = rnd.uniform(i);
      unsigned tmp = batch_starting_segment[i];
      batch_starting_segment[i] = batch_starting_segment[j];
      batch_starting_segment[j] = tmp;
    }
  }

  // Returns a vector containing the segment numbers for the next batch
  void getBatch(vector<unsigned> &batch) {
    generateBatch(batch_starting_segment[num_batches_dispensed], batch);
  }

};

#endif // GMTK_SHUFFLESEGMENTSCHEDULE
