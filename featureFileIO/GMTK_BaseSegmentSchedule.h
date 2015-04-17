/*
 * GMTK_BaseSegmentSchedule.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (C) 2015 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

#ifndef GMTK_BASESEGMENTSCHEDULE_H
#define GMTK_BASESEGMENTSCHEDULE_H

#include "range.h"
#include "GMTK_SegmentSchedule.h"
#include "GMTK_FileSource.h"

// Partial implementation to serve as base class for SegmentSchedule subclasses
class BaseSegmentSchedule : public SegmentSchedule {
 protected:
  unsigned num_training_units;
  unsigned batch_size, short_batch_size;
  unsigned num_batches_dispensed;
  Range *trrng;

  // Returns a vector containing the segment numbers for the next batch
  // The batch consists of a continuous sequence of upto batch_size segments 
  // in the -trrng starting from first_segment

  // Note! this also updates num_batches_dispensed

  virtual void generateBatch(unsigned first_segment, vector<unsigned> &batch) {
    // last batch might be short...
    unsigned cur_batch_size = (num_training_units - 1 > num_batches_dispensed) ?
      batch_size : short_batch_size;
    batch.resize(cur_batch_size);

    Range::iterator s = trrng->begin();
    s.step_by(first_segment);
    for (unsigned i=0; i < cur_batch_size; ++i, ++s) {
      batch[i] = *s;
    }
    num_batches_dispensed = (num_batches_dispensed + 1) % num_training_units;
  }

 public:

  BaseSegmentSchedule(FileSource *obs_source, char const *trrng_str, unsigned batch_size) 
    : batch_size(batch_size), num_batches_dispensed(0)
  {
    assert(obs_source);
    trrng = new Range(trrng_str,0,obs_source->numSegments());
    if (trrng->length() <= 0) {
      error("Training range '%s' specifies empty set.\n", trrng_str);
    }
    unsigned num_segments = trrng->length();
    num_training_units = num_segments / batch_size;
    // last batch might be short...
    if (num_segments % batch_size > 0) {
      num_training_units +=  1;
      short_batch_size = num_segments % batch_size;
    } else {
      short_batch_size = batch_size;
    }
  }

  virtual ~BaseSegmentSchedule() {
    delete trrng;
  }

  // Returns the number of viable training units in the input data.
  unsigned numViableUnits() { return num_training_units; }

};

#endif // GMTK_BASESEGMENTSCHEDULE
