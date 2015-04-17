/*
 * GMTK_SegmentSchedule.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (C) 2015 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

#ifndef GMTK_SEGMENTSCHEDULE_H
#define GMTK_SEGMENTSCHEDULE_H

#include <vector>

// The SegmentSchedule class abstracts various training schedules (order
// which the segments of training data are presented to a learning algorithm).
// SegmentSchedule classes just provide an order for segments, whereas 
// the TrainingSchedule classes schedule over both segments and frames.
// Subclasses implement specific schedules, for example:
//   RandomSegmentSchedule implements uniform random sampling with replacement
//   PermutationSegmentSchedule cycles through a random permuation (selected via cubic residue)
//   ShuffleSegmentSchedule cycles through a random permutation (selected via Knuth shuffle)
//   LinearSegmentSchedule just presents the training data in observation file order

// Note that SegmentSchedule classes assume that a training unit will consist of
// a contiguous sequence of of segments. I.e., if the training unit size is n,
// the training unit consists of segments [i,i+n-1]. SegmentSchedule also assumes
// training units do not overlap. Subclasses can override this last assumption as needed.

class SegmentSchedule {

 public:

  virtual ~SegmentSchedule() {}

  // Returns the number of viable training units in the input data.
  virtual unsigned numViableUnits() = 0;

  // Returns a vector containing the segment numbers for the next batch
  virtual void getBatch(vector<unsigned> &batch) = 0;
};

#endif // GMTK_SEGMENTSCHEDULE
