/*
 * GMTK_RandomSampleSchedule.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (C) 2013 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

#ifndef GMTK_RANDOMSAMPLESCHEDULE_H
#define GMTK_RANDOMSAMPLESCHEDULE_H

#include <string.h>

#include "GMTK_TrainingSchedule.h"
#include "GMTK_ObservationSource.h"
#include "rand.h"
#include "error.h"


class RandomSampleSchedule : public TrainingSchedule {

  float *segment_dist;  // segment_dist[i] = probability of sampling unit from segment i
               
 public:

  RandomSampleSchedule(unsigned feature_offset, unsigned features_per_frame,
		       unsigned label_offset,  unsigned label_domain_size,
		       bool one_hot, unsigned window_radius, unsigned unit_size, 
		       FileSource *obs_source, char const *trrng_str)
    : TrainingSchedule(feature_offset, features_per_frame, label_offset, 
		       label_domain_size, one_hot, window_radius, unit_size,
		       obs_source, trrng_str)
  {
    segment_dist = new float[num_segments];
    unsigned segment_count[num_segments];
    memset((void *)segment_count, 0, num_segments * sizeof(unsigned));
    for (unsigned i=0; i < num_segments; i+=1) {
      if (!obs_source->openSegment(trrng->index(i))) {
	error("ERROR: Unable to open observation file segment %u\n", trrng->index(i));
      }
      segment_count[i] = obs_source->numFrames() - unit_size + 1;
    }
    for (unsigned i=0; i < num_segments; i+=1) {
      segment_dist[i] = (double)segment_count[i] / (double)num_segments;
    }
  } 
  

  ~RandomSampleSchedule() {
    if (segment_dist) delete[] segment_dist;
  }


  // RAND's state is global, and other GMTK functions consume random variables,
  // so RandomSampleSchedule would either need a private PRNG implementation or
  // else resetting it would stomp on all the other RAND clients. Moreover, since
  // the RAND output stream is used by other clients in unpredictable ways, we 
  // couldn't guarantee the same sample sequence even if we did reset the seed.
  TrainingScheduleState const *getState() {
    return NULL;
  }


  void reset(TrainingScheduleState const *state) {
    TrainingSchedule::reset(state);
  }


  bool nextTrainingUnit(unsigned &segment, unsigned &frame) { 
    if (num_units_dispensed >= num_viable_units) return false;

    segment = (unsigned) rnd.sample(num_segments, segment_dist);
    if (!obs_source->openSegment(trrng->index(segment))) {
	error("ERROR: Unable to open observation file segment %u\n", segment);
    } 
    frame = (unsigned) rnd.uniformOpen( obs_source->numFrames() - unit_size + 1 );
    (void) TrainingSchedule::nextTrainingUnit(segment, frame);
    return true;
  }
};

#endif
