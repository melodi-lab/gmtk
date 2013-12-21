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

// Random sample with replacement. Note that sampled training units may overlap.

class RandomSampleSchedule : public TrainingSchedule {

  float *segment_dist;    // segment_dist[i] = probability of sampling unit from segment i

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
    unsigned segment_count[num_segments];  // # of viable units in each segment
    memset((void *)segment_count, 0, num_segments * sizeof(unsigned));
    for (unsigned i=0; i < num_segments; i+=1) {
      if (!obs_source->openSegment(trrng->index(i))) {
	error("ERROR: Unable to open observation file segment %u\n", trrng->index(i));
      }
      unsigned num_frames = obs_source->numFrames();
      if (num_frames >= unit_size) {
	segment_count[i] = num_frames - unit_size + 1;
      } else if (num_frames > 0) {
	segment_count[i] = 1;
      } else {
	segment_count[i] = 0;
      }
    }
    for (unsigned i=0; i < num_segments; i+=1) {
      segment_dist[i] = (double)segment_count[i] / (double)num_segments;
    }
  } 
  

  ~RandomSampleSchedule() {
    if (segment_dist) delete[] segment_dist;
  }


  void nextTrainingUnit(unsigned &segment, unsigned &frame) { 
    segment = (unsigned) rnd.sample(num_segments, segment_dist);
    if (!obs_source->openSegment(trrng->index(segment))) {
	error("ERROR: Unable to open observation file segment %u\n", segment);
    }
    unsigned num_frames = obs_source->numFrames();
    unsigned length;
    if (num_frames >= unit_size) {
      frame = (unsigned) rnd.uniformOpen( num_frames - unit_size + 1 );
    } else if (num_frames > 0) {
      frame = 0; // it's going to be short
    } else {
      assert(0); // this segment should have had 0 probability to be chosen!
    }
    TrainingSchedule::nextTrainingUnit(segment, frame);
  }
};

#endif
