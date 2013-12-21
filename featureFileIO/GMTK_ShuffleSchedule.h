/*
 * GMTK_ShuffleSchedule.h
 *   Return training units according to a random permutation.
 *
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (C) 2013 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

#ifndef GMTK_SHUFFLESCHEDULE_H
#define GMTK_SHUFFLESCHEDULE_H

#include <string.h>

#include "GMTK_TrainingSchedule.h"
#include "GMTK_ObservationSource.h"
#include "rand.h"
#include "error.h"

// Return non-overlapping training units of requested size in an order determined
// by random sampling without replacement (Knuth shuffle random permutation).

class ShuffleSchedule : public TrainingSchedule {

  unsigned *segmentPerm;     // segment # of training unit permutation
  unsigned *framePerm;       // frame # of training unit permutation
  unsigned  idx;             // position in permutation

 public:

  ShuffleSchedule(unsigned feature_offset, unsigned features_per_frame,
		 unsigned label_offset,  unsigned label_domain_size,
		 bool one_hot, unsigned window_radius, unsigned unit_size, 
		 FileSource *obs_source, char const *trrng_str)
    : TrainingSchedule(feature_offset, features_per_frame, label_offset, 
		       label_domain_size, one_hot, window_radius, unit_size,
		       obs_source, trrng_str), idx(0)
  {
    // count viable training units
    num_viable_units = 0;
    Range::iterator* trrng_it = new Range::iterator(trrng->begin());
    while (!trrng_it->at_end()) {
      unsigned i = (unsigned)(*(*trrng_it));
      if (!obs_source->openSegment(i)) {
	error("ERROR: Unable to open observation file segment %u\n", i);
      }
      unsigned num_frames = obs_source->numFrames();
      unsigned units = num_frames / unit_size;
      if (num_frames % unit_size) units += 1;
      if (units == 0) {
	warning("WARNING: segment %u contains no frames\n", i);
      }
      num_viable_units += units;
      (*trrng_it)++;
    }
    if (num_viable_units == 0) {
      error("ERROR: observation files contain no viable training instances\n");
    }

    // initialize segment & frame permutation arrays 
    segmentPerm = new unsigned[num_viable_units];
    framePerm = new unsigned[num_viable_units];
    trrng_it->reset();
    unsigned j=0;
    while (!trrng_it->at_end()) {
      unsigned i = (unsigned)(*(*trrng_it));
      if (!obs_source->openSegment(i)) {
	error("ERROR: Unable to open observation file segment %u\n", i);
      }
      unsigned num_frames = obs_source->numFrames();
      unsigned units = num_frames / unit_size;
      if (num_frames % unit_size) units += 1;
      for (unsigned u = 0, f=0; u < units; u+=1, f+=unit_size, j+=1) {
	// Knuth shuffle - http://en.wikipedia.org/wiki/Random_permutation
	unsigned k = rnd.uniform(j);
	segmentPerm[j] = segmentPerm[k];
	framePerm[j] = framePerm[k];
	segmentPerm[k] = i;
	framePerm[k] = f;
      }
      (*trrng_it)++;
    }
  } 

  
  ~ShuffleSchedule() {
    delete[] segmentPerm;
    delete[] framePerm;
  }

  
  void nextTrainingUnit(unsigned &segment, unsigned &frame) { 
    segment = segmentPerm[idx];
    frame = framePerm[idx];
    idx = (idx + 1) % num_viable_units;
    TrainingSchedule::nextTrainingUnit(segment, frame);
  }
};

#endif
