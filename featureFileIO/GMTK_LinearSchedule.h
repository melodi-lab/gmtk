/*
 * GMTK_LinearSchedule.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (C) 2013 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

#ifndef GMTK_LINEARSCHEDULE_H
#define GMTK_LINEARSCHEDULE_H

#include <string.h>

#include "GMTK_TrainingSchedule.h"
#include "GMTK_ObservationSource.h"
#include "rand.h"
#include "error.h"

// Return non-overlapping training units of requested size in observation source order.
// If a segment's length is not a multiple of the unit size, the last unit from that
// segment will be short.

class LinearSchedule : public TrainingSchedule {

  unsigned curSegment, curFrame;
  unsigned maxFrame;         // last possible unit starting frame in current segment
  Range::iterator* trrng_it; // segment iterator

 public:

  LinearSchedule(unsigned feature_offset, unsigned features_per_frame,
		 unsigned label_offset,  unsigned label_domain_size,
		 bool one_hot, unsigned window_radius, unsigned unit_size, 
		 FileSource *obs_source, char const *trrng_str)
    : TrainingSchedule(feature_offset, features_per_frame, label_offset, 
		       label_domain_size, one_hot, window_radius, unit_size,
		       obs_source, trrng_str)
  {
    num_viable_units = 0;
    trrng_it = new Range::iterator(trrng->begin());
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
    trrng_it->reset();
    unsigned num_frames;
    do {
      curSegment = *(*trrng_it);
      if (!obs_source->openSegment(curSegment)) {
	error("ERROR: Unable to open observation file segment %u\n", curSegment);
      }
      num_frames = obs_source->numFrames();
      if (num_frames == 0) (*trrng_it)++;
    } while (num_frames == 0);
    maxFrame = ((num_frames-1) / unit_size) * unit_size;
    curFrame = 0;
  } 

  
  ~LinearSchedule() {
    if (trrng_it) delete trrng_it;
  }

  
  void nextTrainingUnit(unsigned &segment, unsigned &frame) { 
    if (curFrame > maxFrame) {
      unsigned num_frames;
      do {
	(*trrng_it)++;
	if (trrng_it->at_end()) {
	  trrng_it->reset();
	}
	curSegment = *(*trrng_it);
	if (!obs_source->openSegment(curSegment)) {
	  error("ERROR: Unable to open observation file segment %u\n", curSegment);
	}
	num_frames = obs_source->numFrames();
      } while (num_frames == 0); // must be a viable unit somewhere or the ctor would fail
      maxFrame = ((num_frames-1) / unit_size) * unit_size;
      curFrame = 0;
    } 
    segment = curSegment;
    frame = curFrame;
    curFrame += unit_size;
    TrainingSchedule::nextTrainingUnit(segment, frame);
  }
};

#endif
