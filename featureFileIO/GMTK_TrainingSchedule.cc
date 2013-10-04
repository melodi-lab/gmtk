/*
 * GMTK_TrainingSchedule.cc
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (C) 2013 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

#include <string.h>

#include "GMTK_TrainingSchedule.h"

TrainingSchedule::TrainingSchedule(unsigned feature_offset, unsigned features_per_frame,
				   unsigned label_offset,  unsigned label_domain_size,
				   bool one_hot, unsigned window_radius, unsigned unit_size, 
				   FileSource *obs_source, char const *trrng_str)
   :  feature_offset(feature_offset), features_per_frame(features_per_frame),
      label_offset(label_offset), label_domain_size(label_domain_size),
      one_hot(one_hot), window_radius(window_radius),
      unit_size(unit_size), obs_source(obs_source), num_units_dispensed(0)
{
  assert(obs_source);
  if (window_radius > obs_source->startSkip()) {
    error("ERROR: -startSkip must be >= %u to support requested window radius\n", window_radius);
  }
  if (window_radius > obs_source->endSkip()) {
    error("ERROR: -endSkip must be >= %u to support requested window radius\n", window_radius);
  }
  // these should be error checked before calling ctor
  assert(feature_offset < obs_source->numContinuous());
  assert(feature_offset + features_per_frame <= obs_source->numContinuous());
  if (one_hot) {
    assert(obs_source->numContinuous() <= label_offset);
    assert(label_offset < obs_source->numFeatures());
  } else {
    assert(label_offset < obs_source->numContinuous());
    assert(label_offset + label_domain_size <= obs_source->numContinuous());
  }
  trrng = new Range(trrng_str,0,obs_source->numSegments());
  stride = obs_source->stride();
  total_frames = 0;
  num_viable_units = 0;
  num_segments = trrng->length();
  
  Range::iterator* trrng_it = new Range::iterator(trrng->begin());
  while (!trrng_it->at_end()) {
    unsigned i = (unsigned)(*(*trrng_it));
    if (!obs_source->openSegment(i)) {
      error("ERROR: Unable to open observation file segment %u\n", i);
    }
    total_frames += obs_source->numFrames();

    // Because startSkip & endSkip frames exist, all frames counted by
    // numFrames() are viable for unit size == 1 frame. Larger units
    // exclude the last unit_size - 1 frames.
    num_viable_units += obs_source->numFrames() - unit_size + 1;
    (*trrng_it)++;
  }
  delete trrng_it;
  unit_data = new float[unit_size * features_per_frame * (2 * window_radius + 1)];
  if (one_hot) {
    heated_labels = new float[unit_size * label_domain_size];
  } else {
    heated_labels = NULL;
  }
}


float *
TrainingSchedule::getLabels(unsigned segment, unsigned frame) {
  segment = trrng->index(segment);
  if (!obs_source->openSegment(segment)) {
    error("ERROR: Unable to open observation file segment %u\n", segment);
  } 
  (void) obs_source->loadFrames(frame, unit_size); // ensure necessary data is in memory
  if (one_hot) {
    memset(heated_labels, 0, label_domain_size * unit_size * sizeof(float));
    unsigned *labels = obs_source->unsignedVecAtFrame(frame) + 
      label_offset - obs_source->numContinuous();
    float *dest = heated_labels;
    for (unsigned i=0; i < unit_size; i+=1, labels += stride, dest += label_domain_size) {
      unsigned label = *labels;
      if (label >= label_domain_size) {
	error("ERROR: observed label %u at frame %u in segment %u is too large for label domain size %u\n", label, frame, segment, label_domain_size);
      }
      dest[label] = 1.0;
    }
    return heated_labels;
  } else {
    return obs_source->floatVecAtFrame(frame) + label_offset;
  }
}
