/*
 * GMTK_RandomSampleSchedule.cc
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
#include "GMTK_ObservationSource.h"
#include "rand.h"
#include "error.h"


class RandomSampleSchedule : public TrainingSchedule {
  
  unsigned feature_offset;
  unsigned num_features; 

  unsigned label_offset;
  unsigned label_domain_size;
  bool     one_hot;           // true iff we must expand single int -> vector

  unsigned stride;
  unsigned window_radius;
  unsigned unit_size;

  // Can't use StreamSource because it doesn't know how many segments/frames
  // there are until the end of the input is reached.
  FileSource *obs_source;

  unsigned num_units;
  unsigned num_segments; // in the observation file

  float    *heated_labels; // hold label matrix if one_hot is true

  float    *segment_dist;  // segment_dist[i] = probability of sampling unit from segment i
                    
 public:

  RandomSampleSchedule(unsigned feature_offset, unsigned num_features,
		       unsigned label_offset,  unsigned label_domain_size,
		       bool one_hot, unsigned window_radius,
		       unsigned unit_size, FileSource *obs_source)
    : feature_offset(feature_offset), num_features(num_features),
      label_offset(label_offset), label_domain_size(label_domain_size),
      one_hot(one_hot), window_radius(window_radius),
      unit_size(unit_size), obs_source(obs_source)
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
    assert(feature_offset + num_features <= obs_source->numContinuous());
    if (one_hot) {
      assert(obs_source->numContinuous() < label_offset);
      assert(label_offset + label_domain_size <= obs_source->numFeatures());
    } else {
      assert(label_offset < obs_source->numContinuous());
      assert(label_offset + label_domain_size <= obs_source->numContinuous());
    }

    stride = obs_source->stride();
    num_units = 0;
    num_segments = obs_source->numSegments();
    segment_dist = new float[num_segments];

    for (unsigned i=0; i < num_segments; i+=1) {
      if (!obs_source->openSegment(i)) {
	error("ERROR: Unable to open observation file segment %u\n", i);
      }
      
      // Because startSkip & endSkip frames exist, all frames counted by
      // numFrames() are viable for unit size == 1 frame. Larger units
      // exclude the last unit_size - 1 frames.
      num_units += obs_source->numFrames() - unit_size + 1;
      segment_dist[i] = (float)(obs_source->numFrames() - unit_size + 1);
    }
    for (unsigned i=0; i < num_segments; i+=1) {
      segment_dist[i] /= num_units;
    }
    if (one_hot) {
      heated_labels = new float[unit_size * label_domain_size];
    } else {
      heated_labels = NULL;
    }
  } 


  ~RandomSampleSchedule() {
    if (segment_dist) delete[] segment_dist;
    if (heated_labels) delete[] heated_labels;
  }


  // RAND's state is global, and other GMTK functions consume random variables,
  // so RandomSampleSchedule would either need a private PRNG implementation or
  // else resetting it would stomp on all the other RAND clients. Moreover, since
  // the RAND output stream is used by other clients in unpredictable ways, we 
  // couldn't guarantee the same sample sequence even if we did reset the seed.
  TrainingScheduleState const *getState() {
    return new TrainingScheduleState();
  }

  void reset(TrainingScheduleState const *state) {
  }


  unsigned numTrainingUnits() { return num_units; }


  void describeFeatures(unsigned &numFeatures, unsigned &numInstances, unsigned &stride) {
    numFeatures = num_features; numInstances = unit_size; stride = this->stride;
  }


  bool nextTrainingUnit(unsigned &segment, unsigned &frame) { 
    segment = (unsigned) rnd.sample(num_segments, segment_dist);
    if (!obs_source->openSegment(segment)) {
	error("ERROR: Unable to open observation file segment %u\n", segment);
    } 
    frame = (unsigned) rnd.uniformOpen( obs_source->numFrames() - unit_size + 1 );
    return true;
  }


  float *getFeatures(unsigned segment, unsigned frame) {
    if (!obs_source->openSegment(segment)) {
	error("ERROR: Unable to open observation file segment %u\n", segment);
    } 
    (void) obs_source->loadFrames(frame, unit_size); // ensure necessary data is in memory
    return obs_source->floatVecAtFrame(frame) - window_radius * stride + feature_offset;
  }
    

  float *getLabels(unsigned segment, unsigned frame) {
    if (!obs_source->openSegment(segment)) {
	error("ERROR: Unable to open observation file segment %u\n", segment);
    } 
    (void) obs_source->loadFrames(frame, unit_size); // ensure necessary data is in memory
    if (one_hot) {
      memset(heated_labels, 0, label_domain_size * unit_size * sizeof(float));
      unsigned *labels = obs_source->unsignedVecAtFrame(frame) + label_offset;
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


  void describeLabels(unsigned &domainSize, unsigned &numInstances, unsigned &stride) {
    domainSize = label_domain_size; numInstances = unit_size; 
    if (one_hot) {
      stride = label_domain_size;
    } else {
      stride = this->stride;
    }
  }

};

