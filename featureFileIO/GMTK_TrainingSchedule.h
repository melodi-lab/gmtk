/*
 * GMTK_TrainingSchedule.h
 * 
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (C) 2013 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

#ifndef GMTK_TRAININGSCHEDULE_H
#define GMTK_TRAININGSCHEDULE_H

#include "GMTK_FileSource.h"
#include "range.h"
#include "error.h"

// The TrainingSchedule class abstracts various training schedules (order
// which the training examples are presented to a learning algorithm).
// TrainingSchedule subclasses schedule over both segments and frames, 
// whereas the SegmentSchedule subclasses schedule only over segments.
// Subclasses implement specific schedules, for example:
//   RandomSampleSchedule implements uniform random sampling with replacement
//   PermutationSchedule cycles through a random permuation (selected via cubic residue)
//   ShuffleSchedule cycles through a random permutation (selected via Knuth shuffle)
//   LinearSchedule just presents the training data in observation file order

// Note that TrainingSchedule assumes that a training unit will consist of
// a contiguous sequence of of intances. I.e., if the training unit size is n,
// the n instances come from windows centered on frames [i,i+n-1]. Subclasses
// can override more of the methods if they want non-contiguous units. TrainingSchedule
// also assumes training units may overlap. Again, subclasses can override as needed.


class TrainingSchedule {
 protected:
  unsigned feature_offset;       // offset of first input feature 
  unsigned features_per_frame;   // # features to take from each frame in the window

  unsigned label_offset;         // offset of first label feature 
                                 //   this will be >= # continuous features if 1-hot
  unsigned label_domain_size;
  bool     one_hot;              // true iff we must expand single int -> vector

  unsigned stride;               // obs file frame stride
  unsigned window_radius;        // 2r+1 frames constitute a single training instance
  unsigned unit_size;            // # of instances in a training unit (minibatch size)

  // Can't use StreamSource because it doesn't know how many segments/frames
  // there are until the end of the input is reached.
  FileSource *obs_source;

  unsigned  num_viable_units;    // total # of frames across all segments that
                                 //   can be the start of a training unit (units
                                 //   may overlap, so this may be > total_frames / unit_size)

  unsigned  total_frames;        // total frames between startSkip & endSkip 
                                 //   across all segments
 
  Range    *trrng;               // observation file segment range
  unsigned  num_segments;        // in the observation file segment range

  float    *unit_data;           // hold the training data for a unit
  float    *heated_labels;       // hold label matrix if one_hot is true

  unsigned num_units_dispensed;  // # of calls to nextTrainingUnit 



  // Returns the observation source matrix of the training unit starting at segment, frame.
  // length is set to the number of instances in the training unit (some subclasses may 
  // return short units)

  // The data is organized as in ObservationSource::loadFrames() except
  // the returned pointer is to the first feature of the first instance (rather
  // than the 0th feature of the center frame of the first instance's window).

  // Returns floatVecAtFrame(frame) - window_radius * stride + feature_offset
  // Ensures data for frames [frame - window_radius, frame + length + window_radius) are in memory
  virtual float *getObservations(unsigned segment, unsigned frame, unsigned &length) {
    segment = trrng->index(segment);
    if (!obs_source->openSegment(segment)) {
	error("ERROR: Unable to open observation file segment %u\n", segment);
    }
    unsigned num_frames = obs_source->numFrames();
    length = (num_frames > frame + unit_size - 1) ? unit_size : num_frames - frame;
    (void) obs_source->loadFrames(frame, length); // ensure necessary data is in memory
    return obs_source->floatVecAtFrame(frame) - window_radius * stride + feature_offset;
  }


  // Describes the observation source matrix returned by getObservations()

  // Returns the number of features per frame to use in numFeatures
  // Returns the number of frames in a single training instance in numFrames (window diameter)
  // Returns the number of instances in the training unit in numInstances.
  // Returns the number of floats between each frame in stride (note this is not the stride between training instances).
  virtual void describeObservations(unsigned &numFeatures,
				    unsigned &numFrames,
				    unsigned &numInstances,
				    unsigned &stride)
  {
    numFeatures = features_per_frame;
    numFrames = 2 * window_radius + 1; 
    numInstances = unit_size; 
    stride = this->stride;
  }

 public:

  TrainingSchedule(unsigned feature_offset, unsigned features_per_frame,
		   unsigned label_offset,  unsigned label_domain_size,
		   bool one_hot, unsigned window_radius, unsigned unit_size, 
		   FileSource *obs_source, char const *trrng_str);

  virtual ~TrainingSchedule() {
    if (trrng) delete trrng;
    if (unit_data) delete[] unit_data;
    if (heated_labels) delete[] heated_labels;
  }


  // Returns the number of training instances in the input data.
  virtual unsigned numInstances() { return total_frames; }


  // Returns the number of viable training units in the input data.
  virtual unsigned numViableUnits() { return num_viable_units; }


  // Some training schedules may sample the same training unit more than 
  // once (or not sample a training unit at all), so this is not necessarily 
  // the number of training units the client will process.
  virtual unsigned numTrainingUnitsPerEpoch() { return total_frames / unit_size +
                                                  ( (total_frames % unit_size > 0) ? 1 : 0 ); }

  // Sets segment and frame index to the first instance of the next training unit.
  // frame = center frame of first instance
  virtual void nextTrainingUnit(unsigned &segment, unsigned &frame) {
    num_units_dispensed += 1;
  }


  // Returns the feature matrix of the training unit starting at segment, frame.
  // length is set to the number of training instances it contains.
  // Note that some training schedules may return units shorter than unit_size.

  virtual float *getFeatures(unsigned segment, unsigned frame, unsigned &length) {
    unsigned frames_per_instance = 2 * window_radius + 1;
    float *obsData = getObservations(segment, frame, length);
    float *dest = unit_data, *src;
    for (unsigned b=0; b < length; b+=1) {    // loop over instances in unit
      src = obsData + b * stride;             //   first frame of current instance
      for (unsigned w=0; w < frames_per_instance; w+=1) { // loop over frames in instance
	memcpy((void *)dest, (void const *)src, features_per_frame * sizeof(float));
	dest += features_per_frame;
	src  += stride;
      }
    }
    return unit_data;
  }


  // Describes the feature matrix returned by getFeatures()

  // Returns the number of features in a training instance in numFeatures.
  // Returns the maximum number of instances in the training unit in maxInstances.

  // The features of a particular instance are contiguous in memory.
  virtual void describeFeatures(unsigned &numFeatures,
				unsigned &maxInstances)
  {
    numFeatures = features_per_frame * (2 * window_radius + 1); 
    maxInstances = unit_size; 
  }


  // Returns the label matrix corresponding to getFeatures(segment,frame).

  // The labels may be "1-hot" where all but one element of the domain
  // have value zero and the element of the domain corresponding to the
  // correct label has value 1, or there may be a posterior-like score 
  // for every element of the domain.
  // length is set to the number of training instances it contains.
  // Note that some training schedules may return units shorter than unit_size.
  virtual float *getLabels(unsigned segment, unsigned frame, unsigned &length);


  // Describes the label matrix returned by getLabels()

  // Returns the size of the label's domain in domainSize.
  // Returns the maximum number of instances in the training unit in maxInstances.
  // Returns the number of floats between each instance in stride.
  virtual void describeLabels(unsigned &domainSize,
			      unsigned &maxInstances,
			      unsigned &stride)
  {
    domainSize = label_domain_size; maxInstances = unit_size; 
    if (one_hot) {
      stride = label_domain_size;
    } else {
      stride = this->stride;
    }
  }

};

#endif
