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

// The TrainingSchedule class abstracts various training schedules.
// Subclasses implement specific schedules, for example:
//   RandomSampleSchedule implements uniform random sampling with replacement.
//   PermutationSchedule cycles through a random permuation


// Used to represent the state of a TrainingSchedule object
class TrainingScheduleState {};


class TrainingSchedule {
 protected:
  unsigned feature_offset;       // offset of first input feature 
  unsigned features_per_frame;   // # features to take from each frame in the window

  unsigned label_offset;         // offset of first label feature 
                                 //   this will be >= # continuous features if 1-hot
  unsigned label_domain_size;
  bool     one_hot;              // true iff we must expand single int -> vector

  unsigned stride;               // obs file stride
  unsigned window_radius;
  unsigned unit_size;            // # of instances in a training unit (minibatch size)

  // Can't use StreamSource because it doesn't know how many segments/frames
  // there are until the end of the input is reached.
  FileSource *obs_source;

  unsigned  num_viable_units;    // total # of frames across all segments that
                                 //   can be the start of a training unit

  unsigned  total_frames;        // total frames between startSkip & endSkip 
                                 //   across all segments
 
  Range    *trrng;               // observation file segment range
  unsigned  num_segments;        // in the observation file segment range

  float    *unit_data;           // hold the training data for a unit
  float    *heated_labels;       // hold label matrix if one_hot is true

  unsigned num_units_dispensed;  // # of calls to nextTrainingUnit since last reset()

 public:

  TrainingSchedule(unsigned feature_offset, unsigned features_per_frame,
		   unsigned label_offset,  unsigned label_domain_size,
		   bool one_hot, unsigned window_radius, unsigned unit_size, 
		   FileSource *obs_source, char const *trrng_str);

  ~TrainingSchedule() {
    if (trrng) delete trrng;
    if (unit_data) delete[] unit_data;
    if (heated_labels) delete[] heated_labels;
  }


  // Returns the current state of the TrainingSchedule object
  virtual TrainingScheduleState const *getState() = 0;


  // Sets the current state of the TrainingSchedule object & resets
  // num_units_dispensed to zero.
  virtual void reset(TrainingScheduleState const *state) {
    num_units_dispensed = 0;
  }

  // Returns the number of viable training units in the input data.
  virtual unsigned numViableUnits() { return num_viable_units; }


  // Some training schedules may sample the same training unit more than 
  // once (or not sample a training unit at all), so this is not necessarily 
  // the number of training units the client will process.
  virtual unsigned numTrainingUnits() { return total_frames / unit_size; }

  // Returns true and sets segment and frame index to the first instance of the
  // next training unit, or returns false if the end of the schedule has been reached. 
  // Note that some training schedules never end, and so never return false.
  virtual bool nextTrainingUnit(unsigned &segment, unsigned &frame) {
    num_units_dispensed += 1;
    return false;
  }


  // Returns the observation source matrix of the training unit starting at segment, frame.

  // The data is organized as in ObservationSource::loadFrames() except
  // the returned pointer is to the first feature of the first instance (rather
  // than the first feature of the center frame of the first instance's window).

  // returns floatVecAtFrame(frame) - window_radius * stride + feature_offset
  // ensures data for frames [frame - window_radius, frame + unit_size + window_radius) are in memory
  virtual float *getObservations(unsigned segment, unsigned frame) {
    segment = trrng->index(segment);
    if (!obs_source->openSegment(segment)) {
	error("ERROR: Unable to open observation file segment %u\n", segment);
    } 
    (void) obs_source->loadFrames(frame, unit_size); // ensure necessary data is in memory
    return obs_source->floatVecAtFrame(frame) - window_radius * stride + feature_offset;
  }


  // Describes the observation source matrix returned by getObservations()

  // Returns the number of features per frame to use in numFeatures
  // Returns the number of frames in a training instance in numFrames
  // Returns the number of instances in the training unit in numInstances.
  // Returns the number of floats between each instance in stride.
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


  // Returns the feature matrix of the training unit starting at segment, frame.

  virtual float *getFeatures(unsigned segment, unsigned frame) {
    unsigned frames_per_instance = 2 * window_radius + 1;
    float *obsData = getObservations(segment, frame);
    float *dest = unit_data, *src;
    for (unsigned b=0; b < unit_size; b+=1) { // loop over instances in batch
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
  // Returns the number of instances in the training unit in numInstances.

  // The features of a particular instance are contiguous in memory.
  virtual void describeFeatures(unsigned &numFeatures,
				unsigned &numInstances)
  {
    numFeatures = features_per_frame * (2 * window_radius + 1); 
    numInstances = unit_size; 
  }


  // Returns the label matrix corresponding to getFeatures(segment,frame).

  // The labels may be "1-hot" where all but one element of the domain
  // have value zero and the element of the domain corresponding to the
  // correct label has value 1, or there may be a posterior-like score 
  // for every element of the domain.
  virtual float *getLabels(unsigned segment, unsigned frame);


  // Describes the label matrix returned by getLabels()

  // Returns the size of the label's domain in domainSize.
  // Returns the number of instances in the training unit in numInstances.
  // Returns the number of floats between each instance in stride.
  virtual void describeLabels(unsigned &domainSize,
			      unsigned &numInstances,
			      unsigned &stride)
  {
    domainSize = label_domain_size; numInstances = unit_size; 
    if (one_hot) {
      stride = label_domain_size;
    } else {
      stride = this->stride;
    }
  }
  
};

#endif
