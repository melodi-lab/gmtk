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

  unsigned  num_units;
  Range    *trrng;        // observation file segment range
  unsigned  num_segments; // in the observation file segment range

  float    *heated_labels; // hold label matrix if one_hot is true

  unsigned num_units_dispensed; // # of calls to nextTrainingUnit

 public:

  TrainingSchedule(unsigned feature_offset, unsigned num_features,
		   unsigned label_offset,  unsigned label_domain_size,
		   bool one_hot, unsigned window_radius, unsigned unit_size, 
		   FileSource *obs_source, char const *trrng_str);

  ~TrainingSchedule() {
    if (trrng) delete trrng;
    if (heated_labels) delete[] heated_labels;
  }


  // Returns the current state of the TrainingSchedule object
  virtual TrainingScheduleState const *getState() = 0;


  // Sets the current state of the TrainingSchedule object & resets
  // num_units_dispensed to zero.
  virtual void reset(TrainingScheduleState const *state) {
    num_units_dispensed = 0;
  }

  // Returns the number of training units in the input data.

  // Some training schedules may sample the same training unit more than 
  // once (or not sample a training unit at all), so this is not necessarily 
  // the number of training units the client will process.
  virtual unsigned numTrainingUnits() { return num_units; }


  // Describes the feature matrix returned by getFeatures()

  // Returns the number of features in a training instance in numFeatures.
  // Returns the number of instances in the training unit in numInstances.
  // Returns the number of floats between each instance in stride.
  // The features of a particular instance are contiguous in memory.
  virtual void describeFeatures(unsigned &numFeatures,
				unsigned &numInstances,
				unsigned &stride)
  {
    numFeatures = num_features; numInstances = unit_size; stride = this->stride;
  }


  // Returns true and sets segment and frame index to the first instance of the
  // next training unit, or returns false if the end of the schedule has been reached. 
  // Note that some training schedules never end, and so never return false.
  virtual bool nextTrainingUnit(unsigned &segment, unsigned &frame) {
    num_units_dispensed += 1;
    return false;
  }


  // Returns the feature matrix of the training unit starting at semgnet, frame.

  // The data is organized as in ObservationSource::loadFrames() except
  // the returned pointer is to the first feature of the first instance (rather
  // than the first feature of the center frame of the first instance's window).

  // returns floatVecAtFrame(frame) - window_radius * stride + feature_offset
  // ensures data for frames [frame - window_radius, frame + unit_size + window_radius) are in memory
  virtual float *getFeatures(unsigned segment, unsigned frame) {
    segment = trrng->index(segment);
    if (!obs_source->openSegment(segment)) {
	error("ERROR: Unable to open observation file segment %u\n", segment);
    } 
    (void) obs_source->loadFrames(frame, unit_size); // ensure necessary data is in memory
    return obs_source->floatVecAtFrame(frame) - window_radius * stride + feature_offset;
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
