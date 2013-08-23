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

// The TrainingSchedule class abstracts various training schedules.
// Subclasses implement specific schedules, for example:
//   RandomSampleSchedule implements uniform random sampling with replacement.
//   PermutationSchedule cycles through a random permuation


// Used to represent the state of a TrainingSchedule object
class TrainingScheduleState {};


class TrainingSchedule {

 public:

  // Returns the current state of the TrainingSchedule object
  virtual TrainingScheduleState const *getState() = 0;


  // Sets the current state of the TrainingSchedule object
  virtual void reset(TrainingScheduleState const *state) = 0;


  // Returns the number of training units in the input data.

  // Some training schedules may sample the same training unit more than 
  // once (or not sample a training unit at all), so this is not necessarily 
  // the number of training units the client will process.
  virtual unsigned numTrainingUnits() = 0;


  // Describes the feature matrix returned by getFeatures()

  // Returns the number of features in a training instance in numFeatures.
  // Returns the number of instances in the training unit in numInstances.
  // Returns the number of floats between each instance in stride.
  // The features of a particular instance are contiguous in memory.
  virtual void describeFeatures(unsigned &numFeatures,
				unsigned &numInstances,
				unsigned &stride) = 0;


  // Returns true and sets segment and frame index to the first instance of the
  // next training unit, or returns false if the end of the schedule has been reached. 
  // Note that some training schedules never end, and so never return false.
  virtual bool nextTrainingUnit(unsigned &segment, unsigned &frame) = 0;


  // Returns the feature matrix of the training unit starting at semgnet, frame.

  // The data is organized as in ObservationSource::loadFrames() except
  // the returned pointer is to the first feature of the first instance (rather
  // than the first feature of the center frame of the first instance's window).

  // returns floatVecAtFrame(frame) - window_radius * stride + feature_offset
  // ensures data for frames [frame - window_radius, frame + unit_size + window_radius) are in memory
  virtual float *getFeatures(unsigned segment, unsigned frame) = 0;


  // Returns the label matrix corresponding to getFeatures(segment,frame).

  // The labels may be "1-hot" where all but one element of the domain
  // have value zero and the element of the domain corresponding to the
  // correct label has value 1, or there may be a posterior-like score 
  // for every element of the domain.
  virtual float *getLabels(unsigned segment, unsigned frame) = 0;


  // Describes the label matrix returned by getLabels()

  // Returns the size of the label's domain in domainSize.
  // Returns the number of instances in the training unit in numInstances.
  // Returns the number of floats between each instance in stride.
  virtual void describeLabels(unsigned &domainSize,
			      unsigned &numInstances,
			      unsigned &stride) = 0;
};

#endif
