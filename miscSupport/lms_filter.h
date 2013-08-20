/*-
 * GMTK_LMSFilter.h
 *
 *  Written by Jeff Bilmes <bilmes@ee.washington.edu>
 * 
 *  $Header$
 * 
 * Copyright (C) 2008 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 *
 */


#ifndef LMS_FILTER
#define LMS_FILTER

#include "sArray.h"
#include "debug.h"

#include "adaptive_filter.h"

class LMSFilter : public AdaptiveFilter  {

protected:
public:

  // we use the weights array as a circular buffer, so we
  // store the position of the most recently added sample
  unsigned mostRecentSamplePosition;

  // the x values, as a circular buffer. Most recent position is
  // mostRecentSamplePosition and least recent is
  // (mostRecentSamplePosition + 1) % order.
  sArray<double> history;
  
  double learning_rate;
  unsigned numSamplesLoaded;


  // constructor
  LMSFilter(unsigned _order, double _learning_rate);
  
  virtual bool readyToMakePrediction()  { return (numSamplesLoaded >= order) ; }

  virtual double makePrediction();
  virtual void addNextSampleAndUpdate(double val);
  virtual void init();

};

#endif




