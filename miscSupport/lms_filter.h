/*-
 * GMTK_LMSFilter.h
 *
 *  Written by Jeff Bilmes <bilmes@ee.washington.edu>
 * 
 *  $Header$
 * 
 * Copyright (c) 2008, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about
 * the suitability of this software for any purpose.  It is provided
 * "as is" without express or implied warranty.
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




