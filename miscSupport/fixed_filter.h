/*-
 * GMTK_FixedFilter.h
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


#ifndef FIXED_FILTER
#define FIXED_FILTER

#include "sArray.h"
#include "debug.h"

#include "adaptive_filter.h"

class FixedFilter : public AdaptiveFilter  {

protected:
public:

  double prevValue;
  double prevPrevValue;
  unsigned numSamplesLoaded;

  // constructor
  FixedFilter() : AdaptiveFilter(2) { prevValue = prevPrevValue = 0.0; numSamplesLoaded = 0; }
  
  virtual bool readyToMakePrediction()  { return (numSamplesLoaded >= 2); }

  virtual double makePrediction() {
    if (!readyToMakePrediction()) return 0;
    // return 2*prevPrevValue - prevValue;
    return 2*prevValue - prevPrevValue;
  }
  virtual void addNextSampleAndUpdate(double val) {
    prevPrevValue = prevValue;
    prevValue = val;
    numSamplesLoaded ++;
  }
  virtual void init() {
    numSamplesLoaded = 0;
  }

};

#endif




