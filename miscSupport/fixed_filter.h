/*-
 * GMTK_FixedFilter.h
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
    return 2*prevPrevValue - prevValue;
  }
  virtual void addNextSampleAndUpdate(double val) {
    prevPrevValue = prevValue;
    prevValue = val;
    numSamplesLoaded ++;
  }

};

#endif




