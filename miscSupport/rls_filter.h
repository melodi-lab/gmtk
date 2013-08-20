/*-
 * GMTK_RLSFilter.h
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


#ifndef RLS_FILTER
#define RLS_FILTER

#include "sArray.h"
#include "debug.h"

#include "adaptive_filter.h"

class RLSFilter : public AdaptiveFilter  {
  
  // recursive least squares filter.

protected:
public:


  sArray<double> tmpVec;

  // a order x order concentration matrix.
  sArray<double> concMat;


  // we use the weights array as a circular buffer, so we
  // store the position of the most recently added sample
  unsigned mostRecentSamplePosition;

  // the x values, as a circular buffer. Most recent position is
  // mostRecentSamplePosition and least recent is
  // (mostRecentSamplePosition + 1) % order.
  sArray<double> history;

  // forgetting coefficient = fc, should be set such
  // that 0 < fc <= 1. At 1 is the least amount of forgetting,
  // while close to 0 means that it forgets very rapidly.
  double forgetting_coef; 
  unsigned numSamplesLoaded;


  // constructor
  RLSFilter(unsigned _order, double _forgetting_coef);
  
  bool readyToMakePrediction()  { return (numSamplesLoaded >= order) ; }
  double makePrediction();
  void addNextSampleAndUpdate(double val);
  void init();

};

#endif




