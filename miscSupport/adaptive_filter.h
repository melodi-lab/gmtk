/*-
 * GMTK_AdaptiveFilter
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


#ifndef ADAPTIVE_FILTER
#define ADAPTIVE_FILTER

#include "sArray.h"
#include "debug.h"

class AdaptiveFilter  {

protected:
public:
  
  unsigned order;
  // the weights of the filter. The least recent position corresponds
  // to weights[0], and the most recent to weights[order-1].
  sArray <double> weights;

  AdaptiveFilter(unsigned _order) : order(_order) {
    if (order == 0) 
      error("Error: an adaptive filter can't have an order of %d.",order);
  }
  virtual ~AdaptiveFilter() {}

  virtual bool readyToMakePrediction() = 0;
  virtual double makePrediction() = 0;
  virtual void addNextSampleAndUpdate(double val) = 0;

};

#endif




