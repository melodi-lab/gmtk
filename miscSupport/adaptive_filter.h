/*-
 * GMTK_AdaptiveFilter
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


#ifndef ADAPTIVE_FILTER
#define ADAPTIVE_FILTER

#include "sArray.h"
#include "debug.h"

class AdaptiveFilter  {

protected:
public:
  
  unsigned order;
  // The weights of the filter. The least recent position corresponds
  // to weights[0], and the most recent to weights[order-1].
  sArray <double> weights;

  AdaptiveFilter(unsigned _order) : order(_order) {
    if (order == 0) 
      error("Error: an adaptive filter can't have an order of %d.",order);
  }
  virtual ~AdaptiveFilter() {}

  virtual bool readyToMakePrediction() { return 1; }
  virtual double makePrediction() { return 0; }
  virtual void addNextSampleAndUpdate(double val) {}
  virtual void init() { fprintf(stderr,"abstract function\n"); }

};

#endif




