/*
 * GMTK_SpaceManager.h
 *   GMTK Space Manager, manages space size allocation requests
 *   by various classes. I.e., It keeps track of how much has been
 *   allocated, what the current average allocation is, how do
 *   increase to the next size (e.g., multiply by 2 or something else), etc.
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2003, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about the suitability of this software
 * for any purpose. It is provided "as is" without express or implied warranty.
 *
 *
 * $Header$
 *
 */

#ifndef GMTK_SPACEMANAGER_H
#define GMTK_SPACEMANAGER_H

#include <vector>
#include <string>
#include <map>

#include <stdio.h>
#include <stdlib.h>



#include "debug.h"

class SpaceManager : public IM {

  // When another allocation request comes in, how much we currently
  // should ask to allocate.
  unsigned currentAllocationSize;

  // The initial size of an allocation. Should be >= 1.
  // should be const (but compiler complains)
  unsigned startingSize;

  // The smallest size used by anyone using this space manager
  // but that is >= startingSize.
  unsigned minAllocationSize;

  // Rate of exponential growth for memory allocation
  // should be const (but compiler complains)
  float growthFactor;

  // Rate of additive growth for memory allocation
  // should be const (but compiler complains)
  unsigned growthAddition;

  // Between instances of an object requesting via the same
  // spacemanager, how quickly should we decay. 
  // Must have 0 < decayRate
  // should be const (but compiler complains)
  float decayRate;

  // statistics on memory allocation.
  double sumRequested;
  unsigned numberRequests;

public:


  SpaceManager(const unsigned _startingSize=2,
	       const float _growthFactor=2.0,
	       const unsigned _growthAddition=1,
	       const float _decayRate=1.0)
    : startingSize(_startingSize),
      growthFactor(_growthFactor),
      growthAddition(_growthAddition),
      decayRate(_decayRate)
  {
    if (decayRate <= 0.0)
      error("ERROR: SpaceManager, decay rate %f must be > 0", decayRate);
    if (growthFactor  < 1.0)
      error("ERROR: SpaceManager, growth factor %f must be >= 1", growthFactor);
    if (growthAddition  < 1)
      error("ERROR: SpaceManager, growth addition %d must be > 0", growthAddition);
    if (startingSize < 1)
      error("ERROR: SpaceManager, starting size %d must be >= 1", startingSize);
    reset();
  }

  // Routine to reset the state to right after the constructor was
  // called.
  void reset() {
    minAllocationSize = currentAllocationSize = startingSize;
    sumRequested = 0.0;
    numberRequests = 0;
  }

  // After each object allocates a series of chunks of memory, it 
  // might be advantages to decay the amount allocated rather than to
  // keep growing indefinitely whenever one of the objects happens to require
  // a large amount. This function decreases the currently allocated amount
  // by the decay rate.
  void decay() {
    currentAllocationSize = (int)((double)currentAllocationSize*decayRate);
    if (currentAllocationSize < startingSize)
      currentAllocationSize = startingSize;
  }

  // Return the current recommended size of memory to allocate.
  unsigned currentSize() {
    return currentAllocationSize;
  }

  // Return the min recommended size of memory to allocate.
  unsigned minSize() {
    return minAllocationSize;
  }

  void recomputeMinSize(unsigned arg) {
  }

  // Bump up to the next new size
  unsigned advanceToNextSize() {
    currentAllocationSize = (int)((double)currentAllocationSize*growthFactor) + 
      growthAddition;
    return currentAllocationSize;
  }


}; 


#endif

