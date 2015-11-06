/*
 * GMTK_CliqueValueHolder.h
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *  $Header$
 *
 * Copyright (C) 2001 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 *
 */

#ifndef GMTK_CLIQUEVALUEHOLDER_H
#define GMTK_CLIQUEVALUEHOLDER_H

#if HAVE_CONFIG_H
#include <config.h>
#endif

#if 0
#include "general.h"
#include "vhash_set.h"
#include "vhash_map.h"
#include "logp.h"
#endif

#include "cArray.h"
#include "sArray.h"
#include "debug.h"

#if 0
#include "fixed_filter.h"
#include "lms_filter.h"
#include "rls_filter.h"
#include "counted_ptr.h"

#include "GMTK_RV.h"
#include "GMTK_DiscRV.h"
#include "GMTK_PackCliqueValue.h"
#include "GMTK_SpaceManager.h"
#include "GMTK_FactorInfo.h"
#include "GMTK_ObservationFile.h"

#include <stdio.h>
#include <stdlib.h>

#include <vector>
#include <string>
#include <map>
#endif


class CliqueValueHolder  {

  // The packed clique value size (in unsigned).
  // @@@ sould be const
  // const unsigned cliqueValueSize;
  unsigned cliqueValueSize;

  // The amount that successive allocation chunks grow in size.  Must
  // be >= 1.0.
  // Ideally, should be a const.
  float growthFactor;
  
  // The size of the initial allocation unit. A chunk is an
  // array of packed clique values. When we allocate the n'th
  // chunk, we allocate k(n) new packed clique values, 
  // where k(n) = allocationUnitChunkSize*growthFactor^(n-1)
  // Ideally, should be a const.
  unsigned allocationUnitChunkSize;

  // manage the space here.
  // SpaceManager spaceManager;

  // Current capacity, total number of unsigned values that this
  // object can currently potentially hold max without a resize.
  unsigned capacity;

  // a chunk, i.e., matrix of unsigned numbers constituting
  // allocationUnitChunkSize*growthFactor^(n-1) packed clique values
  // for some n.
  typedef sArray< unsigned > AllocationChunk;

  // Array of chunks.
  cArray< AllocationChunk >  values;

  // The the pointernext position in the current chunk to obtain a clique value
  // to use.
  unsigned* curAllocationPosition;

  // the end position in the current chunk, meaning
  // that we need to reallocate
  unsigned* curAllocationEnd;

  // total number of values currently used (not really
  // needed but kept anyway for debugging).
  unsigned numAllocated;


public:

  // this is the default value for growthFactor, set by -memoryGrowth
  static float defaultGrowthFactor;

  // this is the default value for allocationUnitChunkSize, set by -memoryGrowth
  static unsigned defaultAllocationUnitChunkSize;

  // create an empty object to re-construct later
  CliqueValueHolder() {}
  
  // real constructor
  CliqueValueHolder(unsigned cliqueValueSize);

  CliqueValueHolder(unsigned cliqueValueSize,
		    unsigned allocationUnitChunkSize,
		    float growthFactor=1.25);

  ~CliqueValueHolder() { makeEmpty();  }

  // clear out all existing memory, and get ready for next use.
  void prepare();

  // Empty out and free up all memory, and reset to having
  // nothing added. Make invalid as well.
  void makeEmpty();

  // return a pointer to the next unused clique value
  // for scratch (etc.) without actually allocating it.
  unsigned* curCliqueValuePtr() { return curAllocationPosition; }

  // Actually "allocate" the space pointed to by curCliqueValuePtr()
  // and then advance curCliqueValuePtr() to the next position.
  void allocateCurCliqueValue();

  // return the total number of bytes requested to the OS memory system by this structure.
  unsigned long bytesRequested() { return (unsigned long)capacity*sizeof(unsigned); }

};


#endif
