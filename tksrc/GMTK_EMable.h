/*-
 * GMTK_EMable
 *        Objects that inherit from this are "EMable" in the sense
 *        that we can use them in EM. This class provides some support
 *        for such objects.
 *
 *  Written by Jeff Bilmes <bilmes@ee.washington.edu>
 * 
 *  $Header$
 * 
 * Copyright (c) 2001, < fill in later >
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


#ifndef GMTK_EMABLE_H
#define GMTK_EMABLE_H

#include "fileParser.h"
#include "logp.h"

#include "GMTK_RandomVariable.h"

class EMable {


protected:

  //////////////////////////////////////////////////////////////////
  // A bitmask giving the "state" of the object. Used for
  //   1) error checking
  //   2) swapping parameters when EMables are shared (so
  //      that parameters aren't swapped twice)
  //   3) Loading/Storing accumulators, to make sure they aren't
  //      loaded/stored more than one time.
  enum {
    // Basic data structures allocated (after a read)
    bm_basicAllocated  = (1 << 0),
    // Are the em structures allocated
    bm_emAllocated     = (1 << 1),
    // True when an EM epoch is ongoing, so that we don't call
    // startEmEpoch multiple times when things are tied together.
    bm_emEpochOnGoing  = (1 << 2),
    // Have the current and next parameters been swapped.
    bm_swapped         = (1 << 3),
    // Has the accumulators been loaded/stored.
    bm_accLoadStore    = (1 << 4),

    // Initial State, where no data structures have been allocated.
    bm_initState      = 0x0
  };

  unsigned int bitmask;
 
public:

  
  EMable() { bitmask = 0x0; }
  virtual ~EMable() {}


  //////////////////////////////////
  // EM training                  //
  //////////////////////////////////

  // an "iteration" is an entire pass through the data.

  ////////////////////////////////////////////////////////////////////
  // begins a new epoch. Also ensures that data for EM is allocated
  // and possibly changes the alocated bit, and ongoing bit
  virtual void emStartIteration() = 0;

  ////////////////////////////////////////////////////////////////////
  // Accumulate new data into the internal structures for eam.
  // assumes that the ongoing bit it set.
  virtual void emIncrement(RandomVariable* rv, logpr prob) = 0;

  ////////////////////////////////////////////////////////////////////
  // Accumulate new data into the internal structures for this 
  // em iteration, clears the ongoing bit.
  virtual void emEndIteration() = 0;

  ////////////////////////////////////////////////////////////////////////////
  // if swap bit not set, swaps the current and new parameters, set swap bit.
  // otherwise does nothing.
  virtual void emSwapCurAndNew() = 0;

  /////////////////////////////////////////////////////////////////
  // clear the swap bit, needed for sharing.
  void emClearSwappedBit() { bitmask &= ~bm_swapped; }
  void emSetSwappedBit() { bitmask |= bm_swapped; }

  /////////////////////////////////////////////////////////////////
  // clear the swap bit, needed for sharing.
  void emClearAllocatedBit() { bitmask &= ~bm_emAllocated; }
  void emSetAllocatedBit() { bitmask |= bm_emAllocated; }


  /////////////////////////////////////////////////////////////////
  // clear the swap bit, needed for sharing.
  void emClearOnGoingBit() { bitmask &= ~bm_emEpochOnGoing; }
  void emSetOnGoingBit() { bitmask |= bm_emEpochOnGoing; }

  //////////////////////////////////////////////
  // For parallel EM training.

  ///////////////////////////////////////////////////////////////
  // store the current set of accumulators to a file.
  virtual void emStoreAccumulators(oDataStreamFile& ofile) = 0;

  ///////////////////////////////////////////////////////////////
  // load the current set of accumulators from a file.
  virtual void emLoadAccumulators(iDataStreamFile& ifile) = 0;

  //////////////////////////////////////////////////////////////////////
  // accumulate (add to) the current set of accumulators from a file.
  virtual void emAccumulateAccumulators(iDataStreamFile& ifile) = 0;


  //////////////////////////////////////////////



};


#endif
