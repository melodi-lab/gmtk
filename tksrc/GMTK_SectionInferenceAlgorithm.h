/*
 * GMTK_SectionInferenceAlgorithm.h
 *   Base class for within-section inference algorithms.
 *
 * Written by Richard Rogers <rprogers@uw.edu>
 *
 * Copyright (C) 2015 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

#ifndef GMTK_SECTIONINFERENCEALGORITHM_H
#define GMTK_SECTIONINFERENCEALGORITHM_H


// Specific section inference algorithm implementations:

//   Pedagogical 
//   Sparse join

#include "fileParser.h"

#include "GMTK_InterfaceSeparator.h"

class SectionInferenceAlgorithm {
 public:

  virtual ~SectionInferenceAlgorithm() {}

  // Prepare to do inference on segment of length T. Note that this
  // is the analog of the current JunctionTree::unroll() which sets up the
  // PartitionStructureArray for at most  P' C' C' E' (4 sections, so
  // O(1) memory), and the PartitionTableArray for between 0 and T
  // sections depending on ZeroTable, ShortTable, or LongTable. 
  // JT::unroll() also allocates O(T) memory (optionally memory mapped)
  // to store the the Viterbi values if they aren't being written to
  // a file.

  // Since we now have the ability to write Viterbi values to files, do
  // we want to drop the ability to store them in memory? This would
  // break things, as a Viterbi file would then be a required argument.
  // It would also be slower for applications that do fit in memory
  // because of the higher overhead of file I/O operations.

  // Maybe rename this to prepareForSegment() or something to avoid confusion
  // with O(T) space graph unrolling?
  virtual unsigned unroll(unsigned T) = 0; // returns number of usable frames


  // All message actions are named from the perspective of C_t.

  // compute forward message for C'_t -> C'_{t+1} (aka gather into root)
  virtual InterfaceSeparator computeForwardInterfaceSeparator(unsigned t) = 0; 

  // recieve forward message for C'_{t-1} -> C'_t (sendForwardsCrossPartitions)
  virtual void receiveForwardInterfaceSeparator(unsigned t, InterfaceSeparator const &msg) = 0;


  // compute backward message for C'_{t-1} <- C'_t (aka scatter out of root)
  virtual InterfaceSeparator computeBackwardsInterfaceSeparator(unsigned t) = 0;

  // recieve backward message for C'_t <- C'_{t+1} (sendBackwardCrossPartitions)
  virtual void receiveBackwardInterfaceSeparator(unsigned t, InterfaceSeparator const &msg) = 0;


  // return P(Q_t | X_{?}), where ? depends on the messages C_t has seen so far:
  //        P(Q_t | X_{0:t}) in the forward pass
  //        P(Q_t | X_{0:T-1}) in a (full) backward pass
  //        P(Q_t | X_{0:t+\tau}) in a smoothing backward pass
  virtual logpr probEvidence(unsigned t) = 0;


  // Section subclasses can manage their own message ordering w/in a section.
  // Read/write the section's inference plan (JT, msg orders, etc)

  void createInferencePlan(); // called by gmtkTriangulate ?
  void readInferencePlan(iDataStreamFile *f);
  void writeInferencePlan(ioDataStreamFile *f);
};

#endif
