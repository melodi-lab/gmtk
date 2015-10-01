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
//   LBP

#include "fileParser.h"

#include "GMTK_PtPsIterator.h"
#include "GMTK_SectionSeparator.h"
#include "GMTK_PartitionTables.h"

class SectionInferenceAlgorithm {
 public:

  virtual ~SectionInferenceAlgorithm() {}

  virtual void setJT(JunctionTree *jt) = 0; // BOGUS
  
  // All message actions are named from the perspective of C_t.

  // compute forward message for C'_t -> C'_{t+1} (aka gather into root)
  virtual SectionSeparator *computeForwardInterfaceSeparator(unsigned t, PartitionTables *sectionPosterior) = 0; 

  // recieve forward message for C'_{t-1} -> C'_t (sendForwardsCrossPartitions)
  virtual void receiveForwardInterfaceSeparator(PtPsIterator &t, SectionSeparator *msg, PartitionTables *sectionPosterior) = 0;


  // compute backward message for C'_{t-1} <- C'_t (aka scatter out of root)
  virtual SectionSeparator computeBackwardsInterfaceSeparator(unsigned t) = 0;

  // recieve backward message for C'_t <- C'_{t+1} (sendBackwardCrossPartitions)
  virtual void receiveBackwardInterfaceSeparator(unsigned t, SectionSeparator const &msg) = 0;


  // return P(Q_t | X_{?}), where ? depends on the messages C_t has seen so far:
  //        P(Q_t | X_{0:t}) in the forward pass
  //        P(Q_t | X_{0:T-1}) in a (full) backward pass
  //        P(Q_t | X_{0:t+\tau}) in a smoothing backward pass
  virtual logpr probEvidence(unsigned t, PartitionTables *sectionPosterior) = 0;


  // Section subclasses can manage their own message ordering w/in a section.
  // Read/write the section's inference plan (JT, msg orders, etc)

  void createInferencePlan(); // called by gmtkTriangulate ?
  void readInferencePlan(iDataStreamFile *f);
  void writeInferencePlan(ioDataStreamFile *f);
};

#endif
