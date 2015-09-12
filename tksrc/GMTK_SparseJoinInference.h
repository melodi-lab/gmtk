/*
 * GMTK_SparseJoinInference.h
 *   Efficient Hugin-style message passing inference algorithm within sections.
 *
 * Written by Richard Rogers <rprogers@uw.edu>
 *
 * Copyright (C) 2015 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */


#ifndef GMTK_SPARSEJOININFERENCE_H
#define GMTK_SPARSEJOININFERENCE_H

#include "GMTK_InterfaceSeparator.h"
#include "GMTK_SectionInferenceAlgorithm.h"

class SparseJoinInference : public SectionInferenceAlgorithm {

  // Prepare to do inference on segment of length T.
  virtual unsigned unroll(unsigned T) {return 0;} // returns number of usable frames


  // All message actions are named from the perspective of C_t.

  // compute forward message for C'_t -> C'_{t+1} (aka gather into root)
  virtual InterfaceSeparator computeForwardInterfaceSeparator(unsigned t) {
    InterfaceSeparator is;
    return is;
  } 

  // recieve forward message for C'_{t-1} -> C'_t (sendForwardsCrossPartitions)
  virtual void receiveForwardInterfaceSeparator(unsigned t, InterfaceSeparator const &msg) {}


  // compute backward message for C'_{t-1} <- C'_t (aka scatter out of root)
  virtual InterfaceSeparator computeBackwardsInterfaceSeparator(unsigned t) {
    InterfaceSeparator is;
    return is;
  } 

  // recieve backward message for C'_t <- C'_{t+1} (sendBackwardCrossPartitions)
  virtual void receiveBackwardInterfaceSeparator(unsigned t, InterfaceSeparator const &msg) {}


  // return P(Q_t | X_{?}), where ? depends on the messages C_t has seen so far:
  //        P(Q_t | X_{0:t}) in the forward pass
  //        P(Q_t | X_{0:T-1}) in a (full) backward pass
  //        P(Q_t | X_{0:t+\tau}) in a smoothing backward pass
  virtual logpr probEvidence(unsigned t) {
    logpr p;
    return p;
  }

};

#endif

