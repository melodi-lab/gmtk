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

#include "GMTK_SectionSeparator.h"
#include "GMTK_SectionInferenceAlgorithm.h"
#include "GMTK_PtPsIterator.h"

#include "GMTK_JunctionTree.h"
#include "GMTK_PartitionTables.h"

class SparseJoinInference : public SectionInferenceAlgorithm {

  // temporary bogusness
  void setJT(JunctionTree *jt) {myjt = jt;} // BOGUS

  // All message actions are named from the perspective of C_t.

  // compute forward message for C'_t -> C'_{t+1} (aka gather into root)
  SectionSeparator *computeForwardInterfaceSeparator(unsigned t, PartitionTables *sectionPosterior) {
    assert(myjt);
    return myjt->computeForwardInterfaceSeparator(t, sectionPosterior);
  } 

  // recieve forward message for C'_{t-1} -> C'_t (sendForwardsCrossPartitions)
  void receiveForwardInterfaceSeparator(PtPsIterator &iit, SectionSeparator *msg, PartitionTables *sectionPosterior) {
    myjt->recieveForwardInterfaceSeparator(iit, msg, sectionPosterior);
  }


  // compute backward message for C'_{t-1} <- C'_t (aka scatter out of root)
  SectionSeparator computeBackwardsInterfaceSeparator(unsigned t) {
    SectionSeparator is;
    return is;
  } 

  // recieve backward message for C'_t <- C'_{t+1} (sendBackwardCrossPartitions)
  void receiveBackwardInterfaceSeparator(unsigned t, SectionSeparator const &msg) {}


  // return P(Q_t | X_{?}), where ? depends on the messages C_t has seen so far:
  //        P(Q_t | X_{0:t}) in the forward pass
  //        P(Q_t | X_{0:T-1}) in a (full) backward pass
  //        P(Q_t | X_{0:t+\tau}) in a smoothing backward pass
  logpr probEvidence(unsigned t, PartitionTables *sectionPosterior) {
    return myjt->computeProbEvidence(t, sectionPosterior);
  }

 private:
  JunctionTree *myjt;

};

#endif

