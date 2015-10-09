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
#include "GMTK_SectionScheduler.h"
#include "GMTK_SectionIterator.h"
#include "GMTK_SectionInferenceAlgorithm.h"

#include "GMTK_PartitionTables.h"

class SparseJoinInference : public SectionInferenceAlgorithm {
 public:

  SparseJoinInference(SectionScheduler *jt) : SectionInferenceAlgorithm(jt) {}

  // All message actions are named from the perspective of C_t.

  // compute forward message for C'_t -> C'_{t+1} (aka gather into root)
  SectionSeparator *computeForwardInterfaceSeparator(SectionIterator &t, PartitionTables *section_posterior);

  // recieve forward message for C'_{t-1} -> C'_t (sendForwardsCrossPartitions)
  void receiveForwardInterfaceSeparator(SectionIterator &t, SectionSeparator *msg, PartitionTables *section_posterior);


  // compute backward message for C'_{t-1} <- C'_t (aka scatter out of root)
  SectionSeparator *computeBackwardsInterfaceSeparator(SectionIterator &t);

  // recieve backward message for C'_t <- C'_{t+1} (sendBackwardCrossPartitions)
  void receiveBackwardInterfaceSeparator(SectionIterator &t, SectionSeparator const &msg);

 private:

};

#endif

