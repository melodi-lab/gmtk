/*
 * GMTK_SparseJoinInference.cc
 *   Efficient Hugin-style message passing inference algorithm within sections.
 *
 * Written by Richard Rogers <rprogers@uw.edu>
 *
 * Copyright (C) 2015 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

#include "GMTK_SparseJoinInference.h"


// compute forward message for C'_t -> C'_{t+1} (aka gather into root)
SectionSeparator *
SparseJoinInference::computeForwardInterfaceSeparator(SectionIterator &t, PartitionTables *section_posterior) {
  return section_posterior;
} 

// recieve forward message for C'_{t-1} -> C'_t (sendForwardsCrossPartitions)
void 
SparseJoinInference::receiveForwardInterfaceSeparator(SectionIterator &t, SectionSeparator *msg, PartitionTables *section_posterior) {
}


// compute backward message for C'_{t-1} <- C'_t (aka scatter out of root)
SectionSeparator *
SparseJoinInference::computeBackwardsInterfaceSeparator(SectionIterator &t) {
  return NULL;
} 


// recieve backward message for C'_t <- C'_{t+1} (sendBackwardCrossPartitions)
void 
SparseJoinInference::receiveBackwardInterfaceSeparator(SectionIterator &t, SectionSeparator const &msg) {
}

