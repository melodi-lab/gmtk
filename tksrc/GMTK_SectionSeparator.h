/*
 * GMTK_SectionSeparator.h
 *   The "message" passed between sections
 *
 * Eventially subclasses may include:
 *  - completed single clique (current implementation)
 *  - junction tree
 *  - general factor graph
 *  - no edges (mean field approximation) 
 *
 * SectionInferenceAlgorithm classes must know how to send/recieve
 * the necessary messages to/from a Section.
 *
 * Written by Richard Rogers <rprogers@uw.edu>
 *
 * Copyright (C) 2015 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */


#ifndef GMTK_SECTIONSEPARATOR_H
#define GMTK_SECTIONSEPARATOR_H

#include "GMTK_PartitionTables.h"

typedef PartitionTables SectionSeparator;

#if 0
class SectionSeparator {

  // What's the API for a SectionInferenceAlgorithm to construct an SectionSeparator?
  // What's the API for a SectionInferenceAlgorithm to receive an SectionSeparator?
};
#endif

#endif
