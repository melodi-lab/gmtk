/*
 * GMTK_SectionTablesBase.h
 *   
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (C) 2009 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */


#ifndef GMTK_SECTIONTABLESBASE_H
#define GMTK_SECTIONTABLESBASE_H

#include <vector>
#include <string>
#include <map>

#include <stdio.h>
#include <stdlib.h>

#include "bp_range.h"

#include "GMTK_RV.h"
#include "GMTK_FileParser.h"
#include "GMTK_GMTemplate.h"
#include "GMTK_PartitionStructures.h"
#include "GMTK_JT_Partition.h"
#include "GMTK_MaxCliqueTable.h"

#include "debug.h"

// class mention for forward references.
class GraphicalModel;
class BoundaryTriangulate;
class Partition;
class GMTemplate;
class SectionScheduler;
class SectionIterator;

//
// This object contains only the data tables that are expanded upon
// during inference but that reside within a section. Basically,
// this object contains only the unique information that is relevat to
// a given section as it is expanded during inference, and contains
// no information that is redundant and/or can be deduced by the
// section, clique, and separator structures.  An instance of a
// SectionTable along with a PartitionStructures object together
// constitute all of the information in a section, so these can be
// considered paired, although there will in general be many more
// instances of a SectionTable (the length of the segment) then a
// PartitionStructures (the minimal amount one needs to unroll).  The
// redundant information for all cliques/separators are obtainable via
// corresponding PartitionStructures object.
// 
// This object should be as small as possible to conserve memory since
// we might make an array of these that are very long. We also use
// pointers to make it such that part of a SectionTable (e.g., the
// separatorCliques) may be deallocated and freed and re-constructed
// from the corresponding clique when needed. Also, clique table
// storage is currently (at least as of Jan 2009) more dense (and
// probably more efficient) than separator table storage.
class SectionTablesBase {
 public:

  ConditionalSeparatorTable* separatorCliques;

  // Create a dummy/invalid table for arrays that can be written over.
  SectionTablesBase() : separatorCliques(NULL) {}

  SectionTablesBase(JT_Partition& origin) {
    // first allocate space with empty (and unusable) entries
    separatorCliques = new ConditionalSeparatorTable[origin.separators.size()];
    
    // then actually re-construct the objects in the array appropriately.
    for (unsigned i=0;i<origin.separators.size();i++) {
      new (&separatorCliques[i]) ConditionalSeparatorTable(origin.separators[i]);
    }
  }

  virtual ~SectionTablesBase() { clear(); }

  virtual MaxCliqueTable *getMaxCliques() = 0;

  virtual logpr probEvidence(SectionIterator &inference_it, SectionScheduler &myjt) = 0;

  virtual void projectToOutgoingSeparators(SectionIterator &stss_it,
					   PartitionStructures &sourceSectionStructures, 
					   ConditionalSeparatorTable *separatorTableArray,
					   ConditionalSeparatorTable::SharedLocalStructure *sepSharedStructureArray) = 0;

  virtual void receiveBackwardsSeparators(SectionIterator &stss_it,
					  PartitionStructures &sourceSectionStructures, 
					  ConditionalSeparatorTable *separatorTableArray,
					  ConditionalSeparatorTable::SharedLocalStructure &sepSharedStructure) = 0;

  virtual void printAllCliques(FILE *f, BP_Range *clique_print_range,
			       SectionIterator &stss_it, PartitionStructures &ss,
			       const bool normalize, const bool unlog,
			       const bool justPrintEntropy,
			       ObservationFile *obs_file = NULL) = 0;

  // EM updating.
  virtual void emIncrement(PartitionStructures& ps,
			   const logpr probE, 
			   const bool localCliqueNormalization = false,
			   const double emTrainingBeam = -LZERO) = 0;

  virtual void clear() {
    delete [] separatorCliques;
    separatorCliques = NULL;
  }

  virtual void init(PartitionStructures& ps) = 0;

  // memory use reporting
  virtual void reportMemoryUsageTo(PartitionStructures& ps,FILE *f) = 0;

};


#endif

