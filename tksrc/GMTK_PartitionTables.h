/*
 * GMTK_PartitionTables.h
 *   
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2009, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about the suitability of this software
 * for any purpose. It is provided "as is" without express or implied warranty.
 *
 *
 * $Header$
 *
 */


#ifndef GMTK_PARTITIONTABLES_H
#define GMTK_PARTITIONTABLES_H

#include <vector>
#include <string>
#include <map>

#include <stdio.h>
#include <stdlib.h>

#include "bp_range.h"

#include "GMTK_RV.h"
#include "GMTK_FileParser.h"
#include "GMTK_GMTemplate.h"
#include "GMTK_MaxClique.h"
#include "GMTK_JT_Partition.h"

#include "debug.h"

// class mention for forward references.
class GraphicalModel;
class BoundaryTriangulate;
class Partition;
class GMTemplate;
class JunctionTree;


//
// This object contains only the data tables that are expanded upon
// during inference but that reside within a partition. Basically,
// this object contains only the unique information that is relevat to
// a given partition as it is expanded during inference, and contains
// no information that is redundant and/or can be deduced by the
// partition, clique, and separator structures.  An instance of a
// PartitionTable along with a PartitionStructures object together
// constitute all of the information in a partition, so these can be
// considered paired, although there will in general be many more
// instances of a PartitionTable (the length of the segment) then a
// PartitionStructures (the minimal amount one needs to unroll).  The
// redundant information for all cliques/separators are obtainable via
// corresponding PartitionStructures object.
// 
// This object should be as small as possible to conserve memory since
// we might make an array of these that are very long. We also use
// pointers to make it such that part of a PartitionTable (e.g., the
// separatorCliques) may be deallocated and freed and re-constructed
// from the corresponding clique when needed. Also, clique table
// storage is currently (at least as of Jan 2009) more dense (and
// probably more efficient) than separator table storage.
struct PartitionTables {

  MaxCliqueTable* maxCliques;
  ConditionalSeparatorTable* separatorCliques;

  // Create a dummy/invalid table for arrays that can be written over.
  PartitionTables() 
  : maxCliques(NULL), separatorCliques(NULL) {}

  PartitionTables(JT_Partition& origin);
  ~PartitionTables() { clear(); }

  // EM updating.
  void emIncrement(PartitionStructures& ps,
		   const logpr probE, 
		   const bool localCliqueNormalization = false,
		   const double emTrainingBeam = -LZERO);

  void clear() {
    delete [] maxCliques;
    maxCliques = NULL;
    delete [] separatorCliques;
    separatorCliques = NULL;
  }
  void init(PartitionStructures& ps);

  // memory use reporting
  void reportMemoryUsageTo(PartitionStructures& ps,FILE *f);

};


#endif

