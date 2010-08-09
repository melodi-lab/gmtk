/*
 * GMTK_PartitionStructures.h
 *   GMTK Junction Tree. Exact inference support for GMTK.
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2003, < fill in later >
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


#ifndef GMTK_PARTITIONSTRUCTURES_H
#define GMTK_PARTITIONSTRUCTURES_H

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


// A PartitionStructures partition is used to store partition
// information for the final unrolled set of random variables. This
// does not contain information about the clique tables.  
class PartitionStructures {
  friend class JunctionTree;
public:

  // original partition that this has been cloned from.
  JT_Partition& origin;

  // The MaxClique and Separator table's define what shared structure
  // they need, but we keep it here since it is something that will be
  // reused for many different data instances of a MaxCliqueTable and
  // a ConditionalSeparatorTable.
  sArray< MaxCliqueTable::SharedLocalStructure > maxCliquesSharedStructure;
  sArray< ConditionalSeparatorTable::SharedLocalStructure > separatorCliquesSharedStructure;

  // Store factor cliques somewhere that can be shared. (factorCliques
  // don't currently have data associated with them).
  // sArray< InferenceFactorClique > factorCliques;

  // WARNING: constructor hack to create a very broken object with
  // non-functional reference objects (in order to create an array of
  // these objects and then initialize them later with appropriate
  // references). Do not use until after proper re-constructor.
  PartitionStructures() : origin(*((JT_Partition*)NULL)) {}
  // normal (or re-)constructor
  PartitionStructures(JT_Partition& _origin,
		      vector <RV*>& newRvs,
		      map < RVInfo::rvParent, unsigned >& ppf,
		      const unsigned int frameDelta,
		      const bool has_li_separator = true);

  // destructor
  // ~PartitionStructures() {}


  // Return as a set all RVs (and any of their observed parents) that
  // are contained in this partition structure.
  set <RV*> returnRVsAndTheirObservedParentsAsSet();

  void clearCliqueAndIncommingSeparatorMemory();

  // all the random variables in this partition, other than the left
  // interface (i.e., so this is the partition's innovation).
  set<RV*> allrvs;
  // Information needed to compute viterbi values.
  vector <RV*> hidRVVector;
  // direct pointers to the values of all discrete hidden variables in
  // this strucure, ordered.
  sArray <DiscRVType*> hrvValuePtrs;
  // a packer for this structure.
  PackCliqueValue packer; 



#if 0
  // frame delta might be positive or negative
  static void adjustFramesBy(PartitionStructures& part1,
			     PartitionStructures& part2,
			     const int frameDelta);
#endif

};



#endif

