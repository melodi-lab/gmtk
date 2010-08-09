/*
 * GMTK_JunctionTree.h
 *
 *   JT_Partition is a partition that has JT-related and inference members.
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

#ifndef GMTK_JT_PARTITION_H
#define GMTK_JT_PARTITION_H

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


#include "debug.h"

// class mention for forward references.
class GraphicalModel;
class BoundaryTriangulate;
class Partition;
class GMTemplate;
class JunctionTree;

// TODO: perhaps create a subclass of maxClique at some point, rather than
// adding everything for exact inference to the base class.

// Child class of partition that includes support for doing exact
// inference.
class JT_Partition : public Partition {

  friend class JunctionTree;

  void findInterfaceCliques(const set <RV*>& iNodes,
			    unsigned& iClique,
			    bool& iCliqueSameAsInterface,
			    const string priorityStr);
public:


  // Interface nodes on the "left" of this partition. I.e., To find
  // the left interface clique, find a clique that is a superset of
  // these nodes. Empty if there is no such set (e.g., for a P
  // partition)
  set <RV*> liNodes;

  // Interface nodes on the "right" of this partition. I.e., to
  // compute the root clique of this partition, we find a clique that
  // is a superset of these nodes. Empty if there is no such set
  // (e.g., for an E partition).
  set <RV*> riNodes;

  // Nodes that are not assigned in this partition. If all nodes are
  // forward-time directed, we are guaranteed that they will be
  // assigned in the left adjacent partition. Similarly, if all nodes
  // are backward-time directed, the nodes are assigned in the right
  // adjacent partition. With a bi-directional graph, the nodes could
  // be assigned in either the left or right adjacent partition.
  set <RV*> unassignedInPartition;
  
  // The separators for this partition.  If this is a P partition,
  // then all of the separators in this partition are between cliques
  // that live entirely within this partition.  If this is a C or an E
  // partition, then most of the sperators are between cliques that
  // live entirely within this partition.  The last separator
  // (left-interface or LI separator) in this vector is guaranteed to
  // be the seperator between the right interface (RI) clique of the
  // adjacent partition on the left (if it exists), and the left
  // interface (LI) clique of this partition. This last separator may
  // or may not exist because if this is a P partition, there is no
  // left interface separator, but there is one for a C or an E
  // partition.
  // 
  // Created in: JunctionTree::createSeparators(); 
  // 
  // This final separator is called the LI separator.  Note there also
  // might be VE-separators here as well. The order that the
  // separators are: 1) normal separators, 2) VE-separators, and 3)
  // the final LI separator.
  vector<SeparatorClique> separators;

  // The set of factor cliques (i.e., hard and soft constraints) that
  // live in this partition, corresponding to the 'factor' constructs
  // in the .str file.
  vector<FactorClique> factorCliques;

  void useLISeparator()  { separators[separators.size()-1].skipMe = false; }
  void skipLISeparator() { separators[separators.size()-1].skipMe = true; }

  // number of VE separators among the separators in this partition.
  unsigned numVEseps;

  // create an empty one to be filled in later.
  JT_Partition() {}

  // constructor
  JT_Partition(Partition& from_part,
	       const unsigned int frameDelta,
	       // the left and right interface variables for
	       // this JT partition Empty if doesn't exist
	       // (say for an P or E partition). These have
	       // their own frame deltas since they might be
	       // different.
	       const set <RV*>& from_liVars,
	       const unsigned int liFrameDelta,
	       const set <RV*>& from_riVars,
	       const unsigned int riFrameDelta,
	       // Information todo the mapping.
	       vector <RV*>& newRvs,
	       map < RVInfo::rvParent, unsigned >& ppf);

  JT_Partition(Partition& from_part,
	       const set <RV*>& from_liVars,
	       const set <RV*>& from_riVars);
  

  // returns the left and right interface clique. If not defined,
  // sets the variable to ~0x0.
  void findLInterfaceClique(unsigned& liClique,bool& liCliqueSameAsInterface,
			    const string priorityStr);
  void findRInterfaceClique(unsigned& riClique,bool& riCliqueSameAsInterface,
			    const string priorityStr);

  // return the index of the clique with max/min weight.
  unsigned cliqueWithMaxWeight();
  unsigned cliqueWithMinWeight();


  void clearCliqueSepValueCache(bool force = false) {
    for (unsigned i=0;i<cliques.size();i++)
      cliques[i].clearCliqueValueCache(force);
    for (unsigned i=0;i<separators.size();i++) 
      separators[i].clearSeparatorValueCache(force);
  }

  void clearCliqueValueCache(bool force=false) {
    for (unsigned i=0;i<cliques.size();i++)
      cliques[i].clearCliqueValueCache(force);
  }

  void clearSeparatorValueCache(bool force=false) {
    for (unsigned i=0;i<separators.size();i++) 
      separators[i].clearSeparatorValueCache(force);
  }

  void clearCliqueAndIncommingSeparatorMemoryForClique(unsigned cliqueNo);

};


#endif

