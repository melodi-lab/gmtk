/*
 * GMTK_FACTORCLIQUE.h
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *  $Header$
 *
 *
 *   A clique class.
 *   Note: some texts define a 'clique' as any complete set
 *   while other texts define a 'clique' as a maximally
 *   complete set with respect to the subset operator (i.e., a
 *   clique is one such that no proper superset of the set
 *   of nodes is a clique). In order to avoid confusion,
 *   I adopt here the term 'maxclique' which corresponds
 *   to a maximally complete set. Note, however, that in this
 *   program, the concepts are such that 
 *
 *               'clique == maxclique != complete set'
 *
 *   meaning that cliques are taken to be max cliques.
 *
 * Copyright (C) 2001 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 *
 */

// TODO: perhaps create a subclass or member of maxClique at some point, rather than
// adding everything for exact inference to the base class.


#ifndef GMTK_FACTORCLIQUE_H
#define GMTK_FACTORCLIQUE_H

#if HAVE_CONFIG_H
#include <config.h>
#endif


#include "debug.h"

#include "GMTK_RV.h"
#include "GMTK_DiscRV.h"
#include "GMTK_FactorInfo.h"


#include <stdio.h>
#include <stdlib.h>

#include <vector>
#include <set>
#include <map>

class PartitionStructures;
class PartitionTables;
class SeparatorClique;
class ConditionalSeparatorTable;
// class ConditionalSeparatorTable::SharedLocalStructure;
class MaxClique;
class MaxCliqueTable;




// A factor clique is a (not-necessarily max) clique that implements
// some form of directed "factor" or "constraint" within a regular max
// clique. This coresponds to the 'factor' construct in the .str file.
// A factor clique's nodes will, like a separator, necessarily be a
// subset of some MaxClique, but a factor is implemented quite
// differently. Factors mixed with DAG cpts allow the specificatio of
// true hybrid directed/undirected graphical models in GMTK.
class FactorClique : public IM
{
  friend class FileParser;
  friend class GraphicalModel;
  friend class GMTemplate;
  friend class ConditionalSeparatorTable;
  friend class JunctionTree;

public:


  // information about the factor corresponding
  // to this factor
  FactorInfo* factorInfo;

  // the set of nodes that are involved in the factor/constraint.
  set<RV*> nodes;
  // vector of the same nodes, in order that they
  // appear in the .str file, and the order used to index into
  // the factor functions.
  vector<RV*> orderedNodes;

  // copy constructor 
  FactorClique(const FactorClique& factor)
  { 
    factorInfo = factor.factorInfo;
    nodes = factor.nodes; 
    orderedNodes = factor.orderedNodes;
  }
  FactorClique(FactorInfo& factorInfo,
	       vector <RV*>& unrolled_rvs,
	       map < RVInfo::rvParent, unsigned > ppf,
	       const unsigned offset);


  ~FactorClique() {}

  // TODO: These appear to be undefined?

  // prepare the last set of data structures so that clones of this
  // can be unrolled and inference can occur.
  void prepareForUnrolling();

  // print out everything in this clique to a file.
  void printAllJTInfo(FILE* f);

};



#endif
