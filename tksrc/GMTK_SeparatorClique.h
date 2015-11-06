/*
 * GMTK_MAXCLIQUE.h
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


#ifndef GMTK_SEPARATORCLIQUE_H
#define GMTK_SEPARATORCLIQUE_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "general.h"
#include "vhash_set.h"
#include "vhash_map.h"

#include "cArray.h"
#include "sArray.h"
#include "debug.h"

#include "GMTK_RV.h"
#include "GMTK_DiscRV.h"
#include "GMTK_PackCliqueValue.h"
#include "GMTK_SpaceManager.h"



#include "GMTK_MaxClique.h"
#include "GMTK_ConditionalSeparatorTable.h"

#include <stdio.h>
#include <stdlib.h>

#include <vector>
#include <map>


class SeparatorClique : public IM
{
  friend class FileParser;
  friend class GraphicalModel;
  friend class GMTemplate;
  friend class ConditionalSeparatorTable;
  friend class JunctionTree;

public:

  // beam width for separator-based beam pruning.
  static double separatorBeam;

  // the set of nodes which forms the separator
  set<RV*> nodes;

  //////////////////////////////////////////////////////////////////////
  // For a given iteration order, the intersection of 'nodes' with all
  // nodes in previous seperators according to the iteration.
  set<RV*> accumulatedIntersection;
  // hidden variables of the above
  vector<RV*> hAccumulatedIntersection;
  // structure used to pack and unpack clique values
  PackCliqueValue accPacker;
  // structure to allocate clique values from.  These things are
  // allocated and deleted in bulk (rather than individually) and that
  // is handled by this object, thereby acting effectively as a
  // customized memory management unit.
  CliqueValueHolder accValueHolder;
  // hash table that holds the clique values (a packed vector of RV
  // values). The actual values are pointers to within valueHolder.
  vhash_set< unsigned > accSepValHashSet;


  //////////////////////////////////////////////////////////////////////
  // remainder = nodes - accumulatedIntersection which is precomputed
  // in the non-unrolled version of a separator.
  set<RV*> remainder;
  // hidden variables of the above
  vector<RV*> hRemainder;
  // structure used to pack and unpack clique values
  PackCliqueValue remPacker;
  // structure to allocate clique values from.  These things are
  // allocated and deleted in bulk (rather than individually) and that
  // is handled by this object, thereby acting effectively as a
  // customized memory management unit.
  CliqueValueHolder remValueHolder;
  // hash table that holds the clique values (a packed vector of RV
  // values). The actual values are pointers to within valueHolder.
  vhash_set< unsigned > remSepValHashSet;

  ////////////////////////////////////////////////////////////
  // Data structures for when this is a VE separator. 
  // set to true if this is a virtual evidence (VE) separator.
  bool veSeparator;
  // information about this VE separator
  MaxClique::VESepInfo veSepInfo;
  // a pointer to the ve sep clique that contains the actual probabily table.
  // This is computed in prepareForUnrolling(). Once it is computed,
  // we are not allowed to make any additional copies of this object,
  // or otherwise this pointer might get deleted twice.
  ConditionalSeparatorTable* veSepClique;
  ///////////////////////////////////////////////
  // VE separator files information.
  ///////////////////////////////////////////////
  // Command line: recompute the VE separator tables and save to disk in all cases.
  static bool recomputeVESeparatorTables;
  // File name to read/write VE separator table.
  static const char* veSeparatorFileName;
  // set to true if we are (re-)generating the VE tables,
  // or set to false if we are just reading them in from disk.
  static bool generatingVESeparatorTables;
  // actual file to get ve sep stuff.
  static FILE* veSeparatorFile;
  // The log (base 10) upper limit on a VE sep variable cardinality
  // product. I.e., if the number of parents that need to be iterated
  // over to produce the VE sep table has a prod. of cardinalties
  // greater than this (in log base 10), we don't use this VE sep.
  static float veSeparatorLogProdCardLimit;


  // Memory management parameters set by -memoryGrowth
  static unsigned aiStartingSize;
  static float    aiGrowthFactor;
  static unsigned remStartingSize;
  static float    remGrowthFactor;
  static unsigned sepSpaceMgrStartingSize;
  static float    sepSpaceMgrGrowthRate;
  static float    sepSpaceMgrDecayRate;
  static unsigned remSpaceMgrStartingSize;
  static float    remSpaceMgrGrowthRate;
  static float    remSpaceMgrDecayRate;

  // A boolean flag that the inference code uses to determine if it
  // should skip this separator. This is used when a P partition is
  // empty and the C partition needs to skip its incomming
  // interface separator (which is empty and if it wasn't
  // skipped would lead to a zero probability).
  bool skipMe;

  // copy constructor 
  SeparatorClique(const SeparatorClique& sep)
    : veSeparator(sep.veSeparator)
  { 
    // this constructor only copies the non-filled out information
    // (nodes and veSep status and information) since the other stuff
    // isn't needed for this type of constructor. Note that if other
    // code changes, we might need to add more here.
    nodes = sep.nodes; 
    veSepInfo = sep.veSepInfo;
    // make sure this is NULL as if we copy this in, it will get
    // deleted twice.
    assert ( sep.veSepClique == NULL );;
    veSepClique = NULL;
    skipMe = sep.skipMe;
  }

  // constructor for VE separators.
  SeparatorClique(const MaxClique::VESepInfo& _veSepInfo)
    : veSeparator(true)
  { 
    veSepInfo = _veSepInfo;
    // need nodes to reflect union, to sort, etc.
    for (unsigned i=0;i<veSepInfo.parents.size();i++) {
      nodes.insert(veSepInfo.parents[i]);
    }
    // child is guaranteed not to be NULL.
    assert ( veSepInfo.child != NULL );
    if (veSepInfo.grandChild == NULL) {
      // then this is a PC case.
      // Here, we do not insert the child since
      // it is not officially part of the separator.
    } else {
      // the PCG case. Here we insert the child, since the separator
      // is relative to the grandchild, so we iterate over both the
      // parents and the child. 
      nodes.insert(veSepInfo.child);
      // We do not insert the grandchild since it is not officially
      // part of the separator.
    }
    veSepClique = NULL;
    skipMe = false;
  }

  // construct a separator between two cliques
  SeparatorClique(MaxClique& c1, MaxClique& c2);

#if 0
  // not used any longer. Keep here in case we 
  // need it again in future.
  SeparatorClique(SeparatorClique& from_sep,
		  vector <RV*>& newRvs,
		  map < RVInfo::rvParent, unsigned >& ppf,
		  const unsigned int frameDelta = 0);
#endif
  
  // Create an empty separator, used for re-construction later.
  // Hopefully STL doesn't allocate anything when using default
  // constructor, as if it does, it will be lost. Note that we do not
  // repeatedly construct one of these objects, so while we might have
  // a bit of lost memory as a result of this, it won't constitute an
  // ever-growing memory leak.
  SeparatorClique() : veSeparator(false),veSepClique(NULL) {}

  ~SeparatorClique();

  // compute the weight (log10 state space) of this separator clique.
  float weight(const bool useDeterminism = true) const { 
    return MaxClique::computeWeight(nodes,NULL,useDeterminism); 
  }

  // prepare the last set of data structures so that clones of this
  // can be unrolled and inference can occur.
  void prepareForUnrolling();

  // print out everything in this clique to a file.
  void printAllJTInfo(FILE* f);

  // memory reporting
  void reportMemoryUsageTo(FILE *f);

  // Manages and memorizes the size and space requests made by all
  // corresponding ConditionalSeparatorTables regarding the size of
  // 'separatorValues' array. I.e., the next time a
  // ConditionalSeparatorTable asks for an initial amount of memory,
  // we'll be able to give it something closer to what it previously
  // asked for rather than something too small.
  SpaceManager separatorValueSpaceManager;

  // TODO: Manages and memorizes space for the remainder portion of a
  // separator. Note that this is used quite differently as the
  // 'separatorValueSpaceManager' variable, as a potentially different
  // sized remainder exists for all accumulated intersection
  // size. This therefore only keeps track of the maximum size.
  SpaceManager remainderValueSpaceManager; 

  // used to clear out hash table memory between segments
  void clearSeparatorValueCache(bool force=false);

  void clearInferenceMemory() {
    accSepValHashSet.clear();
    remSepValHashSet.clear();
    // could change to makeEmpty if running on a static graph only 
    if (accPacker.packedLen() > ISC_NWWOH_AI)
      accValueHolder.prepare();
    if (remPacker.packedLen() > ISC_NWWOH_RM)
      remValueHolder.prepare();
  }


  // @@@ need to take out, here for now to satisify STL call of vector.clear().
#if 0
  SeparatorClique& operator=(const SeparatorClique& f) {
    return *this;
  }
#endif

};


#endif
