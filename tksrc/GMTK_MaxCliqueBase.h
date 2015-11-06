/*
 * GMTK_MAXCLIQUEBASE.h
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


#ifndef GMTK_MAXCLIQUEBASE_H
#define GMTK_MAXCLIQUEBASE_H

#if HAVE_CONFIG_H
#include <config.h>
#endif


#include "general.h"
#include "vhash_set.h"
#include "vhash_map.h"
#include "logp.h"
#include "cArray.h"
#include "sArray.h"
#include "debug.h"
#include "fixed_filter.h"
#include "lms_filter.h"
#include "rls_filter.h"
#include "counted_ptr.h"

#include "GMTK_RV.h"
#include "GMTK_DiscRV.h"
#include "GMTK_PackCliqueValue.h"
#include "GMTK_SpaceManager.h"
#include "GMTK_FactorInfo.h"
#include "GMTK_ObservationFile.h"

#include <stdio.h>
#include <stdlib.h>

#include <vector>
#include <string>
#include <map>

class PartitionStructures;
class PartitionTables;
class SeparatorClique;
class ConditionalSeparatorTable;
// class ConditionalSeparatorTable::SharedLocalStructure;
class MaxClique;
class MaxCliqueTable;

#include "GMTK_CliqueValueHolder.h"
#include "GMTK_VHashMapUnsignedUnsignedKeyUpdatable.h"

//////////////////////////////////////////////////////////////////////////////
// Number of words to exist in a packed clique without a hash. These
// may be set to zero to turn on hashing even when packed clique
// values take up less than 1 machine word. To turn off hashing
// always, set to something larger than the largest packed clique
// value in machine words (e.g., 2, 3, etc. but this might consume
// lots of memory). The number of bits required for a packed clique
// value is equal to:
//
//    num_bits_required = \sum_{v \in h(C)} ceil(log2(card(v)))
//
// where C is a clique, h(C) is all *hidden* variables in the clique C
// (i.e., vars with card > 1), and card(v) is the cardinality of the
// variable v (the card is taken from the structure file). Note, we
// use h(C) so that observed variables (or vars with card = 1) are not
// needlessly stored in the clique value.
// --
// MaxCliqueTable Number Words WithOut a Hash (IMC_NWWOH): Namely,
// the number of words that can be stored directly as a packed clique
// value before we resort to using a shared hash table for all
// instances of this origin clique.
#define IMC_NWWOH (1)
// ConditionalSeparatorTable Number Words WithOut a Hash (ISC_NWWOH):
// Namely, the number of words that can be stored directly as a packed
// clique value before we resort to using a shared hash table for all
// instances of this origin clique.  One for the accumulated
// Intersection packed values
#define ISC_NWWOH_AI (1)
// And the same for the remainder in a Separator.
#define ISC_NWWOH_RM (1)
// -- 
//////////////////////////////////////////////////////////////////////////////





/*
 *
 * Static memory for use for packed clique values. We do this here
 * rather than placing things on the stack since some routines might
 * do it repeatly and recursively. Also, making this global will make
 * the inference code below non-reentrant. TODO: change this when
 * getting working with POSIX threads.
 *
 */


namespace CliqueBuffer {
  // TODO: change this for multi-threading.

  // allocate some temporary storage for packed separator values.
  // 128 words is *much* bigger than any possible packed clique value
  // will take on, but it is easy/fast to allocate on the stack right now.
  extern unsigned packedVal[128];
}




class MaxCliqueBase : public IM {

  friend class FileParser;
  friend class GraphicalModel;
  friend class GMTemplate;
  friend class SeparatorClique;

 public:

  // memory management options set by -memoryGrowth
  static unsigned spaceMgrStartingSize;
  static float    spaceMgrGrowthRate;
  static float    spaceMgrDecayRate;

  // Thresholds for -cpbeam pruning. We keep previous clique max value
  // and previous previous clique max value around.
  logpr prevMaxCEValue;
  logpr prevPrevMaxCEValue;
  counted_ptr<AdaptiveFilter> maxCEValuePredictor;

  double prevMaxCEValPrediction;


  // beam width for clique-based beam pruning.
  static double cliqueBeam;
  // beam width to use while building cliques for partial clique pruning.
  static double cliqueBeamBuildBeam;
  // properties of the filter used for partial clique pruning.
  static char* cliqueBeamBuildFilter;
  // use clique cont heuristic
  static bool cliqueBeamContinuationHeuristic;
  // clique build beam expansion factor 
  static double cliqueBeamBuildExpansionFactor;
  // clique build beam maximum number of expansions
  static unsigned cliqueBeamBuildMaxExpansions;
  // clique beam cluster pruning number of clsuters, 0 or 1 to turn off.
  static unsigned cliqueBeamClusterPruningNumClusters;
  // beam width for cluster-based beam pruning.
  static double cliqueBeamClusterBeam;
  // state beam width for state-based cluster pruning
  static unsigned cliqueBeamClusterMaxNumStates;

  // TODO: @@@@ remove this variable, it is here just to gather a few stats for
  // Wed Dec 10 10:47:27 2008 Friday's talk. This will *NOT* work
  // with multiple cliques per partition.
  double prevFixedPrediction;


#define NO_PRUNING_CLIQUEBEAMCLUSTERPRUNINGMAXSTATES (UINT_MAX)

  // forced max number of states in a clique. Set to 0 to turn it off.
  static unsigned cliqueBeamMaxNumStates;
  // fraction of clique to retain, forcibly pruning away everything else. Must be

  // between 0 and 1 (i.e., 0 < v <= 1).
  static float cliqueBeamRetainFraction;
  // a version for clustered states
  static float cliqueBeamClusterRetainFraction;

  // between 0 and 1 (i.e., 0 < v <= 1).
  static double cliqueBeamMassRetainFraction;
  // exponentiate the clique scores before pruning by this amount.
  static double cliqueBeamMassExponentiate;
  // min possible resulting size.
  static unsigned cliqueBeamMassMinSize;
  // additional normal beam to include once mass is covered.
  static double cliqueBeamMassFurtherBeam;

  // a version of the above for, but for cluster/diversity pruning
  static double cliqueBeamClusterMassRetainFraction;
  static double cliqueBeamClusterMassExponentiate;
  static unsigned cliqueBeamClusterMassMinSize;
  static double cliqueBeamClusterMassFurtherBeam;


  // amount to re-sample from pruned clique. =0.0 to turn off.
  static double cliqueBeamUniformSampleAmount;


  // set to true to store all deterministic children that exist in
  // this clique with its parents in the clique sorage. Otherwise, set
  // to false so that the det. children are not stored (which saves
  // space and might but does not necessarily slow things down).
  static bool storeDeterministicChildrenInClique;


  // if this is turned on, GMTK will normalize the score of each clique so
  //  that we do not run out of numeric dynamic range. This will
  //  render prob(evidence) meaningless, but will give the same
  //  results for training or decoding. The special values
  //  f the variabale are:
  //      val == 1 ==> don't do any score normalize (default)
  //      val == 0 ==> normalize (divide) by the maximum clique value
  //      val > 1  ==> normalize (divide) by the given value. A good
  //                value to use woudl be exp(- log(prob(evidence)) / T ) 
  //                and the per-frame log-likelyhood can be estimaed
  //                by a short-length run of the program.
  static double normalizeScoreEachClique;


  // if true, zero cliques abort GMTK. if false, only the segment is
  // aborted and inference continues for the next segment
  static bool failOnZeroClique;


  // the set of nodes which form a max clique, in arbitrary order.
  set<RV*> nodes;

  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////

  // basic constructor with a set of nodes
  MaxCliqueBase(set<RV*> arg) {
    nodes = arg;
  }

  MaxCliqueBase() {}

  virtual ~MaxCliqueBase() {
    // TODO: do this right so that it works with tmp values in these objects.
    //   right now, not deleting these causes a memory leak.
    // if (maxCEValuePredictor != NULL) 
    // delete maxCEValuePredictor; 
  }



  // Complete the max clique by adjusting the random variables
  // themselves.
  virtual void makeComplete() { makeComplete(nodes); }
  // Static version of variable set completion, to complete set of
  // random variables passed in.
  static void makeComplete(const set<RV*> &rvs);



  // The penalty per feature that we pay to include a continous
  // observation in a clique. We include this since if we always
  // charge penalty of a factor of one, the continuous observation
  // might end up in a clique with many more than its own parents,
  // something that might possibly slow things down. The
  // default value is 0, meaning no additional penalty.
  static double continuousObservationPerFeaturePenalty;

  // Static version of various compute weight routines, useful in
  // certain places outside this class where we just have a collection
  // of nodes to compute the weight of.
  static float computeWeight(const set<RV*>& nodes,
			     const RV* node = NULL,
			     const bool useDeterminism = true);
  static float computeWeightWithExclusion(const set<RV*>& nodes,
					  const set<RV*>& unassignedIteratedNodes,
					  const set<RV*>& unionSepNodes,
					  const bool useDeterminism = true);
  static float computeWeightInJunctionTree(const set<RV*>& nodes,
					   const set<RV*>& assignedNodes,
					   const set<RV*>& cumulativeAssignedNodes,
					   const set<RV*>& unassignedIteratedNodes,
					   const set<RV*>& cumulativeUnassignedIteratedNodes,
					   const set<RV*>& separatorNodes,
					   const set<RV*>& unassignedInPartition,
					   set<RV*>* lp_nodes,
					   set<RV*>* rp_nodes,
					   const bool upperBound,
					   const bool moreConservative,
					   const bool useDeterminism);

  // compute the weight (log10 state space) of this clique.
  virtual float weight(const bool useDeterminism = true) const { 
    return computeWeight(nodes,NULL,useDeterminism); 
  }

  // print just the clique nodes
  virtual  void printCliqueNodes(FILE* f);
  
  // memory reporting
  virtual void reportMemoryUsageTo(FILE *f) = 0;
  

  // score reporting
  virtual void reportScoreStats();

  //////////////////////////////////////////////
  // TODO: figure out a way so that the member variables below exist only in
  // a subclass since all of the below is not needed for basic
  // triangulation.


  // TODO: write a destructor for this (not a problem for now since
  // we only create a fixed set of maxcliques per GMTK run). 


  virtual void numParentsSatisfyingChild(unsigned& num,unsigned par,vector <RV*> & parents, RV* child);



};


#endif
