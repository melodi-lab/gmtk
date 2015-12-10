/*-
 * GMTK_MaxClique.cc
 *
 *     maxClique support and JT probabilistic inference. Includes the
 *     implementation of a message in the message passing algorithm.
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (C) 2003 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 *
 */





/*
 * TODO:
 *   1) Maybe: In some of the structures the choice of when
 *   to use a hash table vs. an inlined value is done
 *   by either using an unsigned* or a unsigned. Create
 *   options where we do 'unsigned val[LEN]' and where
 *   only if packed value is > LEN words do we resort to
 *   the hash table.
 *
 *   2) Figure out way to remove the (unsigned**) type casts
 *   in front of call to packer object (see below).   
 *
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif
#if HAVE_HG_H
#  include "hgstamp.h"
#endif


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <float.h>
#include <assert.h>
#include <ctype.h>

#include <iterator>
#include <map>
#include <set>
#include <algorithm>

#include "general.h"
#include "error.h"
#include "debug.h"
#include "rand.h"
#include "psp.h"

#include "GMTK_FileParser.h"
#include "GMTK_RV.h"
#include "GMTK_DiscRV.h"
#include "GMTK_HidDiscRV.h"
#include "GMTK_ObsDiscRV.h"
#include "GMTK_ContRV.h"
#include "GMTK_GMTemplate.h"
#include "GMTK_GMParms.h"
#include "GMTK_MaxCliqueBase.h"
#include "GMTK_MaxClique.h"
#include "GMTK_PedagogicalCliqueTable.h"
#include "GMTK_MDCPT.h"
#include "GMTK_MSCPT.h"
#include "GMTK_MTCPT.h"
#include "GMTK_Mixture.h"
#include "GMTK_ObservationSource.h"
#include "GMTK_SectionScheduler.h"
#include "GMTK_ZeroCliqueException.h"

VCID(HGID)



////////////////////////////////////////////////////////////////////////
// Comment/Uncomment to optimize for speed/reducing memory usage.
#ifndef OPTIMIZE_FOR_MEMORY_USAGE
#  define OPTIMIZE_FOR_MEMORY_USAGE
#endif
////////////////////////////////////////////////////////////////////////


// for sorting an array of CliqueValue descending based on the contained logpr 
struct CliqueValueDescendingProbCompare
{  
  bool operator() (const PedagogicalCliqueTable::CliqueValue& cv1,
		   const PedagogicalCliqueTable::CliqueValue& cv2)
  {
    return (cv1.p > cv2.p);
  }
};





////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Static variables and functions
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

/*
 * integer value to keep track of indenting when running in trace mode.
 */
int PedagogicalCliqueTable::traceIndent = -1;

/*
 * number of spaces per indent. Could be a #define
 *
 */
const unsigned PedagogicalCliqueTable::spi = 1;





float PedagogicalCliqueTable::valuePoolGrowthRate;

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        PedagogicalCliqueTable support
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////




/*-
 *-----------------------------------------------------------------------
 * PedagogicalCliqueTable::PedagogicalCliqueTable()
 *    Clone constructor with frame delta to create a clone but under an unrolling.
 *    I.e., this isn't really a normal constructor, this is a contructor that
 *    sets up a clone of the MaxClique given by the argument from_clique. The
 *    clone is adjusted so that:
 *        1) it is shifted in time by frameDelta
 *        2) not all data structures are retained, only the ones 
 *           necessary to do exact inference.
 *
 * Preconditions:
 *
 * Postconditions:
 *
 * Side Effects:
 *
 * Results:
 *     nothing
 *
 *-----------------------------------------------------------------------
 */
PedagogicalCliqueTable::PedagogicalCliqueTable(MaxClique& origin)
{
  init(origin);
}

// get clique ready for use.
void PedagogicalCliqueTable::init(MaxClique& origin)
{
  numCliqueValuesUsed = 0;
#ifdef TRACK_NUM_CLIQUE_VALS_SHARED
  numCliqueValuesShared = 0;
#endif

  // TODO: optimize this and make depend on if clique is all hidden, has observed, etc.
  // NOTE: This must be set to something greater than 0.
  // cliqueValues.resize(3); // 10000
  cliqueValues.resize(origin.cliqueValueSpaceManager.currentSize());

}


/*-
 *-----------------------------------------------------------------------
 * PedagogicalCliqueTable::SharedLocalStructure::PedagogicalCliqueTable() constructor.
 *    INitialize the shared part of a PedagogicalCliqueTable
 *
 * Preconditions:
 *
 * Postconditions:
 *
 * Side Effects:
 *
 * Results:
 *     nothing
 *
 *-----------------------------------------------------------------------
 */
PedagogicalCliqueTable::SharedLocalStructure::
SharedLocalStructure(MaxClique& _origin,
		     vector <RV*>& newRvs,
		     map < RVInfo::rvParent, unsigned >& ppf,
		     const unsigned int frameDelta)
{

  origin = &_origin;

  set<RV*>::iterator it;

  rv_w_max_frame_num = rv_w_min_frame_num = NULL;

  // clone over nodes RVs.
  fNodes.resize(origin->nodes.size());
  unsigned i=0;
  for (it = origin->nodes.begin();
       it != origin->nodes.end();
       it++) {
    RV* rv = (*it);
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    

    if ( ppf.find(rvp) == ppf.end() ) {
      coredump("ERROR: can't find rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	       rv->name().c_str(),rv->frame(),frameDelta,
	       rvp.first.c_str(),rvp.second);
    }

    RV* nrv = newRvs[ppf[rvp]];
    fNodes[i++] = nrv;

    // compute and score max/min frame RVs.
    if (rv_w_max_frame_num == NULL || nrv->frame() > rv_w_max_frame_num->frame())
      rv_w_max_frame_num = nrv;
    if (rv_w_min_frame_num == NULL || nrv->frame() < rv_w_min_frame_num->frame())
      rv_w_min_frame_num = nrv;

  }

  // and clone over assigned nodes and sorted assigned nodes
  fSortedAssignedNodes.resize(origin->sortedAssignedNodes.size());
  for (i=0;i<origin->sortedAssignedNodes.size();i++) {
    RV* rv = origin->sortedAssignedNodes[i];
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    

    if ( ppf.find(rvp) == ppf.end() ) {
      coredump("ERROR: can't find assigned rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	       rv->name().c_str(),rv->frame(),frameDelta,
	       rvp.first.c_str(),rvp.second);
    }

    RV* nrv = newRvs[ppf[rvp]];
    fSortedAssignedNodes[i] = nrv;
  }

  // do unassignedIteratedNodes
  i=0;
  fUnassignedIteratedNodes.resize(origin->unassignedIteratedNodes.size());
  for (it = origin->unassignedIteratedNodes.begin();
       it != origin->unassignedIteratedNodes.end();
       it++) {
    RV* rv = (*it);
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    

    if ( ppf.find(rvp) == ppf.end() ) {
      coredump("ERROR: can't find assigned rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	       rv->name().c_str(),rv->frame(),frameDelta,
	       rvp.first.c_str(),rvp.second);
    }
    RV* nrv = newRvs[ppf[rvp]];
    fUnassignedIteratedNodes[i++] = nrv;
  }

  // Clique values only store/hash values of hidden (thus necessarily
  // discrete) variables since they are the only thing that change.
  discreteValuePtrs.resize(origin->hashableNodes.size());
  for (i=0;i<discreteValuePtrs.size();i++) {
    // get the unrolled rv for this hidden node
    RV* rv = origin->hashableNodes[i];
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    
    if ( ppf.find(rvp) == ppf.end() ) {
      warning("ERROR: can't find assigned rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	      rv->name().c_str(),rv->frame(),frameDelta,
	      rvp.first.c_str(),rvp.second);
      assert ( ppf.find(rvp) != ppf.end() );
    }
    RV* nrv = newRvs[ppf[rvp]];

    // hidden nodes are always discrete (in this version).
    DiscRV* drv = 
      (DiscRV*)nrv;
    // grab a pointer directly to its value for easy access later.
    discreteValuePtrs[i] = &(drv->val);
  }

  fDeterminableNodes.resize(origin->determinableNodes.size());
  for (i=0;i<fDeterminableNodes.size();i++) {
    // get the unrolled rv for this hidden node
    RV* rv = origin->determinableNodes[i];
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    
    if ( ppf.find(rvp) == ppf.end() ) {
      warning("ERROR: can't find assigned rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	      rv->name().c_str(),rv->frame(),frameDelta,
	      rvp.first.c_str(),rvp.second);
      assert ( ppf.find(rvp) != ppf.end() );
    }
    RV* nrv = newRvs[ppf[rvp]];

    // hidden nodes are always discrete (in this version).
    fDeterminableNodes[i] = RV2DRV(nrv);
  }

}


set <RV*> 
PedagogicalCliqueTable::SharedLocalStructure::returnRVsAndTheirObservedParentsAsSet()
{
  set<RV*> rc;
  for (unsigned i=0;i<fNodes.size();i++) {
    rc.insert(fNodes[i]);
    set <RV*> tmp = fNodes[i]->observedParents();
    // add in observed parents
    std::copy(tmp.begin(), tmp.end(),
	      inserter(rc,rc.end()));
  }
  return rc;
}

set <RV*> 
PedagogicalCliqueTable::SharedLocalStructure::returnRVsAsSet()
{
  set<RV*> rc;
  for (unsigned i=0;i<fNodes.size();i++) {
    rc.insert(fNodes[i]);
  }
  return rc;
}


vector <RV*> 
PedagogicalCliqueTable::SharedLocalStructure::returnRVsAsVector()
{
  vector<RV*> rc;
  for (unsigned i=0;i<fNodes.size();i++) {
    rc.push_back(fNodes[i]);
  }
  assert(rc.size() == fNodes.size());
  return rc;
}



//
// TODO: make proper comments to all of the functions below.
//

/*-
 *-----------------------------------------------------------------------
 * MaxClique::ceGatherFromIncommingSeparators()
 *
 *    Collect Evidence, Gather from Incomming Separators: This routine
 *    is the main driver for creating a clique table during the
 *    collect evidence stage of inference. It assumes that all
 *    separators have been created, based on that the clique will be created.
 *    There are two forms of clique table creation here, the one that
 *    is used is decided at the top of the routine.
 *    
 *        1) clique driven clique creation: Here, we iterate through
 *        all clique values of a cliuqe (based on the unassigned and
 *        assigned nodes in that clique), and for each one, check all
 *        incomming separators to make sure that the separators have
 *        entries for the corresponding intersected variables at the
 *        current clique value. The separator scores are then
 *        multiplied together with the scoring assigned variables,
 *        all of which becomes the clique value score. The clique
 *        value and its score are assigned to the clique.
 *
 *        2) separator driven clique instantiation: Here, the
 *        separators are iterated directly, and later separators are
 *        iterated over only those values which are "compatible" with
 *        earlier separators (meaning any zero probability entries are
 *        skipped entirely, this is done using the separator data
 *        structures which were set up prior to this routine). Once an
 *        entry survives each separator (meaning the separator
 *        intersection is non zero for that entry), it is subjected to
 *        assigned clique nodes, and if it survies them as well
 *        (meaning not zero prob), then the score is multiplied and
 *        added to the clique.
 *
 *        Note that this version first iterates through separators,
 *        then through unassigned clique nodes, and finally assigned
 *        clique nodes.
 *
 *    Inference method 1 is good because normal weight tells you
 *    exactly how computationally costly it is. Inference method 1 is
 *    bad because it is apparentlly *much* slower then method 2 when
 *    we have either many deterministic variables, or 2 when there is
 *    much pruning. Unfortunately, the computational cost of method 2
 *    is much more difficult to predict a-priori since it depends on
 *    pruning, and the degree of sparse/deterministic CPTs. Junction
 *    Tree weight (jtweight) is an attempt to predict this cost.
 *
 * Preconditions:
 *
 *   All the cliques incomming separators *must* have been created and are ready
 *   to be used to produce the clique table.
 *     
 *
 * Postconditions:
 *    Clique table has been created.
 *
 * Side Effects:
 *    Will significantly affect member variables in cliques.
 *
 * Results:
 *     nothing
 *
 *-----------------------------------------------------------------------
 */
void
PedagogicalCliqueTable::
ceGatherFromIncommingSeparators(PedagogicalCliqueTable::SharedLocalStructure& sharedStructure,
				ConditionalSeparatorTable* separatorTableArray,
				ConditionalSeparatorTable::SharedLocalStructure* sepSharedStructureArray)
{

  // syntactic convenience variables.
  MaxClique& origin = *(sharedStructure.origin);

  // we should never try more than 1x in this case.
  assert ( origin.cliqueBeamBuildBeam != (-LZERO) || origin.cliqueBeamBuildMaxExpansions == 1 );

  unsigned cliqueExpansionTry = 0;

  // Max collect-evidence probability for this clique. Used for beam
  // pruning.
  logpr maxCEValue;

  // An estimate of the beam threshold (i.e., est ~ maxCeValue -
  // threshold) for the current clique based on predictions from a
  // history of previous clique's maxCE values. If clique values fall
  // below this estimated threshold, they are pruned.
  logpr cliqueBeamThresholdEstimate;

  while (cliqueExpansionTry < origin.cliqueBeamBuildMaxExpansions) {

    traceIndent=-1; 
    // this is like the sub-main() for collect evidence.
    if (origin.cliqueBeamBuildBeam != (-LZERO)
	&& origin.maxCEValuePredictor.ptr() != NULL
	&& origin.maxCEValuePredictor->readyToMakePrediction()) {

      double currentCliqueBeamBuildBeam = 
	origin.cliqueBeamBuildBeam * (::pow(origin.cliqueBeamBuildExpansionFactor,cliqueExpansionTry));
      double maxCEValPrediction = origin.maxCEValuePredictor->makePrediction();
      double fixedPrediction = 2*origin.prevMaxCEValue.valref() - origin.prevPrevMaxCEValue.valref();

      cliqueBeamThresholdEstimate.valref() = maxCEValPrediction - 
	currentCliqueBeamBuildBeam;
      if (cliqueBeamThresholdEstimate.essentially_zero())
	cliqueBeamThresholdEstimate.set_to_almost_zero();
      // report the various values, 'p' = previous or predicted.

      double denom = origin.prevMaxCEValue.valref();
      // report absolute error if we can't compute relative error.
      if (denom == 0.0) denom = 1.0;
      infoMsg(IM::Inference,IM::Med,"Partial clique beam pruning, ppmax= %f, pmax= %f, ppred= %f, pfpred= %f, p_rel_%%err = %f, p_frel_%%err= %f, pred= %f, pthres= %f.\n",
	      origin.prevPrevMaxCEValue.valref(),
	      origin.prevMaxCEValue.valref(),
	      origin.prevMaxCEValPrediction,
	      origin.prevFixedPrediction,
	      100*::fabs(origin.prevMaxCEValPrediction - origin.prevMaxCEValue.valref())/
	      ::fabs(denom),
	      100*::fabs(origin.prevFixedPrediction - origin.prevMaxCEValue.valref())/
	      ::fabs(denom),
	      maxCEValPrediction,
	      cliqueBeamThresholdEstimate.valref());
      if (cliqueExpansionTry == 0) {
	origin.prevMaxCEValPrediction = maxCEValPrediction;
	origin.prevFixedPrediction = fixedPrediction;
      }
    } else {

      // always prune if we fall below or equal to almost zero.
      cliqueBeamThresholdEstimate.set_to_almost_zero(); 
      origin.prevFixedPrediction = 0;
      origin.prevMaxCEValPrediction = 0;
    }

    maxCEValue.set_to_zero();
    logpr p = 1.0;
    // next, do the actual collect message.
    if (origin.hashableNodes.size() == 0) {
      ceGatherFromIncommingSeparatorsCliqueObserved(sharedStructure,
						    separatorTableArray,
						    sepSharedStructureArray,
						    maxCEValue);
    } else if (true /* !origin.ceSeparatorDrivenInference */) { 
      ceGatherFromIncommingSeparatorsCliqueDriven(sharedStructure,
						  separatorTableArray,
						  sepSharedStructureArray,
						  cliqueBeamThresholdEstimate,
						  maxCEValue); // max value that is returned
    } else {
      // if we're still here, we do regular separator driven inference.
      if (origin.ceReceiveSeparators.size() == 0) {
	if (origin.unassignedIteratedNodes.size() == 0) {
	  ceIterateAssignedNodes(sharedStructure,
				 cliqueBeamThresholdEstimate,
				 maxCEValue, // max value that is returned
				 0, // initial node number
				 p);
	} else {
	  ceIterateUnassignedIteratedNodes(sharedStructure,
					   cliqueBeamThresholdEstimate,
					   maxCEValue,
					   0,
					   p);
	}
      } else {
	ceIterateSeparators(sharedStructure,
			    separatorTableArray,
			    sepSharedStructureArray,
			    cliqueBeamThresholdEstimate,
			    maxCEValue,
			    0,
			    p);
      }
    }

    if (numCliqueValuesUsed == 0) {
      // if we have a zero clique, print message and possibly continue
      // with expanded clique.
      if (message(IM::Inference, IM::Med)) {
	printf("WARNING: ZERO CLIQUE: clique with no entries, try %d out of %d.\n",
	       cliqueExpansionTry+1,origin.cliqueBeamBuildMaxExpansions);
      }
    } else {
      // current pruning level worked.
      break;
    }
    cliqueExpansionTry ++;
  }

  // check if we have a zero clique, and if we do, print message and exit.
  // TODO: rather than exit, pop back to the top and allow continuation and/or
  // beam expansion.
  if (numCliqueValuesUsed == 0) {
    if (MaxClique::failOnZeroClique) 
        error("ERROR: ZERO CLIQUE: clique with no entries. Final probability will be zero.\n");

    // It looks like there's no cleanup to do here - the loop above just breaks
    // without doing anything to cleanup the previous clique expansion tries. - RR
    warning("ZERO CLIQUE: clique with no entries. Final probability will be zero.\n");
    throw ZeroCliqueException();
  }

  // We have some clique entries, so we store new previous max CE
  // values, before any pruning.
  if (origin.cliqueBeamBuildBeam != (-LZERO)) {
    origin.prevPrevMaxCEValue.valref() = origin.prevMaxCEValue.valref();
    origin.prevMaxCEValue.valref() = maxCEValue.valref();
    if (origin.maxCEValuePredictor.ptr() != NULL) {
      origin.maxCEValuePredictor->addNextSampleAndUpdate(maxCEValue.valref());
    }
  }

  // In case that USE_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL is defined:
  // Note that if you want to print the entries here before pruning,
  // and if you want to call printing in the form of say:
  //
  //     printCliqueEntries(sharedStructure,stdout,NULL,false,false);
  // 
  // then you will need to change the printCliqueEntries routine to
  // not use the perment shared structures.


  // Now, we do pruning here. If the temporary clique value pool is
  // being used, we prune here, *before* we copy things out of the
  // temporary pool so that pruned entries are not inserted into
  // permanent locations.
  ceDoAllPruning(origin,maxCEValue);


#ifdef USE_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL
  // Now that pruning has occured, we copy things out of the pool,
  // into the permanent shared pool locaion.  
  // IMPORTANT: once this code is run, it means
  // that we cannot call ceDoAllPruning() any further
  //  (is this true? Indeed, if we prune again, the state
  //  won't be removed from the shared arrays, but is there
  //  an inherent harm in doing so? Perhaps not). 
  if (origin.packer.packedLen() > IMC_NWWOH) {
    // finally, insert surviving entries into global shared pool.
    insertLocalCliqueValuesIntoSharedPool(origin);
    // and free up the local buffer.
    origin.temporaryCliqueValuePool.resize(MaxClique::spaceMgrStartingSize*origin.packer.packedLen());
  }
#endif

  if (origin.normalizeScoreEachClique != 1.0)
    ceDoCliqueScoreNormalization(sharedStructure);

}


/*-
 *-----------------------------------------------------------------------
 * PedagogicalCliqueTable::ceGatherFromIncommingSeparatorsCliqueObserved()
 *
 *    A version of ceGatherFromIncommingSeparators that is is
 *    specifically for cliques that are all observed which means that
 *    all surrounding separators are also all observed.  observed
 *    clique observedclique
 *
 * Preconditions:
 *      Same as the separator driven case.
 *
 * Postconditions:
 *      Same as the separator driven case.
 *
 * Side Effects:
 *      Same as the separator driven case.
 *
 * Results:
 *     nothing
 *
 *-----------------------------------------------------------------------
 */

void
PedagogicalCliqueTable::
ceGatherFromIncommingSeparatorsCliqueObserved(PedagogicalCliqueTable::SharedLocalStructure& sharedStructure,
					      ConditionalSeparatorTable* separatorTableArray,
					      ConditionalSeparatorTable::SharedLocalStructure* sepSharedStructureArray,
					      logpr& maxCEValue)
{

  MaxClique& origin = *(sharedStructure.origin);

  // just do the assigned nodes as the unassigned observed nodes have no effect here.
  unsigned nodeNumber;

  logpr p = 1.0;
  bool zero = false;
  for (nodeNumber = 0; nodeNumber < sharedStructure.fSortedAssignedNodes.size(); nodeNumber ++ ) {
    RV* rv = sharedStructure.fSortedAssignedNodes[nodeNumber];
    if (origin.dispositionSortedAssignedNodes[nodeNumber] == MaxClique::AN_NOTSEP_PROB_SPARSEDENSE || 
	origin.dispositionSortedAssignedNodes[nodeNumber] == MaxClique::AN_NOTSEP_PROB_SPARSEDENSE) {
      // apply the probabiltiy
      logpr cur_p = rv->probGivenParents();

      if (message(Inference,Huge)) {
	psp2(stdout,spi*(traceIndent+1+nodeNumber));
	printf("%d:assigned obs/prob app, Pr[",nodeNumber);
	rv->printNameFrameValue(stdout,false);
	if (message(Inference, Mega)) {
	  if (rv->allParents.size() > 0) {
	    printf("|");
	    printRVSetAndValues(stdout,rv->allParents,false);
	  }
	} else {
	  if (rv->allParents.size() > 0) {
	    printf("|parents");
	  }
	}
	printf("]=%f,crClqPr=%f\n",cur_p.val(),p.val());
      }

      // if at any step, we get zero, then back out.
      if (cur_p.essentially_zero()) {
	// we've got a zero, so might as well stop here now. 
	zero = true;
	break;
      } else {
	p *= cur_p;
      }
    } else { 
      // in none of the other cases do we apply the probability. We
      // still check for zeros though.
      logpr cur_p = rv->probGivenParents();

      if (message(Inference,Huge)) {
	psp2(stdout,spi*(traceIndent+1+nodeNumber));
	printf("%d:assigned obs/zero rmv, Pr[",nodeNumber);
	rv->printNameFrameValue(stdout,false);
	if (message(Inference, Mega)) {
	  if (rv->allParents.size() > 0) {
	    printf("|");
	    printRVSetAndValues(stdout,rv->allParents,false);
	  }
	} else {
	  if (rv->allParents.size() > 0) {
	    printf("|parents");
	  } 
	}
	printf("]=%f,crClqPr=%f\n",cur_p.val(),p.val());
      }

      if (cur_p.essentially_zero()) {
	zero = true;
	// we've got a zero, so might as well stop here now. 
	break;
      }
    }
  }

  // make sure we have an entry.
  if (cliqueValues.size() == 0) {
    cliqueValues.resize(1);
  }
  // we will only be using one value.
  numCliqueValuesUsed = 1;  

  if (zero) {
    // don't bother checking separators.
    maxCEValue.set_to_zero();
    cliqueValues.ptr[0].p.set_to_zero();

  } else {

    // Now, we check all incoming CE separators, make sure the entry for
    // the current clique value it exists, and if it does, multiply by
    // separator probability checking for zeros. Since the clique
    // is all observed, we know that all separators are observed as well.

    for (unsigned sepNumber=0;sepNumber<origin.ceReceiveSeparators.size();sepNumber++) {
      // get a handy reference to the current separator table and its origin
      ConditionalSeparatorTable& sep = 
	separatorTableArray[origin.ceReceiveSeparators[sepNumber]];
      SeparatorClique& sepOrigin = 
	*(sepSharedStructureArray[origin.ceReceiveSeparators[sepNumber]].origin);

      // if the separator is currently not active, we should not be
      // using it for the additional probability here.
      if (sepOrigin.skipMe)
	continue;

      // keep a local variable copy of this around to avoid potential dereferencing.
      ConditionalSeparatorTable::AISeparatorValue * const
	sepSeparatorValuesPtr = sep.separatorValues->ptr; 


      unsigned accIndex = 0;
      // separator consists of all observed values. We just multiply
      // in the one entry the separator, that is guaranteed to be at
      // postiion 0,0. The key thing is that we must not do ANY hash
      // lookup as there are no hash tables in a separator that
      // has no hidden variables.

      ConditionalSeparatorTable::AISeparatorValue& sv
	= sepSeparatorValuesPtr[accIndex];

      // If anyone is almost zero, stop right now. 
      // Where was this allocated?  Search in file for key string
      // "ALLOCATE_REMVALUES_ALL_OBSERVED'
      if (sv.remValues.ptr[0].p.essentially_zero()) {
	maxCEValue.set_to_zero();
	cliqueValues.ptr[0].p.set_to_zero();
	return;
      }

      p *= sv.remValues.ptr[0].p;
    }

    // store max (and only) value.
    maxCEValue = p;

    // finally, save the probability
    cliqueValues.ptr[0].p = p;
  }

  if (message(Inference,High)) {
    // see https://j.ee.washington.edu/trac/gmtk/ticket/214#comment:14
    if (message(Inference,High+5))
      psp2(stdout,spi*(traceIndent+1+sharedStructure.fSortedAssignedNodes.size()));
    infoMsg(IM::Inference, IM::High,"CI:Inserting Observed %d-clique ent #0,pr=%f,sm=%f:",
	    sharedStructure.fNodes.size(),
	    cliqueValues.ptr[0].p.val(),sumProbabilities().val());
    printRVSetAndValues(stdout,sharedStructure.fNodes);
  }

}





void
PedagogicalCliqueTable::
ceGatherFromIncommingSeparatorsCliqueDriven(PedagogicalCliqueTable::SharedLocalStructure& sharedStructure,
					    ConditionalSeparatorTable* separatorTableArray,
					    ConditionalSeparatorTable::SharedLocalStructure* sepSharedStructureArray,
					    logpr cliqueBeamThresholdEstimate,
					    logpr& maxCEValue)
{
  logpr p = 1.0;
  // syntactic convenience variables.
  MaxClique& origin = *(sharedStructure.origin);
  if (origin.unassignedIteratedNodes.size() == 0) {
    ceIterateAssignedNodesCliqueDriven(sharedStructure, cliqueBeamThresholdEstimate, maxCEValue, 0, p);
  } else {
    //    ceIterateUnassignedNodesCliqueDriven(sharedStructure, ,0,p);
  }
}



/*
  A note on trace printing and debugging values to control it for the
  following few routines: 
  Rough range of debug values:

  Med  = 50,
  High = 60,
  Huge = 70,
  Mega = 80,
  Giga = 90,

  We have three routines for iterating a clique.
  ceIterateSeparators
  ceIterateUnassignedIteratedNodes
  ceIterateAssignedNodesRecurse
  - High: print just final clique insertions, nothign else. No indentation. 
  done just in ceIterateAssignedNodesRecurse
  - High+5: add all starts (Separator, Unassigned, and RV iteration starts)
  - Huge: Print all iters (separator, unassigned, & RV iters), but not parent values in RV case.
  - Mega: Also print all parent values at all iterations.
  - Mega+5: also prints continuous observation values (rather than just "=C").

*/


/*-
 *-----------------------------------------------------------------------
 * MaxClique::ceIterateSeparators()
 *
 *    Collect Evidence, Iterate Separators: This routine
 *    is part of separator driven clique instantiation. It 
 *    iterates through each separator and moves on to the next
 *    separator checking it based on the accumulated set of variables
 *    that have been assigned in previous separators. This routine
 *    might therefore be called a 'sparse join', since it does a join
 *    of a bunch of separator tables each of which might be very sparse. 
 *
 *    If we have reached the last separator in this clique, we move
 *    directly on to the unassigned nodes (if any). 
 *
 * Preconditions:
 *
 *   All the cliques incomming separators *must* have been created and
 *   are ready to be used to produce the clique table. 
 *     
 *
 * Postconditions:
 *    We move on to the unasssigned nodes with a probabilty value that
 *    includes the multiplication of all compatible separator entries for this var
 *    value set.
 *
 * Side Effects:
 *    none
 *
 * Results:
 *     nothing
 *
 *-----------------------------------------------------------------------
 */
void
PedagogicalCliqueTable::ceIterateSeparatorsRecurse(PedagogicalCliqueTable::SharedLocalStructure& sharedStructure,
					   ConditionalSeparatorTable* separatorTableArray,
					   ConditionalSeparatorTable::SharedLocalStructure* sepSharedStructureArray,
					   logpr cliqueBeamThresholdEstimate,
					   logpr& maxCEValue,
					   const unsigned sepNumber,
					   const logpr p)
{

  if (p.essentially_zero())
    return;

// TODO: at some point, also include some sort of continuation and pruning heuristic
// at the separator level.
//   if (p*separator_continuation[sepNum]*clique_continuation <= cliqueBeamThresholdEstimate) {
//     fprintf(stderr,"Returning from separator recurse, p.v = %f, thresest = %f\n",
// 	    p.val(),cliqueBeamThresholdEstimate.val());
//     return;
//   }
// or should have a separator threshold which is precomputed cliqueBeamThresholdEstimate-cliuqe_contin



  // syntactic simplicity and cached variables ...
  MaxClique& origin = *(sharedStructure.origin);
  ConditionalSeparatorTable& sep = 
    separatorTableArray[origin.ceReceiveSeparators[sepNumber]];
  ConditionalSeparatorTable::SharedLocalStructure& sepSharedStructure = 
    sepSharedStructureArray[origin.ceReceiveSeparators[sepNumber]];

  DiscRVType** accDiscreteValuePtrs = 
    sepSharedStructure.accDiscreteValuePtrs.ptr;
  DiscRVType** remDiscreteValuePtrs = 
    sepSharedStructure.remDiscreteValuePtrs.ptr;

  SeparatorClique& sepOrigin = 
    *(sepSharedStructure.origin);
  ConditionalSeparatorTable::AISeparatorValue * const
    sepSeparatorValuesPtr = sep.separatorValues->ptr; 


  traceIndent++;
  if (message(Inference,High+5)) {
    psp2(stdout,spi*traceIndent);    
    infoMsg(Inference,High+5,"S%d:Starting separator iter,partSepNo=%d,p=%f,nodes:",
	    sepNumber,origin.ceReceiveSeparators[sepNumber],p.val());
    printRVSet(stdout,sepSharedStructure.fNodes);
  }

  if (sepOrigin.skipMe) {
    // then we completely skip this separator passing onto the next
    // one.
    ceIterateSeparators(sharedStructure,
			separatorTableArray,
			sepSharedStructureArray,
			cliqueBeamThresholdEstimate,
			maxCEValue,
			sepNumber+1,
			p);
    goto ceIterateSeparatorsFinished;
  }


  unsigned sepValueNumber;  
  if (sepOrigin.hAccumulatedIntersection.size() > 0) {
    // look up existing intersected values to see if we have a match
    // and only proceed if we do.

    unsigned *key = &CliqueBuffer::packedVal[0];
    sepOrigin.accPacker.pack((unsigned**)accDiscreteValuePtrs,
			     (unsigned*)key);
    unsigned* indexp = sep.iAccHashMap->find(key);
    if (indexp == NULL) {
      // Then not found in this separator, so it must have (pruned) zero
      // probability. We continue with the next value of the previous
      // separator.
      if (message(Inference,Huge)) {
	psp2(stdout,spi*traceIndent);
	infoMsg(Inference,Huge,"S%d:Separator iter accumulated intersection prune\n",
		sepNumber);
	// TODO: @@@ figure out why we can't do: 
	// printRVSet(stdout,sep.fAccumulatedIntersection);
      }
      goto ceIterateSeparatorsFinished;
    } else {
      // need to further iterate.
      sepValueNumber = *indexp;
    }

  } else {
    // This condition would fail if we've completely pruned away
    // parent clique, but that can only happen when
    // we've got a negative beam. non-negative beam means
    // that we will always have at least one entry.
    // TODO: probably ok to remove this assertion.
    assert ( sep.separatorValues->size() == 1);
    sepValueNumber = 0;
  }

  // Iterate through remainder of separator. 
  // NOTE: we could do some online pruning here as well, but instead
  // we do it in a special separator prune routine, called ceSeparatorPrune().
  if (sepOrigin.hRemainder.size() == 0) {
    // Only one remainder entry (in position 0) and also no need to
    // unpack since all has been covered by accumulated intersection
    // set above in a previous separator. Just continue on with single
    // value. 
    // - 
    // Note that this case also includes the case when the entire
    // separator is observed, i.e., when we'll get an accum inter size
    // of 0 and a hreminder size of 0 -- nothing to unpack here either
    // (no hash tables even exist), so we just continue along.

    if (message(Inference, Huge)) {
      psp2(stdout,spi*traceIndent);
      infoMsg(Inference, Huge,"S%d:Separator iter no-unpack %d,%d,partSepNo=%d,p=%f,sp=%f,nodes:",
	      sepNumber,
	      sepSeparatorValuesPtr[sepValueNumber].remValues.size(),
	      sepSeparatorValuesPtr[sepValueNumber].numRemValuesUsed,
	      origin.ceReceiveSeparators[sepNumber],
	      p.val(),
	      sepSeparatorValuesPtr[sepValueNumber].remValues.ptr[0].p.val());
      printRVSetAndValues(stdout,sepSharedStructure.fNodes);
    }

    // We should have either one or zero entries. The only way zero
    // entries could arrise is if we have done severe separator pruning.
    assert ( (sepSeparatorValuesPtr[sepValueNumber].remValues.size() & ~0x1)
	     == 0x0 );

    // assert ( sepSeparatorValuesPtr[sepValueNumber].remValues.size() == 1 );
    if (sepSeparatorValuesPtr[sepValueNumber].numRemValuesUsed == 1) {
      // Continue down with new probability value.
      // Search for tag 'ALLOCATE_REMVALUES_OPTION' in this file for
      // more info why remValues.ptr[0] exists.
      // Note: could do more separator pruning here.
      ceIterateSeparators(sharedStructure,
			  separatorTableArray,
			  sepSharedStructureArray,
			  cliqueBeamThresholdEstimate,
			  maxCEValue,
			  sepNumber+1,
			  p*
			  sepSeparatorValuesPtr[sepValueNumber].remValues.ptr[0].p);
    }
  } else {

    // TODO: this assertion should be redundant (check above)
    assert ( sepOrigin.remPacker.packedLen() > 0 );

    // TODO: perhaps special case for VE seps, since all probs are == 1, so no need to multiply.

    if (sepOrigin.remPacker.packedLen() <= ISC_NWWOH_RM) {
      for (unsigned i=0;i< sepSeparatorValuesPtr[sepValueNumber].numRemValuesUsed; i++) {

	// TODO: optimize this, pre-compute base array outside of loop.
	sepOrigin.remPacker.unpack(
				   (unsigned*)&(sepSeparatorValuesPtr[sepValueNumber].remValues.ptr[i].val[0]),
				   (unsigned**)remDiscreteValuePtrs);


	if (message(Inference, Huge)) {
	  psp2(stdout,spi*traceIndent);
	  infoMsg(Inference, Huge,"S%d:Separator iter %d of %d,partSepNo=%d,p=%f,sp=%f,nodes:",
		  sepNumber,
		  i,sepSeparatorValuesPtr[sepValueNumber].numRemValuesUsed,
		  origin.ceReceiveSeparators[sepNumber],
		  p.val(),
		  sepSeparatorValuesPtr[sepValueNumber].remValues.ptr[i].p.val());
	  printRVSetAndValues(stdout,sepSharedStructure.fNodes);
	}

	// continue down with new probability value.
	ceIterateSeparators(sharedStructure,
			    separatorTableArray,
			    sepSharedStructureArray,			    
			    cliqueBeamThresholdEstimate,
			    maxCEValue,
			    sepNumber+1,
			    p*sepSeparatorValuesPtr[sepValueNumber].remValues.ptr[i].p);
      }
    } else {
      for (unsigned i=0;i< sepSeparatorValuesPtr[sepValueNumber].numRemValuesUsed; i++) {

	// TODO: optimize this, pre-compute base array outside of loop.
	sepOrigin.remPacker.unpack(
				   (unsigned*)sepSeparatorValuesPtr[sepValueNumber].remValues.ptr[i].ptr,
				   (unsigned**)remDiscreteValuePtrs);

	if (message(Inference, Huge)) {
	  psp2(stdout,spi*traceIndent);
	  infoMsg(Inference, Huge,"S%d:Separator iter %d of %d,partSepNo=%d,p=%f,sp=%f,nodes:",
		  sepNumber,
		  i,sepSeparatorValuesPtr[sepValueNumber].numRemValuesUsed,
		  origin.ceReceiveSeparators[sepNumber],
		  p.val(),
		  sepSeparatorValuesPtr[sepValueNumber].remValues.ptr[i].p.val());
	  printRVSetAndValues(stdout,sepSharedStructure.fNodes);
	}

	// continue down with new probability value.
	ceIterateSeparators(sharedStructure,
			    separatorTableArray,
			    sepSharedStructureArray,
			    cliqueBeamThresholdEstimate,
			    maxCEValue,
			    sepNumber+1,
			    p*sepSeparatorValuesPtr[sepValueNumber].remValues.ptr[i].p);
      }
    }

  }

 ceIterateSeparatorsFinished:
  //  if (message(Inference, High+5))
  traceIndent--;
}


/*-
 *-----------------------------------------------------------------------
 * MaxClique::ceIterateUnassignedIteratedNodes()
 *
 *    Collect Evidence, Iterate Separators: This routine
 *    is part of separator driven clique instantiation. Once
 *    we have a partial clique value that has survied the separators,
 *    the next thing we need to do is iterate through any unassigned
 *    nodes in this cliuqe (if any, hopefully not since these can be costly
 *    and might indicate a poor triangulation, depending on pruning and/or
 *    the sparse/deterministic variables).
 *
 *    If we have reached the last unassigned node in this clique, we
 *    move directly on to the assigned nodes.
 *
 * Preconditions:
 *
 *   All the cliques incomming separators *must* have been created and are ready
 *   to be used to produce the clique table.
 *     
 *
 * Postconditions:
 *    We move on to the asssigned nodes with an probabilty. The prob. has
 *    not been updated from what is passed in since these nodes are unassigned
 *    and so never contribute score to this clique.
 *
 * Side Effects:
 *    none
 *
 * Results:
 *     nothing
 *
 * TODO:
 *     should not force iterating unassigned at the top of the recurssion since it might 
 *     no be relevant or useful for the clique expansion (e.g., switching parents might
 *     mean that this variable isn't needed for most cases).
 *
 *-----------------------------------------------------------------------
 */

void
PedagogicalCliqueTable::ceIterateUnassignedIteratedNodesRecurse(PedagogicalCliqueTable::SharedLocalStructure& sharedStructure,
							logpr cliqueBeamThresholdEstimate,
							logpr& maxCEValue,
							const unsigned nodeNumber,
							const logpr p)
{

  RV* rv = sharedStructure.fUnassignedIteratedNodes[nodeNumber];
  // TODO: update comments here to match others.
  traceIndent++;
  if (message(Inference, High+5)) {
    psp2(stdout,spi*traceIndent);
    infoMsg(Inference, High+5,"U%d:Starting Unassigned iteration of rv %s(%d),p=%f\n",
	    nodeNumber,
	    rv->name().c_str(),rv->frame(),p.val());
  }
  if (rv->hidden()) {
    // only discrete RVs can be hidden for now.
    DiscRV* drv = (DiscRV*)rv;
    // do the loop right here
    drv->val = 0;
    do {
      if (message(Inference, Huge)) {
	psp2(stdout,spi*traceIndent);
	infoMsg(Inference, Huge,"U%d:Unassigned iter of rv %s(%d)=%d,p=%f\n",
		nodeNumber,
		rv->name().c_str(),rv->frame(),drv->val,p.val());
      }
      // continue on, effectively multiplying p by unity.
      ceIterateUnassignedIteratedNodes(sharedStructure,
				       cliqueBeamThresholdEstimate,
				       maxCEValue,
				       nodeNumber+1,p);
    } while (++drv->val < drv->cardinality);
  } else {
    // TODO: Perhaps unassignedIteratedNodes should contain
    // no observed nodes at all since we are not updating 
    // probability here anyway.

    if (message(Inference, Huge)) {
      psp2(stdout,spi*traceIndent);
      // observed, either discrete or continuous
      if (rv->discrete()) {
	infoMsg(Inference, Huge,"U%d:Unassigned pass through observed rv %s(%d)=%d,p=%f\n",
		nodeNumber,
		rv->name().c_str(),rv->frame(),RV2DRV(rv)->val,p.val());
      } else {
	// nothing to do since we get continuous observed value
	// indirectly
	infoMsg(Inference, Huge,"U%d:Unassigned pass through of observed rv %s(%d)=C,p=%f\n",
		nodeNumber,
		rv->name().c_str(),rv->frame(),p.val());
      }
    }
    // continue on, effectively multiplying p by unity.
    ceIterateUnassignedIteratedNodes(sharedStructure,
				     cliqueBeamThresholdEstimate,
				     maxCEValue,
				     nodeNumber+1,p);
  }
  //  if (message(Inference, High+5))
  traceIndent--;
}


/*-
 *-----------------------------------------------------------------------
 * MaxClique::ceIterateAssignedNodesRecurse()
 *
 *    Collect Evidence, Iterate Separators: This routine
 *    is part of separator driven clique instantiation. 
 *  
 *    This routine is a main workhorse of inference. Once we have a
 *    partial clique value that has survied the separators, and has
 *    (possibly) unassigned variable values, the next thing we need to
 *    do is iterate through the assigned nodes in this cliuqe.
 *
 *
 *    If we have reached the last unassigned node in this clique, we
 *    go ahead and assign the value to this clique if it is non-zero.
 *
 *    If we have not reached the last, we 'iterate' over the assigned
 *    node depending on its type (see MaxClique.h for a definition of
 *    the different types of assigned nodes).
 *    
 *    Note: this routine will *never* knowingly add a zero-probability
 *    clique value to the clique. I.e., all current zeros are not
 *    added. This way, during backward pass if anyway, we never need
 *    to check for divide by zero.
 *
 *
 * Preconditions:
 *
 *   All the cliques incomming separators *must* have been iterated, and
 *   same for the unassigned nodes (if any). The clique is assumed
 *   to have at least one hidden node in it.
 *
 * Postconditions:
 *    We move on to the asssigned nodes with an probabilty. The prob. has
 *    not been updated from what is passed in since these nodes are unassigned
 *    and so never contribute score to this clique.
 *
 * Side Effects:
 *    none
 *
 * Results:
 *     nothing
 *
 *-----------------------------------------------------------------------
 */
void
PedagogicalCliqueTable::ceIterateAssignedNodesRecurse(PedagogicalCliqueTable::SharedLocalStructure& sharedStructure,
					      logpr cliqueBeamThresholdEstimate,
					      logpr& maxCEValue,
					      const unsigned nodeNumber,
					      const logpr p)
{

  // syntactic convenience variables.
  MaxClique& origin = *(sharedStructure.origin);

  // First, do (potentially) partial clique beam pruning here.
  // 
  // Note, we may wish, if (nodeNumber == fSortedAssignedNodes.size())
  // to just go ahead and add the clique entry since it is already
  // here, but we still prune since the very last variable that was
  // assigned might have given us a value that falls below threshold.
  if (p*origin.sortedAssignedContinuationScores[nodeNumber] <= cliqueBeamThresholdEstimate)
    return;

  // No pruning, so we go ahead with the continued expansion.

  if (nodeNumber == sharedStructure.fSortedAssignedNodes.size()) {
    // time to store clique value and total probability, p is
    // current clique probability.

    // printf("ceIterateAssignedNodesRecurse: nodeNumber = %d, p = %f,",nodeNumber,p.val());
    // printRVSet(fNodes);

    // keep track of the max clique probability right here.
    if (p > maxCEValue)
      maxCEValue = p;

    if (numCliqueValuesUsed >= cliqueValues.size()) {
      // TODO: optimize this.
      if (numCliqueValuesUsed >= origin.cliqueValueSpaceManager.currentSize())
	origin.cliqueValueSpaceManager.advanceToNextSize();
      // cliqueValues.resizeAndCopy(cliqueValues.size()*2);
      cliqueValues.resizeAndCopy(origin.cliqueValueSpaceManager.currentSize());
    }

    // TODO: figure out if it is possible to get around doing this
    // check (or to help branch prediction predict it, since it will
    // be different for differnet cliques). Answer: We can do this
    // check once at beginning of iteration of assigned nodes, and
    // have two versions of this code.
    if (origin.packer.packedLen() <= IMC_NWWOH) {
      // pack the clique values directly into place
      origin.packer.pack(
			 (unsigned**)sharedStructure.discreteValuePtrs.ptr,
			 (unsigned*)&(cliqueValues.ptr[numCliqueValuesUsed].val[0]));
    } else {

#ifdef USE_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL      
      const unsigned long lindex = numCliqueValuesUsed*origin.packer.packedLen();
      if (lindex >= origin.temporaryCliqueValuePool.size()) {
	// use aggressive growth factor for now to avoid expensive copies.
	origin.temporaryCliqueValuePool.resizeAndCopy(
						      origin.packer.packedLen()*
						      int(1.5+(double)origin.temporaryCliqueValuePool.size()*valuePoolGrowthRate));
      }
      unsigned *pcv = 
	&origin.temporaryCliqueValuePool.ptr[lindex];
      origin.packer.pack((unsigned**)sharedStructure.discreteValuePtrs.ptr,(unsigned*)pcv);
      // store integer value of the location.
      cliqueValues.ptr[numCliqueValuesUsed].ival = lindex;
#else

      // Deal with the hash table to re-use clique values.
      // First, grab pointer to storge where next clique value would
      // be stored if it ends up being used.
      unsigned *pcv = origin.valueHolder.curCliqueValuePtr();
      // Next, pack the clique values into this position.
      origin.packer.pack((unsigned**)sharedStructure.discreteValuePtrs.ptr,(unsigned*)pcv);
      // Look it up in the hash table.
      bool foundp;
      unsigned *key;
      key = origin.cliqueValueHashSet.insert(pcv,foundp);
      if (!foundp) {
	// if it was not found, need to claim this storage that we
	// just used.
	origin.valueHolder.allocateCurCliqueValue();
      }
      // Save the pointer to whatever the hash table decided to use.
      cliqueValues.ptr[numCliqueValuesUsed].ptr = key;
#endif


    }
    // save the probability
    cliqueValues.ptr[numCliqueValuesUsed].p = p;
    numCliqueValuesUsed++;

    if (message(Inference, High)) {
      // see https://j.ee.washington.edu/trac/gmtk/ticket/214#comment:14
      if (message(Inference, High+5))
	psp2(stdout,spi*(traceIndent+1));
      infoMsg(Inference, High,"CI:Inserting %d-clique ent #%d,pr=%f,sm=%f:",
	      sharedStructure.fNodes.size(),
	      (numCliqueValuesUsed-1),
	      cliqueValues.ptr[numCliqueValuesUsed-1].p.val(),sumProbabilities().val());
      printRVSetAndValues(stdout,sharedStructure.fNodes);
    }
    return;
  }
  RV* rv = sharedStructure.fSortedAssignedNodes[nodeNumber];
  // do the loop right here

  traceIndent++;
  if (message(Inference, High+5)) {
    psp2(stdout,spi*traceIndent);
    infoMsg(Inference, High+5,"A%d:Starting assigned iteration of rv %s(%d),crClqPr=%f\n",
	    nodeNumber,
	    rv->name().c_str(),rv->frame(),p.val());
  }

  switch (origin.dispositionSortedAssignedNodes[nodeNumber]) {
  case MaxClique::AN_CPT_ITERATION_COMPUTE_AND_APPLY_PROB: 
    {
      logpr cur_p;
      rv->begin(cur_p);
      do {
	if (message(Inference, Huge)) {
	  psp2(stdout,spi*traceIndent);
	  printf("A%d:assigned iter/prob app, Pr[",nodeNumber);
	  rv->printNameFrameValue(stdout,false);
	  if (message(Inference, Mega)) {
	    if (rv->allParents.size() > 0) {
	      printf("|");
	      printRVSetAndValues(stdout,rv->allParents,false);
	    }
	  } else {
	    if (rv->allParents.size() > 0) {
	      printf("|parents");
	    }
	  }
	  printf("]=%f,crClqPr=%f\n",cur_p.val(),p.val());
	}
	// if at any step, we get zero, then back out.
	if (!cur_p.essentially_zero()) {
	  // Continue, updating probability by cur_p, contributing this
	  // probability to the clique potential.
	  ceIterateAssignedNodesRecurse(sharedStructure,
					cliqueBeamThresholdEstimate,maxCEValue,
					nodeNumber+1,p*cur_p);
	}
	// TODO: backjump
      } while (rv->next(cur_p));
    }
    break;

  case MaxClique::AN_CPT_ITERATION_COMPUTE_PROB_REMOVE_ZEROS:
    {
      logpr cur_p;
      rv->begin(cur_p);
      do {
	// At each step, we compute probability
	if (message(Inference, Huge)) {
	  psp2(stdout,spi*traceIndent);
	  printf("A%d:assigned iter/zero rmv, Pr[",nodeNumber);
	  rv->printNameFrameValue(stdout,false);
	  if (message(Inference, Mega)) {
	    if (rv->allParents.size() > 0) {
	      printf("|");
	      printRVSetAndValues(stdout,rv->allParents,false);
	    }
	  } else {
	    if (rv->allParents.size() > 0) {
	      printf("|parents");
	    } 
	  }
	  printf("]=%f,crClqPr=%f\n",cur_p.val(),p.val());
	}
	// if at any step, we get zero, then back out.
	if (!cur_p.essentially_zero()) {
	  // Continue, do not update probability!!
	  ceIterateAssignedNodesRecurse(sharedStructure,
					cliqueBeamThresholdEstimate,maxCEValue,
					nodeNumber+1,p);
	}
      } while (rv->next(cur_p));
    }
    break;


  case MaxClique::AN_CARD_ITERATION:
    {
      DiscRV* drv = (DiscRV*)rv;
      // do the loop right here
      drv->val = 0;
      do {
	if (message(Inference, Huge)) {
	  psp2(stdout,spi*traceIndent);
	  printf("A%d:assigned card iter, Pr[",nodeNumber);
	  rv->printNameFrameValue(stdout,false);
	  if (message(Inference, Mega)) {
	    if (rv->allParents.size() > 0) {
	      printf("|");
	      printRVSetAndValues(stdout,rv->allParents,false);
	    }
	  } else {
	    if (rv->allParents.size() > 0) {
	      printf("|parents");
	    }
	  }
	  printf("]=???,crClqPr=%f\n",p.val());
	}
	// Continue, do not update probability!!
	ceIterateAssignedNodesRecurse(sharedStructure,
				      cliqueBeamThresholdEstimate,maxCEValue,
				      nodeNumber+1,p);
      } while (++drv->val < drv->cardinality);
    }
    break;

  case MaxClique::AN_COMPUTE_AND_APPLY_PROB:
    {
      // TODO: Make more efficient version of this, based on the type of
      // RV.
      logpr cur_p = rv->probGivenParents();
      // if at any step, we get zero, then back out.
      if (message(Inference, Huge)) {
	psp2(stdout,spi*traceIndent);
	printf("A%d:assigned compute appl prob, Pr[",nodeNumber);
	rv->printNameFrameValue(stdout,false);
	if (message(Inference, Mega)) {
	  if (rv->allParents.size() > 0) {
	    printf("|");
	    printRVSetAndValues(stdout,rv->allParents,false);
	  }
	} else {
	  if (rv->allParents.size() > 0) {
	    printf("|parents");
	  }
	}
	printf("]=%f,crClqPr=%f\n",cur_p.val(),p.val());
      }
      if (!cur_p.essentially_zero()) {
	// Continue, updating probability by cur_p.
	ceIterateAssignedNodesRecurse(sharedStructure,
				      cliqueBeamThresholdEstimate,maxCEValue,				      
				      nodeNumber+1,p*cur_p);
      }
    }
    break;

  case MaxClique::AN_CONTINUE:
    if (message(Inference, Huge)) {
      psp2(stdout,spi*traceIndent);
      printf("A%d:sep cont, non prob, Pr[",nodeNumber);
      rv->printNameFrameValue(stdout,false);
      if (message(Inference, Mega)) {
	if (rv->allParents.size() > 0) {
	  printf("|");
	  printRVSetAndValues(stdout,rv->allParents,false);
	}
      } else {
	if (rv->allParents.size() > 0) {
	  printf("|parents");
	}
      }
      printf("]=???,crClqPr=%f\n",p.val());
    }
    ceIterateAssignedNodesRecurse(sharedStructure,
				  cliqueBeamThresholdEstimate,maxCEValue,
				  nodeNumber+1,p);
    break;

  case MaxClique::AN_CONTINUE_COMPUTE_PROB_REMOVE_ZEROS:
    {
      // TODO: Make more efficient version of this, based on the type of
      // RV.
      logpr cur_p = rv->probGivenParents();
      if (message(Inference, Huge)) {
	psp2(stdout,spi*traceIndent);
	printf("A%d:assigned compute continue, Pr[",nodeNumber);
	rv->printNameFrameValue(stdout,false);
	if (message(Inference, Mega)) {
	  if (rv->allParents.size() > 0) {
	    printf("|");
	    printRVSetAndValues(stdout,rv->allParents,false);
	  }
	} else {
	  if (rv->allParents.size() > 0) {
	    printf("|parents");
	  }
	}
	printf("]=%f,crClqPr=%f\n",cur_p.val(),p.val());
      }
      if (!cur_p.essentially_zero()) {
	// Continue, do not update probability!!
	ceIterateAssignedNodesRecurse(sharedStructure,
				      cliqueBeamThresholdEstimate,
				      maxCEValue,
				      nodeNumber+1,p);
      }
    }
    break;

  default:
    assert(0);
    break;
  }
  //  if (message(Inference, High+5))
  traceIndent--;

}






void 
PedagogicalCliqueTable::ceIterateAssignedNodesCliqueDriven(PedagogicalCliqueTable::SharedLocalStructure& sharedStructure,
							   logpr cliqueBeamThresholdEstimate,
							   logpr& maxCEValue,
							   const unsigned nodeNumber,
							   const logpr p)
{
}




/*-
 *-----------------------------------------------------------------------
 * MaxClique::ceIterateAssignedNodesNoRecurse()
 *
 *    This is a non-recursive version of
 *    ceIterateAssignedNodesRecurse(). The reason for this routine is
 *    to try to produce code that has better branch prediction
 *    behavior and that also creates fewer temporary variables on
 *    routine calls. Timings show that GMTK is about 10-15% faster
 *    using this routine.
 *
 *    This routine is a bit more involved than the recursive version
 *    above, as it uses customize loops (with much use of the evil
 *    goto statements, so the code is a bit ugly). It is essentially
 *    assembly code but expressed in C. Note that the recursive
 *    version is called as a backup when there are no assigned nodes
 *    in a clique or when verbose printing is turned on. For that
 *    reason, and also for pedagogical reasons, the recursive version
 *    is still in place. Any modifications not associated with verbose
 *    printing will need to be done in both places.
 *
 *
 *  TODO: when we find that all values of a variable cause a prune, don't back up to the previous variable (which 
 *        might not be a parent),
 *        rather back up to the latest parent of that variable and change its value. We might even
 *        back out to the separators.
 *
 * Preconditions:
 *
 *   Same as ceIterateAssignedNodesRecurse()
 *
 * Postconditions:
 *
 *   Same as ceIterateAssignedNodesRecurse()
 *
 * Side Effects:
 *
 *   Same as ceIterateAssignedNodesRecurse()
 *
 * Results:
 *
 *   Same as ceIterateAssignedNodesRecurse()
 *
 *-----------------------------------------------------------------------
 */
void
PedagogicalCliqueTable::ceIterateAssignedNodesNoRecurse(PedagogicalCliqueTable::SharedLocalStructure& sharedStructure,
						logpr cliqueBeamThresholdEstimate,
						logpr& maxCEValue,
						const logpr p)
{
  // syntactic convenience variables.
  MaxClique& origin = *(sharedStructure.origin);

  if (p*origin.sortedAssignedContinuationScores[0]  <= cliqueBeamThresholdEstimate)
    return;


  // parray has to be 1 offset, storing p in entry -1
  logpr* parray = origin.probArrayStorage.ptr + 1;
  parray[-1] = p;

  // Another option: compute a local clique beam estimate that can be used
  //    to compare directly against cur_p.
  // set as follows:
  //    cliqueBeamThresholdEstimate = cliqueBeamThresholdEstimate/p;
  //    cliqueBeamThresholdEstimate = cliqueBeamThresholdEstimate/parray[nodeNumber];
  //
  // or could subtract off continuation scores from an array estimate.

  int nodeNumber;
  logpr cur_p;
  for (nodeNumber=0;nodeNumber<(int)sharedStructure.fSortedAssignedNodes.size();nodeNumber++) {
    RV* rv = sharedStructure.fSortedAssignedNodes.ptr[nodeNumber];
    switch (origin.dispositionSortedAssignedNodes.ptr[nodeNumber]) {

    case MaxClique::AN_CPT_ITERATION_COMPUTE_AND_APPLY_PROB: 
      rv->begin(cur_p);
      goto applyProbTag;
    case MaxClique::AN_COMPUTE_AND_APPLY_PROB:
      rv->probGivenParents(cur_p);
      {
      applyProbTag:
	// check for possible zero (could occur with
	// zero score or observations).
	parray[nodeNumber] = parray[nodeNumber-1]*cur_p;
	if (parray[nodeNumber]*origin.sortedAssignedContinuationScores[nodeNumber] 
	    <= cliqueBeamThresholdEstimate) {
	  // Since we have small here, we cancel iterations of all
	  // subsequent clique variables right now, rather than
	  // iterate them with what will end up being below threshold
	  // probability (assuming the Gaussians are never > 1).
	  parray[nodeNumber].set_to_zero(); // @@@ REMOVE
	  nodeNumber++;
	  goto end_of_initial_clique_value;
	}
      }
      break;

    case MaxClique::AN_CPT_ITERATION_COMPUTE_PROB_REMOVE_ZEROS:
      rv->begin(cur_p);
      goto removeZeroTag;
    case MaxClique::AN_CONTINUE_COMPUTE_PROB_REMOVE_ZEROS:
      rv->probGivenParents(cur_p);
      {
      removeZeroTag:
	// check for possible zero (could occur with
	// zero score or observations).
	// TODO: could do cpbeam pruning here and see if 
	//       parray[nodeNumber-1]*cur_p falls below threshold pruning if it does.
	if (cur_p.essentially_zero()) {
	  // Since we have zero here, we cancel iterations of all
	  // subsequent variables right now, rather than iterate them
	  // with what will end up being zero probability.
	  parray[nodeNumber].set_to_zero();
	  nodeNumber++;
	  goto end_of_initial_clique_value;
	} else 
	  parray[nodeNumber] = parray[nodeNumber-1];
      }
      break;

    case MaxClique::AN_CARD_ITERATION:
      {
	DiscRV* drv = (DiscRV*)rv;
	drv->val = 0;
	parray[nodeNumber] = parray[nodeNumber-1];
      }
      break;

    case MaxClique::AN_CONTINUE:
      parray[nodeNumber] = parray[nodeNumber-1];
      break;

    }
  }

 end_of_initial_clique_value:

  // get ready for main loop
  const bool imc_nwwoh_p = (origin.packer.packedLen() <= IMC_NWWOH);
  const int maxNodeNumber = (int)sharedStructure.fSortedAssignedNodes.size()-1;
  nodeNumber--;

  MaxClique::AssignedNodeDisposition cur_disp = origin.dispositionSortedAssignedNodes.ptr[nodeNumber];
  RV* cur_rv = sharedStructure.fSortedAssignedNodes.ptr[nodeNumber];

  // check if zero probability, and if so, skip the first one and
  // continue on.
  if (parray[nodeNumber]*origin.sortedAssignedContinuationScores[nodeNumber]  
      <= cliqueBeamThresholdEstimate)
    goto next_iteration;

  // main loop, iterate through all assigned nodes in this clique.
  do {

    // add a clique value to the clique.
    {
      const logpr final_p = parray[nodeNumber];

      // time to store clique value and total probability, p is current
      // clique probability.
      // keep track of the max clique probability right here.

      if (final_p > maxCEValue)
	maxCEValue = final_p;

      if (numCliqueValuesUsed >= cliqueValues.size()) {
	// TODO: optimize this.
	// cliqueValues.resizeAndCopy(cliqueValues.size()*2);
	if (numCliqueValuesUsed >= origin.cliqueValueSpaceManager.currentSize())
	  origin.cliqueValueSpaceManager.advanceToNextSize();
	cliqueValues.resizeAndCopy(origin.cliqueValueSpaceManager.currentSize());
      }

      // Possibly remove this check if possible. It's probably ok though
      // since as long as this loop runs several times, branch
      // predictions should handle it (so 1 cycle latency).
      if (imc_nwwoh_p) {
	// pack the clique values directly into place
	origin.packer.pack(
			   (unsigned**)sharedStructure.discreteValuePtrs.ptr,
			   (unsigned*)&(cliqueValues.ptr[numCliqueValuesUsed].val[0]));
      } else {


#ifdef USE_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL
	const unsigned long lindex = numCliqueValuesUsed*origin.packer.packedLen();
	if (lindex >= origin.temporaryCliqueValuePool.size()) {
	  // use aggressive growth factor for now to avoid expensive copies.
	  origin.temporaryCliqueValuePool.resizeAndCopy(
							origin.packer.packedLen()*
							int(1.5+(double)origin.temporaryCliqueValuePool.size()*valuePoolGrowthRate));
	}
	unsigned *pcv = 
	  &origin.temporaryCliqueValuePool.ptr[lindex];
	origin.packer.pack((unsigned**)sharedStructure.discreteValuePtrs.ptr,(unsigned*)pcv);
	// store integer value of the location.
	cliqueValues.ptr[numCliqueValuesUsed].ival = lindex;

#if 0
	// grab pointer to storage.
	unsigned *pcv = origin.localValueHolder.curCliqueValuePtr();
	origin.packer.pack((unsigned**)discreteValuePtrs.ptr,(unsigned*)pcv);
	// claim stored value.
	origin.localValueHolder.allocateCurCliqueValue();
	// store pointer to appropriate location.
	cliqueValues.ptr[numCliqueValuesUsed].ptr = pcv;
#endif

#else

	// Deal with the hash table to re-use clique values.
	// First, grab pointer to storge where next clique value would
	// be stored if it ends up being used.
	unsigned *pcv = origin.valueHolder.curCliqueValuePtr();
	// Next, pack the clique values into this position.
	origin.packer.pack((unsigned**)sharedStructure.discreteValuePtrs.ptr,(unsigned*)pcv);
	// Look it up in the hash table.
	bool foundp;
	unsigned *key;
	key = origin.cliqueValueHashSet.insert(pcv,foundp);
	if (!foundp) {
	  // if it was not found, need to claim this storage that we
	  // just used.
	  origin.valueHolder.allocateCurCliqueValue();
	} 
#ifdef TRACK_NUM_CLIQUE_VALS_SHARED	
	else
	  numCliqueValuesShared++;
#endif
	// Save the pointer to whatever the hash table decided to use.
	cliqueValues.ptr[numCliqueValuesUsed].ptr = key;


#endif

      }
      // save the probability
      cliqueValues.ptr[numCliqueValuesUsed].p = final_p;
      numCliqueValuesUsed++;

      /*
       * uncomment this next code to produce lots of messages.
       * note above, if high verbosity is on, recursive version of this
       * routine will print this information.
       */
      /*
	if (message(Inference, Mega)) {
	// psp2(stdout,spi*traceIndent);
	infoMsg(Inference, Mega,"Inserting New Clique Val,pr=%f,sm=%f: ",
	cliqueValues.ptr[numCliqueValuesUsed-1].p.val(),sumProbabilities().val());
	printRVSetAndValues(stdout,fNodes);
	}
      */

    }

  next_iteration:
    // Now we need to move to next clique value. This bunch of next
    // code is like an inlined iterator over clique values. It is done
    // as "inline" in an attempt to avoid branch mis-predicts and since
    // this code appears only one time.
    do {

      bool unfinished;
      switch (cur_disp) {
      case MaxClique::AN_CPT_ITERATION_COMPUTE_AND_APPLY_PROB: 
	{
	  do {
	    // continue going until we get one that is unfinished and within beam.
	    unfinished = cur_rv->next(cur_p);
	    if (unfinished) {
	      REGISTER logpr tmp = parray[nodeNumber-1]*cur_p;
	      if (tmp*origin.sortedAssignedContinuationScores[nodeNumber]  
		  > cliqueBeamThresholdEstimate) {
		parray[nodeNumber] = tmp;
		goto next_node_number;
	      }
	    } else
	      break;
	  } while (1);
	}
	break;

      case MaxClique::AN_CPT_ITERATION_COMPUTE_PROB_REMOVE_ZEROS:
	{
	  // guaranteed not to get zero prob here. Note that
	  // we get the probability into cur_p, but we do not
	  // use it in this case. Also, parray[nodeNumber] is
	  // already up to date from the "begin()" at the beginning
	  // of this variables iteration.
	  unfinished = cur_rv->next(cur_p);
	  // we could include:
	  // assert ( !cur_p.essentially_zero() );
	}
	break;

      case MaxClique::AN_CARD_ITERATION:
	{
	  DiscRV* drv = (DiscRV*)cur_rv;
	  unfinished = (++drv->val < drv->cardinality);
	}
	break;

      default:
	// these cases are handeled by 'default':
	// case MaxClique::AN_COMPUTE_AND_APPLY_PROB:
	// case MaxClique::AN_CONTINUE:
	// case MaxClique::AN_CONTINUE_COMPUTE_PROB_REMOVE_ZEROS:
	{
	  unfinished = false;
	}
	break;
      }

      if (unfinished) {
	// Best case, we continue on to next outer iteration filling
	// in next clique value.
	break;
      } else {
	// we're finished with the current variable.
	if (nodeNumber == 0) {
	  // then we're really done since we're finished with the
	  // first clique variable.
	  return;
	} else {
	  // we're finished with the current variable, need to
	  // continue on with previous variable.
	  nodeNumber--;
	  cur_disp = origin.dispositionSortedAssignedNodes.ptr[nodeNumber];
	  cur_rv = sharedStructure.fSortedAssignedNodes.ptr[nodeNumber];
	}
      }

      // still here? Continue on incrementing previous variable.
    } while (1);

  next_node_number:
    while (nodeNumber < maxNodeNumber) {
      nodeNumber++;
      cur_disp = origin.dispositionSortedAssignedNodes.ptr[nodeNumber];
      cur_rv = sharedStructure.fSortedAssignedNodes.ptr[nodeNumber];

      // need to fill up the rest of the table.
      switch (cur_disp) {
      case MaxClique::AN_CPT_ITERATION_COMPUTE_AND_APPLY_PROB: 
	{
	  cur_rv->begin(cur_p);
	  // might get a zero probability, check condition here.
	  // TODO: keep result in temp and only write back out if need be. @@@
	  REGISTER logpr tmp =  parray[nodeNumber-1]*cur_p;
	  if (tmp*origin.sortedAssignedContinuationScores[nodeNumber] 
	      <= cliqueBeamThresholdEstimate) {
	    // We just did a begin and got zero on the first try. 
	    // we need to continue on with this variable until we finish.
	    goto next_iteration;
	  } else {
	    // only write when we have to.
	    parray[nodeNumber] = tmp;
	  }
	}
	break;

      case MaxClique::AN_COMPUTE_AND_APPLY_PROB:
	{
	  cur_rv->probGivenParents(cur_p);
	  // might get a zero probability, check condition here.
	  REGISTER logpr tmp = parray[nodeNumber-1]*cur_p;
	  if (tmp*origin.sortedAssignedContinuationScores[nodeNumber] 
	      <= cliqueBeamThresholdEstimate) {
	    // Since we have zero here, we cancel iterations of all
	    // subsequent variables right now, rather than iterate
	    // them with what will end up being zero probability.
	    nodeNumber--;
	    cur_disp = origin.dispositionSortedAssignedNodes.ptr[nodeNumber];
	    cur_rv = sharedStructure.fSortedAssignedNodes.ptr[nodeNumber];
	    goto next_iteration;
	  } else
	    parray[nodeNumber] = tmp;
	}
	break;

      case MaxClique::AN_CPT_ITERATION_COMPUTE_PROB_REMOVE_ZEROS:
	{
	  cur_rv->begin(cur_p);
	  // might get a zero probability, check condition here.
	  if (cur_p.essentially_zero()) {
	    goto next_iteration;
	  } else {
	    // TODO: update post-condition comments in RV and CPT iterators as
	    // to when a zero RV can occur and when not.
	    parray[nodeNumber] = parray[nodeNumber-1];
	  }
	}
	break;

      case MaxClique::AN_CONTINUE_COMPUTE_PROB_REMOVE_ZEROS:
	{
	  cur_rv->probGivenParents(cur_p);
	  // might get a zero probability, check condition here.
	  if (cur_p.essentially_zero()) {
	    // Since we have zero here, we cancel iterations of all
	    // subsequent variables right now, rather than iterate them
	    // with what will end up being zero probability.
	    nodeNumber--;
	    cur_disp = origin.dispositionSortedAssignedNodes.ptr[nodeNumber];
	    cur_rv = sharedStructure.fSortedAssignedNodes.ptr[nodeNumber];
	    goto next_iteration;
	  } else {
	    // TODO: update post-condition comments in RV and CPT iterators as
	    // to when a zero RV can occur and when not.
	    parray[nodeNumber] = parray[nodeNumber-1];
	  }
	}
	break;


      case MaxClique::AN_CARD_ITERATION:
	{
	  DiscRV* drv = (DiscRV*)cur_rv;
	  drv->val = 0;
	  parray[nodeNumber] = parray[nodeNumber-1];
	}
	break;

      case MaxClique::AN_CONTINUE:
	parray[nodeNumber] = parray[nodeNumber-1];
	break;

      }
    }

  } while (1);

}




/*-
 *-----------------------------------------------------------------------
 * PedagogicalCliqueTable::clearCliqueAndIncommingSeparatorMemory()
 *
 *    memory clearing routine, this routine clears all significant memory
 *    associated with this clique and all its incomming separators.
 *    It should only be called when using this clique in a collect-evidence
 *    form, where the clique memory is never going to be used again (such as
 *    during prob(Evidence) form of inference, where the goal is just to
 *    compute the probability of evidence.
 *
 *
 * Preconditions:
 *
 *   The basic data structures should be set up.
 *
 * Postconditions:
 *    The memory has been completely freed. It should be possible to reconstruct the
 *    separator again though. Also, if the 'alsoClearOrigins=true', then
 *    we can no longer use an instance of this clique again without constructing
 *    it from scratch (i.e., it would always do a hash insert).
 *    
 *
 * Side Effects:
 *    Changes 
 *      1) all of the memory associated with this clique and its incomming separators
 *      2) if 'alsoClearOrigins' is true, it will delete all of the memory associated with
 *         the origin of the clique.
 *
 * Results:
 *     nothing
 *
 *-----------------------------------------------------------------------
 */
void 
PedagogicalCliqueTable::clearCliqueAndIncommingSeparatorMemory(PedagogicalCliqueTable::SharedLocalStructure& sharedStructure,
						       ConditionalSeparatorTable* separatorTableArray,
						       ConditionalSeparatorTable::SharedLocalStructure* sepSharedStructureArray)
{
  MaxClique& origin = *(sharedStructure.origin);
  // first do the separators
  for (unsigned sepNumber=0;sepNumber<origin.ceReceiveSeparators.size();sepNumber++) {
    ConditionalSeparatorTable& sep = 
      separatorTableArray[origin.ceReceiveSeparators[sepNumber]];
    sep.clearInferenceMemory();
  }
  // and next do self.
  clearInferenceMemory();
}


/*-
 *-----------------------------------------------------------------------
 * MaxClique::ceSendToOutgoingSeparator()
 *
 *    Collect Evidence, Send to outgoing separator.
 *
 *    We have now a fully instantiated clique. This routine Iterates
 *    through the values that are above beam in the clique table, and
 *    instantiates the outgoing separator with those values. If a separator
 *    value has multiple clique values, then values are accumulated into
 *    the separator (in a Veterbi approach, we would take the max, see
 *    veterbi code). 
 *
 *    This code also combine clique pruning right here, rather than needing
 *    to do a separate pruning stage by calling ceCliqueBeamPrune().
 * 
 *
 * Preconditions:
 *
 *   The clique table must be fully instantiated, and contain entries
 *   for this clique.  The outgoing separator must have been created,
 *   setup, and ready to accept a message (projection down) from a
 *   clique
 *
 * Postconditions:
 *    The outgoing separator has been instantiated.
 *
 * Side Effects:
 *    Changes 
 *      1) the clique table, since with pruning, it will potentially shrink the clique table
 *      2) the outgoing separator table, since it is being created.
 *
 * Results:
 *     nothing
 *
 *-----------------------------------------------------------------------
 */
void 
PedagogicalCliqueTable::
ceSendToOutgoingSeparator(PedagogicalCliqueTable::SharedLocalStructure& sharedStructure,
			  ConditionalSeparatorTable* separatorTableArray,
			  ConditionalSeparatorTable::SharedLocalStructure* sepSharedStructureArray)
{
  // syntactic convenience variables.
  MaxClique& origin = *(sharedStructure.origin);
  ceSendToOutgoingSeparator(sharedStructure,
			    separatorTableArray[origin.ceSendSeparator],
			    sepSharedStructureArray[origin.ceSendSeparator]
			    );
}
void 
PedagogicalCliqueTable::
ceSendToOutgoingSeparator(PedagogicalCliqueTable::SharedLocalStructure& sharedStructure,
			  ConditionalSeparatorTable& sep,
			  ConditionalSeparatorTable::SharedLocalStructure& sepSharedStructure)
{

  // printf("At start of sending to outgoing separator\n");

  // keep a local variable copy of this around to avoid potential
  // dereferencing.  This one cannot be const since it might change
  // during a resize, in which case we need to reassign this variable.
  ConditionalSeparatorTable::AISeparatorValue * 
    sepSeparatorValuesPtr = sep.separatorValues->ptr; 

#ifdef TRACK_NUM_CLIQUE_VALS_SHARED  
  infoMsg(IM::Inference, IM::High-2,"ceSendToOutgoingSep: from clique w state space = %d. NumShared = %d, %2.2f percent\n",
	  numCliqueValuesUsed,numCliqueValuesShared, 100*(float)numCliqueValuesShared/(float)numCliqueValuesUsed);
#else
  infoMsg(IM::Inference, IM::High-2,"ceSendToOutgoingSep: from clique w state space = %d.\n",numCliqueValuesUsed);
#endif

  // syntactic convenience variables.
  MaxClique& origin = *(sharedStructure.origin);
  SeparatorClique& sepOrigin = 
    *(sepSharedStructure.origin);  

  // first check if this is an all "observed" clique
  if (origin.hashableNodes.size() == 0) {
    // Everything in this clique is observed or deterministic.  Therefore, there should
    // be one and only one clique value in this clique which is the
    // probability of the assigned probability nodes in this clique.

    // Since this clique is a superset of the separator, it means also
    // that the seperator will consist only of observed nodes. We
    // therefore create one accumulated entry and one rem entry and
    // insert our current probability.

    // printf("=====>>>> Observed clique in forward pass\n");

    // first do some sanity checks.
    assert (sepOrigin.hAccumulatedIntersection.size() == 0);
    assert (sepOrigin.hRemainder.size() == 0);
    assert (sepSharedStructure.remDiscreteValuePtrs.size() == 0);

    const unsigned accIndex = 0;

    // We here keep handy reference for readability.
    // To find out where has this been allocated, 
    // search for tag: ALLOCATE_REMVALUES_ALL_OBSERVED.
    ConditionalSeparatorTable::AISeparatorValue& sv
      = sepSeparatorValuesPtr[accIndex];

    // This must be first time for this entry.  We never need more
    // than one rem value.  Search for tag
    // 'ALLOCATE_REMVALUES_ALL_OBSERVED' in this file for where else
    // this could be done, but we allocate it here.
    if (sv.remValues.size() == 0)
      sv.remValues.resize(1); 
    // in all cases
    sv.numRemValuesUsed = 1;	  

    sv.remValues.ptr[0].p = cliqueValues.ptr[0].p;

    // and we're done already. This was easy!
    return;
  }

  // Note, we could optionally do pruning here rather than where we do
  // it now, namely when we construct the clique. There may be some
  // utility in doing it here (delaying the pruning) since the values
  // at the end of the pruned clique table could be used again before
  // they are deallocated.

  // TODO: do sampling code here, i.e., optionally sample from
  //       remaining portion of what otherwise would be pruned away
  //       portion of the clique.

  // Now, we send the clique to the outgoing separator.

  // next check if the outgoing separator has only obseved values.
  if (sepOrigin.hAccumulatedIntersection.size() == 0 
      && sepOrigin.hRemainder.size() == 0) {

    // printf("Observed separator\n");

    // Then indeed, outgoing separator is all observed values. We do
    // this special separately from the general case since that case
    // is already getting a bit unwieldy.

    // TODO: this can also probably handle the case of disconnected networks.
    // 

    ConditionalSeparatorTable::AISeparatorValue& sv
      = sepSeparatorValuesPtr[0];

    // This must be first time for this entry.  We never need more
    // than one rem value.  Search for tag
    // 'ALLOCATE_REMVALUES_ALL_OBSERVED' in this file for where else
    // this could be done, but we allocate it here.
    if (sv.remValues.size() == 0)
      sv.remValues.resize(1); 
    // in all cases, we use only one entry when the separator is observed.
    sv.numRemValuesUsed = 1;

    // Just sum up the clique entries projecting down into the
    // observed separator, but only do the ones that pass the beam
    // threshold. Arguably, we might not want or need to do beam
    // pruning here since when the separator is all observed, any
    // sparseness that we introduce by clique pruning isn't going to
    // change things later on (since everything is projected to the
    // same point), but we do it here anyway for numerical consistency
    // with the general case.
    for (unsigned cvn=0;cvn<numCliqueValuesUsed;) {
      if (SectionScheduler::viterbiScore) {
	// sv.remValues.ptr[0].p.assign_if_greater(cliqueValues.ptr[cvn].p);
	// TODO: add k-best
	if (cliqueValues.ptr[cvn].p > sv.remValues.ptr[0].p) {
	  sv.remValues.ptr[0].p = cliqueValues.ptr[cvn].p;
	  sv.remValues.ptr[0].backPointer = cvn;
	}
      } else {
	sv.remValues.ptr[0].p += cliqueValues.ptr[cvn].p;
      }
      cvn++;
    }

  } else {
    // We are guaranteed that either we have an accumulated
    // intersection, a reminder, or both (but not neither).  Go
    // through clique values and accumulate into appropriate place
    // within separator.

    // pre-compute a few values that the compiler might spend more
    // time than necessary to re-check.
    const bool imc_nwwoh_p = (origin.packer.packedLen() <= IMC_NWWOH);
    const bool isc_nwwoh_ai_p = (sepOrigin.accPacker.packedLen() <= ISC_NWWOH_AI);
    const bool isc_nwwoh_rm_p = (sepOrigin.remPacker.packedLen() <= ISC_NWWOH_RM);
    const bool sep_origin_hAccumulatedIntersection_exists_p =
      (sepOrigin.hAccumulatedIntersection.size() > 0);
    const bool sep_remDiscreteValuePtrs_exists_p = 
      (sepSharedStructure.remDiscreteValuePtrs.size() > 0);

    for (unsigned cvn=0;cvn<numCliqueValuesUsed;) {{

	// printf("Iteration through clique, iter = %d\n",cvn);

	// TODO: optimize away this conditional check. (and/or use const
	// local variable to indicate it wont change)
	if (imc_nwwoh_p) {
	  origin.packer.unpack((unsigned*)&(cliqueValues.ptr[cvn].val[0]),
			       (unsigned**)sharedStructure.discreteValuePtrs.ptr);
	} else {
	  origin.packer.unpack((unsigned*)cliqueValues.ptr[cvn].ptr,
			       (unsigned**)sharedStructure.discreteValuePtrs.ptr);
	}
	for (unsigned j=0;j<sharedStructure.fDeterminableNodes.size();j++) {
	  RV* rv = sharedStructure.fDeterminableNodes.ptr[j];
	  RV2DRV(rv)->assignDeterministicChild();
	}


	// All hidden random variables now have their discrete
	// value. Accumulate this probability into the given separator.

	/*
	 * There are 3 cases.
	 * 1) Both AI exists and REM exist.
	 * 2) AI exists but REM doesnt exist.
	 * 3) AI does not exist, but REM exists.
	 * x) AI not exist and REM not exist can't occur (covered above).
	 */

	unsigned accIndex;
	// TODO: optimize this check away out of loop.
	if (sep_origin_hAccumulatedIntersection_exists_p) { 
	  // an accumulated intersection exists.

	  // make sure there is at least one available accumulated intersection entry
	  assert ( sep.numSeparatorValuesUsed <= sep.separatorValues->size());
	  if (sep.numSeparatorValuesUsed >= sep.separatorValues->size()) {
	    const unsigned old_size = sep.separatorValues->size();
	    // TODO: optimize this size re-allocation.
	    if (sep.numSeparatorValuesUsed >= sepOrigin.separatorValueSpaceManager.currentSize()) 
	      sepOrigin.separatorValueSpaceManager.advanceToNextSize();
	    sep.separatorValues->resizeAndCopy(sepOrigin.separatorValueSpaceManager.currentSize()); 
	    sepSeparatorValuesPtr = sep.separatorValues->ptr;
	    if (isc_nwwoh_ai_p) {
	      // Then the above resize just invalided all our pointers to keys,
	      // but it did not invalidate the array indices. Go through
	      // and correct the keys within the hash table.
	      // TODO: think of a better way to do this that also looses no efficiency.
	      for (unsigned i=0;i<sep.iAccHashMap->tableSize();i++) {
		if (!sep.iAccHashMap->tableEmpty(i)) {
		  sep.iAccHashMap->tableKey(i)
		    = &(sepSeparatorValuesPtr[sep.iAccHashMap->tableItem(i)].val[0]);
		}
	      }
	    }
	    const unsigned new_size = sep.separatorValues->size();
	    // if (sep.remDiscreteValuePtrs.size() > 0) {
	    if (sep_remDiscreteValuePtrs_exists_p) {
	      for (unsigned i=old_size;i<new_size;i++) {
		// re-construct hash tables only for new entries.
		new (&sepSeparatorValuesPtr[i].iRemHashMap)
		  VHashMapUnsignedUnsignedKeyUpdatable
		  (sepOrigin.remPacker.packedLen(),ConditionalSeparatorTable::remHashMapStartingSize);
		// TODO: potentially preallocate default size of  
		// separatorValues->ptr[i].remValues.resize(default);
		// TODO: potentially create zero size here, and only
		//       grow bigger when we start adding things.
	      }
	    }
	  }

	  unsigned *accKey;
	  // TODO: optimize this check out of loop.
	  if (isc_nwwoh_ai_p) {
	    accKey = &(sepSeparatorValuesPtr[sep.numSeparatorValuesUsed].val[0]);
	    sepOrigin.accPacker.pack((unsigned**)sepSharedStructure.accDiscreteValuePtrs.ptr,
				     accKey);
	  } else {
	    accKey = sepOrigin.accValueHolder.curCliqueValuePtr();
	    sepOrigin.accPacker.pack((unsigned**)sepSharedStructure.accDiscreteValuePtrs.ptr,
				     accKey);
	    // check if this value combination already lives in
	    // origin's value holder hash table and if so, use that.
	    bool foundp;
	    accKey = sepOrigin.accSepValHashSet.insert(accKey,foundp);
	    if (!foundp) {
	      // only allocate a new value if it was inserted.
	      sepOrigin.accValueHolder.allocateCurCliqueValue();
	    }
	    // store the pointer in case we use it.
	    sepSeparatorValuesPtr[sep.numSeparatorValuesUsed].ptr = accKey;
	  }

	  bool foundp;
	  unsigned* accIndexp =
	    sep.iAccHashMap->insert(accKey,
				    sep.numSeparatorValuesUsed,
				    foundp);

	  if (!foundp) {
	    //  add the values we just used. 
	    sep.numSeparatorValuesUsed++;
	  }
	  accIndex = *accIndexp;

	  // TODO: optimize this check out of loop.
	  // if (sep.remDiscreteValuePtrs.size() == 0) {
	  if (!sep_remDiscreteValuePtrs_exists_p) {
	    // 2) AI exists and REM doesnt exist
	    // Then this separator is entirely covered by one or 
	    // more other separators earlier in the order.

	    // go ahead and insert it here to the 1st entry (entry 0).

	    // handy reference for readability.
	    ConditionalSeparatorTable::AISeparatorValue& sv
	      = sepSeparatorValuesPtr[accIndex];

	    // Accumulate the clique's
	    // probability into this separator's probability.
	    if (sv.remValues.size() < 1) {
	      // This must be first time for this entry.
	      // Search for tag 'ALLOCATE_REMVALUES_OPTION' in this file for where else this could
	      // be done.
	      sv.remValues.resize(1);
	      sv.numRemValuesUsed = 1;	  
	      // initialize and assign.
	      sv.remValues.ptr[0].p = cliqueValues.ptr[cvn].p;
	      if (SectionScheduler::viterbiScore)
		sv.remValues.ptr[0].backPointer = cvn;
	    } else {
	      // already there so must have hit before.
	      // we thus accumulate.
	      if (SectionScheduler::viterbiScore) {
		// sv.remValues.ptr[0].p.assign_if_greater(cliqueValues.ptr[cvn].p);
		if (cliqueValues.ptr[cvn].p > sv.remValues.ptr[0].p) {
		  sv.remValues.ptr[0].p = cliqueValues.ptr[cvn].p;
		  sv.remValues.ptr[0].backPointer = cvn;
		}
	      } else {
		sv.remValues.ptr[0].p += cliqueValues.ptr[cvn].p;
	      }
	    }

	    goto next_iteration;
	  }

	} else {
	  accIndex = 0;
	}

	// If we're here, then we are guaranteed must have some remainder
	// pointers, i.e., we could do:
	//    assert (sep.remDiscreteValuePtrs.size() > 0);
	// So, the remainder exists in this separator.
	// either:
	//   1) AI exists and REM exist
	//     or
	//   3) AI does not exist (accIndex == 0), but REM exists
	// 

	// keep handy reference for readability.
	ConditionalSeparatorTable::AISeparatorValue& sv
	  = sepSeparatorValuesPtr[accIndex];

	// make sure there is at least one available entry
	assert (sv.numRemValuesUsed <= sv.remValues.size());
	if (sv.numRemValuesUsed >= sv.remValues.size()) {
	  // TODO: optimize this growth rate.
	  // start small but grow fast.
	  // sv.remValues.resizeAndCopy(1+sv.remValues.size()*2); // *3
	  sv.remValues.resizeAndCopy(sepOrigin.remainderValueSpaceManager.nextSizeFrom(sv.remValues.size()));
	  sepOrigin.remainderValueSpaceManager.setCurrentAllocationSizeIfLarger(sv.remValues.size());
	  if (isc_nwwoh_rm_p) {
	    // Then the above resize just invalided all sv.iRemHashMap's pointers to keys,
	    // but it did not invalidate its array indices. Go through
	    // and correct the keys within the hash table.
	    // TODO: think of a better way to do this that looses no efficiency.
	    for (unsigned i=0;i<sv.iRemHashMap.tableSize();i++) {
	      if (!sv.iRemHashMap.tableEmpty(i)) {
		sv.iRemHashMap.tableKey(i)
		  = &(sv.remValues.ptr[sv.iRemHashMap.tableItem(i)].val[0]);
	      }
	    }
	  }
	}

	unsigned *remKey;
	// pack relevant variable values
	// TODO: optimize away this check.
	if (isc_nwwoh_rm_p) {
	  // grab pointer to next location to be used in this case.
	  remKey = &(sv.remValues.ptr[sv.numRemValuesUsed].val[0]);
	  // pack the remainder pointers
	  sepOrigin.remPacker.pack((unsigned**)sepSharedStructure.remDiscreteValuePtrs.ptr,
				   remKey);
	} else {
	  // grab pointer to next packed clique value to be used.
	  remKey = sepOrigin.remValueHolder.curCliqueValuePtr();
	  sepOrigin.remPacker.pack((unsigned**)sepSharedStructure.remDiscreteValuePtrs.ptr,
				   remKey);
	  // check if this value combination already lives in
	  // origin's value holder hash table and if so, use that.
	  bool foundp;
	  remKey = sepOrigin.remSepValHashSet.insert(remKey,foundp);
	  if (!foundp) {
	    // only allocate a new value if it was inserted.
	    sepOrigin.remValueHolder.allocateCurCliqueValue();
	  }
	  // store the pointer in case we use it.
	  sv.remValues.ptr[sv.numRemValuesUsed].ptr = remKey;
	}

	bool foundp;
	unsigned* remIndexp =
	  sv.iRemHashMap.insert(remKey,
				sv.numRemValuesUsed,
				foundp);
	if (!foundp) {
	  // add the values we just used. 
	  sv.numRemValuesUsed++;
	}

	// We've finally got the entry, so accumulate the clique's
	// probability into this separator's probability.
	if (SectionScheduler::viterbiScore) {
	  // sv.remValues.ptr[*remIndexp].p.assign_if_greater(cliqueValues.ptr[cvn].p);
	  if (cliqueValues.ptr[cvn].p > sv.remValues.ptr[*remIndexp].p) {
	    sv.remValues.ptr[*remIndexp].p = cliqueValues.ptr[cvn].p;
	    sv.remValues.ptr[*remIndexp].backPointer = cvn;
	  }
	} else {
	  sv.remValues.ptr[*remIndexp].p += cliqueValues.ptr[cvn].p;
	}

	// printf("Inserted sep value, iter = %d\n",cvn);

      }
    next_iteration:
      cvn++;
    }
  }

  // To reallocate or not to reallocate, that is the question.  here,
  // we just reallocate for now. The reason why we might want to reallocate is that
  // pruning might have adjusted numCliqueValuesUsed in this routine above.
  // 
  // TODO: reallocate only if change is > some percentage (say 5%),
  // and export to command line.
  // TODO: possibly tell origin.cliqueValueSpaceManager about this pruning event.
  // TODO: if -probE option or gmtkDecode is running, no need to re-allocate here since this
  //       will soon be deleted anyway.
  // TODO: if we do a backward pass, we might not want to resize, and if an entry ends up
  //       becoming probable, re-insert the entry as a valid clique entry.
  if (numCliqueValuesUsed < cliqueValues.size())
    cliqueValues.resizeAndCopy(numCliqueValuesUsed);

  // printCliqueEntries(stdout);

  /////////////////////////////////////
  // And prune the separator as well.
  sep.ceSeparatorPrune(sepOrigin);

}


/*-
 *-----------------------------------------------------------------------
 * PedagogicalCliqueTable::ceCliqueBeamPrune()
 *
 *    Collect Evidence, Clique Prune: This routine will prune away
 *    part of a previously instantiated clique based on the current
 *    clique beam width.
 *    
 *    Note that PedagogicalCliqueTable::ceSendToOutgoingSeparator() does
 *    its own pruning, so when using ceSendToOutgoingSeparator(), this
 *    pruning routine does not need to be called (at least with the
 *    same beam width). This routine, however, is kept here since 1)
 *    at the very right clique of the right most partition, we need to
 *    explicitely prune, and 2) the island algorithm sometimes also
 *    needs to explicitly call pruning.
 *
 * Preconditions:
 *   1) the value of the max clique 'maxCEValue' must have been
 *      computed already.
 *
 *   2) clique table must be created, meaning that either:
 *
 *        PedagogicalCliqueTable::ceIterateAssignedNodesRecurse()
 *
 * Postconditions:
 *    Clique table has been pruned, and memory for it has NOT been re-allocated to
 *    fit the smaller size.
 *
 * Side Effects:
 *    changes the clique size
 *
 * Results:
 *     nothing
 *
 *-----------------------------------------------------------------------
 */

void 
PedagogicalCliqueTable::ceCliqueBeamPrune(MaxClique& origin,
				  logpr maxCEValue)
{
  // return immediately if beam pruning is turned off.
  if (origin.cliqueBeam == (-LZERO))
    return;

  // create an ininitialized variable
  logpr beamThreshold((void*)0);
  if (origin.cliqueBeam != (-LZERO)) {
    // then we do clique table pruning right here rather
    // than a separate call to ceCliqueBeamPrune().
    // break into the logp to avoid unnecessary zero checking.
    beamThreshold.valref() = maxCEValue.valref() - origin.cliqueBeam;
  } else {
    // set beam threshold to a value that will never cause pruning.
    beamThreshold.set_to_zero();
  }

  const unsigned origNumCliqueValuesUsed = numCliqueValuesUsed;
  for (unsigned cvn=0;cvn<numCliqueValuesUsed;) {
    if (cliqueValues.ptr[cvn].p < beamThreshold) {
      // swap with last entry, and decrease numCliqueValuesUsed by one. We
      // swap so that entries at the end can be added back in by a future stage.
      swap(cliqueValues.ptr[cvn],cliqueValues.ptr[numCliqueValuesUsed-1]);
      numCliqueValuesUsed--;
    } else {
      cvn++;
    }
  }

  infoMsg(IM::Inference, IM::Med,"Clique beam pruning: Max cv = %f, thres = %f. Original clique state space = %d, new clique state space = %d, %2.2f%% reduction\n",
	  maxCEValue.valref(),
	  beamThreshold.valref(),
	  origNumCliqueValuesUsed,
	  numCliqueValuesUsed,
	  100*(1.0 - (double)numCliqueValuesUsed/(double)(origNumCliqueValuesUsed>0?origNumCliqueValuesUsed:1)) );

#if 0
  // A version with a bit more information printed.
  infoMsg(IM::Inference, IM::Med,"Clique beam pruning, Max cv = %f, thres = %f. Original clique state space = %d, orig sum = %f, new clique state space = %d, new sum = %f\n",
	  maxCEValue.valref(),
	  beamThreshold.valref(),
	  origNumCliqueValuesUsed,
	  numCliqueValuesUsed,
	  sumProbabilities().valref());
#endif


}



/*-
 *-----------------------------------------------------------------------
 * PedagogicalCliqueTable::ceDoAllPruning()
 *
 *    Clique Prune: This routine will do both k-pruning and beam pruning just
 *    like the routine that projects down to the next separator, but this routine does
 *    just the pruning part.
 *
 * Preconditions:
 *   same as other pruning routines
 *
 * Postconditions:
 *    Clique table has been pruned, and memory for it has been re-allocated to
 *    fit the smaller size.
 *
 * Side Effects:
 *    changes the clique size.
 *
 * Results:
 *     nothing
 *
 *-----------------------------------------------------------------------
 */
void 
PedagogicalCliqueTable::ceDoAllPruning(MaxClique& origin,
			       logpr maxCEValue)
{


  // if all observed and/or deterministic clique, then only one state,
  // so nothing to prune.
  if (origin.hashableNodes.size() == 0)
    return;

  // We seed the random number generator specifically for this clique
  // since it might get called again during island algorithm. 
  // Search for string K-BEAM-SEED elsewhere in this file for further information.
  rnd.seed(&(cliqueValues[0].p.valref()));

  const unsigned long origNumCliqueValuesUsed = numCliqueValuesUsed;
  // printf("DEBUG: orig = %lu, allo = %lu\n",origNumCliqueValuesUsed,cliqueValues.size());

  // TODO: give a command line option to change the order of these pruning
  // options.

  // First, do k-pruning (ideally we would do this after
  // beam pruning, but since beam pruning is integrated
  // into the code above, we do max state pruning first).
  // Prune the minimum of the fixed K size and the percentage size.
  unsigned k;
  k = 2 + (unsigned)((origin.cliqueBeamRetainFraction)*(double)numCliqueValuesUsed);
  if (origin.cliqueBeamMaxNumStates > 0) {
    k = min(k,origin.cliqueBeamMaxNumStates);
  }
  //   printf("nms = %d, pf = %f, ncv = %d, k = %d\n",origin.cliqueBeamMaxNumStates,
  // origin.cliqueBeamRetainFraction,numCliqueValuesUsed,k);
  // printf("starting k pruning with state space %d\n",numCliqueValuesUsed); fflush(stdout);

  if (k < numCliqueValuesUsed) {
    infoMsg(IM::Inference, IM::Med,"Clique k-beam pruning with k=%d: Original clique state space = %d\n",k,
	    numCliqueValuesUsed);
    numCliqueValuesUsed = ceCliqueStatePrune(k,cliqueValues.ptr,numCliqueValuesUsed);
  }

  // printf("ending k pruning\n"); fflush(stdout);

  // next do mass pruning.
  numCliqueValuesUsed = ceCliqueMassPrune(1.0 - origin.cliqueBeamMassRetainFraction,
					  origin.cliqueBeamMassExponentiate,
					  origin.cliqueBeamMassFurtherBeam,
					  origin.cliqueBeamMassMinSize,
					  cliqueValues.ptr,
					  numCliqueValuesUsed);

  // next, do normal beam pruning.
  ceCliqueBeamPrune(origin,maxCEValue);

  // do diversity pruning.
  // printf("starting diversity pruning with state space %d\n",numCliqueValuesUsed); fflush(stdout);
  ceCliqueDiversityPrune(origin,origin.cliqueBeamClusterPruningNumClusters);
  // printf("ending diversity pruning\n"); fflush(stdout);

  // last, add random entries back in.
  if (origin.cliqueBeamUniformSampleAmount != 0) {
    ceCliqueUniformSamplePrunedCliquePortion(origin,origNumCliqueValuesUsed);
  }


  // lastly, resize.
  // To reallocate or not to reallocate, that is the question.  here,
  // we just reallocate for now.
  // 
  // TODO: reallocate only if change is > some percentage (say 5%),
  // and export to command line.
  // e.g., if ((origNumCliqueValuesUsed - numCliqueValuesUsed) > 0.05*origNumCliqueValuesUsed)
  // TODO: if -probE option or gmtkDecode is running, no need to re-allocate here since this
  //       will soon be deleted anyway.
  // TODO: if we do a backward pass, we might not want to resize, and if an entry ends up
  //       becoming probable, re-insert the entry as a valid clique entry.

  // resize only if we've shrunk by more than about 1.6%
  // relative. Assume strength reduction will take care of fast
  // divide.


  //   if ((origNumCliqueValuesUsed - numCliqueValuesUsed) > origNumCliqueValuesUsed/64) {
  //     cliqueValues.resizeAndCopy(numCliqueValuesUsed);
  //   }
  // printf("DEBUG: old cond = %d, new cond = %d\n",((origNumCliqueValuesUsed - numCliqueValuesUsed) > origNumCliqueValuesUsed/64),((cliqueValues.size() - numCliqueValuesUsed) > cliqueValues.size()/64));
  if ((cliqueValues.size() - numCliqueValuesUsed) > cliqueValues.size()/64) {
    cliqueValues.resizeAndCopy(numCliqueValuesUsed);
  }

}

void 
PedagogicalCliqueTable::ceDoCliqueScoreNormalization(PedagogicalCliqueTable::SharedLocalStructure& 
					     sharedStructure)
{
  // syntactic convenience variables.
  MaxClique& origin = *(sharedStructure.origin);

  assert (origin.normalizeScoreEachClique != 1.0);
  
  logpr normValue;
  if (origin.normalizeScoreEachClique == 0.0) {
    // find max score and take inverse
    normValue = maxProb().inverse();
  } else {
    normValue = origin.normalizeScoreEachClique;
  }
  for (unsigned cvn=0;cvn<numCliqueValuesUsed;cvn++) {
    cliqueValues.ptr[cvn].p *= normValue;
  }
}




#ifdef USE_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL
/*-
 *-----------------------------------------------------------------------
 * PedagogicalCliqueTable::insertLocalCliqueValuesIntoSharedPool();
 *
 *    This routine takes all the clique entries and values, and inserts
 *    the values that live currently in the local clique value pool into
 *    the globally shared clique value pool.
 *
 * Preconditions:
 *   The clique must currently be set up so that the clique entries indicate 
 *   an index into the local clique value pool. Note that the clique value pointers
 *   are such that they give the actual integer index into the exact location
 *   of the current local clique value pool (i.e., it is the entry location
 *   multiplied in by the number of words per clique value).
 *
 *   Also, this routine must not be called unless (origin.packer.packedLen() > IMC_NWWOH) is true.
 *
 * Postconditions:
 *    The clique entries are now set so that the values are pointers into
 *    the potentially shared value location in the globally shared pool.
 *
 * Side Effects:
 *    changes the clique entry value pointers.
 *
 * Results:
 *     nothing
 *
 *-----------------------------------------------------------------------
 */
void 
PedagogicalCliqueTable::insertLocalCliqueValuesIntoSharedPool(MaxClique& origin)
{

  // clique value length in bytes
  const unsigned cvlb = origin.packer.packedLenInBytes();

  unsigned k;
  for (k=0;k<numCliqueValuesUsed;k++) {
    // First, grab pointer to storge where next clique value would
    // be stored if it ends up being used.
    unsigned *pcv = origin.valueHolder.curCliqueValuePtr();
    // next, grab pointer to storge where the local clique value currently lives.
    unsigned *lcv = &origin.temporaryCliqueValuePool.ptr[cliqueValues.ptr[k].ival];
    // copy from the local location to the permanent storage pool.
    // TODO: have a packer optimized copy routine here.
    ::memcpy((void*)pcv,(void*)lcv,cvlb);
    // Look it up in the hash table.
    bool foundp;
    unsigned *key;
    key = origin.cliqueValueHashSet.insert(pcv,foundp);
    if (!foundp) {
      // if it was not found, need to claim this storage that we
      // just used.
      origin.valueHolder.allocateCurCliqueValue();
    } 
#ifdef TRACK_NUM_CLIQUE_VALS_SHARED	
    else
      numCliqueValuesShared++;
#endif
    // Save the pointer to whatever the hash table decided to use.
    cliqueValues.ptr[k].ptr = key;
  }

}
#endif


/*-
 *-----------------------------------------------------------------------
 * PedagogicalCliqueTable::ceCliqueStatePrune(unsigned k)
 *
 *    Collect Evidence, Clique Prune: This routine will organize a 
 *    reviously instantiated clique table so that the top k entries
 *    are at the beginning of the table (in arbitrary order) and the
 *    rest of the entries start at the k+1'st position. It does
 *    this using a fast O(N) algorithm.
 *    It returns k, so that pruning can be implemented by calling
 *    this routine and taking/using only the first k entries on return.
 *    It does not do any change of array allocation, so that can
 *    be done by the caller, and only does array organization
 *    within the bounds of the array so given by the arguments.
 *
 *    We can also quickly find the k top entries for Viterbi decoding
 *    by calling this, and then sorting the top k entries.
 *
 * Preconditions:
 *   1) the value of the max clique 'maxCEValue' must have been
 *      computed already.
 *
 *   2) clique table must be created, meaning that either:
 *
 *        PedagogicalCliqueTable::ceIterateAssignedNodesRecurse()
 *
 * Postconditions:
 *
 *    Clique table has been reorganized so that the k top entries are
 *    at the beginning of the table. Memory for it has NOT been
 *    re-allocated and no entries have been removed.
 *
 * Side Effects:
 *    reorganizes the clique table. Note that this routine does
 *    nothing if k >= the size of the table. 
 *
 * Results:
 *     min(k,length of the given array)
 *
 *-----------------------------------------------------------------------
 */

unsigned
PedagogicalCliqueTable::ceCliqueStatePrune(const unsigned k,
				   CliqueValue* curCliqueVals,
				   const unsigned curNumCliqueValuesUsed)
{
  // k can't be larger than the number of clique entries.
  if (k == 0 || k >= curNumCliqueValuesUsed)
    return curNumCliqueValuesUsed;

  // inclusive range to process.
  unsigned lower = 0;
  unsigned upper = curNumCliqueValuesUsed-1;


  // We find the k max elements using a quick sort-like algorithm, and
  // start by processing the elements between lower and upper
  // inclusive. Note that unlike quicksort (which is O(NlogN) average
  // case and O(N^2) worse case), this algorithm is O(N) average case
  // (yes, not dep. on k) since we're only finding the top k values,
  // and those top k values need not be sorted.

  do {
    // k'th largest lives within [lower,upper] inclusive.
    // printf("lower = %d, upper = %d\n",lower,upper);

    if (lower == upper)
      break;

    // Choose a random pivot and then divide the array so that
    // elements >= pivot are on the left, and elements < pivot are on
    // the right. We choose a pivot by taking the median of 3 randomly
    // choosen values.
    unsigned pl1 = rnd.uniform(lower,upper);
    unsigned pl2 = rnd.uniform(lower,upper);
    unsigned pl3 = rnd.uniform(lower,upper);
    unsigned pl;
    // find rough median
    if (curCliqueVals[pl1].p < curCliqueVals[pl2].p) {
      if (curCliqueVals[pl2].p < curCliqueVals[pl3].p) {
	pl = pl2;
      } else {
	pl = pl3;
      }
    } else {
      if (curCliqueVals[pl1].p < curCliqueVals[pl3].p) {
	pl = pl1;
      } else {
	pl = pl3;
      }
    }

    // Uncomment to give a valid but fixed deterministic pivot just for testing.
    // pl = (lower+upper)/2;

    logpr pivot = curCliqueVals[pl].p;
    // printf("pivot location = %d, pivot = %f\n",pl,pivot);

    // swap pivot into first position, which is a valid position for
    // pivot, since pivot >= pivot.
    swap(curCliqueVals[pl],curCliqueVals[lower]);

    unsigned l = lower+1;
    unsigned u = upper;
    while (l < u) {
      if (curCliqueVals[l].p >= pivot) {
	l++;
      } else {
	swap(curCliqueVals[l],curCliqueVals[u]);
	u--;
      }
    }

    // 1) Now all entries with index < l are >= pivot,
    //   which means at this point we know we have
    //   the top l entries at indices [0,l]
    // 2) All entries with index > u are < pivot
    // 3) we have l == u
    // 4) The entry at index u==l is unknown however.


    if (curCliqueVals[l].p >= pivot) {
      l++;
      u++;
    }
    // 1) Now all entries with index < l are >= pivot,
    //   which means at this point we know we have
    //   the top l entries at indices [0,l]
    // 2) All entries with index >= u are < pivot
    // 3) we still have l == u

    // put pivot in its appropriate place.
    l--;
    swap(curCliqueVals[lower],curCliqueVals[l]);
    // 1) Now all entries with index <= l are >= pivot,
    //   which means at this point we know we have
    //   the top l entries at indices [0,l]
    // 2) All entries with index >= u are < pivot
    // 3) we have l +1 == u

    // Now, deal with potential ties.
    // Move any other entries that are equal to pivot to the center.
    // -- Note that this next step helps to significantly speed up the
    // -- algorithm in the case when there are lots of ties --- this
    // -- is done here since one of the main reasons for having this
    // -- form of pruning is that normal beam pruning doesn't work
    // -- well on very "uniform"-like unnormalized distributions,
    // -- i.e., ones for which many ties might be present.
    unsigned ll;
    if (l > lower) {
      ll = l-1; // left of left-most known pivot value.
      unsigned i = lower;
      while (i < ll) {
	if (curCliqueVals[i].p == pivot) {
	  swap(curCliqueVals[i],curCliqueVals[ll]);
	  ll--;
	} else {
	  i++;
	}
      }
    } else
      ll = l;

    // 1) Now all entries with index <= l are >= pivot,
    //   which means at this point we know we have
    //   the top l entries at indices [0,l]
    // 2) All entries with index >= u are < pivot
    // 3) we have l +1 == u
    // 4) either
    //    a) all entries with index <= ll are > pivot, or
    //         (meaning that ll+1 is the number of entries
    //          that are strictly greater than pivot)
    //    b) ll = lower and arr[ll] == pivot
    //       (in which case, ll is the number of entries
    //          that are strictly greater than pivot).

    //     printf("after swapping, l = %d, u = %d,array now:\n",l,u);
    //     for (unsigned j=0;j<len;j++)
    //       printf("arr[%d] = %f\n",j,arr[j]);

    if (ll > k) {
      upper = ll;
    } else if (l > k) {
      // then we're done since we must
      // have a bunch of ties, and pivot value
      // must be correct k'th value.
      break;
    } else if (l == k) {
      // then we have l+1 entries
      // printf("finished, we have l == k == %d\n",l);
      break;
    } else { // l < k, need to search to the right.
      // printf("setting lower to %d\n",u);
      lower = u;
    }

  } while(1);

  return k;

}


/*
 * structure used only for diversity pruning
 */
struct DistClust {
  // the distance of this point to its cluster rep (its nearest center)
  unsigned distance;
  // the cluster number of this point's cluster (i.e., its nearest center).
  unsigned cluster;
};


/*
 * TODO: this next macro is probably generally useful, so ultimately
 * put it somewhere else.
 */

#ifdef USE_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL
#define EXTRACT_KEY_FROM_CVN_INTO_KEY_P(__cvn,__key_p) \
    if (origin.packer.packedLen() <= IMC_NWWOH) { \
      (__key_p) = &(cliqueValues.ptr[(__cvn)].val[0]); \
    } else { \
      (__key_p) = \
	&origin.temporaryCliqueValuePool.ptr[cliqueValues.ptr[(__cvn)].ival]; \
    }
#else
#define EXTRACT_KEY_FROM_CVN_INTO_KEY_P(__cvn,__key_p) \
    if (origin.packer.packedLen() <= IMC_NWWOH) { \
      (__key_p) = &(cliqueValues.ptr[(__cvn)].val[0]); \
    } else { \
      (__key_p) = (unsigned*)cliqueValues.ptr[(__cvn)].ptr; \
    }
#endif





void 
PedagogicalCliqueTable::ceCliqueDiversityPrune(MaxClique& origin,const unsigned numClusters)
{
  /*
   * there are two forms of clustering algorithm. 
   * 
   * 1) The first one clusters only by diversity (i.e., it chooses for
   * the next center, the point that maximizes the min distance to any
   * cluster. That is, each point has a cluster center that is
   * assigned based on how close it is to that center (each point is
   * assigned to the nearest cluster center), and the next center is
   * the point that is maximally distant from its center. 
   *
   * 2) the next cluster algorithm is quite similar except that when
   * there are multiple points that are maximally distant from their
   * centers, rather than breaking ties arbitrarily, we choose the
   * point that has has the highest current score.
   *
   * To get the first algorithm, do not #define the macro DIVERSITY_PRUNE_SCORE_BASED
   * To get the second algorithm, #define the macro DIVERSITY_PRUNE_SCORE_BASED
   */

#define DIVERSITY_PRUNE_SCORE_BASED


  // k can't be larger than the number of clique entries.
  if ((origin.cliqueBeamClusterBeam == (-LZERO)  && 
       (MaxClique::cliqueBeamClusterMaxNumStates 
	== NO_PRUNING_CLIQUEBEAMCLUSTERPRUNINGMAXSTATES)
       && 
       origin.cliqueBeamClusterRetainFraction == 1.0
       && 
       origin.cliqueBeamClusterMassRetainFraction == 1.0)
      || numClusters >= numCliqueValuesUsed) {
    return;
    // note that setting the number of clusters to 1 should produce
    // results identical to the no-clustering case, but we allow
    // that case to pass through here for consistency, sanity checking
    // completeness, and since the k=1 case does not complicate the
    // code below (i.e., it need not be a special case).
  }

  infoMsg(IM::Inference, IM::Med+9,"Diversity/cluster pruning with state space = %d\n",numCliqueValuesUsed);

  /*
    - ideas to speed up: 
    -- factor malloc outside of this routine.
    -- use smaller data types for arrays
    -- use a faster distance function
    -- use a log k data structure ala Feder&Greene (but then constants might
    be worse, and this would be only for relatively large k)
    -- compile with -fno-exceptions (to possibly make routine calls faster)
    -- 

  */

  unsigned long* centers = new unsigned long[numClusters];

  // printf("begin allocate\n");fflush(stdout);
  DistClust* distClusts = new DistClust[numCliqueValuesUsed];
  // printf("end allocate\n");fflush(stdout);


  // find the entry with the max score. The first one is arbitrary.
  centers[0] = 0;
  logpr maxScore = cliqueValues.ptr[0].p;
  // We could also pick a random point as the first cluster rather
  // than the first one, but this makes each run of the inference different.
  // centers[0] = rnd.uniformOpen(numCliqueValuesUsed);

#ifdef DIVERSITY_PRUNE_SCORE_BASED
  // we need to start with the one that has the highest score.  I.e.,
  // the first center is one having the highest score.
  for (unsigned cvn=1;cvn<numCliqueValuesUsed;cvn++) {
    // TODO: unroll this.
    if (cliqueValues.ptr[cvn].p > maxScore) {
      maxScore = cliqueValues.ptr[cvn].p;
      centers[0] = cvn;       
    }
  }
#else
  // nothing to do here.
#endif

  // get the clique entry corresponding to centers[0] which at
  // this point is the clique entry with the max score.
  unsigned* center_key_p;
  EXTRACT_KEY_FROM_CVN_INTO_KEY_P(centers[0],center_key_p);

  unsigned maxDist = 0;
  unsigned maxIndx = 0;
  maxScore.set_to_zero();

  for (unsigned cvn=0;cvn<numCliqueValuesUsed;cvn++) {
    unsigned* current_key_p;
    EXTRACT_KEY_FROM_CVN_INTO_KEY_P(cvn,current_key_p);
    distClusts[cvn].cluster = 0;
    distClusts[cvn].distance = 
      origin.packer.use_distance(center_key_p,current_key_p);
    // find max while we're at it.
#ifdef DIVERSITY_PRUNE_SCORE_BASED
    if (distClusts[cvn].distance > maxDist
	|| (((distClusts[cvn].distance == maxDist && cliqueValues.ptr[cvn].p > maxScore)))) {
      maxIndx = cvn;
      maxDist = distClusts[cvn].distance;
      maxScore = cliqueValues.ptr[cvn].p;
    }
#else
    if (distClusts[cvn].distance > maxDist) {
      maxIndx = cvn;
      maxDist = distClusts[cvn].distance;
      maxScore = cliqueValues.ptr[cvn].p;
    }
#endif
  }

  unsigned k = 1;
  while (k < numClusters) {

    // assign new center using previously computed max index.
    centers[k] = maxIndx;

    // now we need to update all distances. Each min distance is such
    // that it is either the same (i.e., it's cluster center has not
    // changed) or it has decreased (it has a new cluster center). So
    // this means that the distance can never increase (the distances
    // are monotonically decreasing). Therefore, we only need to check
    // against the newly assigned cluster, not all the rest of them.

    EXTRACT_KEY_FROM_CVN_INTO_KEY_P(centers[k],center_key_p);

    // now compute the new max distance and cluster membership.
    maxDist = 0;
    maxIndx = 0;
    maxScore.set_to_zero();
    for (unsigned cvn=0;cvn<numCliqueValuesUsed;cvn++) {

      unsigned* current_key_p;
      EXTRACT_KEY_FROM_CVN_INTO_KEY_P(cvn,current_key_p);

      unsigned cur_dist = 
	origin.packer.use_distance(center_key_p,current_key_p);

      // now check if the distance from the current point to the new
      // center is lower than its current distance, and if it is then
      // this new point gets a new cluster. We break ties by keeping
      // the point assigned to the old cluster, but it might be better
      // to break ties randomly, or to assign it to the cluster based
      // on the cluster center's score.
      if (cur_dist < distClusts[cvn].distance) {
	distClusts[cvn].distance = cur_dist;
	distClusts[cvn].cluster = k;
      }

      // we need unfortunately to re-compute the new max from scratch
      // since the previous max could have been overwritten. One idea:
      //  1) if each cluster keeps track of its own max, and if a cluster
      //     isn't touched by the finding of a new cluster (i.e., the
      //     cluster stays intact), then the points in that cluster
      //     need not be max searched again, only the max point

#ifdef DIVERSITY_PRUNE_SCORE_BASED
      if (distClusts[cvn].distance > maxDist
	  || (((distClusts[cvn].distance == maxDist && cliqueValues.ptr[cvn].p > maxScore)))) {
	maxIndx = cvn;
	maxDist = distClusts[cvn].distance;
	maxScore = cliqueValues.ptr[cvn].p;
      }
#else
      if (distClusts[cvn].distance > maxDist) {
	maxIndx = cvn;
	maxDist = distClusts[cvn].distance;
	maxScore = cliqueValues.ptr[cvn].p;
      }
#endif
    }
    k++;
  }

  // so now we have all the cluster centers and the assignment of all
  // points to their clusters based on the nearest center.

  // Initialize cluster sizes.
  unsigned*  orig_cluster_sizes = new unsigned[numClusters];
  for (k=0;k<numClusters;k++) {
    orig_cluster_sizes[k] = 0;
  }


  const unsigned origNumCliqueValuesUsed = numCliqueValuesUsed;

  // printf("Trying beam pruning\n");

  // Next do normal beam pruning, this next step costs O(n).
  if (MaxClique::cliqueBeamClusterBeam != (-LZERO)) {

    // First, calculate the max score value within each cluster.
    // Note, default values of logpr are set to zero.
    // We also here calculate the cluster sizes.
    logpr* intra_cluster_max_values = new logpr[numClusters];
    for (unsigned cvn=0;cvn<numCliqueValuesUsed;cvn++) {
      const unsigned clust = distClusts[cvn].cluster;
      orig_cluster_sizes[clust]++;
      if (cliqueValues.ptr[cvn].p > intra_cluster_max_values[clust]) {
	intra_cluster_max_values[clust] = cliqueValues.ptr[cvn].p;
      }
    }


    // now we have the clustering do the pruning. 
    // printf("now we have the clustering do the pruning. \n"); fflush(stdout);
    if (IM::messageGlb(IM::Med+9)) {
      printf("Clique cluster beam pruning: Orig cluster-%d sizes: ",numClusters);
      for (k=0;k<numClusters;k++) {
	printf("%u ",orig_cluster_sizes[k]);
	// turn max values into the needed threshold
      }
      printf("\n");
    }


    // Next, we need to calculate the per-cluster pruning threshold:
    for (k=0;k<numClusters;k++) {
      // turn max values into the needed threshold
      intra_cluster_max_values[k].valref() = 
	intra_cluster_max_values[k].valref() - MaxClique::cliqueBeamClusterBeam;
    }

    // Next, we do the actual pruning, and we do this without reordering the
    // clique table entries.
    for (unsigned cvn=0;cvn<numCliqueValuesUsed;) {
      const unsigned clust = distClusts[cvn].cluster;
      if (cliqueValues.ptr[cvn].p < intra_cluster_max_values[clust]) {
	// then we do a prune of this entry.

	assert ( clust < numClusters );

	orig_cluster_sizes[clust]--;

	// swap with last entry, and decrease numCliqueValuesUsed by
	// one.  We swap rather than assign just in case some other
	// method wants to add them back in (but this is slower). If
	// we are assured that this clustering method is the last one
	// called, we can turn the swap into an assignment.
	swap(cliqueValues.ptr[cvn],cliqueValues.ptr[numCliqueValuesUsed-1]);

	// WARNING!! distance values for current entry (cvn) will be invalid
	// after this next step since we only assign the cluster id
	// member in a DistClust object rather than the entire
	// structure. We do this for the sake of speed, and since the
	// distance value is no longer needed at this point and 
	// presumably below.
	distClusts[cvn].cluster = 
	  distClusts[numCliqueValuesUsed-1].cluster;

	numCliqueValuesUsed--;
      } else {
	// the entry 'cvn' survived, so no pruning.
	cvn++;
      }
    }
    // optionally print information done after beam pruning.
    infoMsg(IM::Inference, IM::Med,"Clique cluster beam pruning: Original clique state space = %d, new clique state space = %d, %2.2f%% reduction\n",
	    origNumCliqueValuesUsed,
	    numCliqueValuesUsed,
	    100*(1.0 - (double)numCliqueValuesUsed/(double)(origNumCliqueValuesUsed>0?origNumCliqueValuesUsed:1)) );
    
    if (IM::messageGlb(IM::Med+9)) {
      printf("Clique cluster beam pruning: Post cluster-%d sizes: ",numClusters);
      for (k=0;k<numClusters;k++) {
	printf("%u ",orig_cluster_sizes[k]);
      }
      printf("\n");
    }
    delete [] intra_cluster_max_values;
  } else {
    // still need to calculate cluster sizes
    for (unsigned cvn=0;cvn<numCliqueValuesUsed;cvn++) {
      const unsigned clust = distClusts[cvn].cluster;
      orig_cluster_sizes[clust]++;
    }
  }
  
  // printf("************* trying div state pruning, a = %ul, b = %ul\n",MaxClique::cliqueBeamClusterMaxNumStates,NO_PRUNING_CLIQUEBEAMCLUSTERPRUNINGMAXSTATES);

  if (origin.cliqueBeamClusterMaxNumStates != NO_PRUNING_CLIQUEBEAMCLUSTERPRUNINGMAXSTATES
      ||  
      origin.cliqueBeamClusterRetainFraction < 1.0
      || 
      origin.cliqueBeamClusterMassRetainFraction < 1.0) {

    // Now we do k-beam pruning. The algorithm is to:
    //  1) 'sort' entire list by cluster number, O(n)
    //  2) run k pruning within each set, amortized O(n) cost.
    //  3) reorganize and put pruned guys at the bottom.
    // total pruning cost is O(n log n) (which is probably ok).

    // - need to modify k-prune to take an arguemnt to pointer to clique
    //  valus + length and return new length (and have that array is
    //  modified to have top k at the top).

    // we can use the STL to create a container class that sorts
    // one array using keys from the other.
    // Also note, if for given s and k, sk >= current size, then just return.


    // First, we need to sort the clique values using the clusters,
    // rather than the values, as the key. This will organize the table
    // so that cluster 0 elements will come first, cluster 1 elements will come next, and
    // so on. This can be done in O(n) time rather than doing an O(nlogn) sort.
    // Once done, the entries will be in arbitrary order in each cluster.

    // now calculate the cluster sizes
    unsigned*  cluster_starts = new unsigned[numClusters+1];
    unsigned*  cluster_endp = new unsigned[numClusters]; // current ends
    cluster_endp[0] = cluster_starts[0] = 0;
    for (k=1;k<numClusters;k++) {
      cluster_endp[k] = 
	cluster_starts[k] = cluster_starts[k-1] + orig_cluster_sizes[k-1];
    }
    // include this one for convenience to the code below.
    cluster_starts[numClusters] = numCliqueValuesUsed;

    unsigned cvn=0;
    // pos_cvn_k is the cluster that the position 'cvn' *should* be.
    // Since we are sweeping increasing in cvn, cvn_pos_k starts
    // at zero and increments up to 'numClusters' and switches
    // at teh boundaries given by the 'cluster_starts' array.
    unsigned pos_cvn_k = 0; // this corresponds to cvn == 0.

    // we use a separate index since the loop never needs to occur
    // more than numCliqueValuesUsed times.
    for (unsigned i=0;i<numCliqueValuesUsed;i++) {

      // what is the cluster of current entry at position cvn
      const unsigned cur_cvn_k = distClusts[cvn].cluster;

      if (cur_cvn_k == pos_cvn_k) {
	// then the guy at postion cvn is already in the right partition
	// of the array.
	if (cvn == cluster_endp[cur_cvn_k]) {
	  // then we haven't accounted for this cluster member
	  // yet and we must do so now.
	  cluster_endp[cur_cvn_k]++;
	} else {
	  assert ( cvn < cluster_endp[cur_cvn_k] );
	  // otherwise, this cluster member must be here because
	  // we swaped it here before.
	}

	// move to next entry
	cvn++;

	// deduce and update what the cluster should be of the next cvn.
	if (cvn >= cluster_starts[pos_cvn_k+1])
	  pos_cvn_k++;
      } else {
	// swap the curent one to where it belongs and continue
	swap(cliqueValues.ptr[cvn],cliqueValues.ptr[cluster_endp[cur_cvn_k]]);

	// this next step does not swap the distances only the cluster
	// member items, so leaves the distances totally invalid. This
	// is ok, presumably since the distances are no longer being
	// used and may be a bit faster since there is less swapping
	// going on.
	swap(distClusts[cvn].cluster,
	     distClusts[cluster_endp[cur_cvn_k]].cluster);
	// we've got a new end placement
	cluster_endp[cur_cvn_k]++;
      }
    }

    // Next, we do state-beam pruning within each cluster. We call the state-pruning
    // routine which will leave the top m=maxStateSize clique entries at the beginning
    // of each cluster and the rest will be after the top m=maxStateSize (which will
    // be useful for other forms of clustering later). We place
    // the new sizes of the clusters in the 'cluster_endp' array -- note
    // that the sizes are at most m but could be < m in which case
    // no pruning is done of course. The entire amortized process is O(n)
    // since each sub-step is O(size of cluster k) and 
    // O(n) = \sum_k O(size of cluster k). Note, that we calcluate
    // a local maxStateSize(k) for cluste k based on the global constraint
    // on the max state size within each cluster, and the percentage
    // reduction to retain, taking the minimum of the two.
    unsigned newNumCliqueValuesUsed = numCliqueValuesUsed;
    if (origin.cliqueBeamClusterMaxNumStates != NO_PRUNING_CLIQUEBEAMCLUSTERPRUNINGMAXSTATES
	||  
	origin.cliqueBeamClusterRetainFraction < 1.0) {
      newNumCliqueValuesUsed = 0;
      for (k=0;k<numClusters;k++) {
	
	unsigned maxStateSize = 
	  2 + (unsigned)((origin.cliqueBeamClusterRetainFraction)
			 *(double)orig_cluster_sizes[k]);
	
	if (MaxClique::cliqueBeamClusterMaxNumStates != NO_PRUNING_CLIQUEBEAMCLUSTERPRUNINGMAXSTATES)
	  maxStateSize = min(MaxClique::cliqueBeamClusterMaxNumStates,
			   maxStateSize);

	// what is returned is the new cluster size, which is
	// the min of the original size or the within-cluster state space.
	// We re-use the cluster_endp array that was used as temporary
	// storage above to store the new cluster sizes.
	cluster_endp[k] = 
	  ceCliqueStatePrune(maxStateSize,
			     cliqueValues.ptr + cluster_starts[k],
			     orig_cluster_sizes[k]);
	newNumCliqueValuesUsed += cluster_endp[k];
      }

      assert ( newNumCliqueValuesUsed <= numCliqueValuesUsed );
    } else {
      // we still need to update the sizes
      for (k=0;k<numClusters;k++) {
	cluster_endp[k] = orig_cluster_sizes[k];
      }
    }

    //////////////////////////////////////////////////////////////////////
    // next, we do mass percentage based pruning within each
    // cluster. This is pretty easy given that the clique is organized
    // as it is.
    if (origin.cliqueBeamClusterMassRetainFraction < 1.0) {
      newNumCliqueValuesUsed = 0;
      for (k=0;k<numClusters;k++) {
	cluster_endp[k] =
	  ceCliqueMassPrune(1.0 - origin.cliqueBeamClusterMassRetainFraction,
			    origin.cliqueBeamClusterMassExponentiate,
			    origin.cliqueBeamClusterMassFurtherBeam,
			    origin.cliqueBeamClusterMassMinSize,
			    cliqueValues.ptr + cluster_starts[k],
			    cluster_endp[k]);
	newNumCliqueValuesUsed += cluster_endp[k];
      }
      assert ( newNumCliqueValuesUsed <= numCliqueValuesUsed );
    }


    // Last, we need to re-organize the entire clique table so that
    // all of the to-be pruned entries are at the end. We do the
    // sweep in place, leading to another O(n) algorithm.
    
    // For better cache behavior, we have two pointers scan the array,
    // one that is ahead of the other. The earlier pointer points to
    // the next element that needs to be pruned and the later pointer
    // points to the next element that should not be pruned.  We
    // always maintain that earlier_p < later_p as we scan the
    // array. We also maintain the invariant that everything between
    // earlier_p and later_p should be pruned, so that as soon as
    // earlier_p has enough entries before it, we can stop.
    // Note also that later_p could touch (and prefetch into cache)
    // some memory items that are later used by earlier_p (this
    // wouldn't be the case if we used an alternative approach where
    // we had two pointers, one initalized at the start and one at the
    // end of the array, and where the pointers appraoch each other
    // (respectively incrementing/decrementing). 

    if (newNumCliqueValuesUsed < numCliqueValuesUsed) {
      // We first find the first "hole", i.e., location where
      // something should be moved.
      for (k=0;k<numClusters;k++) {
	if (cluster_endp[k] < orig_cluster_sizes[k])
	  break;
      }
      assert (k < numClusters);

      // next, do the reorganization.
      unsigned early_loc = cluster_starts[k] + cluster_endp[k];
      k++;
      unsigned later_loc = cluster_starts[k];
      unsigned later_loc_end;
      if (k < numClusters) 
	later_loc_end = cluster_starts[k] + cluster_endp[k];
      else 
	later_loc_end = numCliqueValuesUsed;
      while (early_loc < newNumCliqueValuesUsed) {
	// swap earlier (to be pruned) and later (not to be pruned)
	// locations.
	swap(cliqueValues.ptr[early_loc],
	     cliqueValues.ptr[later_loc]);

	// update earlier pointer
	early_loc++;

	// update later pointer
	later_loc++;
	if (later_loc == later_loc_end) {
	  k++;
	  if (k < numClusters) {
	    later_loc = cluster_starts[k];
	    later_loc_end = cluster_starts[k] + cluster_endp[k];
	  } else {
	    // this should only happen at the last iteration.
	    assert ( early_loc  ==  newNumCliqueValuesUsed );
	  }
	}
      }
      numCliqueValuesUsed = newNumCliqueValuesUsed;
    }

    infoMsg(IM::Inference, IM::Med,"Clique cluster state pruning: Original clique state space = %d, new clique state space = %d, %2.2f%% reduction\n",
	    origNumCliqueValuesUsed,
	    numCliqueValuesUsed,
	    100*(1.0 - (double)numCliqueValuesUsed/(double)(origNumCliqueValuesUsed>0?origNumCliqueValuesUsed:1)) );

    if (IM::messageGlb(IM::Med+9)) {
      printf("Clique cluster state pruning: Post cluster-%d sizes: ",numClusters);
      for (k=0;k<numClusters;k++) {
	printf("%u ",cluster_endp[k]);
      }
      printf("\n");
    }

    // free memory for k-beam cluster pruning
    delete [] cluster_starts;
    delete [] cluster_endp;

  }


  // free memory for general stuff.
  delete [] centers;
  delete [] distClusts;
  delete [] orig_cluster_sizes;

}


/*-
 *-----------------------------------------------------------------------
 * PedagogicalCliqueTable::ceCliqueMassPrune(double removeFraction, unsigned minSize)
 *
 *    Collect Evidence, Clique Mass Prune: This routine will prune away
 *    part of a previously instantiated clique so that it has only
 *    'fraction' of the clique mass left. It does this by sorting
 *    and choosing the top entries such that the mass is retained. This
 *    is based on an idea by Andrew McCallum (2005).
 *
 * Preconditions:
 *   1) clique table must be created, meaning that either:
 *
 *        PedagogicalCliqueTable::ceIterateAssignedNodesRecurse()
 *
 * Postconditions:
 *    Clique table has been pruned, and memory for it has NOT been re-allocated to
 *    fit the smaller size.
 *
 * Side Effects:
 *    changes the clique size.
 *
 * Results:
 *     nothing
 *
 *-----------------------------------------------------------------------
 */

unsigned
PedagogicalCliqueTable::ceCliqueMassPrune(const double removeFraction,
				  const double exponentiate,
				  const double furtherBeam,
				  const unsigned minSize,
				  CliqueValue* curCliqueVals,
				  const unsigned curNumCliqueValuesUsed)
{

  if (removeFraction <= 0.0 || curNumCliqueValuesUsed <= minSize)
    return curNumCliqueValuesUsed;

  // sort all current clique values descending based on clique mass value.
  // Note, sorting won't be correct if the exponent is negative.
  if (exponentiate < 0) {
    error("ERROR: trying to do exponentiated mass clique pruning with a negative exponent (%e). Exponent must be non-negative for sensible pruning.\n",exponentiate);
  }
  sort(curCliqueVals,curCliqueVals + curNumCliqueValuesUsed,CliqueValueDescendingProbCompare());

  // logpr loc_maxCEValue = curCliqueVals[0].p;
  // printf("mass pruning: maxVal %.18e, minVal %.18e\n",
  // loc_maxCEValue.val(),curCliqueVals[curNumCliqueValuesUsed-1].p.val());

  logpr origSum = sumExponentiatedProbabilities(exponentiate,
						curCliqueVals,
						curNumCliqueValuesUsed);

  if (origSum.zero())
    return curNumCliqueValuesUsed;

  // logpr desiredSum = origSum* (1.0 - removeFraction);
  logpr desiredSum = (origSum - (origSum*removeFraction));
  logpr actualSum;

#if 0
  printf("DEBUG: mass pruning: origSum %.18e, desiredSum %.18e, diff %.18e\n",
   	 origSum.val(),
   	 desiredSum.val(),
   	 (origSum - desiredSum).val());
#endif
  
  unsigned k;
  for (k=0;k<curNumCliqueValuesUsed;k++) {

    actualSum += curCliqueVals[k].p.pow(exponentiate); // /loc_maxCEValue;

    // printf("k=%d: origSum = %.16e, desiredSum = %.16e, actualSum = %.16e\n",k,origSum.valref(),desiredSum.valref(),actualSum.valref());

    // could use either ">" or ">=" here.
    //   - use ">=" if you want more aggressive pruning (will probably
    //       want to use a bigger minSize in this case).
    //   - use ">" if you want less agressive pruning.
    // This can actually have a big effect since due to the dynamic
    // range of the elements in the clique, the rest of the clique
    // when added to the current sum might not cause any increment at
    // all (all it needs to be is different by less than the min
    // difference, see logp.h). If we were to use ">" here, then
    // we would continue to add all the rest of the clique while
    // not changing actualSum, so no pruning would be done. We
    // therefore use ">=".
    if (actualSum >= desiredSum) {
      // increment k so that it becomes a count rather than an index.
      k++;
      break;
    }
  }
  
  if (furtherBeam != 0.0 && k < curNumCliqueValuesUsed ) {
    logpr curMax = curCliqueVals[k].p;
    logpr threshold = curMax/logpr((void*)0,furtherBeam);
    while (++k < curNumCliqueValuesUsed) {
      if (curCliqueVals[k].p  < threshold)
	break;
    }
  }

  //   // optional code to calculate and print residual (stuff that is pruned away)
  //   unsigned r=k;
  //   logpr residualSum;
  //   for (r=k+1;r<curNumCliqueValuesUsed;r++) {
  //     residualSum += curCliqueVals[r].p; // /loc_maxCEValue;
  //   }
  //   printf("mass pruning: removeFraction %.18e, maxCV = %.18e, actualSum = %.18e, residualSum = %0.18e, actual+residual = %0.18e\n",
  // 	 removeFraction,loc_maxCEValue.val(),actualSum.val(),residualSum.val(),
  // 	 (actualSum+residualSum).val());
  
  
  if (k<minSize)
    k=min(minSize,curNumCliqueValuesUsed);

  const unsigned newStateSpace = min(k,curNumCliqueValuesUsed);
  infoMsg(IM::Inference, IM::Med,"Clique mass-beam pruning: Original state space = %d (exp-mass=%e), new state space = %d, reduction %2.2f%%, desired exp-mass = %e, actual exp-mass = %e\n",
	  curNumCliqueValuesUsed,
	  origSum.val(),
	  newStateSpace,
	  100.0*(1.0 - (double)newStateSpace/(double)curNumCliqueValuesUsed),
	  desiredSum.val(),
	  actualSum.val());
  return newStateSpace;

}


/*-
 *-----------------------------------------------------------------------
 * PedagogicalCliqueTable::ceCliqueUniformSamplePrunedCliquePortion()
 *
 *    After pruning has occured, this routine will uniformly sample a
 *    certain portion of the entries that have been pruned away, and
 *    place them back in the clique. It does this uniformly, i.e.,
 *    irrespectively of the mass/score on each clique entry. The
 *    reason for this is that if pruning by score (which is what all
 *    of the pruning methods so far do) leads us to get a zero
 *    probability, then we add a bit of noise irrespective of the
 *    probabilities back in to the inference so that hopefully one of
 *    the entries will allow us to get a non-zero score at the very
 *    end of inference.
 *  
 *    TODO: select entries rather than uniformly at random, but based also on
 *          some form of "diversity", so that for the number of re-samples, we
 *          get as different a set of entries as possible.
 *    
 *
 * Preconditions:
 *   1) clique table must be created, meaning that either:
 *
 *        PedagogicalCliqueTable::ceIterateAssignedNodesRecurse()
 *
 *   We assume that the clique has been pruned in the past, but that
 *   the memory for the clique still exists in the clique table (cliqueValues), meaning
 *   that no memory reallocation has yet been done based on prunning. We also
 *   assume that the argument 'origNumCliqueValuesUsed' contains the size of the
 *   clique before any pruning has been done, and that 'numCliqueValuesUsed' gives
 *   the number of current clique values used, and that all entries between entry
 *   numCliqueValuesUsed and origNumCliqueValuesUsed-1 inclusive are previously pruned
 *   cliques entries.
 *   
 *
 * Postconditions:
 *    Clique table has been "un-pruned", and memory for it has NOT been re-allocated to
 *    fit the smaller size.
 *
 * Side Effects:
 *    changes the clique size.
 *
 * Results:
 *     nothing
 *
 *-----------------------------------------------------------------------
 */
void 
PedagogicalCliqueTable::ceCliqueUniformSamplePrunedCliquePortion(MaxClique& origin,
							 const unsigned origNumCliqueValuesUsed)
{

  const unsigned numCliqueValuesUsedBeforeSampling = numCliqueValuesUsed;

  if (origin.cliqueBeamUniformSampleAmount == 0.0)
    return;
  else if (origin.cliqueBeamUniformSampleAmount == 1.0) {
    numCliqueValuesUsed = origNumCliqueValuesUsed;
  } else {

    unsigned numEntriesPruned = origNumCliqueValuesUsed - numCliqueValuesUsed;
    if (numEntriesPruned != 0) {
      unsigned numToSample;
      if (origin.cliqueBeamUniformSampleAmount < 1.0) {
	numToSample = (unsigned)(origin.cliqueBeamUniformSampleAmount*(double)numEntriesPruned);
      } else { // > 1.0
	numToSample = (unsigned)origin.cliqueBeamUniformSampleAmount;
      }

      numToSample = min(numToSample,numEntriesPruned);
      if (numToSample == 0) {
	; // do nothing
      } else if (numToSample == numEntriesPruned) {
	numCliqueValuesUsed = origNumCliqueValuesUsed;
      } else {
	while (numToSample > 0) {

	  unsigned entry = rnd.uniform(--numEntriesPruned);

	  // swap the entry to the end of the current clique.
	  swap(cliqueValues[numCliqueValuesUsed],cliqueValues[numCliqueValuesUsed + entry]);
    
	  numToSample --;
	  numCliqueValuesUsed++;
	}
      }
    }
  }

  infoMsg(IM::Inference, IM::Med,"Clique uniform sampling: Upped state space from %d to %d, before pruning state space was %d\n",
	  numCliqueValuesUsedBeforeSampling,numCliqueValuesUsed,origNumCliqueValuesUsed);

}




/*-
 *-----------------------------------------------------------------------
 * PedagogicalCliqueTable::sumProbabilities()
 *
 *    Simply sum up the probabilities of all elements in this clique
 *    and return the results.
 *
 * Preconditions:
 *      Clique data structures must be created.
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      None
 *
 * Results:
 *     sum of probabilities of all elements in the clique.
 *
 *-----------------------------------------------------------------------
 */
logpr
PedagogicalCliqueTable::
sumProbabilities()
{
  logpr p;
  if (numCliqueValuesUsed > 0) {
    // We directly assign first one rather than adding to initialized
    // zero so that logpr's log(0) floating point value is preserved
    // (see logp.h).
    p = cliqueValues.ptr[0].p;
    for (unsigned i=1;i<numCliqueValuesUsed;i++)
      p += cliqueValues.ptr[i].p;
  }
  return p;
}



/*-
 *-----------------------------------------------------------------------
 * PedagogicalCliqueTable::sumExponentiatedProbabilities()
 *
 *    Simply sum up the exponentiated probabilities of all elements in this clique
 *    and return the results.
 *
 * Preconditions:
 *      Clique data structures must be created.
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      None
 *
 * Results:
 *     sum of probabilities of all elements in the clique.
 *
 *-----------------------------------------------------------------------
 */
logpr
PedagogicalCliqueTable::
sumExponentiatedProbabilities(double exponent,
			      CliqueValue* curCliqueVals,
			      const unsigned curNumCliqueValuesUsed)
{
  logpr p;
  if (curNumCliqueValuesUsed > 0) {
    // We directly assign first one rather than adding to initialized
    // zero so that logpr's log(0) floating point value is preserved.
    p = curCliqueVals[0].p.pow(exponent);
    for (unsigned i=1;i<curNumCliqueValuesUsed;i++)
      p += curCliqueVals[i].p.pow(exponent);
  }
  return p;
}



/*-
 *-----------------------------------------------------------------------
 * PedagogicalCliqueTable::cliqueEntropy()
 *
 *    Compute the clique entropy (base 2) of the forcibly normalized clique.
 *
 * Preconditions:
 *      Clique data structures must be created.
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      None
 *
 * Results:
 *     sum of probabilities of all elements in the clique.
 *
 *-----------------------------------------------------------------------
 */
double
PedagogicalCliqueTable::
cliqueEntropy()
{
  logpr sum = sumProbabilities();
  double H = 0.0;
  if (numCliqueValuesUsed > 0) {
    logpr tmp = cliqueValues.ptr[0].p/sum;
    H = tmp.unlog() * tmp.val();
    for (unsigned i=1;i<numCliqueValuesUsed;i++) {
      logpr tmp = cliqueValues.ptr[i].p/sum;
      H += tmp.unlog() * tmp.val();
    }
  }
  // convert to entropy and log base 2.
  return - H / logpr::internal_log(2.0);
}






/*-
 *-----------------------------------------------------------------------
 * PedagogicalCliqueTable::reportMemoryUsageTo()
 *
 *    Report current memory usage of this inference clique *and* the origin clique to the file in
 *    units of MBs.
 *
 * Preconditions:
 *      Clique data structures must be created.
 *
 * Postconditions:
 *      results printed.
 *
 * Side Effects:
 *      None
 *
 * Results:
 *     none
 *
 *-----------------------------------------------------------------------
 */
void
PedagogicalCliqueTable::
reportMemoryUsageTo(MaxClique& origin,FILE *f)
{
  // Memory: IC=Inference Clique, CT=Clique Table
  fprintf(f,"*MEM:IC CT(used=%lu,all=%lu=%luMB), ",
	  (unsigned long)numCliqueValuesUsed,
	  (unsigned long)cliqueValues.size(),
	  (unsigned long)((1+(sizeof(CliqueValue)*(unsigned long)cliqueValues.size())/(1024ul*1024ul))));
  origin.reportMemoryUsageTo(f);
  fprintf(f,"\n");
}




/*-
 *-----------------------------------------------------------------------
 * PedagogicalCliqueTable::maxProbability()
 *
 *    This routine just computes the max clique probabiltiy and
 *    optionally sets the clique to the clique value to the one
 *    that has maximum score. It sets all RVs to the values associated
 *    with the clique value that has maximum score.
 *    
 * Preconditions:
 *
 *     Clique table must be at least partially instantiated.
 *
 * Postconditions:
 *
 *     Hidden RVs now have values of clique value having maximum score
 *     in clique table.
 *
 * Side Effects:
 *
 *    Changes values of RVs associated with this clique
 *
 * Results:
 *    nothing
 *
 *-----------------------------------------------------------------------
 */
logpr
PedagogicalCliqueTable::
maxProbability(PedagogicalCliqueTable::SharedLocalStructure& sharedStructure,
	       bool setCliqueToMaxValue)
{
  // syntactic convenience variables.
  MaxClique& origin = *(sharedStructure.origin);

  // check for empty clique and if so, return zero.
  if (numCliqueValuesUsed == 0)
    return logpr();
  
  if (origin.hashableNodes.size() == 0) {
    // The observed clique case requires no action since this
    // means that the cliuqe (and therefore all its separators)
    // are all observed and already set to their max prob (and only) values.
    return cliqueValues.ptr[0].p;
  } else {
    unsigned max_cvn = 0;
    logpr max_cvn_score = cliqueValues.ptr[0].p;
    
    // find the max score clique entry
    for (unsigned cvn=1;cvn<numCliqueValuesUsed;cvn++) {
      if (cliqueValues.ptr[cvn].p > max_cvn_score) {
	max_cvn_score = cliqueValues.ptr[cvn].p;
	max_cvn = cvn;
      }
    }

    if (setCliqueToMaxValue) {
      // store the current table entry for the max clique.
      back_max_cvn = max_cvn;

      const bool imc_nwwoh_p = (origin.packer.packedLen() <= IMC_NWWOH); 
      if (imc_nwwoh_p) {
	origin.packer.unpack((unsigned*)&(cliqueValues.ptr[max_cvn].val[0]),
			     (unsigned**)sharedStructure.discreteValuePtrs.ptr);
      } else {
	origin.packer.unpack((unsigned*)cliqueValues.ptr[max_cvn].ptr,
			     (unsigned**)sharedStructure.discreteValuePtrs.ptr);
      }
      for (unsigned j=0;j<sharedStructure.fDeterminableNodes.size();j++) {
	RV* rv = sharedStructure.fDeterminableNodes.ptr[j];
	RV2DRV(rv)->assignDeterministicChild();
      }

    }

    return max_cvn_score;
  }
}
logpr
PedagogicalCliqueTable::
maxProb()
{
  if (numCliqueValuesUsed == 0)
    return logpr();

#if 1
  // check for empty clique and if so, return zero.
  logpr mx = cliqueValues.ptr[0].p;
  // find the max score clique entry
  for (unsigned cvn=1;cvn<numCliqueValuesUsed;cvn++) {
    if (cliqueValues.ptr[cvn].p > mx) {
      mx = cliqueValues.ptr[cvn].p;
    }
  }
  return mx;
#else
  // pipeline version exposing independent ops,
  // since max is bad for ilp.

  // check for empty clique and if so, return zero.
  REGISTER logpr mx0 = cliqueValues.ptr[0].p;
  if (numCliqueValuesUsed == 1)
    return mx0;
  REGISTER logpr mx1 = cliqueValues.ptr[1].p;
  // find the max score clique entry
  unsigned end_loc = numCliqueValuesUsed & ~0x1;
  for (unsigned cvn=2;cvn<end_loc;cvn+=2) {
    if (cliqueValues.ptr[cvn].p > mx0) {
      mx0 = cliqueValues.ptr[cvn].p;
    }
    if (cliqueValues.ptr[cvn+1].p > mx1) {
      mx1 = cliqueValues.ptr[cvn+1].p;
    }
  }
  if (numCliqueValuesUsed & 0x1) {
    if (cliqueValues.ptr[numCliqueValuesUsed-1].p > mx1) {
      mx1 = cliqueValues.ptr[numCliqueValuesUsed-1].p;
    }
  }
  if (mx0 > mx1)
    return mx0;
  else
    return mx1;
#endif

}



/*-
 *-----------------------------------------------------------------------
 * PedagogicalCliqueTable::cliqueDiversity()
 *
 *    Compute the clique diversity of the clique.
 *
 * Preconditions:
 *      Clique data structures must be created. In the default version
 *      of this routine, all of the clique entries must have been
 *      removed from the temporary local clique value pool and placed
 *      in the permanent arrays before calling this routine. But you
 *      can, by change the #define below if you want to, use this
 *      routine for debugging before the tmp entries have been copied
 *      out (which might be useful to print before pruning has
 *      occured).
 *
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      None
 *
 * Results:
 *     sum of probabilities of all elements in the clique.
 *
 *-----------------------------------------------------------------------
 */
double
PedagogicalCliqueTable::
cliqueDiversity(MaxClique& origin)
{
  double D = 0.0;
  if (numCliqueValuesUsed > 0) {
    for (unsigned i=1;i<numCliqueValuesUsed;i++) {
      for (unsigned j=i;j<numCliqueValuesUsed;j++) {
	unsigned* i_key_p; 
	unsigned* j_key_p; 
#if defined(USE_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL) && defined(CLIQUE_DIVERSITY_USING_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL)
	// then we print the entries that are currently stored in the temporary clique value pool.
	// This might be useful for debugging to print the entries before they are copied into their
	// permanent locations (i.e., if we wish to print before pruning has occured).
	i_key_p =
	  &origin.temporaryCliqueValuePool.ptr[cliqueValues.ptr[i].ival];
	j_key_p =
	  &origin.temporaryCliqueValuePool.ptr[cliqueValues.ptr[j].ival];
#else
	i_key_p = (unsigned*)cliqueValues.ptr[i].ptr;
	j_key_p = (unsigned*)cliqueValues.ptr[j].ptr;
#endif	
	double cur_dist = 
	  origin.packer.use_distance(i_key_p,j_key_p);
	D += cur_dist;
      }
    }
  }
  D = D /(  ((double)numCliqueValuesUsed)* ((double)(numCliqueValuesUsed+1))/2.0 );
  return D;
}



/*-
 *-----------------------------------------------------------------------
 * PedagogicalCliqueTable::printCliqueEntries()
 *
 *    Simply prints out all elements in the clique, giving both the RV
 *    values and the (not necessarily a probability) score for the
 *    corresponding clique entry. If the normalize option is given,
 *    then the scores can rightfully be interpreted as probabilities,
 *    and if CE/DE has been called, then this will produce true marginal
 *    probabilities over the variables in the clique.
 *
 * Preconditions:
 *
 *      Clique data structures must be created. In the default version
 *      of this routine, all of the clique entries must have been
 *      removed from the temporary local clique value pool and placed
 *      in the permanent arrays before calling this routine. But you
 *      can, by change the #define below if you want to, use this
 *      routine for debugging before the tmp entries have been copied
 *      out (which might be useful to print before pruning has
 *      occured).
 *
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      None
 *
 * Results:
 *     none
 *
 *-----------------------------------------------------------------------
 */
void
PedagogicalCliqueTable::
printCliqueEntries(PedagogicalCliqueTable::SharedLocalStructure& sharedStructure,
		   FILE *f,const char*str, 
		   const bool normalize, const bool unlog,
		   const bool justPrintEntropy)
{
  MaxClique& origin = *(sharedStructure.origin);

  fprintf(f,"--------\n");
  if (str != NULL)
    fprintf(f,"%s ",str);
  
  // TODO: might also want to optionally print the cliqueDiversity() of this
  // clique as well as the entropy.

  fprintf(f,"Printing Clique with %d variables, %d entries, H=%e\n",
	  sharedStructure.fNodes.size(),numCliqueValuesUsed,cliqueEntropy());
  if (justPrintEntropy)
    return;
  logpr sum;
  if (normalize)
    sum = sumProbabilities();
  const bool imc_nwwoh_p = (origin.packer.packedLen() <= IMC_NWWOH);
  const bool clique_has_hidden_vars = (origin.hashableNodes.size() > 0);
  for (unsigned cvn=0;cvn<numCliqueValuesUsed;cvn++) {

    if (clique_has_hidden_vars) {
      if (imc_nwwoh_p) {
	origin.packer.unpack((unsigned*)&(cliqueValues.ptr[cvn].val[0]),
			     (unsigned**)sharedStructure.discreteValuePtrs.ptr);
      } else {

#if defined(USE_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL) && defined(PRINT_CLIQUE_ENTRIES_USING_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL)
	// then we print the entries that are currently stored in the temporary clique value pool.
	// This might be useful for debugging to print the entries before they are copied into their
	// permanent locations (i.e., if we wish to print before pruning has occured).
	origin.packer.unpack((unsigned*)&origin.temporaryCliqueValuePool.ptr[cliqueValues.ptr[cvn].ival],
			     (unsigned**)sharedStructure.discreteValuePtrs.ptr);
#else
        // print the entries assuming they are stored in their final locations.
	origin.packer.unpack((unsigned*)cliqueValues.ptr[cvn].ptr,
			     (unsigned**)sharedStructure.discreteValuePtrs.ptr);
#endif
      }
      for (unsigned j=0;j<sharedStructure.fDeterminableNodes.size();j++) {
	RV* rv = sharedStructure.fDeterminableNodes.ptr[j];
	RV2DRV(rv)->assignDeterministicChild();
      }
    }
    if (normalize) {
      if (unlog) {
	// then print the exponentiated probability
	fprintf(f,"%d: %.8e ",cvn,(cliqueValues.ptr[cvn].p/sum).unlog());
      } else {
	fprintf(f,"%d: %.8e ",cvn,(cliqueValues.ptr[cvn].p/sum).valref());
      }
    } else {
      if (unlog) {
	fprintf(f,"%d: %.8e ",cvn,cliqueValues.ptr[cvn].p.unlog());
      } else {
	// print the log value directly
	fprintf(f,"%d: %.8e ",cvn,cliqueValues.ptr[cvn].p.valref());
      }
    }
    printRVSetAndValues(f,sharedStructure.fNodes);
    fflush(f);
  }
}


#define MAX_CLIQUE_DOMAIN_SIZE (0x7FFFFFFF)

unsigned
PedagogicalCliqueTable::
cliqueDomainSize(PedagogicalCliqueTable::SharedLocalStructure& sharedStructure) {
  if (sharedStructure.discreteValuePtrs.size() == 0) 
    return 0;
  unsigned domain = 1;
  for (unsigned i = 0; i < sharedStructure.fNodes.size(); i += 1) {
    if (sharedStructure.fNodes[i]->discrete() && sharedStructure.fNodes[i]->hidden()) {
      domain *= RV2DRV(sharedStructure.fNodes[i])->cardinality;
    }
  }
#if 0
printf("cliqueDomainSize = %u :", domain);
for (unsigned i = 0; i < sharedStructure.fNodes.size(); i += 1) {
  if (sharedStructure.fNodes[i]->discrete() && sharedStructure.fNodes[i]->hidden()) {
    DiscRV *drv = RV2DRV(sharedStructure.fNodes[i]);
    printf(" %s(%u)=%u", drv->name().c_str(), drv->frame(), drv->cardinality);
  }
}
printf("\n");
#endif
  return domain;
}


void 
PedagogicalCliqueTable::
printCliqueOrder(FILE *f, SharedLocalStructure& sharedStructure, int frameDelta) {
  if (sharedStructure.discreteValuePtrs.size() == 0) 
    return ;
  for (unsigned i = 0; i < sharedStructure.fNodes.size(); i += 1) {
    if (sharedStructure.fNodes[i]->discrete() && sharedStructure.fNodes[i]->hidden()) {
      DiscRV *drv = RV2DRV(sharedStructure.fNodes[i]);
      fprintf(f, " %s(%d)", drv->name().c_str(), (int)drv->frame() + frameDelta);
    }
  }
  fprintf(f,"\n");
}


unsigned
PedagogicalCliqueTable::
cliqueValueMagnitude(PedagogicalCliqueTable::SharedLocalStructure& sharedStructure, unsigned cliqueIndex) {

  MaxClique& origin = *(sharedStructure.origin);
  const bool imc_nwwoh_p = (origin.packer.packedLen() <= IMC_NWWOH);
  const bool clique_has_hidden_vars = (origin.hashableNodes.size() > 0);
  unsigned nDigits = sharedStructure.discreteValuePtrs.size();
  unsigned magnitude = 0;

  if (nDigits == 0) return 0;
  assert(cliqueIndex < cliqueValues.size());

  vector<unsigned> base (nDigits);
  vector<DiscRV *> digitRVs(nDigits);

  for (unsigned i = 0, digit = 0; i < sharedStructure.fNodes.size(); i += 1) {
    if (sharedStructure.fNodes[i]->discrete() && sharedStructure.fNodes[i]->hidden()) {
      digitRVs[digit++] = RV2DRV(sharedStructure.fNodes[i]);
    }
  }
  base[nDigits-1] = 1;
  for (int i = nDigits-2; i >= 0; i -= 1) {
    unsigned card = digitRVs[i+1]->cardinality;
    if (card >= MAX_CLIQUE_DOMAIN_SIZE / base[i+1]) {
      error("ERROR: clique state space too large for non-sparse printing\n");
    }
    base[i] = base[i+1] * card;
  }
  if (digitRVs[0]->cardinality >= MAX_CLIQUE_DOMAIN_SIZE / base[0]) {
    error("ERROR: clique state space too large for non-sparse printing\n");
  }

  if (clique_has_hidden_vars) {
    if (imc_nwwoh_p) {
      origin.packer.unpack((unsigned*)&(cliqueValues.ptr[cliqueIndex].val[0]),
			   (unsigned**)sharedStructure.discreteValuePtrs.ptr);
    } else {
#if defined(USE_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL) && defined(PRINT_CLIQUE_ENTRIES_USING_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL)
      // the entries are currently stored in the temporary clique value pool.
      // This might be useful for debugging to print the entries before they are copied into their
      // permanent locations (i.e., if we wish to print before pruning has occured).
      origin.packer.unpack((unsigned*)&origin.temporaryCliqueValuePool.ptr[cliqueValues.ptr[cliqueIndex].ival],
			   (unsigned**)sharedStructure.discreteValuePtrs.ptr);
#else
      // assume entries are stored in their final locations.
      origin.packer.unpack((unsigned*)cliqueValues.ptr[cliqueIndex].ptr,
			   (unsigned**)sharedStructure.discreteValuePtrs.ptr);
#endif
    }
    for (unsigned j=0;j<sharedStructure.fDeterminableNodes.size();j++) {
      RV* rv = sharedStructure.fDeterminableNodes.ptr[j];
      RV2DRV(rv)->assignDeterministicChild();
    }
    for (unsigned i=0; i < nDigits; i+=1) {
      magnitude += digitRVs[i]->val * base[i];
    }
  }
#if 0
printf("cliqueValueMagnitude(%u) = %u :  %s(%u)=%u", cliqueIndex, magnitude,
       digitRVs[nDigits-1]->name().c_str(), digitRVs[nDigits-1]->frame(), digitRVs[nDigits-1]->val);
for (int i = nDigits-2; i >= 0; i -= 1) {
  printf("  %s(%u)=%u", digitRVs[i]->name().c_str(), digitRVs[i]->frame(), digitRVs[i]->val);
}
printf("\n");
#endif
  return magnitude;
}


// a - b; so positive iff a > b       == 0 if a == b      negative iff a < b
int
PedagogicalCliqueTable::
cliqueValueDistance(PedagogicalCliqueTable::SharedLocalStructure& sharedStructure, 
		    unsigned a, unsigned b) 
{
  MaxClique& origin = *(sharedStructure.origin);
  const bool imc_nwwoh_p = (origin.packer.packedLen() <= IMC_NWWOH);
  const bool clique_has_hidden_vars = (origin.hashableNodes.size() > 0);

  unsigned nDigits = sharedStructure.discreteValuePtrs.size();

  if (nDigits == 0) return 0;
  assert(a < cliqueValues.size() && b < cliqueValues.size());

  vector<unsigned> aVals(nDigits);
  vector<unsigned> bVals(nDigits);
  vector<unsigned> base (nDigits);
  vector<DiscRV *> digitRVs(nDigits);

  for (unsigned i = 0, digit = 0; i < sharedStructure.fNodes.size(); i += 1) {
    if (sharedStructure.fNodes[i]->discrete() && sharedStructure.fNodes[i]->hidden()) {
      digitRVs[digit++] = RV2DRV(sharedStructure.fNodes[i]);
    }
  }

  base[nDigits-1] = 1;
  for (int i = nDigits-2; i >= 0; i -= 1) {
    unsigned card = digitRVs[i+1]->cardinality;
    if (card >= MAX_CLIQUE_DOMAIN_SIZE / base[i+1]) {
      error("ERROR: clique state space too large for non-sparse printing\n");
    }
    base[i] = base[i+1] * card;
  }
  if (digitRVs[0]->cardinality >= MAX_CLIQUE_DOMAIN_SIZE / base[0]) {
    error("ERROR: clique state space too large for non-sparse printing\n");
  }

#if 0
  for (unsigned i=0; i < nDigits; i+=1) {
    printf("  |%s(%u)|=%u", digitRVs[i]->name().c_str(), digitRVs[i]->frame(), digitRVs[i]->cardinality);
  }
  printf("\n");

printf(" < ");
#endif

  if (clique_has_hidden_vars) {
    if (imc_nwwoh_p) {
      origin.packer.unpack((unsigned*)&(cliqueValues.ptr[a].val[0]),
			   (unsigned**)sharedStructure.discreteValuePtrs.ptr);
    } else {
#if defined(USE_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL) && defined(PRINT_CLIQUE_ENTRIES_USING_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL)
      // the entries are currently stored in the temporary clique value pool.
      // This might be useful for debugging to print the entries before they are copied into their
      // permanent locations (i.e., if we wish to print before pruning has occured).
      origin.packer.unpack((unsigned*)&origin.temporaryCliqueValuePool.ptr[cliqueValues.ptr[a].ival],
			   (unsigned**)sharedStructure.discreteValuePtrs.ptr);
#else
      // assume entries are stored in their final locations.
      origin.packer.unpack((unsigned*)cliqueValues.ptr[a].ptr,
			   (unsigned**)sharedStructure.discreteValuePtrs.ptr);
#endif
    }
    for (unsigned j=0;j<sharedStructure.fDeterminableNodes.size();j++) {
      RV* rv = sharedStructure.fDeterminableNodes.ptr[j];
      RV2DRV(rv)->assignDeterministicChild();
    }
    for (unsigned i=0; i < nDigits; i+=1) {
      aVals[i] = digitRVs[i]->val; //*(sharedStructure.discreteValuePtrs[i]);
    }

#if 0
  for (unsigned i=0; i < digitRVs.size(); i+=1) {
    printf(" %s(%u)=%u ", digitRVs[i]->name().c_str(), digitRVs[i]->frame(), digitRVs[i]->val);
  }
  printf(" > - < ");
#endif

#if 1
    if (imc_nwwoh_p) {
      origin.packer.unpack((unsigned*)&(cliqueValues.ptr[b].val[0]),
			   (unsigned**)sharedStructure.discreteValuePtrs.ptr);
    } else {
#if defined(USE_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL) && defined(PRINT_CLIQUE_ENTRIES_USING_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL)
      // the entries are currently stored in the temporary clique value pool.
      // This might be useful for debugging to print the entries before they are copied into their
      // permanent locations (i.e., if we wish to print before pruning has occured).
      origin.packer.unpack((unsigned*)&origin.temporaryCliqueValuePool.ptr[cliqueValues.ptr[b].ival],
			   (unsigned**)sharedStructure.discreteValuePtrs.ptr);
#else
      // assume entries are stored in their final locations.
      origin.packer.unpack((unsigned*)cliqueValues.ptr[b].ptr,
			   (unsigned**)sharedStructure.discreteValuePtrs.ptr);
#endif
    }
    for (unsigned j=0;j<sharedStructure.fDeterminableNodes.size();j++) {
      RV* rv = sharedStructure.fDeterminableNodes.ptr[j];
      RV2DRV(rv)->assignDeterministicChild();
    }
    for (unsigned i=0; i < nDigits; i+=1) {
      bVals[i] = digitRVs[i]->val; //*(sharedStructure.discreteValuePtrs[i]);
    }
#endif
  }

  int diff = 0;
  for (unsigned i=0; i < nDigits; i+=1) {
    diff += (aVals[i] - bVals[i]) * base[i];
  }

#if 0
  for (unsigned i=0; i < digitRVs.size(); i+=1) {
    printf(" %s(%u)=%u ", digitRVs[i]->name().c_str(), digitRVs[i]->frame(), digitRVs[i]->val);
  }
  printf(" > = %d\n", diff);
#endif

#if 0
printf("cliqueValueDistance(%u,%u) = %d\n", a, b, diff);
#endif

  return diff;
}


void
PedagogicalCliqueTable::
printCliqueEntries(PedagogicalCliqueTable::SharedLocalStructure& sharedStructure,
		   ObservationFile *f, const bool normalize, const bool unlog)
{
#if 0
  printf("clique has %d variables (%d hidden discrete), %d entries\n", 
	 sharedStructure.fNodes.size(),sharedStructure.discreteValuePtrs.size(),numCliqueValuesUsed);
#endif

  logpr sum;
  if (normalize)
    sum = sumProbabilities();
#if 0
  const bool imc_nwwoh_p = (origin.packer.packedLen() <= IMC_NWWOH);
  const bool clique_has_hidden_vars = (origin.hashableNodes.size() > 0);
#endif
  sArray<CliqueValueIndex> index(numCliqueValuesUsed);

  for (unsigned cvn=0;cvn<numCliqueValuesUsed;cvn++) {
    index[cvn] = CliqueValueIndex(&sharedStructure, this, cvn);
  }


#if 0
printf("  before:    ");
for (unsigned i=0; i < numCliqueValuesUsed; i+=1){
  printf(" %d", cliqueValueMagnitude(sharedStructure, index[i].index));
}
#endif

  index.qsort();

#if 0
printf("    sort:    ");
for (unsigned i=0; i < numCliqueValuesUsed; i+=1){
  printf(" %d", cliqueValueMagnitude(sharedStructure, index[i].index));
}
printf("\n");
#endif

  float x;

  union {
    Data32  d;
    float   R;
  } zero;
  zero.R = LZERO;
  
  union {
    Data32  d;
    float   R;
  } feature;

  // We don't need to handle printing zero cliques, since they either
  // throw an exception or terminate the program. Thus it should be
  // safe to assume at least 1 value in the clique at this point.

  unsigned skip = cliqueValueMagnitude(sharedStructure, index[0].index);
  assert(skip < cliqueDomainSize(sharedStructure));  // there must be at least 1 non-zero clique value
  for (unsigned i=0; i < skip; i+=1) {
    if (unlog) {
      f->writeFeature(0);
    } else {
      f->writeFeature(zero.d);
    }
  }
  if (normalize) {
    if (unlog) {
      // then print the exponentiated probability
      x = (cliqueValues.ptr[index[0].index].p/sum).unlog();
    } else {
      x = (cliqueValues.ptr[index[0].index].p/sum).valref();
    }
  } else {
    if (unlog) {
      x = cliqueValues.ptr[index[0].index].p.unlog();
    } else {
      // print the log value directly
      x = cliqueValues.ptr[index[0].index].p.valref();
    }
  }
  feature.R = x;
  f->writeFeature(feature.d);
  unsigned cvn = index[0].index;
  for (unsigned cvnIdx = 1; cvnIdx < numCliqueValuesUsed; cvnIdx++) {
    unsigned prev = cvn;
    cvn = index[cvnIdx].index;
    skip = (int)cliqueValueDistance(sharedStructure, cvn, prev) - 1;
    assert(0 <= skip && skip < cliqueDomainSize(sharedStructure));
    for (unsigned i=0; i < skip; i+=1) {
      if (unlog) {
        f->writeFeature(0);
      } else {
        f->writeFeature(zero.d);
      }
    }
    if (normalize) {
      if (unlog) {
	// then print the exponentiated probability
	x = (cliqueValues.ptr[cvn].p/sum).unlog();
      } else {
	x = (cliqueValues.ptr[cvn].p/sum).valref();
      }
    } else {
      if (unlog) {
	x = cliqueValues.ptr[cvn].p.unlog();
      } else {
	// print the log value directly
	x = cliqueValues.ptr[cvn].p.valref();
      }
    }
    feature.R = x;
    f->writeFeature(feature.d);
  }
  skip = cliqueDomainSize(sharedStructure) - 
    cliqueValueMagnitude(sharedStructure,index[numCliqueValuesUsed-1].index) - 1;
  assert(0 <= skip && skip < cliqueDomainSize(sharedStructure));
  for (unsigned i=0; i < skip; i+=1) {
    if (unlog) {
      f->writeFeature(0);
    } else {
      f->writeFeature(zero.d);
    }
  }
}



/*-
 *-----------------------------------------------------------------------
 * PedagogicalCliqueTable::emIncrement()
 *
 *      Do the work of EM increment. That is, for each clique element,
 *      load the values into the variables and then for each assigned prob node
 *      call that nodes emIncrement() routine thus realizing EM training.
 *      The routine has the option of using either a local normalization
 *      (useful when doing double precision 64-bit floating point) or a local
 *      normalization (perhaps better for signle precision, especially on
 *      longer utterances).
 *
 * Preconditions:
 *      Clique data structures must be created.
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      Potentially all assigned probability nodes accumulators are changed.
 *
 * Results:
 *     none
 *
 *-----------------------------------------------------------------------
 */
void
PedagogicalCliqueTable::
emIncrement(PedagogicalCliqueTable::SharedLocalStructure& sharedStructure,
	    const logpr probE,
	    const bool localCliqueNormalization,
	    const double emTrainingBeam)
{
  MaxClique& origin = *(sharedStructure.origin);

  // recompute here each time, shouldn't take too long
  // TODO: re-compute this once for each inference clique.
  sArray< RV*> fAssignedProbNodes;
  unsigned numAssignedProbNodes = 0;
  for (unsigned nodeNumber = 0; 
       nodeNumber < sharedStructure.fSortedAssignedNodes.size(); 
       nodeNumber ++ ) {
    if (origin.dispositionSortedAssignedNodes[nodeNumber] == MaxClique::AN_CPT_ITERATION_COMPUTE_AND_APPLY_PROB ||
	origin.dispositionSortedAssignedNodes[nodeNumber] == MaxClique::AN_COMPUTE_AND_APPLY_PROB) {
      numAssignedProbNodes++;
    }
  }
  if (numAssignedProbNodes == 0)
    return; // nothing to do for this clique

  infoMsg(Huge,"EM Incrementing for clique\n");

  fAssignedProbNodes.resize(numAssignedProbNodes);
  numAssignedProbNodes = 0;
  for (unsigned nodeNumber = 0; nodeNumber < 
	 sharedStructure.fSortedAssignedNodes.size(); nodeNumber ++ ) {
    if (origin.dispositionSortedAssignedNodes[nodeNumber] == MaxClique::AN_CPT_ITERATION_COMPUTE_AND_APPLY_PROB ||
	origin.dispositionSortedAssignedNodes[nodeNumber] == MaxClique::AN_COMPUTE_AND_APPLY_PROB) {
      RV* rv = sharedStructure.fSortedAssignedNodes[nodeNumber];
      fAssignedProbNodes[numAssignedProbNodes++] = rv;
    }
  }

  logpr locProbE((void*)0);
  if (localCliqueNormalization) {
    // this case is probably better/safer when using 32-bit single precision IEEE floating point for logpr.
    locProbE = sumProbabilities();
  } else {
    // probably ok to do this when using double precision IEEE floating point.
    locProbE = probE;
  }

  const bool imc_nwwoh_p = (origin.packer.packedLen() <= IMC_NWWOH);
  const bool clique_has_hidden_vars = (origin.hashableNodes.size() > 0);
  
  logpr beamThreshold((void*)0);
  if (emTrainingBeam != (-LZERO)) {
    // then we do clique table pruning right here rather
    // than a separate call to ceCliqueBeamPrune().
    // break into the logp to avoid unnecessary zero checking.
    beamThreshold.valref() = - emTrainingBeam;
  } else {
    // set beam threshold to a value that will never cause real
    // pruning (i.e., we always prune posteriors that are almost
    // zero).
    beamThreshold.set_to_almost_zero();
  }
  // create unnormalized beam threshold.
  beamThreshold *= locProbE;

  // now go through updating each thing
  for (unsigned cvn=0;cvn<numCliqueValuesUsed;cvn++) {

    // EM pruning here based on unnormalized posterior. Don't bother
    // with things that are below threshold.
    if (cliqueValues.ptr[cvn].p <= beamThreshold)
      continue;

    // if still here, then create the posterior to update the
    // parameters.
    logpr posterior = cliqueValues.ptr[cvn].p/locProbE;

    // printf("EM training, cvn=%d, log(posterior) = %f\n",cvn,posterior.valref());

    if (clique_has_hidden_vars) {
      if (imc_nwwoh_p) {
	origin.packer.unpack((unsigned*)&(cliqueValues.ptr[cvn].val[0]),
			     (unsigned**)sharedStructure.discreteValuePtrs.ptr);
      } else {
	origin.packer.unpack((unsigned*)cliqueValues.ptr[cvn].ptr,
			     (unsigned**)sharedStructure.discreteValuePtrs.ptr);
      }
      for (unsigned j=0;j<sharedStructure.fDeterminableNodes.size();j++) {
	RV* rv = sharedStructure.fDeterminableNodes.ptr[j];
	RV2DRV(rv)->assignDeterministicChild();
      }

    }
    // Increment all assigned probability nodes.  
    // 
    // TODO: integerate out all but cont. variables parents so we
    // don't multiply increment those varialbes for the same parent
    // values (waisting time). Alternatively, this can be done by
    // creating a sub-clique hanging off of this clique which contains
    // only the observed variables and its parents.
    for (unsigned nodeNumber = 0; nodeNumber < fAssignedProbNodes.size(); nodeNumber ++ ) {
      RV* rv = fAssignedProbNodes[nodeNumber];
      rv->emIncrement(posterior);
    }
  }

}



/*-
 *-----------------------------------------------------------------------
 * PedagogicalCliqueTable::deReceiveFromIncommingSeparator()
 *
 *      We are in the DE phase, and the separator that during the CE phase
 *      we sent a message out to is now ready with a backwards message back
 *      to this clique. We iterate through all clique entries and lookup
 *      the corresponding "incomming separator" entry, multiplying in its score
 *      into the corresponding clique entry.
 *
 *      Note that some of the incomming separator values might be zero (for
 *      reason based on higher up in the JT). Therefore, at the end of this
 *      routine we do another DE pruning pass to remove clique entries
 *      that have become zero. Since one separator entry might correspond to
 *      many clique entries, a zero separator entry might cause a lot of backwards
 *      clique pruning. 
 *
 *      Note: A condition might occur where we've got a clique entry without
 *      a corresponding incomming separator entry. On first thought, this
 *      shouldn't occur, but in actuallity it can occur for one of two reasons:
 *
 *     1) if sbeam is turned on, it's possible that the separator entry corresponding to this
 *        clique entry was pruned away. This can happen since separator pruning is later
 *        then the clique which created (projected down into) the separator. The solution
 *        here is to increase or entirely turn off sbeam pruning.
 *     
 *     2) During the island algorithm and using k-pruning (fixed k clique size), there
 *        were ties in the clique probability and since we use a random median to find the top
 *        k entries, we didn't get the same top k entries the 2nd (or 3rd, 4th, etc. depending 
 *        on the island algorithm's 'base' and 'lst' parameters) time that we created the 
 *        clique and then pruned it down to k entries. The solution we employ is to
 *        seed the random number generator in k-pruning to something that is clique dependent.
 *        This is now being done in PedagogicalCliqueTable::ceCliqueStatePrune() above. Search for string 
 *        K-BEAM-SEED elsewhere in this file for further information. Therefore, this 
 *        case should not occur.
 *
 *        The solution (Marked by the keyword "SEPCLIQUEZERO" below) is to just
 *        prune this clique entry away, for case 1 since case 2 will not
 *        happen. 
 *
 *        Another solution would be to remove the clique-dependent seed for case
 *        2, but pruning would then get ugly. This is because if the clique was
 *        all ties, this would mean that we prune all but the intersection of
 *        the survivors of both k-prunings, and this might be small or empty.
 *        In the best of cases, this would probably be much smaller than k, so
 *        we instead just do the clique-dependent seed which is much easier.
 *
 * Preconditions:
 *      Clique data structures and separator must be created.
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      Potentially all assigned probability nodes accumulators are changed.
 *
 * Results:
 *     none
 *
 *-----------------------------------------------------------------------
 */
void 
PedagogicalCliqueTable::
deReceiveFromIncommingSeparator(PedagogicalCliqueTable::SharedLocalStructure& sharedStructure,
				ConditionalSeparatorTable* separatorTableArray,
				ConditionalSeparatorTable::SharedLocalStructure* sepSharedStructureArray)
{
  MaxClique& origin = *(sharedStructure.origin);
  deReceiveFromIncommingSeparator(sharedStructure,
				  separatorTableArray[origin.ceSendSeparator],
				  sepSharedStructureArray[origin.ceSendSeparator]
				  );
}
void 
PedagogicalCliqueTable::
deReceiveFromIncommingSeparator(PedagogicalCliqueTable::SharedLocalStructure& sharedStructure,
				ConditionalSeparatorTable& sep,
				ConditionalSeparatorTable::SharedLocalStructure& sepSharedStructure)
{
  if (SectionScheduler::viterbiScore) {
    return deReceiveFromIncommingSeparatorViterbi(sharedStructure,
						  sep,
						  sepSharedStructure);
  }

  // syntactic convenience variables.
  MaxClique& origin = *(sharedStructure.origin);
  SeparatorClique& sepOrigin = 
    *(sepSharedStructure.origin);  

  // keep a local variable copy of this around to avoid potential dereferencing.
  ConditionalSeparatorTable::AISeparatorValue * const
    sepSeparatorValuesPtr = sep.separatorValues->ptr; 

  if (origin.hashableNodes.size() == 0) {
    // do the observed clique case up front right here so we don't
    // need to keep checking below.
    ConditionalSeparatorTable::AISeparatorValue& sv
      = sepSeparatorValuesPtr[0];
    cliqueValues.ptr[0].p *= sv.remValues.ptr[0].bp();
    return;
  }

  // allocate some temporary storage for packed separator values.
  // 128 words is *much* bigger than any possible packed clique value
  // will take on, but it is easy/fast to allocate on the stack right now.
  // unsigned packedVal[128];
  // but just in case, we assert.
  assert ((sepOrigin.hAccumulatedIntersection.size() == 0)
	  ||
	  (sepOrigin.accPacker.packedLen() < 128)
	  );
  assert ((sepOrigin.hRemainder.size() == 0) 
	  ||
	  (sepOrigin.remPacker.packedLen() < 128 )
	  );
  // If this assertion fails (at some time in the future, probably in
  // the year 2150), then it is fine to increase 128 to something larger.

  // cache check here.
  const bool imc_nwwoh_p = (origin.packer.packedLen() <= IMC_NWWOH);
  for (unsigned cvn=0;cvn<numCliqueValuesUsed;cvn++) {{

      // TODO: optimize away this conditional check.
      if (imc_nwwoh_p) {
	origin.packer.unpack((unsigned*)&(cliqueValues.ptr[cvn].val[0]),
			     (unsigned**)sharedStructure.discreteValuePtrs.ptr);
      } else {
	origin.packer.unpack((unsigned*)cliqueValues.ptr[cvn].ptr,
			     (unsigned**)sharedStructure.discreteValuePtrs.ptr);
      }
      for (unsigned j=0;j<sharedStructure.fDeterminableNodes.size();j++) {
	RV* rv = sharedStructure.fDeterminableNodes.ptr[j];
	RV2DRV(rv)->assignDeterministicChild();
      }

      // All hidden random variables now have their discrete
      // value. Get appropriate entry in the separator
      // for this clique entry.

      /*
       * There are 3 cases.
       * 1) AI exists and REM exist
       * 2) AI exists and REM doesnt exist
       * 3) AI does not exist, but REM exists
       * AI not exist and REM not exist can't occur.
       */

      unsigned accIndex;
      // TODO: optimize this check away out of loop.
      if (sepOrigin.hAccumulatedIntersection.size() > 0) {
	// an accumulated intersection exists.

	sepOrigin.accPacker.pack((unsigned**)sepSharedStructure.accDiscreteValuePtrs.ptr,
				 &CliqueBuffer::packedVal[0]);
	unsigned* accIndexp =
	  sep.iAccHashMap->find(&CliqueBuffer::packedVal[0]);

	// case SEPCLIQUEZERO, see comments in routine heading
	if ( accIndexp == NULL ) {
	  // Then separator entry got pruned away. Force prune of clique
	  // entry as well.
	  cliqueValues.ptr[cvn].p.set_to_zero();
	  goto next_iteration;
	}

	accIndex = *accIndexp;

	// TODO: optimize this check out of loop.
	if (sepSharedStructure.remDiscreteValuePtrs.size() == 0) {
	  // 2) AI exists and REM doesnt exist
	  // Then this separator is entirely covered by one or 
	  // more other separators earlier in the order.

	  // go ahead and insert it here to the 1st entry (entry 0).

	  // handy reference for readability.
	  ConditionalSeparatorTable::AISeparatorValue& sv
	    = sepSeparatorValuesPtr[accIndex];

	  // Multiply in this separator value's probability.
	  cliqueValues.ptr[cvn].p *= sv.remValues.ptr[0].bp();
	  // done
	  goto next_iteration;
	}
      } else {
	// no accumulated intersection exists, everything
	// is in the remainder.
	accIndex = 0;
      }

      if (sepSharedStructure.remDiscreteValuePtrs.size() > 0) {
	// if we're here, then we must have some remainder pointers.
	// Do the remainder exists in this separator.
	// either:
	//   1) AI exists and REM exist
	//     or
	//   3) AI does not exist (accIndex == 0), but REM exists
	// 

	// keep handy reference for readability.
	ConditionalSeparatorTable::AISeparatorValue& sv
	  = sepSeparatorValuesPtr[accIndex];

	sepOrigin.remPacker.pack((unsigned**)sepSharedStructure.remDiscreteValuePtrs.ptr,
				 &CliqueBuffer::packedVal[0]);

	unsigned* remIndexp =
	  sv.iRemHashMap.find(&CliqueBuffer::packedVal[0]);

	// case SEPCLIQUEZERO, see comments in routine heading
	if ( remIndexp == NULL ) {
	  // Then separator entry got pruned away. Force prune of clique
	  // entry as well.
	  cliqueValues.ptr[cvn].p.set_to_zero();
	  goto next_iteration;
	  // warning("ERROR: During distribute evidence, found clique entry without corresponding incomming separator entry. Try increasing -sbeam beam, or just use -cbeam without the -sbeam option.\nClique contains:");
	  // printRVSetAndValues(stdout,fNodes);

	}

	// We've finally got the sep entry. Multiply it it into the
	// current clique value.
	cliqueValues.ptr[cvn].p *= sv.remValues.ptr[*remIndexp].bp();
      } else {
	// Either separator is all observed, or the separator
	// is completely contained in the accumulated intersection.
	// In either case, we multiply in its one value.

	ConditionalSeparatorTable::AISeparatorValue& sv
	  = sepSeparatorValuesPtr[accIndex];

	if (sv.numRemValuesUsed == 1) {
	  // We've finally got the sep entry. Multiply it it into the
	  // current clique value.
	  cliqueValues.ptr[cvn].p *= sv.remValues.ptr[0].bp();
	} else {
	  // case SEPCLIQUEZERO, see comments in routine heading
	  // Then separator entry got pruned away. Force prune of clique
	  // entry as well.
	  cliqueValues.ptr[cvn].p.set_to_zero();
	  goto next_iteration;
	}
      }

    }
  next_iteration:
    ;    
  }

  // Backwards pruning: Fixed backwards/distribute evidence beam
  // pruning here. Since we know log(prob(E)), we can do fairly
  // accurate pruning now. This will be useful particulalry when many
  // of the bp() values are zero (i.e., in this case, we have zero
  // compression). We just prune out the zeros for now.  
  // TODO: integrate this pruning into the above loop instead.  
  // TODO: export backwards beam width to command line & integrate
  //    with -ebeam
  {
    const unsigned origNumCliqueValuesUsed = numCliqueValuesUsed;
    for (unsigned cvn=0;cvn<numCliqueValuesUsed;) {
      if (cliqueValues.ptr[cvn].p.essentially_zero()) {

	// copy the last entry (that might be good) to the current
	// position (which has just been pruned). Do an assignment
	// since the pruned entries are never going to be needed again
	// (they're deallocated below).
	cliqueValues.ptr[cvn] = cliqueValues.ptr[--numCliqueValuesUsed];

	// alternatively, we could swap with last entry, and decrease
	// numCliqueValuesUsed by one (if for some reason we want to
	// possibly use these zero entries someday).
	// swap(cliqueValues.ptr[cvn],cliqueValues.ptr[numCliqueValuesUsed-1]);
	// numCliqueValuesUsed--;

      } else {
	cvn++;
      }
    }
    infoMsg(IM::Inference, IM::High-1,"DE Clique Receive: (old,new) clique state space = (%d,%d).\n",
	    origNumCliqueValuesUsed,numCliqueValuesUsed);
    // TODO: resize only if size difference is large.
    if (numCliqueValuesUsed < origNumCliqueValuesUsed)
      cliqueValues.resizeAndCopy(numCliqueValuesUsed);
  }

}




/*-
 *-----------------------------------------------------------------------
 * PedagogicalCliqueTable::deReceiveFromIncommingSeparatorViterbi()
 *
 *      We are in the DE phase, and the separator that during the CE
 *      phase we sent a message out to is now ready with a backwards
 *      message back to this clique. 
 * 
 *      This routine here is the Viterbi version of
 *      deReceiveFromIncommingSeparator() above.  In this case, we
 *      simply look up the separators back pointer index and use it to
 *      choose the current clique entry. Moreover, we unpack that
 *      clique entry leaving the RVs assigned to what is given in the
 *      clique entry.
 *
 *      An alternative strategy would be to assume that the RVs of the
 *      separator currently point to (i.e., hold the value of) the
 *      viterbi entry containing the back pointer index. We would then
 *      look up the separator entry corresponding to the currently set
 *      RV values via hashing, and choose the current clique entry
 *      accordingly. This, however, interacts poorly with the island
 *      algorithm, and would not generalize to n-best, so we do not do
 *      this here.
 *
 *      Yet another strategy would be to look up the separator
 *      backpointer index and use it to choose a clique entry and then
 *      store the index of that clique entry in the clique (rather
 *      than using the random variable values) since this would
 *      greatly simplify things like the island algorithm and shifting
 *      (but would require a bit more memory).
 *
 *
 * Preconditions:
 *      Clique data structures and separator must be created.
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      Potentially all assigned probability nodes accumulators are changed.
 *
 * Results:
 *     none
 *
 *-----------------------------------------------------------------------
 */
void 
PedagogicalCliqueTable::
deReceiveFromIncommingSeparatorViterbi(PedagogicalCliqueTable::SharedLocalStructure& sharedStructure,
				       ConditionalSeparatorTable& sep,
				       ConditionalSeparatorTable::SharedLocalStructure& sepSharedStructure)
{
  MaxClique& origin = *(sharedStructure.origin);

  if (origin.hashableNodes.size() == 0) {
    // All clique values already set to their max (and only) settings,
    // so no need to do anything.
    return;
  }

  // keep a local variable copy of this around to avoid potential dereferencing.
  ConditionalSeparatorTable::AISeparatorValue * const
    sepSeparatorValuesPtr = sep.separatorValues->ptr; 

  // cache check here.
  const bool imc_nwwoh_p = (origin.packer.packedLen() <= IMC_NWWOH);
  // grab backpointer from seps forward pointer entry directly.
  ConditionalSeparatorTable::AISeparatorValue& sv
    = sepSeparatorValuesPtr[sep.forwPointer.viterbiAccIndex];
  unsigned cvn = sv.remValues.ptr[sep.forwPointer.viterbiRemIndex].backPointer;

  // store the current table entry for the max clique.
  back_max_cvn = cvn;

  // unpack clique value 'cvn' into corresponding random variables and expand
  // any deterministic values.
  if (imc_nwwoh_p) {
    origin.packer.unpack((unsigned*)&(cliqueValues.ptr[cvn].val[0]),
			 (unsigned**)sharedStructure.discreteValuePtrs.ptr);
  } else {
    origin.packer.unpack((unsigned*)cliqueValues.ptr[cvn].ptr,
			 (unsigned**)sharedStructure.discreteValuePtrs.ptr);
  }
  for (unsigned j=0;j<sharedStructure.fDeterminableNodes.size();j++) {
    RV* rv = sharedStructure.fDeterminableNodes.ptr[j];
    RV2DRV(rv)->assignDeterministicChild();
  }

  // printf("*** max RV of clique values set from parent separator ***\n");
  // printRVSetAndValues(stdout,sharedStructure.fNodes,true);
  // fflush(stdout);
}



/*-
 *-----------------------------------------------------------------------
 * PedagogicalCliqueTable::deScatterToOutgoingSeparators()
 *
 *      We are in the DE phase, and we've received the separator
 *      message fror the separator we sent a message out during CE.
 *      Now it is our turn to send (project down to) to the separators 
 *      that we received incomming messages during the CE phase.
 *
 *      Specifically, we "scatter" out to the now outgoing separators which are the same as the
 *      "incomming" separators in the collect evidence stage, so we use that
 *      array here directly here. We do this by iterating through all clique entries
 *      and sending it back out to all CE-incomming separators (unless the separator
 *      is a VE separator in which case it is just a constaint that
 *      we have already acounted for).
 *
 * Preconditions:
 *      Clique data structures and separator must be created.
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      Potentially all assigned probability nodes accumulators are changed.
 *
 * Results:
 *     none
 *
 *-----------------------------------------------------------------------
 */
void 
PedagogicalCliqueTable::
deScatterToOutgoingSeparators(PedagogicalCliqueTable::SharedLocalStructure& sharedStructure,
			      ConditionalSeparatorTable* separatorTableArray,
			      ConditionalSeparatorTable::SharedLocalStructure* sepSharedStructureArray)
{

  if (SectionScheduler::viterbiScore) {
    // While we might think that there is nothing to do in this case
    // since if the RVs associated with the current clique have been
    // set to the appropriate clique table entry, the associated
    // separator RVs have also been set. But in some cases (namely the
    // island algorithm) it is possible for the separator to get
    // changed by an additional forward pass in which case we would be
    // using the wrong value.
    // 
    // Specifically, to get island decoding working, we need to have
    // separator keep track of which entry is current max rather than
    // having it assume that its values are set from the 
    // clique that, on CE, gathers in that separator and, on DE,
    // scatters out.
    // 
    // Consider the following example:
    //  
    //    The vertical bars mark the partition boundaries 
    //    of partitions P0, P1, ...
    //
    //        |           |           |   island   |  ...              
    //  P0    |     P1    |     P2    |     P3     |  ...              
    //        |           |           |            |  ...
    //    C0 -- s0 -- C1 -- s1 -- C2 -- s2 -- C3 --
    //             1>    2>    3>    <4          <b  
    //                         <5  
    // 
    //  P3 is an island in an island algorithm. We have gone through
    //  and created and deleted partitions P0, P1, P2, and P3, but we
    //  store and save P3 as an island.  We then move to the right of
    //  P3, and come back, but we then need to reconstruct partitions
    //  P0 - P3. In other words:
    // 
    //  1> and 3> are ceGatherIntoRoot msgs
    //  2> is a ceSendForwardsCrossPartitions msg
    //  <4 is a deReceveFromIncommingSep msg
    //      and <5 is a deScatterOutofRoot, but has not yet happened.
    //  <b is a deReceveFromIncommingSep msg from before (thus
    //     marked with a 'b').
    // 
    // The potential problem is that <4 might assume that it's 's2' is currently
    // set to appropriate values (max of C3 done by <b), but message 3>
    // changes clique C2 and thus changes s2 so that it no longer
    // necessarily holds max value, since s2 holds whever is 
    // currently set by 3>.
    // 
    // One solution (that adds memory) is to add two variables in
    // InferenceSep, i.e., two indices to get at the Viterbi values
    // for current separator.
    //
    //      unsigned viterbiAccIndex;
    //      unsigned viterbiRemIndex;
    //
    // to keep track of which entries in seps (such as s2) are max.
    // Then <4 would restore those values for each separator.  and
    // deScatterOutofRoot would compute them (since deScatterOutofRoot
    // is called at the time C3 is at its true max value). Note that
    // this works between partitions since the separator between
    // two partitions is contained in the right partition.
    //
    // Note also that this can probably be easily extended to the
    // N-best list case by having multiple index pairs.
    deScatterToOutgoingSeparatorsViterbi(sharedStructure,separatorTableArray,sepSharedStructureArray);
    return; 
  }
  MaxClique& origin = *(sharedStructure.origin);

  if (origin.ceReceiveSeparators.size() == 0)
    return;


  // Note. All separator .bp values have already been initialized to
  // zero when the structure containing them 'a RemainderValue' was
  // constructed. All memory reallocations will have preserved these
  // initializations, so there is no need to scan through initializing
  // bp to zero here.

  if (origin.hashableNodes.size() == 0) {
    // printf("<<<<<==== Observed clique in backwards pass\n");

    // Do the observed clique case up front right here so we don't
    // need to keep checking below. Here, the clique is observed which
    // means that all connecting separators are also observed. We just
    // sweep through all separators updating the values.
    for (unsigned sepNumber=0;sepNumber<origin.ceReceiveSeparators.size();sepNumber++) {
      ConditionalSeparatorTable& sep = 
	separatorTableArray[origin.ceReceiveSeparators[sepNumber]];
      ConditionalSeparatorTable::AISeparatorValue& sv
	= sep.separatorValues->ptr[0];
      SeparatorClique& sepOrigin = 
	*(sepSharedStructureArray[origin.ceReceiveSeparators[sepNumber]].origin);

      // don't distribute to VE separators or to one that is being skipped.
      if (sep.veSeparator() || sepOrigin.skipMe)
	continue;


      // can use assignment rather than += here since there is only one value.
      sv.remValues.ptr[0].bp() = cliqueValues.ptr[0].p;      
      sv.numRemValuesUsed = 1;
    }
  } else {

    // Allocate some temporary storage for packed separator values.
    // 128 words is *much* bigger than any possible packed clique value
    // will take on, but it is easy/fast to allocate on the stack right now.
    // unsigned packedVal[128];

    infoMsg(IM::Inference, IM::High-1,"DE Clique state space = %d.\n",numCliqueValuesUsed);

    // cache check here.
    const bool imc_nwwoh_p = (origin.packer.packedLen() <= IMC_NWWOH); 
    for (unsigned cvn=0;cvn<numCliqueValuesUsed;cvn++) {

      // TODO: optimize away this conditional check.
      if (imc_nwwoh_p) {
	origin.packer.unpack((unsigned*)&(cliqueValues.ptr[cvn].val[0]),
			     (unsigned**)sharedStructure.discreteValuePtrs.ptr);
      } else {
	origin.packer.unpack((unsigned*)cliqueValues.ptr[cvn].ptr,
			     (unsigned**)sharedStructure.discreteValuePtrs.ptr);
      }
      for (unsigned j=0;j<sharedStructure.fDeterminableNodes.size();j++) {
	RV* rv = sharedStructure.fDeterminableNodes.ptr[j];
	RV2DRV(rv)->assignDeterministicChild();
      }

      // now we iterate through all the separators.
      for (unsigned sepNumber=0;sepNumber<origin.ceReceiveSeparators.size();sepNumber++) {

	// get a handy reference to the current separator
	ConditionalSeparatorTable& sep = 
	  separatorTableArray[origin.ceReceiveSeparators[sepNumber]];
	ConditionalSeparatorTable::SharedLocalStructure& sepSharedStructure = 
	  sepSharedStructureArray[origin.ceReceiveSeparators[sepNumber]];
	SeparatorClique& sepOrigin = 
	  *(sepSharedStructure.origin);

	// don't distribute to VE separators or to one that is being skipped.
	if (sep.veSeparator() || sepOrigin.skipMe)
	  continue;

	// keep a local variable copy of this around to avoid potential dereferencing.
	ConditionalSeparatorTable::AISeparatorValue * const
	  sepSeparatorValuesPtr = sep.separatorValues->ptr; 

	// If these assertions fail (at some time in the future, probably in
	// the year 2150), then it is fine to increase 128 to something larger.
	// In fact, 128 is so large, lets not even do the assert.
	// assert ( sepOrigin.accPacker.packedLen() < 128 );
	// assert ( sepOrigin.remPacker.packedLen() < 128 );


	/*
	 * There are 3 cases.
	 * 1) AI exists and REM exist
	 * 2) AI exists and REM doesnt exist
	 * 3) AI does not exist, but REM exists
	 * AI not exist and REM not exist can't occur.
	 */

	unsigned accIndex;
	// TODO: optimize this check away out of loop.
	if (sepOrigin.hAccumulatedIntersection.size() > 0) {
	  // an accumulated intersection exists.

	  sepOrigin.accPacker.pack((unsigned**)sepSharedStructure.accDiscreteValuePtrs.ptr,
				   &CliqueBuffer::packedVal[0]);
	  unsigned* accIndexp =
	    sep.iAccHashMap->find(&CliqueBuffer::packedVal[0]);

	  // we should always find something or else something is wrong.
	  assert ( accIndexp != NULL ); 
	  accIndex = *accIndexp;

	  // TODO: optimize this check out of loop.
	  if (sepSharedStructure.remDiscreteValuePtrs.size() == 0) {
	    // 2) AI exists and REM doesnt exist
	    // Then this separator is entirely covered by one or 
	    // more other separators earlier in the order.

	    // go ahead and insert it here to the 1st entry (entry 0).

	    // handy reference for readability.
	    ConditionalSeparatorTable::AISeparatorValue& sv
	      = sepSeparatorValuesPtr[accIndex];

	    // Add in this clique value's probability.  Note that bp was
	    // initialized during forward pass.
	    sv.remValues.ptr[0].bp() += cliqueValues.ptr[cvn].p;
	    // done, move on to next separator.
	    continue; 
	  } // else, we continue on below.
	} else {
	  // no accumulated intersection exists, everything
	  // is in the remainder.
	  accIndex = 0;
	}

	if (sepSharedStructure.remDiscreteValuePtrs.size() > 0) {
	  // if we're here, then we must have some remainder
	  // pointers.

	  // Do the remainder exists in this separator.
	  // 
	  // either:
	  //   1) AI exists and REM exist
	  //     or
	  //   3) AI does not exist (accIndex == 0), but REM exists
	  // 
	
	  // keep handy reference for readability.
	  ConditionalSeparatorTable::AISeparatorValue& sv
	    = sepSeparatorValuesPtr[accIndex];
	
	  sepOrigin.remPacker.pack((unsigned**)sepSharedStructure.remDiscreteValuePtrs.ptr,
				   &CliqueBuffer::packedVal[0]);

	  unsigned* remIndexp =
	    sv.iRemHashMap.find(&CliqueBuffer::packedVal[0]);

	  // it must exist
	  assert ( remIndexp != NULL );
	
	  // We've finally got the sep entry.  Add in this clique value's
	  // probability.  Note that bp was initialized during forward
	  // pass.
	  sv.remValues.ptr[*remIndexp].bp() += cliqueValues.ptr[cvn].p;
	} else {
	  // Either separator is all observed, or the separator
	  // is completely contained in the accumulated intersection.
	
	  // keep handy reference for readability.
	  ConditionalSeparatorTable::AISeparatorValue& sv
	    = sepSeparatorValuesPtr[accIndex];

	  // We've finally got the sep entry.  Add in this clique value's
	  // probability.  Note that bp was initialized during forward
	  // pass.
	  sv.remValues.ptr[0].bp() += cliqueValues.ptr[cvn].p;
	}
      }
    }
  }

  // lastly iterate through all separators, and all entries in
  // each separator and do the "divide" (subtraction)
  for (unsigned sepNumber=0;sepNumber<origin.ceReceiveSeparators.size();sepNumber++) {
    // get a handy reference to the current separator
    ConditionalSeparatorTable& sep = 
      separatorTableArray[origin.ceReceiveSeparators[sepNumber]];
    SeparatorClique& sepOrigin = 
      *(sepSharedStructureArray[origin.ceReceiveSeparators[sepNumber]].origin);


    // don't distribute to VE separators or to one that is being skipped.
    if (sep.veSeparator() || sepOrigin.skipMe)
      continue;

    // keep a local variable copy of this around to avoid potential dereferencing.
    ConditionalSeparatorTable::AISeparatorValue * const
      sepSeparatorValuesPtr = sep.separatorValues->ptr; 

    for (unsigned aiNo=0;aiNo < sep.numSeparatorValuesUsed; aiNo ++) {
      ConditionalSeparatorTable::AISeparatorValue* aisep = &(sepSeparatorValuesPtr[aiNo]);
      for (unsigned remNo=0; remNo < aisep->numRemValuesUsed; remNo++) {
	ConditionalSeparatorTable::RemainderValue* sep_entry = &(aisep->remValues.ptr[remNo]);
	// We remove p from bp since bp will already have a factor of
	// p in it. We do this by dividing it out.
	// -
	//
	// We must make sure that if CE stage is entirely zero (i.e.,
	// zero decoding), we do not run DE stage, as in that case it
	// might be the case that sep_entry->p == 0.
	//
	// -
	// In the normal case (CE != 0), we could do direct value
	// reference subtraction in log domain (corresponding to
	// divison in original domain) to ensure that compiler creates
	// no temporaries. In other words, this operation could be either:
	//
	//        sep_entry->bp = sep_entry->bp / sep_entry->p;
	// or 
	//        sep_entry->bp.valref() = sep_entry->bp.valref() - sep_entry->p.valref(); 
	// 
	// We do slower version for now until we are certain this is
	// debugged: We assume here that (!sep_entry->p.zero()) is
	// true since we pruned all zero p's above. If we didn't
	// prune, then 'sep_entry->p == zero' would imply that
	// 'sep_entry->bp == zero', and we would need to do a
	// check. Note that this pruning always occurs, regardless of
	// beam.
	sep_entry->bp() /= sep_entry->p;
      }
    }
  }
}



/*-
 *-----------------------------------------------------------------------
 * PedagogicalCliqueTable::deScatterToOutgoingSeparatorsViterbi()
 *
 *      Viterbi version of deScatterToOutgoingSeparators(). We assume
 *      here that all RV values in the separator are currently set to
 *      appropriate max value the current clique. This is done by
 *      actually assuming that the parent clique of the separators is
 *      currently assigned to the max value.
 *
 *      What we do here is find the index (two indices actually) in
 *      each CE-incomming separator corresponding to the current max
 *      clique assignment, and store this index in the separator
 *      itself. These indices can thus at a later time be used to look
 *      up the separator entries corresponding to the currently (at
 *      the time of this routine call) assigned clique RV values, but
 *      regardless of what the RV values happen to be assigned to at
 *      that later time (and they might be different, say, during the
 *      island algorithm).
 *
 * Preconditions:
 *      Clique data structures and separator must be created. RVs within clique 
 *      are assumed to be set to their maximum (or DE Viterbi max) values.
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      changes all CE-incomming separators.
 *
 * Results:
 *     none
 *
 *-----------------------------------------------------------------------
 */
void 
PedagogicalCliqueTable::
deScatterToOutgoingSeparatorsViterbi(PedagogicalCliqueTable::SharedLocalStructure& sharedStructure,
				     ConditionalSeparatorTable* separatorTableArray,
				     ConditionalSeparatorTable::SharedLocalStructure* sepSharedStructureArray)
{

  // allocate some temporary storage for packed separator values.
  // 128 words is *much* bigger than any possible packed clique value
  // will take on, but it is easy/fast to allocate on the stack right now.
  // unsigned packedVal[128];

  MaxClique& origin = *(sharedStructure.origin);

  if (origin.ceReceiveSeparators.size() == 0)
    return;

  if (origin.hashableNodes.size() == 0) {
    // Do the observed clique case up front right here so we don't
    // need to keep checking below. Here, the clique is observed which
    // means that all connecting separators are also observed. We just
    // sweep through all separators updating the values.
    for (unsigned sepNumber=0;sepNumber<origin.ceReceiveSeparators.size();sepNumber++) {
      ConditionalSeparatorTable& sep = 
	separatorTableArray[origin.ceReceiveSeparators[sepNumber]];
      SeparatorClique& sepOrigin = 
	*(sepSharedStructureArray[origin.ceReceiveSeparators[sepNumber]].origin);

      // don't distribute to VE separators or to one that is being skipped.
      if (sep.veSeparator() || sepOrigin.skipMe)
	continue;

      // in this case, since the cliquueis observed, the separators
      // are also observed and so the zero entries correspond to the
      // max (and the only) entries.
      sep.forwPointer.viterbiAccIndex = 0;
      sep.forwPointer.viterbiRemIndex = 0;
    }

  } else {

    // restore parent clique to the cvn stored in the clique table.

    // unpack clique value 'back_max_cvn' into corresponding random variables and expand
    // any deterministic values.
    const bool imc_nwwoh_p = (origin.packer.packedLen() <= IMC_NWWOH); 
    if (imc_nwwoh_p) {
      origin.packer.unpack((unsigned*)&(cliqueValues.ptr[back_max_cvn].val[0]),
			   (unsigned**)sharedStructure.discreteValuePtrs.ptr);
    } else {
      origin.packer.unpack((unsigned*)cliqueValues.ptr[back_max_cvn].ptr,
			   (unsigned**)sharedStructure.discreteValuePtrs.ptr);
    }
    for (unsigned j=0;j<sharedStructure.fDeterminableNodes.size();j++) {
      RV* rv = sharedStructure.fDeterminableNodes.ptr[j];
      RV2DRV(rv)->assignDeterministicChild();
    }

    // now we iterate through all the separators.
    for (unsigned sepNumber=0;sepNumber<origin.ceReceiveSeparators.size();sepNumber++) {
      // get a handy reference to the current separator
      ConditionalSeparatorTable& sep = 
	separatorTableArray[origin.ceReceiveSeparators[sepNumber]];
      ConditionalSeparatorTable::SharedLocalStructure& sepSharedStructure = 
	sepSharedStructureArray[origin.ceReceiveSeparators[sepNumber]];
      SeparatorClique& sepOrigin = 
	*(sepSharedStructure.origin);

      // don't distribute to VE separators or to one that is being skipped.
      if (sepOrigin.veSeparator || sepOrigin.skipMe)
	continue;

      // keep a local variable copy of this around to avoid potential dereferencing.
      ConditionalSeparatorTable::AISeparatorValue * const
	sepSeparatorValuesPtr = sep.separatorValues->ptr; 
    
      /*
       * There are 3 cases.
       * 1) AI exists and REM exist
       * 2) AI exists and REM doesnt exist
       * 3) AI does not exist, but REM exists
       * AI not exist and REM not exist can't occur.
       */

      unsigned accIndex;
      // TODO: optimize this check away out of loop.
      if (sepOrigin.hAccumulatedIntersection.size() > 0) {
	// an accumulated intersection exists.

	sepOrigin.accPacker.pack((unsigned**)sepSharedStructure.accDiscreteValuePtrs.ptr,
				 &CliqueBuffer::packedVal[0]);
	unsigned* accIndexp =
	  sep.iAccHashMap->find(&CliqueBuffer::packedVal[0]);

	// we should always find something or else something is wrong.
	assert ( accIndexp != NULL ); 
	accIndex = *accIndexp;

      } else {
	// no accumulated intersection exists, everything
	// is in the remainder.
	accIndex = 0;
      }

      if (sepSharedStructure.remDiscreteValuePtrs.size() == 0) {
	// 2) AI exists and REM doesnt exist
	// Then this separator is entirely covered by one or 
	// more other separators earlier in the order.
	sep.forwPointer.viterbiAccIndex = accIndex;
	sep.forwPointer.viterbiRemIndex = 0;
      } else {
	// if we're here, then we must have some remainder
	// pointers.

	// Do the remainder exists in this separator.
	// 
	// either:
	//   1) AI exists and REM exist
	//     or
	//   3) AI does not exist (accIndex == 0), but REM exists
	// 
	
	// keep handy reference for readability.
	ConditionalSeparatorTable::AISeparatorValue& sv
	  = sepSeparatorValuesPtr[accIndex];
	
	sepOrigin.remPacker.pack((unsigned**)sepSharedStructure.remDiscreteValuePtrs.ptr,
				 &CliqueBuffer::packedVal[0]);

	unsigned* remIndexp =
	  sv.iRemHashMap.find(&CliqueBuffer::packedVal[0]);

	if (remIndexp == NULL ) {
	  // print out the rvs.
	  fprintf(stderr,"ERROR: can't find separator rvs values from parent clique in fwrd hash table. Separator rv values follow.\n");
	  printRVSetAndValues(stderr,sepSharedStructure.fNodes,true);
	  fprintf(stderr,"Clique random variables follow:\n");
	  printRVSetAndValues(stderr,sharedStructure.fNodes,true);
	  assert ( remIndexp != NULL );
	}
	
	// We've finally got the sep entry.  Store the sep entry's id.
	sep.forwPointer.viterbiAccIndex = accIndex;
	sep.forwPointer.viterbiRemIndex = *remIndexp;
      }
    }
  }
}



/////////////////////////////////////////////	
/// END OF FILE
/////////////////////////////////////////////
