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

#include "GMTK_FileParser.h"
#include "GMTK_RV.h"
#include "GMTK_DiscRV.h"
#include "GMTK_HidDiscRV.h"
#include "GMTK_ObsDiscRV.h"
#include "GMTK_ContRV.h"
#include "GMTK_GMTemplate.h"
#include "GMTK_GMParms.h"
#include "GMTK_MaxClique.h"
#include "GMTK_MDCPT.h"
#include "GMTK_MSCPT.h"
#include "GMTK_MTCPT.h"
#include "GMTK_Mixture.h"
#include "GMTK_ObservationSource.h"
#include "GMTK_SeparatorClique.h"
#include "GMTK_ConditionalSeparatorTable.h"

VCID(HGID)



unsigned ConditionalSeparatorTable::remHashMapStartingSize;

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        ConditionalSeparatorTable support
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////



/*-
 *-----------------------------------------------------------------------
 * ConditionalSeparatorTable::ConditionalSeparatorTable()
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
ConditionalSeparatorTable::
ConditionalSeparatorTable(SeparatorClique& origin)
  : separatorValues(NULL),iAccHashMap(NULL),preserve(false)
{
  init(origin);
}

void ConditionalSeparatorTable::init(SeparatorClique& origin) 
{
  preserve=false;
  if (origin.veSeparator) {
    // For VE separators, our origin contains the separator tables
    // already pre-generated and constant accross all instances of
    // this separator. Therefore, we special case here and just copy
    // over the pointer values tha are stored in the VE separator that
    // is contained in the origin.
    //
    // We need to make sure not to delete these values however when
    // this object is deleted since normally the separator tables will
    // be reused many times.
    separatorValues = origin.veSepClique->separatorValues;
    iAccHashMap =  origin.veSepClique->iAccHashMap;
    numSeparatorValuesUsed = origin.veSepClique->numSeparatorValuesUsed;
    setToVeSeparatorId();
  } else {

    clearInferenceMemory();

    separatorValues = new cArray< AISeparatorValue >;
    // allocate at one value for now.
    if (origin.hAccumulatedIntersection.size() == 0) {
      // in this case, we'll only need one and never more.
      separatorValues->resize(1);
      // there will always be one used value here.
      numSeparatorValuesUsed = 1;
      if (origin.hRemainder.size() > 0) {
	// So we have no accumulated intersection, and a remainder which
	// means we are in a good position to predict the size of the
	// (necessarily single) remainder vectors from the previous
	// times we used it. Therefore, we do just that, but only in
	// this case.
	separatorValues->ptr[0].remValues.resize(origin.remainderValueSpaceManager.currentSize());
	new (&separatorValues->ptr[0].iRemHashMap)VHashMapUnsignedUnsignedKeyUpdatable
	  (origin.remPacker.packedLen(),origin.remainderValueSpaceManager.currentSize());
	// new (&separatorValues->ptr[0].iRemHashMap)VHashMapUnsignedUnsignedKeyUpdatable(origin.remPacker.packedLen(),REM_HASH_MAP_STARTING_SIZE);
      } else {
	// The separator consists of all observed nodes. 
	// Search in file for key string
	// "ALLOCATE_REMVALUES_ALL_OBSERVED'
	// to find where the nec. single entry is allocated.
	// @@@ try allocating it here
	separatorValues->ptr[0].remValues.resize(1);
	separatorValues->ptr[0].numRemValuesUsed = 0;	
      }
    } else {
      // start with something a bit larger
      // TODO: optimize this.
      const unsigned starting_size = origin.separatorValueSpaceManager.currentSize(); // 3,2000;
      separatorValues->resize(starting_size);
      if (origin.hRemainder.size() > 0) {
	for (unsigned i=0;i<starting_size;i++) {
	  // need to re-construct individual hash tables.
	  new (&separatorValues->ptr[i].iRemHashMap)VHashMapUnsignedUnsignedKeyUpdatable
	    (origin.remPacker.packedLen(),remHashMapStartingSize);
	  // TODO: while we potentially could preallocate default size
	  // of separatorValues->ptr[i].remValues.resize(default); here,
	  // we don't really know what it should be. Since there are
	  // multiple remainders here, and each might be drastically
	  // different in size (some even being zero length), it is not
	  // a good idea to allocate anything at all and we let it
	  // lazily be sized as needed. The TODO is to come up with a
	  // better scheme (e.g., keep a counter array to keep track of
	  // the number at each allocation level, and then allocate this
	  // to the min size of the previous time that that a non-zero
	  // number of cases.
	}
      } else {
	// things such as array separatorValues.ptr[i].remValues will be sized as needed later.
	// Search for tag 'ALLOCATE_REMVALUES_OPTION' in this file for where it is allocated.
      }
      // need to re-construct the hash table.
      iAccHashMap = new VHashMapUnsignedUnsignedKeyUpdatable
	(origin.accPacker.packedLen(),starting_size); // 2
      numSeparatorValuesUsed = 0;
    }
  }
}

ConditionalSeparatorTable::SharedLocalStructure::
SharedLocalStructure(SeparatorClique& _origin,
		     vector <RV*>& newRvs,
		     map < RVInfo::rvParent, unsigned >& ppf,
		     const unsigned int frameDelta)
{

  origin = &_origin;

  set<RV*>::iterator it;

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
  }

  i=0;
  fAccumulatedIntersection.resize(origin->accumulatedIntersection.size());
  for (it = origin->accumulatedIntersection.begin();
       it != origin->accumulatedIntersection.end();
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
    // BUG!! (i needs to be incremented, check other cases too, but then this variable is never used).
    fAccumulatedIntersection[i] = nrv;
  }

  i=0;
  fRemainder.resize(origin->remainder.size());
  for (it = origin->remainder.begin();
       it != origin->remainder.end();
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
    fRemainder[i++] = nrv;
  }

  // Separator accumulated intersection values only store/hash values
  // of hidden (thus necessarily discrete) variables since they are
  // the only thing that change.
  accDiscreteValuePtrs.resize(origin->hAccumulatedIntersection.size());
  for (i=0;i<accDiscreteValuePtrs.size();i++) {
    // get the hidden rv for this location
    RV* rv = origin->hAccumulatedIntersection[i];;
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    

    if ( ppf.find(rvp) == ppf.end() ) {
      coredump("ERROR: can't find assigned rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	       rv->name().c_str(),rv->frame(),frameDelta,
	       rvp.first.c_str(),rvp.second);
    }
    RV* nrv = newRvs[ppf[rvp]];

    // hidden nodes are always discrete (in this version).
    // TODO: add hidden continuous variable
    DiscRV* drv = 
      (DiscRV*)nrv;

    // grab a pointer directly to its value for easy access later.
    accDiscreteValuePtrs[i] = &(drv->val);
  }

  // Separator remainder values only store/hash values of hidden (thus
  // necessarily discrete) variables since they are the only thing
  // that change.
  remDiscreteValuePtrs.resize(origin->hRemainder.size());
  for (i=0;i<remDiscreteValuePtrs.size();i++) {
    // get the hidden rv for this location
    RV* rv = origin->hRemainder[i];
    RVInfo::rvParent rvp;
    rvp.first = rv->name();
    rvp.second = rv->frame()+frameDelta;    

    if ( ppf.find(rvp) == ppf.end() ) {
      coredump("ERROR: can't find assigned rv %s(%d+%d)=%s(%d) in unrolled RV set\n",
	       rv->name().c_str(),rv->frame(),frameDelta,
	       rvp.first.c_str(),rvp.second);
    }
    RV* nrv = newRvs[ppf[rvp]];

    // hidden nodes are always discrete (in this version).
    DiscRV* drv = 
      (DiscRV*)nrv;

    // grab a pointer directly to its value for easy access later.
    remDiscreteValuePtrs[i] = &(drv->val);
  }
}



// Next, version of the above constructor to create specifically for a
// VE separator that only lives and stays in its origin separator, and
// used in SeparatorClique::prepareForUnrolling(), where the set of
// RVs are the same as the origin separator, and there is no frame
// delta. This separator is then used to create copies of the
// pointers in the table data structures when we need
// a ve separator that is iterated against a set of random variables.
// 
// A inference vec separator is used just for its values at every
// possible time frame.  See comments above.
//
// this routine assumes that the sepSharedStructure argument is empty
// and that member variables can be removed/resized.

ConditionalSeparatorTable
::ConditionalSeparatorTable(SeparatorClique& origin,
			    ConditionalSeparatorTable::SharedLocalStructure& sepSharedStructure)
  : separatorValues(NULL),iAccHashMap(NULL)
{

  set<RV*>::iterator it;

  sepSharedStructure.origin = &origin;

  // this table is a ve separator
  assert ( origin.veSeparator );
  setToVeSeparatorId();

  // clone over nodes RVs.
  sepSharedStructure.fNodes.resize(origin.nodes.size());
  unsigned i=0;
  for (it = origin.nodes.begin();
       it != origin.nodes.end();
       it++) {
    RV* rv = (*it);
    sepSharedStructure.fNodes[i++] = rv;
  }

  i=0;
  sepSharedStructure.fAccumulatedIntersection.resize(origin.accumulatedIntersection.size());
  for (it = origin.accumulatedIntersection.begin();
       it != origin.accumulatedIntersection.end();
       it++) {
    RV* rv = (*it);
    sepSharedStructure.fAccumulatedIntersection[i] = rv;
  }

  i=0;
  sepSharedStructure.fRemainder.resize(origin.remainder.size());
  for (it = origin.remainder.begin();
       it != origin.remainder.end();
       it++) {
    RV* rv = (*it);
    sepSharedStructure.fRemainder[i++] = rv;
  }

  // Separator accumulated intersection values only store/hash values
  // of hidden (thus necessarily discrete) variables since they are
  // the only thing that change.
  sepSharedStructure.accDiscreteValuePtrs.resize(origin.hAccumulatedIntersection.size());
  for (i=0;i<sepSharedStructure.accDiscreteValuePtrs.size();i++) {
    // get the hidden rv for this location
    RV* rv = origin.hAccumulatedIntersection[i];;
    DiscRV* drv = (DiscRV*)rv;
    // grab a pointer directly to its value for easy access later.
    sepSharedStructure.accDiscreteValuePtrs[i] = &(drv->val);
  }

  // Separator remainder values only store/hash values of hidden (thus
  // necessarily discrete) variables since they are the only thing
  // that change.
  sepSharedStructure.remDiscreteValuePtrs.resize(origin.hRemainder.size());
  for (i=0;i<sepSharedStructure.remDiscreteValuePtrs.size();i++) {
    // get the hidden rv for this location
    RV* rv = origin.hRemainder[i];
    // hidden nodes are always discrete (in this version).
    DiscRV* drv = (DiscRV*)rv;
    // grab a pointer directly to its value for easy access later.
    sepSharedStructure.remDiscreteValuePtrs[i] = &(drv->val);
  }

  separatorValues = new cArray< AISeparatorValue >;
  // allocate at one value for now.
  if (origin.hAccumulatedIntersection.size() == 0) {
    // in this case, we'll only need one and never more.
    separatorValues->resize(1);
    // there will always be one used value here.
    numSeparatorValuesUsed = 1;
    if (origin.hRemainder.size() > 0) {
      // So we have no accumulated intersection, and a remainder which
      // means we are in a good position to predict the size of the
      // (necessarily single) remainder vectors from the previous
      // times we used it. Therefore, we do just that, but only in
      // this case.
      separatorValues->ptr[0].remValues.resize(origin.remainderValueSpaceManager.currentSize());
      new (&separatorValues->ptr[0].iRemHashMap)VHashMapUnsignedUnsignedKeyUpdatable
	(origin.remPacker.packedLen(),origin.remainderValueSpaceManager.currentSize());
      // new (&separatorValues->ptr[0].iRemHashMap)VHashMapUnsignedUnsignedKeyUpdatable(origin.remPacker.packedLen(),REM_HASH_MAP_STARTING_SIZE);
    } else {
      // The separator consists of all observed nodes. 
      // Search in file for key string
      // "ALLOCATE_REMVALUES_ALL_OBSERVED'
      // to find where the nec. single entry is allocated.
      // @@@ try allocating it here
      separatorValues->ptr[0].remValues.resize(1);
      separatorValues->ptr[0].numRemValuesUsed = 0;

    }
  } else {
    // start with something a bit larger
    // TODO: optimize this.
    const unsigned starting_size = origin.separatorValueSpaceManager.currentSize(); // 3,2000;
    separatorValues->resize(starting_size);
    if (origin.hRemainder.size() > 0) {
      for (unsigned i=0;i<starting_size;i++) {
	// need to re-construct individual hash tables.
	new (&separatorValues->ptr[i].iRemHashMap)VHashMapUnsignedUnsignedKeyUpdatable
	  (origin.remPacker.packedLen(),remHashMapStartingSize);
	// TODO: while we potentially could preallocate default size
	// of separatorValues->ptr[i].remValues.resize(default); here,
	// we don't really know what it should be. Since there are
	// multiple remainders here, and each might be drastically
	// different in size (some even being zero length), it is not
	// a good idea to allocate anything at all and we let it
	// lazily be sized as needed. The TODO is to come up with a
	// better scheme (e.g., keep a counter array to keep track of
	// the number at each allocation level, and then allocate this
	// to the min size of the previous time that that a non-zero
	// number of cases.
      }
    } else {
      // things such as array separatorValues.ptr[i].remValues will be sized as needed later.
      // Search for tag 'ALLOCATE_REMVALUES_OPTION' in this file for where it is allocated.
    }
    // need to re-construct the hash table.
    iAccHashMap = new VHashMapUnsignedUnsignedKeyUpdatable
      (origin.accPacker.packedLen(),starting_size); // 2
    numSeparatorValuesUsed = 0;
  }

}


/*-
 *-----------------------------------------------------------------------
 * ConditionalSeparatorTable::insert()
 *
 *    Insert whatever the current RV values are set to into the current
 *    inference separator (based in their RV values).
 *
 * Preconditions:
 *   1) separator tables must be created, meaning that the same
 *      conditions must holed as just before ceSendToOutgoingSeparator() is
 *      called.
 *
 * Postconditions:
 *    Separator table has entry added, according to current accumulated intersection
 *    and remainder, and current RV values.
 *
 * Side Effects:
 *    potential memory allocations and hash table adjustments.
 *
 * Results:
 *     nothing
 *
 *-----------------------------------------------------------------------
 */
void ConditionalSeparatorTable::insert(SeparatorClique& origin,
				       ConditionalSeparatorTable::SharedLocalStructure& sepSharedStructure)
{

  // keep a local variable copy of this around to avoid potential
  // dereferencing.  This one cannot be const since it might change
  // during a resize, in which case we need to reassign this variable.
  ConditionalSeparatorTable::AISeparatorValue * 
    sepSeparatorValuesPtr = separatorValues->ptr; 

  // there must be someplace to insert in this case.
  assert (!(origin.hAccumulatedIntersection.size() == 0 && origin.hRemainder.size() == 0));

  // precompute some constants.
  const bool isc_nwwoh_ai_p = (origin.accPacker.packedLen() <= ISC_NWWOH_AI);
  const bool isc_nwwoh_rm_p = (origin.remPacker.packedLen() <= ISC_NWWOH_RM);
  const bool sep_origin_hAccumulatedIntersection_exists_p =
    (origin.hAccumulatedIntersection.size() > 0);
  const bool sep_remDiscreteValuePtrs_exists_p = 
    (sepSharedStructure.remDiscreteValuePtrs.size() > 0);
  
  unsigned accIndex;
  if (sep_origin_hAccumulatedIntersection_exists_p) { 
    // an accumulated intersection exists.

    // make sure there is at least one available accumulated intersection entry
    assert ( numSeparatorValuesUsed <= separatorValues->size());
    if (numSeparatorValuesUsed >= separatorValues->size()) {
      
      infoMsg(Inference, Max+5,"ac-rsz,");

      const unsigned old_size = separatorValues->size();
      // TODO: optimize this size re-allocation.
      if (numSeparatorValuesUsed >= origin.separatorValueSpaceManager.currentSize()) 
	origin.separatorValueSpaceManager.advanceToNextSize();
      separatorValues->resizeAndCopy(origin.separatorValueSpaceManager.currentSize()); 
      sepSeparatorValuesPtr = separatorValues->ptr;
      if (isc_nwwoh_ai_p) {
	// Then the above resize just invalided all our pointers to
	// keys (which in this case are compressed RV values for the
	// acc inter), but it did not invalidate the hash items (which
	// in this case are the array indices in the accumulated
	// intersection corresponding to a given compressed acc intr
	// RV values). We thus go through and correct the key pointers
	// within the hash table.  TODO: think of a better way to do
	// this that also looses no efficiency.
	for (unsigned i=0;i<iAccHashMap->tableSize();i++) {
	  if (!iAccHashMap->tableEmpty(i)) {
	    iAccHashMap->tableKey(i)
	      = &(sepSeparatorValuesPtr[iAccHashMap->tableItem(i)].val[0]);
	  }
	}
      }
      const unsigned new_size = separatorValues->size();
      // if (remDiscreteValuePtrs.size() > 0) {
      if (sep_remDiscreteValuePtrs_exists_p) {
	for (unsigned i=old_size;i<new_size;i++) {
	  // re-construct hash tables only for new entries.
	  new (&sepSeparatorValuesPtr[i].iRemHashMap)
	    VHashMapUnsignedUnsignedKeyUpdatable
	    (origin.remPacker.packedLen(),remHashMapStartingSize);
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
      accKey = &(sepSeparatorValuesPtr[numSeparatorValuesUsed].val[0]);
      origin.accPacker.pack((unsigned**)sepSharedStructure.accDiscreteValuePtrs.ptr,
			    accKey);
    } else {
      accKey = origin.accValueHolder.curCliqueValuePtr();
      origin.accPacker.pack((unsigned**)sepSharedStructure.accDiscreteValuePtrs.ptr,
			    accKey);
      // check if this value combination already lives in
      // origin's value holder hash table and if so, use that.
      bool foundp;
      accKey = origin.accSepValHashSet.insert(accKey,foundp);
      if (!foundp) {
	// only allocate a new value if it was inserted.
	origin.accValueHolder.allocateCurCliqueValue();
      }
      // store the pointer in case we use it.
      sepSeparatorValuesPtr[numSeparatorValuesUsed].ptr = accKey;
    }

      
    bool foundp;
    unsigned* accIndexp =
      iAccHashMap->insert(accKey,
			  numSeparatorValuesUsed,
			  foundp);

    if (!foundp) {
      //  add the values we just used. 
      numSeparatorValuesUsed++;
    }
    accIndex = *accIndexp;

    infoMsg(IM::Inference, Max+5,"inst:ai=%d,ky=%X,",accIndex,*accKey);

    // TODO: optimize this check out of loop.
    // if (remDiscreteValuePtrs.size() == 0) {
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
	sv.remValues.ptr[0].p = 1.0; // probability being assigned is unity
	// if (JunctionTree::viterbiScore)
	// sv.remValues.ptr[0].backPointer = cvn;
      } else {
	// already there so must have hit before.
	// for now, die with an assertion if this case occurs (can't insert twice)
	assert ( 0 );
      }
      return;
    }

  } else {
    accIndex = 0;
    infoMsg(IM::Inference, Max+5,"inst:ai=%d,",accIndex);
  }




  // If we're here, then we are guaranteed must have some remainder
  // pointers, i.e., we could do:
  //    assert (remDiscreteValuePtrs.size() > 0);
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

    infoMsg(IM::Inference, Max+5,"rm-rsz,u=%d,f=%d,",sv.numRemValuesUsed,sv.remValues.size());

    // TODO: optimize this growth rate.
    // start small but grow fast.
    // sv.remValues.resizeAndCopy(1+sv.remValues.size()*2); // *3
    sv.remValues.resizeAndCopy(origin.remainderValueSpaceManager.nextSizeFrom(sv.remValues.size()));
    origin.remainderValueSpaceManager.setCurrentAllocationSizeIfLarger(sv.remValues.size());

    infoMsg(IM::Inference, Max+5,"=%d,",sv.remValues.size());

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
    origin.remPacker.pack((unsigned**)sepSharedStructure.remDiscreteValuePtrs.ptr,
			  remKey);
  } else {
    // grab pointer to next packed clique value to be used.
    remKey = origin.remValueHolder.curCliqueValuePtr();
    origin.remPacker.pack((unsigned**)sepSharedStructure.remDiscreteValuePtrs.ptr,
			  remKey);
    // check if this value combination already lives in
    // origin's value holder hash table and if so, use that.
    bool foundp;
    remKey = origin.remSepValHashSet.insert(remKey,foundp);
    if (!foundp) {
      // only allocate a new value if it was inserted.
      origin.remValueHolder.allocateCurCliqueValue();
    }
    // store the pointer in case we use it.
    sv.remValues.ptr[sv.numRemValuesUsed].ptr = remKey;
  }

  bool foundp;
  unsigned* remIndexp =
    sv.iRemHashMap.insert(remKey,
			  sv.numRemValuesUsed,
			  foundp);

  infoMsg(IM::Inference, Max+5,"inst:ky=%X,rii=%d\n",*remKey,*remIndexp);

  if (!foundp) {
    // add the values we just used. 
    sv.numRemValuesUsed++;
  } else {
    // already found
    infoMsg(IM::Inference, Max+5,"already found\n");
    fflush(stdout);
    fflush(stderr);
    assert ( 0 );
  }

  // We've finally got the entry, assign the probability as unity.
  sv.remValues.ptr[*remIndexp].p = 1.0;


}


set <RV*> 
ConditionalSeparatorTable::SharedLocalStructure::returnRVsAsSet()
{
  set<RV*> rc;
  for (unsigned i=0;i<fNodes.size();i++) {
    rc.insert(fNodes[i]);
  }
  return rc;
}



/*-
 *-----------------------------------------------------------------------
 * ConditionalSeparatorTable::ceSeparatorPrune()
 *
 *    Collect Evidence, Separator Prune: This routine will prune away
 *    part of a previously instantiated separator based on the current
 *    separator beam width.
 *
 * Preconditions:
 *   1) separator table must be created, meaning that either:
 *
 *        MaxCliqueTable::ceSendToOutgoingSeparator()
 *      must have been called sending a message (projection downto) this separator.
 *      Separator must not be a VE separator.
 *
 * Postconditions:
 *    Separator table has been pruned, and memory for it has been re-allocated to
 *    fit the smaller size. All hash tables adjusted.
 *
 * Side Effects:
 *    potential memory allocations and hash table adjustments.
 *
 * Results:
 *     nothing
 *
 *-----------------------------------------------------------------------
 */

void 
ConditionalSeparatorTable::ceSeparatorPrune(SeparatorClique& origin)
{

  // keep a local variable copy of this around to avoid potential dereferencing.
  AISeparatorValue * const
    separatorValuesPtr = separatorValues->ptr; 

  if (origin.separatorBeam != (-LZERO)) {
    // only do this if separator beam pruning is not turned off.

    // we shouldn't have to check this since we should never be pruning
    // a VE separator.
    // if (veSeparator())
    // return;

    // compute max and current state space.
    logpr maxCEsepValue;
    unsigned originalTotalStateSpace = 0;
    for (unsigned asv=0;asv<numSeparatorValuesUsed;asv++) {
      originalTotalStateSpace += separatorValuesPtr[asv].numRemValuesUsed;
      for (unsigned rsv=0;rsv<separatorValuesPtr[asv].numRemValuesUsed;rsv++) {
	if (separatorValuesPtr[asv].remValues.ptr[rsv].p > maxCEsepValue)
	  maxCEsepValue = separatorValuesPtr[asv].remValues.ptr[rsv].p;
      }
    }

    // Check for all observed case, or a case where there is only one
    // entry, in which case we never do anything.
    if (originalTotalStateSpace == 1)
      return;

    // check if we have a zero separator, and if we do, print message.
    // Do this before separator pruning.
    if (originalTotalStateSpace == 0) {
      infoMsg(IM::Inference, IM::Mod,"WARNING: ZERO SEPARATOR: separator with no entries. Final probability will be zero.\n");
      // nothing to prune, so we return.
      return;
    }

    // create an ininitialized variable
    logpr beamThreshold((void*)0);
    // break into the logp to avoid unnecessary zero checking.
    beamThreshold.valref() = maxCEsepValue.valref() - origin.separatorBeam;

    // pointers to the ht keys for the two entries.
    unsigned** ht_prune_key_p=NULL;
    unsigned** ht_swap_key_p=NULL;

    // go through and shrink guys less than maximum.
    unsigned newTotalStateSpace = 0;  
    for (unsigned asv=0;asv<numSeparatorValuesUsed;asv++) {
      const unsigned origNumRemValuesUsed = separatorValuesPtr[asv].numRemValuesUsed;
      for (unsigned rsv=0;rsv<separatorValuesPtr[asv].numRemValuesUsed;) {
	if (separatorValuesPtr[asv].remValues.ptr[rsv].p < beamThreshold) {

	  if (separatorValuesPtr[asv].numRemValuesUsed > 1) {


	    // We prune away entry for rsv, by swapping it in last
	    // position. Here, however, it is not as easy as with a clique
	    // separator as we have also to deal with the hash
	    // table. Specifically, we need to swap index entries in hash
	    // table as well. Note that we can not remove the hash table
	    // entry for the one that got pruned away without re-hashing
	    // the entire hash table. The reason is that if the entry that
	    // got removed was a collision for another entry that is in
	    // the table, then removing the collision will make the other
	    // entry inaccessible. Therefore, for now, the hash table
	    // does not shrink while the table does.
	    // TODO: test if it is better to just prune here and just
	    // rehash everything.

	    // the index of the entry being swapped with the
	    // one that is being pruned.
	    const unsigned swap_index = separatorValuesPtr[asv].numRemValuesUsed-1;

	    // First, get pointers to hash table index values for the two
	    // entries corresponding to the one we are prunning
	    // and the one ware swapping it with.
	    unsigned* prune_index_p;
	    unsigned* swap_index_p;

	    // the keys for the two entries.
	    unsigned* prune_key_p;
	    unsigned* swap_key_p;


	    if (origin.remPacker.packedLen() <= ISC_NWWOH_RM) {
	      prune_key_p = &(separatorValuesPtr[asv].remValues.ptr[rsv].val[0]);
	      swap_key_p = &(separatorValuesPtr[asv].remValues.ptr[swap_index].val[0]);
	    } else {
	      prune_key_p = separatorValuesPtr[asv].remValues.ptr[rsv].ptr;
	      swap_key_p = separatorValuesPtr[asv].remValues.ptr[swap_index].ptr;
	    }
	
	    prune_index_p =  separatorValuesPtr[asv].iRemHashMap.find(prune_key_p,ht_prune_key_p);
	    // it must exist
	    assert ( prune_index_p != NULL );
	    swap_index_p =  separatorValuesPtr[asv].iRemHashMap.find(swap_key_p,ht_swap_key_p);
	    // it must exist
	    assert ( swap_index_p != NULL );

	    // swap the entries in the separator remainder (rv_val,
	    // prob) table. We can't do this any earlier than here
	    // since the hash finding above uses pointers to these
	    // entries.
	    swap(separatorValuesPtr[asv].remValues.ptr[rsv],
		 separatorValuesPtr[asv].remValues.ptr[swap_index]);

	    // and swap the hash table item index values (i.e., the ht items are integer indices
	    // into the separator remainder (rv_val, prob) table. We want the hash table entry for the
	    // item that we are not pruning to now point to the separator remainder entry that is not
	    // being pruned away (rather than its old entry, which is no longer being used). This
	    // does not change the position in the hash table of this entry, rather it only changes
	    // what the hash table entry is pointing back to in the separator remainder table.
	    swap((*prune_index_p),(*swap_index_p));

	    // and swap the hash table keys if they are pointers to
	    // the arrays which just got swapped. In other words, when
	    // the separator remainder value is small enough to fit in
	    // one machine word, then the hash table key values will
	    // point directly into the slot in the separator remainder
	    // (rv_val, prob) table rather than pointing to some
	    // globally shared value pool. In this case, since the
	    // separator remainder table (containing rv values
	    // directly rather than pointers) has changed, we need to
	    // adjust the hash table so that its pointer to key is
	    // appropriate.
	    if (origin.remPacker.packedLen() <= ISC_NWWOH_RM) {
	      // printf("foobarbaz");
	      swap((*ht_prune_key_p),(*ht_swap_key_p));
	    }

	    // TODO: we don't need to swap above, rather we just need
	    // to update the slot that is not getting pruned away.

	    // decrease values
	    separatorValuesPtr[asv].numRemValuesUsed--;

	    // TODO: the above needs to be looked at for 64bit (64-bit) case.

	  } else {
	    // then separatorValuesPtr[asv].numRemValuesUsed == 1. This
	    // will happen under two conditiosn.
	    //  1) There is no remainder, meaning (sepOrigin.hRemainder.size() == 0), 
	    //     and the entire separator is
	    //     in the accumulated intersection with other separators.
	    //     In this case there is no hash table at all.
	    //  2) All other entries in the remainder for this particular
	    //     accumulated intersection slot have been pruned away. In this
	    //     case, the hash table does exist, but it'll still be deleted
	    //     upon calling the destructor.

	    // so, what we do is just set the the number of values to zero.
	    // Other code will need to thus check for empty separator acc. intersection
	    // values.
	    separatorValuesPtr[asv].numRemValuesUsed = 0;

	  }

	} else {
	  rsv++;
	}
      }
      newTotalStateSpace += separatorValuesPtr[asv].numRemValuesUsed;
      if (separatorValuesPtr[asv].numRemValuesUsed < origNumRemValuesUsed) {
	if (separatorValuesPtr[asv].numRemValuesUsed == 0 ) {
	  // should/could remove accumulator entry here as well. 
	}
	// - re-allocate memory & adjust hash table.
	// - separatorValuesPtr[asv].remValues
	// - possibly re-hash hash tables if necessary.
      }
    }

    infoMsg(IM::Inference, IM::Med,"Separator beam pruning, Max cv = %f, thres = %f. Original sep state space = %d, new sep state space = %d\n",
	    maxCEsepValue.valref(),
	    beamThreshold.valref(),
	    originalTotalStateSpace,newTotalStateSpace);

  }

#if 0

  // reallocate memory so that nothing is wasted.
  for (unsigned asv=0;asv<numSeparatorValuesUsed;asv++) {
    // shrink down if we can gain more than about 1.6% of size.
    if ((separatorValuesPtr[asv].remValues.size() - separatorValuesPtr[asv].numRemValuesUsed) > separatorValuesPtr[asv].remValues.size()/64) {
      separatorValuesPtr[asv].remValues.resizeAndCopy(separatorValuesPtr[asv].numRemValuesUsed);
    }
  }

  // TODO: get the rest of this code working, will need in some cases
  // to re-hash the above hash tables. Not sure this is worth it.

#endif 

}



/*-
 *-----------------------------------------------------------------------
 * ConditionalSeparatorTable::reportMemoryUsageTo()
 *
 *    Report current memory usage of this inference sep clique *and*
 *    the origin sep clique to the file in units of MBs.
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
ConditionalSeparatorTable::
reportMemoryUsageTo(SeparatorClique& origin, FILE *f)
{
  // Memory: IC=Inference Separator  Clique, AI=accumulated intersection
  fprintf(f,"*MEM:ISC AI(used=%lu,all=%lu=%luMB),AIH(%luMB)",
	  (unsigned long)numSeparatorValuesUsed,
	  (unsigned long)separatorValues->size(),
	  (unsigned long)((1+(sizeof(AISeparatorValue)*(unsigned long)separatorValues->size())/(1024ul*1024ul))),
	  (iAccHashMap != NULL)? (unsigned long)(1 + iAccHashMap->bytesRequested()/(1024ul*1024ul)) : 0 );

  // sum up stats from the remainders
  unsigned long remUsed = 0;
  unsigned long remAllocated = 0;
  unsigned long remHashAllocated = 0;

  for (unsigned long i =0; i< numSeparatorValuesUsed; i++) {
    remUsed += separatorValues->ptr[i].numRemValuesUsed;
    remAllocated += separatorValues->ptr[i].remValues.size();
    remHashAllocated += separatorValues->ptr[i].iRemHashMap.bytesRequested();
  }
  // note: remUsed is also equal to the state space size of the separator.

  // REM = remainder
  fprintf(f,"REM(used=%lu,all=%lu=%luMB),RH(%luMB),",
	  (unsigned long)remUsed,
	  (unsigned long)remAllocated,
	  (unsigned long)((1+(sizeof(RemainderValue)*(unsigned long)remAllocated)/(1024ul*1024ul))),
	  (unsigned long)((1+remHashAllocated/(1024ul*1024ul))));
  origin.reportMemoryUsageTo(f);
  fprintf(f,"\n");

}





/////////////////////////////////////////////	
/// END OF FILE
/////////////////////////////////////////////
