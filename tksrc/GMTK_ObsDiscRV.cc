/*
 * GMTK_ObsDiscRV.cc
 *
 * Observed discrete random variable.
 * 
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *  $Header$
 *
 * Copyright (c) 2001, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about the suitability of this software
 * for any purpose. It is provided "as is" without express or implied warranty.
 *
 *
 * The top level GMTK random variable object for the RV class hierarchy.
 *
 *
 *
 */



#include "general.h"
#if HAVE_CONFIG_H
#include <config.h>
#endif
#if HAVE_HG_H
#include "hgstamp.h"
#endif
VCID(HGID)


#include <iostream>
#include <fstream>
#include <typeinfo>

#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <string.h>

#include "GMTK_ObsDiscRV.h"
#include "GMTK_MTCPT.h"


/*-
 *-----------------------------------------------------------------------
 * printSelf()
 *      prints a one-line summary of the detailed information about this RV.
 *
 * Preconditions:
 *      none
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      self is printed.
 *
 *-----------------------------------------------------------------------
 */
void ObsDiscRV::printSelf(FILE *f,bool nl)
{
  printNameFrameValue(f,false);
  fprintf(f," observed discrete cardinality = %d%s",cardinality,nls(nl));
}



/*-
 *-----------------------------------------------------------------------
 * printSelfVerbose()
 *      prints a multi-line verbose description of this RV.
 *
 * Preconditions:
 *      none
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      self is printed.
 *
 *-----------------------------------------------------------------------
 */
void ObsDiscRV::printSelfVerbose(FILE *f)
{
  fprintf(f,"Observed Discrete Random variable:\n");
  printNameFrameValue(f,true);
  fprintf(f," From line %d in file %s\n",rv_info.fileLineNumber,rv_info.rvFileName.c_str());
  fprintf(f,"RV has cardinality = %d\n",cardinality);
}



/*-
 *-----------------------------------------------------------------------
 * cloneRVShell()
 *      clones a shell of the current random variable (see GMTK_RV.h for docs)
 *
 * Preconditions:
 *      RV must be instantiated and with parameters (i.e., what lives in the template RVs).
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      self is printed.
 *
 *-----------------------------------------------------------------------
 */
ObsDiscRV* ObsDiscRV::cloneRVShell()
{
  ObsDiscRV*rv = (ObsDiscRV*)DiscRV::cloneRVShell();
  // make sure to also set value since it might be an inline 'value'
  // observation which needs to be retained.
  rv->val = val;
  return rv;
}



/*-
 *-----------------------------------------------------------------------
 * setToObservedValue()
 *
 *   set the RV (which must be observed) to the observed value in a
 *   file (this should be done once outside of the inference inner
 *   loops).
 *
 * Preconditions:
 *      RV must be instantiated
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      val is changed.
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
void ObsDiscRV::setToObservedValue() 
{
    // observed, so set value from observation matrix

    if (rv_info.rvFeatureRange.filled == RVInfo::FeatureRange::fr_FirstIsValue) {
      val = rv_info.rvFeatureRange.firstFeatureElement;
    } else if (rv_info.rvFeatureRange.filled == RVInfo::FeatureRange::fr_FrameNumIsValue) {
      assert (globalObservationMatrix->active());
      if (frame() >= cardinality) 
	error("ERROR: RV '%s(%d)' has cardinality %d, but current frame value %d too large to store in RV with this cardinality (in segment %d of frame length %d).\n",
	      name().c_str(),
	      frame(),
	      cardinality,
	      frame(),
	      globalObservationMatrix->segmentNumber(),
	      globalObservationMatrix->numFrames());
      val = frame();
    } else if (rv_info.rvFeatureRange.filled == RVInfo::FeatureRange::fr_EmarfNumIsValue) {
      assert (globalObservationMatrix->active());
      FileSource fs;
      if ( typeid(*globalObservationMatrix) != typeid(fs) ) {
	error("ERROR: RV '%s(%d)' online inference cannot support emarf as value",name().c_str(),frame());
      }
      unsigned emarf = globalObservationMatrix->numFrames() - frame();
      if (emarf >= cardinality) 
	error("ERROR: RV '%s(%d)' has cardinality %d, but current emarf value %d too large to store in RV with this cardinality (in segment %d of frame length %d).\n",
	      name().c_str(),
	      frame(),
	      cardinality,
	      emarf,
	      globalObservationMatrix->segmentNumber(),
	      globalObservationMatrix->numFrames());
      val = emarf;
    } else if (rv_info.rvFeatureRange.filled == RVInfo::FeatureRange::fr_NumFramesIsValue) {
      assert (globalObservationMatrix->active());
      FileSource fs;
      if ( typeid(*globalObservationMatrix) != typeid(fs) ) {
	error("ERROR: RV '%s(%d)' online inference cannot support number of frames as value",name().c_str(),frame());
      }
      if (globalObservationMatrix->numFrames() >= cardinality) 
	error("ERROR: RV '%s(%d)' has cardinality %d, but current num frames %d too large to store in RV with this cardinality (in segment %d of frame length %d).\n",
	      name().c_str(),
	      frame(),
	      cardinality,
	      globalObservationMatrix->numFrames(),
	      globalObservationMatrix->segmentNumber(),
	      globalObservationMatrix->numFrames());
      val = globalObservationMatrix->numFrames();
    } else if (rv_info.rvFeatureRange.filled == RVInfo::FeatureRange::fr_SegmentNumIsValue) {
      assert (globalObservationMatrix->active());
      if (globalObservationMatrix->segmentNumber() >= cardinality) 
	error("ERROR: RV '%s(%d)' has cardinality %d, but current segment number %d too large to store in RV with this cardinality (in segment %d of frame length %d).\n",
	      name().c_str(),
	      frame(),
	      cardinality,
	      globalObservationMatrix->segmentNumber(),
	      globalObservationMatrix->segmentNumber(),
	      globalObservationMatrix->numFrames());
      val = globalObservationMatrix->segmentNumber();
    } else if (rv_info.rvFeatureRange.filled == RVInfo::FeatureRange::fr_NumSegmentsIsValue) {
      assert (globalObservationMatrix->active());
      FileSource *fs = static_cast<FileSource *>(globalObservationMatrix);
      FileSource dummyFS;
      if ( typeid(*globalObservationMatrix) != typeid(dummyFS) ) {
	error("ERROR: RV '%s(%d)' online inference cannot support number of segments as value",name().c_str(),frame());
      }
      if (fs->numSegments() >= cardinality) 
	error("ERROR: RV '%s(%d)' has cardinality %d, but current number segments %d too large to store in RV with this cardinality (in segment %d of frame length %d).\n",
	      name().c_str(),
	      frame(),
	      cardinality,
	      fs->numSegments(),
	      globalObservationMatrix->segmentNumber(),
	      fs->numFrames());
      val = fs->numSegments();
    } else if (rv_info.rvFeatureRange.filled == RVInfo::FeatureRange::fr_Range) {
      assert (globalObservationMatrix->active());
      // printf("getting value of random variable '%s', time index %d, el %d\n",
      // label.c_str(),timeIndex,featureElement);
      unsigned tmp = globalObservationMatrix->unsignedAtFrame(frame(),featureElement());
      if (tmp >= (unsigned)cardinality) 
	error("ERROR: RV '%s' (file:line '%s:%d') at time index %d has cardinality %d, "
	      "but feature element position %d in observation file (time %d of segment %d) "
	      "has value %u = 0x%X.\n",
	      name().c_str(),
	      rv_info.rvFileName.c_str(),
	      rv_info.fileLineNumber,
	      frame(),
	      cardinality,
	      featureElement(),
	      frame(),
	      globalObservationMatrix->segmentNumber(),
	      tmp, tmp);
      val = tmp;
    } else {
      // shouldn't happen.
      assert (0);
    }
    // otherwise, we keep the value set to what it was before.
    return;
}



void ObsDiscRV::computeParentsSatisfyingChild(
	    // input arguments
	    unsigned par, // parent number
	    vector <RV*> & parents, 
	    vector <RV*> & hiddenParents,
	    PackCliqueValue& hiddenParentPacker,
	    sArray < DiscRVType*>& hiddenNodeValPtrs,
	    RV* child,
	    // output arguments
	    sArray < unsigned >& packedParentVals,
	    unsigned& num)
{
  assert ( !switching() && deterministic() && curCPT->cptType == CPT::di_MTCPT );
  MTCPT* mtcpt = (MTCPT*) curCPT;
  return mtcpt->computeParentsSatisfyingChild(par,parents,hiddenParents,hiddenParentPacker,
					      hiddenNodeValPtrs,child,packedParentVals,num);
}





/*-
 *-----------------------------------------------------------------------
 * setObservedRVs()
 *   sets the observed RVs to their values, either taking values from
 *   the global observation matrix, or taking the values from the files.
 *
 * Preconditions:
 *   The given RVs must come from the result of unroll, and the observation matrix
 *   *must* be set up and ready to be used.
 *
 * Postconditions:
 *   All discrete observed random variables have a value that is their appropriate observed
 *   value.
 *
 * Side Effects:
 *   Modifies values of discrete observed random variables.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
setObservedRVs(vector <RV*>& rvs)
{
  // Set all discrete observed variables to their values here. Continuous
  // observed rvs get their values from elsehwere.
  vector <RV*>::iterator it;
  for (it = rvs.begin(); it!= rvs.end(); it++) {
    RV *rv = (*it);
    if (rv->discrete() && !rv->hidden()) {
      DiscRV* drv = (DiscRV*)rv;
      drv->setToObservedValue();
    }
  }
}
void
setObservedRVs(set <RV*>& rvs)
{
  // Set all discrete observed variables to their values here. Continuous
  // observed rvs get their values from elsehwere.
  set <RV*>::iterator it;
  for (it = rvs.begin(); it!= rvs.end(); it++) {
    RV *rv = (*it);
    if (rv->discrete() && !rv->hidden()) {
      DiscRV* drv = (DiscRV*)rv;
      drv->setToObservedValue();
    }
  }
}

