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
VCID("$Header$");

#include <iostream>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <string.h>

#include "GMTK_ObsDiscRV.h"


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
  fprintf(f,"observed discrete cardinality = %d%s",cardinality,nls(nl));
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
  fprintf(f,"From line %d in file %s\n",rv_info.fileLineNumber,rv_info.rvFileName.c_str());
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
      assert (globalObservationMatrix.active());
      if (frame() >= cardinality) 
	error("ERROR: RV '%s(%d)' has cardinality %d, but current frame value %d too large to store in RV with this cardinality (in segment %d of frame length %d).\n",
	      name().c_str(),
	      frame(),
	      cardinality,
	      frame(),
	      globalObservationMatrix.segmentNumber(),
	      globalObservationMatrix.numFrames());
      val = frame();
    } else if (rv_info.rvFeatureRange.filled == RVInfo::FeatureRange::fr_NumFramesIsValue) {
      assert (globalObservationMatrix.active());
      if (globalObservationMatrix.numFrames() >= cardinality) 
	error("ERROR: RV '%s(%d)' has cardinality %d, but current num frames %d too large to store in RV with this cardinality (in segment %d of frame length %d).\n",
	      name().c_str(),
	      frame(),
	      cardinality,
	      globalObservationMatrix.numFrames(),
	      globalObservationMatrix.segmentNumber(),
	      globalObservationMatrix.numFrames());
      val = globalObservationMatrix.numFrames();
    } else if (rv_info.rvFeatureRange.filled == RVInfo::FeatureRange::fr_SegmentNumIsValue) {
      assert (globalObservationMatrix.active());
      if (globalObservationMatrix.segmentNumber() >= cardinality) 
	error("ERROR: RV '%s(%d)' has cardinality %d, but current segment number %d too large to store in RV with this cardinality (in segment %d of frame length %d).\n",
	      name().c_str(),
	      frame(),
	      cardinality,
	      globalObservationMatrix.segmentNumber(),
	      globalObservationMatrix.segmentNumber(),
	      globalObservationMatrix.numFrames());
      val = globalObservationMatrix.segmentNumber();
    } else if (rv_info.rvFeatureRange.filled == RVInfo::FeatureRange::fr_NumSegmentsIsValue) {
      assert (globalObservationMatrix.active());
      if (globalObservationMatrix.numSegments() >= cardinality) 
	error("ERROR: RV '%s(%d)' has cardinality %d, but current number segments %d too large to store in RV with this cardinality (in segment %d of frame length %d).\n",
	      name().c_str(),
	      frame(),
	      cardinality,
	      globalObservationMatrix.numSegments(),
	      globalObservationMatrix.segmentNumber(),
	      globalObservationMatrix.numFrames());
      val = globalObservationMatrix.numSegments();
    } else if (rv_info.rvFeatureRange.filled == RVInfo::FeatureRange::fr_Range) {
      assert (globalObservationMatrix.active());
      // printf("getting value of random variable '%s', time index %d, el %d\n",
      // label.c_str(),timeIndex,featureElement);
      unsigned tmp = globalObservationMatrix.unsignedAtFrame(frame(),featureElement());
      if (tmp >= (unsigned)cardinality) 
	error("ERROR: RV '%s' at time index %d has cardinality %d, but feature element position %d in observation file (time %d of segment %d) has value %u.\n",
	      name().c_str(),
	      frame(),
	      cardinality,
	      featureElement(),
	      frame(),
	      globalObservationMatrix.segmentNumber(),
	      tmp);
      val = tmp;
    } else {
      // shouldn't happen.
      assert (0);
    }
    // otherwise, we keep the value set to what it was before.
    return;
}
