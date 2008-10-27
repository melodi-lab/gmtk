/*
 * GMTK_ObsContRV.cc
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
VCID("$Header$")

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <float.h>
#include <assert.h>
#include <set>


#include "error.h"
#include "rand.h"

#include "GMTK_ObsContRV.h"
#include "GMTK_GMParms.h"
#include "GMTK_ObservationMatrix.h"
#include "GMTK_MixtureCommon.h"
#include "GMTK_Mixture.h"
#include "GMTK_RngDecisionTree.h"



/*-
 *-----------------------------------------------------------------------
 * printNameFrameValue()
 *      prints a one-line summary of the detailed information about this RV, optionally including
 *      the continuous observation (which might be long).
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
void
ObsContRV::printNameFrameValue(FILE *f,bool nl)
{
  RV::printNameFrame(f,false);
  // the current global debug level changes the way observed variables are printed.
  if (IM::messageGlb(IM::Mega+5)) {
    // then print out the observation.
    fprintf(f,"=");
    for (unsigned i=firstFeatureElement();i<=lastFeatureElement();i++) {
      // print only 1+3 significant digits for now.
      fprintf(f,"%.3e%s",
	      (*globalObservationMatrix.floatVecAtFrame(frame(), i)),
	      ((i<lastFeatureElement())?",":""));
    }
    if (nl)
      fprintf(f,"\n");
  } else
    fprintf(f,"=C%s",nls(nl));
}


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
void ObsContRV::printSelf(FILE *f,bool nl)
{
  printNameFrameValue(f,false);
  fprintf(f," observed continuous%s",nls(nl));
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
void ObsContRV::printSelfVerbose(FILE *f)
{
  fprintf(f,"Observed Continuous Random variable:\n");
  RV::printNameFrame(f,true);
  fprintf(f,"From line %d in file %s\n",rv_info.fileLineNumber,rv_info.rvFileName.c_str());
  fprintf(f,"Vector Value: (");
  for (unsigned i=firstFeatureElement();i<=lastFeatureElement();i++) {
    fprintf(f,"%d:%f,",
	    i,
	    (*globalObservationMatrix.floatVecAtFrame(frame(), i))
	    );
  }
  fprintf(f,"\n");
}


/*-
 *-----------------------------------------------------------------------
 * probGivenParents()
 *      Computes the probability given the current values of the parents.
 * 
 * Preconditions:
 *      function must be filled in.
 *
 * Postconditions:
 *      What is true after the function is called.
 *
 * Side Effects:
 *      none.
 *
 * Results:
 *      What does the function return, if anything. 
 *
 *-----------------------------------------------------------------------
 */
void
ObsContRV::probGivenParents(logpr& p)
{
  ///////////////////
  // We assume here that the resulting mixture is the correct
  // dimensionality (this is checked in GMTK_FileParser.cc, 
  // in function FileParser::associateWithDataParams(bool)

  if (curMappingOrDirect->direct) {
    p = curMappingOrDirect->mixture->log_p
      (frame(),firstFeatureElement());
    return;
  }

  // need to find which gaussian this will be.
  const unsigned gaussianIndex =
    curMappingOrDirect->mapping.dtMapper->query(allParents,this);

  ///////////////////////////////////////////////////////////
  // Dynamic error checking:
  // the following check needs to be here because DTs might
  // have formulas in their leaves and there is no way
  // to check this statically w/o enumerating through all possible
  // values of the parents of this RV.
  if (!curMappingOrDirect->mapping.collection->validMxIndex(gaussianIndex)) {
    warning("ERROR: random variable '%s' (time frame %d) using decision tree '%s' wants mixture "
	    "with index %d but there are only %d mixtures in collection '%s'",
	    name().c_str(),frame(),curMappingOrDirect->mapping.dtMapper->name().c_str(),
	    gaussianIndex,
	    curMappingOrDirect->mapping.collection->mxSize(),
	    curMappingOrDirect->mapping.collection->name().c_str());
    fprintf(stderr,"Parents configuration :");
    printRVSetAndValues(stderr,allParents);
    error("");
  }
  ////////////////////////////////////////////////////////////
  
  // printf("CRV: '%s', par val %d, gi = %d\n",
  // label.c_str(),(*curConditionalParents)[0]->val,gaussianIndex);
  p = curMappingOrDirect->mapping.collection->mx(gaussianIndex)->log_p
    (frame(),firstFeatureElement());
}




/////////////////
// EM Support //
/////////////////


void
ObsContRV::emIncrement(logpr posterior)
{


  if (curMappingOrDirect->direct) {
    curMappingOrDirect->mixture->emIncrement
      (posterior,frame(),firstFeatureElement());
  } else {

    // need to find which mixture this will be.
    const unsigned mixtureIndex =
      curMappingOrDirect->mapping.dtMapper->query(allParents,this);

    ///////////////////////////////////////////////////////////
    // Dynamic error checking:
    // the following check needs to be here because DTs might
    // have formulas in their leaves and there is no way
    // to check this statically w/o enumerating through all possible
    // values of the parents of this RV.
    if (!curMappingOrDirect->mapping.collection->validMxIndex(mixtureIndex)) {
      error("ERROR: random variable '%s' (time frame %d) using decision tree '%s' wants mixture "
	    "with index %d but there are only %d mixtures in collection '%s'",
	    name().c_str(),frame(),curMappingOrDirect->mapping.dtMapper->name().c_str(),
	    mixtureIndex,
	    curMappingOrDirect->mapping.collection->mxSize(),
	    curMappingOrDirect->mapping.collection->name().c_str());
      fprintf(stderr,"Parents configuration :");
      printRVSetAndValues(stderr,allParents);
      error("");
    }

    // TODO: this needs to be changed when we have
    // different types of mixtures.
    curMappingOrDirect->mapping.collection->mx(mixtureIndex)->emIncrement
      (posterior,frame(),firstFeatureElement());
  }
}




/*-
 *-----------------------------------------------------------------------
 * cloneRVShell()
 *      clones a shell of the current random variable (see< GMTK_RV.h for docs)
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
ObsContRV* ObsContRV::cloneRVShell()
{
  ObsContRV*rv = (ObsContRV*)ContRV::cloneRVShell();
  rv->conditionalMixtures = conditionalMixtures;
  rv->curMappingOrDirect = curMappingOrDirect;
  return rv;
}

