/*-
 * GMTK_ContinuousRandomVariable.cc
 *     Support code for continuous random variables.

 * Written by Jeff Bilmes<bilmes@ee.washington.edu>
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
 */


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <float.h>
#include <assert.h>

#include "general.h"
#include "error.h"
#include "rand.h"

#include "GMTK_ContinuousRandomVariable.h"
#include "GMTK_GMParms.h"
#include "GMTK_ObservationMatrix.h"
#include "GMTK_MixGaussiansCommon.h"
#include "GMTK_MixGaussians.h"


VCID("$Header$");

ContinuousRandomVariable::ContinuousRandomVariable(string _label)
  : RandomVariable(_label,Continuous) 
{

}

/*-
 *-----------------------------------------------------------------------
 * Function
 *      findConditionalParents()
 *     Set up conditional parents pointers and other tables.
 * 
 * Preconditions:
 *      variable must be filled in.
 *
 * Postconditions:
 *      What is true after the function is called.
 *
 * Side Effects:
 *      Changes some internal object structures such as;
 *
 * Results:
 *      What does the function return, if anything. 
 *
 *-----------------------------------------------------------------------
 */
void
ContinuousRandomVariable::findConditionalParents()
{
  cachedIntFromSwitchingState = intFromSwitchingState();

  if ( cachedIntFromSwitchingState < 0 ||
       cachedIntFromSwitchingState >= conditionalParentsList.size()) {
    error("ERROR: CRV %s:%d using DT '%s' got invliad switching position %d\n",
	  label.c_str(),timeIndex, (dtMapper == NULL?"NULL":dtMapper->name().c_str()),cachedIntFromSwitchingState);
  }

  curConditionalParents = & conditionalParentsList[cachedIntFromSwitchingState];
  curMappingOrDirect = &conditionalGaussians[cachedIntFromSwitchingState];
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
logpr
ContinuousRandomVariable::probGivenParents()
{
  // TODO: there should be some check to make sure
  // that the resulting gaussian is the correct dimensionality.

  logpr _cachedProb;
  if (curMappingOrDirect->direct) {
    _cachedProb = 
      curMappingOrDirect->gaussian->log_p
      ((unsigned)timeIndex,firstFeatureElement);
  } else {
    // need to find which gaussian this will be.
    const unsigned gaussianIndex =
      curMappingOrDirect->dtMapper->query(*curConditionalParents);

    ///////////////////////////////////////////////////////////
    // Dynamic error checking:
    // the following check needs to be here because DTs might
    // have formulas in their leaves and there is no way
    // to check this statically w/o enumerating through all possible
    // values of the parents of this RV.
    if ( gaussianIndex >= GM_Parms.mixGaussians.size()) {
      error("ERROR: random variable '%s' (time frame %d) using decision tree '%s' wants GM "
            "with index %d but there are only %d GMs",
	    label.c_str(),timeIndex,curMappingOrDirect->dtMapper->name().c_str(),
	    gaussianIndex,GM_Parms.mixGaussians.size());
    }
    ////////////////////////////////////////////////////////////

    //
    // TODO: this needs to be changed when we have
    // different types of mixtures of Gaussians.
    // printf("CRV: '%s', par val %d, gi = %d\n",
    // label.c_str(),(*curConditionalParents)[0]->val,gaussianIndex);
    _cachedProb = GM_Parms.mixGaussians[gaussianIndex]->log_p
      ((unsigned)timeIndex,firstFeatureElement);
  }
  return _cachedProb;
}

void
ContinuousRandomVariable::makeRandom()
{
  error("not implemented, this should be called to somewhere else");
}

void
ContinuousRandomVariable::makeUniform()
{
  error("not implemented, this should be called to somewhere else");
}



/*-
 *-----------------------------------------------------------------------
 * tieParametersWith()
 *      Ties the parameters of 'this' with whatever those of 'other' are. 
 *      'other' and 'this' must be identical structuraly.
 * 
 * Preconditions:
 *      other must be a fully instantiated RV with parameters, and 'this'
 *      and 'other' must be structurally identical.
 *
 * Postconditions:
 *      'this' has the identical _tied_ parameters with 'other'
 *
 * Side Effects:
 *      Changes the internal parameter data structures, but does not delete anything.
 *
 * Results:
 *      returns nothing.
 *
 *-----------------------------------------------------------------------
 */
void
ContinuousRandomVariable::tieParametersWith(RandomVariable*const _other)
{
  assert ( !(_other -> discrete) );
  ContinuousRandomVariable* other = (ContinuousRandomVariable*)_other;

  if (!identicalStructureWith(*other))
    error("Error, trying to tie parameters of RV '%s' with RV '%s' but they have different structure.",
	  label.c_str(),other->label.c_str());

  conditionalGaussians = other->conditionalGaussians;
  curMappingOrDirect = other->curMappingOrDirect;
}


/*-
 *-----------------------------------------------------------------------
 * clone()
 *      Returns a clone of self.
 * 
 * Preconditions:
 *      self must be filled in.
 *
 * Postconditions:
 *      same as preconditions.
 *
 * Side Effects:
 *      No internal effects.
 *
 * Results:
 *      returns a new random variable.
 *
 *-----------------------------------------------------------------------
 */
RandomVariable*
ContinuousRandomVariable::clone()
{
  ContinuousRandomVariable* rv = 
    (ContinuousRandomVariable*) RandomVariable::clone();
  // might as well set val, although probably won't be useful.
  rv->tieParametersWith(this);
  rv->firstFeatureElement = firstFeatureElement;
  rv->lastFeatureElement = lastFeatureElement;
  return rv;
}





/////////////////
// EM Support //
/////////////////

void
ContinuousRandomVariable::emIncrement(logpr posterior)
{
  findConditionalParents();
  if (curMappingOrDirect->direct) {
    curMappingOrDirect->gaussian->emIncrement
      (posterior,(unsigned)timeIndex,firstFeatureElement);
  } else {
    // need to find which gaussian this will be.
    const unsigned gaussianIndex =
      curMappingOrDirect->dtMapper->query(*curConditionalParents);

    ///////////////////////////////////////////////////////////
    // Dynamic error checking:
    // the following check needs to be here because DTs might
    // have formulas in their leaves and there is no way
    // to check this statically w/o enumerating through all possible
    // values of the parents of this RV.
    if ( gaussianIndex >= GM_Parms.mixGaussians.size()) {
      error("ERROR: random variable '%s' (time frame %d) using decision tree '%s' wants GM "
            "with index %d but there are only %d GMs",
	    label.c_str(),timeIndex,curMappingOrDirect->dtMapper->name().c_str(),
	    gaussianIndex,GM_Parms.mixGaussians.size());
    }

    // TODO: this needs to be changed when we have
    // different types of mixtures of Gaussians.
    GM_Parms.mixGaussians[gaussianIndex]->emIncrement
      (posterior,(unsigned)timeIndex,firstFeatureElement);
  }
}

#if 0
void
ContinuousRandomVariable::emStoreAccumulators(oDataStreamFile& ofile)
{
  error("not implemented");
}

void
ContinuousRandomVariable::emLoadAccumulators(iDataStreamFile& ifile)
{
  error("not implemented");
}


void
ContinuousRandomVariable::emAccumulateAccumulators(iDataStreamFile& ifile)
{
  error("not implemented");
}

#endif
