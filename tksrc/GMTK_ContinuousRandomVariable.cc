
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
#include "GMTK_MixtureCommon.h"
#include "GMTK_Mixture.h"


VCID("$Header$");

ContinuousRandomVariable::ContinuousRandomVariable(RVInfo& _rv_info,
						   string _label)
  : RandomVariable(_rv_info,_label,Continuous) 
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
    error("ERROR: CRV %s:%d using DT '%s' yielded an invalid switching position %d. Must be between 0 and %d.\n",
	  label.c_str(),timeIndex, (dtMapper == NULL?"NULL":dtMapper->name().c_str()),cachedIntFromSwitchingState,conditionalParentsList.size());
  }

  curConditionalParents = & conditionalParentsList[cachedIntFromSwitchingState];
  curMappingOrDirect = &conditionalMixtures[cachedIntFromSwitchingState];
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
  ///////////////////
  // We assume here that the resulting mixture is the correct
  // dimensionality (this is checked in GMTK_FileParser.cc, 
  // in function FileParser::associateWithDataParams(bool)

  logpr _cachedProb;
  if (curMappingOrDirect->direct) {
    _cachedProb = 
      curMappingOrDirect->mixture->log_p
      ((unsigned)timeIndex,firstFeatureElement);
  } else {
    // need to find which gaussian this will be.
    const unsigned gaussianIndex =
      curMappingOrDirect->mapping.dtMapper->query(*curConditionalParents);

    ///////////////////////////////////////////////////////////
    // Dynamic error checking:
    // the following check needs to be here because DTs might
    // have formulas in their leaves and there is no way
    // to check this statically w/o enumerating through all possible
    // values of the parents of this RV.
    if (!curMappingOrDirect->mapping.collection->validMxIndex(gaussianIndex)) {
      error("ERROR: random variable '%s' (time frame %d) using decision tree '%s' wants GM "
            "with index %d but there are only %d GMs in collection '%s'",
	    label.c_str(),timeIndex,curMappingOrDirect->mapping.dtMapper->name().c_str(),
	    gaussianIndex,
	    curMappingOrDirect->mapping.collection->mxSize(),
	    curMappingOrDirect->mapping.collection->name().c_str());
    }
    ////////////////////////////////////////////////////////////

    //
    // TODO: this needs to be changed when we have
    // different types of mixtures.
    // printf("CRV: '%s', par val %d, gi = %d\n",
    // label.c_str(),(*curConditionalParents)[0]->val,gaussianIndex);
    _cachedProb = 
      curMappingOrDirect->mapping.collection->mx(gaussianIndex)->log_p
      ((unsigned)timeIndex,firstFeatureElement);
  }
  if (wtStatus != wt_NoWeight) {
    if (wtStatus == wt_Constant)
      _cachedProb.valref() *= wtWeight;
    else // get weight from observation matrix at current time frame
      _cachedProb.valref() *= 
	(*globalObservationMatrix.floatVecAtFrame((unsigned)timeIndex, wtFeatureElement));
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
ContinuousRandomVariable::tieParametersWith(RandomVariable*const _other,
					    bool checkStructure)
{
  assert ( !(_other -> discrete) );
  ContinuousRandomVariable* other = (ContinuousRandomVariable*)_other;

  if (checkStructure && !identicalStructureWith(*other))
    error("Error, trying to tie parameters of RV '%s' with RV '%s' but they have different structure.",
	  label.c_str(),other->label.c_str());

  conditionalMixtures = other->conditionalMixtures;
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



/*-
 *-----------------------------------------------------------------------
 * cloneWithoutParents()
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
ContinuousRandomVariable::cloneWithoutParents()
{
  ContinuousRandomVariable* rv = 
    (ContinuousRandomVariable*) RandomVariable::cloneWithoutParents();
  // might as well set val, although probably won't be useful.
  rv->tieParametersWith(this,false);
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
    curMappingOrDirect->mixture->emIncrement
      (posterior,(unsigned)timeIndex,firstFeatureElement);
  } else {
    // need to find which mixture this will be.
    const unsigned mixtureIndex =
      curMappingOrDirect->mapping.dtMapper->query(*curConditionalParents);

    ///////////////////////////////////////////////////////////
    // Dynamic error checking:
    // the following check needs to be here because DTs might
    // have formulas in their leaves and there is no way
    // to check this statically w/o enumerating through all possible
    // values of the parents of this RV.
    if (!curMappingOrDirect->mapping.collection->validMxIndex(mixtureIndex)) {
      error("ERROR: random variable '%s' (time frame %d) using decision tree '%s' wants GM "
            "with index %d but there are only %d GMs in collection '%s'",
	    label.c_str(),timeIndex,curMappingOrDirect->mapping.dtMapper->name().c_str(),
	    mixtureIndex,
	    curMappingOrDirect->mapping.collection->mxSize(),
	    curMappingOrDirect->mapping.collection->name().c_str());
    }

    // TODO: this needs to be changed when we have
    // different types of mixtures.
    curMappingOrDirect->mapping.collection->mx(mixtureIndex)->emIncrement
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
