/*-
 * GMTK_MSCPT.cc
 *     A Multi-Dimensional dense Conditional Probability Table class.
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2001, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle, and Jeff Bilmes make no representations about
 * the suitability of this software for any purpose.  It is provided
 * "as is" without express or implied warranty.
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

#include "GMTK_MSCPT.h"
#include "GMTK_GMParms.h"
#include "GMTK_DiscreteRandomVariable.h"
#include "GMTK_GMParms.h"


VCID("$Header$");




////////////////////////////////////////////////////////////////////
//        General create, read, destroy routines 
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * MSCPT::MSCPT()
 *      Constructor
 *
 * Results:
 *      Constructs the object.
 *
 * Side Effects:
 *      None so far.
 *
 *-----------------------------------------------------------------------
 */
MSCPT::MSCPT()
  : CPT(di_MSCPT)
{
}


/*-
 *-----------------------------------------------------------------------
 * MSCPT::setNumParents()
 *      Just sets the number of parents and resizes arrays as appropriate.
 *
 * Results:
 *      no results.
 *
 * Side Effects:
 *      Will change internal arrays of this object.
 *
 *-----------------------------------------------------------------------
 */
void MSCPT::setNumParents(const int _nParents)
{
  CPT::setNumParents(_nParents);
  bitmask &= ~bm_basicAllocated;
}


/*-
 *-----------------------------------------------------------------------
 * MSCPT::setNumCardinality(var,card)
 *      sets the cardinality of var to card
 *
 * Results:
 *      no results.
 *
 * Side Effects:
 *      Will change internal array content of this object.
 *
 *-----------------------------------------------------------------------
 */
void MSCPT::setNumCardinality(const int var, const int card)
{

  CPT::setNumCardinality(var,card);
  // assume that the basic stuff is not allocated.
  bitmask &= ~bm_basicAllocated;
}


/*-
 *-----------------------------------------------------------------------
 * MSCPT::allocateBasicInternalStructures()
 *      Allocates remainder of internal data structures assuming
 *      that numParents and cardinalities are called.
 *
 * Results:
 *      no results.
 *
 * Side Effects:
 *      Will change internal content of this object.
 *
 *-----------------------------------------------------------------------
 */
void MSCPT::allocateBasicInternalStructures()
{
  error("MSCPT::allocateBasicInternalStructures() not implemented");
  // basic stuff is now allocated.
  bitmask |= bm_basicAllocated;
}




/*-
 *-----------------------------------------------------------------------
 * MSCPT::read(is)
 *      read in the table.
 * 
 * Results:
 *      No results.
 *
 * Side Effects:
 *      Changes the 'mdcpt' member function in the object.
 *
 *-----------------------------------------------------------------------
 */

void
MSCPT::read(iDataStreamFile& is)
{

  NamedObject::read(is);
  is.read(_numParents,"MSCPT::read numParents");

  if (_numParents < 0) 
    error("ERROR: reading file '%s', MSCPT '%s' trying to use negative (%d) num parents.",is.fileName(),name().c_str(),_numParents);
  if (_numParents >= warningNumParents)
    warning("WARNING: creating MSCPT with %d parents in file '%s'",_numParents,
	    is.fileName());

  cardinalities.resize(_numParents);
  // read the cardinalities
  for (unsigned i=0;i<_numParents;i++) {
    is.read(cardinalities[i],"MSCPT::read cardinality");
    if (cardinalities[i] <= 0)
      error("ERROR: reading file '%s', MDCPT '%s' trying to use 0 or negative (%d) cardinality table, position %d.",
	    is.fileName(),name().c_str(),cardinalities[i],i);
  }

  // read the self cardinalities
  is.read(_card,"MSCPT::read cardinality");
  if (_card <= 0)
    error("ERROR: reading file '%s', MSCPT '%s' trying to use 0 or negative (%d) cardinality table, position %d.",
	  is.fileName(),name().c_str(),_card,_numParents);


  // read the name of the decision tree
  string str;
  is.read(str);
  if (GM_Parms.dtsMap.find(str) ==  GM_Parms.dtsMap.end()) 
      error("Error: MTCPT '%s' specifies DT name '%s' that does not exist",
	    _name.c_str(),str.c_str());
  dtIndex = GM_Parms.dtsMap[str];

  dt = GM_Parms.dts[dtIndex];
  
  if (_numParents != dt->numFeatures())
      error("ERROR: reading file '%s', MTCPT '%s' with %d parents specifies DT '%s' with %d features that does not match",
	    is.fileName(),
	    _name.c_str(),_numParents,str.c_str(),dt->numFeatures());

  // read the name of the collection giving a list of SPMFs
  string clstr;
  is.read(clstr);
  if (GM_Parms.nclsMap.find(clstr) ==  GM_Parms.nclsMap.end()) 
      error("Error: MTCPT '%s' specifies collection name '%s' that does not exist",
	    _name.c_str(),clstr.c_str());
  ncl = GM_Parms.ncls[GM_Parms.nclsMap[clstr]];
  ncl->fillSpmfTable();

  //////////////////////////////////////////////////////////
  // now go through the dt and make sure each dt leaf (to
  // the extent possible) to a valid collection entry
  RngDecisionTree::iterator it = dt->begin();
  do {
    if (!it.valueNode()) {
      // run-time check will check for this condition.
    } else {
      // if it is a constant value, might as well
      // do the check here since that might inform user
      // of problems earlier on than during run-time.
      const unsigned v = it.value();
      if (!ncl->validSpmfIndex(v))
	error("ERROR: reading file '%s', MSCPT '%s' has DT '%s' with leaf value %d is out of range of collection '%s'",
	      is.fileName(),
	      name().c_str(),
	      str.c_str(),
	      v,
	      ncl->name().c_str());
    }
  } while (++it != dt->end());

  // Check that each Sparse1DPMF pointed to by the colleciton has the same
  // cardinality as self. This will eliminate the need for run-time checks of
  // this condition since all DT leaf values, even if they have integer expressions,
  // must go through a collection.
  // also compute _averageCardinality using a heuristic
  // under an assumption of uniform distribution of hits in the collection.
  _averageCardinality = _maxCardinality = 0;
  for (unsigned u=0;u<ncl->spmfSize();u++) {
    if (ncl->spmf(u)->card() != card()) {
      error("ERROR: reading file '%s', MSCPT '%s' uses collection '%s' referring to SPMF '%s' (position %d within collection) with card %d but MSCPT needs card %d",
	    is.fileName(),
	    name().c_str(),
	    ncl->name().c_str(),
	    ncl->spmf(u)->name().c_str(),
	    u,
	    ncl->spmf(u)->card(),
	    card());
    }
    _averageCardinality += ncl->spmf(u)->length();
    if ((unsigned)ncl->spmf(u)->length() > _maxCardinality)
      _maxCardinality = ncl->spmf(u)->length();

  }
  _averageCardinality /= ncl->spmfSize();


  setBasicAllocatedBit();
}


/*-
 *-----------------------------------------------------------------------
 * MSCPT::write(os)
 *      write out data to file 'os'. 
 * 
 * Results:
 *      No results.
 *
 * Side Effects:
 *      No effects other than  moving the file pointer of os.
 *
 *-----------------------------------------------------------------------
 */
void
MSCPT::write(oDataStreamFile& os)
{
  NamedObject::write(os);
  os.write(_numParents,"MSCPT::write numParents");
  os.writeComment("number parents");os.nl();
  for (unsigned i=0;i<_numParents;i++) {
    os.write(cardinalities[i],"MSCPT::write cardinality");
  }
  os.write(card(),"MSCPT::write cardinality");
  os.writeComment("cardinalities");
  os.nl();
  os.write(dt->name());
  os.nl();
  os.write(ncl->name());
  os.nl();
}


////////////////////////////////////////////////////////////////////
//        Probability Evaluation
////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////
//        Misc Support
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * MSCPT::randomSample()
 *      Takes a random sample given current parent values.
 *  
 * Results:
 *      the sample
 *
 * Side Effects:
 *      none.
 *
 *-----------------------------------------------------------------------
 */
int
MSCPT::randomSample()
{
  assert ( bitmask & bm_basicAllocated );
  
  iterator it = begin();
  logpr uniform = rnd.drand48();
  logpr sum = 0.0;
  do {
    sum += it.probVal;
    if (uniform <= sum)
      break;
    it++;
  } while (!end(it));
  
  return it.val();
}


/*-
 *-----------------------------------------------------------------------
 * MSCPT::normalize()
 *      Re-normalize all the distributions
 *  
 * Results:
 *      No results.
 *
 * Side Effects:
 *      Changes the values of all tables that this uses.
 *      NOTE: normalizing here will also normalize any
 *      Sparse1DPMFs that are shared with this.
 *
 *-----------------------------------------------------------------------
 */
void
MSCPT::normalize()
{
  assert ( bitmask & bm_basicAllocated );
  RngDecisionTree::iterator it = dt->begin();
  do {
    const int v = it.value();

    // TODO: this runtime check isn't needed if
    // we do not allow integer expression leaf nodes for MSCPT
    // DTs and if we check at read time above.
    if (!ncl->validSpmfIndex(v)) {
      error("ERROR: MSCPT '%s' uses a DT named '%s' that resulted in an invalid index (%d) in collection '%s' to a SPMF, of which there are only %d.\n",
	    name().c_str(),
	    dt->name().c_str(),
	    v,
	    ncl->name().c_str(),
	    ncl->spmfSize());
    }

    ncl->spmf(v)->normalize();
    it++;
  } while (it != dt->end());
}


/*-
 *-----------------------------------------------------------------------
 * MSCPT::makeRandom()
 *      Assign random but valid values to the CPT distribution.
 *  
 * Results:
 *      No results.
 *
 * Side Effects:
 *      Changes the values of all tables that this uses.
 *      NOTE: normalizing here will also normalize any
 *      Sparse1DPMFs that are shared with this.
 *
 *-----------------------------------------------------------------------
 */
void
MSCPT::makeRandom()
{
  assert ( bitmask & bm_basicAllocated );
  if (!emAmTrainingBitIsSet())
    return;

  RngDecisionTree::iterator it = dt->begin();
  do {
    const int v = it.value();
    if (!ncl->validSpmfIndex(v)) {
      error("ERROR: MSCPT '%s' uses a DT named '%s' that resulted in an invalid index (%d) in collection '%s' to a SPMF, of which there are only %d.\n",
	    name().c_str(),dt->name().c_str(),v,
	    ncl->name().c_str(),ncl->spmfSize());
    }
    ncl->spmf(v)->normalize();
    it++;
  } while (it != dt->end());
}




/*-
 *-----------------------------------------------------------------------
 * MSCPT::makeUnifrom()
 *      Have distribution be entirely uniform.
 *  
 * Results:
 *      No results.
 *
 * Side Effects:
 *      Changes the values of all tables that this uses.
 *      NOTE: normalizing here will also normalize any
 *      Sparse1DPMFs that are shared with this.
 *
 *-----------------------------------------------------------------------
 */
void
MSCPT::makeUniform()
{
  assert ( bitmask & bm_basicAllocated );
  if (!emAmTrainingBitIsSet())
    return;

  RngDecisionTree::iterator it = dt->begin();
  do {
    const int v = it.value();
    if (!ncl->validSpmfIndex(v)) {
      error("ERROR: MSCPT '%s' uses a DT named '%s' that resulted in an invalid index (%d) in collection '%s' to a SPMF, of which there are only %d.\n",
	    name().c_str(),dt->name().c_str(),v,
	    ncl->name().c_str(),
	    ncl->spmfSize());
    }
    ncl->spmf(v)->normalize();
    it++;
  } while (it != dt->end());
}



////////////////////////////////////////////////////////////////////
//        EM Routines
////////////////////////////////////////////////////////////////////

 
void
MSCPT::emStartIteration()
{
  if (!emAmTrainingBitIsSet())
    return;

  if(emOnGoingBitIsSet())
    return; // already done

  if (!emEmAllocatedBitIsSet()) {
    emSetEmAllocatedBit();
  }
  // EM iteration is now going.
  emSetOnGoingBit();
  emSetSwappableBit();

  // spmf->emStartIteration();

  accumulatedProbability = 0.0;  
}


void
MSCPT::emIncrement(logpr prob,RandomVariable* rv)
{
  if (!emAmTrainingBitIsSet())
    return;

  if (!emOnGoingBitIsSet())
    emStartIteration();

  // this is an MSCPT, so rv must be discrete.a
  assert ( rv -> discrete );

  DiscreteRandomVariable* drv = (DiscreteRandomVariable*)rv;
  // make sure, by checking that drv's curCPT points to this.
  assert ( drv -> curCPT == this );

  // 
  // TODO: This needs to be factored out of the inner most
  // loop!
  becomeAwareOfParentValues(*(drv->curConditionalParents));

  spmf->emIncrement(prob,drv->val);

  accumulatedProbability += prob;
}

void
MSCPT::emEndIteration()
{
  if (!emAmTrainingBitIsSet())
    return;

  if (!emOnGoingBitIsSet())
    return;

  accumulatedProbability.floor();
  if (accumulatedProbability.zero()) {
    warning("WARNING: MSCPT named '%s' did not receive any accumulated probability in EM iteration",name().c_str());
  }

  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  // Now we need to go through all Sparse1DPMFs that this
  // MSCPT might use and call emEndIteration(). We assume
  // however that this is done somewhere else from the
  // global object. The reason we are not able to do it here
  // is that our decision tree can have integer formulas
  // as leaf nodes, and we do not know at this point who all 
  // of the Sparse1DPMFs that are being used by this MDCPT
  // are. Instead, we make the assumption above about
  // the global object (GMTK_GMParms) doing the work.
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////

  // stop EM
  emClearOnGoingBit();
}

void
MSCPT::emSwapCurAndNew()
{
  if (!emAmTrainingBitIsSet())
    return;

  if (!emSwappableBitIsSet())
    return;
  
  // We don't call our member's swap function since
  // they are avaiable only as leaves of a decision tree
  // which depend on parent values. We assume that this
  // happens via the global object.

  emClearSwappableBit();
}



////////////////////////////////////////////////////////////////////
//        Test Driver
////////////////////////////////////////////////////////////////////

#ifdef MAIN

#include "fileParser.h"

int
main()
{

}


#endif
