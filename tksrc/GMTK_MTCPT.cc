/*-
 * GMTK_MTCPT.cc
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

#include "GMTK_MTCPT.h"
#include "GMTK_GMParms.h"
#include "GMTK_DiscreteRandomVariable.h"

VCID("$Header$");




////////////////////////////////////////////////////////////////////
//        General create, read, destroy routines 
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * MTCPT::MTCPT()
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
MTCPT::MTCPT()
  : CPT(di_MTCPT)
{
}


/*-
 *-----------------------------------------------------------------------
 * MTCPT::setNumParents()
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
void MTCPT::setNumParents(const int _nParents)
{
  CPT::setNumParents(_nParents);
  bitmask &= ~bm_basicAllocated;
}


/*-
 *-----------------------------------------------------------------------
 * MTCPT::setNumCardinality(var,card)
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
void MTCPT::setNumCardinality(const int var, const int card)
{

  CPT::setNumCardinality(var,card);
  // assume that the basic stuff is not allocated.
  bitmask &= ~bm_basicAllocated;
}


/*-
 *-----------------------------------------------------------------------
 * MTCPT::allocateBasicInternalStructures()
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
void MTCPT::allocateBasicInternalStructures()
{
  error("MTCPT::allocateBasicInternalStructures() not implemented");
  // basic stuff is now allocated.
  bitmask |= bm_basicAllocated;
}




/*-
 *-----------------------------------------------------------------------
 * MTCPT::read(is)
 *      read in the table.
 * 
 * Results:
 *      No results.
 *
 * Side Effects:
 *      Changes the 'mscpt' member function in the object.
 *
 *-----------------------------------------------------------------------
 */

void
MTCPT::read(iDataStreamFile& is)
{

  NamedObject::read(is);
  is.read(_numParents,"MTCPT::read numParents");

  if (_numParents < 0) 
    error("ERROR: reading file '%s', MTCPT '%s' trying to use negative (%d) num parents.",is.fileName(),name().c_str(),_numParents);
  if (_numParents >= warningNumParents)
    warning("WARNING: creating MTCPT '%s' with %d parents in file '%s'",
	    _numParents,name().c_str(),is.fileName());

  cardinalities.resize(_numParents);
  // read the cardinalities
  for (unsigned i=0;i<_numParents;i++) {
    is.read(cardinalities[i],"MTCPT::read cardinality");
    if (cardinalities[i] <= 0)
      error("ERROR: reading file '%s', MDCPT '%s' trying to use 0 or negative (%d) cardinality table, position %d.",
	    is.fileName(),name().c_str(),cardinalities[i],i);
  }

  // read the self cardinalities
  is.read(_card,"MTCPT::read cardinality");
  if (_card <= 0)
    error("ERROR: reading file '%s', MTCPT '%s' trying to use 0 or negative (%d) cardinality table, position %d.",
	  is.fileName(),name().c_str(),_card,_numParents);


  // Finally read in the ID of the decision tree
  // that maps from parent values to an integer specifying
  // the sparse CPT. 
  string str;
  is.read(str);
  if (GM_Parms.dtsMap.find(str) ==  GM_Parms.dtsMap.end()) 
    error("ERROR: reading file '%s', MTCPT '%s' specifies DT name '%s' that does not exist",is.fileName(),_name.c_str(),str.c_str());
  dtIndex = GM_Parms.dtsMap[str];

  dt = GM_Parms.dts[dtIndex];
  
  if (_numParents != dt->numFeatures())
    error("ERROR: reading file '%s', MTCPT '%s' with %d parents specifies DT '%s' with %d features that does not match",is.fileName(),
	    _name.c_str(),_numParents,str.c_str(),dt->numFeatures());

  // 
  setBasicAllocatedBit();
}


/*-
 *-----------------------------------------------------------------------
 * MTCPT::write(os)
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
MTCPT::write(oDataStreamFile& os)
{
  NamedObject::write(os);
  os.write(_numParents,"MTCPT::write numParents");
  os.writeComment("number parents");os.nl();
  for (unsigned i=0;i<_numParents;i++) {
    os.write(cardinalities[i],"MTCPT::write cardinality");
  }
  os.write(card(),"MTCPT::write cardinality");
  os.writeComment("cardinalities");
  os.nl();
  os.write(dt->name());
  os.nl();
}



////////////////////////////////////////////////////////////////////
//        EM Routines
////////////////////////////////////////////////////////////////////

 
void
MTCPT::emStartIteration()
{
  if(emOnGoingBitIsSet())
    return; // already done

  if (!emEmAllocatedBitIsSet()) {
    emSetEmAllocatedBit();
  }
  // EM iteration is now going.
  emSetOnGoingBit();
  emSetSwappableBit();

  accumulatedProbability = 0.0;  
}


void
MTCPT::emIncrement(logpr prob,RandomVariable* rv)
{
  emStartIteration();

  // this is an MTCPT, so rv must be discrete.a
  assert ( rv -> discrete );

  DiscreteRandomVariable* drv = (DiscreteRandomVariable*)rv;
  // make sure, by checking that drv's curCPT points to this.
  assert ( drv -> curCPT == this );

  // TODO: This needs to be factored out of the inner most
  // loop! (in an MTCPT can safely remove this).
  becomeAwareOfParentValues(*(drv->curConditionalParents));

  accumulatedProbability += prob;
}

void
MTCPT::emEndIteration()
{

  if (!emOnGoingBitIsSet())
    return;

  accumulatedProbability.floor();
  if (accumulatedProbability.zero()) {
    warning("WARNING: MTCPT named '%s' did not receive any accumulated probability in EM iteration",name().c_str());
  }

  // stop EM
  emClearOnGoingBit();

}

void
MTCPT::emSwapCurAndNew()
{
  if (!emSwappableBitIsSet())
    return;
  emClearSwappableBit();
}


void
MTCPT::emStoreAccumulators(oDataStreamFile& ofile)
{
  assert ( basicAllocatedBitIsSet() );
  if ( !emEmAllocatedBitIsSet() ) {
    warning("WARNING: storing zero accumulators for MTCPT '%s'\n",
	    name().c_str());
    emStoreZeroAccumulators(ofile);
    return;
  }
  EMable::emStoreAccumulators(ofile);
}


void
MTCPT::emStoreZeroAccumulators(oDataStreamFile& ofile)
{
  assert ( basicAllocatedBitIsSet() );
  EMable::emStoreZeroAccumulators(ofile);
}

void
MTCPT::emLoadAccumulators(iDataStreamFile& ifile)
{
  assert (basicAllocatedBitIsSet());
  assert (emEmAllocatedBitIsSet());
  EMable::emLoadAccumulators(ifile);
}


void
MTCPT::emAccumulateAccumulators(iDataStreamFile& ifile)
{
  assert ( basicAllocatedBitIsSet() );
  assert ( emEmAllocatedBitIsSet() );
  EMable::emAccumulateAccumulators(ifile);
}


////////////////////////////////////////////////////////////////////
//        Test Driver
////////////////////////////////////////////////////////////////////

#ifdef MAIN

#include "fileParser.h"

GMParms GM_Parms;

int
main()
{

  MTCPT mtcpt;

  // this is a binary variable with three parents
  // the first one is binary, the second one is ternary,
  // and the third one is binary.

  mtcpt.setNumParents(3);
  mtcpt.setNumCardinality(0,2);
  mtcpt.setNumCardinality(1,3);
  mtcpt.setNumCardinality(2,2);
  mtcpt.setNumCardinality(3,3);
  mtcpt.allocateBasicInternalStructures();

  // set to random values
  mtcpt.makeRandom();
  
  // write values to a data file in ASCII.
  oDataStreamFile od ("/tmp/foo.mtcpt",false);

  mtcpt.write(od);

  // now print out some probabilities.
  vector < int > parentVals;
  vector < int > cards;

  parentVals.resize(3);
  cards.resize(3);
  
  cards[0] = cards[1] = cards[2] = 5;

  parentVals[0] = 0;
  parentVals[1] = 0;
  parentVals[2] = 0;

  for (int i =0; i<3;i++) {
    printf("Prob(%d) Given cur Par = %f\n",
	   i,mtcpt.probGivenParents(parentVals,cards,i).unlog());
  }

  parentVals[0] = 0;
  parentVals[1] = 0;
  parentVals[2] = 1;

  for (int i =0; i<3;i++) {
    printf("Prob(%d) Given cur Par = %f\n",
	   i,mtcpt.probGivenParents(parentVals,cards,i).unlog());
  }

  parentVals[0] = 0;
  parentVals[1] = 1;
  parentVals[2] = 1;

  for (int i =0; i<3;i++) {
    printf("Prob(%d) Given cur Par = %f\n",
	   i,mtcpt.probGivenParents(parentVals,cards,i).unlog());
  }

  parentVals[0] = 1;
  parentVals[1] = 2;
  parentVals[2] = 1;

  for (int i =0; i<3;i++) {
    printf("Prob(%d) Given cur Par = %f\n",
	   i,mtcpt.probGivenParents(parentVals,cards,i).unlog());
  }

  // Now iterate over valid values.
  MTCPT::iterator it = mtcpt.begin();
  do {
    printf("Prob of %d is %f\n",
	   it.val(),it.probVal.unlog());
  } while (it != mtcpt.end());


  parentVals[0] = 0;
  parentVals[1] = 0;
  parentVals[2] = 1;
  mtcpt.becomeAwareOfParentValues(parentVals,cards);

  it = mtcpt.begin();
  do {
    printf("Prob of %d is %f\n",
	   it.val(),it.probVal.unlog());
  } while (it != mtcpt.end());


}


#endif
