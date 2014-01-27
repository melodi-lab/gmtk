/*-
 * GMTK_MTCPT.cc
 *     A Multi-Dimensional dense Conditional Probability Table class.
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (C) 2001 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
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

#include "GMTK_GMParms.h"
#include "GMTK_MTCPT.h"
#include "GMTK_DiscRV.h"

#if HAVE_CONFIG_H
#include <config.h>
#endif
#if HAVE_HG_H
#include "hgstamp.h"
#endif
VCID(HGID)





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
  // this CPT is not trained and its maximum possible value is always unity.
  cachedMaxValue.set_to_one();
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
void MTCPT::setNumParents(const unsigned _nParents)
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
void MTCPT::setNumCardinality(const unsigned var, const int card)
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
  error("DeterministicCPT::allocateBasicInternalStructures() not implemented");
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
  is.read(_numParents,"Can't read DeterministicCPT numParents");

  if (_numParents >= warningNumParents)
    warning("WARNING: creating DeterministicCPT '%s' with %d parents in file '%s' line %d",
	    _numParents,name().c_str(),is.fileName(),is.lineNo());

  cardinalities.resize(_numParents);
  // read the cardinalities
  for (unsigned i=0;i<_numParents;i++) {
    is.read(cardinalities[i],"Can't read DeterministicCPT parent cardinality");
    if (cardinalities[i] <= 0)
      error("ERROR: reading file '%s' line %d, DeterministicCPT '%s' trying to use 0 or negative (%d) cardinality table, position %d.",
	    is.fileName(),is.lineNo(),name().c_str(),cardinalities[i],i);
  }

  // read the self cardinalities
  is.read(_card,"Can't read DeterministicCPT self cardinality");
  if (_card <= 0)
    error("ERROR: reading file '%s' line %d, DeterministicCPT '%s' trying to use 0 or negative (%d) cardinality table, position %d.",
	  is.fileName(),is.lineNo(),name().c_str(),_card,_numParents);


  // Finally read in the ID of the decision tree
  // that maps from parent values to an integer specifying
  // the sparse CPT. 
  string str;
  is.read(str);
  if (GM_Parms.dtsMap.find(str) ==  GM_Parms.dtsMap.end()) 
    error("ERROR: reading file '%s' line %d, DeterministicCPT '%s' specifies DT name '%s' that does not exist",
	  is.fileName(),is.lineNo(),_name.c_str(),str.c_str());
  dtIndex = GM_Parms.dtsMap[str];

  dt = GM_Parms.dts[dtIndex];
  
  if (_numParents != dt->numFeatures())
    error("ERROR: reading file '%s' line %d, DeterministicCPT '%s' with %d parents specifies DT '%s' with %d features that does not match",
	  is.fileName(),is.lineNo(),_name.c_str(),_numParents,str.c_str(),dt->numFeatures());

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
  os.write(_numParents,"DeterministicCPT::write numParents");
  os.writeComment("number parents");os.nl();
  for (unsigned i=0;i<_numParents;i++) {
    os.write(cardinalities[i],"DeterministicCPT::write cardinality");
  }
  os.write(card(),"DeterministicCPT::write cardinality");
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

  accumulatedProbability = 0.0;  
}


void
MTCPT::emIncrement(logpr prob, vector <RV*>& parents,RV* rv)
{
  if (!emAmTrainingBitIsSet())
    return;

  if(!emOnGoingBitIsSet())
    emStartIteration();

  // this is an MTCPT, so rv must be discrete.a
  assert ( rv -> discrete() );

  DiscRV* drv = (DiscRV*)rv;
  // make sure, by checking that drv's curCPT points to this.
  assert ( drv -> curCPT == this );

  accumulatedProbability += prob;
}

void
MTCPT::emEndIteration()
{
  if (!emAmTrainingBitIsSet())
    return;

  if (!emOnGoingBitIsSet())
    return;

  accumulatedProbability.floor();
  if (accumulatedProbability.zero()) {
    warning("WARNING: DeterministicCPT named '%s' did not receive any accumulated probability in EM iteration",name().c_str());
  }

  // stop EM
  emClearOnGoingBit();

}

void
MTCPT::emSwapCurAndNew()
{
  if (!emAmTrainingBitIsSet())
    return;

  if (!emSwappableBitIsSet())
    return;
  emClearSwappableBit();
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
