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
    error("MSCPT: read, trying to use negative (%d) num parents.",_numParents);
  if (_numParents >= warningNumParents)
    warning("MSCPT: read, creating MSCPT with %d parents",_numParents);
  cardinalities.resize(_numParents+1);
  // read the cardinalities
  for (unsigned i=0;i<=_numParents;i++) {
    is.read(cardinalities[i],"MSCPT::read cardinality");
    if (cardinalities[i] <= 0)
      error("MSCPT: read, trying to use 0 or negative (%d) cardinality table.",cardinalities[i]);
  }

  // Finally read in the integer ID of the decision tree
  // that maps from parent values to an integer specifying
  // the sparse CPT. 
  is.read(dtIndex);
  if (dtIndex < 0 || dtIndex >= GM_Parms.dts.size())
    error("MSCPT::read, invalid DT index %d\n",dtIndex);

  // TODO: check that the cardinalities of self match
  // with that of the dt.

  //////////////////////////////////////////////////////////
  // set the DT
  dt = GM_Parms.dts[dtIndex];

  //////////////////////////////////////////////////////////
  // now go through the dt and make sure each dt leaf
  // node points to a Sparse1DPMF that has the same
  // cardinality as self.
  RngDecisionTree<unsigned>::iterator it = dt->begin();
  do {
    const unsigned v = it.value();
    if ( v < 0 || v >= GM_Parms.sPmfs.size() )
      error("MSCPT::read, dt leaf value %d refers to invalid Sparse 1D PMF",
	    v);
    if (GM_Parms.sPmfs[v]->card() != card()) {
      error("MSCPT::read, Sparse CPT %d has card %d but CPT has card %d",
	    v,GM_Parms.sPmfs[v]->card(),card());
    }
    it++;
  } while (it != dt->end());

  bitmask |= bm_basicAllocated;
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
  for (unsigned i=0;i<=_numParents;i++) {
    os.write(cardinalities[i],"MSCPT::write cardinality");
  }
  os.writeComment("cardinalities");
  os.nl();
  os.write(dtIndex);
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
  } while (it != end());
  
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
  RngDecisionTree<unsigned>::iterator it = dt->begin();
  do {
    const int v = it.value();
    GM_Parms.msCpts[v]->normalize();
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
  RngDecisionTree<unsigned>::iterator it = dt->begin();
  do {
    const int v = it.value();
    GM_Parms.msCpts[v]->makeRandom();
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
  RngDecisionTree<unsigned>::iterator it = dt->begin();
  do {
    const int v = it.value();
    GM_Parms.msCpts[v]->makeUniform();
    it++;
  } while (it != dt->end());
}





////////////////////////////////////////////////////////////////////
//        EM Routines
////////////////////////////////////////////////////////////////////

void
MSCPT::emStartIteration()
{
  error("Not implemented");
}


void
MSCPT::emIncrement(RandomVariable* rv, logpr prob)
{
  error("Not implemented");
}

void
MSCPT::emEndIteration()
{
  error("Not implemented");
}

void
MSCPT::emSwapCurAndNew()
{
  error("Not implemented"); 
}



////////////////////////////////////////////////////////////////////
//        Test Driver
////////////////////////////////////////////////////////////////////

#ifdef MAIN

#include "fileParser.h"

int
main()
{

  MSCPT mdcpt;

  // this is a binary variable with three parents
  // the first one is binary, the second one is ternary,
  // and the third one is binary.

  mdcpt.setNumParents(3);
  mdcpt.setNumCardinality(0,2);
  mdcpt.setNumCardinality(1,3);
  mdcpt.setNumCardinality(2,2);
  mdcpt.setNumCardinality(3,3);
  mdcpt.allocateBasicInternalStructures();

  // set to random values
  mdcpt.makeRandom();
  
  // write values to a data file in ASCII.
  oDataStreamFile od ("/tmp/foo.mdcpt",false);

  mdcpt.write(od);

  // now print out some probabilities.
  vector < int > parentVals;
  parentVals.resize(3);

  parentVals[0] = 0;
  parentVals[1] = 0;
  parentVals[2] = 0;

  for (int i =0; i<3;i++) {
    printf("Prob(%d) Given cur Par = %f\n",
	   i,mdcpt.probGivenParents(parentVals,i).unlog());
  }

  parentVals[0] = 0;
  parentVals[1] = 0;
  parentVals[2] = 1;

  for (int i =0; i<3;i++) {
    printf("Prob(%d) Given cur Par = %f\n",
	   i,mdcpt.probGivenParents(parentVals,i).unlog());
  }

  parentVals[0] = 0;
  parentVals[1] = 1;
  parentVals[2] = 1;

  for (int i =0; i<3;i++) {
    printf("Prob(%d) Given cur Par = %f\n",
	   i,mdcpt.probGivenParents(parentVals,i).unlog());
  }

  parentVals[0] = 1;
  parentVals[1] = 2;
  parentVals[2] = 1;

  for (int i =0; i<3;i++) {
    printf("Prob(%d) Given cur Par = %f\n",
	   i,mdcpt.probGivenParents(parentVals,i).unlog());
  }

  // Now iterate over valid values.
  MSCPT::iterator it = mdcpt.first();
  do {
    printf("Prob of %d is %f\n",
	   it.val,it.probVal.unlog());
  } while (mdcpt.next(it));


  parentVals[0] = 0;
  parentVals[1] = 0;
  parentVals[2] = 1;
  mdcpt.becomeAwareOfParentValues(parentVals);

  it = mdcpt.first();
  do {
    printf("Prob of %d is %f\n",
	   it.val,it.probVal.unlog());
  } while (mdcpt.next(it));


}


#endif
