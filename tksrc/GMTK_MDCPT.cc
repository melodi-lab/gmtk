/*-
 * GMTK_MDCPT.cc
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


#include "GMTK_MDCPT.h"

VCID("$Header$");




////////////////////////////////////////////////////////////////////
//        General create, read, destroy routines 
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * MDCPT::MDCPT()
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
MDCPT::MDCPT()
{
}


/*-
 *-----------------------------------------------------------------------
 * MDCPT::setNumParents()
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
void MDCPT::setNumParents(const int _nParents)
{

  if (_nParents < 0) 
    error("MDCPT: setNumParents, trying to use negative (%d) num parents.",
	  _nParents);

  // assume that the basic stuff is no longer allocated.
  bitmask &= ~bm_basicAllocated;

  numParents = _nParents;
  cardinalities.resize(numParents+1);
  cumulativeCardinalities.resize(numParents);

}


/*-
 *-----------------------------------------------------------------------
 * MDCPT::setNumCardinality(var,card)
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
void MDCPT::setNumCardinality(const int var, const int card)
{

  if (var < 0)
    error("MDCPT: setNumCardinality, trying to use negative (%d) var.",
	  var);
  if (var > numParents) 
    error("MDCPT: setNumCardinality, trying to use illegal (%d) var.",
	  var);
  if (card <= 0)
    error("MDCPT: setNumCardinality, trying to use illegal (%d) card.",
	  card);

  // assertion should be satisifed by the way that cardinalities
  // is allocated allong with setting num parents.
  assert ( var < cardinalities.len() );

  cardinalities[var] = card;

  // assume that the basic stuff is not allocated.
  bitmask &= ~bm_basicAllocated;
}


/*-
 *-----------------------------------------------------------------------
 * MDCPT::allocateBasicInternalStructures()
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
void MDCPT::allocateBasicInternalStructures()
{
  int numValues = 1;
  for (int i=0;i<=numParents;i++) {
    numValues *= cardinalities[i];
  }

  mdcpt.resize(numValues);

  if (numParents > 0) {
    cumulativeCardinalities[numParents-1] = cardinalities[numParents];
    for (int i=numParents-2; i>=0; i--) {
      cumulativeCardinalities[i] = 
	cumulativeCardinalities[i+1]*cardinalities[i+1];
    }
  }

  // basic stuff is now allocated.
  bitmask |= bm_basicAllocated;
}




/*-
 *-----------------------------------------------------------------------
 * MDCPT::read(is)
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
MDCPT::read(iDataStreamFile& is)
{

  is.read(numParents,"MDCPT::read numParents");

  if (numParents < 0) 
    error("MDCPT: read, trying to use negative (%d) num parents.",numParents);
  if (numParents >= warningNumParents)
    warning("MDCPT: read, creating MDCPT with %d parents",numParents);

  cardinalities.resize(numParents+1);
  cumulativeCardinalities.resize(numParents);

  // read the cardinalities
  int numValues = 1;
  for (int i=0;i<=numParents;i++) {
    is.read(cardinalities[i],"MDCPT::read cardinality");
    if (cardinalities[i] <= 0)
      error("MDCPT: read, trying to use 0 or negative (%d) cardinality table.",cardinalities[i]);
    numValues *= cardinalities[i];
  }

  // cumulativeCardinalities gets set to the
  // reverse cumulative cardinalities of the random
  // variables.
  if (numParents > 0) {
    cumulativeCardinalities[numParents-1] = cardinalities[numParents];
    for (int i=numParents-2;i>=0;i--) {
      cumulativeCardinalities[i] = 
	cumulativeCardinalities[i+1]*cardinalities[i+1];
    }
  }

  // Finally read in the probability values (stored as doubles).
  // NOTE: We could check that things sum to approximately 1 here, if
  // we didn't use a large 1D loop. 
  mdcpt.resize(numValues);
  for (int i=0;i<numValues;i++) {
    double val;
    is.readDouble(val,"MDCPT::read, reading value");
    if (val < 0 || val > 1)
      error("MDCPT: read, invalid pmf value (%g)",val);
    mdcpt[i] = val;
  }
  bitmask |= bm_basicAllocated;
}


/*-
 *-----------------------------------------------------------------------
 * MDCPT::write(os)
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
MDCPT::write(oDataStreamFile& os)
{
  os.write(numParents,"MDCPT::write numParents"); 
  os.writeComment("number parents");os.nl();
  for (int i=0;i<=numParents;i++) {
    os.write(cardinalities[i],"MDCPT::write cardinality");
  }
  os.writeComment("cardinalities");
  os.nl();

  // Finally write in the probability values (stored as doubles).
  // NOTE: We could check that things sum to approximately 1 here, if
  // we didn't use a large 1D loop. 
  int childCard = cardinalities[numParents];
  for (int i=0;i<mdcpt.len();i++) {
    os.writeDouble(mdcpt[i].unlog(),"MDCPT::write, writing value");
    childCard --;
    if (childCard == 0) {
      os.nl();
      childCard = cardinalities[numParents];
    }
  }
}


////////////////////////////////////////////////////////////////////
//        Probability Evaluation
////////////////////////////////////////////////////////////////////



/*-
 *-----------------------------------------------------------------------
 * MDCPT::setParentValues()
 *      Adjusts the current structure so that subsequent calls of
 *      probability routines will be conditioned on the given
 *      assigment to parent values.
 *  
 * Results:
 *      No results.
 *
 * Side Effects:
 *      Changes the mdcpt_ptr
 *
 *-----------------------------------------------------------------------
 */
void
MDCPT::setParentValues( sArray<int>& parentValues)
{

  assert ( parentValues.len() == numParents );
  assert ( bitmask & bm_basicAllocated );
  
  int offset = 0;
  for (int i = 0; i < numParents; i++) {
    if (parentValues[i] < 0 || parentValues[i] >= cardinalities[i]) 
      error("MDCPT:setParentValues: Invalid parent value for parent %d, parentValue = %d but card = %d\n",i,parentValues[i],cardinalities[i]);
    offset += parentValues[i]*cumulativeCardinalities[i];
  }
  mdcpt_ptr = mdcpt.ptr + offset;

}

////////////////////////////////////////////////////////////////////
//        Misc Support
////////////////////////////////////////////////////////////////////




/*-
 *-----------------------------------------------------------------------
 * MDCPT::makeRandom()
 *      Assign random but valid values to the CPT distribution.
 *  
 * Results:
 *      No results.
 *
 * Side Effects:
 *      Changes the mdcpt_ptr
 *
 *-----------------------------------------------------------------------
 */
void
MDCPT::makeRandom()
{
  assert ( bitmask & bm_basicAllocated );

  // Use the inherent structure of the multi-D array
  // so to loop over the final distributions on the child.

  const int child_card = cardinalities[numParents];
  const int num_parent_assignments = mdcpt.len()/child_card;
  logpr *loc_ptr = mdcpt.ptr;
  for (int parent_assignment =0; 
       parent_assignment < num_parent_assignments; 
       parent_assignment ++) {
    logpr sum = 0.0;
    logpr *tmp_loc_ptr = loc_ptr;
    for (int i=0;i<child_card;i++) {
      logpr tmp = rnd.drand48();
      sum += tmp;
      *tmp_loc_ptr ++ = tmp;
    }
    tmp_loc_ptr = loc_ptr;
    for (int i=0;i<child_card;i++) {
      *tmp_loc_ptr /= sum;
      tmp_loc_ptr++;
    }
    loc_ptr += child_card;
  }
}


////////////////////////////////////////////////////////////////////
//        Test Driver
////////////////////////////////////////////////////////////////////

#ifdef MAIN

#include "fileParser.h"

int
main()
{

  MDCPT mdcpt;

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
  sArray < int > parentVals;
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
  MDCPT::iterator it = mdcpt.first();
  do {
    printf("Prob of %d is %f\n",
	   it.val,it.probVal.unlog());
  } while (mdcpt.next(it));


  parentVals[0] = 0;
  parentVals[1] = 0;
  parentVals[2] = 1;
  mdcpt.setParentValues(parentVals);

  it = mdcpt.first();
  do {
    printf("Prob of %d is %f\n",
	   it.val,it.probVal.unlog());
  } while (mdcpt.next(it));


}


#endif
