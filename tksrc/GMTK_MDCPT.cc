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
#include "GMTK_RandomVariable.h"
#include "GMTK_DiscreteRandomVariable.h"
#include "GMTK_GMParms.h"

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
  : CPT(di_MDCPT)
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
  CPT::setNumParents(_nParents);

  // assume that the basic stuff is no longer allocated.
  bitmask &= ~bm_basicAllocated;
  cumulativeCardinalities.resize(_numParents);
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
  CPT::setNumCardinality(var,card);
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
  for (unsigned i=0;i<_numParents;i++) {
    numValues *= cardinalities[i];
  }
  numValues *= card();

  mdcpt.resize(numValues);

  if (_numParents > 0) {
    cumulativeCardinalities[_numParents-1] = card();
    for (int i=_numParents-2; i>=0; i--) {
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

  NamedObject::read(is);

  is.read(_numParents,"MDCPT::read numParents");

  if (_numParents < 0) 
    error("ERROR: reading file '%s', MDCPT '%s' trying to use negative (%d) num parents.",is.fileName(),name().c_str(),_numParents);
  if (_numParents >= warningNumParents)
    warning("WARNING: creating MDCPT with %d parents in file '%s'",_numParents,
	    is.fileName());

  cardinalities.resize(_numParents);
  cumulativeCardinalities.resize(_numParents);

  // read the parent cardinalities
  int numValues = 1;
  for (unsigned i=0;i<_numParents;i++) {
    is.read(cardinalities[i],"MDCPT::read cardinality");
    if (cardinalities[i] <= 0)
      error("ERROR: reading file '%s', MDCPT '%s' trying to use 0 or negative (%d) cardinality table, position %d.",
	    is.fileName(),name().c_str(),cardinalities[i],i);
    numValues *= cardinalities[i];
  }

  // read the self cardinalities
  is.read(_card,"MDCPT::read cardinality");
  if (_card <= 0)
    error("ERROR: reading file '%s', MDCPT '%s' trying to use 0 or negative (%d) cardinality table, position %d.",
	  is.fileName(),name().c_str(),_card,_numParents);
  numValues *= _card;

  // cumulativeCardinalities gets set to the
  // reverse cumulative cardinalities of the random
  // variables.
  if (_numParents > 0) {
    cumulativeCardinalities[_numParents-1] = card();
    for (int i=_numParents-2;i>=0;i--) {
      cumulativeCardinalities[i] = 
	cumulativeCardinalities[i+1]*cardinalities[i+1];
    }
  }

  // Finally read in the probability values (stored as doubles).
  // NOTE: We could check that things sum to approximately 1.0 here
  // (if we didn't use a large 1D loop).
  mdcpt.resize(numValues);
  for (int i=0;i<numValues;i++) {
    double val;
    is.readDouble(val,"MDCPT::read, reading value");
    if (val < 0 || val > 1)
      error("ERROR: reading file '%s', MDCPT '%s' has invalid probability value (%e), table entry number %d",
	    is.fileName(),
	    name().c_str(),
	    val,
	    i);
    mdcpt[i] = val;
  }
  setBasicAllocatedBit();
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
  NamedObject::write(os);
  os.nl();
  os.write(_numParents,"MDCPT::write numParents"); 
  os.writeComment("number parents");os.nl();
  for (unsigned i=0;i<_numParents;i++) {
    os.write(cardinalities[i],"MDCPT::write cardinality");
  }
  os.write(card(),"MDCPT::write cardinality");
  os.writeComment("cardinalities");
  os.nl();

  // Finally write in the probability values (stored as doubles).
  // NOTE: We could check that things sum to approximately 1 here, if
  // we didn't use a large 1D loop. 
  int childCard = card();
  for (int i=0;i<mdcpt.len();i++) {
    os.writeDouble(mdcpt[i].unlog(),"MDCPT::write, writing value");
    childCard --;
    if (childCard == 0) {
      os.nl();
      childCard = card();
    }
  }
}


////////////////////////////////////////////////////////////////////
//        Probability Evaluation
////////////////////////////////////////////////////////////////////



/*-
 *-----------------------------------------------------------------------
 * MDCPT::becomeAwareOfParentValues()
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
MDCPT::becomeAwareOfParentValues( vector<int>& parentValues,
				  vector<int>& cards)
{

  assert ( parentValues.size() == _numParents );
  assert ( bitmask & bm_basicAllocated );
  
  int offset = 0;
  for (unsigned i = 0; i < _numParents; i++) {
    if (parentValues[i] < 0 || parentValues[i] >= cardinalities[i]) 
      error("MDCPT:becomeAwareOfParentValues: Invalid parent value for parent %d, parentValue = %d but card = %d\n",i,parentValues[i],cardinalities[i]);
    offset += parentValues[i]*cumulativeCardinalities[i];
  }
  mdcpt_ptr = mdcpt.ptr + offset;

}



/*-
 *-----------------------------------------------------------------------
 * MDCPT::becomeAwareOfParentValues()
 *      Like above, but uses the explicit array of parents.
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
MDCPT::becomeAwareOfParentValues( vector< RandomVariable * >& parents)
{

  assert ( parents.size() == _numParents );
  assert ( bitmask & bm_basicAllocated );
  
  int offset = 0;
  for (unsigned i = 0; i < _numParents; i++) {
    if ( parents[i]->val < 0 || parents[i]->val >= cardinalities[i])
      error("MDCPT:becomeAwareOfParentValues: Invalid parent value for parent %d, parentValue = %d but card = %d\n",i,parents[i]->val,cardinalities[i]);
    offset += parents[i]->val*cumulativeCardinalities[i];
  }
  mdcpt_ptr = mdcpt.ptr + offset;
}



////////////////////////////////////////////////////////////////////
//        Misc Support
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * MDCPT::randomSample()
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
MDCPT::randomSample()
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
 * MDCPT::normalize()
 *      Re-normalize all the distributions
 *  
 * Results:
 *      No results.
 *
 * Side Effects:
 *      Changes the values of all tables.
 *
 *-----------------------------------------------------------------------
 */
void
MDCPT::normalize()
{
  assert ( bitmask & bm_basicAllocated );

  // Use the inherent structure of the multi-D array
  // so to loop over the final distributions on the child.

  const int child_card = card();
  const int num_parent_assignments = mdcpt.len()/child_card;
  logpr *loc_ptr = mdcpt.ptr;
  for (int parent_assignment =0; 
       parent_assignment < num_parent_assignments; 
       parent_assignment ++) {
    logpr sum = 0.0;
    logpr *tmp_loc_ptr = loc_ptr;
    for (int i=0;i<child_card;i++) {
      sum += *tmp_loc_ptr++;
    }
    tmp_loc_ptr = loc_ptr;
    for (int i=0;i<child_card;i++) {
      *tmp_loc_ptr /= sum;
      tmp_loc_ptr++;
    }
    loc_ptr += child_card;
  }
}


/*-
 *-----------------------------------------------------------------------
 * MDCPT::makeRandom()
 *      Assign random but valid values to the CPT distribution.
 *  
 * Results:
 *      No results.
 *
 * Side Effects:
 *      Changes the values of all tables.
 *
 *-----------------------------------------------------------------------
 */
void
MDCPT::makeRandom()
{
  assert ( bitmask & bm_basicAllocated );

  // Use the inherent structure of the multi-D array
  // so to loop over the final distributions on the child.

  const int child_card = card();
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




/*-
 *-----------------------------------------------------------------------
 * MDCPT::makeUnifrom()
 *      Have distribution be entirely uniform.
 *  
 * Results:
 *      No results.
 *
 * Side Effects:
 *      Changes the values of all tables.
 *
 *-----------------------------------------------------------------------
 */
void
MDCPT::makeUniform()
{
  assert ( bitmask & bm_basicAllocated );

  // Use the inherent structure of the multi-D array
  // so to loop over the final distributions on the child.

  const int child_card = card();
  double u_val = 1.0/(double)child_card;
  const int num_parent_assignments = mdcpt.len()/child_card;
  logpr *loc_ptr = mdcpt.ptr;
  for (int parent_assignment =0; 
       parent_assignment < num_parent_assignments; 
       parent_assignment ++) {
    logpr *tmp_loc_ptr = loc_ptr;
    for (int i=0;i<child_card;i++) {
      *tmp_loc_ptr ++ = u_val;
    }
    loc_ptr += child_card;
  }
}


////////////////////////////////////////////////////////////////////
//        EM Routines
////////////////////////////////////////////////////////////////////

void
MDCPT::emStartIteration()
{
  if (!GM_Parms.amTrainingMDCPTs())
    return;

  if(emOnGoingBitIsSet())
    return; // already done

  if (!emEmAllocatedBitIsSet()) {
    nextMdcpt.resize(mdcpt.len());
    emSetEmAllocatedBit();
  }
  // EM iteration is now going.
  emSetOnGoingBit();
  emSetSwappableBit();

  accumulatedProbability = 0.0;  
  // zero the accumulators
  // or if we want to add priors here, we can do that at this point.
  for (int i=0;i<nextMdcpt.len();i++) {
    nextMdcpt[i].set_to_zero();
  }
}


void
MDCPT::emIncrement(logpr prob,RandomVariable* rv)
{
  if (!GM_Parms.amTrainingMDCPTs())
    return;

  emStartIteration();


  // this is an MDCPT, so rv must be discrete.a
  assert ( rv -> discrete );

  DiscreteRandomVariable* drv = (DiscreteRandomVariable*)rv;
  // make sure, by checking that drv's curCPT points to this.
  assert ( drv -> curCPT == this );

  // 
  // TODO: This needs to be factored out of the inner most
  // loop!
  becomeAwareOfParentValues(*(drv->curConditionalParents));

  // Grab the current offset ...
  int offset = mdcpt_ptr-mdcpt.ptr;

  // ... and use it for the next cpt 
  *(nextMdcpt.ptr + offset + drv->val) += prob;

  accumulatedProbability += prob;
}


void
MDCPT::emEndIteration()
{
  if (!GM_Parms.amTrainingMDCPTs())
    return;

  if (!emOnGoingBitIsSet())
    return;

  if (accumulatedProbability.zero()) {
    warning("WARNING: MDCPT named '%s' did not receive any accumulated probability in EM iteration",name().c_str());
  }

  // now normalize the next ones

  const int child_card = card();
  const int num_parent_assignments = mdcpt.len()/child_card;
  logpr *loc_ptr = nextMdcpt.ptr;
  for (int parent_assignment =0; 
       parent_assignment < num_parent_assignments; 
       parent_assignment ++) {
    logpr sum = 0.0;
    logpr *tmp_loc_ptr = loc_ptr;
    for (int i=0;i<child_card;i++) {
      sum += *tmp_loc_ptr++;
    }

    if (sum == 0.0) {
      error("Ending EM iteration but a row of CPT %s has zero probability of occurance\n",
	  _name.c_str());      
    }
    tmp_loc_ptr = loc_ptr;
    for (int i=0;i<child_card;i++) {
      *tmp_loc_ptr /= sum;
      tmp_loc_ptr++;
    }
    loc_ptr += child_card;
  }

  // stop EM
  emClearOnGoingBit();

}

void
MDCPT::emSwapCurAndNew()
{
  if (!GM_Parms.amTrainingMDCPTs())
    return;

  if (!emSwappableBitIsSet())
    return;

  mdcpt.swapPtrs(nextMdcpt);
  emClearSwappableBit();
}


void
MDCPT::emStoreAccumulators(oDataStreamFile& ofile)
{
  assert ( basicAllocatedBitIsSet() );
  if ( !emEmAllocatedBitIsSet() ) {
    warning("WARNING: storing zero accumulators for MDCPT '%s'\n",
	    name().c_str());
    emStoreZeroAccumulators(ofile);
    return;
  }
  EMable::emStoreAccumulators(ofile);
  for (int i=0;i<nextMdcpt.len();i++) {
    ofile.write(nextMdcpt[i].val(),"mdcw");
  }
}

void
MDCPT::emStoreZeroAccumulators(oDataStreamFile& ofile)
{
  assert ( basicAllocatedBitIsSet() );
  EMable::emStoreZeroAccumulators(ofile);
  for (int i=0;i<mdcpt.len();i++) {
    const logpr p; 
    ofile.write(p.val(),"mdcw");
  }
}

void
MDCPT::emLoadAccumulators(iDataStreamFile& ifile)
{
  assert (basicAllocatedBitIsSet());
  assert (emEmAllocatedBitIsSet());
  EMable::emLoadAccumulators(ifile);
  for (int i=0;i<nextMdcpt.len();i++) {
    ifile.read(nextMdcpt[i].valref(),"mdcr");
  }
}


void
MDCPT::emAccumulateAccumulators(iDataStreamFile& ifile)
{
  assert ( basicAllocatedBitIsSet() );
  assert ( emEmAllocatedBitIsSet() );
  EMable::emAccumulateAccumulators(ifile);
  for (int i=0;i<nextMdcpt.len();i++) {
    logpr tmp;
    ifile.read(tmp.valref(),"mdcr");
    nextMdcpt[i] += tmp;
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

  printf("mdcpt name = %s\n",mdcpt.name().c_str());

  // now print out some probabilities.
  vector < int > parentVals;
  parentVals.resize(3);

  parentVals[0] = 0;
  parentVals[1] = 0;
  parentVals[2] = 0;

  printf("parentVals:");
  for (unsigned i=0;i<3;i++) 
    printf("%d ",parentVals[i]);
  printf("\n");
  for (int i =0; i<3;i++) {
    printf("Prob(%d) Given cur Par = %f\n",
	   i,mdcpt.probGivenParents(parentVals,i).unlog());
  }

  parentVals[0] = 0;
  parentVals[1] = 0;
  parentVals[2] = 1;

  printf("parentVals:");
  for (unsigned i=0;i<3;i++) 
    printf("%d ",parentVals[i]);
  printf("\n");
  for (int i =0; i<3;i++) {
    printf("Prob(%d) Given cur Par = %f\n",
	   i,mdcpt.probGivenParents(parentVals,i).unlog());
  }

  parentVals[0] = 0;
  parentVals[1] = 1;
  parentVals[2] = 1;

  printf("parentVals:");
  for (unsigned i=0;i<3;i++) 
    printf("%d ",parentVals[i]);
  printf("\n");
  for (int i =0; i<3;i++) {
    printf("Prob(%d) Given cur Par = %f\n",
	   i,mdcpt.probGivenParents(parentVals,i).unlog());
  }

  parentVals[0] = 1;
  parentVals[1] = 2;
  parentVals[2] = 1;

  printf("parentVals:");
  for (unsigned i=0;i<3;i++) 
    printf("%d ",parentVals[i]);
  printf("\n");
  for (int i =0; i<3;i++) {
    printf("Prob(%d) Given cur Par = %f\n",
	   i,mdcpt.probGivenParents(parentVals,i).unlog());
  }

  // Now iterate over valid values.
  MDCPT::iterator it = mdcpt.begin();
  do {
    printf("Prob of %d is %f\n",
	   it.val(),it.probVal.unlog());
    it++;
  } while (it != mdcpt.end());


  parentVals[0] = 0;
  parentVals[1] = 0;
  parentVals[2] = 1;
  printf("parentVals:");
  for (unsigned i=0;i<3;i++) 
    printf("%d ",parentVals[i]);
  printf("\n");
  mdcpt.becomeAwareOfParentValues(parentVals);

  it = mdcpt.begin();
  do {
    printf("Prob of %d is %f\n",
	   it.val(),it.probVal.unlog());
    it++;
  } while (it != mdcpt.end());

  mdcpt.makeRandom();
  printf("After randomization\n");
  it = mdcpt.begin();
  do {
    printf("Prob of %d is %f\n",
	   it.val(),it.probVal.unlog());
    it++;
  } while (it != mdcpt.end());


}


#endif
