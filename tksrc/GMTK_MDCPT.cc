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
#include "GMTK_DiscRV.h"
#include "GMTK_GMParms.h"


VCID("$Header$")

/*
 * This routine copies the magnitude of x, the sign of y, and returns the result, i.e.,
 *  copysign(x,y) = fabs(x)*sign(y)
 */      
extern "C" double copysign(double x, double y) __THROW;


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

  if (numValues < 0) {
    fprintf(stderr,"ERROR: DenseCPT '%s' with %d parents, CPT table size too large.  Combined cardinalities of parents & child too large: ",name().c_str(),_numParents);
    fprintf(stderr,"child card = %d, parent cardinalities are:",card());
    for (unsigned i=0;i<_numParents;i++) {
      fprintf(stderr," %d",cardinalities[i]);
    }
    error(" exiting ...");
  }

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

  is.read(_numParents,"Can't read DenseCPT numParents");

  if (_numParents < 0) 
    error("ERROR: reading file '%s' line %d, DenseCPT '%s' trying to use negative (%d) num parents.",
	  is.fileName(),is.lineNo(),name().c_str(),_numParents);
  if (_numParents >= warningNumParents)
    warning("WARNING: creating DenseCPT with %d parents in file '%s' line %d",
	    _numParents,is.fileName(),is.lineNo());

  cardinalities.resize(_numParents);
  cumulativeCardinalities.resize(_numParents);

  // read the parent cardinalities
  int numValues = 1;
  for (unsigned i=0;i<_numParents;i++) {
    is.read(cardinalities[i],"Can't read DenseCPT parent cardinality");
    if (cardinalities[i] <= 0)
      error("ERROR: reading file '%s' line %d, DenseCPT '%s' trying to use 0 or negative (%d) cardinality table, position %d.",
	    is.fileName(),is.lineNo(),name().c_str(),cardinalities[i],i);
    numValues *= cardinalities[i];
  }

  // read the self cardinalities
  is.read(_card,"Can't read DenseCPT self cardinality");
  if (_card <= 0)
    error("ERROR: reading file '%s' line %d, DenseCPT '%s' trying to use 0 or negative (%d) cardinality table, position %d.",
	  is.fileName(),is.lineNo(),name().c_str(),_card,_numParents);
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

  // read optional smoothing parameters. We support either
  //   1) constant dirichlet priors, where a constant alpha
  //      is given right here, where we have a Dirichlet 
  //      with hyperparameter alpha >= 0 which is constant for
  //      all rv values. The syntax for this is:
  //         
  //         DirichletConst <alpha>
  // 
  //      and where the accumulators now have E[counts] + alpha rather
  //      than just alpha. Note that \alpha is given as a real fractional
  //      count value, so it is not a log prob.
  //   2) specify a counts object, where we have a full set of 'numValues'
  //      counts for all values of the random variable for all possible
  //      parent values. This is much more general than the above, as
  //      the counts object can specify a different Dirichlet hyperparameter
  //      for each value of the RV. Syntax for this is:
  //   
  //        DirichletTable table-object
  //
  //      where table-object is a previously defined Dirichlet Table object in the
  //      master file which is compatible with this table.
  
  if (is.readIfMatch(DirichletConstStr,"DenseCPT double value or Dirichlet const spec")) {
    // so we should have a single constant alpha value next.
    is.read(dirichletAlpha,"Can't read DenseCPT Dirichlet hyperparameter");
    smoothingType = DirichletConstVal;
  } else if (is.readIfMatch(DirichletTableStr,"DenseCPT double value or Dirichlet table spec")) {
    // so we should have a pointer to a previously existing count table.
    string dirichletTableName;
    is.read(dirichletTableName);
    if (GM_Parms.dirichletTabsMap.find(dirichletTableName) == GM_Parms.dirichletTabsMap.end()) {
	error("ERROR: reading file '%s' line %d, DenseCPT '%s' specified Dirichlet Table (%s) that does not exist",
	      is.fileName(),is.lineNo(),
	      name().c_str(),
	      dirichletTableName.c_str());

    }
    dirichletTable = GM_Parms.dirichletTabs[GM_Parms.dirichletTabsMap[dirichletTableName]];
    smoothingType = DirichletTableVal;
    // next check that the table matches the CPT.
    if (_numParents + 1 != dirichletTable->numDimensions()) {
	error("ERROR: reading file '%s' line %d, DenseCPT '%s' has %d parents (a %d-D table), but Dirichlet Table '%s' has dimensionality %d",
	      is.fileName(),is.lineNo(),
	      name().c_str(),
	      _numParents,_numParents+1,
	      dirichletTable->name().c_str(),
	      dirichletTable->numDimensions());
    }

    // check that parents match
    for (unsigned i=0;i<_numParents;i++) {    
      if (cardinalities[i] != dirichletTable->dimension(i)) {
	error("ERROR: reading file '%s' line %d, in DenseCPT '%s', %d'th parent has cardinality %d, but Dirichlet Table '%s' has its %d'th dimension of size %d",
	      is.fileName(),is.lineNo(),
	      name().c_str(),
	      i,cardinalities[i],
	      dirichletTable->name().c_str(),
	      i,dirichletTable->dimension(i));
      }
    }

    // check self cardinality
    if (_card != dirichletTable->lastDimension()) {
      error("ERROR: reading file '%s' line %d, in DenseCPT '%s', child has cardinality %d, but Dirichlet Table '%s' has its last dimension of size %d",
	    is.fileName(),is.lineNo(),
	    name().c_str(),
	    _card,
	    dirichletTable->name().c_str(),
	    dirichletTable->lastDimension());
    }

    // everything matches up, but include last sanity check
    assert ( dirichletTable->tableSize() == (unsigned)numValues );

  }

  // Finally read in the probability values (stored as doubles).
  mdcpt.resize(numValues);
  logpr child_sum;
  child_sum.set_to_zero();
  int row=0;;
  // be more forgiving as cardinality increases
  const double threshold = _card*normalizationThreshold;
  for (int i=0;i<numValues;) {

    double val;  // sign bit below needs to be changed if we change this type.
    is.readDouble(val,"Can't read DenseCPT double value");


    // we support reading in both regular probability values
    // (in the range [+0,1] inclusive) and log probability 
    // values (in the range (-infty,-0] inclusive. These
    // ranges give distinct values for probabilties, except for
    // the value 0 which can either be real probability zero (impossible
    // event) or it could be log(1) = 0 (the certain event). Since
    // the IEEE FP standard supports both +0 and -0, and since the
    // ASCII read routines preserve ASCII string '-0.0' to be negative zero,
    // we consider -0.0 as log(1) , and +0.0 as real zero.
    if (val > 1)
      error("ERROR: reading file '%s' line %d, DenseCPT '%s' has invalid probability value (%e), table entry number %d",
	    is.fileName(),is.lineNo(),
	    name().c_str(),
	    val,
	    i);
    if (val > 0) { 
      // regular probability
      mdcpt[i] = val;
    } else if (val < 0) {
      // log base e probability
      mdcpt[i].setFromLogP(val);
    } else {
      // is zero, so need to check sign bit for
      // either -0 (log(1)) or +0 (true zero prob)
      if (copysign(1.0,val)==1.0) {
	// regular zero probability
	mdcpt[i].set_to_zero();
      } else {
	// val == -0, so set to log(1)
	mdcpt[i].set_to_one();	
      }
    }
    child_sum += mdcpt[i];

    i++;
    if (i % _card == 0 && (normalizationThreshold != 0)) {
      // check that child sum is approximately one if (normalizationThreshold != 0)
      // which otherwise would turn it off.
      double abs_diff = fabs(child_sum.unlog() - 1.0);
      if (abs_diff > threshold) 
	error("ERROR: reading file '%s' line %d, row %d of DenseCPT '%s' has probabilities that sum to %e but should sum to unity, absolute difference = %e, current normalization threshold = %f.",
	      is.fileName(),is.lineNo(),
	      row,
	      name().c_str(),
	      child_sum.unlog(),
	      abs_diff,
	      normalizationThreshold);
      // reset
      child_sum.set_to_zero();
      row++;
    }
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
  assert ( basicAllocatedBitIsSet() );
  NamedObject::write(os);
  os.nl();
  os.write(_numParents,"DenseCPT::write numParents"); 
  os.writeComment("number parents");os.nl();
  for (unsigned i=0;i<_numParents;i++) {
    os.write(cardinalities[i],"DenseCPT::write cardinality");
  }
  os.write(card(),"DenseCPT::write cardinality");
  os.writeComment("cardinalities");
  os.nl();

  if (smoothingType == DirichletConstVal) {
    os.write(DirichletConstStr);
    os.write(dirichletAlpha);
    os.nl();
  } else if (smoothingType == DirichletTableVal) {
    os.write(DirichletTableStr);
    os.write(dirichletTable->name());
    os.nl();
  }

  // Finally write in the probability values (stored as doubles).
  normalize();
  int childCard = card();
  for (int i=0;i<mdcpt.len();i++) {
    os.writeDouble(mdcpt[i].unlog(),"DenseCPT::write, writing value");
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
 *
 * Preconditions:
 *      parents must be an array of parents, and actually must
 *      really be an array of DiscRVs. This condition is not
 *      checked, and if it is not true, random results will occur.
 *
 * Postconditions:
 *      we are aware of the parent values.
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
MDCPT::becomeAwareOfParentValues( vector< RV * >& parents,
				  const RV* rv)
{

  assert ( basicAllocatedBitIsSet() );
  assert ( parents.size() == _numParents );
  
  int offset = 0;
  for (unsigned i = 0; i < _numParents; i++) {
    if (((DiscRV*)parents[i])->val >= cardinalities[i])
      error("ERROR: DenseCPT %s used by RV %s(%d), invalid parent value for parent %s(%d) (parent number %d), parentValue = %d but parent RV cardinality = %d\n",
	    name().c_str(),
	    rv->name().c_str(),rv->frame(),
	    parents[i]->name().c_str(),parents[i]->frame(),
	    i,
	    ((DiscRV*)parents[i])->val,cardinalities[i]);
    offset += ((DiscRV*)parents[i])->val*cumulativeCardinalities[i];
  }
  mdcpt_ptr = mdcpt.ptr + offset;
}
void
MDCPT::becomeAwareOfParentValuesAndIterBegin( vector< RV * >& parents,
					      iterator & it,
					      DiscRV* drv,
					      logpr& p)
{

  assert ( basicAllocatedBitIsSet() );
  assert ( parents.size() == _numParents );
  
  int offset = 0;
  for (unsigned i = 0; i < _numParents; i++) {
    if ( ((DiscRV*)parents[i])->val >= cardinalities[i])
      error("ERROR:becomeAwareOfParentValues. DenseCPT %s, invalid parent value for parent %s(%d) (parent number %d), parentValue = %d but RV cardinality = %d\n",
	    name().c_str(),
	    parents[i]->name().c_str(),parents[i]->frame(),
	    i,
	    ((DiscRV*)parents[i])->val,cardinalities[i]);
    offset += ((DiscRV*)parents[i])->val*cumulativeCardinalities[i];
  }
  register logpr* const mdcpt_ptr = mdcpt.ptr + offset;

  it.setCPT(this);
  it.internalStatePtr = (void*)mdcpt_ptr;
  it.drv = drv;

  register DiscRVType value = 0;
  while (mdcpt_ptr[value].essentially_zero()) {
    value++;
    // We keep the following check as we must have that at least one
    // entry is non-zero.  The read code of the MDCPT should ensure
    // this as sure all parameter update procedures, as long as
    // normalizationThreshold is not set to large.
    if (value >= card()) {
      fprintf(stderr,"ERROR: DenseCPT '%s' used for RV '%s(%d)' has a row with all zeros. Parents and values are: ",
	      name().c_str(),drv->name().c_str(),drv->frame());
      printRVSetAndValues(stderr,parents);
      error("Program Exiting...\n"); 
    }
  }
  p = mdcpt_ptr[value];    
  drv->val = value;
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
MDCPT::randomSample(DiscRV*drv)
{
  assert ( basicAllocatedBitIsSet() );
  
  logpr* mdcpt_ptr_p = mdcpt_ptr;
  logpr* mdcpt_ptr_endp = mdcpt_ptr + card();
  logpr uniform = rnd.drand48();
  logpr sum = 0.0;
  do {
    sum += (*mdcpt_ptr_p);
    if (uniform <= sum)
      break;
    mdcpt_ptr_p++;
  } while (mdcpt_ptr_p != mdcpt_ptr_endp);
  
  return (drv->val = (mdcpt_ptr_p - mdcpt_ptr));
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
  assert ( basicAllocatedBitIsSet() );

  assert ( card() > 0 );

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
    int i=0; do {
      sum += *tmp_loc_ptr++;
    } while (++i < child_card);

    tmp_loc_ptr = loc_ptr;
    i=0; do {
      *tmp_loc_ptr /= sum;
      (*tmp_loc_ptr).floor();
      tmp_loc_ptr++;   
    } while (++i < child_card);

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
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;

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
      (*tmp_loc_ptr).floor();
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
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;


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
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
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
  if (smoothingType == NoneVal || !useDirichletPriors) {
    for (int i=0;i<nextMdcpt.len();i++) {
      nextMdcpt[i].set_to_zero();
    }
  } else if  (smoothingType == DirichletConstVal) {
    logpr alpha(dirichletAlpha);
    for (int i=0;i<nextMdcpt.len();i++) {
      nextMdcpt[i] = alpha;
    }
  } else {
    // dirichlet table
    for (int i=0;i<nextMdcpt.len();i++) {
      nextMdcpt[i] = dirichletTable->tableValue(i);
    }
  }

}


void
MDCPT::emIncrement(logpr prob,vector <RV*>& parents,RV* rv)
{
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;

  if (!emOnGoingBitIsSet())
    emStartIteration();

  // this is an MDCPT, so rv must be discrete.a
  assert ( rv -> discrete() );

  DiscRV* drv = RV2DRV(rv);
  // make sure, by checking that drv's curCPT points to this.
  assert ( drv -> curCPT == this );

  MDCPT::becomeAwareOfParentValues(parents,rv);

  // Grab the current offset ...
  int offset = mdcpt_ptr-mdcpt.ptr;

  // ... and use it for the next cpt 
  *(nextMdcpt.ptr + offset + drv->val) += prob;

  accumulatedProbability += prob;
}


void
MDCPT::emEndIteration()
{
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;

  if (!emOnGoingBitIsSet())
    return;

  accumulatedProbability.floor();
  if (accumulatedProbability < minDiscAccumulatedProbability()) {
    warning("WARNING: DenseCPT named '%s' received only %e accumulated probability in EM iteration. Using previous iteraton values.",name().c_str(),accumulatedProbability.val());
    for (int i=0;i<nextMdcpt.len();i++) {
      nextMdcpt[i] = mdcpt[i];
    }
  } else {

    // now normalize the next ones
    const int child_card = card();
    const int num_parent_assignments = mdcpt.len()/child_card;
    logpr *loc_ptr = nextMdcpt.ptr;
    int num_rows_with_zero_counts = 0;
    for (int parent_assignment =0; 
	 parent_assignment < num_parent_assignments; 
	 parent_assignment ++) {
      logpr sum = 0.0;
      logpr *tmp_loc_ptr = loc_ptr;
      for (int i=0;i<child_card;i++) {
	sum += *tmp_loc_ptr++;
      }

      sum.floor();
      if (sum == 0.0) {
	num_rows_with_zero_counts ++;
	logpr *mdcpt_p = mdcpt.ptr + (loc_ptr - nextMdcpt.ptr);
	tmp_loc_ptr = loc_ptr;
	for (int i=0;i<child_card;i++) {
	  *tmp_loc_ptr++ = *mdcpt_p++; 
	}
      } else {
	tmp_loc_ptr = loc_ptr;
	for (int i=0;i<child_card;i++) {
	  *tmp_loc_ptr /= sum;
	  (*tmp_loc_ptr).floor();
	  tmp_loc_ptr++;
	}
      }

      loc_ptr += child_card;
    }

    if (num_rows_with_zero_counts > 0) 
      warning("WARNING: Ending EM iteration but %d rows of DenseCPT '%s' had zero counts. Using previous values for those rows.\n",
	      num_rows_with_zero_counts,
	      _name.c_str());

  }
  // stop EM
  emClearOnGoingBit();

}

void
MDCPT::emSwapCurAndNew()
{
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;

  if (!emSwappableBitIsSet())
    return;

  mdcpt.swapPtrs(nextMdcpt);
  emClearSwappableBit();
}



/*-
 *-----------------------------------------------------------------------
 *
 * Accumulator loading/storing routines for parallel training support.
 * These routines are virtual, and are called from the EMable object
 * which has the code containing the logic and checks about which one
 * to call when.
 *
 *----------------------------------------------------------------------- 
*/



void
MDCPT::emStoreObjectsAccumulators(oDataStreamFile& ofile)
{
  for (int i=0;i<nextMdcpt.len();i++) {
    ofile.write(nextMdcpt[i].val(),"DenseCPT store accums");
  }
}


void
MDCPT::emLoadObjectsDummyAccumulators(iDataStreamFile& ifile)
{
  logpr tmp;
  for (int i=0;i<mdcpt.len();i++) {
    ifile.read(tmp.valref(),"DenseCPT load accums");
  }
}


void
MDCPT::emZeroOutObjectsAccumulators()
{
  for (int i=0;i<nextMdcpt.len();i++) {
    nextMdcpt[i].set_to_zero();
  }
}

void
MDCPT::emLoadObjectsAccumulators(iDataStreamFile& ifile)
{
  for (int i=0;i<nextMdcpt.len();i++) {
    ifile.read(nextMdcpt[i].valref(),"DenseCPT load accums");
  }
}


void
MDCPT::emAccumulateObjectsAccumulators(iDataStreamFile& ifile)
{
  for (int i=0;i<nextMdcpt.len();i++) {
    logpr tmp;
    ifile.read(tmp.valref(),"DenseCPT accumulate accums");
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
