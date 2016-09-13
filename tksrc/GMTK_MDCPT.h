/*-
 * GMTK_MDCPT.h
 *
 *  Written by Jeff Bilmes <bilmes@ee.washington.edu>
 * 
 * Copyright (C) 2001 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */


#ifndef GMTK_MDCPT_H
#define GMTK_MDCPT_H

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include "fileParser.h"
#include "logp.h"

#include "sArray.h"

#include "GMTK_CPT.h"
#include "GMTK_EMable.h"
#include "GMTK_DiscRV.h"
#include "GMTK_NamedObject.h"
#include "GMTK_DirichletTable.h"
#include "GMTK_DirichletPrior.h"

class MDCPT : public CPT, public DirichletPrior  {

  //////////////////////////////////
  // The acutal cpt. This is the table for
  // Pr( variable at mdcpt.len()-1 |  variables from 0 to mdcpt.len()-1 )
  sArray < logpr > mdcpt;

  //////////////////////////////////
  // Support for computing probabilities as below.
  // This indices the section of the table mdcpt above
  // for a particular set of parent values. 
  logpr* mdcpt_ptr; 
  // The accumulative cardinalities to help index into
  // the table above.
  // TODO: make this an sArray.
  vector <int> cumulativeCardinalities;

  //////////////////////////////////
  // Data structures support for EM
  sArray < logpr > nextMdcpt;

  // the overall expected occurence of this CPT
  logpr accumulator;

protected:

  // special constructor for subclasses who want to 
  // use a different type.
  MDCPT(DiscreteImplementaton _cptType) : CPT(_cptType) {};

public:

  ///////////////////////////////////////////////////////////  
  // General constructor
  MDCPT();
  ~MDCPT() { }

    sArray<logpr>& getMdcpt() {return mdcpt;}
    sArray<logpr>& getNextMdcpt() {return nextMdcpt;}

  ///////////////////////////////////////////////////////////    
  // Semi-constructors: useful for debugging.
  // See parent class for further documention.
  void setNumParents(const unsigned _nParents);
  void setNumCardinality(const unsigned var, const int card);
  void allocateBasicInternalStructures();


  ///////////////////////////////////////////////////////////  
  // Probability evaluation, compute Pr( child | parents ), and
  // iterator support. See GMTK_CPT.h for documentation.
  void becomeAwareOfParentValues( vector <RV*>& parents,
				  const RV* rv);
  void begin(iterator& it,DiscRV* drv,logpr& p) {
    assert ( bitmask & bm_basicAllocated );
    it.setCPT(this);
    it.internalStatePtr = (void*)mdcpt_ptr;
    it.drv = drv;

    REGISTER DiscRVType value = 0;
    while (mdcpt_ptr[value].essentially_zero()) {
      value++;
      // We keep the following check as we must have that at least one
      // entry is non-zero.  The read code of the MDCPT should ensure
      // this as sure all parameter update procedures, as long as
      // normalizationThreshold is not set to large.
      if (value >= card()) {
	error("ERROR: DenseCPT '%s' used for RV '%s(%d)' has a row with all zeros. Program Exiting ...\n",
	      name().c_str(),drv->name().c_str(),drv->frame());
      }
    }
    p = mdcpt_ptr[value];    
    drv->val = value;
  }
  void becomeAwareOfParentValuesAndIterBegin(vector < RV *>& parents, 
					     iterator &it, 
					     DiscRV* drv,
					     logpr& p);
  logpr probGivenParents(vector < RV* >& parents,
			 DiscRV* drv) {
    assert ( bitmask & bm_basicAllocated );
    MDCPT::becomeAwareOfParentValues(parents,drv);
    REGISTER DiscRVType val = drv->val;
    assert ( val < card() );
    return *(mdcpt_ptr + val);
  }
  bool next(iterator &it,logpr& p) {
    REGISTER logpr* const loc_mdcpt_ptr = (logpr*)it.internalStatePtr;
    REGISTER DiscRVType value = it.drv->val;
    do {
      value++;
      // don't increment past the last value.
      if (value == card()) {
	return false;
      }
    } while (loc_mdcpt_ptr[value].essentially_zero());
    p = loc_mdcpt_ptr[value];
    it.drv->val = value;    
    return true;
  }

  ///////////////////////////////////
  int randomSample(DiscRV* drv);

  ///////////////////////////////////
  unsigned totalNumberParameters() { return mdcpt.len(); }

  ///////////////////////////////////////////////////////////  
  // Re-normalize the output distributions
  void normalize();
  // set all values to random values.
  void makeRandom();
  // set all values to uniform values.
  void makeUniform();

  //////////////////////////////////////////////
  // read/write basic parameters
  void read(iDataStreamFile& is);
  void write(oDataStreamFile& os);

  //////////////////////////////////
  // Public interface support for EM
  //////////////////////////////////
  void emStartIteration();
  void emIncrement(logpr p,vector <RV*>& parents,RV*rv);
  void emEndIteration();
  void emSwapCurAndNew();

  // parallel training
  void emStoreObjectsAccumulators(oDataStreamFile& ofile,
				  bool writeLogVals = true,
				  bool writeZeros = false);
  void emLoadObjectsDummyAccumulators(iDataStreamFile& ifile);
  void emZeroOutObjectsAccumulators();
  void emLoadObjectsAccumulators(iDataStreamFile& ifile);
  void emAccumulateObjectsAccumulators(iDataStreamFile& ifile);
  const string typeName() { return "DenseCPT"; }
  //////////////////////////////////


};



#endif // defined MDCPT
