/*-
 * GMTK_MTCPT.h
 *
 *  Written by Jeff Bilmes <bilmes@ee.washington.edu>
 * 
 *  $Header$
 * 
 * Copyright (c) 2001, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about
 * the suitability of this software for any purpose.  It is provided
 * "as is" without express or implied warranty.
 *
 */


#ifndef GMTK_MTCPT_H
#define GMTK_MTCPT_H

#include <vector>

#include "fileParser.h"
#include "logp.h"

#include "GMTK_RngDecisionTree.h"
#include "GMTK_DiscRV.h"
#include "GMTK_CPT.h"
#include "GMTK_Sparse1DPMF.h"

#include "GMTK_EMable.h"
#include "GMTK_GMParms.h"
#include "GMTK_NamedObject.h"


class MTCPT : public CPT  {

  //////////////////////////////////
  // Index into the world structure
  // of the decision tree
  unsigned dtIndex; 

  ///////////////////////////////////////
  // Direct pointer to the decision tree.
  RngDecisionTree* dt;

  ////////////////
  // the value that has prob one for current parent
  // values. I.e., P(child = _val| Parents = par_vals) = 1.
  unsigned _val;

public:

  ///////////////////////////////////////////////////////////  
  // General constructor
  MTCPT();
  ~MTCPT() { }

  ///////////////////////////////////////////////////////////    
  void setNumParents(const int _nParents);
  void setNumCardinality(const int var, const int card);
  void allocateBasicInternalStructures();
  ///////////////////////////////////////////////////////////    

  //////////////////////////////////
  bool iterable() { return dt->iterable(); } 

  //////////////////////////////////
  // various forms of probability calculation

  ///////////////////////////////////////////////////////////  
  // Probability evaluation, compute Pr( child | parents ), and
  // iterator support.  See GMTK_CPT.h for documentation.
  void becomeAwareOfParentValues( vector < RV *>& parents,
				  const RV* rv ) {
    _val = dt->query(parents,rv);
    if (_val >= card()) {
      warning("ERROR: Deterministic CPT '%s' of card %d used by RV %s(%d) querying DT '%s' received value %d",
	      name().c_str(),
	      card(),
	      rv->name().c_str(),rv->frame(),
	      dt->name().c_str(),
	      _val);
      fprintf(stderr,"Parents configuration :");
      printRVSetAndValues(stderr,parents);
      error("");
    }
  }
  void begin(iterator& it,DiscRV* drv,logpr& p) {
    assert ( bitmask & bm_basicAllocated );
    it.setCPT(this);
    it.drv = drv;
    p  = 1.0;
    drv->val = _val;
  }
  void becomeAwareOfParentValuesAndIterBegin( vector < RV* >& parents,
					      iterator & it,
					      DiscRV* drv,
					      logpr& p)
  {
    _val = dt->query(parents,drv);
    if (_val >= card()) {
      warning("ERROR: Deterministic CPT '%s' of card %d querying DT '%s' received value %d",
	      name().c_str(),
	      card(),
	      dt->name().c_str(),
	      _val);
      fprintf(stderr,"Parents configuration :");
      printRVSetAndValues(stderr,parents);
      error("");
    }
    assert ( bitmask & bm_basicAllocated );
    it.setCPT(this);
    it.drv = drv;
    p = 1.0;
    drv->val = _val;
  }
  logpr probGivenParents(vector < RV* >& parents,
			 DiscRV* drv) {
    assert ( bitmask & bm_basicAllocated );
    becomeAwareOfParentValues(parents,drv);
    register DiscRVType val = drv->val;
    assert ( val <= card() );
    if (val == _val)
      return 1.0;
    else
      return 0.0;
  }
  inline bool next(iterator &it,logpr& p) {
    // this is an MTCPT so we end here immediately.  We need not set
    // anything more of the iterators state nor adjust 'p' since we
    // return false, so the caller should do no more with this cpt and
    // the random variable values.
    return false;
  }

  void assignDeterministicChild( vector < RV* >& parents,
				 DiscRV* drv )
  {
    drv->val = dt->query(parents,drv);
    if (drv->val >= card()) {
      warning("ERROR: Deterministic CPT '%s' of card %d querying DT '%s' received value %d",
	      name().c_str(),
	      card(),
	      dt->name().c_str(),
	      drv->val);
      fprintf(stderr,"Parents configuration :");
      printRVSetAndValues(stderr,parents);
      error("");
    }
  }


  // used for elimination/triangulation
  virtual unsigned averageCardinality() { return 1; }
  virtual unsigned maxCardinality() { return 1; }


  // random sample given current parents value
  int randomSample(DiscRV* drv) { return (drv->val = _val); }
  
  //////////
  // these routines do nothing for an MTCPT since
  // there is only one possible value for a DT given the parents.
  ///////////////////////////////////////////////////////////  
  // Re-normalize the output distributions
  void normalize() {}
  // set all values to random values.
  void makeRandom() {}
  // set all values to uniform values.
  void makeUniform() {}


  ///////////////////////////////////
  // get parameters from dense pmfs
  unsigned totalNumberParameters() { return 0; }

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
  void emStoreObjectsAccumulators(oDataStreamFile& ofile) {}
  void emLoadObjectsDummyAccumulators(iDataStreamFile& ifile) {}
  void emZeroOutObjectsAccumulators() {}
  void emLoadObjectsAccumulators(iDataStreamFile& ifile) {}
  void emAccumulateObjectsAccumulators(iDataStreamFile& ifile) {}
  const string typeName() { return "DeterministicCPT"; }
  //////////////////////////////////

  void computeParentsSatisfyingChild(
	    // input arguments
	    unsigned par, // parent number
	    vector <RV*> & parents, 
	    vector <RV*> & hiddenParents,
	    PackCliqueValue& hiddenParentPacker,
	    sArray < DiscRVType*>& hiddenNodeValPtrs,
	    RV* child,
	    // output arguments
	    sArray < unsigned >& packedParentVals,
	    unsigned& num)
  {
    return dt->computeParentsSatisfyingChild(par,parents,hiddenParents,hiddenParentPacker,
					     hiddenNodeValPtrs,child,packedParentVals,num);
  }

  void computeParentsChildSatisfyingGrandChild(
	    // input arguments
	    unsigned par, // parent number
	    vector <RV*> & parents, 
	    vector <RV*> & hiddenParents,
	    PackCliqueValue& hiddenParentPacker,
	    sArray < DiscRVType*>& hiddenNodeValPtrs,
	    RV* child,
	    RV* grandChild,
	    // output arguments
	    sArray < unsigned >& packedParentVals,
	    unsigned& num)
  {
    return dt->computeParentsChildSatisfyingGrandChild(par,parents,hiddenParents,hiddenParentPacker,
						       hiddenNodeValPtrs,child,grandChild,
						       packedParentVals,num);
  }



};



#endif // defined MTCPT
