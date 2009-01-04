/*-
 * GMTK_MSCPT.h
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


#ifndef GMTK_MSCPT_H
#define GMTK_MSCPT_H

#include <vector>

#include "fileParser.h"
#include "logp.h"

#include "GMTK_RngDecisionTree.h"
#include "GMTK_DiscRV.h"
#include "GMTK_CPT.h"
#include "GMTK_Sparse1DPMF.h"
#include "GMTK_NameCollection.h"

#include "GMTK_EMable.h"
#include "GMTK_GMParms.h"
#include "GMTK_NamedObject.h"
#include "GMTK_NameCollection.h"

class MSCPT : public CPT {

  //////////////////////////////////
  // Index into the world structure
  // of the decision tree
  unsigned dtIndex; 

  ///////////////////////////////////////
  // Direct pointer to the decision tree.
  RngDecisionTree* dt;

  ///////////////////////////////////////
  // Direct pointer to collection (indirect
  // mapping to spmfs).
  NameCollection* ncl;

  ///////////////////////////////////////
  // mapping from DT leaves to SPMFs
  NameCollection* spmfCollection;

  ///////////////////////////////////////
  // Index of world's sparse mass function,
  // cached for current value of parents.
  unsigned spmfIndex;

  ///////////////////////////////////////
  // Direct index to the Sparse PMF
  Sparse1DPMF* spmf;

  unsigned _averageCardinality;
  unsigned _maxCardinality;

public:

  ///////////////////////////////////////////////////////////  
  // General constructor
  MSCPT();
  ~MSCPT() { }

  ///////////////////////////////////////////////////////////    
  void setNumParents(const int _nParents);
  void setNumCardinality(const int var, const int card);
  void allocateBasicInternalStructures();

  //////////////////////////////////
  bool iterable() { return dt->iterable(); } 

  virtual unsigned averageCardinality() { return _averageCardinality; }
  virtual unsigned maxCardinality() { return _maxCardinality; }


  ///////////////////////////////////////////////////////////    

  ///////////////////////////////////////////////////////////  
  // Probability evaluation, compute Pr( child | parents ), and
  // iterator support. See GMTK_CPT.h for documentation.
  void becomeAwareOfParentValues( vector < RV *>& parents,
				  const RV* rv ) {
    spmfIndex = dt->query(parents,rv);
    if (!ncl->validSpmfIndex(spmfIndex)) {
      error("ERROR: Sparse CPT '%s' uses DT '%s' with invalid SPMF index '%d' in collection '%s' of size %d\n",
	    name().c_str(),dt->name().c_str(),spmfIndex,
	    ncl->name().c_str(),ncl->spmfSize());
    }
    spmf = ncl->spmf(spmfIndex);
    if (spmf->card() != card()) {
      warning("ERROR: Sparse CPT '%s' of card %d querying DT '%s' received index %d of SPMF '%s' (offset %d in collection '%s') having card %d",
	      name().c_str(),
	      card(),
	      dt->name().c_str(),
	      spmfIndex,
	      spmf->name().c_str(),
	      spmfIndex,ncl->name().c_str(),
	      spmf->card());
      fprintf(stderr,"Parents configuration :");
      printRVSetAndValues(stderr,parents);
      error("");
    }
  }
  void begin(iterator& it, DiscRV* drv,logpr& p) {
    assert ( bitmask & bm_basicAllocated );
    it.setCPT(this);
    it.internalState = -1;
    // store the current spmf
    it.internalStatePtr = (void*) spmf;
    do {
      it.internalState++;
      // We keep the following assertion as we
      // must have that at least one entry is non-zero.
      // The read code of the MDCPT should ensure this
      // as sure all parameter update procedures.
      assert (it.internalState < (int)spmf->length());

      p = spmf->probAtEntry(it.internalState);
      // NOTE: this is a sparse CPT, so the user really shouldn't
      // be placing zeros in a sparse CPT. It might be the case,
      // however, that during parameter training, some CPT entry
      // converges to zero, so we keep this check here, and once
      // it essentially becomes zero, we just skip it over, never
      // placing it in a clique entry.
    } while (p.essentially_zero());
    drv->val = spmf->valueAtEntry(it.internalState); 
    it.drv = drv;
  }
  void becomeAwareOfParentValuesAndIterBegin(vector < RV *>& parents ,
					     iterator &it,
					     DiscRV* drv,
					     logpr& p) {
    spmfIndex = dt->query(parents,drv);
    if (!ncl->validSpmfIndex(spmfIndex)) {
      error("ERROR: Sparse CPT '%s' uses DT '%s' with invalid SPMF index '%d' in collection '%s' of size %d\n",
	    name().c_str(),dt->name().c_str(),spmfIndex,
	    ncl->name().c_str(),ncl->spmfSize());
    }
    Sparse1DPMF* const spmf = ncl->spmf(spmfIndex);
    if (spmf->card() != card()) {
      warning("ERROR: Sparse CPT '%s' of card %d querying DT '%s' received index %d of SPMF '%s' (offset %d in collection '%s') having card %d",
	      name().c_str(),
	      card(),
	      dt->name().c_str(),
	      spmfIndex,
	      spmf->name().c_str(),
	      spmfIndex,ncl->name().c_str(),
	      spmf->card());
      fprintf(stderr,"Parents configuration :");
      printRVSetAndValues(stderr,parents);
      error("");
    }
    assert ( bitmask & bm_basicAllocated );
    it.setCPT(this);
    it.internalState = -1;
    // store the current spmf
    it.internalStatePtr = (void*) spmf;
    do {
      it.internalState++;
      // We keep the following assertion as we
      // must have that at least one entry is non-zero.
      // The read code of the MDCPT should ensure this
      // as sure all parameter update procedures.
      assert (it.internalState < (int)spmf->length());

      p = spmf->probAtEntry(it.internalState);
      // NOTE: this is a sparse CPT, so the user really shouldn't
      // be placing zeros in a sparse CPT. It might be the case,
      // however, that during parameter training, some CPT entry
      // converges to zero, so we keep this check here, and once
      // it essentially becomes zero, we just skip it over, never
      // placing it in a clique entry.
    } while (p.essentially_zero());
    drv->val = spmf->valueAtEntry(it.internalState); 
    it.drv = drv;
  }
  logpr probGivenParents(vector < RV* >& parents,
			 DiscRV* drv) {
    assert ( bitmask & bm_basicAllocated );
    becomeAwareOfParentValues(parents,drv);
    register DiscRVType val = drv->val;
    assert ( val <= card() );
    return spmf->prob(val);
  }
  bool next(iterator &it,logpr& p) {
    Sparse1DPMF* const cur_spmf = (Sparse1DPMF*)it.internalStatePtr;
    do {
      it.internalState++;
      // don't increment past the last value.
      if (it.internalState >= (int)cur_spmf->length()) {
	it.internalState = cur_spmf->length();
	return false;
      }
      p = cur_spmf->probAtEntry(it.internalState);
      // NOTE: this is a sparse CPT, so the user really shouldn't
      // be placing zeros in a sparse CPT. It might be the case,
      // however, that during parameter training, some CPT entry
      // converges to zero, so we keep this check here, and once
      // it essentially becomes zero, we just skip it over, never
      // placing it in a clique entry.
    } while (p.essentially_zero());
    it.drv->val = cur_spmf->valueAtEntry(it.internalState); 
    return true;
  }

  // go through all possible probabilities for this CPT and recompute
  // the max.
  void recomputeCachedMaxValue();


  int randomSample(DiscRV*drv);

  ///////////////////////////////////////////////////////////  
  // Re-normalize the output distributions
  void normalize();
  // set all values to random values.
  void makeRandom();
  // set all values to uniform values.
  void makeUniform();

  ///////////////////////////////////////////////////////////
  // get parameters from dense pmfs, since this heavily
  // depends on the DT mapping.
  unsigned totalNumberParameters() { return 0;  }

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
				  bool writeZeros = false) {}
  void emLoadObjectsDummyAccumulators(iDataStreamFile& ifile) {}
  void emZeroOutObjectsAccumulators() {}
  void emLoadObjectsAccumulators(iDataStreamFile& ifile) {}
  void emAccumulateObjectsAccumulators(iDataStreamFile& ifile) {}
  const string typeName() { return "SparseCPT"; }
  //////////////////////////////////

};



#endif // defined MSCPT
