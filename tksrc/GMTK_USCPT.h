/*-
 * GMTK_USCPT.h
 *
 * Special "unity score" cpt, which may only be used with 
 * discrete observations variables with no parents.
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


#ifndef GMTK_USCPT_H
#define GMTK_USCPT_H

#include "fileParser.h"
#include "logp.h"

#include "sArray.h"

#include "GMTK_MDCPT.h"
#include "GMTK_EMable.h"
#include "GMTK_DiscRV.h"
#include "GMTK_NamedObject.h"


// name to use in collections when refering to one of these objects.
#define USMDCPT_NAME "internal:UnityScore"

class USCPT : public MDCPT {

public:

  ///////////////////////////////////////////////////////////  
  // General constructor
  USCPT()
    : MDCPT(di_USCPT)
  {
    _numParents = 0;
    setName(USMDCPT_NAME);
    // we use a special cardinality here just for this
    // type of CPT.
    _card = 0;
    setBasicAllocatedBit();
  }
  ~USCPT() { }

  ///////////////////////////////////////////////////////////    
  // Semi-constructors: useful for debugging.
  // See parent class for further documention.
  void setNumParents(const int _nParents) {}
  void setNumCardinality(const int var, const int card) {}
  void allocateBasicInternalStructures() {}

  ///////////////////////////////////////////////////////////  
  // Probability evaluation, compute Pr( child | parents ), and
  // iterator support.  See GMTK_CPT.h for documentation.
  void becomeAwareOfParentValues( vector <RV *>& parents, const RV* rv ) 
  { /* a USCPT has no parents so do nothing */ }
  void begin(iterator& it,DiscRV* drv,logpr* p) {
    // this routine should never be called for this object, since the variable
    // is always observed.
    assert ( 0 );
    // include code to keep compiler happy
    it.setCPT(this);
  }
  void becomeAwareOfParentValuesAndIterBegin(vector < RV *>& parents, 
					     iterator &it, 
					     DiscRV* drv,
					     logpr& p)
  {
    // this routine should never be called for this object, since the variable
    // is always observed.
    assert ( 0 );
    // include code to keep compiler happy
    it.drv = drv;
    it.setCPT(this);
  }
  logpr probGivenParents(vector < RV*>& parents,
			 DiscRV* drv) 
  {
    logpr val((void*)NULL);
    val.set_to_one();
    return val; 
  }
  bool next(iterator &_it,logpr& p) {
    // this routine should never be called for this object, since the variable
    // is always observed.
    assert ( 0 );
    // include code to keep compiler happy
    return false;
  }


  ///////////////////////////////////
  // don't set drv's value since it is observed.
  int randomSample(DiscRV*drv) { return 0; }

  ///////////////////////////////////
  unsigned totalNumberParameters() { return 0; }

  ///////////////////////////////////////////////////////////  
  // Re-normalize the output distributions
  void normalize() {}
  // set all values to random values.
  void makeRandom() {}
  // set all values to uniform values.
  void makeUniform() {}

  //////////////////////////////////////////////
  // read/write basic parameters
  void read(iDataStreamFile& is) {}
  void write(oDataStreamFile& os) {}

  //////////////////////////////////
  // Public interface support for EM
  //////////////////////////////////
  void emStartIteration() {}
  void emIncrement(logpr p,vector <RV*>& parents,RV*rv) {}
  void emEndIteration() {}
  void emSwapCurAndNew() {}


  // parallel training
  void emStoreAccumulators(oDataStreamFile& ofile) {}
  void emLoadAccumulators(iDataStreamFile& ifile) {}
  void emAccumulateAccumulators(iDataStreamFile& ifile) {}

  void emStoreObjectsAccumulators(oDataStreamFile& ofile) {}
  void emLoadObjectsDummyAccumulators(iDataStreamFile& ifile) {}
  void emZeroOutObjectsAccumulators() {}
  void emLoadObjectsAccumulators(iDataStreamFile& ifile) {}
  void emAccumulateObjectsAccumulators(iDataStreamFile& ifile) {}
  const string typeName() { return "UnityScoreCPT"; }
  //////////////////////////////////

};



#endif // defined USCPT
