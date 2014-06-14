/*-
 * GMTK_DeepCPT.h
 *
 *   This CPT gives GMTK the ability to use a deep neural netowrk
 *   to compute probabilities.
 *
 *   This CPT does P(X | \Pi_X = j) = NN(I(j)) where
 *     - \Pi_X are the parents of RV X
 *     - NN is a DeepNN with |Domain(X)| outputs and \Sum_{Y \in \Pi_X} |Domain(Y)| inputs
 *     - I produces an |\Pi_X|-hot vector (0s everywhere except 1s in the positions
 *       corresponding to the value(s) observed for \Pi_X) of length \Sum_{Y \in \Pi_X} |Domain(Y)| 
 *
 *  Written by Richard Rogers <rprogers@ee.washington.edu>
 * 
 * 
 * Copyright (C) 2013 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */


#ifndef GMTK_DEEPCPT_H
#define GMTK_DEEPCPT_H

#include "fileParser.h"
#include "logp.h"

#include "sArray.h"

#include "GMTK_CPT.h"
#include "GMTK_RV.h"
#include "GMTK_DiscRVType.h"
#include "GMTK_NamedObject.h"
#include "GMTK_DeepNN.h"


class DeepCPT : public CPT {


  ////////////////
  // current parent values
  vector <DiscRVType> curParentValues;

  ////////////////
  // the Deep Multi-Layer Perceptron
  DeepNN *dmlp;

  ////////////////
  // remember the computed CPT so we don't have to recompute it
  vector<DiscRVType> cachedParentValues;
  double   *cached_CPT;

  // Apply the deep neural network to get the probability
  logpr applyNN(DiscRV * drv);

public:

  ///////////////////////////////////////////////////////////  
  // General constructor, 
  DeepCPT() : CPT(di_DeepCPT), dmlp(NULL), cached_CPT(NULL)
  { }

  ~DeepCPT() { if (cached_CPT) delete[] cached_CPT; }


  DeepNN *getDeepNN() { return dmlp; }

  // a DeepCPT is considered iterable since its implementation can
  // change not only from segment to segment, but even within a
  // segment.
  bool iterable() { return true; }


  ///////////////////////////////////////////////////////////    
  // Semi-constructors: useful for debugging.
  // See parent class for further documention.
  void setNumParents(const int _nParents) {
    assert ( _nParents == 1 );
    CPT::setNumParents(_nParents);
    curParentValues.resize(_nParents);
    bitmask &= ~bm_basicAllocated;
  }
  void setNumCardinality(const int var, const int card) {
    CPT::setNumCardinality(var,card);
    bitmask &= ~bm_basicAllocated;
  }
  void allocateBasicInternalStructures() {
    error("DeepVECPT::allocateBasicInternalStructures() not implemented");
  }

  ///////////////////////////////////////////////////////////  
  // Probability evaluation, compute Pr( child | parents ), and
  // iterator support. See GMTK_CPT.h for documentation.
  void becomeAwareOfParentValues( vector <RV *>& parents,
				  const RV* rv );
  void begin(iterator& it,DiscRV* drv, logpr& p);
  void becomeAwareOfParentValuesAndIterBegin(vector < RV *>& parents, 
					     iterator &it, 
					     DiscRV* drv,
					     logpr& p);
  logpr probGivenParents(vector <RV *>& parents,
			 DiscRV * drv);
  bool next(iterator &it,logpr& p);

  ///////////////////////////////////
  // random sample given current parent value.
  int randomSample(DiscRV* drv);

  ////////////////////////////////////////////////////////////////////
  // How many parameters does this consume? We return 0 since the
  // virtual evidence does not constitute parameters in the normal
  // sense.
  unsigned totalNumberParameters() { return 0; }


  ///////////////////////////////////////////////////////////  
  // These routines are no-ops in this case since all
  // distributions come from DeepNN.
  void normalize() {}
  // set all values to random values.
  void makeRandom() {}
  // set all values to uniform values.
  void makeUniform() {}

  //////////////////////////////////////////////
  // read/write basic parameters
  void read(iDataStreamFile& is);
  void write(oDataStreamFile& os);

  //////////////////////////////////
  // Public interface support for EM
  //////////////////////////////////
  void emStartIteration() {}
  void emIncrement(logpr p,vector <RV*>& parents,RV*rv) {}
  void emEndIteration() {}
  void emSwapCurAndNew() {}

  // parallel training (DeepVECPTs are never EM trained)
  void emStoreObjectsAccumulators(oDataStreamFile& ofile,
				  bool writeLogVals = true,
				  bool writeZeros = false) {}
  void emLoadObjectsDummyAccumulators(iDataStreamFile& ifile) {}
  void emZeroOutObjectsAccumulators() {}
  void emLoadObjectsAccumulators(iDataStreamFile& ifile) {}
  void emAccumulateObjectsAccumulators(iDataStreamFile& ifile) {}
  const string typeName() { return "DeepCPT"; }
  //////////////////////////////////
};


#endif // defined DeepCPT
