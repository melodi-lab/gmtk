/*-
 * GMTK_LatticeNodeCPT.h
 *     CPT for lattice node transitions
 *
 *  Written by Gang Ji <gang@ee.washington.edu>
 * 
 *  $Header$
 * 
 * Copyright (C) 2001 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 *
 */


#ifndef GMTK_LATTICE_NODE_CPT_H
#define GMTK_LATTICE_NODE_CPT_H


#include "GMTK_CPT.h"
#include "GMTK_LatticeADT.h"


class LatticeNodeCPT : public CPT {

 public:
  LatticeNodeCPT();
  ~LatticeNodeCPT();

  ///////////////////////////////////////////////////////////  
  // virtual functions from class CPT

  // Allocate memory, etc. for the internal data structures
  // for this CPT, depending on current _numParents & cardinalities.
  virtual void allocateBasicInternalStructures() {}

  ///////////////////////////////////////////////////////////
  // Probability evaluation, compute Pr( child | parents ), and
  // iterator support. See GMTK_CPT.h for documentation.
  virtual void becomeAwareOfParentValues(vector< RV* >& parents, const RV* rv) { assert(0); }
  virtual void begin(iterator& it, DiscRV* drv,logpr& p) { assert(0); }
  virtual void becomeAwareOfParentValuesAndIterBegin(vector< RV* >& parents, iterator &it, DiscRV* drv, logpr& p);
  virtual logpr probGivenParents(vector< RV* >& parents, DiscRV* drv);
  virtual bool next(iterator &it,logpr& p);

  // Include here an extra routine that returns the probability
  // of the child 'val' given the parents are the assigned to
  // the set of values set during the most previous call to
  // becomeAwareOfParentValues.
  logpr probGivenParents(const int val) { assert(0); return logpr(0.0); }

  ///////////////////////////////////////////////////////////
  // Given the current parent values, generate a random sample.
  virtual int randomSample(DiscRV*drv) { assert(0); return 0; }

  ///////////////////////////////////////////////////////////
  // Re-normalize the output distributions
  virtual void normalize() {}
  // set all values to random values.
  virtual void makeRandom() {}
  // set all values to uniform values.
  virtual void makeUniform() {}

  ///////////////////////////////////////////////////////////
  // read in the basic parameters, assuming file pointer
  // is located at the correct position.
  virtual void read(iDataStreamFile& is) {}
  ///////////////////////////////////////////////////////////
  // Do nothing.
  virtual void write(oDataStreamFile& os) {}
	
  void setLatticeADT(const LatticeADT &latticeAdt);

  inline bool useTimeParent() const { return _latticeAdt->useTimeParent(); }

  ////////////////////////////////////////////////////////////////////////////
  // from base class EMable
  void emStartIteration() {}
  void emIncrement(logpr prob,vector < RV* >& parents, RV*r) {}
  void emEndIteration() {}
  void emSwapCurAndNew() {}

  // return the number of parameters for object.
  virtual unsigned totalNumberParameters() {return 0;}

  ///////////////////////////////////////////////////////////////
  // virtual functions for objects to do the actual work.
  virtual void emStoreObjectsAccumulators(oDataStreamFile& ofile,
					  bool writeLogVals = true,
					  bool writeZeros = false) {}
  virtual void emLoadObjectsDummyAccumulators(iDataStreamFile& ifile) {}
  virtual void emZeroOutObjectsAccumulators() {}
  virtual void emLoadObjectsAccumulators(iDataStreamFile& ifile) {}
  virtual void emAccumulateObjectsAccumulators(iDataStreamFile& ifile) {}

 protected:
  // returns the type of the sub-object in string
  // form that is suitable for printing and identifying
  // the type of the object.
  virtual const string typeName() {return std::string("LatticeNodeCPT");}

    void printLargerCPT(shash_map_iter<unsigned, LatticeADT::LatticeEdgeList> & cpt, LatticeADT::LatticeNode & cur_node, unsigned cur_node_id);


  const LatticeADT *_latticeAdt;
};


#endif // end GMTK_LATTICE_NODE_CPT_H
