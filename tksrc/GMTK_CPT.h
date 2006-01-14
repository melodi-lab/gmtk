/*-
 * GMTK_CPT
 *      .h file the .cc file.
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


#ifndef GMTK_CPT_H
#define GMTK_CPT_H

#include <vector> 

#include "logp.h"
#include "sArray.h"
#include "fileParser.h"

#include "GMTK_NamedObject.h"
#include "GMTK_EMable.h"

class RV;
class DiscRV;

class DiscRV;

/*
 * Generic interface class to all CPT random variables.
 * This gives the probability values for a discrete child with
 * a number of discrete parents. 
 * 
 * Expected class structure.
 * 
 * CPT
 *  |
 *  +-- MDCPT - Multi-d Dense CPT
 *  |
 *  +-- MSCPT - Multi-d Sparse (decision tree based) CPT
 *  |
 *  +-- MTCPT - Multi-d deTerminisitc (decision tree based) CPT
 *  |     |
 *  |     +-- USCPT - special unity (1) score CPT for discrete features.
 *  |
 *  +-- NGramCPT - ARPA ngram language model with backing-off support
 *  |
 *  +-- FNGramCPT - factored language model as in SRILM 1.4+
 *  |
 *  +-- VECPT - virtual evidence CPT
 *  | 
 *  +-- LatticeEdgeCPT - implementation of an HTK lattice edge as a GMTK CPT
 *  +-- LatticeNodeCPT - implementation of an HTK lattice node as a GMTK CPT
 *
 *  Still TODO:
 *      Neural Network/MLP CPT
 *
 */

class MDCPT;
class MSCPT;
class MTCPT;
class USCPT;

class CPT : public EMable {

protected:

  ///////////////////////////////////////////////////////////  
  // The number of "parents" of this CPT, so if we were
  // to expand this out, the dimensionality of this discrete scalar variable
  // and its parents 
  unsigned _numParents;

  ///////////////////////////////////////////////////////////  
  // issue a warning if the number of parents becomes greater than this.
  static unsigned warningNumParents;

  ///////////////////////////////////////////////////////////  
  // The cardinality of each variable, this array is
  // of size (_numParents).
  // cardinality[0] = the cardinality of the first parent
  // cardinality[1] = the cardinality of the 2nd parent
  //    ...
  // cardinality[_nuParents-1] = cardinality of the last parent
  vector < unsigned > cardinalities;

  // cardinality of self (the child)
  unsigned _card;

public:

  /////////////////////////////////////////////////////////////////
  // threshold to ensure input data sums to approximately unity
  static double normalizationThreshold;

  enum DiscreteImplementaton {
    di_MDCPT, // Dense CPT
    di_MSCPT, // Sparse CPT
    di_MTCPT, // Deterministic CPT
    di_USCPT, // Unity Score CPT 
    di_NGramCPT, // Ngram "language model" CPT
    di_FNGramCPT, // factored ngram "language model" CPT
    di_VECPT,    // Virtual Evidence CPT
    di_LatticeNodeCPT,	// lattice node CPT
    di_LatticeEdgeCPT,	// lattice edge CPT
    di_unknown
  };
  const DiscreteImplementaton cptType;


  ///////////////////////////////////////////////////////////  
  // General constructor, does nothing actually.
  CPT(const DiscreteImplementaton _cptType) : cptType(_cptType) { _card = 0; }
  virtual ~CPT() {}

  ////////////////////////////////////////////////
  unsigned numParents() { return _numParents; }

  // Returns the cardinality of this CPT (i.e.,
  // the number of possible values (either with zero or non zero
  // probab) that a RV with this CPT may take on.
  unsigned card() {
    return _card;
  }

  // used for elimination/triangulation
  virtual unsigned averageCardinality() { return (unsigned)card(); }
  virtual unsigned maxCardinality() { return (unsigned)card(); }

  unsigned parentCardinality(const unsigned int par) {
    assert ( par < numParents() );
    return cardinalities[par];
  }

  ///////////////////////////////////////////////////////////    
  // Semi-constructors: useful for debugging.
  // Functions to force the internal structures to be particular values.
  // Force the number of parents to be such.
  virtual void setNumParents(const unsigned _nParents);
  // Set the cardinality. If var = _numParents, this sets
  // the cardinality of the child. Otherwise, it sets the
  // cardinality of the parent. 
  virtual void setNumCardinality(const unsigned var, const int card);
  // Allocate memory, etc. for the internal data structures
  // for this CPT, depending on current _numParents & cardinalities.
  virtual void allocateBasicInternalStructures() = 0;
  // compare the cardinalities of this CPT with that of an other. REturn
  // true if they are equil false otherwise.
  bool compareCardinalities(CPT& cpt);
  // return true if it is possible for this CPT implementation to
  // change from one utterance to the next. Child classes can override
  // this default definition if they so choose.
  virtual bool iterable() { return false; }


  class iterator {
    friend class CPT;
    friend class MDCPT;
    friend class MSCPT;
    friend class MTCPT;
    friend class USCPT;
    friend class NGramCPT;
    friend class FNGramCPT;
    friend class VECPT;


    // The cpt for this iterator.
    CPT* cpt;
  protected:
    void setCPT(CPT* _cpt) { cpt = _cpt; }
  public:
    // Three fields that a CPT subclass can
    // modify as desired to implement a CPT.
    // The actual discrete random variable that this
    // iterator affects. It will use this to set the RVs
    // value and get its cardinality, etc. 
    DiscRV* drv;
    // some internal state pointer to be used.
    void* internalStatePtr;
    // some general scratch internal state, usable
    // as an int or unsigned.
    union {
      int internalState;
      unsigned uInternalState;
    };

    iterator(CPT* _cpt) : cpt(_cpt) {}
    iterator(const iterator& it) :cpt(it.cpt) 
    { internalState = it.internalState; 
      drv = it.drv;
      internalStatePtr = it.internalStatePtr; 
    }
    // allow to construct an empty iterator, to be filled in later
    iterator() : cpt(NULL) {}

  };


  ///////////////////////////////////////////////////////////  
  // Probability evaluation, compute Pr( child | parents ), and
  // iterator support.

  // becomeAwareOfParentValues: sets the parent values to a particular
  // assignment. All subsequent calls to to begin/next will
  // return the probability of the RV given that the parents are at
  // the particular value. This routine assumes (as does similar ones
  // below) that all RVs are really DiscRVs.
  virtual void becomeAwareOfParentValues( vector < RV *>& parents,
					  const RV* rv) = 0;
  // Creates an iterator for the first value of a random variable
  // using this CPT. Note that becomeAwareOfParentValues() MUST have
  // been called before calling the begin iterator.
  virtual void begin(iterator& it,DiscRV* drv,logpr& p) = 0;
  // The above two together in one call with only one virtual dispatch.
  virtual void becomeAwareOfParentValuesAndIterBegin(vector < RV *>& parents, 
						     iterator &it, 
						     DiscRV* drv,
						     logpr& p) = 0;


  // Returns the probability of 'drv->val' given the values of
  // the given parent RVs.  Note that this routine re-evaluates parent
  // values (and does not need becomeAwareOfParentValues to be
  // called), so it is convenient for one time probability evaluation.
  virtual logpr probGivenParents(vector < RV* >& parents,
				 DiscRV* drv) = 0;

  // Given a current iterator, return true if it is a valid next
  // value, otherwise return false so a loop can terminate. If next()
  // returns true, then it also sets the argument p to the next
  // probability value, and otherwise does nothing to p. If next()
  // returns false, then nothing more should be assumed about either
  // the CPT or the value of the caller RV. Note that if next returns
  // true, then p is guaranteed not to be zero probability (meaning it
  // is up to the CPT iterator to advance to the next non-zero
  // probability value).
  virtual bool next(iterator &,logpr& p) = 0;


  // If this CPT is deterministic function, then this routine sets the
  // child as a function of the current values of the parents.
  virtual void assignDeterministicChild(vector < RV *>& parents,
					DiscRV* drv) {
    // unless child class re-defines this routine, it is always
    // an error to call it. 
    coredump("INTERNAL ERROR: invalidly called RV::assignDeterministicChild()");
  }


  ///////////////////////////////////////////////////////////  
  // Given the current parent values, generate a random sample.
  virtual int randomSample(DiscRV* drv) = 0;

  //////////////////////////////////
  // Public interface support for EM that subclasses must support.
  //////////////////////////////////
  void emStartIteration() = 0;
  void emIncrement(logpr p,vector <RV*>& parents,RV*rv) = 0;
  void emEndIteration() = 0;
  void emSwapCurAndNew() = 0;

  ///////////////////////////////////////////////////////////  
  // Re-normalize the output distributions
  virtual void normalize() = 0;
  // set all values to random values.
  virtual void makeRandom() = 0;
  // set all values to uniform values.
  virtual void makeUniform() = 0;


  ///////////////////////////////////////////////////////////    
  // read in the basic parameters, assuming file pointer 
  // is located at the correct position.
  virtual void read(iDataStreamFile& is) = 0;
  ///////////////////////////////////////////////////////////    
  // write out the basic parameters, starting at the current
  // file position.
  virtual void write(oDataStreamFile& os) = 0;

};


#endif
