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


class RandomVariable;
class DiscreteRandomVariable;


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
 *  +-- NGramCPT - ngram language model with backing-off support
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
  vector < int > cardinalities;

  // cardinality of self (the child)
  // TODO: turn into unsigned.
  int _card;

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
    di_VECPT,    // Virtual Evidence CPT
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
  int card() {
    return _card;
  }
  unsigned ucard() {
    return (unsigned)_card;
  }

  // used for elimination/triangulation
  virtual unsigned averageCardinality() { return (unsigned)card(); }
  virtual unsigned maxCardinality() { return (unsigned)card(); }

  int parentCardinality(const unsigned int par) {
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


  ///////////////////////////////////////////////////////////  
  // Probability evaluation, compute Pr( child | parents )
  // 
  // becomeAwareOfParentValues: sets the parent values to a particular
  // assignment. All subsequent calls to to probGivenParents
  // will return the probability of the RV given that the
  // parents are at the particular value.
  virtual void becomeAwareOfParentValues( vector < RandomVariable *>& parents ) = 0;
  // A version that doesn't use any random variables, useful for
  // debugging, but guaranteed not to be called by any random variable
  // and/or GMTK inference code.
  virtual void becomeAwareOfParentValues( vector <int>& parentValues,
					  vector <int>& cards) = 0;

  // return the probability of 'val' given the parents are the
  // assigned to the set of values set during the most previous call
  // to becomeAwareOfParentValues.
  virtual logpr probGivenParents(const int val) = 0;
  // Similar to the above, but convenient for one time probability
  // evaluation.
  virtual logpr probGivenParents(vector < RandomVariable *>& parents,
				 const int val) = 0;
  // A version that doesn't use any random variables, useful for debugging,
  // but guaranteed not to be called by any random variable and/or GMTK inference
  // code.
  virtual logpr probGivenParents(vector <int>& parentValues, 
				 vector <int>& cards, 
				 const int val) = 0;


  class iterator {
    friend class CPT;
    friend class MDCPT;
    friend class MSCPT;
    friend class MTCPT;
    friend class USCPT;
    friend class NGramCPT;


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
    DiscreteRandomVariable* drv;
    // some internal state.
    int internalState;
    // some internal state pointer to be used.
    void* internalStatePtr;

    // The probability of the variable being value 'val()'. CPT
    // subclasses must make sure it is set correctly.
    logpr probVal;

    iterator(CPT* _cpt) : cpt(_cpt) {}
    iterator(const iterator& it) :cpt(it.cpt) 
    { internalState = it.internalState; 
      drv = it.drv;
      internalStatePtr = it.internalStatePtr; 
      probVal = it.probVal; }
    // allow to construct an empty iterator, to be filled in later
    iterator() : cpt(NULL) {}

  };


  // Creates an iterator for the first value of a random variable
  // using this CPT.
  virtual iterator begin(DiscreteRandomVariable* drv) = 0;
  virtual void begin(iterator& it,DiscreteRandomVariable* drv) = 0;
  virtual void begin(iterator& it,DiscreteRandomVariable* drv,logpr& p) = 0;
  // TODO: remove this next routine and always use the version with p as argument.  
  virtual void becomeAwareOfParentValuesAndIterBegin
  (  vector < RandomVariable *>& parents , iterator &it, 
     DiscreteRandomVariable* drv) = 0;
  virtual void becomeAwareOfParentValuesAndIterBegin
  (  vector < RandomVariable *>& parents , iterator &it, 
     DiscreteRandomVariable* drv,logpr& p) = 0;

  // Given a current iterator, return true if it is a valid next
  // value, otherwise return false so a loop can terminate. Sets
  // the probability in the iterator state. 
  virtual bool next(iterator &) = 0;
  // a version that also sets an argument probability and
  // does *NOT* set the probability value in the iterator.
  virtual bool next(iterator &,logpr& p) = 0;

  // returns true if iterate is at end state
  virtual bool end(iterator &it) = 0;

  // given an internal state of an iterator, return
  // the value corresponding to this random variable.
  virtual int valueAtIt(const int internalState) { return internalState; }


  ///////////////////////////////////////////////////////////  
  // Given the current parent values, generate a random sample.
  virtual int randomSample(DiscreteRandomVariable* drv) = 0;


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
