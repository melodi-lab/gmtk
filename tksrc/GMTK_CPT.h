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
 *  |
 *  +-- USCPT - special unity score CPT for discrete features.
 * 
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
  int _card;

public:

  /////////////////////////////////////////////////////////////////
  // threshold to ensure input data sums to approximately unity
  static double normalizationThreshold;

  enum DiscreteImplementaton { di_MDCPT, di_MSCPT, di_MTCPT, di_USCPT, di_unknown };
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
  virtual void becomeAwareOfParentValues( vector <int>& parentValues,
					  vector <int>& cards) = 0;
  // Another version of becomeAwareOfParentValues but this
  // one explicitely takes an array of random variable parents.
  virtual void becomeAwareOfParentValues( vector < RandomVariable *>& parents ) = 0;
  // return the probability of 'val' given the parents are the
  // assigned to the set of values set during the most previous call to 
  // becomeAwareOfParentValues.
  virtual logpr probGivenParents(const int val) = 0;
  // Similar to the above two. This is convenient for one time
  // probability evaluation.
  virtual logpr probGivenParents(vector <int>& parentValues, 
		                                                                                                                            		 vector <int>& cards, 
				 const int val) = 0;
  virtual logpr probGivenParents(vector < RandomVariable *>& parents,
				 const int val) = 0;

  

  class iterator {
    friend class CPT;
    friend class MDCPT;
    friend class MSCPT;
    friend class MTCPT;
    friend class USCPT;

    // An integer internal state which "hopefully" will
    // be enough for each derived class.
    CPT* cpt;
  protected:
    void setCPT(CPT* _cpt) { cpt = _cpt; }
  public:
    // hold the value right here.
    unsigned value;
    // two variables to be used by CPTs
    // to hold their internal state.
    int internalState;
    void* internalStatePtr;
    // return the value of the RV, the iterators make sure it is set
    // right.
    unsigned val() { return value; }
    // the probability of the variable being value 'val()'
    logpr probVal;

    iterator(CPT* _cpt) : cpt(_cpt) {}
    iterator(const iterator& it) :cpt(it.cpt) 
    { internalState = it.internalState; 
      internalStatePtr = it.internalStatePtr; 
      probVal = it.probVal; }
    // allow to construct an empty iterator, to be filled in later
    iterator() : cpt(NULL) {}

     // prefix
    iterator& operator ++() { cpt->next(*this); return *this; }
    // postfix
    iterator operator ++(int) {
      iterator tmp=*this; ++*this; return tmp;
    }

    bool operator == (const iterator &it) 
       { return it.internalState == internalState; } 
    bool operator != (const iterator &it)
       { return it.internalState != internalState; } 
  };

  // returns an iterator for the first one.
  virtual iterator begin() = 0;

  // creates an iterator for the first one.
  virtual void begin(iterator& it) = 0;
  // returns true if iterate is at end state
  virtual bool end(iterator &it) = 0;
  // Given a current iterator, return true if it
  // is a valid next value, otherwise return false so
  // a loop can terminate.
  virtual bool next(iterator &) = 0;


  // given an internal state of an iterator, return
  // the value corresponding to this random variable.
  virtual int valueAtIt(const int internalState) { return internalState; }


  ///////////////////////////////////////////////////////////  
  // Given the current parent values, generate a random sample.
  virtual int randomSample() = 0;


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

