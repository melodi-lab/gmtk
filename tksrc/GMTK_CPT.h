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


#ifndef GMTK_CPT
#define GMTK_CPT


#include "logp.h"
#include "sArray.h"
#include "fileParser.h"

#include "GMTK_RandomVariable.h"


/*
 * Generic interface class to all CPT random variables.
 * This gives the probability values for a discrete child with
 * a number of discrete parents. 
 * 
 * Expected class structure.
 * 
 * CPT
 *  |
 *  +-- MDCPT - dense
 *  |
 *  +-- SMDCPT - sparse
 *  |
 *  +-- SD_CPT - decision tree
 * 
 * 
 */
class CPT {

protected:

  ///////////////////////////////////////////////////////////  
  // The number of "parents" of this CPT, so if we were
  // to expand this out, the dimensionality of this discrete scalar variable
  // and its parents 
  int numParents;

  ///////////////////////////////////////////////////////////  
  // issue a warning if the number of parents becomes greater than this.
  static int warningNumParents;

  ///////////////////////////////////////////////////////////  
  // The cardinality of each variable, this array is
  // of size (numParents+1).
  // cardinality[numParents]   = the cardinality of the child
  // cardinality[numParents-1] = the cardinality of the first parent
  // cardinality[numParents-2] = the cardinality of the 2nd parent
  // and so on.
  sArray < int > cardinalities;


public:

  ///////////////////////////////////////////////////////////  
  // General constructor, does nothing actually.
  CPT() {}
  virtual ~CPT() {}

  ///////////////////////////////////////////////////////////    
  // Semi-constructors: useful for debugging.
  // Functions to force the internal structures to be particular values.
  // Force the number of parents to be such.
  virtual void setNumParents(const int _nParents) = 0;
  // Set the cardinality. If var = numParents, this sets
  // the cardinality of the child. Otherwise, it sets the
  // cardinality of the parent. 
  virtual void setNumCardinality(const int var, const int card) = 0;
  // Allocate memory, etc. for the internal data structures
  // for this CPT, depending on current numParents & cardinalities.
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
  virtual void becomeAwareOfParentValues( sArray <int>& parentValues ) = 0;
  // Another version of becomeAwareOfParentValues but this
  // one explicitely takes an array of random variable parents.
  virtual void becomeAwareOfParentValues( sArray < randomVariable *>& parents ) = 0;
  // return the probability of 'val' given the parents are the
  // assigned to the set of values set during the most previous call to 
  // becomeAwareOfParentValues.
  virtual logpr probGivenParents(const int val) = 0;
  // Similar to the above two. This is convenient for one time
  // probability evaluation.
  virtual logpr probGivenParents(sArray <int>& parentValues, 
				 const int val) = 0;
  virtual logpr probGivenParents(sArray < randomVariable *>& parents,
				 const int val) = 0;


  // Returns an integer that gives, for the current parent assignment
  // the number of possible (i.e., with non-zero probability)
  // different values of this random variable.
  virtual int numValsGivenParents() = 0;
  class iterator {
    friend class CPT;
    // An integer internal state which "hopefully" will
    // be enough for each derived class.
    int internalState;
  public:
    iterator() {}
    int val;
    // the probability of the variable being value 'val'
    logpr probVal;
    // change the iterator to be the next valid value.
    void next();

  };

  // returns an iterator for the first one.
  virtual iterator first() = 0;
  // Given a current iterator, return true if it
  // is a valid next value, otherwise return false so
  // a loop can terminate.
  virtual bool next(iterator &) = 0;

  // Given the current parent values, generate a random sample.
  virtual int randomSample();


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


#endif // defined GMTK_CPT

