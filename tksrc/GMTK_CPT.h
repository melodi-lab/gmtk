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


#include "fileParser.h"
#include "logp.h"
#include "sArray.h"


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
class CPT : public EMable {

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
  // General constructor
  CPT();

  ///////////////////////////////////////////////////////////  
  // Probability evaluation, compute Pr( child | parents )
  // 
  // setParentValues: sets the parent values to a particular
  // assignment. All subsequent calls to to probGivenParents
  // will return the probability of the RV given that the
  // parents are at the particular value.
  virtual void setParentValues( sArray <int>& parentValues ) = 0;
  // return the probability of 'val' given the parents are the
  // assigned to the set of values set during the most previous call to 
  // setParentValues.
  virtual logpr probGivenParents(const int val) = 0;
  // Similar to the above two. This is convenient for one time
  // probability evaluation.
  virtual logpr probGivenParents(sArray <int>& parentValues, 
				 const int val) = 0;

  // Returns an integer that gives, for the current parent assignment
  // the number of possible (i.e., with non-zero probability)
  // different values of this random variable.
  virtual int numValsGivenParents() = 0;
  class iterator {
    friend class CPT;
    int internalState;
  public:
    // the value of the variable
    int val;
    // the probability of the variable being value 'val'
    logpr probVal;
  };
  // returns an iterator for the first one.
  virtual iterator first() = 0;
  // Given a current iterator, return the next one in the sequence.
  virtual iterator next(iterator &) = 0;

  
  ///////////////////////////////////////////////////////////  
  // Re-normalize the output distributions
  virtual void normalize() = 0;

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

