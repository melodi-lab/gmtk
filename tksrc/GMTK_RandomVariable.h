/* 
 * GMTK_RandomVariable.h
 * Defines the functions that all random variables must satisfy in order
 * for the inference, sampling, and other generic routines to work.
 *
 * Written by Geoffrey Zweig <gzweig@us.ibm.com>
 * 
 * Copyright (c) 2001, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about the suitability of this software 
 * for any purpose. It is provided "as is" without express or implied warranty.
 *
 */ 

#ifndef GMTK_RANDOMVARIABLE
#define GMTK_RANDOMVARIABLE

#include "error.h"
#include "sArray.h"
#include "logp.h"

/* The random variable class defines the basic functions that random
   variables must implement. 

   Switching parents are guaranteed to be discrete.
   A nodes conditioning parents must be lower indexed than it.
*/

struct RandomVariable
{
    bool discrete;
    // Is the variable discrete?
    // Inference conditions on this; cliques keep track of the values of 
    // their discrete members.

    short val;
    // in the discrete case, the actual value of the variable.
    // Cliques keep track of the values of their discrete members.

    bool hidden;
    // Is the node hidden, or is it an observation.

    // The probability of a variable taking a particular value depends on
    // the values of its parents. There come in two varieties: switching 
    // parents and conditional parents. The values of the switching parents
    // determines an additional set of conditional parents. The variable
    // is conditioned on these conditional parents. Note that a switching 
    // parent may also be a conditional parent.
    // A node's switching and conditional parents must all be lower numbered 
    // than the node itself, thus allowing a topological ordering of 
    // the graph induced by the union of all possible conditioning relationships
    sArray<randomVariable *> switchingParents, conditionalParents;
  
    sArray<randomVariable *> allPossibleParents;
    // allPossibleParents is the union of the switchingParents and the
    // all possible conditionalParents.
    // Used to determine topological orderings

    sArray<randomVariable *> allPossibleChildren;
    // The set of variables that use this variable either as a switching
    // parent or as a (possible) conditional parent.

    void findConditionalParents();
    // Looks at the values of the switching parents, and sets the 
    // conditionalParents array appropriately.

    void findAllPossibleParents();
    // Iterates through all possible instantiations of the conditioning 
    // parents and unions together the parents. Stores the result in the
    // allPossibleParents array.

    int address(sArray<randomVariable *>);
    // Looks up the values of a vector of parents and computes an
    // integer index into a 1 dimensional array.

    virtual logpr probGivenParents() {error("probGivenParents() Undefined\n");}
    // For continuous variables.
    // The inference algorithm guarantees that when this is called, the
    // variable and its parents will be clamped to appropriate values.
    // The virtualness of this function is expected to be negligible comared
    // to the other work done for continuous variables.

    virtual void clampFirstValue() = 0;
    // Sets a variable to the first value that is possible.
    // Values with 0 probability can be ignored, as long as some value is set.
    // This function must always do a clamping.
    // Observation variables (already clamped do nothing).

    virtual bool clampNextValue() = 0;
    // Sets a variable to the next possible value.
    // Returns false when there are no more values or the variable is an
    // observation.

    virtual void makeRandom() = 0;
    // Sets the parameters determining probGivenParents() to random values.

    virtual void makeUniform() = 0;
    // Sets the parameters determining probGivenParents() to uniform values.

    virtual void instantiate() = 0;
    // Sets the variable according to the probability distribution determined
    // by its parent's values.
    // Used for simulation.

    virtual void tieWith(randomVariable *rv) = 0;
    // In a DGM, the analogous occurrences of a variable in different time
    // slices all share the same distributions and accumulators. This function
    // tells a variable to use the data structures associated with just one
    // member of its equivalence class.

    virtual void zeroAccumulators() = 0;
    // Called at the beginning of an EM iteration.

    virtual void increment(logpr posterior) = 0;
    // On the backward pass, this will be called with a posterior probability
    // for the variable and its parents with all possible sets of values.
    // In inference, the GM is used as a global memory, with variable values
    // dynamically changed as inference proceeds. The inference algorithm
    // guarantees that the variable and its parents will have been set to
    // the the appropriate values when increment() is called, so that the
    // variable can just "look" at its parents to see what their values are.
    // Then the count of seeing the variable and its parents with those values
    // is incremented by the posterior amount.
    // For continuous variables, statistics are accumulated.

    virtual void update() = 0;
    // At the end of each EM iteration, this is called to convert the 
    // accumulated statistics into probability distributions.

    int nodeNum;  
    // A unique node index. nodes in a GM are numbered consecutively from 0.

    int timeIndex;
    // What time frame does the node belong to?
    // Counting starts from 0.

    virtual void cacheValue() = 0;
    // It can be useful to tell a variable to cache the value it's currently
    // set to. For example, wheen keeping track of the best values seen
    // so far in enumerativeViterbi(). Calling this function tells a variable
    // to store its current value (wherever it wants).
    // Not expected to be called many times, so just leave it virtual.

    virtual void restoreCachedValue() = 0;
    // Sets the variable's value to its cached value.
};

#endif
