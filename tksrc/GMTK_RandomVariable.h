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

struct RandomVariable
{
    // Each variable can have some predecessors and successors.
    // The values stored in the Parent and Child arrays define the 
    // graph topology.
    sArray<RandomVariable *> Parent, Child;

    /* These next members are used to control the inference loops in a
       clique tree. (The constituents of each clique are random variables.)
    */

    bool hidden;
    // If a variable is not hidden, inference has no choices to make about
    // its value.

    bool discrete;
    // If not hidden then inference has no choices to make.
    // If hidden and discrete, inference will loop over all possible values.
    // If hidden and continuous, something else is done.

    virtual void clampFirstValue() {error("ClampFirstValue() Undefined\n");}
    // Sets a variable to the first of its possible values.
    // Values that have 0 probability may be ignored.
    // Observtaion variables have just one possible value.
    // Instantiated variables are responsible for knowing about their values.

    virtual bool clampNextValue() {error("ClampNextValue() Undefined\n");}
    // returns false for observations.
    // returns false when no values are left to be considered.
    // Undefined for hidden continuous nodes.
    // Used to iterate over the possible values of discrete nodes.

    virtual logp probGivenParents() {error("probGivenParents() Undefined\n");}
    // The inference algorithm guarantees that when this is called, the
    // variable and its parents will be clamped to appropriate values.

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

    virtual void increment(logp posterior) = 0;
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
};

#endif
