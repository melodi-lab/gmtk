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


TODO:
Break definitions into two grand sections,
  discrete and continuous support.



#ifndef GMTK_RANDOMVARIABLE
#define GMTK_RANDOMVARIABLE

#include "error.h"
#include "sArray.h"
#include "logp.h"

/* The random variable class defines the basic functions that random
   variables must implement. For speed, there are some special non-virtual
   functions for discrete and deterministic variables, which can be handled
   with table lookups.

   Switching parents are guaranteed to be discrete.
   A nodes conditioning parents must be lower indexed than it.
*/

struct RandomVariable
{
    // Each variable can have some predecessors and successors.
    // The values stored in the Parent and Child arrays define the 
    // graph topology.
    // When switching parents are present, the graph topology is taken to
    // contain the union of all parent relations.
    // Induced cycles are prevented by requiring that for some given node
    // ordering, all the switching parents of a variable are lower-numbered
    // that the variable itself.
    // This means that a static graph with a predetermined topological order
    // can be compiled, even when there is switching parentage.
    // Parent and Child arrays refer only to hidden variables. 
  
    sArray<randomVariable *> switchingParents, conditionalParents;
  
    sArray<randomVariable *> allPossibleParents;
    // allPossibleParents is the union of the switchingParents and the
    // all possible conditionalParents.
    // Used to determine topological orderings

    sArray<randomVariable *> allPossibleChildren;
    // The set of variables that use this variable either as a switching
    // parent or as a (possible) conditional parent.

    findConditionalParents();
    // Looks at the values of the switching parents, and sets the 
    // conditional parents appropraitely

    findAllPossibleParents();
    // Iterates through all possible instantiations of the conditioning 
    // parents and unions together the parents. Stores the result in the
    // allPossibleParents array.

    sarray<sparseCPT> discreteCPT;
    // entry i of the sArray gives the CPT to use when the switching parents
    // have instantiation i. i is determined by a multidimensional array
    // indexing function.

    /* These next members are used to control the inference loops in a
       clique tree. (The constituents of each clique are random variables.)
    */

    const bool hidden;
    // If a variable is not hidden, inference has no choices to make about
    // its value.

    const bool discrete;
    // If not hidden then inference has no choices to make.
    // If hidden and discrete, inference will loop over all possible values.
    // If hidden and continuous, something else is done.

    // NOT NEEDED ANYMORE??? 
    // bool deterministic;
    // A special case of discrete variables, where the parents' values 
    // uniquely determine the variable's value. No looping required.

    int val;
    // in the discrete case, the actual value of the variable

    virtual logpr probGivenParents() {error("probGivenParents() Undefined\n");}
    // For continuous variables.
    // The inference algorithm guarantees that when this is called, the
    // variable and its parents will be clamped to appropriate values.
    // The virtualness of this function is expected to be negligible comared
    // to the other work done for continuous variables.

    int address(sArray<randomVariable *>);
    // looks up the values of the the switching parents and computes a 
    // integer index into a 1 dimensional array

    inline logpr discreteProbGivenParents() 
    {findConditionalParents(); return discreteCPT[address(conditionalParents)][i];}
    // For discrete variables.
    // The inference algorithm guarantees that when this is called, the
    // variable and its parents will be clamped to appropriate values.
    // This does a table lookup in the discrete case.  

/* Add a pointer to the appropriate GMTK_DiscretePDF table here */
/* Then discreteProbGivenParents() will follow it. */
    discreteVariable *dvp;


    sArray< int > possibleDiscreteValues;
    // The set of possible discrete values for whatever the
    // parents are currently set to. One may just iterate over this
    // array to obtain the set of values that this random variable
    // may be set to.

    virtual void setPossibleDiscreteValues();
    // Fills the table possibleDiscreteValues based on the set of
    // current parent values. Called only once for each 
    // set of parent values.

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

    inline void discreteIncrement();
    // For discrete variables, incrementing can be done with a table lookup.
    // Avoid a virtual function call.

/* Add a pointer to the appropriate GMTK_DiscretePDF table here */
/* Then discreteIncrement() will follow it to do the increment. */

    virtual void update() = 0;
    // At the end of each EM iteration, this is called to convert the 
    // accumulated statistics into probability distributions.

    int nodeNum;  
    // A unique node index. nodes in a GM are numbered consecutively from 0.

    int timeIndex;
    // What time frame does the node belong to?
    // Counting starts from 0.
};

#endif
