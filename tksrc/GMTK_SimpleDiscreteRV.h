/* 
 * GMTK_SimpleDiscreteRV.h
 * A very simple discrete random variable implementation
 *
 * Written by Geoffrey Zweig <gzweig@us.ibm.com> & Jeff Bilmes <bilmes@ee.washington.edu>
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

#ifndef GMTK_SIMPLERANDOMVARIABLE
#define GMTK_SIMPLERANDOMVARIABLE

#include "GMTK_RandomVariable.h"
#include <vector>
#include <map>
#include <stdlib.h>


class SimpleDiscreteRV : public RandomVariable
{
public:

    SimpleDiscreteRV(char * _label, vartype vt, int card)
    : RandomVariable(_label, vt, card) {;}

    // distribution over own values, given parents values
    // fo no switching. If switching, need a bunch of these
    map<vector<int>, vector<logpr> > dist_given_parents;

    ////////////////////////////////////////////////////////////////////////
    // Looks at the values of the switching parents, and sets the 
    // conditionalParents array pointer appropriately.
    virtual void findConditionalParents()
    { if (conditionalParentsList.size()!=1) error("Can't work w/ switching\n");
        curConditionalParents = &(conditionalParentsList[0]);
    }

    ////////////////////////////////////////////////////////////////////////
    // The inference algorithm guarantees that when this is called, the
    // variable and its parents will be clamped to appropriate values.
    // This also assumes that findConditionalParents has been called.
    virtual logpr probGivenParents() 
    { 
        vector<int> pvals;
        for (unsigned i=0; i<curConditionalParents->size(); i++)
            pvals.push_back((*curConditionalParents)[i]->val);
        return dist_given_parents[pvals][val];
    }

    ////////////////////////////////////////////////////////////////////////
    // Sets a variable to the first value that is possible.
    // Values with 0 probability can be ignored, as long as some value is set.
    // This function must always do a clamping.
    // Observation variables (already clamped do nothing).
    // It is assumed that the parent values are clamped at this point.
    virtual void clampFirstValue() {if (hidden) val = 0;}

    ////////////////////////////////////////////////////////////////////////
    // Sets a variable to the next possible value.
    // Returns false when there are no more values or the variable is an
    // observation.
    // It is assumed that the parent values are clamped at this point.
    virtual bool clampNextValue() 
    {
        if (!hidden) return false;
        if (val == cardinality-1) return false;
        val++;
        return true;
    }

    void recMakeRandom(int pos);

    ////////////////////////////////////////////////////////////////////////
    // Sets the parameters determining probGivenParents() to random values.
    // In general, must do this for each set of switching values
    virtual void makeRandom() {findConditionalParents(); recMakeRandom(0);}

    ////////////////////////////////////////////////////////////////////////
    // Sets the parameters determining probGivenParents() to uniform values.
    virtual void makeUniform() {error("makeUniform undefined");}

    ////////////////////////////////////////////////////////////////////////
    // Sets the variable according to the probability distribution determined
    // by its parent's values. (i.e., sample from this distribution)
    // Used for simulation.
    virtual void instantiate() {val=(rand()%cardinality);}

    ////////////////////////////////////////////////////////////////////////
    // In a DGM, the analogous occurrences of a variable in different time
    // slices all share the same distributions and accumulators. This function
    // tells a variable to use the data structures associated with just one
    // member of its equivalence class.
    virtual void tieWith(RandomVariable *rv) {error("tieWith undefined");}

    ////////////////////////////////////////////////////////////////////////
    // It can be useful to tell a variable to cache the value it's currently
    // set to. For example, wheen keeping track of the best values seen
    // so far in enumerativeViterbi(). Calling this function tells a variable
    // to store its current value (wherever it wants).
    // Not expected to be called many times, so just leave it virtual.
    // This like a push operation, onto a one element size stack.
    int cv;
    virtual void cacheValue() {cv=val;}


    ////////////////////////////////////////////////////////////////////////
    // Sets the variable's value to its cached value.
    // This like a pop operation, off of the stack.
    virtual void restoreCachedValue() {val=cv;}


    /////////////////////////////////////////
    /////////////////////////////////////////
    // SUPPORT FOR EM  VARIABLES ////////////
    /////////////////////////////////////////
    /////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////
    // Called at the beginning of an EM iteration.
    virtual void zeroAccumulators(){error("zeroAccumulators undefined\n");} 


    ////////////////////////////////////////////////////////////////////////
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
    virtual void increment(logpr posterior) {error("increment undefined");}


    ////////////////////////////////////////////////////////////////////////
    // At the end of each EM iteration, this is called to convert the 
    // accumulated statistics into probability distributions.
    virtual void update() {error("update undefined");}


    /////////////////////////////////////////
    /////////////////////////////////////////
    // END OF SUPPORT FOR EM  VARIABLES /////
    /////////////////////////////////////////
    /////////////////////////////////////////
};

#endif
