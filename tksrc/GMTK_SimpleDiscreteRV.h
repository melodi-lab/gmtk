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

    SimpleDiscreteRV(string _label, vartype vt, int card, bool make_cpt=false)
    : RandomVariable(_label, vt, card) 
    {if (make_cpt) {
        dist_given_parents=new map<vector<int>, vector<logpr> >;
        counts_given_parents=new map<vector<int>, vector<logpr> >;}
    }

    // distribution over own values, given parents values
    // fo no switching. If switching, need a bunch of these
    map<vector<int>, vector<logpr> > *dist_given_parents;

    // this holds the counts accumulated in EM
    map<vector<int>, vector<logpr> > *counts_given_parents;

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
        return (*dist_given_parents)[pvals][val];
    }

    ////////////////////////////////////////////////////////////////////////
    // Sets a variable to the first value that is possible.
    // Values with 0 probability can be ignored, as long as some value is set.
    // This function must always do a clamping.
    // Observation variables (already clamped do nothing).
    // It is assumed that the parent values are clamped at this point.
    virtual void clampFirstValue() {findConditionalParents();
    if (hidden) val = 0;}

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
    virtual void makeRandom() 
    {findConditionalParents(); recMakeRandom(0);}

    ////////////////////////////////////////////////////////////////////////
    // Sets the parameters determining probGivenParents() to uniform values.
    virtual void makeUniform() {error("makeUniform undefined");}

    ////////////////////////////////////////////////////////////////////////
    // Sets the variable according to the probability distribution determined
    // by its parent's values. (i.e., sample from this distribution)
    // Used for simulation.
    virtual void instantiate();

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

    /////////////////////////////////////////////////////////////////////////
    // stores a variable's value elsewhere
    virtual void storeValue(VariableValue &vv) {vv.ival = val;}

    /////////////////////////////////////////////////////////////////////////
    // sets a variables value as specified
    virtual void setValue(VariableValue &vv) {val = vv.ival;}


    /////////////////////////////////////////
    /////////////////////////////////////////
    // SUPPORT FOR EM  VARIABLES ////////////
    /////////////////////////////////////////
    /////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////
    // Called at the beginning of an EM iteration.
    virtual void emStartIteration();

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
    virtual void emIncrement(logpr posterior);

    ////////////////////////////////////////////////////////////////////////
    // At the end of each EM iteration, this is called to convert the 
    // accumulated statistics into probability distributions.
    virtual void emEndIteration();

    /////////////////////////////////////////
    /////////////////////////////////////////
    // END OF SUPPORT FOR EM  VARIABLES /////
    /////////////////////////////////////////
    /////////////////////////////////////////

    void emClearEmAllocatedBit() {;}
    void emClearSwappedBit() {;}
    void emSwapCurAndNew() {;}

  RandomVariable* create() {
    RandomVariable *rv = new SimpleDiscreteRV(label, 
      ((discrete)?(Discrete):(Continuous)), cardinality);
    return rv;
  }

    RandomVariable * clone() 
    { 
      RandomVariable *rv = RandomVariable::clone();
      rv->tieParametersWith(this);
      return rv;
    } 

    void tieParametersWith(RandomVariable *other) 
    {SimpleDiscreteRV *p = (SimpleDiscreteRV *) other;
     dist_given_parents=p->dist_given_parents;
     counts_given_parents=p->counts_given_parents;}
};

#endif
