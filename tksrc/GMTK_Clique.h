
/* 
 * GMTK_Clique.h
 * The basic clique data structure.
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

#ifndef GMTK_CLIQUE_H
#define GMTK_CLIQUE_H 

#include <vector>
#include <map>
#include <list>
#include "GMTK_RandomVariable.h"

struct CliqueValue
{
    // A clique has a bunch of possible values, depending on the possible
    // combinations of its constituent variables. With just discrete hidden
    // variables, these can be enumerated out. The probabilities
    // associated with clique values will be stored in these structures.

    logpr lambda, pi;

    // For the clique values associated with non-separators, it is convenient 
    // to keep track of the clique values that are consistent in the surrounding
    // separators. Note that there is just one in each of the surrounding
    // separators. 
    CliqueValue *pred, *succ;

    // The underlying variable values corresponding to a particular clique
    // instantiation
    vector<RandomVariable::DiscreteVariableType> values;
};

struct Clique
{
    vector<RandomVariable *> member;
    // vector of pointers to the constituent variables in the underlying GM.

    vector<RandomVariable *> newMember;
    // vector of pointers to the variables that are present in the clique,
    // but not in its parent.
    // The entries must be organized so that a variable's parents occur
    // before the variable.

    vector<RandomVariable *> discreteMember;
    // A list of the discrete members of the clique.

    vector<RandomVariable *> conditionalProbabilityNode;
    // Each clique has a set of nodes assigned to it, that contribute to
    // its conditional probability. This array stores them.
    // it is set on demand by findConditionalProbabilityNodes()

    vector<RandomVariable::DiscreteVariableType> clampedValues;
    // What are the values of the discrete variables (which should be clamped)
    // set in the course of inference

    void cacheClampedValues();
    // Reads the values of the discrete members and stores them in clampedValues

    void findConditionalProbabilityNodes() 
    {; /* current implementation ignores possible switching efficiencies */}
    // With switching parents, the set of variables assigned to a clique 
    // depends on the value of the clique. This function is called after all
    // the variables have been clamped, and stores the appropriate 
    // RandomVariable pointers in conditionalProbabilityNode.

    bool separator;
    // Is the clique a separator?

    list<CliqueValue> instantiation;
    // This stores all the possible instantiations of a clique.
    // Pruning occurs by removing low probability instantiations.

    map<vector<RandomVariable::DiscreteVariableType>, CliqueValue *> instantiationAddress;
    // A separator clique sums over multiple values from its non-separator
    // parent. instantiationAddress keeps track of the address of the 
    // unique CliqueValue that all the parent values with the same 
    // underlying variable values sum into.

    logpr probGivenParents();
    // This computes the conditional probability of the clique.
    // It calls findConditionalProbabilityNodes(), and then for each
    // member of conditionalProbabilityNode it invokes probGivenParents()
    
    void enumerateValues(int new_member_num, CliqueValue *pred_val,
    bool viterbi=false);
    // This is a recursive function that enumerates all the possible
    // instantiations of a clique.  For consistency with the value of
    // the parent clique, only the new members are instantiated.

    void prune(logpr beam);
    // Removes all instantiations whose forwards probability is less than
    // beam*max.

    void reveal();
    // shows who the members and newMembers are
};

#endif
