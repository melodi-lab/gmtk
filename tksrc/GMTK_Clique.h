
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

#include <vector>
#include <map>

struct cliqueValue
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
    cliqueValue *pred, *succ;
};

struct clique
{
    sArray<randomVariable *> member;
    // vector of pointers to the constituent variables in the underlying GM.

    sArray<randomVariable *> conditionalProbabilityNode;
    // Each clique has a set of nodes assigned to it, that contribute to
    // its conditional probability. This array stores them.

    void findConditionalProbabilityNodes();
    // With switching parents, the set of variables assigned to a clique 
    // depends on the value of the clique. This function is called after all
    // the variables have been clamped, and stores the appropriate 
    // randomVariable pointers in conditionalProbabilityNode.

    bool separator;
    // Is the clique a separator?

    clique *parent, *child;
    // In a clique chain, that's all there is.

    map<vector<short>, cliqueValue> instantiation;
    // This stores all the possible instantiations of a clique.
    // Pruning occurs by removing low probability instantiations.
    // The vector<short> stores the hidden variable assignments corresponding
    // to the cliqueValue. 

    logpr probGivenParents();
    // This computes the conditional probability of the clique.
    // It calls findConditionalProbabilityNodes(), and then for each
    // member of conditionalProbabilityNode it invokes probGivenParents()
    
};
