/*
 * GMTK_GM.h
 * Defines the basic Graphical Model functions.
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

#ifndef GMTK_GM
#define GMTK_GM

#include "sArray.h"
#include "error.h"
#include "randomVariable.h"
#include "logp.h"

struct GMTK_GM
{
    sArray<randomVariable *> node;
    // This holds all the variables in the graph.
    // The topology os determined by the Parent and Child arrays associated
    // with each random variable.

    sArray<randomVariable *> topologicalOrder;
    // A topological ordering of the nodes; useful for simulation and
    // enumerative inference.

    int numEquivalenceClasses;
    // The variables in the network may be divided into equivalence classes,
    // with the members of each equivalence class sharing the same parameters.

    sArray<randomVariable *> representativeOfEquivalenceClass;
    // All the members of an equivalence class will use the parameters 
    // associated with this particular member. To save memory, the other
    // member's parameter data structures simply point to this representative.

    void findTopologicalOrder(randomVariable *rv = NULL);
    // Guarantees that all variables at time t-1 occur before any variables
    // at time t.

    void makeRandom();
    // Goes over each variable in the graph and calls its makeRandom function.

    void makeUniform();
    // Goes over each variable in the graph and calls its makeUniform function.

    void simulate();
    // Goes over the variables in topological order and instantiates them 
    // according to the probability distribution given their parent's values.

    logpr enumerateProb(unsigned pos=0, logpr p=1);
    // A recursive function that enumerates all possible values of the 
    // hidden variables and computes the total data prob by brute force.
    // Good for verifying DP and other inference results.

    logpr logDataProb;

    int numNodes;
    // how many nodes are in the graph, i.e. node.len()

    int sliceSize;
    // how many nodes in a time slice

    int numSlices; 
    // how many time slices there are in the network
};

#endif
