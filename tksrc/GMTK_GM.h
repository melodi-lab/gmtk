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

#ifndef GMTK_GM_H
#define GMTK_GM_H

#include "sArray.h"
#include "logp.h"
#include "GMTK_RandomVariable.h"
#include "GMTK_CliqueChain.h"

struct GMTK_GM
{
    sArray<RandomVariable *> node;
    // This holds all the variables in the graph.
    // The topology os determined by the Parent and Child arrays associated
    // with each random variable.

    sArray<RandomVariable *> topologicalOrder;
    // A topological ordering of the nodes; useful for simulation and
    // enumerative inference.

    void findTopologicalOrder(RandomVariable *rv = NULL);
    // Guarantees that all variables at time t-1 occur before any variables
    // at time t.

    void reveal(sArray<RandomVariable *> order);
    // Go through the nodes in the specified order and show them.

    int numEquivalenceClasses;
    // The variables in the network may be divided into equivalence classes,
    // with the members of each equivalence class sharing the same parameters.

    sArray<RandomVariable *> representativeOfEquivalenceClass;
    // All the members of an equivalence class will use the parameters 
    // associated with this particular member. To save memory, the other
    // member's parameter data structures simply point to this representative.

    void makeRandom();
    // Goes over each variable in the graph and calls its makeRandom function.

    void makeUniform();
    // Goes over each variable in the graph and calls its makeUniform function.

    void simulate();
    // Goes over the variables in topological order and instantiates them 
    // according to the probability distribution given their parent's values.

    void enumerateProb(int pos=0, logpr p=1.0);
    // A recursive function that enumerates all possible values of the 
    // hidden variables and computes the total data prob by brute force.
    // Good for verifying DP and other inference results.
    // Call with emMode=false to get the data probability.
    // Call with emMode=true and p=1/dataProb to do EM

    bool emMode;
    // Is enumerateProb() being called as part of EM?

    void cacheValues();
    // Tells all the nodes int the graph to cache the values they are 
    // currently set to.
    // Observation variables can ignore this.

    void restoreCachedValues();
    // Tells all the nodes in the graph to clamp themselves to their cached
    // values.
    // Observation variables can ignore this.

    void enumerateViterbiProb(int pos=0, logpr p=1.0);
    // Computes the probability of the likeliest instantiation of the hidden
    // variables, and stores it in logViterbiProb. 
    // Has the side effect that at termination, the network is clamped to its
    // liekliest value.

    logpr dataProb, viterbiProb;

    int numNodes;
    // how many nodes are in the graph, i.e. node.len()

    int sliceSize;
    // how many nodes in a time slice

    int numSlices; 
    // how many time slices there are in the network

    void emIncrementStatistics(logpr posterior);
    // Tells each variable to increment its counts by the posterior.

    void emInitialize();
    // Tells each variable to zero out its accumulators and get ready to 
    // start EM.

    void emUpdate();
    // Tells each variable to update its parameters at the end of an EM 
    // iteration.

    void clampFirstExample() {error("can't clamp examples yet");}
    // Clamps the observation variables according to the first example.

    bool clampNextExample() {error("can't clamp examples yet"); return false;}
    // Clamps the observation variables according to the next example.

    void enumerativeEM(int iterations);
    // Does EM using brute force inference.

    CliqueChain *chain;
    // A pointer to a clique chain representation of the GM.

    void cliqueChainEM(int iterations, logpr beam=0.0);
    // Does EM using dynamic programming on a clique chain.
};

#endif
