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
    sArray<randomVariable> node;
    // This holds all the variables in the graph.
    // The topology os determined by the Parent and Child arrays associated
    // with each random variable.

    sArray<randomVariable> topologicalOrder;
    // A topological ordering of the nodes; useful for simulation and
    // enumerative inference.

    void findTopologicalOrder();

    void makeRandom();
    // Goes over each variable in the graph and calls its makeRandom function.

    void makeUniform();
    // Goes over each variable in the graph and calls its makeUniform function.

    void simulate();
    // Goes over the variables in topological order and instantiates them 
    // according to the probability distribution given their parent's values.

    logp enumerateProb(unsigned nodenum=0, logp p=1);
    // A recursive function that enumerates all possible values of the 
    // hidden variables and computes the total data prob by brute force.
    // Good for verifying DP and other inference results.

    logp logDataProb;
};

#endif
