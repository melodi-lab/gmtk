
/* 
 * GMTK_CliqueChain.h
 * The clique chain data structure.
 * Supports pruning and on-line decoding.
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

#ifndef CLIQUECHAIN_H
#define CLIQUECHAIN_H

#include <list>
#include "GMTK_Clique.h"

struct CliqueChain
{
    vector<Clique> cliques;
    // The actual cliques in the chain.

    vector<Clique *> preorder, postorder;
    // Pointers to the cliques in pre and post-order.

    bool forwardPass(logpr beam=0.0, bool viterbi=false);
    // Computes the alpha probabilities and/or viterbi clique value pointers.
    // Prunes away entries that are less than beam*max.

    void backwardPass();
    // In the backward pass, lambdas are computed for each of the CliqueValues
    // stored in the forward pass. When these are multiplied with the 
    // corresponding pis, and divided by the data prob, the poserior of each
    // clique instantiation results.

    bool doViterbi(logpr beam=0.0);
    // Computes the likeliest value of each clique, and clamps the variables
    // in the underlying network correspondingly.
    // Prunes away entries that are less than beam*max.
    // returns false of the forwards pass fails (0 prob after pruning)

    bool computePosteriors(logpr beam=0.0);
    // Calculates the lambdas and pis for all the cliques.
    // Prunes away entries that are less than beam*max.
    // returns false of the forwards pass fails (0 prob after pruning)

    void incrementEMStatistics();
    // Multiplies the lambdas and pis for each clique instantiation, clamps
    // the network, and increments the EM statistics for each node assigned
    // to a clique.

    logpr dataProb, viterbiProb, backwardDataProb;
};

#endif
