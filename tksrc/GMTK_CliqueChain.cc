
/* 
 * GMTK_CliqueChain.cc
 * Clique chain functions.
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

#include "GMTK_CliqueChain.h"

/*-
 *-----------------------------------------------------------------------
 * Function
 *     forwardPass computes the pi probabilities, pruning away any 
 *     entries that with probabilities less than beam*max.
 *     In viterbi mode, it does max's instead of sums.
 *
 * Results:
 *     returns true if the last clique has at least one unpruned, non-zero
 *     probability value. returns false otherwise.
 *
 * Side Effects:
 *     The instantiation lists of the cliques are filled up.
 *     dataProb gets the data probability  
 *     viterbiProb gets the viterbi probability
 *
 *-----------------------------------------------------------------------
 */
bool CliqueChain::forwardPass(logpr beam, bool viterbi)
{
    // zero  out the clique instantiations.
    for (int i=0; i<preorder.len(); i++)
        preorder[i]->instantiation.clear();

    preorder[0]->enumerateValues(0, NULL, viterbi);
    for (int i=0; i<preorder.len()-1; i++)
    {
        // we are done enumerating the values for preorder[i], and can
        // free the memory used in its instantiationAddress
        preorder[i]->instantiationAddress.clear();  

        // prune the low probability entries in preorder[i]
        preorder[i]->prune(beam);

        // propagate the surviving entries to preorder[i+1]
        list<CliqueValue>::iterator li;
        for (li=preorder[i]->instantiation.begin(); 
        li!=preorder[i]->instantiation.end(); li++) 
        {
            // clamp the values of the variables in the parent clique
            assert(int(li->values.size()) == preorder[i]->discreteMember.len());
            for (int j=0; j<preorder[i]->discreteMember.len(); j++)
                preorder[i]->discreteMember[j]->val = li->values[j];

            // compute the clique values of the child clique that are 
            // consistent with this parent instantiation.
            preorder[i+1]->enumerateValues(0, &(*li), viterbi);
        }
    }

    // clean up after the last clique
    postorder[0]->prune(beam);  // didn't get a chance yet
    postorder[0]->instantiationAddress.clear();

    // check if nothing had any probability, or pruning was too drastic
    if (postorder[0]->instantiation.size() == 0) 
        return false; 

    // look at the probabilities
    logpr sum=0.0, max=0.0;
    list<CliqueValue>::iterator li;
    for (li=postorder[0]->instantiation.begin(); 
    li!=postorder[0]->instantiation.end(); li++)
    {
        sum += li->pi;
        if (li->pi > max) max = li->pi;
    }
    dataProb=viterbiProb=0.0;
    if (viterbi)
        viterbiProb=max;
    else
        dataProb=sum;

    return true;
}

/*-
 *-----------------------------------------------------------------------
 * Function
 *    backwardPass computes and stores lambdas for each of the CliqueValues
 *    that survived pruning on the forwards pass.  
 *
 * Results:
 *
 * Side Effects:
 *    resets the lambdas in the instantiation objects    
 *
 *-----------------------------------------------------------------------
 */
void CliqueChain::backwardPass()
{
    // zero out the separator lambdas
    list<CliqueValue>::iterator li;
    for (int i=1; i<postorder.len(); i+=2)
    {
        Clique *cl = postorder[i];
        for (li=cl->instantiation.begin(); li!=cl->instantiation.end(); li++)
            li->lambda = 0.0;
    }

    // compute the lambdas 
    // non-separator clique values simply get the lambda of the separator
    // clique value that is descended from them.
    // separator clique values get the sum of the lambdas of the non-separator
    // clique values descended from them, multiplied by the probGivenParents
    // of those descended values. Note that this was cleverly cached on the
    // forward pass.

    // first do the last clique, where the lambdas are all 1.
    Clique *cl = postorder[0];
    for (li=cl->instantiation.begin(); li!=cl->instantiation.end(); li++)
    {
        logpr t = li->lambda;  // retrieve the cached probGivenParents
        li->lambda = 1.0;
        if (postorder.len() > 1)  // watch out for degenerate case
            li->pred->lambda += t;
    }

    // now do the middle cliques, whose instantiations must pull in a lambda 
    // from the instantiations derived from them, and push a lambda to the 
    // instantiations they derive from.
    for (int i=2; i<postorder.len()-2; i++)
    {
        Clique *cl = postorder[i];
        for (li=cl->instantiation.begin(); li!=cl->instantiation.end(); li++)
        {
            logpr t = li->lambda;
            li->lambda = li->succ->lambda;
            li->pred->lambda += li->lambda*t;
        }
    }

    // now do the root, which does not do any pushing
    cl = preorder[0];
    if (preorder.len() > 1)
        for (li=cl->instantiation.begin(); li!=cl->instantiation.end(); li++)
            li->lambda =  li->succ->lambda;
}

/*-
 *-----------------------------------------------------------------------
 * Function
 *    doViterbi does a forward pass in Viterbi mode, and then clamps
 *    all the nodes in the network to their likeliest values. 
 *
 * Results:
 *    Returns true if successful, false if pruning removes all candidates,
 *    or all complete instantiations have 0 probability.
 *
 * Side Effects:
 *    The nodes in the graphical model are clamped to their likeliest values.
 *
 *-----------------------------------------------------------------------
 */

bool CliqueChain::doViterbi(logpr beam)
{
    bool viterbi=true;
    if (!forwardPass(beam, viterbi))  // 
        return false;

    // first find the likeliest instantiation of the last clique
    logpr maxprob = 0.0;
    CliqueValue *best;
    Clique *cl = postorder[0];
    list<CliqueValue>::iterator li;
    for (li=cl->instantiation.begin(); li!=cl->instantiation.end(); li++)
        if (li->pi >= maxprob)
        {
            maxprob = li->pi;
            best = &(*li);
        }
    assert(maxprob==viterbiProb); // already computed on forward pass
            
    // then trace backwards and clamp the best values
    for (int i=0; i<postorder.len(); i++)
    {
        cl = postorder[i];
        // clamp the values
        for (int j=0; j<cl->discreteMember.len(); j++)
            cl->discreteMember[j]->val = best->values[j];
        // set best to the best instantiation of the predecessor clique
        best = best->pred;
    }

    return true;
}

/*-
 *-----------------------------------------------------------------------
 * Function
 *     computePosteriors computes the lambdas and pis for every clique
 *     instantiation. However, it does prune away entries whose pi probability
 *     is less than beam*max 
 *
 * Results:
 *
 * Side Effects:
 *     lambdas and pis set 
 *     hidden nodes left clamped in an undefined configuration
 *
 *-----------------------------------------------------------------------
 */

void CliqueChain::computePosteriors(logpr beam)
{
    if (!forwardPass(beam))
// really need to throw an exception so that the program can continue with the
// next example
        error("Zero probability on the forward pass");
    backwardPass();
}

/*-
 *-----------------------------------------------------------------------
 * Function
 *     incrementEMStatistics() multiplies the lambdas and the pis for 
 *     each clique instantiation, clamps the nodes in the network, and
 *     increments the EM statistics for each node assigned to the clique.
 *
 * Results:
 *
 * Side Effects:
 *
 *-----------------------------------------------------------------------
 */

void CliqueChain::incrementEMStatistics()
{
    list<CliqueValue>::iterator li;

    // go over the cliques and increment the statistics for all the nodes 
    // assigned to each clique
    for (int i=0; i<preorder.len(); i++)
    {
        Clique *cl = preorder[i];
        
        // consider the contribution of each instantiation of the clique
        for (li=cl->instantiation.begin(); li!=cl->instantiation.end(); li++)
        {
            // clamp the instantiation values
            for (int j=0; j<cl->discreteMember.len(); j++)
                cl->discreteMember[j]->val = li->values[j];

            // find the conditional families, given the switching parent values
            cl->findConditionalProbabilityNodes();
            
            // do the updates
            logpr posterior = li->lambda*li->pi/dataProb;
            for (int j=0; j<cl->conditionalProbabilityNode.len(); j++)
                cl->conditionalProbabilityNode[j]->increment(posterior);
        }
    }
}
