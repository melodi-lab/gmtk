
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

bool CliqueChain::forwardPass(logpr beam=0, bool viterbi=false)
{
    // zero  out the clique instantiations.
    for (int i=0; i<preorder.len(); i++)
        preorder[i]->instantiation.clear();

    preorder[0]->enumerateValues(0, NULL, viterbi);
    for (int i=0; i<preorder.len()-1; i++)
    {
        preorder[i]->prune(beam);
        map<vector<DISCRETE_VARIABLE_TYPE>, CliqueValue>::iterator mi;
        for (mi=preorder[i]->instantiation.begin(); 
        mi!=preorder[i]->end(); mi++) 
        {
            // clamp the values of the variables in the parent clique
            vector<DISCRETE_VARIABLE_TYPE> &vals = (*mi).first;
            assert(vals.size() == preorder[i]->discreteMember.len());
            for (int j=0; j<preorder[i]->discreteMember.len(); j++)
                preorder[i]->discreteMember[j]->val = vals[j];

            // compute the clique values of the child clique that are 
            // consistent with this parent instantiation.
            preorder[i+1]->enumerateValues(0, &(*mi).second, viterbi);
        }
    }

    if (postorder[0]->instantiation.size() == 0) return false; 
    // nothing had any probability, or pruning was too drastic

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
    map<vector<DISCRETE_VARIABLE_TYPE>, CliqueValue>::iterator mi;
    for (int i=1; i<postorder.len(); i+=2)
    {
        Clique *cl = postorder[i];
        for (mi=cl->instantiation.begin(); mi!=cl->instantiation.end(); mi++)
            (*mi).second.lambda = 0;
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
    for (mi=cl->instantiation.begin(); mi!=cl->instantiation.end(); mi++)
    {
        CliqueValue &cv = (*mi).second;
        logpr t = cv.lambda;  // retrieve the cached probGivenParents
        cv.lambda = 1;
        if (postorder.len() > 1)  // watch out for degenerate case
            cv.pred->lambda += t;
    }

    // now do the middle cliques, whose instantiations must pull in a lambda 
    // from the instantiations derived from them, and push a lambda to the 
    // instantiations they derive from.
    for (int i=2; i<postorder.len()-2; i++)
    {
        Clique *cl = postorder[i];
        for (mi=cl->instantiation.begin(); mi!=cl->instantiation.end(); mi++)
        {
            CliqueValue &cv = (*mi).second;
            logpr t = cv.lambda;
            cv.lambda = cv->succ.lambda;
            cv.pred->lambda += cv.lambda*t;
        }
    }

    // now do the root, which does not do any pushing
    cl = preorder[0];
    if (preorder.len() > 1)
        for (mi=cl->instantiation.begin(); mi!=cl->instantiation.end); mi++)
        {
            CliqueValue &cv = (*mi).second;
            cv.lambda =  cv.succ->lambda;
        }
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
    logpr maxprob = 0;
    Clique *cl = postorder[0];
    map<vector<DISCRETE_VARIABLE_TYPE>, CliqueValue>::iterator mi, best;
    for (mi=cl->instantiation.begin(); mi!=cl->instantiation.end(); cl++)
        if ((*mi).second.pi >= maxprob)
        {
            maxprob = (*mi).second.pi;
            best = mi;
        }
    gm->viterbiProb = maxprob;
            
    // then trace backwards and clamp the best values
    for (int i=0; i<postorder.len(); i++)
    {
        // clamp the values
        for (int j=0; j<discreteMember.size(); j++)
            discreteMember[j]->val = (*best).first[j];
        // set best to the best instantiation of the predecessor clique
            
}
