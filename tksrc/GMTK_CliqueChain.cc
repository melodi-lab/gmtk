
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
    for (unsigned i=0; i<preorder.size(); i++)
    {
        preorder[i]->recycleAllCliqueValues();
        preorder[i]->instantiation.clear();
    }
    assert(Clique::nextfree == int(Clique::gip.size())-1);

    preorder[0]->enumerateValues(0, -1, viterbi);
    for (unsigned i=0; i<preorder.size()-1; i++)
    {
        // we are done enumerating the values for preorder[i], and can
        // free the memory used in its instantiationAddress
        if (i%2)  // a separator; instantiationAddress only used in separators
            preorder[i]->instantiationAddress.clear();  

        // prune the low probability entries in preorder[i]
        // do not prune separators
        // this ensures that each entry on the instantaition list of a
        // non-separator clique will have a successor, and that li->succ->lambda
        // can  be dereferenced on the backward pass without an extra check.
        if (!preorder[i]->separator)
            preorder[i]->prune(beam);

        // propagate the surviving entries to preorder[i+1]
        for (unsigned j=0; j<preorder[i]->instantiation.size(); j++)
        {
            CliqueValue &cv = Clique::gip[preorder[i]->instantiation[j]];
            // clamp the values of the variables in the parent clique
            assert(cv.values->size()==preorder[i]->discreteMember.size());
            for (unsigned k=0; k<preorder[i]->discreteMember.size(); k++)
                preorder[i]->discreteMember[k]->val = (*cv.values)[k];

            // compute the clique values of the child clique that are 
            // consistent with this parent instantiation.
            preorder[i+1]->enumerateValues(0, preorder[i]->instantiation[j], 
                viterbi);
        }
    }

    // clean up after the last clique
    postorder[0]->prune(beam);  // didn't get a chance yet
    // postorder[0] not a separator; don't need to clear its instantiationAddrs

    // check if nothing had any probability, or pruning was too drastic
    if (postorder[0]->instantiation.size() == 0) 
    {
        dataProb = viterbiProb = 0.0;
        return false; 
    }

    // look at the probabilities
    logpr sum=0.0, max=0.0;
    for (unsigned j=0; j<postorder[0]->instantiation.size(); j++)
    {
        int inst_idx = postorder[0]->instantiation[j];
        CliqueValue &cv = Clique::gip[inst_idx];
        sum += cv.pi;
        if (cv.pi > max) max = cv.pi;
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
    // compute the lambdas:
    // Non-separator clique values simply get the lambda of the separator
    // clique value that is descended from them.
    // Separator clique values get the sum of the lambdas of the non-separator
    // clique values descended from them, multiplied by the probGivenParents
    // of those descended values. Note that this was cleverly cached on the
    // forward pass.
    // separator lambdas are zero already

    backwardDataProb = 0.0;  // lets see how closely it matches the forward DP.

    // first do the last clique, where the lambdas are all 1.
    Clique *cl = postorder[0];
    for (unsigned j=0; j<cl->instantiation.size(); j++)
    {
        CliqueValue &cv = Clique::gip[cl->instantiation[j]];
        logpr t = cv.lambda;  // retrieve the cached probGivenParents
        cv.lambda = 1.0;
        if (postorder.size() > 1)  // watch out for degenerate case
            Clique::gip[cv.pred].lambda += t;
	else 
            backwardDataProb += cv.lambda*cv.pi;
    }

    // now do the middle cliques, whose instantiations must pull in a lambda 
    // from the instantiations derived from them, and push a lambda to the 
    // instantiations they derive from.
    for (int i=2; i<int(postorder.size())-2; i+=2)
    {
        Clique *cl = postorder[i];
        for (unsigned j=0; j<cl->instantiation.size(); j++)
        {
            CliqueValue &cv = Clique::gip[cl->instantiation[j]];
            // retrieve the cached probGivenParents
            logpr t = cv.lambda;
	    // update non-separator clique value
            cv.lambda = Clique::gip[cv.succ].lambda;
	    // update separator clique value
	    Clique::gip[cv.pred].lambda += cv.lambda*t;
        }
    }

    // now do the root, which does not do any pushing
    if (preorder.size() > 1) 
    {
        cl = preorder[0];
        for (unsigned j=0; j<cl->instantiation.size(); j++)    
        {
              CliqueValue &cv = Clique::gip[cl->instantiation[j]];
	      backwardDataProb += (cv.lambda=Clique::gip[cv.succ].lambda)*cv.pi;
        }
    }
    cout << "Forward data prob: " << dataProb.val() << " Backward data prob: "
    << backwardDataProb.val() << endl;
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
    if (!forwardPass(beam, viterbi))  
    {
        cout << "Zero probability on the forward pass\n";
        return false;
    }

    // first find the likeliest instantiation of the last clique
    logpr maxprob = 0.0;
    int best = -1;
    Clique *cl = postorder[0];
    for (unsigned j=0; j<cl->instantiation.size(); j++)
    {
        CliqueValue &cv = Clique::gip[cl->instantiation[j]];
        if (cv.pi >= maxprob)
        {
            maxprob = cv.pi;
            best = cl->instantiation[j];
        }
    }
    assert(maxprob==viterbiProb); // already computed on forward pass
            
    // then trace backwards and clamp the best values
    for (unsigned i=0; i<postorder.size(); i++)
    {
        cl = postorder[i];
        // clamp the values
        for (unsigned j=0; j<cl->discreteMember.size(); j++)
            cl->discreteMember[j]->val = (*Clique::gip[best].values)[j];
        // set best to the best instantiation of the predecessor clique
        best = Clique::gip[best].pred;
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

bool CliqueChain::computePosteriors(logpr beam)
{
    if (!forwardPass(beam))
    {
        cout << "Zero probability on the forward pass\n";
        return false;
    }
    assert(dataProb != 0.0);
    backwardPass();
    return true;
}

/*-
 *-----------------------------------------------------------------------
 * Function
 *     incrementEMStatistics() multiplies the lambdas and the pis for 
 *     each clique instantiation, clamps the nodes in the network, and
 *     increments the EM statistics for each node assigned to the clique.
 *     It assumes that the separators do not have any conditionalProbability
 *     nodes assigned to them. If it 
 *
 * Results:
 *
 * Side Effects:
 *
 *-----------------------------------------------------------------------
 */

void CliqueChain::incrementEMStatistics()
{
    // go over the cliques and increment the statistics for all the nodes 
    // assigned to each clique

    for (unsigned i=0; i<preorder.size(); i+=2)  // no need to do separators
    {
        Clique *cl = preorder[i];

	// Uncomment for debugging, sum up the posteriors
        // logpr sum;

        // If using floats, the forward and backward data probs may
        // not be the same. it might be better to use an average. When
	// using doubles, this shouldn't be necessary (or in fact
	// have any effect).
        logpr aveDataProb = logpr((void*)NULL, 
            (dataProb.val()+backwardDataProb.val())/2);

        // consider the contribution of each instantiation of the clique
        for (unsigned k=0; k<cl->instantiation.size(); k++)
        {
            CliqueValue &cv = Clique::gip[cl->instantiation[k]];
            // clamp the instantiation values
            for (unsigned j=0; j<cl->discreteMember.size(); j++)
                cl->discreteMember[j]->val = (*cv.values)[j];

            // find the conditional families, given the switching parent values
            cl->findConditionalProbabilityNodes();

            // do the updates
            logpr posterior = cv.lambda*cv.pi/aveDataProb;

	    // 
	    // Uncomment for debugging, sum up the posteriors
	    // sum += posterior;

            for (unsigned j=0; j<cl->conditionalProbabilityNode.size(); j++)
                cl->conditionalProbabilityNode[j]->emIncrement(posterior);
        }
	// Uncomment for debugging, print the sum of the posteriors
	// printf("Sum of posteriors = %f\n",sum.unlog());
    }
}
