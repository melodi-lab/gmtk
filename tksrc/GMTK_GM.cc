
/*
 * GMTK_GM.cc
 * Basic Graphical Model functions.
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

#include "GMTK_GM.h"
#include <map>

/*
 *-------------------------------------------------------------------------
 * Function:
 * verifyTopologicalOrder
 *    
 * Results:
 * This verifies that the nodes were declared in topological order
 *
 * Side Effects:
 * The node ordering is copied to the topological array.
*/
void GMTK_GM::verifyTopologicalOrder()
{
    map<RandomVariable *, unsigned> position_of;
    for (unsigned i=0; i<node.size(); i++)
        position_of[node[i]] = i;
    for (unsigned i=0; i<node.size(); i++)
        for (unsigned j=0; j<node[i]->allPossibleParents.size(); j++)
            if (position_of[node[i]->allPossibleParents[j]] >= i)
                error("topological ordering violated");
    topologicalOrder = node;
}

/*
 *---------------------------------------------------------------------------
 * Function:
 * makeRandom
 *
 * Results:
 * The parameters of all the variables in the network are set randomly.
 *
 * Side Effects:
 * None.
*/

void GMTK_GM::makeRandom()
{
    for (int i=0; i<numEquivalenceClasses; i++)
        representativeOfEquivalenceClass[i]->makeRandom();
}

/*
 *---------------------------------------------------------------------------
 * Function:
 * makeUniform
 *
 * Results:
 * The parameters of all the variables in the network are set uniformly.
 *
 * Side Effects:
 * None.
*/

void GMTK_GM::makeUniform()
{
    for (int i=0; i<numEquivalenceClasses; i++)
        representativeOfEquivalenceClass[i]->makeUniform();
}

/*
 *---------------------------------------------------------------------------
 * Function:
 * simulate
 *
 * Results:
 * Each variable in the network sets its value according to the conditional
 * probability distribution appropriate to its parent values;
 *
 * Side Effects:
 * Each variable in the model is clamped to some value.
*/

void GMTK_GM::simulate()
{
    for (unsigned i=0; i<node.size(); i++)
        topologicalOrder[i]->instantiate();
}

/*
 *---------------------------------------------------------------------------
 * Function:
 * enumerateProb 
 *
 * Results:
 * dataProb is computed if emMode=false
 * em statistics are incremented for all variables if emMode=true 
 *
 * Side Effects:
 * The hidden variables in the network are left clamped to some arbitrary value
*/

void GMTK_GM::enumerateProb(int pos, logpr p)
{
    if (pos == 0)             // first call
        if (!emMode)          // computing the data prob
            dataProb = 0.0;     // initialize

    if (unsigned(pos) == node.size())  // all the nodes are instantiated
    {
        if (emMode)
            emIncrementStatistics(p);
        else
            dataProb += p;
        return;
    }

    RandomVariable *rv = topologicalOrder[pos];
    rv->clampFirstValue();
    do
    {
        enumerateProb(pos+1, p*rv->probGivenParents());
    } while (rv->clampNextValue());

    if (pos == 0)  // all done with everything
        if (dataProb == 0.0 && !emMode)
            cout << "Warning: Data probability is 0\n";
}

/*
 *---------------------------------------------------------------------------
 * Function:
 * cacheValues
 *
 * Results:
 * The current value of each variable is stored in a cache associated with
 * the variable.
 *
 * Side Effects:
 * None.
*/

void GMTK_GM::cacheValues()
{
    for (unsigned i=0; i<node.size(); i++)
        node[i]->cacheValue();
}

/*
 *---------------------------------------------------------------------------
 * Function:
 * restoreCachedValues
 *
 * Results:
 * Each variable sets its value to its cached value.
 *
 * Side Effects:
 * None.
*/

void GMTK_GM::restoreCachedValues()
{
    for (unsigned i=0; i<node.size(); i++)
        node[i]->restoreCachedValue();
}

/*
 *---------------------------------------------------------------------------
 * Function:
 * viterbiProb 
 *
 * Results:
 * viterbiProb is computed. 
 * Each variable is left clamped with its likeliest value.
 *
 * Side Effects:
 * Each variable is left clamped with its likeliest value.
*/

void GMTK_GM::enumerateViterbiProb(int pos, logpr p)
{
    if (pos == 0)              // first call
        viterbiProb = 0.0;     // initialize

    if (unsigned(pos) == node.size())  // all the nodes are instantiated
    {
        if (p > viterbiProb)
        {
            viterbiProb = p;
            cacheValues();
        }
        return;
    }

    RandomVariable *rv = topologicalOrder[pos];
    rv->clampFirstValue();
    do
    {
        enumerateViterbiProb(pos+1, p*rv->probGivenParents());
    } while (rv->clampNextValue());

    if (pos == 0)  // all done with everything
        if (viterbiProb == 0.0)
            cout << "Warning: All instantiations have 0 probability.\n"
                 << "Network not clamped\n";
        else
       	    restoreCachedValues();
}

/*
 *---------------------------------------------------------------------------
 * Function:
 * emIncrementStatistics 
 *
 * Results:
 * Each variable increments its accumulators by the posterior.
 *
 * Side Effects:
 * None.
*/

void GMTK_GM::emIncrementStatistics(logpr p)
{
    for (unsigned i=0; i<node.size(); i++)
        node[i]->increment(p);
}

/*
 *---------------------------------------------------------------------------
 * Function:
 * emInitialize 
 *
 * Results:
 * Each variable zeros out its accumulators.
 *
 * Side Effects:
 * None.
*/

void GMTK_GM::emInitialize()
{
    for (unsigned i=0; i<node.size(); i++)
        node[i]->zeroAccumulators();
}

/*
 *---------------------------------------------------------------------------
 * Function:
 * emUpdate
 *
 * Results:
 * Each variable updates its parameters at the end of an EM iteration.
 *
 * Side Effects:
 * None.
*/

void GMTK_GM::emUpdate()
{
    for (unsigned i=0; i<node.size(); i++)
        node[i]->update();
}

/*
 *---------------------------------------------------------------------------
 * Function:
 * enumerativeEM
 *
 * Results:
 * Does EM and updates parameters of the network.
 *
 * Side Effects:
 * Variable parameters are changed in accordance with EM.
 * Hidden variables are left clamped in an arbitrary configuration.
*/

void GMTK_GM::enumerativeEM(int iterations)
{
    emInitialize();
    for (int i=0; i<iterations; i++)
    {
        logpr total_data_prob = 0.0;
        clampFirstExample();
        do
        {
            // first get the total data probability
            emMode=false;
            enumerateProb();  
            total_data_prob *= dataProb;

            // then increment the em counts
            emMode=true;
            enumerateProb(0, 1.0/dataProb);
        } while (clampNextExample());
        cout << "Total data prob is: " << total_data_prob.val() << endl;
        emUpdate();
    }
}

/*-
 *-----------------------------------------------------------------------
 * Function
 *    cliqueChainEM does EM using dynamic programming on a clique chain
 *
 * Results:
 *
 * Side Effects:
 *    parameters are updated
 *    hidden variables are left set in an undifined way.
 *
 *-----------------------------------------------------------------------
 */

void GMTK_GM::cliqueChainEM(int iterations, logpr beam)
{
    emInitialize();
    for (int i=0; i<iterations; i++)
    {
        logpr total_data_prob = 0.0;
        clampFirstExample();
        do
        {
            // first compute the probabilities
            chain->computePosteriors(beam);
            total_data_prob *= chain->dataProb;

            // then increment the em counts
            chain->incrementEMStatistics();
        } while (clampNextExample());
        cout << "Total data prob is: " << total_data_prob.val() << endl;
        emUpdate();
    }
}

void GMTK_GM::reveal(vector<RandomVariable *> order, bool show_vals)
{
    for (unsigned i=0; i<order.size(); i++)
        order[i]->reveal(show_vals);
}
