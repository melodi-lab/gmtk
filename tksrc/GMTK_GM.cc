
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
    for (unsigned i=0; i<node.size(); i++)
        node[i]->makeRandom();
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
    for (unsigned i=0; i<node.size(); i++)
        node[i]->makeUniform();
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
            emIncrement(p);
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
 * storeValues
 *
 * Results:
 *
 * Side Effects:
 * The values of all the observed variables in the network are 
 * put in the vv array
*/

void GMTK_GM::storeValues(vector<VariableValue> &vv)
{
    vv.clear();
    int observations = 0;
    for (unsigned i=0; i<node.size(); i++)
        if (!node[i]->hidden)
            observations++;
    vv.resize(observations);
    int p=0;
    for (unsigned i=0; i<node.size(); i++)
        if (!node[i]->hidden)
            node[i]->storeValue(vv[p++]);
}

/*
 *---------------------------------------------------------------------------
 * Function:
 * setValues
 *
 * Results:
 *
 * Side Effects:
 * The values of all the observed variables in the network are 
 * set according to the values in the vv array
*/

void GMTK_GM::setValues(vector<VariableValue> &vv)
{
    for (unsigned i=0,p=0; i<node.size(); i++)
        if (!node[i]->hidden)
            node[i]->setValue(vv[p++]);
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
 * emIncrement
 *
 * Results:
 * Each variable increments its accumulators by the posterior.
 *
 * Side Effects:
 * None.
*/

void GMTK_GM::emIncrement(logpr p)
{
    for (unsigned i=0; i<node.size(); i++)
        node[i]->emIncrement(p);
}

/*
 *---------------------------------------------------------------------------
 * Function:
 * emStartIteration 
 *
 * Results:
 * Each variable zeros out its accumulators.
 *
 * Side Effects:
 * None.
*/

void GMTK_GM::emStartIteration()
{
    for (unsigned i=0; i<node.size(); i++)
        node[i]->emStartIteration();
}

/*
 *---------------------------------------------------------------------------
 * Function:
 * emEndIteration
 *
 * Results:
 * Each variable updates its parameters at the end of an EM iteration.
 *
 * Side Effects:
 * None.
*/

void GMTK_GM::emEndIteration()
{
    for (unsigned i=0; i<node.size(); i++)
        node[i]->emEndIteration();
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
    for (int i=0; i<iterations; i++)
    {
        emStartIteration();
        logpr total_data_prob = 1.0;
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
        emEndIteration();
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
    for (int i=0; i<iterations; i++)
    {
        emStartIteration();
        logpr total_data_prob = 1.0;
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
        emEndIteration();
    }
}



/*-
 *-----------------------------------------------------------------------
 * Function
 *      reveal 
 * 
 * Preconditions:
 *      none 
 *
 * Postconditions:
 *      none 
 *
 * Side Effects:
 *      prints the state of the variables in the network 
 *
 * Results:
 *      nil 
 *
 *-----------------------------------------------------------------------
 */

void GMTK_GM::reveal(vector<RandomVariable *> order, bool show_vals)
{
    for (unsigned i=0; i<order.size(); i++)
        order[i]->reveal(show_vals);
}

/*-
 *-----------------------------------------------------------------------
 * Function
 *      enumerativeExamplProb 
 * 
 * Preconditions:
 *      network must have initialized conditional probability representations 
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      hidden variable values are in an undefined state
 *      observed variables have the values of the last example in the array 
 *
 * Results:
 *      returns total data log-likelihood
 *
 *-----------------------------------------------------------------------
 */

logpr GMTK_GM::enumerativeExampleProb(vector<vector<VariableValue > > &example)
{
    setExampleStream(&example);
    logpr tp = 1.0;
    emMode=false;
    clampFirstExample();
    do
    {
        enumerateProb();
        tp *= dataProb;
    } while (clampNextExample());
    return tp;
}
    

/*-
 *-----------------------------------------------------------------------
 * Function
 *      cliqueChainExampleProb computes the likelihood of a set of examples
 *      using dynamic programming on a clique chain.
 * 
 * Preconditions:
 *      network conditional probabilities must be initialized
 *      clique chain must be set up
 *      
 * Postconditions:
 *      none. 
 *
 * Side Effects:
 *      hidden variable values are in an undefined state
 *      observed variables have the values of the last example in the array 
 *
 * Results:
 *      returns total data log-likelihood
 *
 *-----------------------------------------------------------------------
 */

logpr GMTK_GM::cliqueChainExampleProb(vector<vector<VariableValue > > &example)
{
    setExampleStream(&example);
    logpr tp = 1.0;
    clampFirstExample();
    do
    {
        chain->computePosteriors();
        tp *= chain->dataProb;
    } while (clampNextExample());
    return tp;
}
