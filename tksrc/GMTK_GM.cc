
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

/*
 *-------------------------------------------------------------------------
 * Function:
 * findTopologicalOrder
 *    
 * Results:
 * The topologicalOrder array is filled with pointers to the random variables
 * in the graph, in topological order.
 * It is guaranteed that all nodes from time slice i appear in this ordering
 * before any nodes from slice i+1.
 *
 * Side Effects:
 * None.
*/
void GMTK_GM::findTopologicalOrder(randomVariabe *rv)
{

    /* Finds a topological order of the random variables.
       The first call is made with no parameters.
    */

    static int *constraints_for;
    static randomVariable **temp_top0_order;
    static int topo_pos;
    if (rv == NULL)  // the first call
    {
        constraints_for = new int[numNodes];
        temp_topo_order = new randomVariable[numNodes];
        topo_pos = 0;
        for (int i=0; i<numNodes; i++)
        {
            node[i]->findAllPossibleParents();
            constraints_for[i] += node[i]->allPossibleParents.len();
        }

        // now start processing all nodes w/o and predecessors
        for (int i=0; i<numNodes; i++)
            if (constraints_for[i] == 0)
                findTopologicalOrder(&node[i]);

        assert(topo_pos == numNodes);
        delete [] constraints_for;

        // Now "pack" the ordering so that all the nodes from time t occur
        // before any nodes from time t+1. 
        // Since variables are only conditioned on variables earlier in 
        // time, this will be possible.
        assert((numNodes % sliceSize) == 0);
        int *pos = new int[numSlices];   // position in the ith slice
        for (int i=0; i<numSlices; i++)
            pos[i] = 0;
        for (int i=0; i<numNodes; i++)
        {
            int slice = temp_topo_order[i]->timeIndex;
            topologicalOrder[slice_size*slice + pos[slice]++] = 
                temp_topo_order[i];
        }

        delete [] pos;
        delete [] temp_topo_order;
        return;
    }

    // not a first call
    // we have a node with no predecessors
    assert(constraints_for[rv->nodeNum] == 0);
    temp_topo_order[topo_pos++] = rv 
    constraints_for[rv->nodeNum] = -1;  // mark as processed
    for (int i=0; i<rv->allPossibleChildren.len(); i++)
    {
        assert(constraints_for[rv->allPossibleChildren[i]] > 0);
        if (--constrains_for[rv->allPossibleChildren[i]] == 0)
            findTopologicalOrder(rv->allPossibleChildren[i]);
    }
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
    for (int i=0; i<numNodes; i++)
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

void GMTK_GM::enumerateProb(unsigned pos, logpr p)
{
    if (pos == 0)             // first call
        if (!emMode)          // computing the data prob
            dataProb = 0;     // initialize

    if (pos == numNodes)  // all the nodes are instantiated
    {
        if (emMode)
            emIncrementStatistics(p);
        else
            dataProb += p;
        return;
    }

    randomVariable *rv = topologicalOrder[pos];
    if (rv->hidden)
    {
        assert(rv->discrete);
	rv->setPossibleDiscreteValues();
        for (int i=0; i<rv->possibleDiscreteValues.len(); i++)
        {
            rv->val = rv->possibleDiscreteValues[i];
            enumerateProb(pos+1, p*rv->discreteProbGivenParents());
        }
    }
    else
    {
        if (rv->discrete)
            enumerateProb(pos+1, p*rv->discreteProbGivenParents();
        else
            enumerateProb(pos+1, p*rv->probGivenParents();
    }

    if (pos == 0)  // all done with everything
        if (dataProb == 0 && !emMode)
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
    for (int i=0; i<numNodes; i++)
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
    for (int i=0; i<numNodes; i++)
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

void GMTK_GM::enumerateViterbiProb(unsigned pos, logpr p)
{
    if (pos == 0)            // first call
        viterbiProb = 0;     // initialize

    if (pos == numNodes)  // all the nodes are instantiated
    {
        if (p > viterbiProb)
        {
            viterbiProb = p;
            cacheValues();
        }
        return;
    }

    randomVariable *rv = topologicalOrder[pos];
    if (rv->hidden)
    {
        assert(rv->discrete);
	rv->findPossibleDiscreteValues();
        for (int i=0; i<rv->possibleDiscreteValues.len(); i++)
        {
            rv->val = rv->possibleDiscreteValues[i];
            enumerateViterbiProb(pos+1, p*rv->discreteProbGivenParents());
        }
    }
    else
    {
        if (rv->discrete)
            enumerateViterbiProb(pos+1, p*rv->discreteProbGivenParents();
        else
            enumerateViterbiProb(pos+1, p*rv->probGivenParents();
    }

    if (pos == 0)  // all done with everything
        if (viterbiProb == 0)
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
    for (int i=0; i<numNodes; i++)
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
    for (int i=0; i<numNodes; i++)
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
    for (int i=0; i<numNodes; i++)
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
        logpr total_data_prob = 0;
        clampFirstExample();
        do
        {
            // first get the total data probability
            emMode=false;
            enumerativeProb();  
            total_data_prob *= dataProb;

            // then increment the em counts
            emMode=true;
            enumerativeProb(0, 1.0/DataProb);
        } while (clampNextExample);
        cout << "Total data prob is: " << total_data_prob << endl;
        emUpdate();
    }
}
