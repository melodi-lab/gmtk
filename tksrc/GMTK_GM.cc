
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
 * logDataProb is computed. 
 *
 * Side Effects:
 * The hidden variables in the network are left clamped to some arbitrary value
*/

void GMTK_GM::enumerateProb(unsigned pos, logpr p)
{
    if (pos == 0)         // first call
        logDataProb = 0;  // initialize

    if (pos == numNodes)  // all the nodes are instantiated
    {
        logDataProb += p;
        return;
    }

    randomVariable *rv = topologicalOrder[pos];
    if (rv->hidden)
    {
        assert(rv->discrete);
        for (int i=0; i<rv->numVals; i++)
        {
            rv->val = i;
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
}
