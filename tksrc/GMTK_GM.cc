
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
#include <set>
#include <algorithm>

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
        node[i]->emClearAllocatedBit();
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
    for (unsigned i=0; i<node.size(); i++)
        node[i]->emClearSwappedBit();
}

void GMTK_GM::emSwapCurAndNew()
{
    for (unsigned i=0; i<node.size(); i++)
        node[i]->emSwapCurAndNew();
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
    logpr last_dp = 0.0;
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
        if (total_data_prob > last_dp)
            emSwapCurAndNew();
        last_dp = total_data_prob;
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
    logpr last_dp = 0.0;
    for (int i=0; i<iterations; i++)
    {
        emStartIteration();
        logpr total_data_prob = 1.0;
        clampFirstExample();
        do
        {
            // first compute the probabilities
            if (!chain->computePosteriors(beam))
            {
                cout << "Skipping example due to 0 probability\n";
                continue;
            }
            total_data_prob *= chain->dataProb;

            // then increment the em counts
            chain->incrementEMStatistics();
        } while (clampNextExample());
        cout << "Total data prob is: " << total_data_prob.val() << endl;
        emEndIteration();
        if (total_data_prob > last_dp)
            emSwapCurAndNew();
        last_dp = total_data_prob;
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

logpr GMTK_GM::cliqueChainExampleProb(vector<vector<VariableValue > > &example,
logpr beam)
{
    setExampleStream(&example);
    logpr tp = 1.0;
    clampFirstExample();
    do
    {
        if (!chain->computePosteriors(beam))
        {
            cout << "Skipping example due to 0 probability\n";
            continue;
        }
        tp *= chain->dataProb;
    } while (clampNextExample());
    return tp;
}

/*-
 *-----------------------------------------------------------------------
 * Function
 *      GM2CliqueChain creates a clique chain corresponding to the GM.
 *      It fully initializes the chain, so that all forms of inference and 
 *      learning can be done.
 * 
 * Preconditions:
 *      GM must be fully initialized
 *
 * Postconditions:
 *      CliqueChain is ready to go
 *
 * Side Effects:
 *      Clique chain created
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */

void GMTK_GM::GM2CliqueChain()
{
    /* A valid clique chain/tree satisfies two constraints (see Zweig thesis):

     * 1)Each variable can be found in at least one clique with all its parents.

     *    - When switching parents are present, the constraint is relaxed: each
     *    variable can be found in at least one clique with 
     *    a) all its switching parents and
     *    b) one occurrence of each set of possible conditional parents

     * 2)The cliques encountered on the (unique) path between two occurrences 
     *   of a variable also contain the variable.
 
     * This routine creates a clique chain according to the following simple
     * algorithm:

     * Maintain a frontier clique that it initialized with a variable that
     * has no parents.
     * Repeat:
     *     - If there are variables in the frontier whose children are all 
     *     present, create a new frontier clique by removing them.
     *     - Otherwise create a new frontier clique by adding a variable whose
     *     parents are all present.
     *     - The new frontier clique is added to the chain as the child of
     *     the old one.

     * This satisfies constraint 1) because a variable is only added
     * when its parents are present. 
     * It satisfies constraint 2) because the variable is added, remains for
     * some time in consecutive cliques, and then is removed.

     * After the clique chain is created, it is condensed by retaining only the 
     * sequence of maximal cliques. Then separator cliques are added.

     * If desired, another routine can be added to split up the cliques so
     * that they satisfy the looser constraints the exist when switching 
     * parents are present. (This should come right before making the
     * separators.)
     */
       
    // This assumes that the variables are present in the node array 
    // according to the order in which they should be added/removed

    // remove a node when all its children have been added
    // use this to keep track of how many have been added.
    map<RandomVariable *, unsigned> num_children_added_for;

    chain = new CliqueChain;

    unsigned maxops = 2*node.size()-1;  // This many add/removes
    vector<Clique> cl(maxops);
    set<Clique *> maximal;  // To keep track of the maximal cliques
    cl[0].member.push_back(node[0]);  // initialize
    cl[0].newMember.push_back(node[0]);
    unsigned num_adds=1;
    assert(node.size() > 0);
    while (num_adds != node.size())
    {
        unsigned i=num_adds;

        // inherit the new nodes of the predecessor
        // accumulate new nodes up to a maximal clique
        // note that they are guaranteed to be in topological order
        cl[i].newMember = cl[i-1].newMember;

        // do all the removes possible
        for (unsigned j=0; j<cl[i-1].member.size(); j++)
        if (num_children_added_for[cl[i-1].member[j]] 
        != cl[i-1].member[j]->allPossibleChildren.size())
        {
            cl[i].member.push_back(cl[i-1].member[j]);
        }
        else
        {
            maximal.insert(&cl[i-1]);
            cl[i].newMember.clear();
        }

        // do an add
        cl[i].member.push_back(node[num_adds]);
        cl[i].newMember.push_back(node[num_adds]);
        for (unsigned j=0; j<node[num_adds]->allPossibleParents.size(); j++)
            num_children_added_for[node[num_adds]->allPossibleParents[j]]++;
        num_adds++;

        // keep the clique members sorted for later intersection
        sort(cl[i].member.begin(), cl[i].member.end());
    }
    maximal.insert(&cl[num_adds-1]);

    // copy the maximal cliques to the clique chain leaving space for separators
    int num_cliques = 2*maximal.size() - 1;  // maximals and separators
    assert(num_cliques>0);
    chain->cliques.resize(num_cliques); 
    set<Clique *>::iterator si;
    int p=0;
    for (si=maximal.begin(); si!=maximal.end(); si++, p+=2)
        chain->cliques[p] = *(*si);

    // make the separators by finding the intersection between successive
    // cliques. 
    for (unsigned i=1; i<chain->cliques.size(); i+=2)
    {
        // make a temporary vector big enough to hold the intersection
        vector<RandomVariable *> tempv(chain->cliques[i-1].member.size());
        vector<RandomVariable *>::iterator se;  // iterator to mark the end

        // get pointers to the ranges to search over
        vector<RandomVariable *>::iterator 
            c1_start= chain->cliques[i-1].member.begin();
        vector<RandomVariable *>::iterator 
            c1_end= chain->cliques[i-1].member.end();
        vector<RandomVariable *>::iterator 
            c2_start= chain->cliques[i+1].member.begin();
        vector<RandomVariable *>::iterator 
            c2_end= chain->cliques[i+1].member.end();

        // do the intersection, allocate memory, and store
        se=set_intersection(c1_start, c1_end, c2_start, c2_end, tempv.begin());
        chain->cliques[i].member.resize(se-tempv.begin());
        copy(tempv.begin(), se, chain->cliques[i].member.begin());

       assert(chain->cliques[i].member.begin()!=chain->cliques[i].member.end());
       
    }
    
    // now initialize all the other clique data members
  
    assert(chain->cliques.size());
    assert(chain->cliques.size()%2);
    for (unsigned i=0; i<chain->cliques.size(); i++)
    {
        Clique *cl = &chain->cliques[i];
        chain->preorder.push_back(cl);
        cl->separator = i%2;

        // the discrete members
        for (unsigned j=0; j<cl->member.size(); j++)
            if (cl->member[j]->discrete)
                cl->discreteMember.push_back(cl->member[j]);

        // the conditional probability nodes
        // set it up conservatively so it will work regardless of switching
        cl->conditionalProbabilityNode = cl->newMember;
    }
    chain->postorder = chain->preorder;
    reverse(chain->postorder.begin(), chain->postorder.end());
}

/*-
 *-----------------------------------------------------------------------
 * Function
 *      showCliques prints out vital information about the cliques
 * 
 * Preconditions:
 *      cliqueChain is set up
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */

void GMTK_GM::showCliques()
{
    for (unsigned i=0; i<chain->preorder.size(); i++)
    {
        cout << "Clique " << i << ":" << endl;
        chain->preorder[i]->reveal();
    }
}

/*-
 *-----------------------------------------------------------------------
 * Function
 *      unroll duplicates one or more time slices a specified number of times 
 * 
 * Preconditions:
 *      # The network up to the template slices must be fully set up
 *      # Variables from a time slice must be contiguous in the node array
 *      # node[i-period] must be analogous to node[i] in the existing network
 *      # timeIndex numbering starts at 0
 *      # Frames from first_frame to last_frame are replicated. For this
 *        to be well-defined, first_frame-1 must be exactly analogous to 
 *        last_frame
 *      # all the prents of variables in the tail must be in the repeating
 *        segments
 *
 * Postconditions:
 *      The network is ready for inference
 *
 * Side Effects:
 *      RandomVariables created and parent pointers set 
 *
 * Results:
 *      nil 
 *
 *-----------------------------------------------------------------------
 */

void GMTK_GM::unroll(int first_frame, int last_frame, int times)
{
    // first frame is the first frame of the repeating segment
    // last frame is the last frame of the repeating segment

    map<int,int> slice_size_for_frame;
    map<int,vector<RandomVariable *>::iterator> start_of_slice, end_of_slice;

    // put the unrolled variables here
    // start with the variables before the tail segment
    vector<RandomVariable *> unrolled;
    for (unsigned i=0; i<node.size(); i++)
        if (node[i]->timeIndex <= last_frame)
            unrolled.push_back(node[i]);

    // find out how big each slice is
    vector<RandomVariable *>::iterator vi;
    for (vi=node.begin(); vi!=node.end(); vi++)
        end_of_slice[(*vi)->timeIndex] = vi; 
    for (vi=node.end()-1; vi>=node.begin(); vi--)
        start_of_slice[(*vi)->timeIndex] = vi;
    for (int i=first_frame; i<=last_frame; i++)
        slice_size_for_frame[i] = end_of_slice[i]-start_of_slice[i]+1;

    // number of slices before unrolling

    // for all existing variables, map their address to the slice and offset
    map<RandomVariable *, pair<int,int> > slice_info_for;
    map<pair<int, int>, RandomVariable *> rv_for_slice;
    int cur_slices = (*node.rbegin())->timeIndex;
    int period = last_frame-first_frame+1;  // over which segments repeat
    repeatPeriod = period;  // store
    for (int i=0; i<=cur_slices; i++)
    {
        int offset = 0;
        for (vi=start_of_slice[i]; vi<=end_of_slice[i]; vi++,offset++)
        {
            int new_slice = i;
            if ((*vi)->timeIndex > last_frame) // part of the tail
               new_slice += times*period; // this is where the tail node will be

            // note where the variable is before unrolling 
            slice_info_for[*vi] = pair<int,int>(i, offset);
 
            // note where it will be after unrolling. 
            // variables in slices up to lasft_frame stay put
            // variables in the tail advance by times*period slices
            rv_for_slice[pair<int,int>(new_slice,offset)] = *vi;
        }
    }
    
    // will need to know which nodes are in the tail segment. record it now
    set<RandomVariable *> in_tail;
    for (unsigned i=0; i<node.size(); i++)
        if (node[i]->timeIndex > last_frame)
            in_tail.insert(node[i]);

    // duplicate the range over and over again
    int cs = last_frame+1;  // slice being created
    unsigned start_of_last_repeated_seg=0, end_of_last_repeated_seg=0;
    for (int i=0; i<times; i++)
    {
        if (i==times-1)
        {
            start_of_last_repeated_seg = unrolled.size();
            startTimeOfLastSegment = cs;
        }
        // make the analogous variables in the new slices
        vi = start_of_slice[first_frame];
        for (int j=first_frame; j<=last_frame; j++, cs++)
            for (int k=0; k<slice_size_for_frame[j]; k++, vi++)    
            {
                RandomVariable *nrv = (*vi)->clone();  // new random variable
                nrv->timeIndex = cs;
                assert(rv_for_slice.find(pair<int,int>(cs, k)) ==
                       rv_for_slice.end());  // should not add twice
                rv_for_slice[pair<int,int>(cs, k)] = nrv; 
                slice_info_for[nrv] = pair<int,int>(cs, k);
                unrolled.push_back(nrv);  
            }
        assert(vi==end_of_slice[last_frame]+1);
        if (i==times-1)
            end_of_last_repeated_seg = unrolled.size()-1;
    }

    // add the tail nodes in the network
    set<RandomVariable *>::iterator si;
    for (si=in_tail.begin(); si!=in_tail.end(); si++)
        unrolled.push_back(*si);

    // will add the parents again to set up the allPossibleChildren and 
    // allPossibleParents arrays, so clear them out now
    for (unsigned i=0; i<unrolled.size(); i++)
    {
        unrolled[i]->allPossibleChildren.clear(); 
        unrolled[i]->allPossibleParents.clear();
    }
   
    // update all the parents in the network
    for (unsigned i=0; i<unrolled.size(); i++)
        if (unrolled[i]->timeIndex <= last_frame)
            continue;  // parents unchanged
        else 
        {
            if (in_tail.count(unrolled[i])==1)
            {
                // the parents in the tail segment point to variables in the 
                // original (unrolled) frames
                // update them to point where they should
                // also update the time index to the right thing
                unrolled[i]->timeIndex += times*period;
            }

            // how many periods ahead are we?
            int periods = (unrolled[i]->timeIndex-last_frame-1)/period + 1; 
            if (in_tail.count(unrolled[i])==1)
                periods = times;  // override

            for (unsigned j=0; j<unrolled[i]->switchingParents.size(); j++)
            {
                pair<int, int> slice_info = 
                    slice_info_for[unrolled[i]->switchingParents[j]];
                int pslice = slice_info.first + periods*period;
                int poffset = slice_info.second;
                unrolled[i]->switchingParents[j] = 
                    rv_for_slice[pair<int,int>(pslice, poffset)];
                if (unrolled[i]->switchingParents[j]==NULL)
                    error("Inconsistency in unrolling. Structure OK?\n");

                // update the data structure for fast splicing
                if (i>=start_of_last_repeated_seg &&i<=end_of_last_repeated_seg)
                {
                    RandomVariable *this_parent = 
                        unrolled[i]->switchingParents[j];
                    parentsToUpdate.push_back(
                        &(unrolled[i]->switchingParents[j]));
                    for (int s=1; s<=times-1; s++)
                    {
                        int slice = this_parent->timeIndex - s*period;
                        pair<RandomVariable *,int> p(this_parent, s);
                        pair<int,int> v(slice, poffset);
                        analogue_k_past[p] = rv_for_slice[v];
                    }
                }
            }

            for (unsigned j=0; j<unrolled[i]->conditionalParentsList.size();j++)
                for (unsigned k=0; 
                k<unrolled[i]->conditionalParentsList[j].size(); k++)
                {
                    pair<int, int> slice_info = 
                    slice_info_for[unrolled[i]->conditionalParentsList[j][k]];
                    int pslice = slice_info.first + periods*period;
                    int poffset = slice_info.second;
                    unrolled[i]->conditionalParentsList[j][k] =
                        rv_for_slice[pair<int,int>(pslice, poffset)];
                    if (unrolled[i]->conditionalParentsList[j][k]==NULL)
                        error("Inconsistency in unrolling. Structure OK?\n");

                    // update the data structure for fast splicing
                    if (i>=start_of_last_repeated_seg 
                    && i<=end_of_last_repeated_seg)
                    {
                        RandomVariable *this_parent = 
                            unrolled[i]->conditionalParentsList[j][k];
                        parentsToUpdate.push_back(
                            &(unrolled[i]->conditionalParentsList[j][k]));
                        for (int s=1; s<=times-1; s++)
                        {
                            int slice = this_parent->timeIndex - s*period;
                            pair<RandomVariable *,int> p(this_parent, s);
                            pair<int,int> v(slice, poffset);
                            analogue_k_past[p] = rv_for_slice[v];
                        }
                    }
                }
        }
    
    // update all the allPossibleChildren arrays -- add the parents 
    // all over again
    for (unsigned i=0; i<unrolled.size(); i++)
        unrolled[i]->setParents(unrolled[i]->switchingParents, 
                            unrolled[i]->conditionalParentsList);

    node = unrolled;
}

void GMTK_GM::spliceOut(int segments)
{
    parentUpdates.clear();
    vector<RandomVariable **>::iterator si;
    for (si=parentsToUpdate.begin(); si!=parentsToUpdate.end(); si++)         
    {
        RandomVariable ** ptr_adr = *si;
        RandomVariable * ptr_val = *ptr_adr;
        pair<RandomVariable **, RandomVariable *> old_values(ptr_adr, ptr_val);
        parentUpdates.push_back(old_values);
        ptr_val = analogue_k_past[pair<RandomVariable *,int>(ptr_val,segments)];
    }
    tempNode = node;
    tempTopo = topologicalOrder;

    int start_purge = startTimeOfLastSegment - segments*repeatPeriod;
    int end_purge = startTimeOfLastSegment-1;
    node.clear();
    for (unsigned i=0; i<tempNode.size(); i++)
        if (tempNode[i]->timeIndex < start_purge 
        || tempNode[i]->timeIndex > end_purge)
            node.push_back(node[i]);
    topologicalOrder.clear();
    for (unsigned i=0; i<tempTopo.size(); i++)
        if (tempTopo[i]->timeIndex < start_purge 
        || tempTopo[i]->timeIndex > end_purge)
            topologicalOrder.push_back(node[i]);
}

void GMTK_GM::restoreNet()
{
    for (unsigned i=0; i<parentUpdates.size(); i++)
        *parentUpdates[i].first = parentUpdates[i].second;
    node = tempNode;
    topologicalOrder = tempTopo;
}
