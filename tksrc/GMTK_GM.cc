
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
#include "GMTK_ObservationMatrix.h"

#include <stdio.h>

#include <map>
#include <set>
#include <algorithm>

#include "GMTK_GMParms.h"
// the observation matrix is referred to in the functions for clamping examples
#include "GMTK_ObservationMatrix.h"


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
        node[i]->instantiate();
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

    RandomVariable *rv = node[pos];
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
        {
            assert(node[i]->val < node[i]->cardinality);
            node[i]->storeValue(vv[p++]);
        }
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
        {
            node[i]->setValue(vv[p++]);
            assert(node[i]->val < node[i]->cardinality);
        }
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

    RandomVariable *rv = node[pos];
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
	printf("Total data prob is: %1.9e\n",total_data_prob.val());
        ::fflush(stdout);
        GM_Parms.emEndIteration();
        if (total_data_prob > last_dp)
	  GM_Parms.emSwapCurAndNew();
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
        logpr total_data_prob = 1.0;
        clampFirstExample();
        do
        {
	     if (globalObservationMatrix.active()) {
	        globalObservationMatrix.printSegmentInfo();
                ::fflush(stdout);
	     }
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
	printf("Total data prob is: %1.9e\n",total_data_prob.val());
        ::fflush(stdout);
        GM_Parms.emEndIteration();
        // if (total_data_prob > last_dp)
	GM_Parms.emSwapCurAndNew();

	/////////////////////////////////////////////////////////
	// write out the basic parameters after each iteration 
	// for debugging purposes. 
	{
	  char buff[1024];
	  ::sprintf(buff,"%d",i);
	  string fname = string("gmparms_it") + buff + string(".gmp");
	  oDataStreamFile of(fname.c_str());
	  GM_Parms.writeAll(of);
	}
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

        // size the cached value array
        cl->clampedValues.resize(cl->discreteMember.size());

        // the conditional probability nodes
        // set it up conservatively so it will work regardless of switching
        // in all cases, separators should never have and conditional 
        // probability nodes assigned to them!!
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

void GMTK_GM::cloneVariables(vector<RandomVariable *> &from,
vector<RandomVariable *> &to)
{
    if (to.size() > 0)
        error("Cloning variables into a non-empty array. Memory Leak.");

    map<RandomVariable *, RandomVariable *> analogue_for;
    for (unsigned i=0; i<from.size(); i++)
    {
        RandomVariable *rv = from[i]->clone();
        analogue_for[from[i]] = rv;
        to.push_back(rv);
    }

    // adjust the internal pointers to point to the new set of variables
    for (unsigned i=0; i<to.size(); i++)
    {
        RandomVariable *rv = to[i];
        for (unsigned j=0; j<rv->switchingParents.size(); j++)
            rv->switchingParents[j] = 
                analogue_for[rv->switchingParents[j]];
        for (unsigned j=0; j<rv->allPossibleParents.size(); j++)
            rv->allPossibleParents[j] = 
                analogue_for[rv->allPossibleParents[j]];
        for (unsigned j=0; j<rv->allPossibleChildren.size(); j++)
            rv->allPossibleChildren[j] = 
                analogue_for[rv->allPossibleChildren[j]];
        for (unsigned j=0; j<rv->conditionalParentsList.size(); j++)
            for (unsigned k=0; k<rv->conditionalParentsList[j].size(); k++)
                rv->conditionalParentsList[j][k] =
                    analogue_for[rv->conditionalParentsList[j][k]];
    } 
}



/*-
 *-----------------------------------------------------------------------
 * Function
 *      SetupForVariableLengthUnrolling() initializes some data structures
 *      that are used to do unrolling. 
 * 
 * Preconditions:
 *      The initial network template must be fully initialized, including
 *      parent arrays, CPTs, gaussians, etc.
 *
 * Postconditions:
 *      gmTemplate array set up.
 *
 *-----------------------------------------------------------------------
 */

void GMTK_GM::setupForVariableLengthUnrolling(int first_frame, int last_frame)
{
    // create the template of the GM to be used for unrolling

    cloneVariables(node, gmTemplate);

    firstChunkFrame = first_frame;
    lastChunkFrame = last_frame;

    obsInTemplate=obsInRepeatSeg=0;
    for (unsigned i=0; i<node.size(); i++)
    {
        if (node[i]->timeIndex>=first_frame 
        && node[i]->timeIndex<=last_frame && !node[i]->hidden)
            obsInRepeatSeg++;
        if (!node[i]->hidden)
            obsInTemplate++;
    }
}


/*-
 *-----------------------------------------------------------------------
 * Function
 *      unroll duplicates one or more time slices a specified number of times 
 * 
 * Preconditions:
 *      # The network up to the template slices must be fully set up
 *      # This includes setting up the CPTs 
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

    if (gmTemplate.size() == 0)  
        error("Must call setupForVariableLengthUnrolling() before unroll");

    // restore the template in preparation for the unrolling
    for (unsigned i=0; i<node.size();i++)
        delete node[i];
    node.clear();
    cloneVariables(gmTemplate, node);

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
    vector<RandomVariable *> tail_node;
    for (unsigned i=0; i<node.size(); i++)
        if (node[i]->timeIndex > last_frame)
        {
            in_tail.insert(node[i]);
            tail_node.push_back(node[i]);
        }

    // duplicate the range over and over again
    int cs = last_frame+1;  // slice being created
    for (int i=0; i<times; i++)
    {
        // make the analogous variables in the new slices
        vi = start_of_slice[first_frame];
        for (int j=first_frame; j<=last_frame; j++, cs++)
            for (int k=0; k<slice_size_for_frame[j]; k++, vi++)    
            {
                RandomVariable *nrv = (*vi)->clone();  // new random variable
                nrv->timeIndex = cs;
                nrv->allPossibleParents.clear();  // will set parents again
                nrv->allPossibleChildren.clear(); // ditto 
                assert(rv_for_slice.find(pair<int,int>(cs, k)) ==
                       rv_for_slice.end());  // should not add twice
                rv_for_slice[pair<int,int>(cs, k)] = nrv; 
                slice_info_for[nrv] = pair<int,int>(cs, k);
                unrolled.push_back(nrv);  
            }
    }

    // add the tail nodes in the network
    for (unsigned i=0; i<tail_node.size(); i++)
        unrolled.push_back(tail_node[i]);

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

                }
        }
    
    // update all the allPossibleChildren arrays -- add the parents 
    // all over again
    for (unsigned i=0; i<unrolled.size(); i++)
        unrolled[i]->setParents(unrolled[i]->switchingParents, 
                            unrolled[i]->conditionalParentsList);

    node = unrolled;
}

void GMTK_GM::setSize(int repeat_segs)
{
    // delete the old clique chain -- we're going to make a new one
    if (chain)
        delete chain;

    // unroll the model to the correct length
    unroll(firstChunkFrame, lastChunkFrame, repeat_segs);

    // make a clique chain to go with the new model
    GM2CliqueChain();
}

void GMTK_GM::setExampleStream(char *obs_file_name)
{
    globalObservationMatrix.openFile(obs_file_name);
    using_files = true;
}

void GMTK_GM::clampFirstExample()
{
    CliqueValue::global_val_set.clear();
    if (using_files)
    {
        if (globalObservationMatrix.numSegments()==0)
            error("no data available");
        globalObservationMatrix.loadSegment(0);
        int frames = globalObservationMatrix.numFrames;
	if ((frames-framesInTemplate)%framesInRepeatSeg != 0) {
	  error("ERROR: segment 0 in observation file does not have a number of frames compatible with unrolling specified structure file");
	}
        setSize((frames-framesInTemplate)/framesInRepeatSeg);
        expos = 0;
	GM_Parms.clampNextDTs();
    }
    else
    {
        if (example==NULL) error("Example array not set.");
        if (example->size()==0) error("No examples.");
        assert(((*example)[0].size()-obsInTemplate)%obsInRepeatSeg==0);
        setSize(((*example)[0].size()-obsInTemplate)/obsInRepeatSeg);
        setValues((*example)[0]);
        expos=0;
    }
}

bool GMTK_GM::clampNextExample()
{
    CliqueValue::global_val_set.clear();
    if (using_files)
    {
        if (expos==globalObservationMatrix.numSegments()-1)
            return false;
        expos++;
        globalObservationMatrix.loadSegment(expos);
        int frames = globalObservationMatrix.numFrames;
	if ((frames-framesInTemplate)%framesInRepeatSeg != 0) {
	  error("ERROR: segment %d in observation file does not have a number of frames compatible with unrolling specified structure file",expos);
	}
        setSize((frames-framesInTemplate)/framesInRepeatSeg);
	GM_Parms.clampNextDTs();
    }
    else
    {
        if (expos==example->size()-1) 
           return false;
        expos++;
        assert( ((*example)[expos].size()-obsInTemplate)%obsInRepeatSeg==0);
        setSize( ((*example)[expos].size()-obsInTemplate)/obsInRepeatSeg);
        setValues((*example)[expos]); 
    }
    return true;
}
