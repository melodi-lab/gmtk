/*
 * GMTK_GM.h
 * Defines the basic Graphical Model functions.
 *
 * Written by Geoffrey Zweig <gzweig@us.ibm.com> & Jeff Bilmes <bilmes@ee.washington.edu>
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

#ifndef GMTK_GM_H
#define GMTK_GM_H

#include <set>

#include "logp.h"
#include "GMTK_RandomVariable.h"
#include "GMTK_CliqueChain.h"

struct GMTK_GM
{
    vector<RandomVariable *> node;
    // This holds all the variables in the graph.
    // The topology os determined by the Parent and Child arrays associated
    // with each random variable.

    vector<RandomVariable *> topologicalOrder;
    // A topological ordering of the nodes; useful for simulation and
    // enumerative inference.

    void verifyTopologicalOrder();
    // verifies that the node array is in topological order
    // copies it to the topological array

    void reveal(vector<RandomVariable *> order, bool show_vals = false);
    // Go through the nodes in the specified order and show them.

    GMTK_GM() {example=NULL;}
    ~GMTK_GM() { delete chain; }

    void makeRandom();
    // Goes over each variable in the graph and calls its makeRandom function.

    void makeUniform();
    // Goes over each variable in the graph and calls its makeUniform function.

    void simulate();
    // This could be called 'sample' since it produces a sample from
    // the set of RVs constituting a GM.
    // Goes over the variables in topological order and instantiates them 
    // according to the probability distribution given their parent's values.

    void enumerateProb(int pos=0, logpr p=1.0);
    // A recursive function that enumerates all possible values of the 
    // hidden variables and computes the total data prob by brute force.
    // Good for verifying DP and other inference results.
    // Call with emMode=false to get the data probability.
    // Call with emMode=true and p=1/dataProb to do EM

    void cliqueChainProb(logpr beam=0.0)
    { chain->computePosteriors(beam); dataProb = chain->dataProb; }
    // Computes the likelihood of the observed variable values with
    // dynamic programming on a clique chain

    bool emMode;
    // when data probability is computed, this controls whether counts
    // are accumulated. If emMode==true, they are accumulated

    void cacheValues();
    // Tells all the nodes int the graph to cache the values they are 
    // currently set to.
    // Observation variables can ignore this.

    void restoreCachedValues();
    // Tells all the nodes in the graph to clamp themselves to their cached
    // values.
    // Observation variables can ignore this.

    void storeValues(vector<VariableValue> &vv);
    // stores the values of all the observed variables

    void setValues(vector<VariableValue> &vv);
    // sets the values of all the observed variables

    void enumerateViterbiProb(int pos=0, logpr p=1.0);
    // Computes the probability of the likeliest instantiation of the hidden
    // variables, and stores it in logViterbiProb. 
    // Has the side effect that at termination, the network is clamped to its
    // likeliest value.

    void cliqueChainViterbiProb(logpr beam=0.0) 
    {chain->doViterbi(beam); viterbiProb=chain->viterbiProb;}
    // Computes the probability of the likeliest instantiation of the hidden
    // variables, and stores it in logViterbiProb. 
    // Has the side effect that at termination, the network is clamped to its
    // likeliest value.
    
    logpr dataProb, viterbiProb;

    void emIncrement(logpr posterior);
    // Tells each variable to increment its counts by the posterior.

    void emStartIteration();
    // Tells each variable to zero out its accumulators and get ready to 
    // start EM.

    void emEndIteration();
    // Tells each variable to update its parameters at the end of an EM 
    // iteration.

    void emSwapCurAndNew();
    // exchanges the newly estimated parameters with the old set, at the
    // end of an em iteration.
    // user has control to do this, or not

    vector<vector<VariableValue> > *example;
    // In EM there must be a way of iterating over the examples.
    // This store the examples to look at.
// need a file-based alternative too

    // how far in the example array have we gotten?
    unsigned expos; 

    void setExampleStream(vector<vector<VariableValue> > *_example)
    {example=_example;}
    // The examples to be iterated over are held in here
// need a file-based analog too

    void clampFirstExample() 
    { if (example==NULL) error("Example array not set.");
      if (example->size()==0) error("No examples.");
      setValues((*example)[0]);
      expos=0;
    }

    // Clamps the observation variables according to the first example.
// need to add a member that specifies if file-based or array-based example
// reading is to be done.

    bool clampNextExample() 
    { if (expos==example->size()-1) return false; 
      setValues((*example)[++expos]); return true; }
    // Clamps the observation variables according to the next example.

    void enumerativeEM(int iterations);
    // Does EM using brute force inference.

    CliqueChain *chain;
    // A pointer to a clique chain representation of the GM.

    void cliqueChainEM(int iterations, logpr beam=0.0);
    // Does EM using dynamic programming on a clique chain.

    logpr enumerativeExampleProb(vector<vector<VariableValue > > &example);
    // computes the likelihood of the examples

    logpr cliqueChainExampleProb(vector<vector<VariableValue > > &example,
          logpr beam=0.0);
    // computes the likelihood of the examples

    void GM2CliqueChain();
    // this creates and initializes a clique chain corresponding to the GM

    void showCliques();
    // shows the cliques in preorder

    void unroll(int first_frame, int last_frame, int times);
    // duplicates the structure from first_frame to last_frame "times" times
    // e.g. unroll(1,2,5) adds 10 frames -- 5 1,2 chunks

    map<pair<RandomVariable *, int>, RandomVariable *> analogue_k_past;
    // To make speedy length adjustments possible, we will keep track of
    // the analogue of every varable that is referenced in the parent lists
    // of the last slice.

    void spliceOut(int segments);
    // splices out "segments" repeating portions of the network.
 
    void restoreNet();
    // undoes the last splice operation
   
    vector<pair<RandomVariable **, RandomVariable *> > parentUpdates;
    // changing the pointers in the first part of the pair to the second part
    // undoes the parent changes made in a splice

    vector<RandomVariable **> parentsToUpdate;
    // which guys to change to do the splice

    vector<RandomVariable *> tempNode, tempTopo;
    // used to store temporary arrays while adjusting length

    int repeatPeriod;   // how many slices is the reepeating period of the net?
  
    int startTimeOfLastSegment;  // last repeating segment
};

#endif
