
/* 
 * GMTK_SimpleDiscreteRV.h
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

#include <numeric>
#include "GMTK_SimpleDiscreteRV.h"

#include "rand.h"

void SimpleDiscreteRV::recMakeRandom(int pos)
{
    if (unsigned(pos)==(*curConditionalParents).size())  // at the end
    {
        // store the parents values
        vector<int> pvals;
        for (unsigned i=0; i<(*curConditionalParents).size(); i++)
            pvals.push_back((*curConditionalParents)[i]->val);

        // make a vector of random self values
        vector<logpr> vals;
        logpr s=0.0;
        for (int i=0; i<cardinality; i++)
        {
	    logpr v = rnd.drand48();
            s += v;
            vals.push_back(v); 
        }
        for (int i=0; i<cardinality; i++)
            vals[i] /= s;

        // associate the two
        (*dist_given_parents)[pvals] = vals;
        return; 
    }

    RandomVariable *par = (*curConditionalParents)[pos];
    assert(par->discrete);
    for (int i=0; i<par->cardinality; i++)
    {
        par->val = i;
        recMakeRandom(pos+1);
    }
}

void SimpleDiscreteRV::instantiate()
{
    findConditionalParents();
    // store the parents values
    vector<int> pvals;
    for (unsigned i=0; i<(*curConditionalParents).size(); i++)
        pvals.push_back((*curConditionalParents)[i]->val);
    vector<logpr> &vals = (*dist_given_parents)[pvals];

    logpr target = 0.9999*logpr(float(rand())/RAND_MAX);
    int pos = 0;
    logpr sum = vals[0]; 
    while (sum < target)
        sum += vals[++pos];
    assert(pos<cardinality);
    val = pos;
}
   
void SimpleDiscreteRV::emStartIteration()
{
    // assume that the CPTs have been set up by now
    assert(dist_given_parents->size());        

    *counts_given_parents = *dist_given_parents;  // get the sizes right
    map<vector<int>, vector<logpr> >::iterator mi;
    for (mi=(*counts_given_parents).begin(); mi!=(*counts_given_parents).end();
    mi++)
    {
        vector<logpr> &v = (*mi).second;
        fill(v.begin(), v.end(), logpr(0.0));  // zero out the counts
    }
}

void SimpleDiscreteRV::emEndIteration()
{
    // inefficient as hell
    *dist_given_parents = *counts_given_parents;
    map<vector<int>, vector<logpr> >::iterator mi;
    for (mi=(*dist_given_parents).begin();mi!=(*dist_given_parents).end();mi++)
    {
        vector<logpr> &v = (*mi).second;
        logpr sum = accumulate(v.begin(), v.end(), logpr(0.0));
        if (sum != 0.0)
            for (unsigned i=0; i<v.size(); i++) v[i] /= sum;
        else
            for (unsigned i=0; i<v.size(); i++) v[i] = 1.0/v.size();
    }
}

void SimpleDiscreteRV::emIncrement(logpr posterior)
{
    findConditionalParents();
    // store the parents values
    vector<int> pvals;
    for (unsigned i=0; i<(*curConditionalParents).size(); i++)
        pvals.push_back((*curConditionalParents)[i]->val);
    (*counts_given_parents)[pvals][val] += posterior;
}

