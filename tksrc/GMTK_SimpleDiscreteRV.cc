
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

#include "GMTK_SimpleDiscreteRV.h"

void GMTK_SimpleDiscreteRV::recMakeRandom(int pos)
{
    if (pos==curConditionalParents.len())  // at the end
    {
        // store the parents values
        vector<int> pvals;
        for (int i=0; i<curConditionalParents.len(); i++)
            pvals.push_back(curConditionalParents[i];

        // make a vector of random self values
        vector<logpr> vals;
        logpr s=0;
        for (int i=0; i<cardinality; i++)
        {
            logpr v = logpr(float(rand())/RAND_MAX);
            s += v;
            vals.push_back(v); 
        }
        for (int i=0; i<cardinality; i++)
            vals[i] /= s;

        // associate the two
        dist_given_parents[pvals] = vals;
        return; 
    }

    RandomVariable *par = curConditionalParents[pos];
    assert(par->discrete);
    for (int i=0; i<par->cardinality; i++)
    {
        par->val = i;
        recMakeRandom(pos+1);
    }
}

