
/* 
 * GMTK_Clique.cc
 * The basic Clique data structure.
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

#include "general.h"
VCID("$Header$");

#include "GMTK_Clique.h"

set<vector<RandomVariable::DiscreteVariableType>, VecCompare > 
CliqueValue::global_val_set;

vector<CliqueValue> Clique::gip;  // the actual global instantiation pool
vector<unsigned> Clique::freelist;     // and the freelist
int Clique::nextfree=-1;

static const float mem_factor=1.2;  // resize by this factor

unsigned Clique::newCliqueValue()
{
    if (nextfree==-1)         // nothing more left
    {
        freelist.clear();     // the previous entries have all been used
        int oldsize = gip.size(), newsize = int(mem_factor*gip.size());
        newsize = max(newsize, 100000);  // don't mess around at the beginning
        gip.resize(newsize);  // make more CliqueValues
        for (int i=oldsize; i<newsize; i++)
            freelist.push_back(i);  // add the indexes of the new CliqueValues
        nextfree = freelist.size()-1;  // initialize the nextfree pointer
    }
    return freelist[nextfree--];
}

void Clique::recycleCliqueValue(unsigned idx)
{
    // resizing the gip changes its size. 
    // after a resize, the freelist's size is equal to the number of 
    // CliqueValues added; not the size of the gip. Therefore, it may be 
    // necessary to make it bigger, e.g. when absolutely everything is freed up
    if (nextfree == int(freelist.size())-1)  
        freelist.resize(int(mem_factor*freelist.size()));

    freelist[++nextfree] = idx;
}

void Clique::cacheClampedValues()
{
    for (unsigned i=0; i<discreteMember.size(); i++)
        clampedValues[i] = discreteMember[i]->val;
}

logpr Clique::probGivenParents()
{
    logpr p = 1.0;
    findConditionalProbabilityNodes();
    for (unsigned i=0; i<conditionalProbabilityNode.size(); i++)
        p *= conditionalProbabilityNode[i]->probGivenParents();
    return p;
}

/*-
 *-----------------------------------------------------------------------
 * Function
 *      prune: prune away all CliqueValues that lie outside of the beam.
 *
 * Results:
 *      nothing
 *
 * Side Effects:
 *      The low probability entries are removed from the instantiation list
 *      Their memory is recycled and placed back on the freelist
 *
 *-----------------------------------------------------------------------
 */
void Clique::prune(logpr beam)
{
    // find the maximum probability instantiation
    logpr maxv = 0.0;
    for (unsigned i=0; i<instantiation.size(); i++)
        if (gip[instantiation[i]].pi > maxv)
            maxv = gip[instantiation[i]].pi;

    // swap back the below-threshold entries
    logpr threshold = maxv*beam;
    int i=0, j=instantiation.size()-1;
    while (i <= j)
        if (gip[instantiation[i]].pi < threshold)
            swap(instantiation[i], instantiation[j--]);
        else
            i++;
        
    // delete the low-probability guys  -- from i to the end
    for (unsigned k=i; k<instantiation.size(); k++)
        recycleCliqueValue(instantiation[k]);
    instantiation.erase(&instantiation[i], instantiation.end()); 
}

void Clique::enumerateValues(int new_member_num, int pred_val, bool viterbi)
{
    if (separator)
    {
        cacheClampedValues();
     
        // Make sure we have a clique value to work with
        // instantiationAddress tells if the instantiation was seen before
        map<vector<RandomVariable::DiscreteVariableType>, int>::iterator mi;
        CliqueValue *cv;
        if ((mi=instantiationAddress.find(clampedValues)) == 
        instantiationAddress.end())               // not seen before
        {
            int ncv = newCliqueValue();
            cv = &gip[ncv];
            instantiation.push_back(ncv);                
            instantiationAddress[clampedValues] = ncv;  
            // instead of storing the underlying variable values, set the pred
            // pointer, and use the previous non-separators values --
            // they are guaranteed to be consistent
            cv->values = NULL;  // will use the parent's 
            cv->lambda = cv->pi = 0.0;
            cv->pred = pred_val;
        }
        else
            cv = &gip[(*mi).second];              // will word with old value

        if (!viterbi)
        {
	    // accumulate in probability
            cv->pi += gip[pred_val].pi;
            gip[pred_val].succ = instantiationAddress[clampedValues];
        }
        else if (gip[pred_val].pi >= cv->pi)
        {
	    // replace value since it is greater
            cv->pi = gip[pred_val].pi;
            cv->pred = pred_val;
        }
    }
    else if (new_member_num == int(newMember.size())) 
    // base case: all members fixed
    {
	// Then all members of this clique have their values clamped,
	// and we are ready to compute the probablity of this clique.
	// Each variable and its parents are guaranteed to be in this
	// clique and all such variables are clamped, so making this
	// possible.

        logpr pi = probGivenParents();
	// copy in the clique value -- if it has a nonzero probability
        // otherwise, discard it to avoid further propagation
        if (pi != 0.0)    
        {
            cacheClampedValues();
            int ncv = newCliqueValue();
            instantiation.push_back(ncv);
            CliqueValue *cv = &gip[ncv];
            cv->pi = cv->lambda = pi;  // cache value in lambda
            cv->pred = pred_val;
            set<vector<RandomVariable::DiscreteVariableType>,
                VecCompare>::iterator si;
            si = CliqueValue::global_val_set.insert(
                CliqueValue::global_val_set.begin(), clampedValues);
            cv->values = &(*si);
            if (pred_val!=-1)   // not doing root
                cv->pi *= gip[pred_val].pi;
        }
    }
    else
    {
	// loop over all values of the current variable number,
	// and recurse.
        newMember[new_member_num]->clampFirstValue();
        do
        {
            enumerateValues(new_member_num+1, pred_val, viterbi);
        } while (newMember[new_member_num]->clampNextValue());
    }
}

void Clique::reveal()
{
    for (unsigned i=0; i<member.size(); i++)
    {
        cout << member[i]->label << "-" << member[i]->timeIndex << " ";
        for (unsigned j=0; j<newMember.size(); j++)
            if (newMember[j] == member[i])
                cout << "(new) ";
    }
    cout << endl;
}
