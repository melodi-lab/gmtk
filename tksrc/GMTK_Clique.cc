
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

#include "GMTK_Clique.h"

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
 *      Totally wipes out the previous instantiation object and 
 *      replaces it with the sparse pruned one.
 *
 *-----------------------------------------------------------------------
 */
void Clique::prune(logpr beam)
{
    // find the maximum probability instantiation
    logpr maxv = 0.0;
    list<CliqueValue>::iterator li;
    for (li=instantiation.begin(); li!=instantiation.end(); li++)
        if (li->pi > maxv)
            maxv = li->pi;
        
    // make a list of all the low probability entries
    logpr threshold = maxv*beam;
    vector<list<CliqueValue>::iterator> kill;
    for (li=instantiation.begin(); li!=instantiation.end(); li++)
        if (li->pi < threshold)
            kill.push_back(li);

    // delete them
    for (unsigned i=0; i<kill.size(); i++)
        instantiation.erase(kill[i]);
}

void Clique::enumerateValues(int new_member_num, CliqueValue *pred_val,
bool viterbi)
{
    if (separator)
    {
        cacheClampedValues();
     
        // Make sure we have a clique value to work with
        // instantiationAddress tells if the instantiation was seen before
        map<vector<RandomVariable::DiscreteVariableType>, CliqueValue *>::iterator mi;
        CliqueValue *cv;
        if ((mi=instantiationAddress.find(clampedValues)) == 
        instantiationAddress.end())                  // not seen before
        {
            instantiation.push_back(CliqueValue());  // auto-init pi to 0.0
            cv = &instantiation.back();              // work with new value
            instantiationAddress[clampedValues] = cv;   // store
            // store the underlying variable values with the instantiation
            cv->values = clampedValues;   
        }
        else
            cv = (*mi).second;                    // will word with old value

        if (!viterbi)
        {
	    // accumulate in probability
            cv->pi += pred_val->pi;
            pred_val->succ = cv;
        }
        else if (pred_val->pi >= cv->pi)
        {
	    // replace value since it is greater
            cv->pi = pred_val->pi;
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
            instantiation.push_back(CliqueValue());
            CliqueValue *cv = &instantiation.back();
            cv->pi = cv->lambda = pi;  // cache value in lambda
            cv->pred = pred_val;
            cv->values = clampedValues;
            if (pred_val)   // not doing root
                cv->pi *= pred_val->pi;
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
