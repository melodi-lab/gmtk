
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

void GMTK_Clique::cacheClampedValues()
{
    for (int i=0; i<discreteMember.len(); i++)
        clampedValues[i] = discreteMember[i]->val;
}

logpr GMTK_Clique::probGivenParents()
{
    logpr p = 1;
    findConditionalProbabilityNodes();
    for (int i=0; i<conditionalProbabilityNode.len(); i++)
        p *= conditionalProbabilityNode[i]->probGivenParents();
    return p;
}

void GMTK_Clique::enumerateValues(int new_member_num, cliqueValue *pred_val,
bool viterbi)
{
    if (separator)
    {
        cacheClampedValues();
     
        // make sure something is in the table.
        if (instantiation.find(campedValues) == instantiation.end())
        {
            cliqueValue cv;
            cv.pi = 0;
            instantiation[clampedValues] = cv;
        }
        
        if (!viterbi)
        {
            instantiation[clampedValues].pi += pred_val->pi;
            pred_val->succ = &instantiation[clampedValues];
        }
        else if (pred_val->pi >= instantiation[clampedValues].pi)
        {
            instantiation[clampedValues].pi = pred_val->pi;
            instantiation[clampedValues].pred = pred_val;
        }
    }
    else if (new_member_num == newMember.len())  // base case: all members fixed
    {
        cacheClampedValues();
        cliqueValue cv;
        cv.pi = cv.lambda = probGivenParents();  // cache value in lambda
        cv.pred = pred_val;
        if (pred_val)   // not doing root
            cv.pi *= pred_val->pi;
        instantiation[clampedValues] = cv;
    }
    else
    {
        newMember[new_member_num]->clampFirstValue();
        do
        {
            enumerateValues(new_member_num+1, pred_val, viterbi);
        } while (member[new_member_num]->clampNextValue();
}

