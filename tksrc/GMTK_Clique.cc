
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

/*-
 *-----------------------------------------------------------------------
 * Function
 *      findConditionalProbabilityNodes, 
 *
 * Results:
 *      What does the function return, if anything. What might it do?
 *
 * Side Effects:
 *      Any exernal side effects
 *
 *-----------------------------------------------------------------------
 */


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
void GMTK_Clique::prune(logpr beam)
{
    // Can't prune separators because their clique entries have pointers
    // pointing to them => can't move.
    if (separator)   
        return;

    logpr maxv = 0;
    map<vector<DISCRETE_VARIABLE_TYPE>, CliqueValue>::iterator mi;
    for (mi=instantiation.begin(); mi!=instantiation.end(); mi++)
        if ((*mi).second.pi > maxv)
            maxv = (*mi).second.pi;
        
    logpr threshold = maxv*beam;
    map<vector<DISCRETE_VARIABLE_TYPE>, CliqueValue> pruned;
    for (mi=instantiation.begin(); mi!=instantiation.end(); mi++)
        if ((*mi).second.pi >= threshold)
            pruned.insert(*mi);

    instantiation = pruned;
}

void GMTK_Clique::enumerateValues(int new_member_num, CliqueValue *pred_val,
bool viterbi)
{
    if (separator)
    {
        cacheClampedValues();
     
        // make sure something is in the table.
	// this occurs here for separators only because all
	// values of all variables are clamped at this point (since
	// a separator is a subset of the previous non-separator clique)
        if (instantiation.find(campedValues) == instantiation.end())
        {
            CliqueValue cv;
            cv.pi = 0;
            instantiation[clampedValues] = cv;
        }
        
        if (!viterbi)
        {
	    // accumulate in probability
            instantiation[clampedValues].pi += pred_val->pi;
// We are relying on the address of instantiation[clampedValues] not to
// change as new entries are added to the map
            pred_val->succ = &instantiation[clampedValues];
        }
        else if (pred_val->pi >= instantiation[clampedValues].pi)
        {
	    // replace value since it is greater
            instantiation[clampedValues].pi = pred_val->pi;
            instantiation[clampedValues].pred = pred_val;
        }
    }
    else if (new_member_num == newMember.len())  // base case: all members fixed
    {
	// Then all members of this clique have their values clamped,
	// and we are ready to compute the probablity of this clique.
	// Each variable and its parents are guaranteed to be in this
	// clique and all such variables are clamped, so making this
	// possible.
        cacheClampedValues();

        CliqueValue cv;
        cv.pi = cv.lambda = probGivenParents();  // cache value in lambda
        cv.pred = pred_val;
        if (pred_val)   // not doing root
            cv.pi *= pred_val->pi;
	// copy in the clique value -- if it has a nonzero probability
        // otherwise, discard it to avoid further propagation
        if (cv.pi != 0)
// does logpr know what 0 is??
            instantiation[clampedValues] = cv;
    }
    else
    {
	// loop over all values of the current variable number,
	// and recurse.
        newMember[new_member_num]->clampFirstValue();
        do
        {
            enumerateValues(new_member_num+1, pred_val, viterbi);
        } while (member[new_member_num]->clampNextValue();
    }
}

