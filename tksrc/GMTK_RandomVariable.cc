
/* 
 * GMTK_RandomVariable.cc
 * Defines the functions that all random variables must satisfy in order
 * for the inference, sampling, and other generic routines to work.
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


/*-
 *-----------------------------------------------------------------------
 *  setParents
 *     sets up the switching and conditional parents.
 *
 * Results:
 *
 * Side Effects:
 *     sets the switchingParents array
 *     sets the conditionalParentsList array
 *     adds to the allPossibleChildren array
 *
 *-----------------------------------------------------------------------
 */

#include <set>

void GMTK_RandomVariable::setParents(sArray<RandomVariable *> &sparents,
         sArray<sArray<RandomVariable *> > &cpl)
{
    switchingParents = sparents;
    conditionalParentsList = cpl;

    // add this to the child list of each of its parents
    // note that a parent may occur twice, e.g. as a switching and conditional
    // parent, or on twqo different conditional parent lists
    // avoid adding this to the child list twice
    set<RandomVariable> parents;  // the set of parents this is a child of
    for (int i=0; i<sparents.len(); i++)
    {
         parents.insert(sparents[i]);
         sparents[i]->allPossibleChildren.push(this);
    }
    
    for (int i=0; i<cpl.len(); i++)
        for (int j=0; j<cpl[i].len(); j++)
            if (parents.find(cpl[i].[j]) == parents.end())
            {
                cpl[i].[j]->allPossibleChildren.push(this);
                parents.push(cpl[i].[j]);
            }
}
