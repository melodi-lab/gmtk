
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

#include <set>
#include "GMTK_RandomVariable.h"

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
 *     sets the allPossibleParents Array
 *     adds to the allPossibleChildren array of the variable's parents
 *
 *-----------------------------------------------------------------------
 */

void RandomVariable::setParents(vector<RandomVariable *> &sparents,
         vector<vector<RandomVariable *> > &cpl)
{
    switchingParents = sparents;
    conditionalParentsList = cpl;

    // add this to the child list of each of its parents
    // note that a parent may occur twice, e.g. as a switching and conditional
    // parent, or on twqo different conditional parent lists
    // avoid adding this to the child list twice
    set<RandomVariable *> parents;  // the set of parents this is a child of
    for (unsigned i=0; i<sparents.size(); i++)
    {
         parents.insert(sparents[i]);
         sparents[i]->allPossibleChildren.push_back(this);
    }
    
    for (unsigned i=0; i<cpl.size(); i++)
        for (unsigned j=0; j<cpl[i].size(); j++)
            if (parents.find(cpl[i][j]) == parents.end())
            {
                cpl[i][j]->allPossibleChildren.push_back(this);
                parents.insert(cpl[i][j]);
            }

    // set up the all possible parents array.
    allPossibleParents.resize(parents.size());
    set<RandomVariable *>::iterator si;
    int p=0;
    for (si = parents.begin(); si != parents.end(); si++)
        allPossibleParents[p++] = *si;
}

void RandomVariable::reveal(bool show_vals)
{
    cout << label << "-" << timeIndex << " : ";
    if (discrete) cout << "discrete (" << cardinality << ") ";
    if (hidden) cout << "hidden "; else cout << "observed ";
    if (!hidden && discrete || show_vals) cout << "val (" << val << ")";

    cout << " possible parents: ";
    for (unsigned i=0; i<allPossibleParents.size(); i++)
        cout << allPossibleParents[i]->label << " ";
    cout << endl;
}
