
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

void GMTK_GM::findTopologicalOrder(randomVariabe *rv)
{

    /* Finds a topological order of the random variables.
       It is assumed that the union of all possible switching dependencies
       are explicitly listed in the parent arrays of the random variables.
       Not that the switching parent(s) as well as actual conditioning 
       parents must be explicitly listed.
       The first call is made with no parameters.
    */

    static 

}
