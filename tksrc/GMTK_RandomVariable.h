/* 
 * GMTK_RandomVariable.h
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


#ifndef GMTK_RANDOMVARIABLE
#define GMTK_RANDOMVARIABLE

struct RandomVariable
{
    // Each variable can have some predecessors and successors.
    //
    sArray<RandomVariable *> Parent, Child;

};

#endif
