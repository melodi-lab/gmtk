/* 
 * GMTK_DiscreteRandomVariable.h
 * A specific case of the RandomVariable class
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

#ifndef GMTK_DISCRETERANDOMVARIABLE
#define GMTK_DISCRETERANDOMVARIABLE

#include "GMTK_RandomVariable.h"
#include "GMTK_CPT.h"

struct DiscreteRandomVariable : public RandomVariable
{
    logpr probGivenParents();

    void clampFirstValue();

    bool clampNextValue();

    void makeRandom();

    void makeUniform();

    void instantiate();

    void tieWith(randomVariable *rv);

    void zeroAccumulators();

    void increment(logpr posterior);

    void update();

    DISCRETE_VARIABLE_TYPE cached_val;
    void cacheValue() = {cached_val=val;}

    void restoreCachedValue() = {val=cached_val;}
};

#endif
