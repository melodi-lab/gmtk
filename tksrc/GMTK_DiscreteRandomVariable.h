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

struct RandomVariable
{
    virtual logpr probGivenParents() = 0;

    virtual void clampFirstValue() = 0;

    virtual bool clampNextValue() = 0;

    virtual void makeRandom() = 0;

    virtual void makeUniform() = 0;

    virtual void instantiate() = 0;

    virtual void tieWith(randomVariable *rv) = 0;

    virtual void zeroAccumulators() = 0;

    virtual void increment(logpr posterior) = 0;

    virtual void update() = 0;

    short cached_val;
    virtual void cacheValue() = {cached_val=val;}

    virtual void restoreCachedValue() = {val=cached_val;}
};

#endif
