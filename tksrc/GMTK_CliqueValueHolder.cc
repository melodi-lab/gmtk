/*-
 * GMTK_CliqueValueHolder.cc
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (C) 2003 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 *
 */



/*
 * TODO:
 *   1) Maybe: In some of the structures the choice of when
 *   to use a hash table vs. an inlined value is done
 *   by either using an unsigned* or a unsigned. Create
 *   options where we do 'unsigned val[LEN]' and where
 *   only if packed value is > LEN words do we resort to
 *   the hash table.
 *
 *   2) Figure out way to remove the (unsigned**) type casts
 *   in front of call to packer object (see below).   
 *
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif
#if HAVE_HG_H
#include "hgstamp.h"
#endif

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <float.h>
#include <assert.h>
#include <ctype.h>

#include <iterator>
#include <map>
#include <set>
#include <algorithm>

#include "general.h"
#include "error.h"
#include "debug.h"

#include "GMTK_CliqueValueHolder.h"

VCID(HGID)



float CliqueValueHolder::defaultGrowthFactor;
unsigned CliqueValueHolder::defaultAllocationUnitChunkSize;


CliqueValueHolder::CliqueValueHolder(unsigned _cliqueValueSize) 
  : cliqueValueSize(_cliqueValueSize),
    growthFactor(defaultGrowthFactor),
    allocationUnitChunkSize(defaultAllocationUnitChunkSize)

{
  assert ( allocationUnitChunkSize > 0 );
  prepare();
}


CliqueValueHolder::CliqueValueHolder(unsigned _cliqueValueSize,
				     unsigned _allocationUnitChunkSize,
				     float _growthFactor)
  : cliqueValueSize(_cliqueValueSize),
    growthFactor(_growthFactor),
    allocationUnitChunkSize(_allocationUnitChunkSize)
{
  assert ( allocationUnitChunkSize > 0 );
  prepare();
}


// make sure there is rom for one element.
void
CliqueValueHolder::prepare()
{
  makeEmpty();
  values.resize(1);
  // newSize *MUST* be a multiple of 'cliqueValueSize' or else
  // this code will fail.
  unsigned newSize = cliqueValueSize*allocationUnitChunkSize;
  values[values.size()-1].resize(newSize);
  capacity = newSize;
  numAllocated = 0;
  curAllocationPosition = values[values.size()-1].ptr;
  curAllocationEnd = values[values.size()-1].ptr + newSize;

}


// free up all memory, postcondition: make an invalid object until
// next prepare is called.
void
CliqueValueHolder::makeEmpty()
{
  for (unsigned i=0;i<values.size();i++) {
    values[i].clear();
  }
  values.clear();
  capacity = 0;
  numAllocated = 0;
}



void
CliqueValueHolder::allocateCurCliqueValue()
{

  curAllocationPosition += cliqueValueSize;
  numAllocated++;

  // first to a cheap and fast allocation of new clique value storage
  if (curAllocationPosition != curAllocationEnd) {
    return;
  } 

  // if here, we need to allocate another chunk add a new chunk so we
  // don't need to re-copy all the existing ones already.
  values.resizeAndCopy(values.size()+1);

  // newSize *MUST* be a multiple of 'cliqueValueSize' or else.
  // this code will fail.
  // TODO: optimize this re-sizing.
  unsigned newSize = cliqueValueSize*
    unsigned(1+allocationUnitChunkSize*
	     ::pow(growthFactor,values.size()-1));

  values[values.size()-1].resize(newSize);

  capacity += newSize;
  curAllocationPosition = values[values.size()-1].ptr;
  curAllocationEnd = values[values.size()-1].ptr + newSize;

}


/////////////////////////////////////////////	
/// END OF FILE
/////////////////////////////////////////////
