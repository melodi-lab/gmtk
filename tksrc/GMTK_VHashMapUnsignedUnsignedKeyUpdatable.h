/*
 * GMTK_MAXCLIQUE.h
 *
 * a special vhash class, for mapping from keys consisting of
 * compressed sets of RV values, to items consisting of array indices.
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (C) 2001 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

#ifndef GMTK_VHASHMAPUNSIGNEDUNSIGNEDKEYUPDATABLE_H
#define GMTK_VHASHMAPUNSIGNEDUNSIGNEDKEYUPDATABLE_H

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include "general.h"
#include "vhash_map.h"
#include "debug.h"

typedef vhash_map < unsigned, unsigned > VHashMapUnsignedUnsigned;
class VHashMapUnsignedUnsignedKeyUpdatable : public VHashMapUnsignedUnsigned {
public:
  //////////////////////
  // constructor for empty invalid object. 
  // WARNING: this will create an invalid object. It is assumed
  // that this object will re-reconstructed later.
  VHashMapUnsignedUnsignedKeyUpdatable() {}
  // constructor
  VHashMapUnsignedUnsignedKeyUpdatable(const unsigned arg_vsize,
				       unsigned approximateStartingSize = 
				       HashTableDefaultApproxStartingSize);

  //
  // Direct access to tables and keys (made available for speed).  Use
  // sparingly, and only if you know what you are doing about the
  // internals of a hash table. Note also that this breaks
  // encapsulation, meaning that if the implementation of the
  // internals of the parent hash table change, this code might break.
  unsigned*& tableKey(const unsigned i) { return table.ptr[i].key; }
  unsigned& tableItem(const unsigned i) { return table.ptr[i].item; }
  bool tableEmpty(const unsigned i) { return table.ptr[i].empty(); }
  unsigned tableSize() { return table.size(); }

};

#endif
