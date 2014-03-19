/*
 * shash_map.h
 *   General data structure for a scalar-based (size 1 key) hash-table
 *   implementation of a set. Each element of the set is a
 *   vector. This hash type is really meant for arrays of words (e.g.,
 *   32 or 64 bit quantities), but could also be used for character
 *   strings.
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (C) 2003 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 *
 *
 * $Header$
 *
*/

#ifndef SHASH_MAP_H
#define SHASH_MAP_H

#include <stdio.h>
#include <stdlib.h>

#include "debug.h"
#include "sArray.h"
#include "hash_abstract.h"

// _Key must be 32 bits.
template <class _Key, class _Data> 
class shash_map : public hash_abstract {
protected:

#ifdef COLLECT_COLLISION_STATISTICS
public:
#endif


  ////////////////////////////////////////////////////////////
  // the set "buckets", which includes only key 
  struct MBucket {
    _Key key;
    _Data item;
    bool active;
    MBucket() :active(false) { }
    inline bool empty() { return (!active); }
    inline bool keyEqual(const _Key k2) {
      return (key == k2);
    }
  };

  //////////////////////////////////////////////////////////
  // the actual hash table, an array of pointers to T's
  sArray < MBucket > table;

  ////////////////////////////////////////////////////////////
  // resize: resizes the table to be new_size.  
  inline void resize(int new_size) {
    // make a new table (nt),  and re-hash everyone in the
    // old table into the new table.

    // the next table, used for table resizing.
    sArray < MBucket > nt;

#if defined(HASH_PRIME_SIZE)
    nt.resize(new_size);
#else
    // In this case, new_size *MUST* be a power of two, or everything
    // will fail.
    nt.resize(new_size);
    // make sure to re-set mask before a rehash.
    sizeMask = new_size-1;
#ifdef HASH_LOC_FOLD
    logSize = ceilLog2(new_size);
#endif
#endif

    for (unsigned i=0;i<table.size();i++) {
      if (!table.ptr[i].empty()) {
	// assume that the key is unique, so collision detection
	// can be done just by a bucket being non-empty.
	unsigned a = entryOfUnique(table.ptr[i].key,nt);
	nt.ptr[a].key = table.ptr[i].key;
	nt.ptr[a].active = true;
	nt.ptr[a].item = table.ptr[i].item;
      }
    }

    // iterate over only non-empty entries in hash table.
    table.swap(nt);
    // clear out nt which is now old table.
    nt.clear();

    numEntriesToCauseResize = (int)(loadFactor*(float)new_size);

  }

public:

  // Dummy constructor to create an invalid object to be re-constructed
  // later.
  // This was commented because it will conflict with the following
  // constructor with default size
  //shash_map() {}

  ////////////////////
  // constructor
  //    All entries in this hash table have the same size (= 1).
  shash_map(unsigned approximateStartingSize = 
	     hash_abstract::HashTableDefaultApproxStartingSize)
  {
    // make sure we hash at least one element, otherwise do {} while()'s won't work.

    numberUniqueEntriesInserted=0;
#if defined(HASH_PRIME_SIZE)
    findPrimesArrayIndex(approximateStartingSize);
    initialPrimesArrayIndex = primesArrayIndex;
    // create the actual hash table here.
    resize(HashTable_PrimesArray[primesArrayIndex]);
#else
    // create the actual hash table here.
    resize(nextPower2(approximateStartingSize));
#endif

#ifdef COLLECT_COLLISION_STATISTICS
    maxCollisions = 0;
    numCollisions = 0;
    numInserts = 0;
#endif

  }


  /////////////////////////////////////////////////////////
  // clear out the table entirely, including deleting
  // all memory pointed to by the T* pointers. 
  void clear(unsigned approximateStartingSize = 
	     hash_abstract::HashTableDefaultApproxStartingSize) 
  {
    table.clear();
    numberUniqueEntriesInserted=0;
#if defined(HASH_PRIME_SIZE)
    primesArrayIndex=initialPrimesArrayIndex;
    resize(HashTable_PrimesArray[primesArrayIndex]);
#else
    resize(nextPower2(approximateStartingSize));
#endif
  }

    //Do not clear the data, only reset the structure
    void reset(unsigned approximateStartingSize = 
	     hash_abstract::HashTableDefaultApproxStartingSize) {

        table.reset(1);
        numberUniqueEntriesInserted=0;
        #if defined(HASH_PRIME_SIZE)
            findPrimesArrayIndex(approximateStartingSize);
            initialPrimesArrayIndex = primesArrayIndex;
            resize(HashTable_PrimesArray[primesArrayIndex]);
        #else
            resize(nextPower2(approximateStartingSize));
        #endif
    }


#ifdef COLLECT_COLLISION_STATISTICS
  // clear out the statistics
  void clearStats() {
    maxCollisions = numCollisions = numInserts = 0;
  }
#endif

  ///////////////////////////////////////////////////////
  // insert an item <key> into the hash table.  Return a pointer to
  // the key in the hash table after it has been inserted.  The foundp
  // argument is set to true when the key has been found. 
  inline _Data* insert(_Key key, _Data val,
		       bool&foundp = hash_abstract::global_foundp) {

#ifdef COLLECT_COLLISION_STATISTICS
    numInserts++;
#endif



    // compute the address
    unsigned a = entryOf(key,table);

    if (table.ptr[a].empty()) {
      foundp = false;

      table.ptr[a].key = key;
      table.ptr[a].active = true;
      table.ptr[a].item = val;

      // time to resize if getting too big.
      if (++numberUniqueEntriesInserted >= numEntriesToCauseResize) {
#if defined(HASH_PRIME_SIZE)
	if (primesArrayIndex == (HashTable_SizePrimesArray-1)) 
	  error("ERROR: Hash table error, table size exceeds max size of %lu",
		HashTable_PrimesArray[primesArrayIndex]);
	resize(HashTable_PrimesArray[++primesArrayIndex]);	
#else
	resize(nextPower2(table.size()+1));
#endif

	// need to re-get location
	a = entryOf(key,table);
	assert (!table.ptr[a].empty());
      }
    } else {
      foundp = true;
    }
    return &(table.ptr[a].item);
  }

  ////////////////////////////////////////////////////////
  // remove an item according to key
  inline void remove(_Key key, bool&foundp = hash_abstract::global_foundp) {
    // compute the address
    unsigned a = entryOf(key,table);
    if ( table.ptr[a].empty() ) {
      foundp = false;
    } else {
      foundp = true;
      table.ptr[a].active = false;
    }
  }


  inline _Data* insertUnique(_Key key, _Data val,
			    bool&foundp = hash_abstract::global_foundp)
  {
    return insert(key,val,foundp);
  }


  ////////////////////////////////////////////////////////
  // search for key returning true if the key is found, otherwise
  // don't change the table. Return true when found.
  inline _Data* find(_Key key) {
    const unsigned a = entryOf(key,table);
    if (!table.ptr[a].empty())
      return &table.ptr[a].item;
    else
      return NULL;
  }

};


#endif // defined SHASH_MAP

