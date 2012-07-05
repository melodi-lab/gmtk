/*
 * shash_set.h
 *   General data structure for a scalar-based (size 1 key) hash-table
 *   implementation of a set. Each element of the set is a
 *   vector. This hash type is really meant for arrays of words (e.g.,
 *   32 or 64 bit quantities), but could also be used for character
 *   strings.
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2003, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about the suitability of this software
 * for any purpose. It is provided "as is" without express or implied warranty.
 *
 *
 * $Header$
 *
*/

#ifndef SHASH_SET_H
#define SHASH_SET_H

#include <stdio.h>
#include <stdlib.h>

#include "debug.h"
#include "sArray.h"
#include "hash_abstract.h"

// _Key must be a basic type. This is a hash set containing vectors
// (_Key*).  All elements of the set have the same length, given by
// object constructor.
template <class _Key> 
class shash_set : public hash_abstract {

#ifdef COLLECT_COLLISION_STATISTICS
public:
#endif


  ////////////////////////////////////////////////////////////
  // the set "buckets", which includes only key 
  struct SBucket {
    _Key key;
    bool active;
    SBucket() :active(false) { }
    inline bool empty() { return (!active); }
    inline bool keyEqual(const _Key k2) {
      return (key == k2);
    }
  };

  //////////////////////////////////////////////////////////
  // the actual hash table, an array of pointers to T's
  sArray < SBucket > table;

  ////////////////////////////////////////////////////////////
  // resize: resizes the table to be new_size.  
  inline void resize(int new_size) {
    // make a new table (nt),  and re-hash everyone in the
    // old table into the new table.

    // the next table, used for table resizing.
    sArray < SBucket > nt;

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
  shash_set() {}

  ////////////////////
  // constructor
  //    All entries in this hash table have the same size (= 1).
  shash_set(unsigned approximateStartingSize = 
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
  inline _Key* insert(_Key key,
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
    return &(table.ptr[a].key);
  }


  inline _Key* insertUnique(_Key key,
			    bool&foundp = hash_abstract::global_foundp)
  {
    return insert(key,foundp);
  }


  ////////////////////////////////////////////////////////
  // search for key returning true if the key is found, otherwise
  // don't change the table. Return true when found.
  inline bool find(_Key key) {
    const unsigned a = entryOf(key,table);
    return !(table.ptr[a].empty());
  }

  ////////////////////////////////////////////////////////
  // another version of find that also returns a pointer to the key in
  // the hash table itself. Returns true if the key is contained in
  // the set, false otherwise.
  inline bool find(_Key key,_Key*& key_pp) {
    const unsigned a = entryOf(key,table);
    if (table.ptr[a].empty()) {
      key_pp = NULL;
      return false;
    } else {
      key_pp = &(table.ptr[a].key);
      return true;
    }
  }

};


#endif // defined SHASH_SET

