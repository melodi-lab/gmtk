/*
 * vhash_map.h
 *   General data structure for a vector-based hash-table implementation
 *   of a set. Each element of the set is a vector.
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

#ifndef VHASH_MAP_H
#define VHASH_MAP_H

#include <stdio.h>
#include <stdlib.h>

#include "debug.h"
#include "sArray.h"
#include "hash_abstract.h"

// _Key must be a basic type. This is a hash set containing vectors
// (_Key*).  All elements of the set have the same length, given by
// object constructor.
template <class _Key, class _Data> 
class vhash_map : public hash_abstract {
protected:

#ifdef COLLECT_COLLISION_STATISTICS
public:
#endif

  //////////////////////////////////////////////////////////
  // The size of all the vectors in this set must be the same. Here
  // we store the size once, so they don't need to be stored in each
  // element. Ideally would be const.
  unsigned ksize;

  ////////////////////////////////////////////////////////////
  // the map "buckets", which includes only key and data item.
  struct MBucket {
    _Key* key;
    _Data item;
    MBucket() :key(NULL) {}

    inline bool empty() { return (key == NULL); }

    inline bool keyEqual(const _Key* k2, const _Key* const k2_endp) {
      const _Key* k1 = key;
      do {
	if (*k1++ != *k2++)
	  return false;
      } while (k2 != k2_endp);
      return true;
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

    // printf("resizing table: current elements = %d, current table size =%d, new table size = %d\n",
    // numberUniqueEntriesInserted,table.size(),new_size);

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
	unsigned a = entryOfUnique(table.ptr[i].key,ksize,nt);
	nt.ptr[a].key = table.ptr[i].key;
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
  vhash_map() {}

  ////////////////////
  // constructor
  //    All entries in this hash table have the same size given
  //    by the argument arg_ksize.
  vhash_map(const unsigned arg_ksize,
	     unsigned approximateStartingSize = 
	     hash_abstract::HashTableDefaultApproxStartingSize)
    : ksize(arg_ksize) 
  {
    // make sure we hash at least one element, otherwise do {} while()'s won't work.
    assert (ksize > 0);

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
  inline _Data* insert(_Key* key,
		       _Data val,
		       bool&foundp = hash_abstract::global_foundp) {

#ifdef COLLECT_COLLISION_STATISTICS
    numInserts++;
#endif

    // compute the address
    unsigned a = entryOf(key,ksize,table);

    // printf("inserting to entry %d, empty = %d\n",a,empty(table[a]));
    if (table.ptr[a].empty()) {
      foundp = false;

      table.ptr[a].key = key;
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
	a = entryOf(key,ksize,table);
	assert (!table.ptr[a].empty());
      }
    } else {
      foundp = true;
    }
    return &table.ptr[a].item;
  }


  ///////////////////////////////////////////////////////
  // This version of insert() assumes that the key is unique, meaning
  // that the caller guarantees that the key has never been inserted
  // before.  WARNING: if it is the case that 'key' has been inserted
  // before, the internal state of the hash table will be
  // inconsistent.
  inline _Data* insertUnique(_Key* key,
			     _Data val,
			     bool&foundp = hash_abstract::global_foundp) {

#ifdef COLLECT_COLLISION_STATISTICS
    numInserts++;
#endif

    // compute the address
    unsigned a = entryOfUnique(key,ksize,table);

    if (table.ptr[a].empty()) {
      foundp = false;

      table.ptr[a].key = key;
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
	a = entryOf(key,ksize,table);
	assert (!table.ptr[a].empty());
      }
    } else {
      // we should never find a unique element that we inserted.
      assert ( 0 );
    }
    return &table.ptr[a].item;
  }


  ////////////////////////////////////////////////////////
  // search for key returning true if the key is found, otherwise
  // don't change the table. Return true when found.
  // search for key returning data item if the key is found, otherwise
  // don't change the table. Return pointer to the data item
  // when found, and NULL when not found.
  _Data* find(_Key* key) {
    const unsigned a = entryOf(key,ksize,table);
    if (!table.ptr[a].empty())
      return &(table.ptr[a].item);
    else 
      return NULL;
  }




  ////////////////////////////////////////////////////////
  // Another version of find that returns not only the data item, but
  // also a pointer to the keyp in the hash table itself, which can be
  // modified if need be. If item is not found, key_pp is *NOT*
  // modified.
  _Data* find(_Key* key,_Key**& key_pp) {
    const unsigned a = entryOf(key,ksize,table);
    // printf("find: entry %d\n",a);
    if (!table.ptr[a].empty()){
      key_pp = &table[a].key;
      return &table[a].item;
    } else {
      return NULL;
    }
  }

};


#endif // defined VHASH_MAP

