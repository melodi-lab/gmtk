/*
 * shash_set.h
 *   General data structure for a scalar-based hash-table implementation
 *   of a set. Each element of the set is a scalar.
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 * Modified by Gang Ji <gang@ee.washington.edu>
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

// _Key must be a basic type.
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
    SBucket() : active(false) {}
  };

  //////////////////////////////////////////////////////////
  // the actual hash table, an array of pointers to T's
  sArray < SBucket > table;

#ifdef COLLECT_COLLISION_STATISTICS
  unsigned maxCollisions;
  unsigned numCollisions;
  unsigned numInserts;
#endif

  ///////////////////////////////////////////////////////////////////////
  // since this is a double hash table, we define two address
  // functions, h1() and h2().  h1() gives the starting position of in
  // the array of the key, and h2() gives the increment when we have a
  // collision.
  unsigned h1(const _Key& key, const unsigned arg_size) {
    return (3367900314ul + key) % arg_size;
  }

  ///////////////////////////////////////////////////////////////////////
  // h2() is the increment of of key when we have a collision.  "The
  // value of h2(key) must be relatively prime to the hash-table m for
  // the entire hash table to be searched" (from Corman, Leiserson,
  // Rivest). Therefore, we have result of the h2() function satisfy
  //         1) it must be greater than zero
  //         2) it must be strictly less than table.size()
  //
  unsigned h2(const _Key& key, const unsigned arg_size) {
    return (key + 1) % (arg_size - 1) + 1;
  }

  ///////////////////////////////////////
  // return true if the bucket is empty
  bool empty(const SBucket* bucket) {
    return bucket->active == false;
  }
  // return true if the bucket is empty
  bool empty(const SBucket& bucket) {
    return bucket.active == false;
  }

  //////////////////////////////////////////////////////////////////
  // return the entry of key in table a_table
  unsigned entryOf(const _Key& key, sArray<SBucket> & a_table) {
    const unsigned size = a_table.size();
    unsigned a = h1(key, size);

#ifdef COLLECT_COLLISION_STATISTICS
    unsigned collisions=0;
#endif

    if ( (! a_table[a].active) || a_table[a].key == key ) {
      return a;
    }
    const unsigned inc = h2(key, size);
    do {
#ifdef COLLECT_COLLISION_STATISTICS
      collisions++;
#endif
      a = (a+inc) % size;
    } while ( a_table[a].active && a_table[a].key != key );

#ifdef COLLECT_COLLISION_STATISTICS
    if (collisions > maxCollisions)
      maxCollisions = collisions;
    numCollisions += collisions;
#endif

    return a;
  }


  ////////////////////////////////////////////////////////////
  // resize: resizes the table to be new_size.  new_size *MUST* be a
  // prime number, or everything will fail.
  void resize(int new_size) {
    // make a new table (nt),  and re-hash everyone in the
    // old table into the new table.

    // the next table, used for table resizing.
    sArray < SBucket > nt;
    nt.resize(new_size);

    for ( unsigned i=0;i<table.size();i++) {
      if ( table[i].active ) {
        unsigned a = entryOf(table[i].key, nt);
        nt[a].active = true;
        nt[a].key = table[i].key;
      }
    }

    // printf("resizing begin, %d current elements, new els = %d\n",_totalNumberEntries,new_size);
    // iterate over only non-empty entries in hash table.
    table.swap(nt);
    // clear out nt which is now old table.
    nt.clear();
  }

public:

  ////////////////////
  // constructor
  //    All entries in this hash table have the same size given
  //    by the argument arg_vsize.
  shash_set(unsigned approximateStartingSize = hash_abstract::HashTableDefaultApproxStartingSize) {
    _totalNumberEntries=0;
    findPrimesArrayIndex(approximateStartingSize);
    initialPrimesArrayIndex = primesArrayIndex;
    // create the actual hash table here.
    resize(HashTable_PrimesArray[primesArrayIndex]);

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
    _totalNumberEntries=0;
    primesArrayIndex=initialPrimesArrayIndex;
    resize(HashTable_PrimesArray[primesArrayIndex]);
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
  _Key* insert(_Key key, bool&foundp = hash_abstract::global_foundp) {

#ifdef COLLECT_COLLISION_STATISTICS
    numInserts++;
#endif

    // compute the address
    unsigned a = entryOf(key, table);
    if ( ! table[a].active ) {
      foundp = false;

      table[a].active = true;
      table[a].key = key;

      // time to resize if getting too big.
      // TODO: probably should resize a bit later than 1/2 entries being used.
      // TODO: precompute size at which we do a resize.
      if ( ++_totalNumberEntries >= table.size() / 2 ) {
	  if ( primesArrayIndex == (HashTable_SizePrimesArray - 1) )
	    error("ERROR: Hash table error, table size exceeds max size of %lu",
		  HashTable_PrimesArray[primesArrayIndex]);
	  resize(HashTable_PrimesArray[++primesArrayIndex]);
	  // need to re-get location
	  a = entryOf(key, table);
	  assert (table[a].active);
      }
    } else {
      foundp = true;
    }
    return &table[a].key;
  }


  ////////////////////////////////////////////////////////
  // search for key returning true if the key is found, otherwise
  // don't change the table. Return true when found.
  bool find(_Key key) {
    const unsigned a = entryOf(key, table);
    return table[a].active;
  }


  ////////////////////////////////////////////////////////
  // another version of find that also returns a pointer to the key in
  // the hash table itself. Returns true if the key is contained in
  // the set, false otherwise.
  bool find(_Key* key,_Key**& key_pp) {
    const unsigned a = entryOf(key,table);
    if (empty(table[a])) {
      key_pp = NULL;
      return false;
    } else {
      key_pp = &(table[a].key);
      return true;
    }
  }
};


#endif // defined VHASH_SET

