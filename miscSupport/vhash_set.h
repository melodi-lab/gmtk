/*
 * vhash_set.h
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

#ifndef VHASH_SET_H
#define VHASH_SET_H

#include <stdio.h>
#include <stdlib.h>

#include "debug.h"
#include "sArray.h"
#include "hash_abstract.h"

// _Key must be a basic type. This is a hash set containing vectors (_Key*).
// All elements of the set have the same length, given by object constructor.
template <class _Key> 
class vhash_set : public hash_abstract {

#ifdef COLLECT_COLLISION_STATISTICS
public:
#endif

  //////////////////////////////////////////////////////////
  // The size of all the vectors in this set must be the
  // same. Here we store the size once, so they don't
  // need to be stored in each element.
  // @@@ should be const
  unsigned vsize;

  ////////////////////////////////////////////////////////////
  // the set "buckets", which includes only key 
  struct SBucket {
    _Key* key;
    SBucket() :key(NULL) {}
  };

  //////////////////////////////////////////////////////////
  // the actual hash table, an array of pointers to T's
  sArray < SBucket > table;

#ifdef COLLECT_COLLISION_STATISTICS
  unsigned maxCollisions;
  unsigned numCollisions;
  unsigned numInserts;
#endif


  //////////////////////////////////////////////////
  // return true if the two keys are equal.
  bool keyEqual(const _Key* k1,const _Key* k2) {
    const _Key* k1_endp = k1+vsize;
    do {
      if (*k1++ != *k2++)
	return false;
    } while (k1 != k1_endp);
    return true;
  }

  ///////////////////////////////////////////////////////////////////////
  // since this is a double hash table, we define two address
  // functions, h1() and h2().  h1() gives the starting position of in
  // the array of the key, and h2() gives the increment when we have a
  // collision.
  unsigned h1(const _Key* key, const unsigned arg_size)
  {
    const _Key* keyp = key + vsize;

    unsigned long a = vsize;
    // unsigned long a = 0;
    do {
      --keyp;
      // a =65599*a + (*keyp) + 1;
      // a = 402223*a + (*keyp) + 1;
      // a = 611953*a + (*keyp) + 1;
      a = 3367900313ul*a + (*keyp) + 1;
    } while (keyp != key);

    // See tk/tcl
    // a = 0;
    // do {
    // --keyp;
    // a += (a <<3) + (*keyp);
    // } while (keyp != key);

    return a % arg_size;
  }

  ///////////////////////////////////////////////////////////////////////
  // h2() is the increment of of key when we have a collision.  "The
  // value of h2(key) must be relatively prime to the hash-table m for
  // the entire hash table to be searched" (from Corman, Leiserson,
  // Rivest). Therefore, we have result of the h2() function satisfy
  //         1) it must be greater than zero 
  //         2) it must be strictly less than table.size()
  // 
  unsigned h2(const _Key* key, const unsigned arg_size)
  {

    unsigned long a=0;
    const _Key* keyp = key + vsize;
    do {
      --keyp;
      // a =65599*a + (*keyp) + 1;
      // a = 402223*a + (*keyp) + 1;
      // a = 611953*a + (*keyp) + 1;
      // a = 1500450271ul*a + (*keyp) + 1;
      a = 3267000013ul*a + (*keyp) + 1;
    } while (keyp != key);
    return (
       (a % (arg_size-1)) // this gives [0 : (vsize-2)]
       +
       1                  // this gives [1 : (vsize-1)] 
       );
  }

  ///////////////////////////////////////
  // return true if the bucket is empty
  bool empty(const SBucket* bucket) {
    return (bucket->key == NULL);
  }
  // return true if the bucket is empty
  bool empty(const SBucket& bucket) {
    return (bucket.key == NULL);
  }

  //////////////////////////////////////////////////////////////////
  // return the entry of key in table a_table
  unsigned entryOf(const _Key* key,
		   sArray<SBucket> & a_table) {
    const unsigned size = a_table.size();
    unsigned a = h1(key,size);

#ifdef COLLECT_COLLISION_STATISTICS
    unsigned collisions=0;
#endif

    // printf("entryOf: at entry %d\n",a);

    if (empty(a_table.ptr[a]) || keyEqual(a_table.ptr[a].key,key)) {
      return a;
    }
    const unsigned inc = h2(key,size);
    do {
#ifdef COLLECT_COLLISION_STATISTICS
      collisions++;
#endif
      a = (a+inc) % size;
      // printf("entryOf: C. now at entry %d\n",a);
    } while ( !empty(a_table.ptr[a])
	      &&
	      (!keyEqual(a_table.ptr[a].key,key)) 
	      );

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

    for (unsigned i=0;i<table.size();i++) {
      if (!empty(table.ptr[i])) {
	  unsigned a = entryOf(table.ptr[i].key,nt);
	  nt[a].key = table.ptr[i].key;
      }
    }

    // printf("resizing begin, %d current elements, new els = %d\n",_totalNumberEntries,new_size);
    // iterate over only non-empty entries in hash table.
    table.swap(nt);
    // clear out nt which is now old table.
    nt.clear();
  }

public:

  // Dummy constructor to create an invalid object to be re-constructed
  // later.
  vhash_set() {}

  ////////////////////
  // constructor
  //    All entries in this hash table have the same size given
  //    by the argument arg_vsize.
  vhash_set(const unsigned arg_vsize,
	    unsigned approximateStartingSize = 
	    hash_abstract::HashTableDefaultApproxStartingSize)
    : vsize(arg_vsize) 
  {
    // make sure we hash at least one element, otherwise do {} while()'s won't work.
    assert (vsize > 0);

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
  _Key* insert(_Key* key,
	       bool&foundp = hash_abstract::global_foundp) {

#ifdef COLLECT_COLLISION_STATISTICS
    numInserts++;
#endif

    // compute the address
    unsigned a = entryOf(key,table);
    // printf("inserting to entry %d, empty = %d\n",a,empty(table[a]));
    if (empty(table.ptr[a])) {
      foundp = false;

      table.ptr[a].key = key;

      // time to resize if getting too big.
      // TODO: probably should resize a bit later than 1/2 entries being used.
      // TODO: precompute size at which we do a resize.
      if (++_totalNumberEntries >= table.size()/2) {
	  if (primesArrayIndex == (HashTable_SizePrimesArray-1)) 
	    error("ERROR: Hash table error, table size exceeds max size of %lu",
		  HashTable_PrimesArray[primesArrayIndex]);
	  resize(HashTable_PrimesArray[++primesArrayIndex]);
	  // need to re-get location
	  a = entryOf(key,table);
	  assert (!empty(table.ptr[a]));
      }
    } else {
      foundp = true;
    }
    return table.ptr[a].key;
  }


  ////////////////////////////////////////////////////////
  // search for key returning true if the key is found, otherwise
  // don't change the table. Return true when found.
  bool find(_Key* key) {
    const unsigned a = entryOf(key,table);
    // printf("find: entry %d\n",a);
    return !empty(table.ptr[a]);
  }


  ////////////////////////////////////////////////////////
  // another version of find that also returns a pointer to the key in
  // the hash table itself. Returns true if the key is contained in
  // the set, false otherwise.
  bool find(_Key* key,_Key**& key_pp) {
    const unsigned a = entryOf(key,table);
    if (empty(table.ptr[a])) {
      key_pp = NULL;
      return false;
    } else {
      key_pp = &(table.ptr[a].key);
      return true;
    }
  }



};


#endif // defined VHASH_SET

