/*
 * hash_map_list.h
 *
 *   A general data structure for a hash-table with efficient (linear)
 *   iterators through all the elements in the hash table. It is a
 *   hash table but where each entry is in a linked list.
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

#ifndef HASH_MAP_LIST_H
#define HASH_MAP_LIST_H

#include <stdio.h>
#include <stdlib.h>

#include "debug.h"
#include "sArray.h"
#include "hash_abstract.h"

template <class _Key, class _Data> 
class hash_map_list : public hash_abstract {

  ////////////////////////////////
  // the "buckets", which include a data item and 
  // a pointer to the next bucket thereby forming
  // a simple linked list.
  struct Bucket {
    _Key key;
    _Data item;
    // a bucket is considered empty when pointer (next == NULL).
    Bucket* next;
    Bucket() : next(NULL) {}
  };

  //////////////////////////////////////////////////////////
  // the actual hash table, an array of pointers to T's
  sArray < Bucket > table;

  // For getting to the singular linked list.
  Bucket* first;
  Bucket* last;


#ifdef COLLECT_COLLISION_STATISTICS
public:
  unsigned maxCollisions;
  unsigned numCollisions;
  unsigned numInserts;
private:
#endif

  //////////////////////////////////////////////////////////////////
  // since this is a double hash table, we define two address
  // functions, h1() and h2().  h1() gives the starting position of in
  // the array of the key, and h2() gives the increment when we have a
  // collision.
  unsigned h1(const _Key key, const unsigned arg_size)
  {
    unsigned long a;
    // a =65599*key + 1;
    // a = 402223*key + 1;
    // a = 611953*key + 1;
    // a = 3367900313ul*key + 1;
    // a = 3367900313ul*key + key + 1;
    a = key+1;
    // TODO: when we know cardinality of a random variable,
    // we can scale a by tableSize/cardinality which might
    // give more uniform distribution over bucket entries.
    return a % arg_size;
  }

  ////////////////////////////////////////////////////////////
  // h2() is the increment of of key when we have a collision.  "The
  // value of h2(key) must be relatively prime to the hash-table m for
  // the entire hash table to be searched" (from Corman, Leiserson,
  // Rivest). Therefore, we have result of the h2() function satisfy
  //         1) it must be greater than zero 
  //         2) it must be strictly less than table.size()
  // 
  unsigned h2(const _Key key, const unsigned arg_size)
  {
    unsigned long a;

    // a =65599*key + 1;
    // a = 402223*key + 1;
    // a = 611953*key + 1;
    // a = 1500450271ul*key + 1;
    // a = 3267000013ul*a + key + 1;
    // a = 3267000013ul*key + key + 1;
    a = key+1;
    return (
       (a % (arg_size-1)) // this gives [0 : (size-2)]
       +
       1                  // this gives [1 : (size-1)] 
       );
  }


  ///////////////////
  // return true when the bucket is empty.
  bool empty(const Bucket* bucket,
	     const Bucket* a_last) {
    // a bucket is empty when
    //   1) it is not a_last, and
    //   2) it's next pointer is NULL
    return (bucket->next == NULL && bucket != a_last); 
  }

  /////////////////////////////////////////////////////////////////////
  // return the entry of key in table a_table with last bucket a_last.
  unsigned entryOf(_Key key,
		   sArray<Bucket> & a_table,
		   const Bucket* a_last) {
    const unsigned size = a_table.size();
    unsigned a = h1(key,size);

#ifdef COLLECT_COLLISION_STATISTICS
    unsigned collisions=0;
#endif

    if ( empty(&a_table[a],a_last) ||
	 (a_table[a].key == key) ) {
      return a;
    }
    const unsigned inc = h2(key,size);
    do {
#ifdef COLLECT_COLLISION_STATISTICS
      collisions++;
#endif
      a = (a+inc) % size;
    } while ( !empty(&a_table[a],a_last)
	      &&
	      (a_table[a].key != key) 
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
    sArray < Bucket > nt;

    Bucket* nfirst=NULL;
    Bucket* nlast=NULL;
    nt.resize(new_size);

    // printf("resizing begin, %d current elements\n",_totalNumberEntries);
    // iterate over only non-empty entries in hash table.
    iterator it = begin();
    for (; it != end(); it++) {
      unsigned a = entryOf(it.b->key,nt,nlast);
      nt[a].key = it.b->key;
      nt[a].item = it.b->item;
      // printf("resizing key = %d, item = %f\n",it.b->key,it.b->item);

      if (nfirst == NULL)
	nfirst = nlast = &nt[a];      
      else {
	nlast->next = &nt[a];
	nlast = &nt[a];
      }
    }
    // finally, swap nt with the existing table.
    table.swap(nt);
    // clear out nt which is now old table.
    nt.clear();
    // set up new first and last pointers.
    first = nfirst;
    last = nlast;
  }


public:

  ////////////////////
  // constructor
  hash_map_list(unsigned approximateStartingSize = 
		hash_abstract::HASH_TABLE_DEFAULT_APPROX_STARTING_SIZE) {
    _totalNumberEntries=0;
    findPrimesArrayIndex(approximateStartingSize);
    initialPrimesArrayIndex = primesArrayIndex;
    first = last = NULL;
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
	     hash_abstract::HASH_TABLE_DEFAULT_APPROX_STARTING_SIZE) {
    table.clear();
    _totalNumberEntries=0;
    primesArrayIndex=initialPrimesArrayIndex;
    first = last = NULL;
    resize(HashTable_PrimesArray[primesArrayIndex]);
  }

#ifdef COLLECT_COLLISION_STATISTICS
  // clear out the statistics
  void clearStats() {
    maxCollisions = numCollisions = numInserts = 0;
  }
#endif

  ///////////////////////////////////////////////////////
  // insert an item <key,val> into the hash table.
  // Return a pointer to the data item in the hash table
  // after it has been inserted. This way, it is possible
  // to insert an empty item (using a _Data() constructor)
  // and then adjust its value later via the returned pointer.
  // The foundp argument is set to true when the key has
  // been found. In this case (when it was found), the existing
  // _Data item is not adjusted to be val.
  _Data* insert(_Key key,
		_Data val,
		bool&foundp = hash_abstract::global_foundp) {

#ifdef COLLECT_COLLISION_STATISTICS
    numInserts++;
#endif

    // compute the address
    unsigned a = entryOf(key,table,last);
    // printf("inserting %d to entry %d, empty = %d\n",key,a,empty(&table[a],last));
    if (empty(&table[a],last)) {
      foundp = false;

      table[a].key = key;
      table[a].item = val;
      table[a].next = NULL;

      if (first == NULL)
	first = last = &table[a];      
      else {
	last->next = &table[a];
	last = &table[a];
      }
      // time to resize if getting too big.
      // TODO: probably should resize a bit later than 1/2 entries being used.
      if (++_totalNumberEntries >= table.size()/2) {
	  if (primesArrayIndex == (HashTable_SizePrimesArray-1)) 
	    error("ERROR: Hash table error, table size exceeds max size of %lu",
		  HashTable_PrimesArray[primesArrayIndex]);
	  resize(HashTable_PrimesArray[++primesArrayIndex]);
	  a = entryOf(key,table,last);
      }
    } else {
      foundp = true;
    }
    return &table[a].item;
  }

  ////////////////////////////////////////////////////////
  // search for key returning true if the key is found, otherwise
  // don't change the table. Return pointer to the data item
  // when found, and NULL when not found.
  _Data* find(_Key key) {
    const unsigned a = entryOf(key,table,last);
    // const bool foundp = ((table[a].next != NULL) || (last == &table[a]));
    const bool foundp = !empty(&table[a],last);
    if (foundp) 
      return &table[a].item;
    else
      return NULL;
  }

  /////////////////////////////////////////////////////////
  // search for key, returning a reference to the data item
  // regardless if the value was found or not. This will
  // insert the key into the hash table (similar to the
  // way STL's operator[] works with the hash_set.h and hash_map.h) 
  _Data& operator[](_Key key) {
    _Data* resp = insert(key,_Data());
    // return a reference to the location regardless
    // of if it was found or not.
    return *resp;
  }

  ///////////////////////////////////////////
  // a dummy type that is used to simulate STL's
  // way of doing iterator end (i.e., it == obj.end()).
  typedef char iterator_end_object;

  //////////////////////////////////////////////////////
  // the iterator class for this object.
  class iterator {

    friend class hash_map_list;
    // current bucket 
    Bucket* b;
  public:
    iterator() : b(NULL) {}

    // move to the next bucket
    iterator& operator ++(int) {
      if (b)
	b = b->next;
      return *this;
    }
    _Key key() { return b->key; }
    _Data& operator*() { return b->item; }

    // include an end() operator here since we know it here.
    bool end() { return (b == NULL); }
    // use operator == only with a dummy end object.
    bool operator == (iterator_end_object eo) {
      return end();
    }
    bool operator != (iterator_end_object eo) {
      return !end();
    }

  };

  // A version that takes an iterator as argument and so
  // can be used with an existing iterator w/o needing to
  // create tmp objects with lots of construction/destruction.
  void begin(iterator& it) {
    it.b = first;
  }

  // STL like iterator over all elements in the hash table
  iterator begin() {
    iterator it;
    begin(it);
    return it;
  }

  // Begin an iterator starting at key if it exists, otherwise 
  // create an iterator which end is true.
  void begin(iterator& it,_Key key) {
    const unsigned a = entryOf(key,table,last);
    const bool foundp = !empty(&table[a],last);
    if (foundp) {
      it.b = &table[a];
    } else {
      it.b = NULL;
    }
  }

  // return an iterator end object, any one will do, just
  // to trigger the end call of the iterator.
  iterator_end_object end() { return (iterator_end_object)0; }

};


#endif // defined HASH_MAP_LIST

