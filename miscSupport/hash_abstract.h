/*
 * abstract_hash.h
 *   Generic stuff for hash tables.
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

#ifndef HASH_ABSTRACT_H
#define HASH_ABSTRACT_H

#include "debug.h"
#include "error.h"

class hash_abstract {

protected:

  static const unsigned HashTableDefaultApproxStartingSize;
  static const unsigned HashTable_SizePrimesArray;
  static const unsigned HashTable_PrimesArray[];
  
  // NOTE: this variable makes the code non-reentrant for threads
  static bool global_foundp;

  ////////////////////////////////
  // total number of entries in the hash table
  unsigned _totalNumberEntries;

  ////////////////////////////////
  // the index into the above prime array of the current size.
  unsigned primesArrayIndex; 

  ////////////////////////////////
  // the initial starting index into the above prime 
  // array based on the contructors approximate starting
  // size.
  unsigned initialPrimesArrayIndex; 


  //////////////////////////////////
  // find the starting index in the array of primes starting
  // at the approximateStartingSize (i.e., finds the next prime
  // greater).
  void findPrimesArrayIndex(unsigned approximateStartingSize) {
    // do a linear search for now, but could do bin-search.
    for (primesArrayIndex=0;
	 primesArrayIndex<HashTable_SizePrimesArray;
	 primesArrayIndex++) {
      if (HashTable_PrimesArray[primesArrayIndex] >= approximateStartingSize)
	break;
    }
    if (primesArrayIndex == HashTable_SizePrimesArray)
      error("ERROR: hash_map_list can't create hash table of such a large approximate size %u\n",approximateStartingSize);
  }

public:

  hash_abstract() : _totalNumberEntries(0) {}

  ////////////////////////////////////////////////////////////////
  // return the total number of entries that have been inserted into
  // the hash table (this is different than the total allocated size
  // of any table).
  unsigned totalNumberEntries() { return _totalNumberEntries; }

};


#endif // defined HASH_ABSTRACT

