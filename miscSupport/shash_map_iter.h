/*
 * shash_map_iter.h
 *
 *   General data structure for a scalar-based hash-table
 *   implementation of a map. Each element of the map has a scalar
 *   key. This hash map one has an iterator. Acutally, this object
 *   should probably be called shash_map_witer (With Iter).
 *
 * Written by Gang Ji <gang@ee.washington.edu>
 * Mods by <bilmes@ee.washington.edu>
 *
 * Copyright (C) 2003 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 *
 *
 *
 * $Header$
 *
 */


#ifndef SHASH_MAP_ITER_H
#define SHASH_MAP_ITER_H

#include "shash_map.h"

// _Key must be a basic type. We map from these scalars of _Key's
// to objects given by _Data.
template <class _Key, class _Data>
class shash_map_iter : public shash_map<_Key, _Data> {

public:

	////////////////////
	// constructor
	//    All entries in this hash table have the same size given
	//    by the argument arg_vsize.
	shash_map_iter(unsigned approximateStartingSize = hash_abstract::HashTableDefaultApproxStartingSize) : shash_map<_Key, _Data>(approximateStartingSize) {}

	shash_map_iter(const shash_map_iter<_Key, _Data> &map) : shash_map<_Key, _Data>((shash_map<_Key, _Data>)map){}

	const shash_map_iter& operator = (const shash_map_iter<_Key, _Data> &map) {
		return *this = (shash_map<_Key, _Data>)map;
	}

	///////////////////////////////////////////////////////
	// the iterator class for this object
	class iterator {
		friend class shash_map_iter;
		// current bucket
		typename shash_map<_Key,_Data>::MBucket* b;
		typename shash_map<_Key,_Data>::MBucket* end_b;
	public:
		iterator() : b(NULL), end_b(NULL) {}
 
		_Key key() { return b->key; }
		_Data& operator*() { return b->item; }

		bool next() {
			do {
				++b;
			} while ( b != end_b && (! b->active) );
			return b != end_b;
		}
	};
 
	// A version that takes an iterator as argument and so
	// can be used with an existing iterator w/o needing to
	// create tmp objects with lots of construction/destruction.
	const bool begin(iterator &it) {
	  it.b = this->table.ptr;
	  it.end_b = it.b + this->table.len();
	  while ( it.b != it.end_b && (! it.b->active) ) {
	    ++it.b;
	  }
	  // return true if there is an item here.
	  return it.b != it.end_b;
	}
};


#endif // defined VHASH_MAP_ITER
