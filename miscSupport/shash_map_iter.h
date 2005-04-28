/*
 * shash_map.h
 *   General data structure for a scalar-based hash-table implementation
 *   of a map. Each element of the map has a scalar key.
 *
 * Written by Gang Ji <gang@ee.washington.edu>
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


#ifndef SHASH_MAP2_H
#define SHASH_MAP2_H

#include "shash_map.h"

// _Key must be a basic type. We map from these scalars of _Key's
// to objects given by _Data.
template <class _Key, class _Data>
class shash_map2 : public shash_map<_Key, _Data> {

public:

	////////////////////
	// constructor
	//    All entries in this hash table have the same size given
	//    by the argument arg_vsize.
	shash_map2(unsigned approximateStartingSize = hash_abstract::HashTableDefaultApproxStartingSize) : shash_map<_Key, _Data>(approximateStartingSize) {}

	shash_map2(const shash_map2<_Key, _Data> &map) : shash_map2((shash_map<_Key, _Data>)map){}

	const shash_map2& operator = (const shash_map2<_Key, _Data> &map) {
		return *this = (shash_map<_Key, _Data>)map;
	}

	///////////////////////////////////////////////////////
	// the iterator class for this object
	class iterator {
		friend class shash_map2;
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
	const iterator begin() {
		iterator it;
		it.b = this->table.ptr;
		it.end_b = it.b + this->table.len();
		while ( it.b != it.end_b && (! it.b->active) ) {
			++it.b;
		}

		return it;
	}
	const iterator end() {
		iterator it;
		it.b = it.end_b = this->table.ptr + this->table.len();
	}
};


#endif // defined VHASH_MAP2

