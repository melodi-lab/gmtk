/*- -*- C++ -*-
 * hash_mtable
 *      .h file
 *
 * This is a data structure for hash multi-table.  The idea is to concatenate
 * a set of hash tables into one contiguous memory chunk.
 *
 *  Written by Gang Ji <gang@ee.washington.edu>
 *
 *  $Header$
 *
 * Copyright (c) 2001, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about
 * the suitability of this software for any purpose.  It is provided
 * "as is" without express or implied warranty.
 *
 */

#ifndef GMTK_HASH_M_TABLE_H
#define GMTK_HASH_M_TABLE_H


#include <cstdlib>
#include <cassert>

#include "fileParser.h"


/*-
 * This is special designed for unsigned hash map.
 *
 * Note: 1. The key can only be non-negative int.
 *       2. There is no boundary checking.
 *       3. There is no checking for insertion to a full table.
 */
template<typename DataT>
class HashMTable {

public:
	/*-
	 *-----------------------------------------------------------------------
	 * HashMTable::HashMTable
	 *      Default constructor
	 *
	 * Results:
	 *      None.
	 *
	 * Side Effects:
	 *      Create an empty table.
	 *
	 *-----------------------------------------------------------------------
	 */
	HashMTable() : _table(NULL), _tableSize(0), _size(0) {}

	/*-
	 *-----------------------------------------------------------------------
	 * HashMTable::HashMTable
	 *      Constructor.
	 *      The size should be precomputed by the calller.
	 *
	 * Results:
	 *      None.
	 *
	 * Side Effects:
	 *      None.
	 *
	 *-----------------------------------------------------------------------
	 */
	HashMTable(unsigned size) : _size(0) {
		_tableSize = size;

		if ( (_table = new MBucket [_tableSize]) == NULL )
			error("out of memory in HashMTable::HashMTable");
	}

	/*-
	 *-----------------------------------------------------------------------
	 * HashMTable::~HashMTable
	 *      Default destructor
	 *
	 * Results:
	 *      None.
	 *
	 * Side Effects:
	 *      Clean up the memory.
	 *
	 *-----------------------------------------------------------------------
	 */
	~HashMTable() {delete [] _table;}

	///////////////////////////////////////////////////////////
	// retrieving table size and number of active entries.
	unsigned tableSize() const {return _tableSize;}
	unsigned size() const {return _size;}

	/*-
	 *-----------------------------------------------------------------------
	 * HashMTable::resize
	 *      Set the size of the table to a new size.
	 *      This new size should be pre-computed by the caller.
	 *
	 * Results:
	 *      None.
	 *
	 * Side Effects:
	 *      There is no copy of the old data.  The original data
	 *      will be lost.
	 *
	 *-----------------------------------------------------------------------
	 */
	void resize(unsigned size) {
		delete [] _table;

		_tableSize = size;
		_size = 0;

		if ( (_table = new MBucket [_tableSize]) == NULL )
			error("out of memory in HashMTable::HashMTable");
	}

	/*-
	 *-----------------------------------------------------------------------
	 * HashMTable::insert
	 *      Insert an entry into the hash multi-table.
	 *      The offset is the offset of the hash sub-tabel block startingSize point.
	 *      The blockSize is the block size of the sub-table.
	 *
	 * Results:
	 *      None.
	 *
	 * Side Effects:
	 *      This will not check whether the caller is trying to insert into
	 *      a full table.  It is the caller's responsibility to maintain
	 *      the block size, offset and filling percentage.
	 *
	 *-----------------------------------------------------------------------
	 */
	void insert(int key, const DataT &item, const unsigned offset, const unsigned blockSize) {
		assert(key >= 0);
		assert(blockSize > 0);

		MBucket *pos = findPos(key, offset, blockSize);

		if ( pos->key < 0 )
			++_size;
		pos->key = key;
		pos->item = item;
	}

	/*-
	 *-----------------------------------------------------------------------
	 * HashMTable::find
	 *      Find an entry according to key value.
	 *      The offset is the offset of the hash sub-tabel block startingSize point.
	 *      The blockSize is the block size of the sub-table.
	 *
	 * Results:
	 *      Return the pointer to the value entry, null if not found.
	 *
	 * Side Effects:
	 *      None.
	 *
	 *-----------------------------------------------------------------------
	 */
	DataT* find(int key, const unsigned offset, const unsigned blockSize) {
		assert(key >= 0);
		if ( blockSize == 0 )
			return NULL;

		MBucket *pos = findPos(key, offset, blockSize);

		if ( pos->key >= 0 )
			return &pos->item;
		else
			return NULL;
	}

	///////////////////////////////////////////////////////////
	// Reading/writing the whole hash multi-table into/from a file.
	/*-
	 *-----------------------------------------------------------------------
	 * HashMTable::readHashMTableForComplexType
	 *      Read the whole hash multi-table from a file.
	 *      In ASCII format, the file looks like this
	 *          table_size size    % table size and number of active entries
	 *          index key value    % index in the table, key, and complex type
	 *
	 * Results:
	 *      None.
	 *
	 * Side Effects:
	 *      None.
	 *
	 *-----------------------------------------------------------------------
	 */
	void readHashMTableForComplexType(iDataStreamFile &ifs) {
		ifs.readUnsigned(_tableSize);
		ifs.readUnsigned(_size);

		// create buffer
		delete [] _table;
		if ( (_table = new MBucket [_tableSize]) == NULL )
			error("out of memory in HashMTable::readObject");

		unsigned index;
		for ( unsigned i = 0; i < _size; ++i ) {
			ifs.readUnsigned(index);
			ifs.readInt(_table[index].key);
			_table[index].item.readObject(ifs);
		}
	}

	/*-
	 *-----------------------------------------------------------------------
	 * HashMTable::writeHashMTableForComplexType
	 *      Write the whole hash multi-table into a file.
	 *      In ASCII format, the file looks like this
	 *          table_size size    % table size and number of active entries
	 *          index key value    % index in the table, key, and complex type
	 *
	 * Results:
	 *      None.
	 *
	 * Side Effects:
	 *      None.
	 *
	 *-----------------------------------------------------------------------
	 */
	void writeHashMTableForComplexType(oDataStreamFile &ofs) const {
		ofs.writeUnsigned(_tableSize);
		ofs.writeUnsigned(_size);
		ofs.nl();

		MBucket *p = _table;
		for ( unsigned i = 0; i < _tableSize; ++i, ++p ) {
			if ( p->key >= 0 ) {
				ofs.writeUnsigned(i);
				ofs.writeInt(p->key);
				p->item.writeObject(ofs);
				ofs.nl();
			}
		}
	}

	/*-
	 *-----------------------------------------------------------------------
	 * HashMTable::readHashMTableForBasicType
	 *      Read the whole hash multi-table from a file.
	 *      Basic types are something like 'float', 'int', ...
	 *      In ASCII format, the file looks like this
	 *          table_size size    % table size and number of active entries
	 *          index key value    % index in the table, key, and basic type
	 *
	 * Results:
	 *      None.
	 *
	 * Side Effects:
	 *      None.
	 *
	 *-----------------------------------------------------------------------
	 */
	void readHashMTableForBasicType(iDataStreamFile &ifs) {
		ifs.readUnsigned(_tableSize);
		ifs.readUnsigned(_size);

    	// create buffer
		delete [] _table;
		if ( (_table = new MBucket [_tableSize]) == NULL )
			error("out of memory in HashMTable::readObject");

		unsigned index;
		for ( unsigned i = 0; i < _size; ++i ) {
			ifs.readUnsigned(index);
			ifs.readInt(_table[index].key);
			ifs.read(_table[index].item);
		}
	}

	/*-
	 *-----------------------------------------------------------------------
	 * HashMTable::readHashMTableForBasicType
	 *      Read the whole hash multi-table from a file.
	 *      Basic types are something like 'float', 'int', ...
	 *      In ASCII format, the file looks like this
	 *          table_size size    % table size and number of active entries
	 *          index key value    % index in the table, key, and basic type
	 *
	 * Results:
	 *      None.
	 *
	 * Side Effects:
	 *      None.
	 *
	 *-----------------------------------------------------------------------
	 */
	void writeHashMTableForBasicType(oDataStreamFile &ofs) {
		ofs.writeUnsigned(_tableSize);
		ofs.writeUnsigned(_size);
		ofs.nl();

		MBucket *p = _table;
		for ( unsigned i = 0; i < _tableSize; ++i, ++p ) {
			if ( p->key >= 0 ) {
				ofs.writeUnsigned(i);
				ofs.writeInt(p->key);
				ofs.write(p->item);
				ofs.nl();
			}
		}
	}

protected:
	/*-
	 * bucket for the hash table
	 * In order to save memory, key < 0 means the entry is empty.
	 */
	struct MBucket {
		int key;
		DataT item;

		// key must be >=0 -1 means empty
		MBucket() : key(-1) {}
	};

	///////////////////////////////////////////////////////////
	// hashing functions.
	inline unsigned h1(const int key, const unsigned blockSize) const {
		return (key + 3367900314ul) % blockSize;
	}

	inline unsigned h2(const int key, const unsigned blockSize) const {
		return (key + 1) % (blockSize - 1) + 1;
	}

	///////////////////////////////////////////////////////////
	/*-
	 *-----------------------------------------------------------------------
	 * HashMTable::find
	 *      Retrieving the pointer according to key.
	 *      The offset is the offset of the hash sub-tabel block startingSize point.
	 *      The blockSize is the block size of the sub-table.
	 *
	 * Results:
	 *      Return the pointer to the hash entry.
	 *
	 * Side Effects:
	 *      None.
	 *
	 *-----------------------------------------------------------------------
	 */
	inline MBucket* findPos(const int key, const unsigned shift, const unsigned blockSize) const {
		unsigned a = h1(key, blockSize);

		MBucket* start = _table + shift;

		if ( ((start + a)->key < 0) || ((start + a)->key == key) )
			return _table + shift + a;

		const unsigned inc = h2(key, blockSize);
		do {
			a = (a + inc) % blockSize;
		} while ( ((start + a)->key >= 0) && ((start + a)->key != key) );

		return start + a;
	}

	MBucket *_table;			// hash multi-table entries
	unsigned _tableSize;		// the overall table size
	unsigned _size;				// number of active entries in the table
};


#endif
