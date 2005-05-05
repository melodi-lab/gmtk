
#ifndef STR_HASH_MAP_H
#define STR_HASH_MAP_H


#include <cstring>


/**
 * string hash table
 */
template< typename DataT >
class strhash_map {

protected:
	/**
	 * bucket used in string hash tables
	 */
	struct HashBucket {
	public:
		/**
		 * construct inactive bucket
		 */
		HashBucket() : _str(NULL) {}
		/**
		 * destructor
		 */
		~HashBucket() { delete [] _str; }

		/**
		 * check whether bucket is active
		 * 
		 * @return true if the bucket is active 
		 */
		inline bool active() const { return _str != NULL; }

		/**
		 * check whether bucket is active
		 * 
		 * @return true if the bucket is inactive 
		 */
		inline bool notActive() const { return _str == NULL; }

		/**
		 * activate a bucket by assigning values
		 * 
		 * @param s key value
		 * @param d data value
		 */
		inline void activate(const char * s, const DataT &d) {
			_str = new char [strlen(s) + 1];
			strcpy(_str, s);
			data = d;
		}

		/**
		 * de-activate a bucket
		 */
		inline void deActivate() { delete [] _str; _str = NULL; }
		
		/**
		 * check equality of key
		 *
		 * @param s the key value
		 * @return true if keys are equal
		 */		
		inline bool keyEqual(const char *s) { return ! strcmp(_str, s); }

		/**
		 * check equality of key
		 *
		 * @param s the key value
		 * @return true if keys are different
		 */		
		inline bool keyNotEqual(const char *s) { return strcmp(_str, s); }

		/**
		 * get the key value
		 * 
		 * @return key value
		 */
		inline const char * getKey() const { return _str; }

		/** data */
		DataT data;

	protected:
		/** key value of bucket */
		char *_str;
	};


public:
	/**
	 * constructor
	 * 
	 * @param numBits number of bits for table size
	 */
	strhash_map(unsigned size = 0) : _size(0) {
		_numBits = numBits(size);
		if ( _numBits < 8 )
			_numBits = 8;
		_tableSize = 0x1 << _numBits;
		_halfFullSize = _tableSize >> 1;
		_mask = _tableSize - 1;
		_table = new HashBucket [_tableSize];
	}

	/**
	 * copy constructor
	 *
	 * @param table the table to be copied
	 */
	strhash_map(const strhash_map<DataT>& table) : _table(NULL), _tableSize(0) {
		*this = table;
	}
		
	/**
	 * destructor
	 */
	virtual ~strhash_map() { delete [] _table; }

	/**
	 * assignment operator
	 *
	 * @param table the table to be copied
	 * @return reference to itself
	 */
	const strhash_map<DataT>& operator = (const strhash_map<DataT>& table) {
		if ( this != &table ) {
			// resize memory
			if ( _tableSize != table._tableSize ) {
				_tableSize = table._tableSize;
				delete [] _table;
				_table = new HashBucket [_tableSize];
			}

			// copy data
			_size = table._size;
			_numBits = table._numBits;
			_halfFullSize = table._halfFullSize;
			_mask = table._mask;

			HashBucket * const end_p = _table + _tableSize;
			HashBucket *p = _table, *tp = table._table;
			while ( p != end_p ) {
				if ( tp->active() ) {
					p->activate(tp->getKey(), tp->data);
				}
			}
		}

		return *this;
	}

	/**
	 * get the number of active items
	 * 
	 * @return the number of data stored
	 */
	inline unsigned long size() const { return _size; }

	/**
	 * get the table size
	 * 
	 * @return the table size
	 */
	inline unsigned long tableSize() const { return _tableSize; }

	/**
	 * change the size of the table
	 *
	 * @param size new size of the table
	 */
	void resize(unsigned size) {
		_numBits = numBits(size);
		_tableSize = 0x1 << _numBits;
		_halfFullSize = _tableSize >> 1;
		_mask = _tableSize - 1;
		delete [] _table;
		_table = new HashBucket [_tableSize];
	}

	/**
	 * insert a data
	 * 
	 * @param key key of entry
	 * @param data data of entry
	 */
	void insert(const char* key, const DataT &data) {
		if ( _size > _halfFullSize )
			rehash();

		HashBucket* pos = findPos(key);

		if ( pos->notActive() ) {
			pos->activate(key, data);
			++_size;
		} else {
			pos->data = data;
		}
	}

	/**
	 * remove a data
	 * 
	 * @param key key of entry
	 */
	void remove(const char* key) {
		HashBucket* pos = findPos(key);

		if ( pos->active() ) {
			pos->deActivate();
			--_size;
		}
	}

	/**
	 * check existance according to key
	 *
	 * @param key the key of entry
	 * @return true if it exists in the hash table
	 */
	inline bool contains(const char* key) const {
		HashBucket* pos = findPos(key);
		return pos->active();
	}

	/**
	 * get data according to key
	 * 
	 * @param key the key of entry
	 * @param data the corresponding data
	 * @return true if key exists
	 */
	inline bool find(const char* key, DataT &data) const {
		HashBucket*pos = findPos(key);
		if ( pos->active() ) {
			data = pos->data;
			return true;
		} else
			return false;
	}

	/**
	 * get data according to key
	 * 
	 * @param key the key of entry
	 * @return coresponding data
	 * @throws Exception if key does not exist
	 */
	inline const DataT get(const char* key) const {
		HashBucket *pos = findPos(key);
		assert( pos->active() );
		return pos->data;
	}

	/**
	 * clear the hash table
	 */
	void clear() {
		HashBucket* const end_p = _table + _tableSize;
		for ( HashBucket *p = _table; p != end_p; ++p ) {
			p->deActivate();
		}
		_size = 0;
	}

protected:
	/**
	 * find the position of bucket according to key
	 * 
	 * @param key the key of entry
	 * @return pointer of bucket in the hash table
	 */
	HashBucket* findPos(const char* key) const {
		register unsigned long hashValue = hashing(key, 0);
		unsigned long pos = ((hashValue >> _numBits) ^ hashValue) & _mask;

		while ( _table[pos].active() && _table[pos].keyNotEqual(key) ) {
			hashValue = hashing(key, hashValue);
			pos = ((hashValue >> _numBits) ^ hashValue) & _mask;
		}

		return _table + pos;
	}

	/**
	 * double the hash table size
	 */
	void rehash() {
		HashBucket *tmp = _table;
		HashBucket * const end_p = tmp + _tableSize;

		++_numBits;
		_halfFullSize = _tableSize;
		_tableSize <<= 1;
		_mask = _tableSize - 1;
		_size = 0;

		_table = new HashBucket [_tableSize];

		for ( HashBucket *p = tmp; p != end_p; p++ ) {
			if ( p->active() )
				insert(p->getKey(), p->data);
		}

		delete [] tmp;
	}

	static inline unsigned long hashing(const char* str, unsigned hashValue) {
		while ( *str ) {
			hashValue += (hashValue <<3) + (*str++);
		}

		return hashValue;
	}

	/**
	 * calcualte the number of bits for hash table given the approximate size
	 *
	 * In this implemention, there is no loop.
	 *   1. calculate the number of leading zeros in x
	 *   2. calculate the number of bits
	 *   3. when x is not power of two, add one more bit
	 *
	 * @param x the approximately number of items in hash table
	 * @return number of bits for hash table
	 */
	static inline unsigned numBits(unsigned long x) {
		if (x == 0)
			return 1;

		unsigned n = 0;
		if ( x <= 0x0000FFFF ) { n += 16; x <<= 16; }
		if ( x <= 0x00FFFFFF ) { n +=  8; x <<=  8; }
		if ( x <= 0x0FFFFFFF ) { n +=  4; x <<=  4; }
		if ( x <= 0x3FFFFFFF ) { n +=  2; x <<=  2; }
		if ( x <= 0x7FFFFFFF ) { n +=  1; }

		return sizeof(unsigned long) * 8 - n + ((x & (x-1)) != 0);
	}

	/** table entries */
	HashBucket *_table;
	/** number of active entries */
	unsigned long _size;

	/** hash table bits */
	unsigned _numBits;
	/** hash table size */
	unsigned long _tableSize;
	/** rehash limit */
	unsigned long _halfFullSize;
	/** mask of hash key */
	unsigned long _mask;
};


#endif // ifndef G_HASH_TABLE_H

