/*- -*- C++ -*-
 * GMTK_Vocab
 *      .h file the .cc file.
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


#ifndef GMTK_VOCAB_H
#define GMTK_VOCAB_H

#include "debug.h"
#include "GMTK_NamedObject.h"
 

/*-
 * ADT for vocabulary
 * This is designed for reading n-grams from a DARPA file.
 */
class Vocab : public NamedObject, public IM {

public:
	///////////////////////////////////////////////////////////
	// constructors and destructor
	Vocab();
	Vocab(unsigned size);
	virtual ~Vocab();

	///////////////////////////////////////////////////////////
	// indexing and retreving string
	unsigned index(const char *word) const;
	const char * word(unsigned index) const;

	void read(const char *filename);
	void read(iDataStreamFile& is);

	///////////////////////////////////////////////////////////
	// Misc methods
        // returns number of vocab items
	unsigned size() const {return _size;}
	void resize(unsigned size);

	///////////////////////////////////////////////////////////
	// Return whether <unk> is a word or not.
	bool unkIsWord() const {
		return _unkIndex < _size;
	}

protected:
	/**
	 * entry in the array
	 */
	struct HashEntry {
		char *key;		/// key for the entry (i.e., the vocab string).
		unsigned wid;	/// data for the entry (i.e., the word or vocab int id).

		// constructors and destructor
		HashEntry() : key(NULL) {}
		HashEntry(const char* theKey, int theWid);
		~HashEntry() {delete [] key;}

		const HashEntry& operator = (const HashEntry& he);
	};

	///////////////////////////////////////////////////////////
	// hash function
	inline unsigned hash(const char* key) const {
		register unsigned hashValue = 0;
		// a simple string hash function.
		while ( *key )
		  hashValue += (hashValue << 7) + *key++;   // this is more efficient

		return hashValue % _tableSize;
	}

	// returns the type of the sub-object in string
	// form that is suitable for printing and identifying
	// the type of the object.
	virtual const string typeName() {return std::string("Vocab");}

	///////////////////////////////////////////////////////////
	// insertion and location
	void insert(const char *word, unsigned wid);
	HashEntry* findPos(const char* key) const;

	unsigned _size;				// number of words
	unsigned _tableSize;		// table size
	HashEntry *_indexTable;		// index hash table
	const char **_stringTable;		// string array table

	// <unk> is a special token in a language model as is used to
	// represent all unknown words. All such words get the same
	// probability (and are considered similar in a context). This
	// is the index for <unk>.
	unsigned _unkIndex;			// index for <unk>
};


#endif
