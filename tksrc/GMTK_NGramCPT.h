/*- -*- C++ -*-
 * GMTK_NGramCPT
 *      .h file the GMTK_NGramCPT.cc file.
 *      Generic CPT for NGram language model
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

#ifndef GMTK_NGRAM_CPT_H
#define GMTK_NGRAM_CPT_H


#include "GMTK_CPT.h"
#include "GMTK_Vocab.h"
#include "GMTK_RandomVariable.h"
#include "fileParser.h"
#include "hash_mtable.h"


/*-
 * CPT for n-gram language models
 */
class NGramCPT : public CPT {

public:

	/*-
	 *-----------------------------------------------------------------------
	 * NGramCPT::NGramCPT
	 *      Default conxtructor
	 *
	 * Results:
	 *      Initialize values
	 *
	 * Side Effects:
	 *      None
	 *
	 *-----------------------------------------------------------------------
	 */
	NGramCPT() : CPT(di_NGramCPT), _contextStartBlockSize(0), _probStartBlockSize(0), _totalNumberOfParameters(0), _lmIndexFile(NULL), _lmIndexFileBin(false) {
		_numParents = _numExistParents = 0;
	}

	/*-
	 *-----------------------------------------------------------------------
	 * ~NGramCPT::~NGramCPT
	 *      Default destructor
	 *
	 * Results:
	 *      Clean up memory
	 *
	 * Side Effects:
	 *      None
	 *
	 *-----------------------------------------------------------------------
	 */
	~NGramCPT() {
		delete [] _lmIndexFile;
	}

	///////////////////////////////////////////////////////////  
	// virtual functions from class CPT

	///////////////////////////////////////////////////////////
	// Semi-constructors: useful for debugging.
	// Functions to force the internal structures to be particular values.
	// Force the number of parents to be such.
	virtual void setNumParents(const unsigned nParents);

	// Allocate memory, etc. for the internal data structures
	// for this CPT, depending on current _numParents & cardinalities.
	virtual void allocateBasicInternalStructures() {}

	///////////////////////////////////////////////////////////  
	// Probability evaluation, compute Pr( child | parents )
	// 
	// becomeAwareOfParentValues: sets the parent values to a particular
	// assignment. All subsequent calls to to probGivenParents
	// will return the probability of the RV given that the
	// parents are at the particular value.
	virtual void becomeAwareOfParentValues(vector<int>& parentValues, vector<int>& cards);
	// Another version of becomeAwareOfParentValues but this
	// one explicitely takes an array of random variable parents.
	virtual void becomeAwareOfParentValues(vector<RandomVariable *>& parents);
	// return the probability of 'val' given the parents are the
	// assigned to the set of values set during the most previous call to 
	// becomeAwareOfParentValues.
	virtual logpr probGivenParents(const int val);
	// Similar to the above two. This is convenient for one time
	// probability evaluation.
	virtual logpr probGivenParents(vector<int>& parentValues, vector<int>& cards, const int val);
	virtual logpr probGivenParents(vector<RandomVariable *>& parents, const int val);

	// returns an iterator for the first one.
	virtual iterator begin();

	// creates an iterator for the first one.
	virtual void begin(iterator& it);
	virtual void becomeAwareOfParentValuesAndIterBegin(vector<RandomVariable *>& parents, iterator &it);
	// returns true if iterate is at end state
	virtual bool end(iterator &it)	{return it.value >= ucard();}
	// Given a current iterator, return true if it
	// is a valid next value, otherwise return false so
	// a loop can terminate.
	virtual bool next(iterator &);

	///////////////////////////////////////////////////////////  
	// Given the current parent values, generate a random sample.
	virtual int randomSample();


	///////////////////////////////////////////////////////////  
	// Re-normalize the output distributions
	virtual void normalize() {}
	// set all values to random values.
	virtual void makeRandom() {}
	// set all values to uniform values.
	virtual void makeUniform() {}

	///////////////////////////////////////////////////////////    
	// read in the basic parameters, assuming file pointer 
	// is located at the correct position.
	virtual void read(iDataStreamFile& is);
	///////////////////////////////////////////////////////////    
	// Do nothing.
	virtual void write(oDataStreamFile& os) {}

	void read(const char *lmFile, const Vocab &vocab);
	void readNGramIndexFile(iDataStreamFile &is);
	void writeNGramIndexFile(oDataStreamFile &os);


	////////////////////////////////////////////////////////////////////////////
	// from base class EMable

	////////////////////////////////////////////////////////////////////////////
	// if swap bit not set, swaps the current and new parameters, set swap bit.
	// otherwise does nothing.
	virtual void emSwapCurAndNew() {}

	// return the number of parameters for object.
	virtual unsigned totalNumberParameters() {return _totalNumberOfParameters;}

	///////////////////////////////////////////////////////////////
	// virtual functions for objects to do the actual work.
	virtual void emStoreObjectsAccumulators(oDataStreamFile& ofile) {}
	virtual void emLoadObjectsDummyAccumulators(iDataStreamFile& ifile) {}
	virtual void emZeroOutObjectsAccumulators() {}
	virtual void emLoadObjectsAccumulators(iDataStreamFile& ifile) {}
	virtual void emAccumulateObjectsAccumulators(iDataStreamFile& ifile) {}

protected:
	/*-
	 * hash entry type for hash multi-table.
	 */
	struct ContextHashEntry {
		double bow;						// backing-off weight

		// info for the probability hash table
		unsigned probOffset;			// offset for probability hash table
		unsigned probBlockSize;			// block size for hash table

		// info to the next level context hash table
		unsigned nextContextOffset;		// offset for next level context hash table
		unsigned nextContextBlockSize;	// block size

		ContextHashEntry() : bow(0), probOffset(0), probBlockSize(0), nextContextOffset(0), nextContextBlockSize(0) {}

		// for serialization
		void readObject(iDataStreamFile &ifs);
		void writeObject(oDataStreamFile &ofs);
	};

	/*-
	 *-----------------------------------------------------------------------
	 * Function
	 *      Find context hash entry from context and length.
	 *
	 * Results:
	 *      Return the pointer to the context hash entry.
	 *
	 * Side Effects:
	 *      For context of zero length, please use _contextStartBlockSize
	 *      with zero shift.
	 *
	 *-----------------------------------------------------------------------
	 */
	inline ContextHashEntry* findContextHashEntry(int *context, unsigned length) {
		if ( length == 0 )
			error("cannot refer to unigram in context search");

		int i = length -1;
		ContextHashEntry* entry = _contextTable.find(context[i], 0, _contextStartBlockSize);
		while ( (entry != NULL) && i > 0 )
			entry = _contextTable.find(context[--i], entry->nextContextOffset, entry->nextContextBlockSize);

		return entry;
	}

	// returns the type of the sub-object in string
	// form that is suitable for printing and identifying
	// the type of the object.
	virtual const string typeName() {return std::string("NGramCPT");}

	///////////////////////////////////////////////////////////////
	// context hash M table and probability hash M table
	HashMTable<ContextHashEntry> _contextTable;
	HashMTable<double> _probTable;
	// the starting block size for the hash M table
	unsigned _contextStartBlockSize;
	unsigned _probStartBlockSize;

	// parents values
	std::vector<ContextHashEntry*> _contextPointers;
	unsigned _numExistParents;		// The size of _parentValues is fixed.  This is added so that lower order ngram can share the table.

	// total number of probabilities and backing-off weights
	unsigned _totalNumberOfParameters;

	char *_lmIndexFile;				// filename for indexing file.  see writeNGramIndexFile for file format.
	bool _lmIndexFileBin;			// indexing file is binary
};


#endif
