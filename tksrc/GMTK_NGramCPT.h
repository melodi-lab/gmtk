/*- -*- C++ -*-
 * GMTK_NGramCPT
 *      .h file the GMTK_NGramCPT.cc file.
 *      Generic CPT for NGram language model
 *
 * Written by Gang Ji <gang@ee.washington.edu>
 * Modifications by J. Bilmes. <bilmes@ee.washington.edu>
 *
 *  This file contains an implementation for ngrams data type and its
 *  associated backoff data-structure. It supports the ARPA format
 *  N-gram language model format. There are actually three ways to use
 *  an ARPA language model in GMTK:
 *
 * 1. use lm file with Vocab object
 *          bigram
 *          1 % number of parents
 *          VOCAB_SIZE VOCAB_SIZE % cards
 *          ./DATA/bigram.arpa vocab % ARPA lm file and vocabulary object
 *
 * 2. use ascii indexing file created by gmtkNGramIndex
 *          bigram
 *          1 % number of parents
 *          VOCAB_SIZE VOCAB_SIZE % cards
 *          ./DATA/bigram.arpa.idx [ascii] % ARPA lm indexing file
 *
 * 3. use binary indexing file created by gmtkNGramIndex
 *          bigram
 *          1 % number of parents
 *          VOCAB_SIZE VOCAB_SIZE % cards
 *          ./DATA/bigram.arpa.idx [binary] % ARPA lm indexing file
 *
 * The latter two read in much faster, but first require the
 * conversion of the ARPA lm file to the GMTK idx file format.
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
#include "GMTK_DiscRV.h"
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
		_numParents = 0;
		_numberOfActiveIterators = 0;
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
	// Probability evaluation, compute Pr( child | parents ), and
	// iterator support. See GMTK_CPT.h for documentation.
	virtual void becomeAwareOfParentValues(vector< RV* >& parents, const RV* rv);
	virtual void begin(iterator& it, DiscRV* drv,logpr& p);
	virtual void becomeAwareOfParentValuesAndIterBegin(vector< RV* >& parents, iterator &it, DiscRV* drv, logpr& p);
	virtual logpr probGivenParents(vector< RV* >& parents, DiscRV* drv);
	virtual bool next(iterator &it,logpr& p);

	// Include here an extra routine that returns the probability
	// of the child 'val' given the parents are the assigned to
	// the set of values set during the most previous call to
	// becomeAwareOfParentValues.
	logpr probGivenParents(const int val);


	///////////////////////////////////////////////////////////
	// Given the current parent values, generate a random sample.
	virtual int randomSample(DiscRV*drv);

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
	void emStartIteration() {}
	void emIncrement(logpr prob,vector < RV* >& parents, RV*r) {}
	void emEndIteration() {}
	void emSwapCurAndNew() {}

	// return the number of parameters for object.
	virtual unsigned totalNumberParameters() {return _totalNumberOfParameters;}

	///////////////////////////////////////////////////////////////
	// virtual functions for objects to do the actual work.
	virtual void emStoreObjectsAccumulators(oDataStreamFile& ofile,
						bool writeLogVals = true,
						bool writeZeros = false) {}
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
	inline ContextHashEntry* findContextHashEntry(unsigned *context, unsigned length) {
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

	logpr probBackingOff(const int val, ContextHashEntry**ptr, unsigned numOfExistsParents);

	///////////////////////////////////////////////////////////////
	// context hash M table and probability hash M table
	HashMTable<ContextHashEntry> _contextTable;
	HashMTable<double> _probTable;
	// the starting block size for the hash M table
	unsigned _contextStartBlockSize;
	unsigned _probStartBlockSize;

	// parents values, TODO: convert this to fast sArray
	std::vector<void *> _contextEntriesStack;
	unsigned _numberOfActiveIterators;

	// total number of probabilities and backing-off weights
	unsigned _totalNumberOfParameters;

	char *_lmIndexFile;				// filename for indexing file.  see writeNGramIndexFile for file format.
	bool _lmIndexFileBin;			// indexing file is binary
};


#endif
