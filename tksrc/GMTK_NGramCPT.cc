/*-
 * GMTK_NGramCPT.cc
 *     N-Gram CPT
 *
 * Written by Gang Ji <gang@ee.washington.edu>
 *
 * Copyright (c) 2001, < fill in later >
 *
 *  $Header$
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle, and Jeff Bilmes make no representations about
 * the suitability of this software for any purpose.  It is provided
 * "as is" without express or implied warranty.
 *
 */


#include <vector>
#include <queue>
#include <cmath>
#include <cassert>


#include "GMTK_NGramCPT.h"
#include "GMTK_GMParms.h"
#include "fileParser.h"
#include "vhash_map.h"
#include "rand.h"
#include "error.h"


#define MAX_LINE_LENGTH 1024


// ARPA format uses log10.
// this should be commented out
//#define M_LN10  1.0


// this function is defined in GMTK_Vocab.cc
extern unsigned nextPrime(unsigned n);


/*-
 * This ADT is used only when loading from a ARPA lm file.
 * A Trie ADT is utilized for loading.
 */
struct ContextTreeEntry {
	double bow;		// backing-off weight

	// information for the prob table
	unsigned probOffset;
	unsigned probBlockSize;

	// information for the context table
	unsigned contextOffset;
	unsigned contextBlockSize;

	// pointer to next level context hash table
	vhash_map<int, ContextTreeEntry> *next;

	// information for the next level context
	// hash table (no iterations for hashmap)
	std::vector<int> keys;

	ContextTreeEntry() : bow(0), probOffset(0), probBlockSize(0), contextOffset(0), contextBlockSize(0), next(NULL) {}
 	~ContextTreeEntry() {delete next;}
};


/*-
 *-----------------------------------------------------------------------
 * NGramCPT::setNumParents
 *      Set the number of parents for ngram (n-1)
 *
 * Results:
 *      Set the _parentValues size.
 *
 * Side Effects:
 *      Set the _parentValues size.
 *
 *-----------------------------------------------------------------------
 */
void NGramCPT::setNumParents(const unsigned nParents) {
	_numParents = nParents;
	_parentValues.resize(_numParents);
	cardinalities.resize(_numParents);
}


/*-
 *-----------------------------------------------------------------------
 * NGramCPT::becomeAwareOfParentValues
 *      Record the value of parents.
 *
 * Results:
 *
 * Side Effects:
 *      Unlike other CPTs, n-gram models requirs backing-off support.
 *      In this case, only the value of parents are recorded.
 *
 *-----------------------------------------------------------------------
 */
void NGramCPT::becomeAwareOfParentValues(vector<int>& parentValues, vector<int>& cards) {
	// TODO: This will vanish.
	_numExistParents = parentValues.size();
	assert(_numExistParents <= _numParents);		// This is not a "==" so that lower order ngram can use it.
	assert(cards.size() == _numExistParents);

	for ( unsigned i = 0; i < _numExistParents; ++i ) {
		if ( parentValues[i] < 0 || parentValues[i] >= cards[i] )
			error("NGramCPT:becomeAwareOfParentValues: Invalid parent value for parent %d, parentValue = %d but card = %d\n", i, parentValues[i], cardinalities[i]);
		_parentValues[i] = parentValues[i];
	}
}


/*-
 *-----------------------------------------------------------------------
 * NGramCPT::becomeAwareOfParentValues
 *      Record the value of parents.
 *
 * Results:
 *
 * Side Effects:
 *      Unlike other CPTs, n-gram models requirs backing-off support.
 *      In this case, only the value of parents are recorded.
 *
 *-----------------------------------------------------------------------
 */
void NGramCPT::becomeAwareOfParentValues(vector<RandomVariable *>& parents) {
	// TODO: instead of copy the parents values, record the hash entry pointers.
	_numExistParents = parents.size();
	assert(_numExistParents <= _numParents);		// This is not a "==" so that lower order ngram can use it.

	for ( unsigned i = 0; i < _numExistParents; ++i ) {
		if ( parents[i]->val < 0 || parents[i]->val >= cardinalities[i])
			error("ERROR:becomeAwareOfParentValues. NGram CPT %s, invalid parent value for parent %s(%d) (parent number %d), parentValue = %d but RV cardinality = %d\n",
				name().c_str(), parents[i]->name().c_str(),parents[i]->frame(), i, parents[i]->val,cardinalities[i]);

		_parentValues[i] = parents[i]->val;
	}
}


/*-
 *-----------------------------------------------------------------------
 * NGramCPT::probGivenParents
 *      Retrieve the probability for a given value.
 *
 * Results:
 *      Return the probability given the parents.
 *
 * Side Effects:
 *
 *-----------------------------------------------------------------------
 */
logpr NGramCPT::probGivenParents(const int val) {
	// TODO: make it faster by not copying the parent values.
	double *p = _probTable.find(val, 0, _probStartBlockSize);

	double prob;

	if ( p != NULL )
		prob = *p;
	else
		return logpr(0.0);

	if ( _numExistParents == 0 )
		return logpr(NULL, prob);

	int i = _numExistParents - 1;

	// the first one is special since the shift is 0
	ContextHashEntry *ce = _contextTable.find(_parentValues[i], 0, _contextStartBlockSize);
	if ( ce == NULL )
		return logpr(NULL, prob);

	p = _probTable.find(val, ce->probOffset, ce->probBlockSize);
	if ( p != NULL )
 		prob = *p;
	else
 		prob += ce->bow;

	// search next context
	while ( --i >= 0 ) {
		ce = _contextTable.find(_parentValues[i], ce->nextContextOffset, ce->nextContextBlockSize);
		if ( ce == NULL )
			return logpr(NULL, prob);

		p = _probTable.find(val, ce->probOffset, ce->probBlockSize);
		if ( p != NULL )
			prob = *p;
		else
			prob += ce->bow;
	}

	return logpr(NULL, prob); // this will create a logp with log value
}


/*-
 *-----------------------------------------------------------------------
 * NGramCPT::probGivenParents
 *      Retrieve the probability for a given value.
 *
 * Results:
 *      Return the probability given the parents.
 *
 * Side Effects:
 *
 *-----------------------------------------------------------------------
 */
logpr NGramCPT::probGivenParents(vector<int>& parentValues, vector<int>& cards, const int val) {
	becomeAwareOfParentValues(parentValues, cards);
	return probGivenParents(val);
}


/*-
 *-----------------------------------------------------------------------
 * NGramCPT::probGivenParents
 *      Retrieve the probability for a given value.
 *
 * Results:
 *      Return the probability given the parents.
 *
 * Side Effects:
 *
 *-----------------------------------------------------------------------
 */
logpr NGramCPT::probGivenParents(vector<RandomVariable *>& parents, const int val) {
	becomeAwareOfParentValues(parents);
	return probGivenParents(val);
}


/*-
 *-----------------------------------------------------------------------
 * NGramCPT::begin
 *      Beginning of an iterator.
 *
 * Results:
 *      Return the iterator from the begin of the ngram.
 *
 * Side Effects:
 *
 *-----------------------------------------------------------------------
 */
CPT::iterator NGramCPT::begin() {
	// TODO: this will vanish.
	iterator it(this);
	*(&it.value) = 0;
	it.probVal = probGivenParents(it.value);

	if ( it.probVal.essentially_zero() ) {
		// go to first entry which is not zero.
		do {
			it.value++;
			// We keep the following assertion as we
			// must have that at least one entry is non-zero.
			// The read code of the MDCPT should ensure this
			// as sure all parameter update procedures.
			assert(it.value < ucard());
			it.probVal = probGivenParents(it.value);
		} while (it.probVal.essentially_zero());
	}

	return it;
}


/*-
 *-----------------------------------------------------------------------
 * NGramCPT::begin
 *      Set an iterator to the begin.
 *
 * Results:
 *      None.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
void NGramCPT::begin(CPT::iterator& it) {
	it.value = 0;
	it.probVal = probGivenParents(it.value);

	if (it.probVal.essentially_zero()) {
		// go to first entry which is not zero.
		do {
			it.value++;
			// We keep the following assertion as we
			// must have that at least one entry is non-zero.
			// The read code of the MDCPT should ensure this
			// as sure all parameter update procedures.
			assert (it.value < ucard());
			it.probVal = probGivenParents(it.value);
		} while (it.probVal.essentially_zero());
	}
}


/*-
 *-----------------------------------------------------------------------
 * NGramCPT::becomeAwareOfParentValuesAndIterBegin
 *      Set the parents and set an iterator to the begin.
 *
 * Results:
 *      None.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
void NGramCPT::becomeAwareOfParentValuesAndIterBegin(vector<RandomVariable *>& parents, iterator &it) {
	becomeAwareOfParentValues(parents);
	begin(it);
}


/*-
 *-----------------------------------------------------------------------
 * NGramCPT::next
 *      Advance an iterator.
 *
 * Results:
 *      Return true if the next exists.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
bool NGramCPT::next(iterator &it) {
	do{
		if ( (++it.value) >= ucard() )
			return false;
		it.probVal = probGivenParents(it.value);
	} while ( it.probVal.essentially_zero() );

	return true;
}


/*-
 *-----------------------------------------------------------------------
 * NGramCPT::randomSample
 *      Generate a random sample according to the distribtution.
 *
 * Results:
 *      Return a random sample.
 *
 * Side Effects:
 *      This is using a inline function NGramCPT::end() rather
 *      than the virtual version due to speed reason.
 *      Please be carefult to override it for sub-classes.
 *
 *-----------------------------------------------------------------------
 */
int NGramCPT::randomSample() {
	iterator it;
	begin(it);		// note that it = begin() will vanish.
	logpr prob = rnd.drand48();
	logpr sum = 0.0;

	do {
		sum += it.probVal;
		if ( prob <= sum )
			break;

		next(it);
	} while ( ! NGramCPT::end(it) );

	return it.value;
}


/*-
 *-----------------------------------------------------------------------
 * NGramCPT::read
 *      Read a language model file form ARPA file or index file.
 *      (1) Read a ARPA file
 *              1 % number of parents (order -1)
 *              10 10 % cardinalities
 *              bigram.lm vocab % ARPA lm file and vocab file
 *      (2) Read a index file
 *              1 % number of parents (order -1)
 *              10 10 % cardinalities
 *              bigram.lm.idx [ascii] % ARPA lm file or [binary]
 *      (3) For standard ARPA file format, please refer to other documents.
 *      (4) The index file is for speed reason because there are two passes
 *          in reading an ARPA file.  Refer writeNGramIndexFile for details.
 *
 *
 * Results:
 *      None.
 *
 * Side Effects:
 *      When reading the index file, it will look for bigram.lm.idx.
 *      The index file can be created by gmtkNGramIndex.
 *
 *-----------------------------------------------------------------------
 */
void NGramCPT::read(iDataStreamFile &is) {
	NamedObject::read(is);
	is.read(_numParents,"NGramCPT::read numParents");

	if ( _numParents < 0 ) 
		error("ERROR: reading file '%s', NGramCPT '%s' trying to use negative (%d) num parents.",is.fileName(),name().c_str(),_numParents);
	if ( _numParents >= warningNumParents )
		warning("WARNING: creating NGramCPT '%s' with %d parents in file '%s'", _numParents,name().c_str(),is.fileName());

	cardinalities.resize(_numParents);
	_parentValues.resize(_numParents);
	// read the cardinalities
	for ( unsigned i = 0; i < _numParents; i++ ) {
		is.read(cardinalities[i], "NGramCPT::read cardinality");
		if ( cardinalities[i] <= 0 )
			error("ERROR: reading file '%s', NGramCPT '%s' trying to use 0 or negative (%d) cardinality table, position %d.", is.fileName(),name().c_str(),cardinalities[i],i);
	}

	// read the self cardinalities
	is.read(_card, "NGramCPT::read cardinality");
	if ( _card <= 0 )
		error("ERROR: reading file '%s', NGramCPT '%s' trying to use 0 or negative (%d) cardinality table, position %d.", is.fileName(),name().c_str(),_card,_numParents);

	// read in the string for ARPA language model file name and vocabulary file name
	char *lmFile;
	string vocabName;
	if ( ! is.readStr(lmFile) )
		error("ERROR: reading file '%s', NGramCPT '%s' trying to find ARPA LM filename", is.fileName(), name().c_str());
	if ( ! is.read(vocabName) )
		error("ERROR: reading file '%s', NGramCPT '%s' trying to find vocab name", is.fileName(), name().c_str());

	delete [] _lmIndexFile;
	//_lmIndexFile = new char [strlen(lmFile) + 10];
	_lmIndexFile = new char [strlen(lmFile) + 1];
	strcpy(_lmIndexFile, lmFile);
	if ( strcmp(vocabName.c_str(), "[ascii]") == 0 ) {
		_lmIndexFileBin = false;
		iDataStreamFile sis(_lmIndexFile, false, false);	// ascii no cpp
		readNGramIndexFile(sis);
	} else if ( strcmp(vocabName.c_str(), "[binary]") == 0 ) {
		_lmIndexFileBin = true;
		iDataStreamFile sis(_lmIndexFile, true);			// binary
		readNGramIndexFile(sis);
	} else {
		if ( GM_Parms.vocabsMap.find(vocabName) ==  GM_Parms.vocabsMap.end())
			error("ERROR: reading file '%s', NGramCPT '%s' specifies Vobab name '%s' that does not exist", is.fileName(), _name.c_str(), vocabName.c_str());

		Vocab *vocab = GM_Parms.vocabs[GM_Parms.vocabsMap[vocabName]];
		read(lmFile, *vocab);
	}

	delete [] lmFile;
}


/*-
 *-----------------------------------------------------------------------
 * NGramCPT::read
 *      Read in ngram from ARPA file with coresponding vocab.
 *
 * Results:
 *      Create all the data structure.
 *
 * Side Effects:
 *      The memory at the loading time will be bigger than after.
 *
 *-----------------------------------------------------------------------
 */
void NGramCPT::read(const char *lmFile, const Vocab &vocab) {
	// the goal of this preprocessing is to
	// collect information for hash table sizes
 
	ContextTreeEntry startEntry;

	// phase I
	// 1. Read in the vocabulary
	// 2. Set up the indices table
	int * indices = new int [vocab.size()];
	if ( indices == NULL )
		error("out of memory in NGramCPT::read");

	int i;
	for ( i = 0; i < (int)vocab.size(); ++i )
		indices[i] = i;

	iDataStreamFile ifs(lmFile, false, false);	// ascii, no cpp

	// Read in the ARPA header
	char *line = new char [MAX_LINE_LENGTH];
	// readin '\data\'
	do {
		ifs.readLine(line, MAX_LINE_LENGTH);
	} while ( strstr(line, "\\data\\") == NULL );

	// read in number of ngrams
	unsigned *numNGrams = new unsigned [_numParents + 2];
	if ( ! numNGrams )
		error("out of memory in NGramCPT::read");

	unsigned index, num;
	char c;
	_totalNumberOfParameters = 0;
	for ( i = 1; i < (int)_numParents + 2; ++i ) {
		{
			char *foo;
			if ( ! ifs.readStr(foo) )
				error("ERROR: cannot read string in %s.", lmFile);
			if ( strcmp(foo, "ngram") != 0 )
				error("ERROR: Wrong ARPA format in %s at %s.", lmFile, foo);
			delete [] foo;
		}
		if ( ! ifs.readUnsigned(index) )
			error("ERROR: Wrong ARPA format in %s.", lmFile);
		ifs.readChar(c);
		if ( c != '=' )
			error("ERROR: Wrong ARPA format in %s at %c.", lmFile, c);
		if ( ! ifs.readUnsigned(num) )
			error("ERROR: Wrong ARPA format in %s.", lmFile);
		numNGrams[index] = num;
		_totalNumberOfParameters += num; // After this, we only need to collect bow
	}

	// phase II
	// 1. accumulate the statistics for probability
	//    tables
	// 2. contruct the context Trie for later usage
	unsigned filePos = ifs.ftell();

	unsigned j, k;
	char seps[] = " \t\r\n";
	char *tok;
	int *context = new int [_numParents];
	int wid;

	// step 1: readin unigram as this is the easy case
	// looking for \1-gram
	do {
		ifs.readLine(line, MAX_LINE_LENGTH);
	} while ( line[0] != '\\');
	index = atoi(line+1);

	num = numNGrams[index];	// number of ngrams for index
	startEntry.probOffset = 0;
	unsigned probTotalSize = nextPrime(num);
	startEntry.probBlockSize = probTotalSize;

	for ( j = 0; j < num; ++j ) {
		ifs.readLine(line, MAX_LINE_LENGTH);

		// skip prob (we don't have prob hash table yet.)
		if ( (tok = strtok(line, seps)) == NULL )
			error("error reading line %s", line);

		// read in word
		if ( (tok = strtok(NULL, seps)) == NULL )
			error("error reading line %s", line);
		if ( (wid = vocab.index(tok)) >= (int)vocab.size() )
			error("index %d exceeds vocab size %d in NGramCPT::read", wid, vocab.size());

		ContextTreeEntry contextEntry;
		// read in bow if any
		if ( (tok = strtok(NULL, seps)) != NULL ) {
			// there is bow
			contextEntry.bow = atof(tok) * M_LN10;
			++_totalNumberOfParameters;

			// insert the context into datebase only when there is bow
			if ( startEntry.next == NULL )
				startEntry.next = new vhash_map<int, ContextTreeEntry>(1, num);
			startEntry.next->insert(indices + wid, contextEntry);
			startEntry.keys.push_back(wid);
		}
	}

	// step 2: read in higher order ngrams
	for ( i = 1; i < (int)_numParents; ++i ) {
		// looking for \n-gram
		do {
			ifs.readLine(line, MAX_LINE_LENGTH);
		} while ( line[0] != '\\');
		index = atoi(line+1);

		num = numNGrams[index];	// number of ngrams for index

		ContextTreeEntry *prevContextEntry = NULL;
		unsigned accumNum = 0;
		for ( j = 0; j < num; ++j ) {
			ifs.readLine(line, MAX_LINE_LENGTH);

			// get prob
			if ( (tok = strtok(line, seps)) == NULL )
				error("error reading line %s", line);

			// read in context
			for ( k = 0; k < index - 1; ++k ) {
				if ( (tok = strtok(NULL, seps)) == NULL )
					error("error reading line %s", line);
				if ( (context[k] = vocab.index(tok)) >= (int)vocab.size() )
					error("index %d exceeds vocab size %d in NGramCPT::read", wid, vocab.size());
			}

			// update the previous level of hashing
			ContextTreeEntry *contextEntry = &startEntry;
			for ( int m = index -2; m >= 0; --m )
				contextEntry = contextEntry->next->find(indices + context[m]);

			if ( prevContextEntry == NULL )
				prevContextEntry = contextEntry;	// first on the n-gram
			else if ( prevContextEntry != contextEntry ) {
				// we need to update the numbers in the entry
				prevContextEntry->probBlockSize = nextPrime(accumNum);
				prevContextEntry->probOffset = probTotalSize;
				probTotalSize += prevContextEntry->probBlockSize;
				accumNum = 0;
				prevContextEntry = contextEntry;
			}
			
			++accumNum;

			// now create next level hashing
			// read in word
			if ( (tok = strtok(NULL, seps)) == NULL )
				error("error reading line %s", line);
			if ( (context[index-1] = vocab.index(tok)) >= (int)vocab.size() )
				error("index %d exceeds vocab size %d in NGramCPT::read", wid, vocab.size());

			ContextTreeEntry entry;
			// read in bow if any
			if ( (tok = strtok(NULL, seps)) != NULL ) {
				// there is bow
				entry.bow = atof(tok) * M_LN10;
				++_totalNumberOfParameters;

				// insert the context into datebase only when there is bow
				contextEntry = &startEntry;
				for ( int m = index -1; m > 0; --m ) {
					if ( contextEntry->next == NULL )
						error("wrong lm: lower order context didn't appear before.");
					contextEntry = contextEntry->next->find(indices + context[m]);
				}
				
				if ( contextEntry == NULL )
					error("error in ngram file");
	
				if ( contextEntry->next == NULL )
					contextEntry->next = new vhash_map<int, ContextTreeEntry>(1, 3);
	 
				contextEntry->next->insert(indices + context[0], entry);
				contextEntry->keys.push_back(context[0]);
			}
		}

		// well, we need to worry about the last one
		if ( prevContextEntry ) {
			// we need to update the numbers in the entry
			prevContextEntry->probBlockSize = nextPrime(accumNum);
			prevContextEntry->probOffset = probTotalSize;
			probTotalSize += prevContextEntry->probBlockSize;
		}
	}

	// step 3: read in highest order ngram
	// looking for \n-gram
	if ( _numParents > 0 ) {	// in case of unigram
		do {
			ifs.readLine(line, MAX_LINE_LENGTH);
		} while ( line[0] != '\\');
		index = atoi(line+1);
	
		num = numNGrams[index];	// number of ngrams for index
	
		ContextTreeEntry *prevContextEntry = NULL;
		unsigned accumNum = 0;
		for ( j = 0; j < num; ++j ) {
			ifs.readLine(line, MAX_LINE_LENGTH);
	
			// get prob
			if ( (tok = strtok(line, seps)) == NULL )
				error("error reading line %s", line);
	
			// read in context
			for ( k = 0; k < index - 1; ++k ) {
				if ( (tok = strtok(NULL, seps)) == NULL )
					error("error reading line %s", line);
				if ( (context[k] = vocab.index(tok)) >= (int)vocab.size() )
					error("index %d exceeds vocab size %d in NGramCPT::read", wid, vocab.size());
			}
	
			// skip the rest
	
			ContextTreeEntry *contextEntry = &startEntry;
			for ( int m = index -2; m >= 0; --m )
				contextEntry = contextEntry->next->find(indices + context[m]);
	
			if ( prevContextEntry == NULL )
				prevContextEntry = contextEntry;	// this is the first of the n-gram
			else if ( prevContextEntry != contextEntry ) {
				// we need to update the numbers in the entry
				prevContextEntry->probBlockSize = nextPrime(accumNum);
				prevContextEntry->probOffset = probTotalSize;
				probTotalSize += prevContextEntry->probBlockSize;
				accumNum = 0;
				prevContextEntry = contextEntry;
			}

			++accumNum;
		}
	
		// well, we need to worry about the last one
		if ( prevContextEntry ) {
			// we need to update the numbers in the entry
			prevContextEntry->probBlockSize = nextPrime(accumNum);
			prevContextEntry->probOffset = probTotalSize;
			probTotalSize += prevContextEntry->probBlockSize;
		}
	}

	// reset the buffer for prob table
	_probTable.resize(probTotalSize);

	// step 4: figure out the statistics for context
	// hash tables.  A breadth-first search is needed.
	unsigned contextTotalSize = 0;

	std::queue<ContextTreeEntry *> contextQueue;
	contextQueue.push(&startEntry);
	ContextTreeEntry *currentContextEntry, *nextContextEntry;
	std::vector<int>::const_iterator iit;
	while ( ! contextQueue.empty() ) {
		currentContextEntry = contextQueue.front();
		contextQueue.pop();
		if ( currentContextEntry->next != NULL && currentContextEntry->keys.size() > 0 ) {
			currentContextEntry->contextOffset = contextTotalSize;
			currentContextEntry->contextBlockSize = nextPrime(currentContextEntry->keys.size());
			contextTotalSize += currentContextEntry->contextBlockSize;

			for ( iit = currentContextEntry->keys.begin(); iit != currentContextEntry->keys.end(); ++iit ) {
				nextContextEntry = currentContextEntry->next->find(indices + *iit);
				contextQueue.push(nextContextEntry);
			}
		}
	}

	// phase III:
	// create the hash M table for context
	_contextTable.resize(contextTotalSize);
	_contextStartBlockSize = startEntry.contextBlockSize;

	contextQueue.push(&startEntry);
	while ( ! contextQueue.empty() ) {
		currentContextEntry = contextQueue.front();
		contextQueue.pop();

		if ( currentContextEntry->next != NULL && currentContextEntry->keys.size() > 0 ) {
			for ( iit = currentContextEntry->keys.begin(); iit != currentContextEntry->keys.end(); ++iit ) {
				nextContextEntry = currentContextEntry->next->find(indices + *iit);

				ContextHashEntry hEntry;
				hEntry.bow = nextContextEntry->bow;
				hEntry.probOffset = nextContextEntry->probOffset;
				hEntry.probBlockSize = nextContextEntry->probBlockSize;
				hEntry.nextContextOffset = nextContextEntry->contextOffset;
				hEntry.nextContextBlockSize = nextContextEntry->contextBlockSize;

				_contextTable.insert(*iit, hEntry, currentContextEntry->contextOffset, currentContextEntry->contextBlockSize);

				contextQueue.push(nextContextEntry);
			}
		}
	}

	// phase IV:
	// load the probablities from the file

	// step 1: set the block size and rewind file position
	_probStartBlockSize = startEntry.probBlockSize;
	// rewind the file pointer
	ifs.fseek(filePos-1,SEEK_SET);

	double prob;

	// step 2: load unigram probabilities
	// looking for \n-gram
	do {
		ifs.readLine(line, MAX_LINE_LENGTH);
	} while ( line[0] != '\\');
	index = atoi(line+1);

	num = numNGrams[index];	// number of ngrams for index

	for ( j = 0; j < num; ++j ) {
		ifs.readLine(line, MAX_LINE_LENGTH);

		// get prob
		if ( (tok = strtok(line, seps)) == NULL )
			error("error reading line %s", line);
		prob = atof(tok) * M_LN10;

		// read in context
		for ( k = 0; k < index - 1; ++k ) {
			if ( (tok = strtok(NULL, seps)) == NULL )
				error("error reading line %s", line);
			if ( (context[k] = vocab.index(tok)) >= (int)vocab.size() )
				error("index %d exceeds vocab size %d in NGramCPT::read", wid, vocab.size());
		}

		// read in word
		if ( (tok = strtok(NULL, seps)) == NULL )
			error("error reading line %s", line);
		if ( (wid = vocab.index(tok)) >= (int)vocab.size() )
			error("index %d exceeds vocab size %d in NGramCPT::read", wid, vocab.size());

		_probTable.insert(wid, prob, 0, _probStartBlockSize);
	}

	// step 3: loading higher order probabilities
	for ( i = 1; i < (int)_numParents + 1; ++i ) {
		// looking for \n-gram
		do {
			ifs.readLine(line, MAX_LINE_LENGTH);
		} while ( line[0] != '\\');
		index = atoi(line+1);

		num = numNGrams[index];	// number of ngrams for index

		for ( j = 0; j < num; ++j ) {
			ifs.readLine(line, MAX_LINE_LENGTH);
	
			// get prob
      		if ( (tok = strtok(line, seps)) == NULL )
				error("error reading line %s", line);
			prob = atof(tok) * M_LN10;

			// read in context
			for ( k = 0; k < index - 1; ++k ) {
				if ( (tok = strtok(NULL, seps)) == NULL )
					error("error reading line %s", line);
				if ( (context[k] = vocab.index(tok)) >= (int)vocab.size() )
					error("index %d exceeds vocab size %d in NGramCPT::read", wid, vocab.size());
			}

			// find the context entry
			ContextHashEntry *contextEntry = findContextHashEntry(context, index - 1);
			if ( contextEntry == NULL )
				error("error in ngram file");

			// read in word
			if ( (tok = strtok(NULL, seps)) == NULL )
				error("error reading line %s", line);
			if ( (wid = vocab.index(tok)) >= (int)vocab.size() )
				error("index %d exceeds vocab size %d in NGramCPT::read", wid, vocab.size());

			_probTable.insert(wid, prob, contextEntry->probOffset, contextEntry->probBlockSize);
		}
	}

	delete [] line;
	delete [] numNGrams;
	delete [] context;
	delete [] indices;
}


/*-
 *-----------------------------------------------------------------------
 * NGramCPT::readNGramIndexFile
 *      Read ngram from a index file.
 *
 * Results:
 *      None.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
void NGramCPT::readNGramIndexFile(iDataStreamFile& ifs) {
	ifs.read(_contextStartBlockSize);
	ifs.read(_probStartBlockSize);
	ifs.read(_totalNumberOfParameters);

	// read in table
	_contextTable.readHashMTableForComplexType(ifs);
	_probTable.readHashMTableForBasicType(ifs);
}


/*-
 *-----------------------------------------------------------------------
 * NGramCPT::writeNGramIndexFile
 *      Write ngram to a serialized index file.
 *      The file looks like
 *          context_start_block_size prob_start_block_size
 *          total_number_of_parameters
 *          hash_multi_table_for_context_search
 *          hash_multi_table_for_prob_search
 *
 * Results:
 *      None.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
void NGramCPT::writeNGramIndexFile(oDataStreamFile& ofs) {
	ofs.write(_contextStartBlockSize);
	ofs.write(_probStartBlockSize);
	ofs.write(_totalNumberOfParameters);

	// writing the table
	_contextTable.writeHashMTableForComplexType(ofs);
	_probTable.writeHashMTableForBasicType(ofs);
}


/*-
 *-----------------------------------------------------------------------
 * NGramCPT::ContextHashEntry::readObject
 *      Read context hash entry from a serialized index file.
 *
 * Results:
 *      None.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
void NGramCPT::ContextHashEntry::readObject(iDataStreamFile &ifs) {
	ifs.readDouble(bow);
	ifs.readUnsigned(probOffset);
	ifs.readUnsigned(probBlockSize);
	ifs.readUnsigned(nextContextOffset);
	ifs.readUnsigned(nextContextBlockSize);
}


/*-
 *-----------------------------------------------------------------------
 * NGramCPT::ContextHashEntry::writeObject
 *      Write context hash entry to a serialized index file.
 *
 * Results:
 *      None.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
void NGramCPT::ContextHashEntry::writeObject(oDataStreamFile &ofs) {
	ofs.writeDouble(bow);
	ofs.writeUnsigned(probOffset);
	ofs.writeUnsigned(probBlockSize);
	ofs.writeUnsigned(nextContextOffset);
	ofs.writeUnsigned(nextContextBlockSize);
}
