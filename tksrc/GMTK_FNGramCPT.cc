/*-
 * GMTK_FNGramCPT.cc
 *
 * Written by Gang Ji <gang@ee.washington.edu>
 * Modifications by J. Bilmes. <bilmes@ee.washington.edu>
 *
 * Portions of this code were borrowed from the FLM extensions to the
 * SRI language model toolkit, but only those parts that were written
 * by Jeff Bilmes <bilmes@ee.washington.edu> as part of and during the
 * JHU CLSP 2002 summer workshop.  All such code is public domain.  NO
 * SRI PORTIONS OF SRILM CODE ARE CONTAINED HERE!! THEREFORE, THIS
 * CODE IS NOT UNDER SRI COPYRIGHT.
 *
 * This part of the code has the implementation of class FNGramCPT for
 * gmtk.  Please see GMTK_FNGramCPT.h for more information.
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


#include <cstring>
#include <string>
#include <vector>

#include "GMTK_FNGramCPT.h"
#include "GMTK_DiscRV.h"
#include "GMTK_GMParms.h"
#include "fileParser.h"
#include "rand.h"
#include "vshash_map.h"
#include "shash_map.h"


#define MAX_LINE_LENGTH 1024


// this function is defined in GMTK_Vocab.cc
extern unsigned nextPrime(unsigned n);


/*-
 * This ADT is used only when loading from a ARPA lm file.
 * A Trie ADT is utilized for loading.
 */
struct FNGramContextTreeHashEntry {
	double bow;		// backing-off weight

	// information for the next level context
	// hash table (no iterations for hashmap)
	std::vector<int> wids;
	std::vector<float> probs;

	FNGramContextTreeHashEntry() : bow(0) {
	}
};


struct FNGramContextTreeNode {
	vshash_map<unsigned, FNGramContextTreeHashEntry> *entries;
	std::vector<unsigned*> keys;

	FNGramContextTreeNode() : entries(NULL) {}
	~FNGramContextTreeNode() { delete entries; }
};


/*******************************
 * methods for class FNGramImp *
 *******************************/


/*-
 *-----------------------------------------------------------------------
 * FNGramImp::FNGramImp
 *      Default constructor for factored ngram internal representation
 *
 * Results:
 *      None.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
FNGramImp::FNGramImp() : _numParents(0), _parents(NULL), _numberOfBGNodes(0), _bgNodes(NULL), _countFileName(NULL),
		_lmFileName(NULL), _probabilities(NULL), _counts(NULL), _probStartBlockSize(0) {
}


/*-
 *-----------------------------------------------------------------------
 * FNGramImp::~FNGramImp
 *      Default destructor for factored ngram CPT implementation
 *
 * Results:
 *      None.
 *
 * Side Effects:
 *      Clean up the memory.
 *
 *-----------------------------------------------------------------------
 */
FNGramImp::~FNGramImp() {
	delete [] _parents;
	delete [] _bgNodes;
	delete [] _countFileName;
	delete [] _lmFileName;
	delete _probabilities;
	delete _counts;
}



/*-
 *-----------------------------------------------------------------------
 * FNGramImp::setNumParents
 *      Set number of parents for this factored ngram
 *
 * Results:
 *      None.
 *
 * Side Effects:
 *      Reset some dynamic memory.
 *
 *-----------------------------------------------------------------------
 */
void FNGramImp::setNumParents(const unsigned nParents) {
	if ( nParents > maxNumParentsPerChild )
		error("Error: number of parents %d should be smaller than %d", nParents, maxNumParentsPerChild);

	_numParents = nParents;

	// create buffers
	if ( _parents != NULL )
		delete [] _parents;
	_parents = new ParentType [_numParents];

	cardinalities.resize(_numParents);

	// this is how many possible nodes in total
	_numberOfBGNodes = 1 << _numParents;

	if ( _bgNodes != NULL )
		delete [] _bgNodes;
	_bgNodes = new BackoffGraphNode [_numberOfBGNodes];
}


/*-
 *-----------------------------------------------------------------------
 * FNGramImp::read
 *      Read in the parameters from a master file
 *
 * Results:
 *      None.
 *
 * Side Effectes:
 *      None.
 *-----------------------------------------------------------------------
 */
void FNGramImp::read(iDataStreamFile &is) {
	// step 1:
	// read in information from master file

	// read in name
	NamedObject::read(is);

	// read in number of parents
	int numParents;
	is.read(numParents, "Can't read FNGramCPT's number of parents");
	if ( numParents < 0 )
		error("Error: reading file '%s' line %d, FNGramCPT '%s' trying to use negative (%d) num parents.", is.fileName(),is.lineNo(), name().c_str(), numParents);
	setNumParents(numParents);

	// read the cardinalities of parents
	for ( unsigned i = 0; i < _numParents; i++ ) {
		is.read(cardinalities[i], "Can't read FNGramCPT's parent cardinality");
		if ( cardinalities[i] <= 0 )
			error("Error: reading file '%s' line %d, FNGramCPT '%s' trying to use 0 or negative (%d) cardinality table, position %d.",
				is.fileName(), is.lineNo(),name().c_str(), cardinalities[i], i);
	}

	// read the self cardinalities
	is.read(_card, "Cant' read FNGramCPT self cardinality");
	if ( _card <= 0 )
		error("Error: reading file '%s' line %d, FNGramCPT '%s' trying to use 0 or negative (%d) cardinality table, position %d.",
			is.fileName(), is.lineNo(),name().c_str(), _card, _numParents);

	// read in flm filename
	char* flmFilename = NULL;
	is.readStr(flmFilename, "Can't read FNGramCPT FLM filename");

	// read in vocab mapping
	char* tag = NULL;
	is.readStr(tag, "Can't read FNGramCPT vocab map");
	const char seps[] = "-:";
	char* tok = strtok(tag, seps);
	char tagChar;
	unsigned tagId = 0;
	std::vector<Vocab*> vocabPtrs;
	while ( tok != NULL ) {
		// get the tag char
		tagChar = tok[0];
		if ( _tagMap.find(tagChar) != NULL )
			error("Warning: tag %c was defined before in file %s, 2nd-time at line %d",
			      tagChar,
			      is.fileName(),is.lineNo());

		// read in the vocab object name
		if ( (tok = strtok(NULL, seps)) == NULL )
			error("Error: tag %c should followed by vocab name in %s line %d",
			      tagChar, is.fileName(),is.lineNo());

		std::string vocabName = std::string(tok);

		if ( GM_Parms.vocabsMap.find(vocabName) ==  GM_Parms.vocabsMap.end())
			error("ERROR: reading file '%s' line %d, NGramCPT '%s' specifies Vobab name '%s' that does not exist", is.fileName(),is.lineNo(), _name.c_str(), vocabName.c_str());
		vocabPtrs.push_back(GM_Parms.vocabs[GM_Parms.vocabsMap[vocabName]]);
		_tagMap.insert(tagChar, tagId++);

		tok = strtok(NULL, seps);
	}

	// step 1:
	// read in specification file
	iDataStreamFile spec_ifs(flmFilename, false, false);	// ascii, no cpp
	readFNGramSpec(spec_ifs);

	// step 2:
	// read in language mdoel file
	readLMFile(vocabPtrs);

	// 3. read in count file if necessary
	bool needToReadCounts = false;

	// walk through the valid nodes to decide whether we need counts
	// descend down the BG, level by level except buttom (0)
	for ( int level = _numParents; level > 0; level-- ) {
		LevelIter liter(_numParents, level);
		for ( unsigned nodeAtLevel; liter.next(nodeAtLevel); ) {
			if ( _bgNodes[nodeAtLevel].valid && _bgNodes[nodeAtLevel].requiresGenBackoff() ) {
				switch ( _bgNodes[nodeAtLevel].backoffStrategy ) {
				case CountsNoNorm:
				case CountsSumCountsNorm:
				case CountsSumNumWordsNorm:
				case CountsProdCardinalityNorm:
				case CountsSumCardinalityNorm:
				case CountsSumLogCardinalityNorm:
					needToReadCounts = true;
					break;
				default:
					;
				}
			}

			// get out of the loop
			if ( needToReadCounts )
				break;
		}

		// get out of the loop
		if ( needToReadCounts )
		break;
	}

	if ( needToReadCounts )
		readCountFile(vocabPtrs);

	// 4. compute statisitcs
	computeCardinalityFunctions();

	delete [] flmFilename;
	delete [] tag;
}


/*-
 *-----------------------------------------------------------------------
 * FNGramImp::readFNGramSpec
 *      Read factored ngram specification.
 *
 * Results:
 *      None.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
void FNGramImp::readFNGramSpec(iDataStreamFile &ifs) {
	char c;

	// read in child node name
	ifs.readChar(c);
	unsigned* it = _tagMap.find(c);
	if ( it == NULL )
		error("Error: tag %c not found in %s", c, ifs.fileName());
	_childIndex = *it;

	ifs.readChar(c);
	if ( c != ':' )
		error("Error: expecting ':' in flm file %s", ifs.fileName());

	// read in parents
	unsigned nParents;
	ifs.readUnsigned(nParents);
	if ( nParents != _numParents )
		error("Error: Number of parents %d doesn't match %d in %s", nParents, _numParents, ifs.fileName());
	if ( _numParents > maxNumParentsPerChild )
		error("Error: number of parents %d should be smaller than %d", _numParents, maxNumParentsPerChild);

	unsigned i = 0;
	while ( i < _numParents ) {
		_parents[i].read(ifs, _tagMap);

		i++;
	}

	// read in word counter filename and language model filename
	ifs.readStr(_countFileName);
	ifs.readStr(_lmFileName);

	// read in how many active backing-off graph node is specified
	unsigned numberOfSpecifiedNodes;
	ifs.readUnsigned(numberOfSpecifiedNodes);

	// now read in the nodes
	// this version only supports one line per node
	char *line = new char [2048];
	for ( i = 0; i < numberOfSpecifiedNodes; i++ ) {
		// read the id of the node
		char *str = NULL;
		ifs.readStr(str);

		// node id and mask (constaint) should have information
		// of parents set.  So instead of putting this into
		// type BackingoffGraphNode, I put the code here.
		unsigned nodeId = parseNodeString(str);
		if ( nodeId >= _numberOfBGNodes )
			error("Error: node specifier must be between 0x0 and 0x%x inclusive in %s", _numberOfBGNodes - 1, ifs.fileName());
		delete [] str;
		str = NULL;

		readBackoffGraphNode(ifs, _bgNodes[nodeId], nodeId);
	}

	delete [] line;

	// poste process the reading of flm specifying file
	_bgNodes[_numberOfBGNodes - 1].valid = true;
	_bgNodes[_numberOfBGNodes - 1].order = numBitsSet(_numberOfBGNodes - 1) + 1;
	// descend down the BG, level by level
	for ( int level = _numParents; level >= 0; level-- ) {
		LevelIter liter(_numParents, level);
		bool allAreNull = true;

		for ( unsigned nodeAtLevel; liter.next(nodeAtLevel); ) {
			if ( ! _bgNodes[nodeAtLevel].valid )
				continue;
			allAreNull = false;

			BGChildIter citer(_numParents, nodeAtLevel, _bgNodes[nodeAtLevel].backoffConstraint);
			unsigned int numChildrenUsed = 0;
			for ( unsigned child; citer.next(child); ) {
				_bgNodes[child].valid = true;
				_bgNodes[child].order = numBitsSet(child) + 1;
				numChildrenUsed++;
				// make sure kn-count-parent has counts itself.
				if ( _bgNodes[child].knCountParent != ~0x0u ) {
					const unsigned kncp = _bgNodes[child].knCountParent;
					if ( kncp >= _numberOfBGNodes || (! _bgNodes[kncp].valid)  )
						error("Error: kn-counts-parent argument %d must specify a parent that exists and is in use in %s", kncp, ifs.fileName());
				}
			}

			// everybody must have a child.
			if ( nodeAtLevel > 0 && numChildrenUsed == 0 )
				error("ERROR: backoff graph node 0x%X has no children with backoff constraint 0x%X. Must have at least one child in %s",
					nodeAtLevel, _bgNodes[nodeAtLevel].backoffConstraint, ifs.fileName());
			_bgNodes[nodeAtLevel].numBGChildren = numChildrenUsed;
		}

		if ( allAreNull ) {
			// no count object was created
			// NOTE: we might not want to consider this an error, if we for example
			// want not to backoff to lower levels in backoff graph. In that case,
			// probabilities would become zero, however.
			error("ERROR: backoff constraints leave level %d of backoff graph entirely unexpanded, lower distribution order never reached in %s", level, ifs.fileName());
		}
	}
}


/*-
 *-----------------------------------------------------------------------
 * FNGramImp::readLMFile
 *      Read in the language model probabilities.
 *
 * Results:
 *      None.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
void FNGramImp::readLMFile(std::vector<Vocab*> vocabs) {
	// first, we create all objects
	// descend down the BG, level by level
	// for the last one, we keep it empty because there will be no context.
	for ( int level = _numParents; level > 0; level-- ) {
		LevelIter liter(_numParents, level);
		for ( unsigned nodeAtLevel; liter.next(nodeAtLevel); ) {
			if ( _bgNodes[nodeAtLevel].valid ) {
				// argv size is order - 1
				_bgNodes[nodeAtLevel].contextTable = new vshash_map<unsigned, BackoffGraphNode::HashEntry>(_bgNodes[nodeAtLevel].order - 1, 3);
			}
		}
	}
	_probabilities = new HashMTable<logpr>();

	unsigned* numberOfNGrams = new unsigned [_numberOfBGNodes];
	memset(numberOfNGrams, 0, sizeof(unsigned) * _numberOfBGNodes);
	_totalNumberOfParameters = 0;

	// open the flm file
	iDataStreamFile ifs(_lmFileName, false, false);		// ascii no cpp
	char *line = new char [MAX_LINE_LENGTH];

	// skip anything before '\data\'
	do {
		ifs.readLine(line, MAX_LINE_LENGTH);
	} while ( strstr(line, "\\data\\") != line );

	// read in number of ngrams
	unsigned nodeId;
	unsigned number;
	for ( unsigned i = 0; i < _numberOfBGNodes; i++ ) {
		ifs.readLine(line, MAX_LINE_LENGTH);
		// TODO: fix compiler warning here about "warning: int format, unsigned int arg (arg 4)"
		if ( sscanf(line, "ngram %x=%u", &nodeId, &number) < 2 )
			error("Error: reading lm file %s in %s", _lmFileName, line);
		if ( nodeId > _numberOfBGNodes )
			error("Error: node id %x should be smaller than %x in %s", nodeId, _numberOfBGNodes, _lmFileName);
		numberOfNGrams[nodeId] = number;
	}

	float prob, bow;
	unsigned wid;
	char *tok;
	const char seps[] = " \t\n";

	// phase I
	// read in context information
	// put all information in FNGramContextTreeEntry

	FNGramContextTreeNode* treeNodes = new FNGramContextTreeNode [_numberOfBGNodes];
	std::vector<unsigned*> contextPool;

	std::vector<unsigned> unigramWids;
	std::vector<float> unigramProbs;

	// looking for '\0x0-grams:'
	ifs.readLine(line, MAX_LINE_LENGTH);
	while ( strstr(line, "\\end\\") != line ) {
		if ( sscanf(line, "\\%x-gram:", &nodeId) < 1 )
			error("Error: reading lm file %s in %s", _lmFileName, line);
		if ( nodeId > _numberOfBGNodes )
			error("Error: node id %x should be smaller than %x in %s", nodeId, _numberOfBGNodes, _lmFileName);

		if ( nodeId == 0 ) {
			// unigram
			for ( unsigned i = 0; i < numberOfNGrams[nodeId]; i++ ) {
				ifs.readLine(line, MAX_LINE_LENGTH);

				// get the probability
				tok = strtok(line, seps);
				if ( tok == NULL )
					error("Error: reading flm %s arround %s", _lmFileName, line);
				prob = atof(tok) * M_LN10;

				// read in word
				if ( (tok = strtok(NULL, seps)) == NULL )
					error("Error: reading flm %s arround %s", _lmFileName, line);
				if ( (wid = vocabs[_childIndex]->index(tok)) == (unsigned)vocabs[_childIndex]->index("<unk>")
						&& strcmp(tok, "<unk>") != 0 )
					error("Error: word %s not in vocab", tok);

				unigramWids.push_back(wid);
				unigramProbs.push_back(prob);
				_totalNumberOfParameters++;

				// unlike normal ARPA files, the bow is not here.
			}
		} else if ( numberOfNGrams[nodeId] != 0 ) {	// we do nothing if no probabilities available
			// figure out which bits are on
			std::vector<unsigned> onPos = bitsOn(nodeId);
			if ( onPos.size() != _bgNodes[nodeId].order - 1 )
				error("Error: something is wrong with bitsOn at node %d", nodeId);


			if ( treeNodes[nodeId].entries == NULL ) {
				// need to create new entry
				treeNodes[nodeId].entries = new vshash_map<unsigned, FNGramContextTreeHashEntry>(_bgNodes[nodeId].order - 1, 3);
			}

			// now read in all the ngram for this node
			for ( unsigned i = 0; i < numberOfNGrams[nodeId]; i++ ) {
				ifs.readLine(line, MAX_LINE_LENGTH);

				// get the probability
				tok = strtok(line, seps);
				if ( tok == NULL )
					error("Error: reading flm %s arround %s", _lmFileName, line);
				prob = atof(tok) * M_LN10;

				// read in the context
				// we will not delete the following here
				// this will be deleted in the destructor of FNGramCPT
				// the reason is in vhash_map, only key pointer is copied over
				// in this implementation, i will also reverse the context order
				// say, if the context is a,b,c for d
				// then it will be c,b,a for d
				// the reason is in fngram specification, it is usually in the other order
				// word(-1), word(-2), tag(-1), ...
				unsigned *context = new unsigned [_bgNodes[nodeId].order - 1];
				for ( int j = _bgNodes[nodeId].order - 2; j >= 0; j-- ) {
					if ( (tok = strtok(NULL, seps)) == NULL )
						error("Error: reading flm %s arround %s", _lmFileName, line);
					if ( (context[j] = vocabs[_parents[onPos[j]].index]->index(tok))
							== (unsigned)vocabs[_parents[onPos[j]].index]->index("<unk>")
							&& strcmp(tok, "<unk>") != 0 )
						error("Error: context %s not in vocab with pos %d and vocab index %d", tok, j, _parents[onPos[j]].index);
				}

				FNGramContextTreeHashEntry* entry = treeNodes[nodeId].entries->find(context);
				if ( entry == NULL ) {
					// we need to insert a new one
					entry = treeNodes[nodeId].entries->insert(context, FNGramContextTreeHashEntry());
					treeNodes[nodeId].keys.push_back(context);
					contextPool.push_back(context);
				} else {
					delete [] context;
					context = NULL;
				}

				// read in word
				if ( (tok = strtok(NULL, seps)) == NULL )
					error("Error: reading flm %s arround %s", _lmFileName, line);
				if ( (wid = vocabs[_childIndex]->index(tok)) == (unsigned)vocabs[_childIndex]->index("<unk>") && strcmp(tok, "<unk>") != 0 )
					error("Error: word %s not in vocab", tok);

				// add word id into keys
				entry->wids.push_back(wid);
				entry->probs.push_back(prob);
				_totalNumberOfParameters++;

				// read in bow if exists
				if ( (tok = strtok(NULL, seps)) != NULL ) {
					bow = atof(tok) * M_LN10;
					entry->bow = bow;
					_totalNumberOfParameters++;
				}
			}
		}

		ifs.readLine(line, MAX_LINE_LENGTH);
	}

	// phase II
	// in this phase, we will get all the statistics of nodes

	// step 1:
	// set up the unigram size
	_probStartBlockSize = nextPrime(numberOfNGrams[0]);
	unsigned offSet = _probStartBlockSize;

	// step 2 take care of higher order fngrams
	for ( unsigned i = 1; i < _numberOfBGNodes; i++ ) {
		if ( numberOfNGrams[i] == 0 )
			continue;

		FNGramContextTreeNode& treeNode = treeNodes[i];
		BackoffGraphNode& bgNode = _bgNodes[i];

		if ( bgNode.contextTable == NULL ) {
			bgNode.contextTable = new vshash_map<unsigned, BackoffGraphNode::HashEntry> (bgNode.order - 1, 3);
		}

		for ( std::vector<unsigned*>::const_iterator it = treeNode.keys.begin(); it != treeNode.keys.end(); it++ ) {
			BackoffGraphNode::HashEntry nodeEntry;
			FNGramContextTreeHashEntry *treeEntry = treeNode.entries->find(*it);
			if ( treeEntry == NULL )
				error("Error: somethins is wrong");
			nodeEntry.backingOffWeight.setFromLogP(treeEntry->bow);
			nodeEntry.probTableOffSet = offSet;
			nodeEntry.probTableBlockSize = nextPrime(treeEntry->wids.size());
			offSet += nodeEntry.probTableBlockSize;

			bgNode.contextTable->insert(*it, nodeEntry);
		}
	}

	// phase III
	// insert all probabilities

	// step 1:
	// initial the hash multi-table
	_probabilities->resize(offSet);		// right now, offSet is the total size

	// step 2:
	// the unigram probabilities
	for ( unsigned i = 0; i < unigramWids.size(); i++ ) {
		logpr logProb(NULL, unigramProbs[i]);
		_probabilities->insert(unigramWids[i], logProb, 0, _probStartBlockSize);
	}

	// step 3:
	// higher order fngrams
	for ( unsigned i = 1; i < _numberOfBGNodes; i++ ) {
		if ( numberOfNGrams[i] == 0 )
			continue;

		FNGramContextTreeNode& treeNode = treeNodes[i];
		BackoffGraphNode& bgNode = _bgNodes[i];

		for ( std::vector<unsigned*>::const_iterator it = treeNode.keys.begin(); it != treeNode.keys.end(); it++ ) {
			FNGramContextTreeHashEntry *treeEntry = treeNode.entries->find(*it);
			if ( treeEntry == NULL )
				error("Error: somethins is wrong");

			BackoffGraphNode::HashEntry *nodeEntry = bgNode.contextTable->find(*it);
			if ( nodeEntry == NULL )
				error("Error: somethins is wrong");

			for ( unsigned j = 0; j < treeEntry->wids.size(); j++ ) {
				logpr logProb(NULL, treeEntry->probs[j]);
				_probabilities->insert(treeEntry->wids[j], logProb, nodeEntry->probTableOffSet, nodeEntry->probTableBlockSize);
			}
		}
	}

	delete [] line;
	delete [] numberOfNGrams;
	delete [] treeNodes;

	for ( std::vector<unsigned*>::iterator it = contextPool.begin(); it != contextPool.end(); it++ ) {
		delete [] *it;
		*it = NULL;
	}
}


/*-
 *-----------------------------------------------------------------------
 * FNGramImp::readCountFile
 *      Read count file (after loading language models).
 *
 * Results:
 *      None.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
void FNGramImp::readCountFile(std::vector<Vocab*> vocabs) {
	assert(_probStartBlockSize > 0 );	// make sure language model is already loaded

	if ( _counts == NULL )
		_counts = new HashMTable<unsigned>();
	_counts->resize(_probabilities->tableSize());

	unsigned nodeId;
	unsigned wid;
	unsigned* context = new unsigned [_numParents];
	unsigned count;

	char* line = new char [MAX_LINE_LENGTH];
	char* tok;
	const char seps[] = " \t\n";

	// open the file in ascii and no cpp mode
	iDataStreamFile ifs(_countFileName, false, false);

	while ( ifs.readLine(line, MAX_LINE_LENGTH) ) {
		if ( (tok = strtok(line, seps)) == NULL ) {
			// empty line
			continue;
		}

		// the first is the index of the node
		sscanf(tok, "%x", &nodeId);
		if ( nodeId >= _numberOfBGNodes )
			error("Error: node id %x should be smaller than %x in count file %s", nodeId, _numberOfBGNodes, _countFileName);

		if ( nodeId == 0 ) {
			// unigram
			if ( (tok = strtok(NULL, seps)) == NULL )
				error("Error: expecting word in %s in %s", line, _countFileName);
			if ( (wid = vocabs[_childIndex]->index(tok)) == (unsigned)vocabs[_childIndex]->index("<unk>") && strcmp(tok, "<unk>") != 0 )
				error("Error: word %s not in vocab", tok);

			if ( (tok = strtok(NULL, seps)) == NULL )
				error("Error: expecting count in %s in %s", line, _countFileName);
			count = atoi(tok);

			_counts->insert(wid, count, 0, _probStartBlockSize);
		} else {
			// we need to make sure that the node is valid
			if ( ! _bgNodes[nodeId].valid )
				continue;		// we skip the count

			// figure out what are the on bits so that vocab can find their word ids
			std::vector<unsigned> onPos = bitsOn(nodeId);

			// read in context
			// remember that we need to reverse the order of context
			for ( int i = _bgNodes[nodeId].order - 2; i >= 0; i-- ) {
				if ( (tok = strtok(NULL, seps)) == NULL )
					error("Error: expecting context in %s in %s", line, _countFileName);
				if ( (context[i] = vocabs[_parents[onPos[i]].index]->index(tok))
						== (unsigned)vocabs[_parents[onPos[i]].index]->index("<unk>")
						&& strcmp(tok, "<unk>") != 0 )
					error("Error: context %s not in vocab with pos %d and vocab index %d at node %x", tok, i, _parents[onPos[i]].index, nodeId);
			}

			// find the context hash entry
			BackoffGraphNode::HashEntry *nodeEntry = _bgNodes[nodeId].contextTable->find(context);
			if ( nodeEntry == NULL )
				continue;		// I found these cases exist. strange though.

			// read in the word
			if ( (tok = strtok(NULL, seps)) == NULL )
				error("Error: expecting word in %s in %s", line, _countFileName);
			if ( (wid = vocabs[_childIndex]->index(tok)) == (unsigned)vocabs[_childIndex]->index("<unk>") && strcmp(tok, "<unk>") != 0 )
				error("Error: word %s not in vocab", tok);

			// well, if the probability doesn't exit, we don't insert
			// otherwise, the hash table will overful. ( infinite loop!)
			if ( _probabilities->find(wid, nodeEntry->probTableOffSet, nodeEntry->probTableBlockSize) == NULL )
				continue;

			// read in count
			if ( (tok = strtok(NULL, seps)) == NULL )
				error("Error: expecting count in %s in %s", line, _countFileName);
			count = atoi(tok);

			// insert the count
			_counts->insert(wid, count, nodeEntry->probTableOffSet, nodeEntry->probTableBlockSize);
		}
	}

	delete [] context;
	delete [] line;
}


/*-
 *-----------------------------------------------------------------------
 * FNGramImp::strtolplusb
 *      John Henderson's extension to strtol with 0b prefix.
 *
 * Results:
 *      Integer correpoding to the bits.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
long FNGramImp::strtolplusb(const char *nptr, char **endptr, int base) {
	const char *i;
	long sign;

	/**
	 * We should only try to be clever if the base 2 is specified, or no
	 * base is specified, just like in the 0x case.
	 */
	if ( base != 2 && base != 0)
		return strtol(nptr, endptr, base);

	i = nptr;

	/* skip white space */
	while ( *i && isspace(*i) )
		i++;

	/* decide what the sign should be */
	sign = 1;
	if ( *i ) {
		if ( *i == '+' ) {
			i++;
		} else if ( *i == '-' ) {
			sign = -1;
			i++;
		}
	}

	/**
	 * If we're not at the end, and we're "0b" or "0B", then return the result
     * base 2.  Let strtol do the work for us.
	 */
	if ( *i && *i == '0' && *(i+1) && (*(i+1) == 'b' || *(i+1)=='B') && *(i+2) && (*(i+2)=='1' || *(i+2)=='0') )
		return sign*strtol(i + 2, endptr, 2);
	else
		/* otherwise use the given base */
		return strtol(nptr, endptr, base);
}


/*-
 *-----------------------------------------------------------------------
 * FNGramImp::numBitsSet
 *      Calcualte how many bits are on.
 *
 * Results:
 *      Return the number of 1's.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
unsigned FNGramImp::numBitsSet(unsigned u) {
	unsigned count = 0;

	while ( u ) {
		count += (u&0x1);
		u >>= 1;
	}

	return count;
}


const std::vector<unsigned> FNGramImp::bitsOn(unsigned n) {
	std::vector<unsigned> ons;
	unsigned pos = 0;
	while ( n != 0 ) {
		if ( n & 0x1u )
			ons.push_back(pos);
		n >>= 1;
		pos++;
	}

	return ons;
}


/*-
 *-----------------------------------------------------------------------
 * FNGramImp::parseNodeString
 *      Read a node string such as 'M1,M0,S1,S2,S+2'
 *
 * Results:
 *      Bits id.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
unsigned FNGramImp::parseNodeString(char *str) {
	if ( str == NULL )
		return 0;

	unsigned bits = 0;

	char *p = str;

	// skip white spaces which should not happen due to prepareNext
	while ( *p != '\0' && isspace(*p) ) {
		p++;
	}

	// get parent
	char *endptr;
	endptr = str;
	bits = (unsigned) strtolplusb(str, &endptr, 0);
	if ( endptr != str ) {
		return bits;
	}

	while ( *p ) {
		// parse tokens of the form TAGNUMBER,TAGNUMBER,TAGNUMBER
		// where TAG is one of the parent names, and number
		// is the negative of the parent position. For example,
		// given parents of the form
		//    W : 3 M(-1) M(-2) M(0) S(-1) S(-2) S(+2)
		// a valid string would be
		//     M1,M0,S1,S2,S+2
		// and which could correspond to parents
		//
		// M(-1), M(0), S(-1), S(-2), S(+2)
		//
		// and be bit vector
		//      0b11101
		// which is returned.

		ParentType parent;
		unsigned* it = _tagMap.find(*p);
		if ( it == NULL )
			error("Error: parent '%c' not exist", *p);
		parent.index = *it;
		p++;

		// be careful about the sign
		// if we see W1 it means W(-1)
		int sign = -1;
		if ( *p == '-' ) {
			sign = -1;
			p++;
		} else if ( *p == '+' ) {
			sign = +1;
			p++;
		}

		// get absolute shift
		parent.offset = strtol(p, &endptr, 0) * sign;
		if ( endptr == p ) {
			error("Can't form integer at parent specifier in string (%s)", str);
		}
		p = endptr;

		// search for parent and position.
		unsigned i;
		for ( i = 0; i < _numParents; i++ ) {
			if ( _parents[i] == parent ) {
				// found
				if ( bits & (1 << i) ) {
					// already set, might be an oversite or error by user, give warning
					warning("WARNING: parent specifier given twice in string (%s)", str);
				}
				bits |= (1 << i);
				break;
			}
		}
		if ( i == _numParents ) {
			error("Can't find a valid parent specifier in string (%s)", str);
		}

		if ( *p == ',' )
			p++;
	}

	return bits;
}


/*-
 *-----------------------------------------------------------------------
 * FNGramImp::readBackoffGraphNode
 *      Read the node properties from a file stream.
 *
 * Results:
 *      None.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
void FNGramImp::readBackoffGraphNode(iDataStreamFile &ifs, BackoffGraphNode &bgNode, unsigned nodeId) {
	// line should have the form
	// NODE_NUM BACKOFFCONSTRAINT <options>
	// options include what ngram-count.cc uses on comand line for
	// discount options. I've given up on any (even slighty) fancy
	// C string parsing for now, so this is just a string of tokens which
	// are parsed in a very simple way.
	// TODO: do a proper multi-line tokenizer here.

	//
	// Current set of Node Options
	//
	// gtmin [num]
	// gtmax [num]
	// gt [fileName string]
	// cdiscount [double]
	// ndiscount []
	// wbdiscount []
	// kndiscount []
	// ukndiscount []
	// kn-counts-modified []
	// kn-counts-modify-at-end []
	// kn [fileName string]
	// interpolate []
	// write [fileName string]
	// strategy [option]
	//    where [option] is one of:
	//            counts_no_norm
	//            counts_sum_counts_norm
	//            counts_sum_num_words_norm
	//            counts_prod_card_norm
	//            counts_sum_card_norm
	//            counts_sum_log_card_norm
	//            bog_node_prob
	//

	char *line = new char [2048];
	ifs.readLine(line, 2048);
	const char seps[] = " \t\n";

	// get bit mask (constraint)
	char *tok = strtok(line, seps);
	bgNode.backoffConstraint = parseNodeString(tok);

	// get other options
	while ( (tok = strtok(NULL, seps)) != NULL ) {
		if ( strcmp(tok, "gtmin") == 0 ) {
			if ( (tok = strtok(NULL, seps)) == NULL )
				error("Error: gtmin argument needs a value reading factored spec file in %s", line);
			if ( sscanf(tok, "%u", &bgNode.gtmin) != 1 )
				error("Error: gtmin argument needs integer value in %s", line);
		} else if ( strcmp(tok, "gtmax") == 0 ) {
			if ( (tok = strtok(NULL, seps)) == NULL )
				error("Error: gtmax argument needs a value reading factored spec file in %s", line);
			if ( sscanf(tok, "%u", &bgNode.gtmax) != 1 )
				error("Error: gtmax argument needs integer value in %s", line);
		} else if ( strcmp(tok, "gt") == 0 ) {
			if ( (tok = strtok(NULL, seps)) == NULL )
				error("Error: gt argument needs a value reading factored spec file in %s", line);
			delete [] bgNode.gtFile;
			bgNode.gtFile = copyToNewStr(tok);
		} else if ( strcmp(tok, "cdiscount") == 0 ) {
			if ( (tok = strtok(NULL, seps)) == NULL )
				error("Error: cdiscount argument needs a value in %s", line);
			double tmp;
			char *endptr;
			tmp = strtod(tok, &endptr);
			if ( endptr == tok )
				error("Error: cdiscount argument (%s) should be floating point value %d in %s", tok, line);
			bgNode.cdiscount = tmp;
		} else if ( strcmp(tok, "ndiscount") == 0 ) {
			bgNode.ndiscount = true;
		} else if ( strcmp(tok, "wbdiscount" ) == 0 ) {
			bgNode.wbdiscount = true;
		} else if ( strcmp(tok, "kndiscount") == 0 ) {
			bgNode.kndiscount = true;
		} else if ( strcmp(tok, "ukndiscount") == 0 ) {
			bgNode.ukndiscount = true;
		} else if ( strcmp(tok, "kn-counts-modified") == 0 ) {
			bgNode.knCountsModified = true;
		} else if ( strcmp(tok, "kn-counts-modify-at-end" ) == 0) {
			bgNode.knCountsModifyAtEnd= true;
		} else if ( strcmp(tok, "kn-count-parent") == 0 ) {
			if ( (tok = strtok(NULL, seps)) == NULL )
				error("Error: kn-count-parent argument needs a parent specifier in %s", line);
			bgNode.knCountParent = parseNodeString(tok);
		} else if ( strcmp(tok, "kn") == 0 ) {
			if ( (tok = strtok(NULL, seps)) == NULL )
				error("Error: kn argument needs a value reading factored spec file in %s", line);
			delete [] bgNode.knFile;
			bgNode.knFile = copyToNewStr(tok);
		} else if ( strcmp(tok, "interpolate") == 0 ) {
			bgNode.interpolate = true;
		} else if ( strcmp(tok, "write") == 0 ) {
			if ( (tok = strtok(NULL, seps)) == NULL )
				error("Error: write argument needs a value reading factored spec file in %s", line);
			delete [] bgNode.countFile;
			bgNode.countFile = copyToNewStr(tok);
		} else if ( strcmp(tok, "strategy") == 0 ) {
			if ( (tok = strtok(NULL, seps)) == NULL )
				error("Error: strategy argument needs a value reading factored spec file in %s", line);
			if ( strcmp(tok, "counts_no_norm") == 0 ) {
				bgNode.backoffStrategy = CountsNoNorm;
			} else if ( strcmp(tok, "counts_sum_counts_norm") == 0 ) {
				bgNode.backoffStrategy = CountsSumCountsNorm;
			} else if ( strcmp(tok, "counts_sum_num_words_norm") == 0 ) {
				bgNode.backoffStrategy = CountsSumNumWordsNorm;
			} else if ( strcmp(tok, "counts_prod_card_norm") == 0 ) {
				bgNode.backoffStrategy = CountsProdCardinalityNorm;
			} else if ( strcmp(tok, "counts_sum_card_norm") == 0 ) {
				bgNode.backoffStrategy = CountsSumCardinalityNorm;
			} else if ( strcmp(tok, "counts_sum_log_card_norm") == 0 ) {
				bgNode.backoffStrategy = CountsSumLogCardinalityNorm;
			} else if ( strcmp(tok, "bog_node_prob") == 0 ) {
				bgNode.backoffStrategy = BogNodeProb;
			} else
				error("Error: unknown strategy argument (%s) when reading factored spec file", tok);
		} else if ( strcmp(tok, "combine") == 0 ) {
			if ( (tok = strtok(NULL, seps)) == NULL )
				error("Error: combine argument needs a value reading factored spec file in %s", line);
			if ( strcmp(tok, "max") == 0 ) {
				bgNode.backoffCombine = MaxBgChild;
			} else if ( strcmp(tok, "min") == 0 ) {
				bgNode.backoffCombine = MinBgChild;
			} else if ( (strcmp(tok, "avg") == 0) || (strcmp(tok, "mean") == 0) ) {
				bgNode.backoffCombine = AvgBgChild;
			} else if ( strcmp(tok, "wmean") == 0 ) {
				bgNode.backoffCombine = WmeanBgChild;
				// next set of tokens must have a combination of (node_spec, weight)
				// for each child.
				// we compute this below, but we compute it here since we don't
				// have the quantity numBGchildren yet.
				BGChildIter citer(_numParents, nodeId, bgNode.backoffConstraint);
				unsigned int numChildrenUsed = 0;
				for ( unsigned child; citer.next(child); ) {
					if ( ~child & (bgNode.backoffConstraint & nodeId)) {
						numChildrenUsed++;
					}
				}

				logpr *wmean = new logpr[numChildrenUsed];
				for ( unsigned cnum = 0; cnum < numChildrenUsed; cnum++ ) {
					wmean[cnum] = 0.0;
				}
				for ( unsigned cnum = 0; cnum < numChildrenUsed; cnum++ ) {
					double value;
					unsigned childSpec = parseNodeString(tok);

					if ( (tok = strtok(NULL, seps)) == NULL )
						error("Error: combine wmean needs a value for weight in %s", line);

					char *endptr;
					value = strtod(tok, &endptr);
					if ( endptr == tok || value < 0.0 )
						error("Error: combine wmean invalid weight value in %s", line);

					citer.init();
					unsigned int cpos = 0;
					for ( unsigned child; citer.next(child); ) {
						if ( ! (~child & (bgNode.backoffConstraint & nodeId)))
							continue;
						if ( child == childSpec )
							break;
						cpos++;
					}
					if ( cpos == numChildrenUsed )
						error("Error: combine wmean, invalid child node given in %s", line);

					// load them in the array in the order that they will
					// be encountered when doing a child iter.
					wmean[cpos].setFromP(value);

					// read the next
					if ( (tok = strtok(NULL, seps)) == NULL )
						error("Error: combine wmean needs more node/value in %s", line);
				}
				logpr sum;
				for ( unsigned cnum = 0; cnum < numChildrenUsed; cnum++ ) {
					sum += wmean[cnum];
				}
				// normalize and convert to logp
				for ( unsigned cnum = 0; cnum < numChildrenUsed; cnum++ ) {
					wmean[cnum] = wmean[cnum] / sum;
				}
				bgNode.wmean = wmean;
			} else if ( strcmp(tok, "sum") == 0 ) {
				bgNode.backoffCombine = SumBgChild;
			} else if ( strcmp(tok, "prod") == 0 ) {
				bgNode.backoffCombine = ProdBgChild;
			} else if ( strcmp(tok, "gmean") == 0 ) {
				bgNode.backoffCombine = GmeanBgChild;
			} else
				error("Error: unknown combine argument (%s) when reading factored spec file", tok);
		} else
			error("Error: unknown argument (%s) when reading factored spec file", tok);
	}

	delete [] line;
}


/*-
 *-----------------------------------------------------------------------
 * FNGramImp::computeCardinalityFunctions
 *      compute statics for calculating scores.
 *
 * Results:
 *      None.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
void FNGramImp::computeCardinalityFunctions() {
	assert(_numParents == cardinalities.size());

	// descend down the BG, level by level
	for ( int level = _numParents; level >= 0; level-- ) {
		LevelIter liter(_numParents, level);
		for ( unsigned nodeAtLevel; liter.next(nodeAtLevel); ) {
			if ( _bgNodes[nodeAtLevel].valid ) {
				// for child
				_bgNodes[nodeAtLevel].prodCardinalities = _bgNodes[nodeAtLevel].sumCardinalities = (double)_card;
				_bgNodes[nodeAtLevel].sumLogCardinalities = log((double)_card);

				// for parents
				std::vector<unsigned> pos = bitsOn(nodeAtLevel);
				for ( std::vector<unsigned>::const_iterator it = pos.begin(); it != pos.end(); it++ ) {
					_bgNodes[nodeAtLevel].prodCardinalities *= cardinalities[*it];
					_bgNodes[nodeAtLevel].sumCardinalities += cardinalities[*it];
					_bgNodes[nodeAtLevel].sumLogCardinalities += log((double)cardinalities[*it]);
				}
			}
		}
	}
}


/*-
 *-----------------------------------------------------------------------
 * FNGramImp::probBackingOff
 *      calculate the probability with backing-off support
 *
 * Results:
 *      Return the probability.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
inline logpr FNGramImp::probBackingOff(unsigned val, unsigned nodeId, BackoffGraphNode::HashEntry **contextEntries, const std::vector<unsigned> &parentsValues) {
	logpr* probPtr;

	if ( nodeId == 0 ) {
		// well, we are at the unigram node
		if ( (probPtr = _probabilities->find(val, 0, _probStartBlockSize)) != NULL )
			return *probPtr;
		return logpr(0.0);	// return zero prob
	}

	// at higher level nodes
	BackoffGraphNode::HashEntry * contextEntry = contextEntries[nodeId];
	if ( contextEntry != NULL && (probPtr = _probabilities->find(val, contextEntry->probTableOffSet, contextEntry->probTableBlockSize)) != NULL ) {
		// easy job, just return the value
		return *probPtr;
	}

	// now we need to back-off
	if ( (contextEntry != NULL) && (! contextEntry->backingOffWeight.essentially_zero()) ) {
		// we already have backing-off weight
		// back-off directly
		return contextEntry->backingOffWeight * bgChildProbBackingOff(val, nodeId, contextEntries, parentsValues);
	}

	// we don't have backing-off weight available
	// do we need generalized backing-off?
	if ( ! _bgNodes[nodeId].requiresGenBackoff() ) {
		return bgChildProbBackingOff(val, nodeId, contextEntries, parentsValues);
	}

	// now we really need to compute the backing-off weights
	logpr sum = bgChildProbSum(nodeId, contextEntries, parentsValues);

	/*
	right now, the following is not supported.
	1. we need to change becomeAwareOfParents so that NULL context entry will not appear
	2. we might need to keep a copy of parent value so that contex tentry can be
	   created on the fly
	3. we might keep a cache of this because in the decoding, the memory will grow too big.
	// for future usage of the same context, we store the value
	if ( contextEntry != NULL )
		contextEntry->backingOffWeight.setFromLogP(- sum.val());
	*/
#ifdef FNGRAM_BOW_GROW_DYNA
	// in this case, we will store the calculated bow in the context table for later usage
	if ( contextEntries[nodeId] == NULL ) {
		// construct the context value
		unsigned *context = new unsigned [_bgNodes[nodeId].order - 1];
		std::vector<unsigned> onPos = bitsOn(nodeId);	// this will be something like 0, 1, 3
		for ( unsigned i = 0; i < _bgNodes[nodeId].order - 1; i++ ) {
			context[i] = parentsValues[onPos[i]];
		}

		contextEntries[nodeId] = _bgNodes[nodeId].contextTable->insert(context, BackoffGraphNode::HashEntry());
	}

	contextEntries[nodeId]->backingOffWeight.setFromLogP(-sum.val());
#endif

	return bgChildProbBackingOff(val, nodeId, contextEntries, parentsValues) / sum;
}


/*-
 *-----------------------------------------------------------------------
 * FNGramImp::bgChildProbBackingOff
 *      general graph backoff algorithm for multiple BG children.
 *
 * Results:
 *      Return the probability.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
inline logpr FNGramImp::bgChildProbBackingOff(unsigned val, unsigned nodeId, BackoffGraphNode::HashEntry **contextEntries, const std::vector<unsigned> &parentsValues) {
	// create a reference to speed up
	BackoffGraphNode& bgNode = _bgNodes[nodeId];

	// special case numBGChildren == 1 for speed.
	if ( bgNode.numBGChildren == 1 ) {
		// Jeff was using iterators here.  But I don't think its necessary
		unsigned bg_child = nodeId & (~bgNode.backoffConstraint);
		if ( _bgNodes[bg_child].valid )
			return probBackingOff(val, bg_child, contextEntries, parentsValues);
	}

	// still here? Do the general case.
	logpr bo_prob;
	if ( bgNode.backoffCombine == ProdBgChild || bgNode.backoffCombine == GmeanBgChild ) {
		bo_prob.set_to_one();
		unsigned bg_child;
		BGChildIter citer(_numParents, nodeId, bgNode.backoffConstraint);
		while ( citer.next(bg_child) ) {
			if ( _bgNodes[bg_child].valid ) {
				// multiply the probs (add the log probs)
				bo_prob *= probBackingOff(val, bg_child, contextEntries, parentsValues);
			}
		}

		if ( bgNode.backoffCombine == GmeanBgChild )
			bo_prob.setFromLogP(bo_prob.val() / bgNode.numBGChildren);
	} else if ( bgNode.backoffCombine == SumBgChild || bgNode.backoffCombine == AvgBgChild) {
		bo_prob.set_to_zero();
		unsigned bg_child;
		BGChildIter citer(_numParents, nodeId, bgNode.backoffConstraint);
		while ( citer.next(bg_child) ) {
			if ( _bgNodes[bg_child].valid ) {
				// add the probs
				bo_prob += probBackingOff(val, bg_child, contextEntries, parentsValues);
			}
		}
		if ( bgNode.backoffCombine == AvgBgChild )
			bo_prob /= (double)bgNode.numBGChildren;
	} else if ( bgNode.backoffCombine == WmeanBgChild ) {
		bo_prob.set_to_zero();
		unsigned bg_child;
		BGChildIter citer(_numParents, nodeId, bgNode.backoffConstraint);
		unsigned cpos = 0;
		while ( citer.next(bg_child) ) {
			if ( _bgNodes[bg_child].valid ) {
				// add the probs by weights
				bo_prob += bgNode.wmean[cpos] * probBackingOff(val, bg_child, contextEntries, parentsValues);
				cpos++;
			}
		}
	} else {
		// choose only one backoff node
		unsigned chosen_descendant = boNode(val, nodeId, contextEntries, parentsValues);
		bo_prob = probBackingOff(val, chosen_descendant, contextEntries, parentsValues);
	}

	return bo_prob;
}


/*-
 *-----------------------------------------------------------------------
 * FNGramImp::bgChildProbSum
 *
 * Results:
 *      Return the sum of all sub nodes probabilities.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
inline logpr FNGramImp::bgChildProbSum(unsigned nodeId, BackoffGraphNode::HashEntry **contextEntries, const std::vector<unsigned> &parentsValues) {
	logpr total;

	for ( unsigned i = 0; i < _card; i++ ) {
		total += bgChildProbBackingOff(i, nodeId, contextEntries, parentsValues);
	}

	return total;
}


/*-
 *-----------------------------------------------------------------------
 * FNGramImp::boNode
 * For a given node in the BG, compute the child node that
 * we should backoff to. The 'context' argument is assumed
 * to be in reverse order with respect to the cound tries.
 *
 * Results:
 *      selected node.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
unsigned FNGramImp::boNode(unsigned val, unsigned nodeId, BackoffGraphNode::HashEntry **contextEntries, const std::vector<unsigned> &parentsValues) {
	// all nodes with one parent back off to the unigram.
	const unsigned nbitsSet = numBitsSet(nodeId);
	if ( nbitsSet == 1 )
		return 0;

	// create a reference to speed up
	BackoffGraphNode& bgNode = _bgNodes[nodeId];

	BGChildIter citer(_numParents, nodeId, bgNode.backoffConstraint);

	unsigned bg_child;						// backoff-graph child
	unsigned chosen_bg_child = ~0x0u;
	unsigned a_bg_child = ~0x0u;					// arbitrary child, if everything fails, we use this
	const bool domax = (bgNode.backoffCombine == MaxBgChild);	// do min if domax == false.
	const double initScore = domax ? -1e220 : 1e220;
	double bestScore = initScore;

	// find child node with the largest counts
	while ( citer.next(bg_child) ) {
		// We max/min over the BO value. The BO value is determined by the
		// BO algorithm (i.e., counts, normalized counts, etc.)  in the
		// specs object for this node.

		double score = backoffValueRSubCtxW(val, bgNode.backoffStrategy, bg_child, contextEntries, parentsValues);

		if ( score == -1e200 )	// TODO: change this to NaN or Inf, and a #define (also see below)
			continue;	// continue presumably because of a NULL counts object
		if ( a_bg_child == ~0x0u )
			a_bg_child = bg_child;

		if ( domax && (score > bestScore) || ! domax && (score < bestScore) ) {
			chosen_bg_child = bg_child;
			bestScore = score;
		}
	}

	// make sure that we have at least one valid child
	assert ( a_bg_child != ~0x0u );

	// if we only have one child, or if we have two children and have
	// not found a best child, just return an arbitrary child node.
	if ( (bgNode.numBGChildren == 1) || (chosen_bg_child == ~0x0u && nbitsSet == 2) )
		return a_bg_child;

	if ( chosen_bg_child == ~0x0u ) {
		// Then we did not found any BG-child with a score for this
		// context. We back off to the child that has the best
		// combined score of its children. We keep
		// doing this, as long as possible.

		unsigned great;
		for ( great = 0; (great + 2) < nbitsSet ; great++ ) {
			citer.init();
			while ( citer.next(bg_child) ) {
				double score = initScore;
				// get score for child
				if ( great == 0 ) {
					// children of child iter
					BGChildIter gciter(_numParents, bg_child, _bgNodes[bg_child].backoffConstraint);
					unsigned bg_grandchild;
					while ( gciter.next(bg_grandchild) ) {
						double tmp = backoffValueRSubCtxW(val, bgNode.backoffStrategy, bg_grandchild, contextEntries, parentsValues);
						// compute local max min of offspring
						if ( domax && (tmp > score) || ! domax && (tmp < score) ) {
							score = tmp;
						}
					}
				} else {
					// grandchildren of child iter
					BGGrandChildIter descendant_iter(_numParents, bg_child, _bgNodes, great - 1);
					unsigned bg_grandchild;
					while ( descendant_iter.next(bg_grandchild) ) {
						double tmp = backoffValueRSubCtxW(val, bgNode.backoffStrategy, bg_grandchild, contextEntries, parentsValues);
						if ( domax && (tmp > score) || ! domax && (tmp < score) ) {
							score = tmp;
						}
					}
				}

				if ( score == -1e200 )		// TODO: change this to NaN or Inf, and a #define (also see above)
					continue;		// presumably because of a NULL counts objects
				if ( domax && (score > bestScore) || ! domax && (score < bestScore) ) {
					chosen_bg_child = bg_child;
					bestScore = score;
				}
			}

			if ( chosen_bg_child != ~0x0u )
				break;
		}

		// still not found, chose an arbitrary child
		if ( chosen_bg_child == ~0x0u )
			chosen_bg_child = a_bg_child;
	}

	return chosen_bg_child;
}


/*-
 *-----------------------------------------------------------------------
 * FNGramImp::backoffValueRSubCtxW
 *      general backoff node selection strategy.
 *
 * Results:
 *      score for the child node.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
double FNGramImp::backoffValueRSubCtxW(unsigned val, BackoffNodeStrategy parentsStrategy, unsigned nodeId, BackoffGraphNode::HashEntry **contextEntries, const std::vector<unsigned> &parentsValues) {
	// create a reference to speed up
	BackoffGraphNode& bgNode = _bgNodes[nodeId];

	// we should never select a path with null counts
	// this could occur because of backoff constraints
	if ( ! bgNode.valid )
		return -1e200;		// TODO: return NaN or Inf to signal this condition (and use #define)

	switch ( parentsStrategy ) {
	case CountsNoNorm: {
			if ( nodeId != 0 && contextEntries[nodeId] == NULL )
				return 0.0;
			unsigned *cnt;
			if ( nodeId == 0 )
				cnt = _counts->find(val, 0, _probStartBlockSize);
			else
				cnt = _counts->find(val, contextEntries[nodeId]->probTableOffSet, contextEntries[nodeId]->probTableBlockSize);
			if ( cnt )
				return (double)*cnt;
			else
				return 0.0;
		}
		break;
	case CountsSumCountsNorm: {
			// choose the node with the max normalized counts,
			// this is equivalent to choosing the largest
			// maximum-likelihod (ML) prob (a max MI criterion).
			if ( nodeId != 0 && contextEntries[nodeId] == NULL )
				return 0.0;
			unsigned *cnt;
			if ( nodeId == 0 )
				cnt = _counts->find(val, 0, _probStartBlockSize);
			else
				cnt = _counts->find(val, contextEntries[nodeId]->probTableOffSet, contextEntries[nodeId]->probTableBlockSize);
			if ( ! cnt )
				return 0.0;

			// get the counts of all entries
			double denominator = 0;
			unsigned key, count;
			HashMTable<unsigned>::iterator it(*_counts, contextEntries[nodeId]->probTableOffSet, contextEntries[nodeId]->probTableBlockSize);
			while ( it.next(key, count) ) {
				denominator += count;
			}

			return ((double)(*cnt)) / denominator;
		}
		break;
	case CountsSumNumWordsNorm: {
			// normalize by the number of words that have occured.
			if ( nodeId != 0 && contextEntries[nodeId] == NULL )
				return 0.0;
			unsigned *cnt;
			if ( nodeId == 0 )
				cnt = _counts->find(val, 0, _probStartBlockSize);
			else
				cnt = _counts->find(val, contextEntries[nodeId]->probTableOffSet, contextEntries[nodeId]->probTableBlockSize);
			if ( ! cnt )
				return 0.0;

			return (double)(*cnt) / ((double)_counts->size(contextEntries[nodeId]->probTableOffSet, contextEntries[nodeId]->probTableBlockSize));
		}
		break;
	case CountsProdCardinalityNorm:
	case CountsSumCardinalityNorm:
	case CountsSumLogCardinalityNorm: {
			if ( nodeId != 0 && contextEntries[nodeId] == NULL )
				return 0.0;
			unsigned *cnt;
			if ( nodeId == 0 )
				cnt = _counts->find(val, 0, _probStartBlockSize);
			else
				cnt = _counts->find(val, contextEntries[nodeId]->probTableOffSet, contextEntries[nodeId]->probTableBlockSize);
			if ( ! cnt )
				return 0.0;

			double numerator = *cnt;
			double denominator;
			if ( parentsStrategy == CountsProdCardinalityNorm )
				denominator = bgNode.prodCardinalities;
			else if ( parentsStrategy == CountsSumLogCardinalityNorm )
				denominator = bgNode.sumCardinalities;
			else
				denominator = bgNode.sumLogCardinalities;
			return numerator / denominator;
		}
		break;
	case BogNodeProb:
		// chose the one with the maximum smoothed probability
		// this is also a max MI criterion.
		return probBackingOff(val, nodeId, contextEntries, parentsValues).unlog();
	default:
		error("Error: unknown backoff strategy, value = %d\n", parentsStrategy);
	}

	return -1.0;
}


/*************************************
 * methods for FNGramCPT::ParentType *
 *************************************/


/*-
 *-----------------------------------------------------------------------
 * FFNGramImp::ParentType::read
 *      Read a parent specification such as W(-1).
 *
 * Results:
 *      None.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
void FNGramImp::ParentType::read(iDataStreamFile &ifs, shash_map<char, unsigned> &tagMap) {
	char c;
	ifs.readChar(c);
	unsigned* it = tagMap.find(c);
	if ( it == NULL )
		error("Error: tag %c not found in %s", c, ifs.fileName());
	index = *it;

	ifs.readChar(c);
	if ( c != '(' )
		error("Error: reading flm. expecting ( in %s", ifs.fileName());

	ifs.readInt(offset);

	ifs.readChar(c);
	if ( c != ')' )
		error("Error: reading flm. expecting ( in %s", ifs.fileName());
}


void FNGramImp::ParentType::parseFromString(char* str, shash_map<char, unsigned> &tagMap) {
	if ( str == NULL )
		error("Error: string is null in FNGramCPT::ParentType::parseFromString");

	unsigned* it = tagMap.find(str[0]);
	if ( it == NULL )
		error("Error: tag %c not found in %s", str[0], str);
	index = *it;

	if ( str[1] != '(' )
		error("Error: expecting ( in %s", str);
	char* end = strchr(str, ')');
	if ( end == NULL )
		error("Error: expecting ( in %s", str);
	*end = '\0';

	offset = atoi(str + 2);
}

/*-
 *-----------------------------------------------------------------------
 * FNGramImp::ParentType::operator ==
 *      Compare two parents.
 *
 * Results:
 *      Return true iff the two are identical.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
bool FNGramImp::ParentType::operator == (const ParentType &parent) const {
	return index == parent.index && offset == parent.offset;
}


/*******************************************
 * methods for FNGramImp::BackoffGraphNode *
 *******************************************/


/*-
 *-----------------------------------------------------------------------
 * FNGramImp::BackoffGraphNode::BackoffGraphNode
 *      Default constructor.
 *
 * Results:
 *      None.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
FNGramImp::BackoffGraphNode::BackoffGraphNode() : gtFile(NULL),  knFile(NULL), countFile(NULL), wmean(NULL), contextTable(NULL) {
	valid = false;
	order = 0;
	numBGChildren = 0;
	backoffConstraint = ~0x0u;
	backoffStrategy = CountsSumCountsNorm;
	backoffCombine = MaxBgChild;
	gtmin = 1;
	gtmax = 7;
	gtFile = 0;
	cdiscount = -1.0;
	ndiscount = false;
	wbdiscount = false;
	kndiscount = false;
	ukndiscount = false;
	knCountsModified = false;
	knCountsModifyAtEnd = false;
	knCountParent = ~0x0u;
	interpolate = false;
}


/*-
 *-----------------------------------------------------------------------
 * FNGramImp::BackoffGraphNode::~BackoffGraphNode
 *      Default destructor.
 *
 * Results:
 *      None.
 *
 * Side Effects:
 *      Clean up memory.
 *
 *-----------------------------------------------------------------------
 */
FNGramImp::BackoffGraphNode::~BackoffGraphNode() {
	delete [] gtFile;
	delete [] knFile;
	delete [] countFile;
	delete [] wmean;
	delete contextTable;
}


inline bool FNGramImp::BackoffGraphNode::requiresGenBackoff() {
	return ! ((numBGChildren <= 1) || (backoffCombine == AvgBgChild) || (backoffCombine == WmeanBgChild));
}


/************************************
 * methods for FNGramImp::LevelIter *
 ************************************/

/*-
 *-----------------------------------------------------------------------
 * FNGramImp::LevelIter::next
 *      Default constructor.
 *
 * Results:
 *      None.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
FNGramImp::LevelIter::LevelIter(const unsigned int _numParents,const unsigned int _level)
: numParents(_numParents), numNodes(1 << _numParents), level(_level) {
	init();
}


/*-
 *-----------------------------------------------------------------------
 * FNGramImp::LevelIter::init
 *      Initialize the iterator.
 *
 * Results:
 *      None.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
void FNGramImp::LevelIter::init() {
	state = 0;
}


/*-
 *-----------------------------------------------------------------------
 * FNGramImp::LevelIter::next
 *      Advance the iterator.
 *
 * Results:
 *      Return true if next exits.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
bool FNGramImp::LevelIter::next(unsigned int&node) {
	for ( ; state < numNodes; state++ ) {
		if ( numBitsSet(state) == level ) {
			node = state++;
			return true;
		}
	}
	return false;
}


/**************************************
 * methods for FNGramImp::BGChildIter *
 **************************************/


/*-
 *-----------------------------------------------------------------------
 * FNGramImp::BGChildIter::BGChildIter
 *      Default constructor.
 *
 * Results:
 *      None.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
FNGramImp::BGChildIter::BGChildIter(const unsigned int _numParents, const unsigned int _homeNode, const unsigned _constraint)
: numParents(_numParents), homeNode(_homeNode), constraint(_constraint & homeNode) {	// remove extra 1's
	init();
}


/*-
 *-----------------------------------------------------------------------
 * FNGramImp::BGChildIter::init
 *      Initialize the iterator.
 *
 * Results:
 *      None.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
void FNGramImp::BGChildIter::init() {
	mask = 1 << (numParents - 1);
}


/*-
 *-----------------------------------------------------------------------
 * FNGramImp::BGChildIter::next
 *      Advance an iterator.
 *
 * Results:
 *      Return true if next exists.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
bool FNGramImp::BGChildIter::next(unsigned int&node) {
	while ( mask != 0 ) {
		if ( (mask & constraint) != 0 ) {
			node = (~mask) & homeNode;
			mask >>= 1;
			return true;
		}
		mask >>= 1;
	}

	return false;
}


/*-
 *-----------------------------------------------------------------------
 * FNGramImp::BGGrandChildIter::BGGrandChildIter
 *      Default constructor.
 *
 * Results:
 *      None.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
FNGramImp::BGGrandChildIter::BGGrandChildIter(const unsigned _numParents, const unsigned _homeNode, BackoffGraphNode *_fngramNodes, const unsigned _great)
: numParents(_numParents), numNodes(1<<_numParents), homeNode(_homeNode), numBitsSetOfHomeNode(numBitsSet(homeNode)), fngramNodes(_fngramNodes), great(_great) {
	init();
}


/*-
 *-----------------------------------------------------------------------
 * FNGramImp::BGGrandChildIter::next
 *      Advance an iterator.
 *
 * Results:
 *      Return true if next exists.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
bool FNGramImp::BGGrandChildIter::next(unsigned int&node) {
	// TODO make this faster
	for ( ; state >= 0; state-- ) {
		// all bits in child=state must also be on in homeNode
		if ( ((state & homeNode) == (unsigned)state) && ((great + 2 + numBitsSet(state)) == numBitsSetOfHomeNode) && fngramNodes[state].valid ) {
			node = state--;
			return true;
		}
	}

	return false;
}


const unsigned FNGramImp::maxNumParentsPerChild = 32;


/*************************
 * methods for FNGramCPT *
 *************************/

/*-
 *-----------------------------------------------------------------------
 * FNGramCPT::FNGramCPT
 *      Default constructor for factored ngram CPT
 *
 * Results:
 *      None.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
FNGramCPT::FNGramCPT() : CPT(di_FNGramCPT), _fngram(NULL), _startNode(0), _numberOfActiveIterators(0) {
	_numParents = 0;
}


/*-
 *-----------------------------------------------------------------------
 * FNGramCPT::~FNGramCPT
 *      Default destructor for factored ngram CPT
 *
 * Results:
 *      None.
 *
 * Side Effects:
 *      Clean up the memory.
 *
 *-----------------------------------------------------------------------
 */
FNGramCPT::~FNGramCPT() {
	for ( unsigned i = 0; i < _contextEntriesStack.size(); i++ ) {
		free(_contextEntriesStack[i]);
	}
}


void FNGramCPT::setNumParents(unsigned nParents) {
	_numParents = nParents;

	// reset size
	cardinalities.resize(_numParents);
	_parentsPositions.resize(_numParents);
}


void FNGramCPT::setFNGramImp(FNGramImp *fngram) {
	if ( _numParents > fngram->_numParents )
		error("number of parents %d is bigger than internal fngram %d", _numParents, fngram->_numParents);

	_fngram = fngram;
	_card = _fngram->_card;
	if ( _numParents == _fngram->_numParents ) {
		for ( unsigned i = 0; i < _numParents; i++ ) {
			_parentsPositions[i] = i;
			cardinalities[i] = _fngram->cardinalities[i];
		}
	}

	_startNode = (1 << _fngram->_numParents) - 1;

#ifdef FNGRAM_BOW_GROW_DYNA
	_parentsValues.resize(_fngram->_numParents);
#endif

	_contextEntriesStack.resize(1);
	_contextEntriesStack[0] = (void *) malloc(sizeof(FNGramImp::BackoffGraphNode::HashEntry*) * 4 * _fngram->_numberOfBGNodes);
	_numberOfActiveIterators = 0;
}


void FNGramCPT::setParentsPositions(const vector<unsigned> &parentsPositions) {
	if ( parentsPositions.size() != _numParents )
		error("number of positions for parents %d is different from number of parents %d in FNGramCPT::setParentsPositions", parentsPositions.size(), _numParents);

	_parentsPositions = parentsPositions;

	// calculate the starting node
	_startNode = 0x0U;

	for ( unsigned i = 0; i < _numParents; i++ ) {
		if ( _parentsPositions[i] > _numParents )
			error("parent %d position %d is bigger or equal to number of parents %d in FNGramCPT::setParentsPositions", i, _parentsPositions[i], _numParents);
		_startNode |= (1U << _parentsPositions[i]);
		cardinalities[i] = _fngram->cardinalities[_parentsPositions[i]];
	}

	if ( ! _fngram->_bgNodes[_startNode].valid )
		error("the node in FNGramCPT %x is not valid in backing-off path of FLM", _startNode);
}


/*-
 *-----------------------------------------------------------------------
 * FNGramCPT::becomeAwareOfParentValues
 *      Set context entry pointers when parents values are known.
 *
 * Results:
 *      None.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
void FNGramCPT::becomeAwareOfParentValues(vector< RV* >& parents, const RV* rv) {
	error("FNGramCPT::becomeAwareOfParentValues shouldn't be used alone");
}


/*-
 *-----------------------------------------------------------------------
 * FNGramCPT::begin
 *      Set an iterator to the begin (and sets the probability directly)
 *
 * Results:
 *      None.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
void FNGramCPT::begin(CPT::iterator& it, DiscRV* drv, logpr& p) {
	error("FNGramCPT::begin shouldn't be used alone");
}


/*-
 *-----------------------------------------------------------------------
 * FNGramCPT::becomeAwareOfParentValuesAndIterBegin
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
void FNGramCPT::becomeAwareOfParentValuesAndIterBegin(vector< RV*>& parents,
						      iterator &it,
						      DiscRV*drv,
						      logpr& p) {
	assert(parents.size() == _numParents);

	FNGramImp::BackoffGraphNode::HashEntry **ptr;

	if ( _numberOfActiveIterators >= _contextEntriesStack.size() * 4 ) {
		void *tmpPtr = (void *) malloc(sizeof(FNGramImp::BackoffGraphNode::HashEntry*) * 4 * _fngram->_numberOfBGNodes);
		_contextEntriesStack.push_back(tmpPtr);
		ptr = (FNGramImp::BackoffGraphNode::HashEntry**)tmpPtr;
	} else {
		ptr = ((FNGramImp::BackoffGraphNode::HashEntry**)_contextEntriesStack[_numberOfActiveIterators / 4]) + (_numberOfActiveIterators % 4) * _fngram->_numberOfBGNodes;
	}
	// increment counter and set inter point for iterator
	_numberOfActiveIterators++;

	memset(ptr, 0, sizeof(FNGramImp::BackoffGraphNode::HashEntry*) * _fngram->_numberOfBGNodes);

	unsigned * context = new unsigned [_numParents];

	// descend down the BG, level by level except buttom (0)
	for ( int level = _numParents; level > 0; level-- ) {
		FNGramImp::LevelIter liter(_numParents, level);
		for ( unsigned nodeAtLevel; liter.next(nodeAtLevel); ) {
			// we make sure it is somewhere below startNode
			if ( _fngram->_bgNodes[nodeAtLevel].valid && (nodeAtLevel | _startNode == _startNode) ) {
				// figure out which bits are on
				std::vector<unsigned> onPos = FNGramImp::bitsOn(nodeAtLevel);	// this will be something like 0, 1, 3
				for ( unsigned i = 0; i < _fngram->_bgNodes[nodeAtLevel].order - 1; i++ ) {
					context[i] = RV2DRV(parents[onPos[i]])->val;
				}
				ptr[nodeAtLevel] = _fngram->_bgNodes[nodeAtLevel].contextTable->find(context);
			}
		}
	}

	delete [] context;

#ifdef FNGRAM_BOW_GROW_DYNA
	for ( unsigned i = 0; i < _numParents; i++ )
		_parentsValues[_parentsPositions[i]] = RV2DRV(parents[i])->val;
#endif

	it.drv = drv;
	it.internalStatePtr = ptr;
	register DiscRVType value = 0;
	p = _fngram->probBackingOff(value, _startNode, ptr, _parentsValues);

	while ( p.essentially_zero() ) {
		value++;
		// We keep the following assertion as we
		// must have that at least one entry is non-zero.
		// The read code of the FNGramCPT should ensure this
		// as sure all parameter update procedures.
		assert(value < card());
		p = _fngram->probBackingOff(value, _startNode, ptr, _parentsValues);
	}
	drv->val = value;
}


/*-
 *-----------------------------------------------------------------------
 * FNGramCPT::probGivenParents
 *      Retrieve the probability for a given value.
 *
 * Results:
 *      Return the probability given the parents.
 *
 * Side Effects:
 *
 *-----------------------------------------------------------------------
 */
logpr FNGramCPT::probGivenParents(vector < RV* >& parents, DiscRV* drv) {
	assert(parents.size() == _numParents);

	FNGramImp::BackoffGraphNode::HashEntry ** contextEntries = new FNGramImp::BackoffGraphNode::HashEntry * [_fngram->_numberOfBGNodes];
	memset(contextEntries, 0, sizeof(FNGramImp::BackoffGraphNode::HashEntry*) * _fngram->_numberOfBGNodes);

	unsigned * context = new unsigned [_numParents];

	// descend down the BG, level by level except buttom (0)
	for ( int level = _numParents; level > 0; level-- ) {
		FNGramImp::LevelIter liter(_numParents, level);
		for ( unsigned nodeAtLevel; liter.next(nodeAtLevel); ) {
			// we make sure it is somewhere below startNode
			if ( _fngram->_bgNodes[nodeAtLevel].valid && (nodeAtLevel | _startNode == _startNode) ) {
				// figure out which bits are on
				std::vector<unsigned> onPos = FNGramImp::bitsOn(nodeAtLevel);	// this will be something like 0, 1, 3
				for ( unsigned i = 0; i < _fngram->_bgNodes[nodeAtLevel].order - 1; i++ ) {
					context[i] = RV2DRV(parents[onPos[i]])->val;
				}
				contextEntries[nodeAtLevel] = _fngram->_bgNodes[nodeAtLevel].contextTable->find(context);
			}
		}
	}

	delete [] context;

#ifdef FNGRAM_BOW_GROW_DYNA
	for ( unsigned i = 0; i < _numParents; i++ )
		_parentsValues[_parentsPositions[i]] = RV2DRV(parents[i])->val;
#endif

	logpr prob = _fngram->probBackingOff(drv->val, _startNode, contextEntries, _parentsValues);

	delete [] contextEntries;
	return prob;
}



/*-
 *-----------------------------------------------------------------------
 * FNGramCPT::next
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
bool FNGramCPT::next(iterator &it, logpr& p) {
	do{
		if ( (++it.drv->val) >= card() ) {
			// need to clean the stack pointer
			--_numberOfActiveIterators;

			return false;
		}
		p = _fngram->probBackingOff(it.drv->val, _startNode, (FNGramImp::BackoffGraphNode::HashEntry **)it.internalStatePtr, _parentsValues);
	} while ( p.essentially_zero() );

	return true;
}


/*-
 *-----------------------------------------------------------------------
 * FNGramCPT::randomSample
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
int FNGramCPT::randomSample(DiscRV* drv)
{
  logpr prob = rnd.drand48();
  iterator it;
  logpr p;
  FNGramCPT::begin(it,drv,p);
  logpr sum;
  do {
    sum += p;
    if ( prob <= sum )
      break;
  } while ( FNGramCPT::next(it,p) );
  return drv->val;
}


/*-
 *-----------------------------------------------------------------------
 * FNGramCPT::read
 *      Read a language model file form ARPA file or index file.
 *
 * Results:
 *      None.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
void FNGramCPT::read(iDataStreamFile &is) {
}
