/*
 * GMTK_FNGramCPT
 *      .h file for the GMTK_FNGramCPT.cc file.
 *
 *    This is the data type for factored langauge models in GMTK.
 *
 * Written by Gang Ji <gang@ee.washington.edu>
 * Borrowed some part from SRI language model toolkit with the factored
 * language model support written by Jeff Bilmes <bilmes@ee.washington.edu>
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


#ifndef GMTK_FNGRAM_CPT_H
#define GMTK_FNGRAM_CPT_H


#include <vector>

#include "GMTK_CPT.h"
#include "GMTK_Vocab.h"
#include "GMTK_RandomVariable.h"
#include "GMTK_DiscreteRandomVariable.h"
#include "fileParser.h"
#include "vshash_map.h"
#include "shash_map.h"
#include "hash_mtable.h"
#include "logp.h"


// with this macro defined, the code will grow the backing-off value
// hash table dynamically.
#define FNGRAM_BOW_GROW_DYNA


/*-
 * data type for factored language model CPT
 */
class FNGramCPT : public CPT {

public:
	// constructors and destructor
	FNGramCPT();
	~FNGramCPT();

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
	// return the probability of 'val' given the parents are the
	// assigned to the set of values set during the most previous call
	// to becomeAwareOfParentValues.
	virtual logpr probGivenParents(DiscreteRandomVariable* drv);
	// Similar to the above, but convenient for one time probability
	// evaluation.
	virtual logpr probGivenParents(vector < RandomVariable *>& parents, DiscreteRandomVariable* drv);

	// returns an iterator for the first one.
	virtual iterator begin(DiscreteRandomVariable* drv);

	virtual void begin(iterator& it, DiscreteRandomVariable* drv);
	virtual void begin(iterator& it, DiscreteRandomVariable* drv, logpr& p);
	virtual void becomeAwareOfParentValuesAndIterBegin(vector<RandomVariable *>& parents, iterator &it, DiscreteRandomVariable* drv);
	virtual void becomeAwareOfParentValuesAndIterBegin(vector<RandomVariable *>& parents, iterator &it, DiscreteRandomVariable* drv, logpr& p);

	// Given a current iterator, return true if it is a valid next
	// value, otherwise return false so a loop can terminate.
	virtual bool next(iterator &it);
	virtual bool next(iterator &it,logpr& p);

	// returns true if iterate is at end state
	virtual bool end(iterator &it) {return it.drv->val >= (int)ucard();}

	///////////////////////////////////////////////////////////
	// Given the current parent values, generate a random sample.
	virtual int randomSample(DiscreteRandomVariable*drv);

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

	void setTagMap(shash_map<char, unsigned> &tagMap) {_tagMap = tagMap;}
	void readFNGramSpec(iDataStreamFile &ifs);
	void readLMFile(std::vector<Vocab*> vocabs);
	void readCountFile(std::vector<Vocab*> vocabs);

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
	// data structures

	/** different backing-off node strategies */
	enum BackoffNodeStrategy {
		CountsNoNorm,				// use absolute counts
		CountsSumCountsNorm,			// norm by sum counts at level (== event MI maximization)
		CountsSumNumWordsNorm,			// norm by num words
		CountsProdCardinalityNorm,		// norma by prod. cardinality of gram
		CountsSumCardinalityNorm,		// norma by sum. cardinality of gram
		CountsSumLogCardinalityNorm,		// norma by sum log cardinality of gram
		BogNodeProb				// backoff graph node probability
	};

	/** combination methods for each node with more than one child */
	enum BackoffNodeCombine {
		MaxBgChild,
		MinBgChild,
		SumBgChild,
		ProdBgChild,
		AvgBgChild,
		GmeanBgChild,
		WmeanBgChild
	};

	/** parent for the factored language mdoel */
	struct ParentType {
		void read(iDataStreamFile &ifs, shash_map<char, unsigned> &tagMap);
		void parseFromString(char* str, shash_map<char, unsigned> &tagMap);
		bool operator == (const ParentType &parent) const;

		// attributes
		unsigned index;		// index for the parent, same as tag id
		int offset;		// offset can be both postive and negative
	};

	/** node in the backing-off graph */
	struct BackoffGraphNode {
		// hash entry
		struct HashEntry {
			logpr backingOffWeight;
			unsigned probTableOffSet;
			unsigned probTableBlockSize;

			HashEntry() : backingOffWeight(0.0), probTableOffSet(0), probTableBlockSize(0) {
			}
		};

		// methods
		BackoffGraphNode();
		~BackoffGraphNode();

		// this will be removed for release
		void print() const;

		inline bool requiresGenBackoff();

		// properties
		bool valid;				// whether this node is used
		unsigned order;
		unsigned numBGChildren;
		unsigned backoffConstraint;		// bit mask for children
		BackoffNodeStrategy backoffStrategy;
		BackoffNodeCombine backoffCombine;
		unsigned gtmin;
		unsigned gtmax;
		char *gtFile;
		double cdiscount;
		bool ndiscount;
		bool wbdiscount;
		bool kndiscount;
		bool ukndiscount;
		char *knFile;
		bool knCountsModified;
		bool knCountsModifyAtEnd;
		unsigned knCountParent;
		bool interpolate;
		char *countFile;			// counts file for this node
		logpr *wmean;

		// statistics
		double prodCardinalities;
		double sumCardinalities;
		double sumLogCardinalities;

		// language model data
		vshash_map<unsigned, HashEntry> *contextTable;
	};

	/**
	 * iterating each levels starting from the top (all parents exist).
	 */
	class LevelIter {
	public:
		LevelIter(const unsigned int _numParents,const unsigned int _level);

		void init();
		inline bool next(unsigned int&node);

	protected:
		const unsigned int numParents;
		const unsigned int numNodes;
		unsigned int state;
		const unsigned int level;
	};

	/**
	 * this is an iterator for all possible children for a node
	 * Note that this iterator doesn't consider the backing-off constaint.
	 * Therefore, it iterates ALL childen.
	 */
	class BGChildIter {
	public:
		BGChildIter(const unsigned _numParents, const unsigned _homeNode, const unsigned constraint);
		void init();
		inline bool next(unsigned &node);

	protected:
		unsigned mask;
		const unsigned numParents;
		const unsigned homeNode;
		const unsigned constraint;
	};

	/**
	 * iterator for grand children
	 */
	class BGGrandChildIter {
	public:
		BGGrandChildIter(const unsigned _numParents, const unsigned _homeNode, BackoffGraphNode *_fngramNodes, const unsigned _great=0);
		void init() { state = ((int)homeNode - ((1 << (great + 1)) - 1)); }
		bool next(unsigned int&node);

	protected:
		const unsigned numParents;
		const unsigned numNodes;
		int state;

		const unsigned homeNode;
		const unsigned numBitsSetOfHomeNode;
		BackoffGraphNode const* fngramNodes;
		const unsigned great;
	};

	// misc helping methods
	static long strtolplusb(const char *nptr, char **endptr, int base);
	static inline unsigned numBitsSet(unsigned u);
	static const std::vector<unsigned> bitsOn(unsigned u);

	// returns the type of the sub-object in string
	// form that is suitable for printing and identifying
	// the type of the object.
	virtual const string typeName() {return std::string("FNGramCPT");}

	// methods for parsing structure file
	unsigned parseNodeString(char *str);
	void readBackoffGraphNode(iDataStreamFile &ifs, BackoffGraphNode &bgNode, unsigned nodeId);

	// statistics
	void computeCardinalityFunctions();

	// following procedures are for calculating probabilities when parents are known
	logpr probBackingOff(unsigned val, unsigned nodeId);
	logpr bgChildProbBackingOff(unsigned val, unsigned nodeId);
	logpr bgChildProbSum(unsigned nodeId);
	unsigned boNode(unsigned val, unsigned nodeId);
	double backoffValueRSubCtxW(unsigned val, BackoffNodeStrategy parentsStrategy, unsigned nodeId);

	// data fields
	unsigned _childIndex;				// child id as W in P(W|F1, F2,...)
	ParentType *_parents;				// parent specifications
	shash_map<char, unsigned> _tagMap;	// parents tag map

	unsigned _numberOfBGNodes;			// total number of possible backing-off graph nodes
	BackoffGraphNode *_bgNodes;			// backing-off graph nodes

	char *_countFileName;				// count file. might be useful for some backing-off strategies
	char *_lmFileName;					// language model file in sudo ARPA format

	HashMTable<logpr> *_probabilities;
	HashMTable<unsigned> *_counts;
	unsigned _probStartBlockSize;

	// total number of probabilities and backing-off weights
	unsigned _totalNumberOfParameters;

	// this will be set when parents value is known
	BackoffGraphNode::HashEntry ** _contextEntries;

#ifdef FNGRAM_BOW_GROW_DYNA
	std::vector<unsigned> _parentsValues;
#endif

	static const unsigned maxNumParentsPerChild;
};


#endif
