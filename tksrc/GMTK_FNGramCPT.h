/*
 * GMTK_FNGramCPT
 *      .h file for the GMTK_FNGramCPT.cc file.
 *
 * The data type for factored langauge models in GMTK.
 *
 * Written by Gang Ji <gang@ee.washington.edu>
 * Modifications by J. Bilmes. <bilmes@ee.washington.edu>
 *
 * Portions of this code were borrowed from the SRI language model
 * toolkit, but only those parts that were written by Jeff Bilmes
 * <bilmes@ee.washington.edu> as part of the JHU CLSP 2002 workshop.  
 * NO OTHER PORTIONS OF SRI CODE ARE CONTAINED HERE!! THEREFORE, THISb
 * CODE IS NOT UNDER SRI COPYRIGHT.
 *
 * This part of the code has the implementation of class FNGramCPT for
 * gmtk.  Please see GMTK_FNGramCPT.h for more information.
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
#include "GMTK_DiscRV.h"
#include "fileParser.h"
#include "vshash_map.h"
#include "shash_map.h"
#include "hash_mtable.h"
#include "logp.h"


// with this macro defined, the code will grow the backing-off value
// hash table dynamically.
#define FNGRAM_BOW_GROW_DYNA


/**
 * This is what's inside FNGramCPT
 * This is seperated from FNGramCPT so that it can support fewer parents and the FLM itself.
 */
class FNGramImp : public NamedObject {
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

public:
	FNGramImp();
	~FNGramImp();

	void setNumParents(const unsigned nParents);

	void setTagMap(shash_map<char, unsigned> &tagMap) {_tagMap = tagMap;}
	void read(iDataStreamFile &ifs);
	void readFNGramSpec(iDataStreamFile &ifs);
	void readLMFile(std::vector<Vocab*> vocabs);
	void readCountFile(std::vector<Vocab*> vocabs);

	logpr probBackingOff(unsigned val, unsigned nodeId, BackoffGraphNode::HashEntry **contextEntries, const std::vector<unsigned> &parentsValues);

	// return the number of parameters for object.
	unsigned totalNumberOfParameters() {return _totalNumberOfParameters;}

protected:
	friend class FNGramCPT;

	// misc helping methods
	static long strtolplusb(const char *nptr, char **endptr, int base);
	static inline unsigned numBitsSet(unsigned u);
	static const std::vector<unsigned> bitsOn(unsigned u);

	// methods for parsing structure file
	unsigned parseNodeString(char *str);
	void readBackoffGraphNode(iDataStreamFile &ifs, BackoffGraphNode &bgNode, unsigned nodeId);

	// statistics
	void computeCardinalityFunctions();

	// following procedures are for calculating probabilities when parents are known
	logpr bgChildProbBackingOff(unsigned val, unsigned nodeId, BackoffGraphNode::HashEntry **contextEntries, const std::vector<unsigned> &parentsValues);
	logpr bgChildProbSum(unsigned nodeId, BackoffGraphNode::HashEntry **contextEntries, const std::vector<unsigned> &parentsValues);
	unsigned boNode(unsigned val, unsigned nodeId, BackoffGraphNode::HashEntry **contextEntries, const std::vector<unsigned> &parentsValues);
	double backoffValueRSubCtxW(unsigned val, BackoffNodeStrategy parentsStrategy, unsigned nodeId,
		BackoffGraphNode::HashEntry **contextEntries, const std::vector<unsigned> &parentsValues);

	// data fields
	unsigned _childIndex;				// child id as W in P(W|F1, F2,...)
	unsigned _card;
	unsigned _numParents;				// number of parents
	ParentType *_parents;				// parent specifications
	std::vector<unsigned> cardinalities;
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

	static const unsigned maxNumParentsPerChild;
};


/*-
 * data type for factored language model CPT
 */
class FNGramCPT : public CPT {

public:
	// constructors and destructor
	FNGramCPT();
	~FNGramCPT();

	void setFNGramImp(FNGramImp *fngram);
	void setParentsPositions(const vector<unsigned> &parentsPositions);

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

	////////////////////////////////////////////////////////////////////////////
	// from base class EMable
	void emStartIteration() {}
	void emIncrement(logpr prob,vector < RV* >& parents, RV*r) {}
	void emEndIteration() {}
	void emSwapCurAndNew() {}

	// return the number of parameters for object.
	virtual unsigned totalNumberParameters() {return _fngram->totalNumberOfParameters();}

	///////////////////////////////////////////////////////////////
	// virtual functions for objects to do the actual work.
	virtual void emStoreObjectsAccumulators(oDataStreamFile& ofile) {}
	virtual void emLoadObjectsDummyAccumulators(iDataStreamFile& ifile) {}
	virtual void emZeroOutObjectsAccumulators() {}
	virtual void emLoadObjectsAccumulators(iDataStreamFile& ifile) {}
	virtual void emAccumulateObjectsAccumulators(iDataStreamFile& ifile) {}

protected:
	// returns the type of the sub-object in string
	// form that is suitable for printing and identifying
	// the type of the object.
	virtual const string typeName() {return std::string("FNGramCPT");}

	FNGramImp *_fngram;

	unsigned _startNode;
	std::vector<unsigned> _parentsPositions;

	// this will be set when parents value is known
#ifdef FNGRAM_BOW_GROW_DYNA
	std::vector<unsigned> _parentsValues;
#endif

	FNGramImp::BackoffGraphNode::HashEntry ** _contextEntries;
};


#endif
