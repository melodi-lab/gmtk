/*-
 * GMTK_LatticeADT.h
 *      .h file for GMTK_LatticeADT.cc, HTK lattice support
 *      distributions.
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


#ifndef GMTK_LATTICE_ADT_H
#define GMTK_LATTICE_ADT_H


#include "GMTK_NamedObject.h"
#include "GMTK_Vocab.h"
#include "fileParser.h"
#include "shash_map_iter.h"
#include "logp.h"


/**
 * HTK lattice support
 */
class LatticeADT : public NamedObject {
public:
	LatticeADT();
	~LatticeADT();

	// read from HTK lattice format file
	void readFromHTKLattice(iDataStreamFile &ifs, const Vocab &vocab);

	// read from GMTK master file
	void read(iDataStreamFile &is);

	// this lattice is iterable
	inline bool iterable() const { return _latticeFile != NULL; }

	void seek(unsigned nmbr);
	void initializeIterableLattice(const string &fileName);
	void beginIterableLattice();
	void nextIterableLattice();

	// reset the frame indices
	void resetFrameIndices(unsigned numFrames);

	friend class LatticeNodeCPT;
	friend class LatticeEdgeCPT;

protected:
	/**
	 * lattice edge
	 */
	struct LatticeEdge {
		/** emission id from current state to next state */
		unsigned emissionId;
		/** acoustic score */
		logpr ac_score;
		/** acoustic score */
		logpr lm_score;
		/** probability score */
		logpr prob_score;

		LatticeEdge() : emissionId(0) {}
	};

	/**
	 * lattice node information
	 */
	struct LatticeNode {
		/** absolute time for this node */
		float time;
		/** starting frame number */
		unsigned startFrame;
		/** ending frame number */
		unsigned endFrame;
		/** possible out-going edges */
		shash_map_iter<unsigned, LatticeEdge> edges;

		LatticeNode() : startFrame(0), endFrame(0) {}
	};

	/** lattice nodes */
	LatticeNode *_latticeNodes;
	/** number of nodes in lattice this should be smaller than node cardinality */
	unsigned _numberOfNodes;
	/** number of links in lattice */
	unsigned _numberOfLinks;
	/** start node id */
	unsigned _start;
	/** end node id */
	unsigned _end;
	/** language model scale */
	double _lmscale;
	/** word penalty */
	double _wdpenalty;
	/** acoustic model scale */
	double _acscale;
	/** ? */
	double _amscale;
	/** log base */
	double _base;

	/** how many frame relaxation is allowed in CPT */
	unsigned _frameRelax;

	/** if this is an iterable cpt */
	iDataStreamFile* _latticeFile; // the file pointer
	string _latticeFileName; // the file name
	unsigned _numLattices; // number of lattices in this file
	int _curNum; // the current lattice number
	string _curName; // the current lattice cpt name

	// the following is GM paramters
	unsigned _nodeCardinality;
	unsigned _wordCardinality;
};


#endif
