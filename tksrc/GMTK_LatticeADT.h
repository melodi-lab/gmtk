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
#include "shash_map2.h"
#include "logp.h"


/**
 * HTK lattice support
 */
class LatticeADT : public NamedObject {
public:
	LatticeADT();
	~LatticeADT();

	void readFromHTKLattice(iDataStreamFile &ifs, const Vocab &vocab);

	void read(iDataStreamFile &is);

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

		LatticeEdge() : emissionId(0) {}
	};

	/**
	 * lattice node information
	 */
	struct LatticeNode {
		/** starting frame number */
		unsigned startFrame;
		/** ending frame number */
		unsigned endFrame;
		/** possible out-going edges */
		shash_map2<unsigned, LatticeEdge> edges;

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
	/** default frame rate this is not yet used */
	double _frameRate;

	// the following is GM paramters
	unsigned _nodeCardinality;
	unsigned _wordCardinality;
};


#endif
