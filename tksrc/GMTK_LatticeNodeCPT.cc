/*-
 * GMTK_LatticeNodeCPT.cc
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


#include "GMTK_LatticeNodeCPT.h"
#include "GMTK_CPT.h"
#include "GMTK_DiscRV.h"

LatticeNodeCPT::LatticeNodeCPT() : CPT(di_LatticeNodeCPT) , _latticeAdt(NULL) {
	// some values are fixed
	_numParents = 1;
	cardinalities.resize(1);
}


LatticeNodeCPT::~LatticeNodeCPT() {
}


void LatticeNodeCPT::becomeAwareOfParentValuesAndIterBegin(vector< RV* >& parents, iterator &it, DiscRV* drv, logpr& p) {
	// check out the out-going edges based on parent value
	shash_map2<unsigned, LatticeADT::LatticeEdge>::iterator *pit = new shash_map2<unsigned, LatticeADT::LatticeEdge>::iterator();
	assert(_latticeAdt->_latticeNodes[RV2DRV(parents[0])->val].edges.totalNumberEntries() != 0);
	*pit = _latticeAdt->_latticeNodes[RV2DRV(parents[0])->val].edges.begin();

	// set up the iternal state for iterators
	it.internalStatePtr = (void*)pit;
	it.drv = drv;

	// set up the values
	drv->val = pit->key();
	p = (**pit).ac_score;		// TODO: instead of acoustic score, implement other scores
}


logpr LatticeNodeCPT::probGivenParents(vector< RV* >& parents, DiscRV* drv) {
	LatticeADT::LatticeEdge* outEdge = _latticeAdt->_latticeNodes[RV2DRV(parents[0])->val].edges.find(drv->val);


	if ( outEdge == NULL )
		return logpr(0.0);
	else
		return outEdge->ac_score;	// TODO: implement other scores
}


bool LatticeNodeCPT::next(iterator &it, logpr& p) {
	shash_map2<unsigned, LatticeADT::LatticeEdge>::iterator* pit = (shash_map2<unsigned, LatticeADT::LatticeEdge>::iterator*) it.internalStatePtr;
	if ( pit->next() ) {
		it.drv->val = pit->key();
		p = (**pit).ac_score;		// TODO: instead of acoustic score, implement other scores
		return true;
	} else {
		delete pit;
		return false;
	}
}


void LatticeNodeCPT::setLatticeADT(const LatticeADT &latticeAdt) {
	_latticeAdt = &latticeAdt;
	_card = cardinalities[0] = _latticeAdt->_nodeCardinality;
}


#ifdef MAIN


#include "rand.h"
#include "GMTK_ObservationMatrix.h"
#include "GMTK_GMParms.h"
#include "GMTK_HidDiscRV.h"
#include "GMTK_LatticeEdgeCPT.h"
#include "GMTK_FileParser.h"
#include "fileParser.h"


RAND rnd;
GMParms GM_Parms;
ObservationMatrix globalObservationMatrix;


int main(int argc, char**argv) {
	iDataStreamFile pf("testLatticeCPT.master", false, true);
	GM_Parms.read(pf);
	GM_Parms.finalizeParameters();

	RVInfo nodeInfo, wordInfo;
	HidDiscRV node(nodeInfo, 0, 50), nodeprnt(nodeInfo, 1, 50), word(wordInfo, 0, 13);
	vector<RV*> parents(1);
	parents[0] = &nodeprnt;

	logpr prob;

	LatticeNodeCPT *nodeCpt = GM_Parms.latticeNodeCpts[GM_Parms.latticeNodeCptsMap[string("sampleLattice")]];

	nodeprnt.val = 0;
	node.val = 1;
	prob = nodeCpt->probGivenParents(parents, &node);
	printf("0->1: %f\n", prob.val());

	nodeprnt.val = 15;
	node.val = 18;
	prob = nodeCpt->probGivenParents(parents, &node);
	printf("15->18: %f\n", prob.val());

	LatticeNodeCPT::iterator it;
	nodeCpt->becomeAwareOfParentValuesAndIterBegin(parents, it, &node, prob);

	do {
		printf("%d->%d: %f\n", nodeprnt.val, node.val, prob.val());
	} while ( nodeCpt->next(it, prob) );


	vector<RV*> edgeParents(2);
	edgeParents[0] = &nodeprnt;
	edgeParents[1] = &node;

	LatticeEdgeCPT *edgeCpt = GM_Parms.latticeEdgeCpts[GM_Parms.latticeEdgeCptsMap[string("sampleLattice")]];

	nodeprnt.val = 0;
	node.val = 1;
	word.val = 11;
	prob = edgeCpt->probGivenParents(edgeParents, &word);
	printf("0->1:11 %f\n", prob.val());

	word.val = 10;
	prob = edgeCpt->probGivenParents(edgeParents, &word);
	printf("0->1:10 %f\n", prob.val());

	nodeprnt.val = 15;
	node.val = 18;
	word.val = 6;
	prob = edgeCpt->probGivenParents(edgeParents, &word);
	printf("15->18:6 %f\n", prob.val());

	LatticeEdgeCPT::iterator eit;
	edgeCpt->becomeAwareOfParentValuesAndIterBegin(edgeParents, eit, &word, prob);
	do {
		printf("%d->%d:%d %f\n", nodeprnt.val, node.val, word.val, prob.val());
	} while ( edgeCpt->next(eit, prob) );

	FileParser fp("testLattice.str");
       	fp.parseGraphicalModel();
        // create the rv variable objects
        fp.createRandomVariableGraph();
        // Make sure that there are no directed loops in the graph.
        fp.ensureValidTemplate();
        fp.associateWithDataParams(FileParser::noAllocate);

}


#endif

