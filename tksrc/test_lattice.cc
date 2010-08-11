/**
 * @file test_lattice.cc
 * @author Gang Ji
 * @author http://welcome.to/rainier
 *
 * $Id$
 */


#include "rand.h"
#include "fileParser.h"
#include "GMTK_FileParser.h"
#include "GMTK_ObservationMatrix.h"
#include "GMTK_GMParms.h"
#include "GMTK_Vocab.h"
#include "GMTK_LatticeADT.h"
#include "GMTK_LatticeEdgeCPT.h"
#include "GMTK_LatticeNodeCPT.h"
#include "GMTK_HidDiscRV.h"


// global variables that must be defined
RAND rnd;
GMParms GM_Parms;
ObservationMatrix globalObservationMatrix;
 

int main() {
	// read in master file
	{
		iDataStreamFile pf("testLattice.master", false, true);
		GM_Parms.read(pf);
	}
	GM_Parms.finalizeParameters();

	// read in structure file
	FileParser fp("testLattice.str");

	// parse the file
	fp.parseGraphicalModel();
	// create the rv variable objects
	fp.createRandomVariableGraph();
	// Make sure that there are no directed loops in the graph.
	fp.ensureValidTemplate();
	fp.associateWithDataParams(FileParser::noAllocate);

	// find the lattice node cpt
	if ( GM_Parms.latticeNodeCptsMap.find(string("sampleLattice")) == GM_Parms.latticeNodeCptsMap.end() ) {
		printf("cannot find fngram in LatticeNodeCPT\n");
		return -1;
	}
	LatticeNodeCPT *nodeCpt = GM_Parms.latticeNodeCpts[GM_Parms.latticeNodeCptsMap[string("sample")]];

	// define misc RVs
	RVInfo nodeInfo;
	HidDiscRV prntNode(nodeInfo, 1, 24), node(nodeInfo, 0, 24);
	logpr prob;
	LatticeNodeCPT::iterator it;
	vector< RV* > parents(1);
	parents[0] = &prntNode;

	prntNode.val = 0;
	nodeCpt->becomeAwareOfParentValuesAndIterBegin(parents, it, &node, prob);
	do {
		printf("from %d to %d with %f\n", prntNode.val, node.val, prob.val());
	} while ( (! prob.zero()) && nodeCpt->next(it, prob) );

	prntNode.val = 1;
	nodeCpt->becomeAwareOfParentValuesAndIterBegin(parents, it, &node, prob);
	do {
		printf("from %d to %d with %f\n", prntNode.val, node.val, prob.val());
	} while ( (! prob.zero()) && nodeCpt->next(it, prob) );

	// find the lattice edge cpt
	if ( GM_Parms.latticeEdgeCptsMap.find(string("sampleLattice")) == GM_Parms.latticeEdgeCptsMap.end() ) {
		printf("cannot find fngram in LatticeNodeCPT\n");
		return -1;
	}
	LatticeEdgeCPT *edgeCpt = GM_Parms.latticeEdgeCpts[GM_Parms.latticeEdgeCptsMap[string("sample")]];

	RVInfo wordInfo;
	HidDiscRV wordNode(wordInfo, 0, 11);
	parents.resize(2);
	parents[0] = &prntNode;
	parents[1] = &node;

	prntNode.val = 0;
	node.val = 1;
	LatticeEdgeCPT::iterator eit;
	edgeCpt->becomeAwareOfParentValuesAndIterBegin(parents, eit, &wordNode, prob);
	do {
		printf("from %d to %d emitting %d with %f\n", prntNode.val, node.val, wordNode.val, prob.val());
	} while ( (! prob.zero()) && edgeCpt->next(eit, prob) );

	return 0;
}
