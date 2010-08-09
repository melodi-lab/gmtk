/**
 *: test_fngram.cc
 */

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "GMTK_Vocab.h"
#include "GMTK_FNGramCPT.h"
#include "GMTK_GMParms.h"
#include "GMTK_FileParser.h"
#include "GMTK_HidDiscRV.h"
#include "fileParser.h"
#include "logp.h"
#include "rand.h"
#include "error.h"


RAND rnd;
GMParms GM_Parms;
ObservationMatrix globalObservationMatrix;


void testFLM0() {
	// read in master file
	{
		iDataStreamFile pf("testFlm.master", false, true);
		GM_Parms.read(pf);
	}
	GM_Parms.finalizeParameters();

	// read in structure file
	FileParser fp("testFlm.str");

	// parse the file
	fp.parseGraphicalModel();
	// create the rv variable objects
	fp.createRandomVariableGraph();
	// Make sure that there are no directed loops in the graph.
	fp.ensureValidTemplate();
	fp.associateWithDataParams(FileParser::noAllocate);

	// get ngram pointer
	if ( GM_Parms.fngramImpsMap.find(string("fngram")) == GM_Parms.fngramImpsMap.end() ) {
		printf("cannot find fngram in FNGramImp\n");
		return;
	}

	if ( GM_Parms.fngramCptsMap.find(string("fngram")) == GM_Parms.fngramCptsMap.end() ) {
		printf("cannot find fngram in FNGramCPT\n");
		return;
	}

	FNGramCPT *fngram = GM_Parms.fngramCpts[GM_Parms.fngramCptsMap[string("fngram")]];

	printf("fngram pointer: %x\n", fngram);

	// read in vocabulary
	Vocab voc_w(27680), voc_a(3213);
	voc_w.read("word_flm.voc.w");
	voc_a.read("word_flm.voc.a");

	RVInfo wordInfo, aInfo;
	HidDiscRV word(wordInfo, 2, 27680), word1(wordInfo, 1, 27680), word2(wordInfo, 0, 27680), a(aInfo, 2, 3213);
	vector<RV*> parents(3);
	parents[0] = &word1;
	parents[1] = &word2;
	parents[2] = &a;

	logpr prob;

	//p( W-do | W-i,<s>,A-[sil])      = [0x7 gram] 0.0193205 [ -1.71398 ]
	word1.val = voc_w.index("W-i");
	word2.val = voc_w.index("<s>");
	a.val = voc_a.index("A-[sil]");
	word.val = voc_w.index("W-do");

	prob = fngram->probGivenParents(parents, &word);
	printf("P(W-do | W-i,<s>,A-[sil]) = %f\n", prob.val() / M_LN10);

	//p( W-i | <s>,<s>,A-head)        = [0x3 gram] 0.0199197 [ -1.70072 ]
	word1.val = voc_w.index("<s>");
	word2.val = voc_w.index("<s>");
	a.val = voc_a.index("A-head");
	word.val = voc_w.index("W-i");

	prob = fngram->probGivenParents(parents, &word);
	printf("P(W-i | <s>,<s>,A-head) = %f\n", prob.val() / M_LN10);

	//p( W-and | W-not,W-am,A-[sil])  = [0x1 gram] 0.00588722 [ -2.23009 ]
	word1.val = voc_w.index("W-not");
	word2.val = voc_w.index("W-am");
	a.val = voc_a.index("A-[sil]");
	word.val = voc_w.index("W-and");

	prob = fngram->probGivenParents(parents, &word);
	printf("P(W-and | W-not,W-am,A-[sil]) = %f\n", prob.val() / M_LN10);

	//p( W-i | <unk>,<s>,A-[sil])     = [0x0 gram] 0.0107713 [ -1.96773 ]
	word1.val = voc_w.index("<unk>");
	word2.val = voc_w.index("<s>");
	a.val = voc_a.index("A-[sil]");
	word.val = voc_w.index("W-i");

	prob = fngram->probGivenParents(parents, &word);
	printf("P(W-i | <unk>,<s>,A-[sil]) = %f\n", prob.val() / M_LN10);

	printf("try testing iterators\n");
	FNGramCPT::iterator it;
	fngram->becomeAwareOfParentValuesAndIterBegin(parents, it, &word, prob);
	printf("%x, word %d, prob %f\n", it.internalStatePtr, word.val, prob.val() / M_LN10);
	while ( fngram->next(it, prob) ) {
		if ( word.val == voc_w.index("W-i") )
			break;
	}
	printf("%x, word %d, prob %f\n", it.internalStatePtr, word.val, prob.val() / M_LN10);

	word1.val = voc_w.index("<s>");
	word2.val = voc_w.index("<s>");
	a.val = voc_a.index("A-head");
	FNGramCPT::iterator it2;
	fngram->becomeAwareOfParentValuesAndIterBegin(parents, it2, &word, prob);
	printf("%x, word %d, prob %f\n", it2.internalStatePtr, word.val, prob.val() / M_LN10);
	fngram->next(it2, prob);
	printf("%x, word %d, prob %f\n", it2.internalStatePtr, word.val, prob.val() / M_LN10);
}


void testFLM1() {
	// read in master file
	{
		iDataStreamFile pf("testFlm.master", false, true);
		GM_Parms.read(pf);
	}
	GM_Parms.finalizeParameters();
	if ( GM_Parms.fngramImpsMap.find(string("fngram")) == GM_Parms.fngramImpsMap.end() ) {
		printf("cannot find fngram in FNGramImp\n");
		return;
	}

	// read in structure file
	FileParser fp("testFlm2.str");

	// parse the file
	fp.parseGraphicalModel();
	// create the rv variable objects
	fp.createRandomVariableGraph();
	// Make sure that there are no directed loops in the graph.
	fp.ensureValidTemplate();
	fp.associateWithDataParams(FileParser::noAllocate);

	// get ngram pointer
	printf("number of fngram imps %d\n", GM_Parms.fngramImps.size());
	printf("number of fngrams %d\n", GM_Parms.fngramCpts.size());
	if ( GM_Parms.fngramCptsMap.find(string("fngram:0,1")) == GM_Parms.fngramCptsMap.end() ) {
		printf("cannot find fngram in FNGramCPT\n");
		return;
	}

	FNGramCPT *fngram = GM_Parms.fngramCpts[GM_Parms.fngramCptsMap[string("fngram:0,1")]];

	printf("fngram pointer: %x\n", fngram);

	// read in vocabulary
	Vocab voc_w(27680);
	voc_w.read("word_flm.voc.w");

	RVInfo wordInfo;
	HidDiscRV word(wordInfo, 2, 27680), word1(wordInfo, 1, 27680), word2(wordInfo, 0, 27680);
	vector<RV*> parents(2);
	parents[0] = &word1;
	parents[1] = &word2;

	logpr prob;

	//p( W-do | W-i,<s>)      = [0x7 gram] [-1.741229]
	word1.val = voc_w.index("W-i");
	word2.val = voc_w.index("<s>");
	word.val = voc_w.index("W-do");

	prob = fngram->probGivenParents(parents, &word);
	printf("P(W-do | W-i,<s>) = %f\n", prob.val() / M_LN10);

	//p( W-i | <s>,<s>)        = [0x3 gram] [-1.367162]
	word1.val = voc_w.index("<s>");
	word2.val = voc_w.index("<s>");
	word.val = voc_w.index("W-i");

	prob = fngram->probGivenParents(parents, &word);
	printf("P(W-i | <s>,<s>) = %f\n", prob.val() / M_LN10);

	//p( W-and | W-not,W-am)  = [0x1 gram] [-2.007548]
	word1.val = voc_w.index("W-not");
	word2.val = voc_w.index("W-am");
	word.val = voc_w.index("W-and");

	prob = fngram->probGivenParents(parents, &word);
	printf("P(W-and | W-not,W-am) = %f\n", prob.val() / M_LN10);

	//p( W-i | <unk>,<s>)     = [0x0 gram] [-1.96773]
	word1.val = voc_w.index("<unk>");
	word2.val = voc_w.index("<s>");
	word.val = voc_w.index("W-i");

	prob = fngram->probGivenParents(parents, &word);
	printf("P(W-i | <unk>,<s>) = %f\n", prob.val() / M_LN10);
}


int main() {
	ieeeFPsetup();
	set_new_handler(memory_error);

	testFLM0();
}
