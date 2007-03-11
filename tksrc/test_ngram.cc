/***************************************************************************
                          main.cpp  -  description
                             -------------------
    begin                : Mon May 17 16:14:21 PDT 2004
    copyright            : (C) 2004 by Gang Ji
    email                : gang@ee.washington.edu
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "GMTK_Vocab.h"
#include "GMTK_NGramCPT.h"
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

/*
void testVocab() {
	Vocab *vocab = new Vocab(13);
	
	vocab->read("w2.dct");
	
	std::cout << "vocab size " << vocab->size() << std::endl;
	std::cout << "index of oh " << vocab->index("oh") << std::endl;
	std::cout << "index of foo " << vocab->index("foo") << std::endl;
	
	delete vocab;
}


void testUnigram() {
	NGramCPT unigram;
	unigram.setNumParents(0);
	unigram.setNumCardinality(0, 13);
	
	Vocab *vocab = new Vocab(13);
	vocab->read("w2.dct");
	
	unigram.read("unigram.lm", *vocab);

	logpr p;
	
	//-1.218495       five
	p = unigram.probGivenParents(vocab->index("five"));
	std::cout << "five " << p.val() << std::endl;
	
	// -infty          foo
	p = unigram.probGivenParents(vocab->index("foo"));
	std::cout << "foo " << p.val() << std::endl;

	delete vocab;
}

void testBigram() {
	NGramCPT bigram;
	bigram.setNumParents(1);
	bigram.setNumCardinality(1, 14);
	
	Vocab *vocab = new Vocab(13);
	vocab->read("w2.dct");
	
	bigram.read("bigram.lm", *vocab);
	
	std::vector<int> context(1), card(1);
	logpr p;
	
	card[0] = 14;
	context[0] = vocab->index("five");
	bigram.setNumCardinality(0, card[0]);
	bigram.becomeAwareOfParentValues(context, card);

	// -infty <s> foo
	p = bigram.probGivenParents(vocab->index("foo"));
	std::cout << "five foo " << p.val() << std::endl;
	
	p = bigram.probGivenParents(vocab->index("four"));
	std::cout << "five four " << p.val() << std::endl;
	
	context[0] = vocab->index("foo");
	bigram.becomeAwareOfParentValues(context, card);
	p = bigram.probGivenParents(vocab->index("</s>"));
	std::cout << "foo </s> " << p.val() << std::endl;

	delete vocab;
}


void testTrigram() {
	Vocab *vocab = new Vocab(13);
	vocab->read("w2.dct");
	
	NGramCPT trigram;
	trigram.setNumParents(2);
	trigram.setNumCardinality(0, 14);
	trigram.setNumCardinality(1, 14);
	trigram.setNumCardinality(2, 14);
	
	trigram.read("trigram.lm", *vocab);
	
	logpr p;
	std::vector<int> context(2), card(2);
	card[0] = 14; card[1] = 14;
	
	context[0] = vocab->index("<s>");
	context[1] = vocab->index("one");
	trigram.becomeAwareOfParentValues(context, card);
	p = trigram.probGivenParents(vocab->index("foo"));
	std::cout << "<s> one foo " << p.val() << std::endl;
	
	context[0] = vocab->index("five");
	context[1] = vocab->index("four");
	trigram.becomeAwareOfParentValues(context, card);
	p = trigram.probGivenParents(vocab->index("oh"));
	std::cout << "five four oh " << p.val() << std::endl;

	context[0] = vocab->index("five");
	context[1] = vocab->index("two");
	trigram.becomeAwareOfParentValues(context, card);
	p = trigram.probGivenParents(vocab->index("eight"));
	std::cout << "five two eight " << p.val() << std::endl;

	delete vocab;
}
*/


void testBigTrigram() {
	// read in master file
	{
		iDataStreamFile pf("testLm.master", false, true);
		GM_Parms.read(pf);
	}
	GM_Parms.finalizeParameters();

	// read in structure file
	FileParser fp("testLm.str");

	// parse the file
	fp.parseGraphicalModel();
	// create the rv variable objects
	fp.createRandomVariableGraph();
	// Make sure that there are no directed loops in the graph.
	fp.ensureValidTemplate();
	fp.associateWithDataParams(FileParser::noAllocate);

	// get ngram pointer
	if ( GM_Parms.ngramCptsMap.find(string("ngram")) == GM_Parms.ngramCptsMap.end() ) {
		printf("cannot find fngram in NGramCPT\n");
		return;
	}

	NGramCPT *ngram = GM_Parms.ngramCpts[GM_Parms.ngramCptsMap[string("ngram")]];

	printf("ngram pointer: %x\n", ngram);

	Vocab *vocab = new Vocab(27682);
	vocab->read("w.dct");

	logpr prob;
	RVInfo wordInfo;
	HidDiscRV word(wordInfo, 2, 27682), word1(wordInfo, 1, 27682), word2(wordInfo, 0, 27682);
	vector<RV*> parents(2);
	parents[0] = &word1;
	parents[1] = &word2;

	//-1.840867       ability to control => -4.23875
	word1.val = vocab->index("ability");
	word2.val = vocab->index("to");

	word.val = vocab->index("control");
	prob = ngram->probGivenParents(parents, &word);
	std::cout << "prob: ability to control " << prob.val() << std::endl;

	//-2.38328        to look -0.7311508
	//-0.4941311      ability to      -0.1190469
	// should be -2.5 => -5.76182
	word.val = vocab->index("look");
	prob = ngram->probGivenParents(parents, &word);
	std::cout << "prob: ability to look " << prob.val() << std::endl;

	//-4.775465       loop    -0.2603517
	//-2.100044       to      -0.8924304
	//-0.4941311      ability to      -0.1190469
	// should be -5.7 => -13.3249
	word.val = vocab->index("loop");
	prob = ngram->probGivenParents(parents, &word);
	std::cout << "prob: ability to loop " << prob.val() << std::endl;

	std::cout << "testing some from the ppl file...\n";

	// p( halloween | it is )  = [1gram] 5.81301e-07 [ -6.2356 ] => -14.358
	word1.val = vocab->index("it");
	word2.val = vocab->index("is");
	word.val = vocab->index("halloween");
	prob = ngram->probGivenParents(parents, &word);
	std::cout << "prob: it is halloween " << prob.val() << std::endl;

	// p( think | (%hesitation) i )       = [2gram] 0.0557179 [ -1.254 ] => -2.88745
	std::cout << "hesitation: " << (word1.val = vocab->index("(%hesitation)")) << std::endl;
	word2.val = vocab->index("i");
	word.val = vocab->index("think");
	prob = ngram->probGivenParents(parents, &word);
	std::cout << "prob: (%hesitation) i think " << prob.val() << std::endl;

	printf("try testing iterators\n");
	NGramCPT::iterator it;
	ngram->becomeAwareOfParentValuesAndIterBegin(parents, it, &word, prob);
	printf("%x, word %d, prob %f\n", it.internalStatePtr, word.val, prob.val());
	while ( ngram->next(it, prob) ) {
		if ( word.val == vocab->index("think") )
			break;
	}
	printf("%x, word %d, prob %f\n", it.internalStatePtr, word.val, prob.val());

	word1.val = vocab->index("<s>");
	word2.val = vocab->index("<s>");
	NGramCPT::iterator it2;
	ngram->becomeAwareOfParentValuesAndIterBegin(parents, it2, &word, prob);
	printf("%x, word %d, prob %f\n", it2.internalStatePtr, word.val, prob.val());
	ngram->next(it2, prob);
	printf("%x, word %d, prob %f\n", it2.internalStatePtr, word.val, prob.val());

	delete vocab;
}


int main(int argc, char *argv[]) {
	rnd.seed();

	std::cout << "testing NGramCPT...\n";
	//testUnigram();
	//testBigram();
	//testTrigram();
	testBigTrigram();
	std::cout << "end testing NGramCPT.\n\n";

	return EXIT_SUCCESS;
}
