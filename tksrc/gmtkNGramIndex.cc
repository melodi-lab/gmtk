/*
 * "Copyright 2001, University of Washington and International Business Machines Corporation. All Rights Reserved
 *
 *    Written by Gang Ji <gang@ee.washington.edu>
 *
 *    The goal of this code is to provide files for fast loading of an ARPA language file.
 * Instead of using
 *          bigram
 *          1 % number of parents
 *          VOCAB_SIZE VOCAB_SIZE % cards
 *          ./DATA/bigram.arpa vocab % ARPA lm file and vocabulary object
 * the user can use
 *          bigram
 *          1 % number of parents
 *          VOCAB_SIZE VOCAB_SIZE % cards
 *          ./DATA/bigram.arpa.idx [ascii] % ARPA lm indexing file
 * or
 *          bigram
 *          1 % number of parents
 *          VOCAB_SIZE VOCAB_SIZE % cards
 *          ./DATA/bigram.arpa.idx [binary] % ARPA lm indexing file
 *
 *
 * Copyright (c) 2001, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about the suitability of this software
 * for any purpose. It is provided "as is" without express or implied warranty.
 *
*/


#include <cstdlib>
#include <cstdio>

#include "GMTK_Vocab.h"
#include "GMTK_NGramCPT.h"
#include "GMTK_GMParms.h"
#include "GMTK_ObservationMatrix.h"
#include "rand.h"
#include "arguments.h"
#include "ieeeFPsetup.h"
#include "error.h"
#include "version.h"

VCID("$Header$");


/*
 * command line arguments
 */
static char * lmFile = NULL;
static char* vocabFile = NULL;
static bool outBin = false;
static bool print_version_and_exit = false;


Arg Arg::Args[] = {

	/////////////////////////////////////////////////////////////
	// input parameter/structure file handling

	Arg("lmFile", Arg::Req, lmFile, "Input ARPA lm file"),
	
	Arg("vocab", Arg::Req, vocabFile, "vocab file"),

	Arg("outBin", Arg::Opt, outBin, "Use binary for output index file"),

	Arg("version", Arg::Opt, print_version_and_exit, "Print GMTK version number and exit."),

	// final one to signal the end of the list
	Arg()

};


/*
 * definition of needed global arguments
 */
RAND rnd;
GMParms GM_Parms;
ObservationMatrix globalObservationMatrix;


int main(int argc, char *argv[]) {
#if 0
	////////////////////////////////////////////
	// set things up so that if an FP exception
	// occurs such as an "invalid" (NaN), overflow
	// or divide by zero, we actually get a FPE
	ieeeFPsetup();

	////////////////////////////////////////////
	// parse arguments
	Arg::parse(argc,argv);

	if (print_version_and_exit) {
		printf("%s\n",gmtk_version_id);
	}
	
	// figure out how many words in the vocab file
	unsigned card = 0;
	unsigned len = 1024;
	char * word = new char[len];
	FILE *fp = fopen(vocabFile, "r");
	if ( fp == NULL )
		error("cannot open file %s", vocabFile);
	while ( ! feof(fp) ) {
		if ( fscanf(fp, "%s", word) > 0 )
			++card;
	}
	fclose(fp);

	// figure out ngram order
	unsigned order = 0;
	if ( (fp = fopen(lmFile, "r")) == NULL )
		qerror("cannot open file %s", lmFile);

	do {
		if ( getline(&word, &len, fp) < 0 )
			error("wrong ARPA format in %s", lmFile);
	} while ( strstr(word, "\\data\\") != word );

	do {
		if ( getline(&word, &len, fp) < 0 )
			error("wrong ARPA format in %s", lmFile);
		if ( strstr(word, "ngram") != NULL && strchr(word, '=') != NULL )
			++order;
	} while ( word[0] != '\\' );

	fclose(fp);
	delete [] word;

	// read in vocab file
	Vocab vocab(card);
	vocab.read(vocabFile);

	// read in lm
	NGramCPT ngram;
	ngram.setNumParents(order - 1);
	for ( int i = 0; i < (int)order; i++ ) {
		ngram.setNumCardinality(i, card);
	}

	fprintf(stderr, "reading %d-order lm file...\n", order);
	ngram.read(lmFile, vocab);

	// save the n-gram into data for next fast reading
	char * indexFile = new char [strlen(lmFile) + 10];
	strcpy(indexFile, lmFile);
	strcat(indexFile, ".idx");

	fprintf(stderr, "saving index file...\n");
	oDataStreamFile ofs(indexFile, outBin);
	ngram.writeNGramIndexFile(ofs);

	delete [] indexFile;

	return 0;
#endif
}
