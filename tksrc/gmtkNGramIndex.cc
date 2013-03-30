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
 * Copyright (C) 2001 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 *
 *
*/

#if HAVE_CONFIG_H
#include <config.h>
#endif
#if HAVE_HG_H
#include "hgstamp.h"
#endif

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

VCID(HGID)



/*
 * command line arguments
 */
static char * lmFile = NULL;
static char* vocabFile = NULL;
static bool outBin = false;


Arg Arg::Args[] = {

	/////////////////////////////////////////////////////////////
	// input parameter/structure file handling

        Arg("\n*** Input files ***\n"),

	Arg("lmFile", Arg::Req, lmFile, "Input ARPA language model file"),
	
	Arg("vocab", Arg::Req, vocabFile, "vocab file"),

        Arg("\n*** Output format ***\n"),

	Arg("outBin", Arg::Opt, outBin, "Use binary for output index file"),

	// final one to signal the end of the list
	Arg()

};


/*
 * definition of needed global arguments
 */
RAND rnd;
GMParms GM_Parms;
ObservationMatrix obsMatrix;
ObservationMatrix *globalObservationMatrix = &obsMatrix;


int main(int argc, char *argv[]) {
	////////////////////////////////////////////
	// set things up so that if an FP exception
	// occurs such as an "invalid" (NaN), overflow
	// or divide by zero, we actually get a FPE
	ieeeFPsetup();
	set_new_handler(memory_error);

	////////////////////////////////////////////
	// parse arguments
	bool parse_was_ok = Arg::parse(argc,argv,
"\nThis program indexes ARPA language model files to make them more efficient\n");

	if(!parse_was_ok) {
	  Arg::usage(); exit(-1);
	}


	// figure out how many words in the vocab file
	unsigned card = 0;
	size_t len = 1024;
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
		error("cannot open file %s", lmFile);

	do {
#if defined(HAVE_GETLINE)
                if ( getline(&word, &len, fp) < 0 )
#else
		if ( fgets(word, len, fp) == NULL )
#endif
			error("wrong ARPA format in %s", lmFile);
	} while ( strstr(word, "\\data\\") != word );

	do {
#if defined(HAVE_GETLINE)
                if ( getline(&word, &len, fp) < 0 )
#else
		if ( fgets(word, len, fp) == NULL )
#endif
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
}
