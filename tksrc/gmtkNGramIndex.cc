/*
 * "Copyright 2001, University of Washington and International Business Machines Corporation. All Rights Reserved
 *
 *    Written by Gang Ji <gang@ee.washington.edu>
 *
 * NO WARRANTY
 * THE PROGRAM IS PROVIDED ON AN "AS IS" BASIS, WITHOUT WARRANTIES OR
 * CONDITIONS OF ANY KIND, EITHER EXPRESS OR IMPLIED INCLUDING, WITHOUT
 * LIMITATION, ANY WARRANTIES OR CONDITIONS OF TITLE, NON-INFRINGEMENT,
 * MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. Each Recipient is
 * solely responsible for determining the appropriateness of using the Program
 * and assumes all risks associated with such use, including but not limited
 * to the risks and costs of program errors, compliance with applicable laws,
 * damage to or loss of data, programs or equipment, and unavailability or
 * interruption of operations.

 * DISCLAIMER OF LIABILITY
 * THE UNIVERSITY OF WASHINGTON, INTERNATIONAL BUSINESS MACHINES CORPORATION,
 * JEFF BILMES AND GEOFFREY ZWEIG SHALL NOT HAVE ANY LIABILITY FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING WITHOUT LIMITATION LOST PROFITS), HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE  OF
 * THE PROGRAM, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGES."
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
		error("cannot open file %s", lmFile);

	do {
		if ( getline(&word, &len, fp) < 0 )
			error("wrong ARPA format in %s", lmFile);
	} while ( word[0] != '\\' );
	if ( strstr(word, "\\data\\") == NULL )
		error("wrong ARPA format in %s", lmFile);
	
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

	fprintf(stderr, "reading lm file...\n");
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
