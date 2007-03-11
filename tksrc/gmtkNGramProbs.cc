/*
 * "Copyright 2001, University of Washington and International Business Machines Corporation. All Rights Reserved
 *
 *    Written by Jeff Bilmes and Geoffrey Zweig
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

/*
 * This program converts from ascii trainable parameters to binary
 * and vice versa.
 *
 *
 */


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <float.h>
#include <assert.h>

#include "general.h"
#include "error.h"
#include "rand.h"
#include "arguments.h"
#include "ieeeFPsetup.h"
//#include "spi.h"
#include "version.h"

VCID("$Header$");

#include "GMTK_FileParser.h"
#include "GMTK_RandomVariable.h"
#include "GMTK_DiscreteRandomVariable.h"
#include "GMTK_ContinuousRandomVariable.h"
#include "GMTK_GM.h"
#include "GMTK_GMParms.h"
#include "GMTK_ObservationMatrix.h"
#include "GMTK_MixtureCommon.h"
#include "GMTK_GaussianComponent.h"
#include "GMTK_MeanVector.h"
#include "GMTK_DiagCovarVector.h"
#include "GMTK_DlinkMatrix.h"
#include "GMTK_NGramCPT.h"


/*
 * command line arguments
 */

char *inputMasterFile = NULL;

char *outputProbFile = NULL;

char *cppCommandOptions = NULL;

bool print_version_and_exit = false;

Arg Arg::Args[] = {

	/////////////////////////////////////////////////////////////
	// input parameter/structure file handling

	Arg("inputMasterFile", Arg::Opt, inputMasterFile, "Input file of multi-level master CPP processed GM input parameters"),

	Arg("outputProbFile", Arg::Opt, outputProbFile, "Output file of probabilities"),

	Arg("cppCommandOptions", Arg::Opt, cppCommandOptions, "Command line options to give to cpp"),

	Arg("version", Arg::Opt, print_version_and_exit, "Print GMTK version number and exit."),

	// final one to signal the end of the list
	Arg()

};

/*
 * definition of needed global arguments
 */
RAND rnd(false);
GMParms GM_Parms;
ObservationMatrix globalObservationMatrix;

int
main(int argc,char*argv[]) {
	////////////////////////////////////////////
	// set things up so that if an FP exception
	// occurs such as an "invalid" (NaN), overflow
	// or divide by zero, we actually get a FPE
	ieeeFPsetup();

	////////////////////////////////////////////
	// parse arguments
	bool parse_was_ok = Arg::parse(argc, (char**)argv);

	if ( print_version_and_exit ) {
		printf("%s\n", gmtk_version_id);
		exit(0);
	}

	if ( ! parse_was_ok ) {
		Arg::usage();
		exit(-1);
	}


	MixtureCommon::checkForValidRatioValues();
	MeanVector::checkForValidValues();
	DiagCovarVector::checkForValidValues();
	DlinkMatrix::checkForValidValues();

	/////////////////////////////////////////////
	if ( inputMasterFile == NULL ) {
		warning("ERROR: need to specify command line parameters inputMasterFile");
		Arg::usage();
		error("");
	} else if ( outputProbFile == NULL ) {
		warning("ERROR: need to specify command line parameters outputProbFile");
		Arg::usage();
		error("");
	}

	////////////////////////////////////////////
	iDataStreamFile pf(inputMasterFile,false,true,cppCommandOptions);
	GM_Parms.read(pf);

	GM_Parms.finalizeParameters();

	printf("Finished reading in all parameters and structures\n");
	printf("Total number of trainable parameters in input files = %u\n", GM_Parms.totalNumberParameters());

	GM_Parms.markUsedMixtureComponents();

	///////////////////////////////////////////
	// find the factored ngrams
	if ( GM_Parms.ngramCpts.size() != 1 )
		error("can only have one ngram in master file");

	NGramCPT *ngram = GM_Parms.ngramCpts[0];
	ngram->printProbs(outputProbFile);

	exit_program_with_status(0);
}
