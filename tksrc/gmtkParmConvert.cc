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
#include "spi.h"
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


/*
 * command line arguments
 */

char *inputMasterFile=NULL;
char *outputMasterFile=NULL;

char *outputTrainableParameters=NULL;
bool binOutputTrainableParameters=false;

char *inputTrainableParameters=NULL;
bool binInputTrainableParameters=false;

double varFloor = 1e-10;
char *cppCommandOptions = NULL;

unsigned allocateDenseCpts=0;
bool seedme = false;

char *strFileName=NULL;

bool print_version_and_exit = false;

Arg Arg::Args[] = {

  /////////////////////////////////////////////////////////////
  // input parameter/structure file handling

  Arg("strFile",Arg::Opt,strFileName,"Optional Graphical Model Structure File"),
  Arg("inputMasterFile",Arg::Opt,inputMasterFile,"Input file of multi-level master CPP processed GM input parameters"),
  Arg("outputMasterFile",Arg::Opt,outputMasterFile,"Output file to place master CPP processed GM output parameters"),

  Arg("inputTrainableParameters",Arg::Opt,inputTrainableParameters,"File of only and all trainable parameters"),
  Arg("binInputTrainableParameters",Arg::Opt,binInputTrainableParameters,"Binary condition of trainable parameters file"),

  Arg("outputTrainableParameters",Arg::Opt,outputTrainableParameters,"File to place only and all trainable output parametes"),
  Arg("binOutputTrainableParameters",Arg::Opt,binOutputTrainableParameters,"Binary condition of output trainable parameters?"),

  Arg("cppCommandOptions",Arg::Opt,cppCommandOptions,"Command line options to give to cpp"),

  Arg("varFloor",Arg::Opt,varFloor,"Variance Floor"),
  Arg("floorVarOnRead",Arg::Opt,DiagCovarVector::floorVariancesWhenReadIn,
       "Floor the variances to varFloor when they are read in"),
  Arg("cptNormThreshold",Arg::Opt,CPT::normalizationThreshold,"Read error if |Sum-1.0|/card > norm_threshold"),


  Arg("seed",Arg::Opt,seedme,"Seed the random number generator"),
  Arg("allocateDenseCpts",Arg::Opt,allocateDenseCpts,"Automatically allocate any undefined CPTs. arg = 1 means use random initial CPT values. arg = 2, use uniform values"),

  Arg("version",Arg::Opt,print_version_and_exit,"Print GMTK version number and exit."),

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
main(int argc,char*argv[])
{
  ////////////////////////////////////////////
  // set things up so that if an FP exception
  // occurs such as an "invalid" (NaN), overflow
  // or divide by zero, we actually get a FPE
  ieeeFPsetup();

  ////////////////////////////////////////////
  // parse arguments
  Arg::parse(argc,argv);

  if (print_version_and_exit)
    printf("%s\n",gmtk_version_id);


  MixtureCommon::checkForValidRatioValues();
  MeanVector::checkForValidValues();
  DiagCovarVector::checkForValidValues();
  DlinkMatrix::checkForValidValues();

  ////////////////////////////////////////////
  // set global variables/change global state from args
  GaussianComponent::setVarianceFloor(varFloor);
  if (seedme)
    rnd.seed();
  /////////////////////////////////////////////

  if (strFileName == NULL || allocateDenseCpts == 0) {
    // check this only if the structure file is not given or
    // if we are not allocating Dense CPTs, since
    // in that case there is no way the user could specify
    // automatic allocation of such CPTs.
    if ((inputMasterFile == NULL) && (inputTrainableParameters == NULL)) {
      warning("ERROR: need to specify command line parameters inputMasterFile or inputTrainableParameters (or both) when no structure file is given");
      Arg::usage();
      error("");
    }
  }

  ////////////////////////////////////////////
  if (inputMasterFile != NULL) {
    iDataStreamFile pf(inputMasterFile,false,true,cppCommandOptions);
    GM_Parms.read(pf);
  }
  if (inputTrainableParameters != NULL) {
    // flat, where everything is contained in one file
    iDataStreamFile pf(inputTrainableParameters,binInputTrainableParameters,true,cppCommandOptions);
    GM_Parms.readTrainable(pf);
  }
  if (strFileName != NULL) {
    // load up the structure file as we might want
    // it to allocate some Dense CPTs.
    FileParser fp(strFileName,cppCommandOptions);

    // parse the file
    fp.parseGraphicalModel();
    // create the rv variable objects
    fp.createRandomVariableGraph();
    // make sure that there are no directed loops in the graph
    // by imposing the S,SE,E,NE constrains
    fp.ensureS_SE_E_NE();
    // link the RVs with the parameters that are contained in
    // the bn1_gm.dt file.
    if (allocateDenseCpts == 0)
      fp.associateWithDataParams(FileParser::noAllocate);
    else if (allocateDenseCpts == 1)
      fp.associateWithDataParams(FileParser::allocateRandom);
    else if (allocateDenseCpts == 2)
      fp.associateWithDataParams(FileParser::allocateUniform);
    else
      error("Error: command line argument '-allocateDenseCpts d', must have d = {0,1,2}\n");
  }

  printf("Finished reading in all parameters and structures\n");
  printf("Total number of trainable parameters in input files = %u\n",
	 GM_Parms.totalNumberParameters());

  GM_Parms.markUsedMixtureComponents();

  if (outputMasterFile != NULL) {
    GM_Parms.write(outputMasterFile);
  }
  if (outputTrainableParameters != NULL) {
    oDataStreamFile of(outputTrainableParameters,binOutputTrainableParameters);
    GM_Parms.writeTrainable(of);
  }

  exit_program_with_status(0);
}

