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

/*-
 * gmtkEMtrain.cc
 *     Train up a GM using EM
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

VCID("$Header$");

#include "GMTK_FileParser.h"
#include "GMTK_RandomVariable.h"
#include "GMTK_DiscreteRandomVariable.h"
#include "GMTK_ContinuousRandomVariable.h"
#include "GMTK_GM.h"
#include "GMTK_GMParms.h"
#include "GMTK_ObservationMatrix.h"
#include "GMTK_MixGaussiansCommon.h"
#include "GMTK_GaussianComponent.h"
#include "GMTK_MeanVector.h"
#include "GMTK_DiagCovarVector.h"


/*
 * command line arguments
 */
bool seedme = false;
float pruneRatio=0.0;
char *obsFileName=NULL;
char *strFileName=NULL;
char *outFileName=NULL;
bool writeParametersBeforeEachEMIteration=true;
bool binOutFile=false;
char *parmsFileName=NULL;
bool binParmsFile=false;
char *parmsPtrFileName=NULL;
unsigned maxEMIterations=3;
bool randomizeParams = true;
bool enem = false;
double mcvr = 1e20;
double mcsr = 1e10;
double varFloor = 1e-10;
char *trrng_str="all";
float lldp = 0.0;
float mnlldp = 0.01;
char *storeAccFile = NULL;
char *progStoreAccFile = NULL;
char *loadAccFile = NULL;
char *loadAccRange = NULL;
char *llStoreFile = NULL;
bool accFileIsBinary = true;
int startSkip = 0;
int endSkip = 0;

ARGS ARGS::Args[] = {

 ARGS("obsFile",ARGS::Req,obsFileName,"File containing observations"),

 ARGS("parmsPtrFile",ARGS::Opt,parmsPtrFileName,"Multi-level GM Parms File"),

 ARGS("parmsFile",ARGS::Opt,parmsFileName,"Single-level GM Parms File"),
 ARGS("binParmsFile",ARGS::Opt,binParmsFile,"Is Single-level GM Parms File binary? (def=false)"),

 ARGS("outFileName",ARGS::Opt,outFileName,"File to place output parametes"),
 ARGS("binOutFile",ARGS::Opt,binOutFile,"Output parametes binary? (def=false)"),
 ARGS("wpbeei",ARGS::Opt,writeParametersBeforeEachEMIteration,
      "Write Paramaeters Before Each EM Iteration? (def=true)"),

 ARGS("strFile",ARGS::Req,strFileName,"GM Structure File"),

 ARGS("seed",ARGS::Opt,seedme,"Seed the RN generator"),

 ARGS("maxEmIters",ARGS::Opt,maxEMIterations,"Max number of EM iterations to do"),

 ARGS("pruneRatio",ARGS::Opt,pruneRatio,"Pruning Ratio, values less than this*max are pruned"),

 // support for splitting and vanishing
 ARGS("mcvr",ARGS::Opt,MixGaussiansCommon::mixCoeffVanishRatio,"Mixture Coefficient Vanishing Ratio"),
 ARGS("mcsr",ARGS::Opt,MixGaussiansCommon::mixCoeffSplitRatio,"Mixture Coefficient Splitting Ratio"),

 ARGS("meanCloneSTDfrac",ARGS::Opt,MeanVector::cloneSTDfrac,"Fraction of mean to use for STD in mean clone"),
 ARGS("covarCloneSTDfrac",ARGS::Opt,DiagCovarVector::cloneSTDfrac,"Fraction of var to use for STD in covar clone"),

 ARGS("varFloor",ARGS::Opt,varFloor,"Variance Floor"),
 ARGS("lldp",ARGS::Opt,lldp,"Log Likelihood difference percentage for termination"),
 ARGS("mnlldp",ARGS::Opt,mnlldp,"Absolute value of max negative Log Likelihood difference percentage for termination"),

 ARGS("trrng",ARGS::Opt,trrng_str,"Range to train over segment file"),

 ARGS("storeAccFile",ARGS::Opt,storeAccFile,"Store accumulators file"),
 ARGS("loadAccFile",ARGS::Opt,loadAccFile,"Load accumulators file"), 
 ARGS("loadAccRange",ARGS::Opt,loadAccRange,"Load accumulators file range"), 
 ARGS("llStoreFile",ARGS::Opt,llStoreFile,"File to store previous sum LL's"), 
 ARGS("accFileIsBinary",ARGS::Opt,accFileIsBinary,"Binary accumulator files (def true)"), 

 ARGS("startSkip",ARGS::Opt,startSkip,"Frames to skip at beginning (i.e., first frame is buff[startSkip])"),
 ARGS("endSkip",ARGS::Opt,endSkip,"Frames to skip at end (i.e., last frame is buff[len-1-endSkip])"),


 ARGS("random",ARGS::Opt,randomizeParams,"Randomize the parameters"),
 ARGS("enem",ARGS::Opt,enem,"Run enumerative EM"),

 // final one to signal the end of the list
 ARGS()
};

/*
 * definition of needed global arguments
 */
RAND rnd(seedme);
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
  ARGS::parse(argc,argv);

  ////////////////////////////////////////////
  // check for valid argument values.
  MixGaussiansCommon::checkForValidRatioValues();
  MeanVector::checkForValidValues();
  DiagCovarVector::checkForValidValues();
  if (lldp < 0.0 || mnlldp < 0.0)
    error("lldp & mnlldp must be >= 0");
  if (pruneRatio < 0.0)
    error("pruneRatio must be >= 0");
  if (startSkip < 0 || endSkip < 0)
    error("startSkip/endSkip must be >= 0");
  if (!parmsFileName && !parmsPtrFileName) 
    error("ERROR: arguments must specify a parameter file.");  

  ////////////////////////////////////////////
  // set global variables/change global state from args
  GaussianComponent::setVarianceFloor(varFloor);
  if (seedme)
    rnd.seed();

  /////////////////////////////////////////////
  // read in all parameters
  if (parmsFileName) {
    // flat, where everything is contained in one file.
    iDataStreamFile parmsFile(parmsFileName,binParmsFile);
    GM_Parms.readAll(parmsFile);
  }
  if (parmsPtrFileName) {
    // this file is always ASCII.
    iDataStreamFile parmsFile(parmsPtrFileName);
    GM_Parms.read(parmsFile);
  }


  /////////////////////////////
  // read in the structure of the GM, this will
  // die if the file does not exist.
  FileParser fp(strFileName);
  // parse the file
  fp.parseGraphicalModel();
  // create the rv variable objects
  fp.createRandomVariableGraph();
  // make sure that there are no directed loops in the graph
  // by imposing the S,SE,E,NE constrains
  fp.ensureS_SE_E_NE();
  // link the RVs with the parameters that are contained in
  // the bn1_gm.dt file.
  fp.associateWithDataParams();
  // now associate the RVs with a GM
  GMTK_GM gm;
  fp.addVariablesToGM(gm);
  gm.verifyTopologicalOrder();

  ////////////////////////////////////
  // set up the observation stream
  gm.setExampleStream(obsFileName);


  gm.GM2CliqueChain();
  // gm.showCliques();

  if (randomizeParams) {
    printf("NOTE: Randomizing initial parameters\n");
    GM_Parms.makeRandom();
    printf("NOTE: Writing iniial randomized parameters to random.gmp\n");
    oDataStreamFile of("random.gmp");
    GM_Parms.writeAll(of);
  }

  /////////////////////////////////////
  // finaly, start training.
  gm.setupForVariableLengthUnrolling(fp.firstChunkFrame(),fp.lastChunkFrame());
  if (enem) {
    warning("**** WARNING: Doing enumerative EM!!! ****");
    gm.enumerativeEM(maxEMIterations);
  } else {
    gm.cliqueChainEM(maxEMIterations, pruneRatio,
		     writeParametersBeforeEachEMIteration,
		     ((outFileName)?(string(outFileName)):("")),
		     binOutFile);
  }

  printf("____ PROGRAM ENDED SUCCESSFULLY AT ");
  print_date_string(stdout);
  printf(" ____\n");
  return 0;
}
