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
#include "GMTK_ProgramDefaultParms.h"

/*
 * command line arguments
 */
bool seedme = false;
float beam=-LZERO;

char *strFileName=NULL;

// char *outputTrainableParameters="outParms%d.gmp";
char *outputTrainableParameters=NULL;
bool binOutputTrainableParameters=false;
bool writeParametersAfterEachEMIteration=true;

char *inputMasterFile=NULL;
char *outputMasterFile=NULL;
char *inputTrainableParameters=NULL;
bool binInputTrainableParameters=false;
char *objsToNotTrainFile=NULL;

unsigned maxEMIterations=3;
bool randomizeParams = false;
bool enem = false;
double mcvr = 1e20;
double mcsr = 1e10;
double varFloor = GMTK_DEFAULT_VARIANCE_FLOOR;
char *trrng_str="all";
float lldp = 0.001;
float mnlldp = 0.01;

char *loadAccFile = NULL;
char *loadAccRange = NULL;
char *storeAccFile = NULL;
bool accFileIsBinary = true;

// file to store log likelihood of this iteration.
char *llStoreFile = NULL;

int startSkip = 0;
int endSkip = 0;

// observation file support

char *obsFileName;

#define MAX_NUM_OBS_FILES (3)
char *ofs[MAX_NUM_OBS_FILES] = { NULL, NULL, NULL }; 
unsigned nfs[MAX_NUM_OBS_FILES] = { 0, 0, 0 };
unsigned nis[MAX_NUM_OBS_FILES] = { 0, 0, 0 };
char *frs[MAX_NUM_OBS_FILES] = { "all", "all", "all" };
char *irs[MAX_NUM_OBS_FILES] = { "all", "all", "all" };
char *fmts[MAX_NUM_OBS_FILES] = { "pfile", "pfile", "pfile" };
bool iswps[MAX_NUM_OBS_FILES] = { false, false, false };

int show_cliques=0;

char *cppCommandOptions = NULL;

int bct=GMTK_DEFAULT_BASECASETHRESHOLD;
int ns=GMTK_DEFAULT_NUM_SPLITS;

unsigned allocateDenseCpts=0;

bool print_version_and_exit = false;

Arg Arg::Args[] = {

  /////////////////////////////////////////////////////////////
  // observation input file handling

  Arg("of1",Arg::Req,ofs[0],"Observation File 1"),
  Arg("nf1",Arg::Opt,nfs[0],"Number of floats in observation file 1"),
  Arg("ni1",Arg::Opt,nis[0],"Number of ints in observation file 1"),
  Arg("fr1",Arg::Opt,frs[0],"Float range for observation file 1"),
  Arg("ir1",Arg::Opt,irs[0],"Int range for observation file 1"),
  Arg("fmt1",Arg::Opt,fmts[0],"Format (htk,bin,asc,pfile) for observation file 1"),
  Arg("iswp1",Arg::Opt,iswps[0],"Endian swap condition for observation file 1"),


  Arg("of2",Arg::Opt,ofs[1],"Observation File 2"),
  Arg("nf2",Arg::Opt,nfs[1],"Number of floats in observation file 2"),
  Arg("ni2",Arg::Opt,nis[1],"Number of ints in observation file 2"),
  Arg("fr2",Arg::Opt,frs[1],"Float range for observation file 2"),
  Arg("ir2",Arg::Opt,irs[1],"Int range for observation file 2"),
  Arg("fmt2",Arg::Opt,fmts[1],"Format (htk,bin,asc,pfile) for observation file 2"),
  Arg("iswp2",Arg::Opt,iswps[1],"Endian swap condition for observation file 2"),


  Arg("of3",Arg::Opt,ofs[2],"Observation File 3"),
  Arg("nf3",Arg::Opt,nfs[2],"Number of floats in observation file 3"),
  Arg("ni3",Arg::Opt,nis[2],"Number of ints in observation file 3"),
  Arg("fr3",Arg::Opt,frs[2],"Float range for observation file 3"),
  Arg("ir3",Arg::Opt,irs[2],"Int range for observation file 3"),
  Arg("fmt3",Arg::Opt,fmts[2],"Format (htk,bin,asc,pfile) for observation file 3"),
  Arg("iswp3",Arg::Opt,iswps[2],"Endian swap condition for observation file 3"),


  /////////////////////////////////////////////////////////////
  // input parameter/structure file handling

  Arg("cppCommandOptions",Arg::Opt,cppCommandOptions,"Additional CPP command line"),

  Arg("inputMasterFile",Arg::Req,inputMasterFile,"Input file of multi-level master CPP processed GM input parameters"),
  Arg("outputMasterFile",Arg::Opt,outputMasterFile,"Output file to place master CPP processed GM output parameters"),

  Arg("inputTrainableParameters",Arg::Opt,inputTrainableParameters,"File of only and all trainable parameters"),
  Arg("binInputTrainableParameters",Arg::Opt,binInputTrainableParameters,"Binary condition of trainable parameters file"),

  Arg("objsNotToTrain",Arg::Opt,objsToNotTrainFile,"File listing trainable parameter objects to not train."),

  Arg("outputTrainableParameters",Arg::Opt,outputTrainableParameters,"File to place only and all trainable output parametes"),
  Arg("binOutputTrainableParameters",Arg::Opt,binOutputTrainableParameters,"Binary condition of output trainable parameters?"),


  Arg("wpaeei",Arg::Opt,writeParametersAfterEachEMIteration,
      "Write Parameters After Each EM Iteration Completes"),


  Arg("strFile",Arg::Req,strFileName,"Graphical Model Structure File"),



  /////////////////////////////////////////////////////////////
  // general files

  Arg("seed",Arg::Opt,seedme,"Seed the random number generator"),
  Arg("maxEmIters",Arg::Opt,maxEMIterations,"Max number of EM iterations to do"),
  Arg("beam",Arg::Opt,beam,"Beam width (less than max*exp(-beam) are pruned away)"),
  Arg("allocateDenseCpts",Arg::Opt,allocateDenseCpts,"Automatically allocate any undefined CPTs. arg = 1 means use random initial CPT values. arg = 2, use uniform values"),

  // support for splitting and vanishing
  Arg("mcvr",Arg::Opt,MixtureCommon::mixCoeffVanishRatio,"Mixture Coefficient Vanishing Ratio"),
  Arg("botForceVanish",Arg::Opt,MixtureCommon::numBottomToForceVanish,"Number of bottom mixture components to force vanish"),
  
  Arg("mcsr",Arg::Opt,MixtureCommon::mixCoeffSplitRatio,"Mixture Coefficient Splitting Ratio"),
  Arg("topForceSplit",Arg::Opt,MixtureCommon::numTopToForceSplit,"Number of top mixture components to force split"),

  Arg("meanCloneSTDfrac",Arg::Opt,MeanVector::cloneSTDfrac,"Fraction of mean to use for STD in mean clone"),
  Arg("covarCloneSTDfrac",Arg::Opt,DiagCovarVector::cloneSTDfrac,"Fraction of var to use for STD in covar clone"),
  Arg("dlinkCloneSTDfrac",Arg::Opt,DlinkMatrix::cloneSTDfrac,"Fraction of var to use for STD in covar clone"),

  Arg("cloneShareMeans",Arg::Opt,GaussianComponent::cloneShareMeans,"Gaussian component clone shares parent mean"),
  Arg("cloneShareCovars",Arg::Opt,GaussianComponent::cloneShareCovars,"Gaussian component clone shares parent covars"),
  Arg("cloneShareDlinks",Arg::Opt,GaussianComponent::cloneShareDlinks,"Gaussian component clone shares parent dlinks"),


  Arg("varFloor",Arg::Opt,varFloor,"Variance Floor"),
  Arg("floorVarOnRead",Arg::Opt,DiagCovarVector::floorVariancesWhenReadIn,
       "Floor the variances to varFloor when they are read in"),



  Arg("lldp",Arg::Opt,lldp,"Log Likelihood difference percentage for termination"),
  Arg("mnlldp",Arg::Opt,mnlldp,"Absolute value of max negative Log Likelihood difference percentage for termination"),

  Arg("trrng",Arg::Opt,trrng_str,"Range to train over segment file"),

  Arg("storeAccFile",Arg::Opt,storeAccFile,"Store accumulators file"),
  Arg("loadAccFile",Arg::Opt,loadAccFile,"Load accumulators file"), 
  Arg("loadAccRange",Arg::Opt,loadAccRange,"Load accumulators file range"), 
  Arg("llStoreFile",Arg::Opt,llStoreFile,"File to store previous sum LL's"), 
  Arg("accFileIsBinary",Arg::Opt,accFileIsBinary,"Binary accumulator files (def true)"), 

  Arg("startSkip",Arg::Opt,startSkip,"Frames to skip at beginning (i.e., first frame is buff[startSkip])"),
  Arg("endSkip",Arg::Opt,endSkip,"Frames to skip at end (i.e., last frame is buff[len-1-endSkip])"),
  
  Arg("cptNormThreshold",Arg::Opt,CPT::normalizationThreshold,"Read error if |Sum-1.0|/card > norm_threshold"),
  Arg("random",Arg::Opt,randomizeParams,"Randomize the parameters"),
  Arg("enem",Arg::Opt,enem,"Run enumerative EM"),

Arg("showCliques",Arg::Opt,show_cliques,"Show the cliques after the network has been unrolled k times."),


  Arg("numSplits",Arg::Opt,ns,"Number of splits to use in logspace recursion (>=2)."),

  Arg("baseCaseThreshold",Arg::Opt,bct,"Base case threshold to end recursion (>=2)."),

  Arg("componentCache",Arg::Opt,MixtureCommon::cacheComponentsInEmTraining,"Cache component probabilities during EM training, speeds things up but uses more memory."),

  Arg("version",Arg::Opt,print_version_and_exit,"Print GMTK version number and exit."),

  // final one to signal the end of the list
  Arg()

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
  set_new_handler(memory_error);

  ////////////////////////////////////////////
  // parse arguments
  Arg::parse(argc,argv);

  if (print_version_and_exit)
    printf("%s\n",gmtk_version_id);

  ////////////////////////////////////////////
  // check for valid argument values.
  int nfiles = 0;
  unsigned ifmts[MAX_NUM_OBS_FILES];
  for (int i=0;i<MAX_NUM_OBS_FILES;i++) {
    if (ofs[i] != NULL && nfs[i] == 0 && nis[i] == 0)
      error("ERROR: command line parameters must specify one of nf%d and ni%d as not zero",
	    i+1,i+1);
    nfiles += (ofs[i] != NULL);
    if (strcmp(fmts[i],"htk") == 0)
      ifmts[i] = HTK;
    else if (strcmp(fmts[i],"binary") == 0)
      ifmts[i] = RAWBIN;
    else if (strcmp(fmts[i],"ascii") == 0)
      ifmts[i] = RAWASC;
    else if (strcmp(fmts[i],"pfile") == 0)
      ifmts[i] = PFILE;
    else
      error("ERROR: Unknown observation file format type: '%s'\n",fmts[i]);
  }

  globalObservationMatrix.openFiles(nfiles,
				    (const char**)&ofs,
				    (const char**)&frs,
				    (const char**)&irs,
				    (unsigned*)&nfs,
				    (unsigned*)&nis,
				    (unsigned*)&ifmts,
				    (bool*)&iswps,
				    startSkip,
				    endSkip);

#if 0
  // for debugging
  for (int i=0;i<globalObservationMatrix.numSegments();i++) {
    printf("loading segment %d\n",i);
    globalObservationMatrix.loadSegment(i);
  }
#endif

  MixtureCommon::checkForValidRatioValues();
  MeanVector::checkForValidValues();
  DiagCovarVector::checkForValidValues();
  DlinkMatrix::checkForValidValues();
  if (lldp < 0.0 || mnlldp < 0.0)
    error("lldp & mnlldp must be >= 0");
  if (beam < 0.0)
    error("beam must be >= 0");
  if (startSkip < 0 || endSkip < 0)
    error("startSkip/endSkip must be >= 0");

  ////////////////////////////////////////////
  // set global variables/change global state from args
  GaussianComponent::setVarianceFloor(varFloor);
  if (seedme)
    rnd.seed();

  /////////////////////////////////////////////
  // read in all the parameters
  if (inputMasterFile) {
    // flat, where everything is contained in one file, always ASCII
    iDataStreamFile pf(inputMasterFile,false,true,cppCommandOptions);
    GM_Parms.read(pf);
  }
  if (inputTrainableParameters) {
    // flat, where everything is contained in one file
    iDataStreamFile pf(inputTrainableParameters,binInputTrainableParameters,true,cppCommandOptions);
    GM_Parms.readTrainable(pf);
  }

  /////////////////////////////
  // read in the structure of the GM, this will
  // die if the file does not exist.
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

  GM_Parms.markObjectsToNotTrain(objsToNotTrainFile,cppCommandOptions);
  // lastly, load the internal objects.
  GM_Parms.loadGlobal();

  printf("Finished reading in all parameters and structures\n");

  // make sure that all observation variables work
  // with the global observation stream.
  fp.checkConsistentWithGlobalObservationStream();
  // now associate the RVs with a GM
  GMTK_GM gm(&fp);
  fp.addVariablesToGM(gm);
  gm.verifyTopologicalOrder();

  ////////////////////////////////////
  // set up the observation stream
  gm.setExampleStream(obsFileName,trrng_str);
  GM_Parms.checkConsistentWithGlobalObservationStream();

  GM_Parms.setFirstUtterance( gm.trrng->min() ); 

  gm.setCliqueChainRecursion(ns, bct);

  gm.GM2CliqueChain();
  gm.setupForVariableLengthUnrolling(fp.firstChunkFrame(),fp.lastChunkFrame());
  if (show_cliques)
  {
    cout << "The cliques in the unrolled network are:\n";
    gm.setSize(show_cliques);
    gm.showCliques();
  }

  if (randomizeParams) {
    printf("### GMTK is randomizing all trainable parameters and writing them to random.gmp ####\n");
    GM_Parms.makeRandom();
    oDataStreamFile of("random.gmp");
    GM_Parms.writeTrainable(of);
  }

  /////////////////////////////////////
  // finaly, start training.
  logpr pruneRatio;
  pruneRatio.valref() = -beam;
  if (enem) {
    warning("******************************************");
    warning("**** WARNING: Doing enumerative EM!!! ****");
    warning("******************************************");
    gm.enumerativeEM(maxEMIterations);
  } else {
    gm.cliqueChainEM(maxEMIterations, 
		     pruneRatio,
		     writeParametersAfterEachEMIteration,
		     outputTrainableParameters,
		     binOutputTrainableParameters,
		     outputMasterFile,
		     loadAccFile,
		     loadAccRange,
		     storeAccFile,
		     accFileIsBinary,
		     llStoreFile,
		     lldp);
  }

  exit_program_with_status(0);
}
