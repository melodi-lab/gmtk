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
#include "GMTK_DlinkMatrix.h"
#include "GMTK_ProgramDefaultParms.h"

/*
 * command line arguments
 */
bool seedme = false;
float beam=-LZERO;

char *strFileName=NULL;

char *prmOutFile="outParms%d.gmp";
bool binPrmOutFile=false;
bool writeParametersAfterEachEMIteration=true;

char *prmMasterFile=NULL;
char *prmTrainableFile=NULL;
bool binPrmTrainableFile=false;
char *objsToNotTrainFile=NULL;

unsigned maxEMIterations=3;
bool randomizeParams = true;
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

char *argsFile = NULL;
char *cppCommandOptions = NULL;

int bct=GMTK_DEFAULT_BASECASETHRESHOLD;
int ns=GMTK_DEFAULT_NUM_SPLITS;

void makeArgs(Argument_List &args)
{
  bool required=1,optional=0;

  /////////////////////////////////////////////////////////////
  // observation input file handling

  args.add("of1",required,&ofs[0],
           "Observation File 1");
  args.add("nf1",optional,&nfs[0],
           "Number of floats in observation file 1");
  args.add("ni1",optional,&nis[0],
           "Number of ints in observation file 1");
  args.add("fr1",optional,&frs[0],
           "Float range for observation file 1");
  args.add("ir1",optional,&irs[0],
           "Int range for observation file 1");
  args.add("fmt1",optional,&fmts[0],
           "Format (htk,bin,asc,pfile) for observation file 1");
  args.add("iswp1",optional,&iswps[0],
           "Endian swap condition for observation file 1");

  args.add("of2",optional,&ofs[1],
           "Observation File 2");
  args.add("nf2",optional,&nfs[1],
           "Number of floats in observation file 2");
  args.add("ni2",optional,&nis[1],
           "Number of ints in observation file 2");
  args.add("fr2",optional,&frs[1],
           "Float range for observation file 2");
  args.add("ir2",optional,&irs[1],
           "Int range for observation file 2");
  args.add("fmt2",optional,&fmts[1],
           "Format (htk,bin,asc,pfile) for observation file 2");
  args.add("iswp2",optional,&iswps[1],
           "Endian swap condition for observation file 2");

  args.add("of3",optional,&ofs[2],
           "Observation File 3");
  args.add("nf3",optional,&nfs[2],
           "Number of floats in observation file 3");
  args.add("ni3",optional,&nis[2],
           "Number of ints in observation file 3");
  args.add("fr3",optional,&frs[2],
           "Float range for observation file 3");
  args.add("ir3",optional,&irs[2],
           "Int range for observation file 3");
  args.add("fmt3",optional,&fmts[2],
           "Format (htk,bin,asc,pfile) for observation file 3");
  args.add("iswp3",optional,&iswps[2],
           "Endian swap condition for observation file 3");

  /////////////////////////////////////////////////////////////
  // input parameter/structure file handling

  args.add("strFile",required,&strFileName,
           "GM Structure File");
  args.add("prmMasterFile",required,&prmMasterFile,
           "Multi-level master CPP processed GM Parms File");
  args.add("prmTrainableFile",optional,&prmTrainableFile,
           "File containing Trainable Parameters");
  args.add("binPrmTrainableFile",optional,&binPrmTrainableFile,
           "Is Binary? File containing Trainable Parameters");
  args.add("cppCommandOptions",optional,&cppCommandOptions,
           "Command line options to give to cpp");

  args.add("objsToNotTrainFile",optional,&objsToNotTrainFile,
           "File list list trainable parm objs not train.");


  args.add("prmOutFile",optional,&prmOutFile,
           "File to place *TRAINABLE* output parametes");
  args.add("binPrmOutFile",optional,&binPrmOutFile,
           "Output parametes binary? (def=false)");


  args.add("wpaeei",optional,&writeParametersAfterEachEMIteration,
           "Write Parameters *After* Each EM Iteration? (def=true)");

  /////////////////////////////////////////////////////////////
  // general files

  args.add("seed",optional,&seedme,
           "Seed the RN generator");
  args.add("maxEmIters",optional,&maxEMIterations,
           "Max number of EM iterations to do");
  args.add("beam",optional,&beam,
           "Beam width (less than max*exp(-beam) are pruned away)");

  // support for splitting and vanishing
  args.add("mcvr",optional,&MixGaussiansCommon::mixCoeffVanishRatio,
           "Mixture Coefficient Vanishing Ratio");
  args.add("botForceVanish",optional,
           &MixGaussiansCommon::numBottomToForceVanish,
           "Number of bottom mixture components to force vanish");
  
  args.add("mcsr",optional,&MixGaussiansCommon::mixCoeffSplitRatio,
           "Mixture Coefficient Splitting Ratio");
  args.add("topForceSplit",optional,
           &MixGaussiansCommon::numTopToForceSplit,
           "Number of top mixture components to force split");

  args.add("meanCloneSTDfrac",optional,&MeanVector::cloneSTDfrac,
           "Fraction of mean to use for STD in mean clone");
  args.add("covarCloneSTDfrac",optional,
           &DiagCovarVector::cloneSTDfrac,
           "Fraction of var to use for STD in covar clone");
  args.add("dlinkCloneSTDfrac",optional,&DlinkMatrix::cloneSTDfrac,
           "Fraction of var to use for STD in covar clone");

  args.add("cloneShareMeans",optional,
           &GaussianComponent::cloneShareMeans,
           "Gaussian component clone shares parent mean");
  args.add("cloneShareCovars",optional,
           &GaussianComponent::cloneShareCovars,
           "Gaussian component clone shares parent covars");
  args.add("cloneShareDlinks",optional,
           &GaussianComponent::cloneShareDlinks,
           "Gaussian component clone shares parent dlinks");


  args.add("varFloor",optional,&varFloor,
           "Variance Floor");
  args.add("floorVarOnRead",optional,
           &DiagCovarVector::floorVariancesWhenReadIn,
           "Floor the variances to varFloor when they are read in");



  args.add("lldp",optional,&lldp,
           "Log Likelihood difference percentage for termination");
  args.add("mnlldp",optional,&mnlldp,
           "Absolute value of max negative Log Likelihood difference percentage for termination");

  args.add("trrng",optional,&trrng_str,
           "Range to train over segment file");

  args.add("storeAccFile",optional,&storeAccFile,
           "Store accumulators file");
  args.add("loadAccFile",optional,&loadAccFile,
           "Load accumulators file"); 
  args.add("loadAccRange",optional,&loadAccRange,
           "Load accumulators file range"); 
  args.add("llStoreFile",optional,&llStoreFile,
           "File to store previous sum LL's"); 
  args.add("accFileIsBinary",optional,&accFileIsBinary,
           "Binary accumulator files (def true)"); 

  args.add("startSkip",optional,&startSkip,
           "Frames to skip at beginning (i.e., first frame is buff[startSkip])");
  args.add("endSkip",optional,&endSkip,
           "Frames to skip at end (i.e., last frame is buff[len-1-endSkip])");
  
  args.add("cptNormThreshold",optional,&CPT::normalizationThreshold,
           "Read error if |Sum-1.0|/card > norm_threshold");
  args.add("random",optional,&randomizeParams,
           "Randomize the parameters");
  args.add("enem",optional,&enem,"Run enumerative EM");

  args.add("showCliques",optional,&show_cliques,
           "Show the cliques after the netwok has been unrolled k times.");

  args.add("argsFile",optional,&argsFile,
           "File to get args from (overrides specified comand line args).");

  args.add("numSplits",optional,&ns,
           "Number of splits to use in logspace recursion (>=2).");

  args.add("baseCaseThreshold",optional,&bct,
           "Base case threshold to end recursion (>=2).");

  args.add("gaussianCache",optional,
           &MixGaussiansCommon::cacheGaussiansInEmTraining,
           "Cache Gaussians evaluations during EM training. true will speeds things up, but uses more memory.");
}

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
  Argument_List args;
  makeArgs(args);
  args.parse(argc, argv);

  ////////////////////////////////////////////
  // check for valid argument values.
  int nfiles = 0;
  unsigned ifmts[MAX_NUM_OBS_FILES];
  for (int i=0;i<MAX_NUM_OBS_FILES;i++) {
    if (ofs[i] != NULL && nfs[i] == 0 && nis[i] == 0)
      error("ERROR: command line must specify one of nf%d and ni%d not zero",
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

  MixGaussiansCommon::checkForValidRatioValues();
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
  if (prmMasterFile) {
    // flat, where everything is contained in one file, always ASCII
    iDataStreamFile pf(prmMasterFile,false,true,cppCommandOptions);
    GM_Parms.read(pf);
  }
  if (prmTrainableFile) {
    // flat, where everything is contained in one file
    iDataStreamFile pf(prmTrainableFile,binPrmTrainableFile,true,cppCommandOptions);
    GM_Parms.readTrainable(pf);
  }
  GM_Parms.markObjectsToNotTrain(objsToNotTrainFile,cppCommandOptions);

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
  fp.associateWithDataParams();
  // make sure that all observation variables work
  // with the global observation stream.
  fp.checkConsistentWithGlobalObservationStream();
  // now associate the RVs with a GM
  GMTK_GM gm;
  fp.addVariablesToGM(gm);
  gm.verifyTopologicalOrder();

  ////////////////////////////////////
  // set up the observation stream
  gm.setExampleStream(obsFileName,trrng_str);
  GM_Parms.checkConsistentWithGlobalObservationStream();

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
    printf("NOTE: Randomizing initial parameters\n");
    GM_Parms.makeRandom();
    printf("NOTE: Writing iniial randomized parameters to random.gmp\n");
    oDataStreamFile of("random.gmp");
    GM_Parms.writeAll(of);
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
		     prmOutFile,
		     binPrmOutFile,
		     loadAccFile,
		     loadAccRange,
		     storeAccFile,
		     accFileIsBinary,
		     llStoreFile,
		     lldp);
  }

  exit_program_with_status(0);
}
