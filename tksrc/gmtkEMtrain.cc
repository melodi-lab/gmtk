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


/*
 * command line arguments
 */
bool seedme = false;
float beam=1000;

char *strFileName=NULL;

char *prmOutFile="outParms%d.gmp";
bool binPrmOutFile=false;
bool writeParametersAfterEachEMIteration=true;

char *prmMasterFile=NULL;
char *prmTrainableFile=NULL;
bool binPrmTrainableFile=false;

unsigned maxEMIterations=3;
bool randomizeParams = true;
bool enem = false;
double mcvr = 1e20;
double mcsr = 1e10;
double varFloor = 1e-10;
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


ARGS ARGS::Args[] = {

  /////////////////////////////////////////////////////////////
  // observation input file handling

  ARGS("of1",ARGS::Req,ofs[0],"Observation File 1"),
  ARGS("nf1",ARGS::Opt,nfs[0],"Number of floats in observation file 1"),
  ARGS("ni1",ARGS::Opt,nis[0],"Number of ints in observation file 1"),
  ARGS("fr1",ARGS::Opt,frs[0],"Float range for observation file 1"),
  ARGS("ir1",ARGS::Opt,irs[0],"Int range for observation file 1"),
  ARGS("fmt1",ARGS::Opt,fmts[0],"Format (htk,bin,asc,pfile) for observation file 1"),
  ARGS("iswp1",ARGS::Opt,iswps[0],"Endian swap condition for observation file 1"),


  ARGS("of2",ARGS::Opt,ofs[1],"Observation File 2"),
  ARGS("nf2",ARGS::Opt,nfs[1],"Number of floats in observation file 2"),
  ARGS("ni2",ARGS::Opt,nis[1],"Number of ints in observation file 2"),
  ARGS("fr2",ARGS::Opt,frs[1],"Float range for observation file 2"),
  ARGS("ir2",ARGS::Opt,irs[1],"Int range for observation file 2"),
  ARGS("fmt2",ARGS::Opt,fmts[1],"Format (htk,bin,asc,pfile) for observation file 2"),
  ARGS("iswp2",ARGS::Opt,iswps[1],"Endian swap condition for observation file 2"),


  ARGS("of3",ARGS::Opt,ofs[2],"Observation File 3"),
  ARGS("nf3",ARGS::Opt,nfs[2],"Number of floats in observation file 3"),
  ARGS("ni3",ARGS::Opt,nis[2],"Number of ints in observation file 3"),
  ARGS("fr3",ARGS::Opt,frs[2],"Float range for observation file 3"),
  ARGS("ir3",ARGS::Opt,irs[2],"Int range for observation file 3"),
  ARGS("fmt3",ARGS::Opt,fmts[2],"Format (htk,bin,asc,pfile) for observation file 3"),
  ARGS("iswp3",ARGS::Opt,iswps[2],"Endian swap condition for observation file 3"),


  /////////////////////////////////////////////////////////////
  // input parameter/structure file handling

  ARGS("prmMasterFile",ARGS::Req,prmMasterFile,"Multi-level master CPP processed GM Parms File"),

  ARGS("prmTrainableFile",ARGS::Opt,prmTrainableFile,"File containing Trainable Parameters"),
  ARGS("binPrmTrainableFile",ARGS::Opt,binPrmTrainableFile,"Is Binary? File containing Trainable Parameters"),


  ARGS("prmOutFile",ARGS::Opt,prmOutFile,"File to place *TRAINABLE* output parametes"),
  ARGS("binPrmOutFile",ARGS::Opt,binPrmOutFile,"Output parametes binary? (def=false)"),


  ARGS("wpaeei",ARGS::Opt,writeParametersAfterEachEMIteration,
      "Write Parameters *After* Each EM Iteration? (def=true)"),

  ARGS("strFile",ARGS::Req,strFileName,"GM Structure File"),

  /////////////////////////////////////////////////////////////
  // general files

  ARGS("seed",ARGS::Opt,seedme,"Seed the RN generator"),
  ARGS("maxEmIters",ARGS::Opt,maxEMIterations,"Max number of EM iterations to do"),
  ARGS("beam",ARGS::Opt,beam,"Beam, values less than this*max are pruned"),

  // support for splitting and vanishing
  ARGS("mcvr",ARGS::Opt,MixGaussiansCommon::mixCoeffVanishRatio,"Mixture Coefficient Vanishing Ratio"),
  ARGS("mcsr",ARGS::Opt,MixGaussiansCommon::mixCoeffSplitRatio,"Mixture Coefficient Splitting Ratio"),

  ARGS("meanCloneSTDfrac",ARGS::Opt,MeanVector::cloneSTDfrac,"Fraction of mean to use for STD in mean clone"),
  ARGS("covarCloneSTDfrac",ARGS::Opt,DiagCovarVector::cloneSTDfrac,"Fraction of var to use for STD in covar clone"),
  ARGS("dlinkCloneSTDfrac",ARGS::Opt,DlinkMatrix::cloneSTDfrac,"Fraction of var to use for STD in covar clone"),

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

    globalObservationMatrix.openFiles(nfiles,
				      (const char**)&ofs,
				      (const char**)&frs,
				      (const char**)&irs,
				      (unsigned*)&nfs,
				      (unsigned*)&nis,
				      (unsigned*)&ifmts,
				      (bool*)&iswps);
  }

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
    iDataStreamFile pf(prmMasterFile,false);
    GM_Parms.read(pf);
  }
  if (prmTrainableFile) {
    // flat, where everything is contained in one file, always ASCII
    iDataStreamFile pf(prmTrainableFile,false);
    GM_Parms.readTrainable(pf);
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
  gm.setExampleStream(obsFileName,trrng_str);


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
