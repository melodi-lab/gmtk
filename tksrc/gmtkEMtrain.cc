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
char *outFileName="outParms%d.gmp";
bool writeParametersAfterEachEMIteration=true;
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

#define MAX_NUM_FILES (3)
char *of[MAX_NUM_FILES] = { NULL, NULL, NULL }; 
int nf[MAX_NUM_FILES] = { 0, 0, 0 };
int ni[MAX_NUM_FILES] = { 0, 0, 0 };
char *fr[MAX_NUM_FILES] = { "all", "all", "all" };
char *ir[MAX_NUM_FILES] = { "all", "all", "all" };
char *fmt[MAX_NUM_FILES] = { "pfile", "pfile", "pfile" };
bool iswp[MAX_NUM_FILES] = { false, false, false };


char *of1=NULL;        // observation file 1
int nf1=0;             // number of floats in observation file 1
int ni1=0;             // number of ints in observation file 1
char *fr1="all";       // float range 1
char *ir1="all";       // int range 1
char *fmt1="pfile";    // format for file 1
bool iswp1=false;      // endian swap flag for file 1

char *of2=NULL;        // observation file 2
int nf2=0;             // number of floats in observation file 2
int ni2=0;             // number of ints in observation file 2
char *fr2="all";       // float range 2
char *ir2="all";       // int range 2
char *fmt2="pfile";    // format for file 2
bool iswp2=false;      // endian swap flag for file 2


char *of3=NULL;        // observation file 3
int nf3=0;             // number of floats in observation file 3
int ni3=0;             // number of ints in observation file 3
char *fr3="all";       // float range 3
char *ir3="all";       // int range 3
char *fmt3="pfile";    // format for file 3
bool iswp3=false;      // endian swap flag for file 3


ARGS ARGS::Args[] = {

 // observation file handling

 ARGS("of1",ARGS::Req,of[0],"Observation File 1"),
 ARGS("nf1",ARGS::Opt,nf[0],"Number of floats in observation file 1"),
 ARGS("ni1",ARGS::Opt,nf[0],"Number of ints in observation file 1"),
 ARGS("fr1",ARGS::Opt,fr[0],"Float range for observation file 1"),
 ARGS("ir1",ARGS::Opt,ir[0],"Int range for observation file 1"),
 ARGS("fmt1",ARGS::Opt,fmt[0],"Format (htk,bin,asc,pfile) for observation file 1"),
 ARGS("iswp1",ARGS::Opt,iswp[0],"Endian swap condition for observation file 1"),


 ARGS("of2",ARGS::Opt,of[1],"Observation File 2"),
 ARGS("nf2",ARGS::Opt,nf[1],"Number of floats in observation file 2"),
 ARGS("ni2",ARGS::Opt,nf[1],"Number of ints in observation file 2"),
 ARGS("fr2",ARGS::Opt,fr[1],"Float range for observation file 2"),
 ARGS("ir2",ARGS::Opt,ir[1],"Int range for observation file 2"),
 ARGS("fmt2",ARGS::Opt,fmt[1],"Format (htk,bin,asc,pfile) for observation file 2"),
 ARGS("iswp2",ARGS::Opt,iswp[1],"Endian swap condition for observation file 2"),


 ARGS("of3",ARGS::Opt,of[2],"Observation File 3"),
 ARGS("nf3",ARGS::Req,nf[2],"Number of floats in observation file 3"),
 ARGS("ni3",ARGS::Opt,nf[2],"Number of ints in observation file 3"),
 ARGS("fr3",ARGS::Opt,fr[2],"Float range for observation file 3"),
 ARGS("ir3",ARGS::Opt,ir[2],"Int range for observation file 3"),
 ARGS("fmt3",ARGS::Opt,fmt[2],"Format (htk,bin,asc,pfile) for observation file 3"),
 ARGS("iswp3",ARGS::Opt,iswp[2],"Endian swap condition for observation file 3"),


 ARGS("obsFile",ARGS::Req,obsFileName,"File containing observations"),

 ARGS("parmsPtrFile",ARGS::Opt,parmsPtrFileName,"Multi-level GM Parms File"),

 ARGS("parmsFile",ARGS::Opt,parmsFileName,"Single-level GM Parms File"),
 ARGS("binParmsFile",ARGS::Opt,binParmsFile,"Is Single-level GM Parms File binary? (def=false)"),

 ARGS("outFileName",ARGS::Opt,outFileName,"File to place output parametes"),
 ARGS("binOutFile",ARGS::Opt,binOutFile,"Output parametes binary? (def=false)"),
 ARGS("wpaeei",ARGS::Opt,writeParametersAfterEachEMIteration,
      "Write Parameters *After* Each EM Iteration? (def=true)"),

 ARGS("strFile",ARGS::Req,strFileName,"GM Structure File"),

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
  if (of1 != NULL & nf1 <= 0 && ni2 <= 0)
    error("ERROR: command line must specify one of nf1 and ni1 not zero");
  if (of2 != NULL & nf2 <= 0 && ni2 <= 0)
    error("ERROR: command line must specify one of nf2 and ni2 not zero");
  if (of3 != NULL & nf3 <= 0 && ni3 <= 0)
    error("ERROR: command line must specify one of nf3 and ni3 not zero");
  int nfiles = 0;
  nfiles += (of2 != NULL);
  nfiles += (of3 != NULL);
  {  
    char *ofs[3] = { of1, of2, of3 };
    char *frs[3] = { fr1, fr2, fr3 };
    char *irs[3] = { ir1, ir2, ir3 };
    int  nfs[3] = { nf1, nf2, nf3 };
    int  nis[3] = { ni1, ni2, ni3 };

    int fmts[3];

    int i = 0;
    if (strcmp(fmt1,"htk") == 0)
      fmts[i] = HTK;
    else if (strcmp(fmt1,"binary") == 0)
      fmts[i] = RAWBIN;
    else if (strcmp(fmt1,"ascii") == 0)
      fmts[i] = RAWASC;
    else if (strcmp(fmt1,"pfile") == 0)
      fmts[i] = PFILE;
    else {
      printf("ERROR: Unknown observation file format type: '%s'\n",fmt1);
      exit(-1); 
    }

    bool iswps[] = { iswp1, iswp2, iswp3 };
    globalObservationMatrix.openFiles(nfiles,
				      (char**)&ofs,
				      (char**)&frs,
				      (char**)&irs,
				      (int*)&nfs,
				      (int*)&nis,
				      (int*)fmts,
				      (bool*)iswps);
  }


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
  // @@@@  gm.setExampleStream(obsFileName,trrng_str);


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
		     outFileName,
		     binOutFile,
		     loadAccFile,
		     loadAccRange,
		     storeAccFile,
		     accFileIsBinary,
		     llStoreFile,
		     lldp);
  }

  exit_program_with_status(0);
}
