/*
 * gmtkUnroll.cc
 * Unroll a graph
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
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
#include "GMTK_GMTemplate.h"
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

unsigned num_times=5;


char *cppCommandOptions = NULL;

int bct=GMTK_DEFAULT_BASECASETHRESHOLD;
int ns=GMTK_DEFAULT_NUM_SPLITS;

Arg Arg::Args[] = {


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

  Arg("varFloor",Arg::Opt,varFloor,"Variance Floor"),
  Arg("floorVarOnRead",Arg::Opt,DiagCovarVector::floorVariancesWhenReadIn,
       "Floor the variances to varFloor when they are read in"),

  Arg("numTimes",Arg::Opt,num_times,"Number of times to unroll"),

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
  GM_Parms.markObjectsToNotTrain(objsToNotTrainFile,cppCommandOptions);

  /////////////////////////////
  // read in the structure of the GM, this will
  // die if the file does not exist.
  FileParser fp(strFileName,cppCommandOptions);

  printf("Finished reading in all parameters and structures\n");

  // parse the file
  fp.parseGraphicalModel();
  // create the rv variable objects
  fp.createRandomVariableGraph();
  fp.ensureValidTemplate();
  // link the RVs with the parameters that are contained in
  // the bn1_gm.dt file.
  fp.associateWithDataParams();
  // make sure that all observation variables work
  // with the global observation stream.
  // fp.checkConsistentWithGlobalObservationStream();

  // GMTemplate gm_template;
  // fp.addVariablesToTemplate(gm_template);
  // gm_template.print();
  vector <RandomVariable*> vars;

  fp.unroll(num_times,vars);

  for (unsigned i=0;i<vars.size();i++) {
    printf("frame %d, variable %s\n",
	   vars[i]->timeIndex,
	   vars[i]->name().c_str());
  }


  exit_program_with_status(0);
}
