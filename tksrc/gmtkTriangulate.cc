/*
 * gmtkTriangulate.cc
 * triangulate a graph
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
#include "GMTK_MixGaussiansCommon.h"
#include "GMTK_GaussianComponent.h"
#include "GMTK_MeanVector.h"
#include "GMTK_DiagCovarVector.h"
#include "GMTK_DlinkMatrix.h"
#include "GMTK_ProgramDefaultParms.h"
#include "GMTK_AnyTimeTriangulation.h"
#include "GMTK_GMTemplate.h"


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

int bct=GMTK_DEFAULT_BASECASETHRESHOLD;
int ns=GMTK_DEFAULT_NUM_SPLITS;

int startSkip = 0;
int endSkip = 0;

int showFrCliques = 0;
int jut = 0;
char* anyTimeTriangulate = NULL;
bool reTriangulate = false;
bool rePartition = false;

int allocateDenseCpts=-1;
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

char *cppCommandOptions = NULL;


char* triangulationHeuristic="WFS";
char* faceHeuristic="SFWC";

bool findBestFace = true;

char* forceLeftRight="";

Arg Arg::Args[] = {


  /////////////////////////////////////////////////////////////
  // input parameter/structure file handling

  Arg("cppCommandOptions",Arg::Opt,cppCommandOptions,"Additional CPP command line"),
  Arg("strFile",Arg::Req,strFileName,"Graphical Model Structure File"),

  /////////////////////////////////////////////////////////////
  // Triangulation Options
  Arg("triangulationHeuristic",Arg::Opt,triangulationHeuristic,"Elimination heuristic, one+ of S=size,T=time,F=fill,W=wght,E=entr,P=pos,H=hint"),

  Arg("findBestFace",Arg::Opt,findBestFace,"Run find-best-face (exponential time) algorithm or not."),
  Arg("forceLeftRight",Arg::Opt,forceLeftRight,"Use only either left (L) or right (R) face heuristic, rather than best of both"),
  Arg("faceHeuristic",Arg::Opt,faceHeuristic,"Face heuristic, one+ of S=size,F=fill,W=wght,C=min(max(C-clique)),M=min(max(clique))"),

  Arg("unroll",Arg::Opt,jut,"Unroll graph & triangulate using heuristics. DON'T use P,C,E constrained triangulation."),
  Arg("anyTimeTriangulate",Arg::Opt,anyTimeTriangulate,"Run the any-time triangulation algorithm for given duration."),
  Arg("reTriangulate",Arg::Opt,reTriangulate,"Re-Run triangluation algorithm even if .str.trifile exists with existing partition elimination ordering"),
  Arg("rePartition",Arg::Opt,rePartition,"Re-Run the exponential partition finding algorithm even if .str.trifile exists with existing partition ordering"),

  // eventually this next option will be removed.
  Arg("showFrCliques",Arg::Opt,showFrCliques,"Show frontier alg. cliques after the network has been unrolled k times and exit."),

  /////////////////////////////////////////////////////////////
  // General Options

  Arg("allocateDenseCpts",Arg::Opt,allocateDenseCpts,"Automatically allocate any undefined CPTs. arg = -1, no read params, arg = 0 noallocate, arg = 1 means use random initial CPT values. arg = 2, use uniform values"),

  Arg("seed",Arg::Opt,seedme,"Seed the random number generator"),
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

#if 0
  // don't read in parameters for now
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
  GM_Parms.loadGlobal();
  GM_Parms.markObjectsToNotTrain(objsToNotTrainFile,cppCommandOptions);
#endif

  /////////////////////////////
  // read in the structure of the GM, this will
  // die if the file does not exist.
  FileParser fp(strFileName,cppCommandOptions);

  printf("Finished reading in structure\n");

  // parse the file
  fp.parseGraphicalModel();
  // create the rv variable objects
  fp.createRandomVariableGraph();

  // Make sure that there are no directed loops in the graph.
  {
    vector <RandomVariable*> vars;
    vector <RandomVariable*> vars2;
    // just unroll it one time to make sure it is valid, and
    // then make sure graph has no loops.
    fp.unroll(1,vars);
    if (!GraphicalModel::topologicalSort(vars,vars2))
      // TODO: fix this error message.
      error("ERROR. Graph is not directed and contains has a directed loop.\n");
  }

  
  // link the RVs with the parameters that are contained in
  // the bn1_gm.dt file.
  if (allocateDenseCpts >= 0) {
    if (allocateDenseCpts == 0)
      fp.associateWithDataParams(FileParser::noAllocate);
    else if (allocateDenseCpts == 1)
      fp.associateWithDataParams(FileParser::allocateRandom);
    else if (allocateDenseCpts == 2)
      fp.associateWithDataParams(FileParser::allocateUniform);
    else
      error("Error: command line argument '-allocateDenseCpts d', must have d = {0,1,2}\n");
  }

  // make sure that all observation variables work
  // with the global observation stream.
  // fp.checkConsistentWithGlobalObservationStream();

  if (showFrCliques) {
    // if this option is set, just run the frontier algorithm
    // and then quite. TODO: remove all of this.
    GMTK_GM gm(&fp);
    fp.addVariablesToGM(gm);
    gm.verifyTopologicalOrder();

    gm.setCliqueChainRecursion(ns, bct);

    gm.GM2CliqueChain();
    gm.setupForVariableLengthUnrolling(fp.firstChunkFrame(),fp.lastChunkFrame());
    printf("Frontier cliques in the unrolled network are:\n");
    gm.setSize(showFrCliques);
    gm.showCliques();
    exit(0);
  }


  GMTemplate gm_template(fp);


  if (jut > 0) {
    gm_template.unrollAndTriangulate(string(triangulationHeuristic),
				     jut);
  } else {
    GMInfo gm_info;

    set<RandomVariable*> P;
    set<RandomVariable*> C;
    set<RandomVariable*> E;
    vector<MaxClique> Pcliques;
    vector<MaxClique> Ccliques;
    vector<MaxClique> Ecliques;
    vector<RandomVariable*> Pordered;
    vector<RandomVariable*> Cordered;
    vector<RandomVariable*> Eordered;

    string tri_file = string(strFileName) + ".trifile";
    if (rePartition && !reTriangulate) {
      warning("Warning: rePartition=T option forces reTriangulate option to be true.\n");
      reTriangulate = true;
    }

    // first check of tri_file exists
    if (rePartition || fsize(tri_file.c_str()) == 0) {
      // then do everything (partition & triangulation)

      oDataStreamFile os(tri_file.c_str());
      fp.writeGMId(os);
      // run partition given options
      gm_template.findPartitions(string(faceHeuristic),
				 string(forceLeftRight),
				 string(triangulationHeuristic),
				 findBestFace,
				 gm_info);
      if (anyTimeTriangulate == NULL) {
	// just run simple triangulation.
	gm_template.triangulatePartitions(string(triangulationHeuristic),
					  gm_info);
      } else {
	// run the anytime algorithm.

	// In this case, here we only run triangulation on the
	// provided new P,C,E partitions until the given amount of time
	// has passed, and we save the best triangulation. This 
	// is seen as a separate step, and would be expected
	// to run for a long while.
	// AnyTimeTriangulation anyTime(string(anyTimeTriangulate),gm_template);
	error("Anytime algorithm not yet implemented\n");
      }
      {
	printf("\n--- Printing final clique set and clique weights---\n");
	float maxWeight = -1.0;
	for (unsigned i=0;i<gm_info.Pcliques.size();i++) {
	  float curWeight = gm_template.computeWeight(gm_info.Pcliques[i].nodes);
	  printf("   --- P curWeight = %f\n",curWeight);
	  if (curWeight > maxWeight) maxWeight = curWeight;
	}
	for (unsigned i=0;i<gm_info.Ccliques.size();i++) {
	  float curWeight = gm_template.computeWeight(gm_info.Ccliques[i].nodes);
	  printf("   --- C curWeight = %f\n",curWeight);
	  if (curWeight > maxWeight) maxWeight = curWeight;
	}
	for (unsigned i=0;i<gm_info.Ecliques.size();i++) {
	  float curWeight = gm_template.computeWeight(gm_info.Ecliques[i].nodes);
	  printf("   --- E curWeight = %f\n",curWeight);
	  if (curWeight > maxWeight) maxWeight = curWeight;
	}
	printf("--- Final set has max clique weight = %f ---\n",maxWeight);
      }
      

      gm_template.storePartitions(os,gm_info);
      gm_template.storePartitionTriangulation(os,
					      gm_info);

    } else if (reTriangulate && !rePartition) {
      // utilize the parition information already there but still
      // run re-triangulation
      
      // first get the id and partition information.
      {
	iDataStreamFile is(tri_file.c_str());
	if (!fp.readAndVerifyGMId(is))
	  error("ERROR: triangulation file '%s' does not match graph given in structure file '%s'\n",tri_file.c_str(),strFileName);
	gm_template.findPartitions(is,
				   gm_info);
      }
      // now using the partition triangulate
      if (anyTimeTriangulate == NULL) {
	// just run simple triangulation.
	gm_template.triangulatePartitions(string(triangulationHeuristic),
					  gm_info);
      } else {
	// In this case, here we only run triangulation on the
	// provided new P,C,E partitions until the given amount of time
	// has passed, and we save the best triangulation. This 
	// is seen as a separate step, and would be expected
	// to run for a long while.

	AnyTimeTriangulation attr(gm_template,string(anyTimeTriangulate));
	attr.triangulatePartitions(gm_info);

      }
      // write everything out anew
      oDataStreamFile os(tri_file.c_str());
      fp.writeGMId(os);
      gm_template.storePartitions(os,gm_info);
      gm_template.storePartitionTriangulation(os,
					      gm_info);

    } else {
      // utilize both the partition information and elimination order information
      // already computed and contained in the file

      iDataStreamFile is(tri_file.c_str());
      if (!fp.readAndVerifyGMId(is))
	error("ERROR: triangulation file '%s' does not match graph given in structure file '%s'\n",tri_file.c_str(),strFileName);

      gm_template.findPartitions(is,
				 gm_info);
      gm_template.triangulatePartitions(is,
					gm_info);
    }
  }

  exit_program_with_status(0);
}
