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
#include "debug.h"
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
#include "GMTK_GMTemplate.h"
#include "GMTK_GMParms.h"
#include "GMTK_ObservationMatrix.h"
#include "GMTK_MixtureCommon.h"
#include "GMTK_GaussianComponent.h"
#include "GMTK_MeanVector.h"
#include "GMTK_DiagCovarVector.h"
#include "GMTK_DlinkMatrix.h"
#include "GMTK_ProgramDefaultParms.h"
#include "GMTK_Signals.h"
#include "GMTK_BoundaryTriangulate.h"
#include "GMTK_JunctionTree.h"

/*
 * command line arguments
 */
static bool seedme = false;
static char *strFileName=NULL;
static double varFloor = GMTK_DEFAULT_VARIANCE_FLOOR;
static int bct=GMTK_DEFAULT_BASECASETHRESHOLD;
static int ns=GMTK_DEFAULT_NUM_SPLITS;
static int showFrCliques = 0;
static int jut = -1;
static char* anyTimeTriangulate = NULL;
static bool reTriangulate = false;
static bool rePartition = false;
static unsigned maxNumChunksInBoundary = 1; 
static unsigned chunkSkip = 1; 
static int allocateDenseCpts=0;
static char *cppCommandOptions = NULL;
static char* triangulationHeuristic="completed";
static char* boundaryHeuristic="SFW";
static bool findBestBoundary = true;
static double traverseFraction = 1.0;
static bool loadParameters = false;
static bool jtWeight = true;

static char *inputMasterFile=NULL;
// static char *outputMasterFile=NULL;
static char *inputTrainableParameters=NULL;
static char *inputTriangulatedFile=NULL;
static char *outputTriangulatedFile=NULL;
static bool binInputTrainableParameters=false;
// static char *objsToNotTrainFile=NULL;

static bool noReTriP = false;
static bool noReTriC = false;
static bool noReTriE = false;
static bool continueTriangulating = false;

static bool noBoundaryMemoize = false;
static unsigned numBackupFiles = 10;
static char* forceLeftRight="";
static unsigned verbosity = IM::Default;
static bool printResults = false;
// uncomment when reading in for sparse CPTs
// static char *inputMasterFile=NULL;
// static char *inputTrainableParameters=NULL;
// static bool binInputTrainableParameters=false;

Arg Arg::Args[] = {

  /////////////////////////////////////////////////////////////
  // input parameter/structure file handling

  Arg("cppCommandOptions",Arg::Opt,cppCommandOptions,"Additional CPP command line"),
  Arg("strFile",Arg::Req,strFileName,"Graphical Model Structure File"),
  Arg("inputMasterFile",Arg::Opt,inputMasterFile,"Input file of multi-level master CPP processed GM input parameters"),
  Arg("inputTrainableParameters",Arg::Opt,inputTrainableParameters,"File of only and all trainable parameters"),
  Arg("inputTriangulatedFile",Arg::Opt,inputTriangulatedFile,"Non-default previous triangulated file to start with"),
  Arg("outputTriangulatedFile",Arg::Opt,outputTriangulatedFile,"File name to write resulting triangulation to"),


  /////////////////////////////////////////////////////////////
  // Triangulation Options
  Arg("triangulationHeuristic",
      Arg::Opt,triangulationHeuristic,
      "Triang. heuristic, >1 of S=size,T=time,F=fill,W=wght,X=rev-time,P=pos,H=hint,R=rnd,N=wght-w/o-det"),

  Arg("jtWeight",
      Arg::Opt,jtWeight,
      "True means use an estimate of the JT weight to score triangulation rather than sum of weight"),

  Arg("jtwUB",
      Arg::Opt,JunctionTree::jtWeightUpperBound,
      "True means jtWeight is allways an upper bound on true JT weight, false means jtWeight is estimate"),

  Arg("pfCobWeight",
      Arg::Opt,MaxClique::continuousObservationPerFeaturePenalty,
      "Per-Feature Dimension Continuous Observation Log penalty to use in clique weight calc"),

  Arg("findBestBoundary",
      Arg::Opt,findBestBoundary,
      "Run the (exponential time) boundary algorithm or not."),

  Arg("traverseFraction",
      Arg::Opt,traverseFraction,
      "Fraction of current interface to traverse in boundary recursion."),

  Arg("noBoundaryMemoize",
      Arg::Opt,noBoundaryMemoize,
      "Do not memoize boundaries (less memory but runs slower)"),

  Arg("forceLeftRight",
      Arg::Opt,forceLeftRight,
      "Run boundary algorithm only for either left (L) or right (R) interface, rather than both"),

  Arg("boundaryHeuristic",
      Arg::Opt,boundaryHeuristic,
      "Boundary heuristic, >1 of S=size,F=fill,W=wght,N=wght-w/o-det,M=max-clique,C=max-C-clique,A=st-spc,Q=C-st-spc"),

  Arg("M",
      Arg::Opt,maxNumChunksInBoundary,
      "Max number simultaneous chunks in which boundary may simultaneously exist"),

  Arg("S",
      Arg::Opt,chunkSkip,
      "Number of chunks that should exist between boundaries"),

  Arg("unroll",
      Arg::Opt,jut,
      "Unroll graph & triangulate using heuristics. DON'T use P,C,E constrained triangulation."),

  Arg("anyTimeTriangulate",
      Arg::Opt,anyTimeTriangulate,
      "Run the any-time triangulation algorithm for given duration."),

  Arg("rePartition",
      Arg::Opt,rePartition,
      "Re-Run the boundary algorithm even if .str.trifile exists to produce new partition and new triangulation."),

  Arg("reTriangulate",
      Arg::Opt,reTriangulate,
      "Re-Run only triangluation using existing partition given in .trifile."),

  Arg("continueTriangulating",
      Arg::Opt,continueTriangulating,
      "When re-triangulating existing .tri file, continue besting existing triangulations"),

  Arg("noReTriP",
      Arg::Opt,noReTriP,
      "When re-triangulating existing .tri file, don't re-triangulate P, keep old"),
  Arg("noReTriC",
      Arg::Opt,noReTriC,
      "When re-triangulating existing .tri file, don't re-triangulate C, keep old"),
  Arg("noReTriE",
      Arg::Opt,noReTriE,
      "When re-triangulating existing .tri file, don't re-triangulate E, keep old"),


  Arg("numBackupFiles",Arg::Opt,numBackupFiles,"Number of backup .trifiles (_bak0,_bak1,etc.) to keep."),

  Arg("printResults",Arg::Opt,printResults,"Print information about result of triangulation."),

  Arg("loadParameters",Arg::Opt,loadParameters,"Also load in all trainable parameters."),

  // eventually this next option will be removed.
  Arg("showFrCliques",Arg::Opt,showFrCliques,"Show frontier alg. cliques after the network has been unrolled k times and exit."),

  /////////////////////////////////////////////////////////////
  // General Options

  Arg("allocateDenseCpts",Arg::Opt,allocateDenseCpts,"Automatically allocate any undefined CPTs. arg = -1, no read params, arg = 0 noallocate, arg = 1 means use random initial CPT values. arg = 2, use uniform values"),

  Arg("seed",Arg::Opt,seedme,"Seed the random number generator"),
  Arg("verbosity",Arg::Opt,verbosity,"Verbosity (0 <= v <= 100) of informational/debugging msgs"),
  // final one to signal the end of the list
  Arg()

};

#define MYBS(x) ((x)?"T":"F")

/*
 * A routine to create a string that contains all
 * relevant command line options to the current triangulation.
 * This string will be saved in the .trifile as a comment
 * so the user will know how the trifile was generated.
 *
 */
void createCommandLineOptionString(string& res)
{
  char buff[2048];

  res.clear();

  sprintf(buff,"triangulationHeuristic: %s, ",triangulationHeuristic);
  res += buff;
  
  sprintf(buff,"jtWeight: %s, ",MYBS(jtWeight));
  res += buff;

  sprintf(buff,"pfCobWeight: %f, ",MaxClique::continuousObservationPerFeaturePenalty);
  res += buff;
  
  sprintf(buff,"findBestBoundary: %s, ",MYBS(findBestBoundary));
  res += buff;
  
  sprintf(buff,"traverseFraction: %f, ",traverseFraction);
  res += buff;

  sprintf(buff,"noBoundaryMemoize: %s, ",MYBS(noBoundaryMemoize));
  res += buff;

  sprintf(buff,"forceLeftRight: %s, ",MYBS(forceLeftRight));  
  res += buff;

  sprintf(buff,"boundaryHeuristic: %s, ",boundaryHeuristic);
  res += buff;

  if (anyTimeTriangulate != NULL) {
    sprintf(buff,"anyTimeTriangulate: %s, ",anyTimeTriangulate);
    res += buff;
  }
}



/*
 *
 * definition of needed global arguments
 *
 */
RAND rnd(seedme);
GMParms GM_Parms;
ObservationMatrix globalObservationMatrix;

/*
 *
 * backupTriFile:
 *    Make a backup copys by renaming the file triFile since it might
 *    be a mistake to delete it and since these file scan take
 *    a while to generate.
 *
 *  TODO: move this to general.{cc,h}
 *
 */
void
backupTriFile(const string &triFile) 
{
  if (numBackupFiles == 0)
    return;
  for (unsigned bk_num=(numBackupFiles-1);bk_num>0;bk_num--) {
    char buff[1024];
    sprintf(buff,"%d",bk_num-1);
    string curFile =  triFile + "_bak" + buff;
    if (fsize(curFile.c_str()) == 0)
      continue;

    sprintf(buff,"%d",bk_num);
    string backupFile = triFile + "_bak" + buff;
    if (rename(curFile.c_str(),backupFile.c_str()) != 0)
      infoMsg(IM::Warning,"Warning: could not create backup triangulation file %s.\n",
	      backupFile.c_str());
  }
  if (fsize(triFile.c_str()) == 0)
    return;
  string backupFile = triFile + "_bak0";
  if (rename(triFile.c_str(),backupFile.c_str()) != 0)
    infoMsg(IM::Warning,"Warning: could not create backup triangulation file %s.\n",
	    backupFile.c_str());
}



int
main(int argc,char*argv[])
{
  ////////////////////////////////////////////
  // set things up so that if an FP exception
  // occurs such as an "invalid" (NaN), overflow
  // or divide by zero, we actually get a FPE
  ieeeFPsetup();
  set_new_handler(memory_error);
  InstallSignalHandlers();

  ////////////////////////////////////////////
  // parse arguments
  Arg::parse(argc,argv);
  (void) IM::setGlbMsgLevel(verbosity);

  if (chunkSkip < 1)
    error("Argument error: chunk skip parameter S must be >= 1\n");
  if (maxNumChunksInBoundary < 1)
    error("Argument error: max number chunks in boundary parameter M must be >= 1\n");

  if (abs(MaxClique::continuousObservationPerFeaturePenalty) > 1.0) {
    infoMsg(IM::Warning,"###\n### !!!DANGER WILL ROBINSON!! LARGE -pfCobWeight VALUE %f MIGHT CAUSE FLOATING POINT EXCEPTION. SUGGEST REDUCE IT IF FPE OCCURS!! ###\n###\n",MaxClique::continuousObservationPerFeaturePenalty);
  }

  MixtureCommon::checkForValidRatioValues();
  MeanVector::checkForValidValues();
  DiagCovarVector::checkForValidValues();
  DlinkMatrix::checkForValidValues();

  ////////////////////////////////////////////
  // set global variables/change global state from args
  GaussianComponent::setVarianceFloor(varFloor);
  if (seedme)
    rnd.seed();

  // if (chunkSkip > maxNumChunksInBoundary)
  //  error("ERROR: Must have S<=M at this time.\n");

  /////////////////////////////////////////////
  if (loadParameters) {
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
  }
  GM_Parms.finalizeParameters();

  /////////////////////////////
  // read in the structure of the GM, this will
  // die if the file does not exist.
  FileParser fp(strFileName,cppCommandOptions);
  infoMsg(IM::Tiny,"Finished reading in structure\n");

  // parse the file
  fp.parseGraphicalModel();
  // create the rv variable objects
  fp.createRandomVariableGraph();
  // Make sure that there are no directed loops in the graph.
  fp.ensureValidTemplate();

  if (loadParameters) {
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
    // do nothing else if this option is given
    exit(0);
  }

  BoundaryTriangulate triangulator(fp,maxNumChunksInBoundary,chunkSkip,traverseFraction);

  triangulator.noBoundaryMemoize = noBoundaryMemoize;

  TimerClass* timer = NULL;
  timer = new TimerClass;
  // Initialize the timer if anyTimeTriangulate is selected
  if (anyTimeTriangulate != NULL) {
    time_t given_time;
    given_time = timer->parseTimeString( string(anyTimeTriangulate) );
    if (given_time == 0) {
      error("ERROR: Must specify a non-zero amount of time for -anyTimeTriangulate"); 
    }
    infoMsg(IM::Low, "Triangulating for %d seconds\n", (int)given_time);
    timer->Reset(given_time);
  }
  else { 
    timer->DisableTimer();
  }
  triangulator.useTimer(timer);

  if (jut >= 0) {
    // then Just Unroll, Triangulate, and report on quality of triangulation.
    triangulator.unrollAndTriangulate(string(triangulationHeuristic),
				     jut);
  } else {

    GMTemplate gm_template(fp,maxNumChunksInBoundary,chunkSkip);
    string input_tri_file, output_tri_file;

    if (inputTriangulatedFile == NULL) {
      input_tri_file = string(strFileName) + GMTemplate::fileExtension;
    }
    else {
      input_tri_file = string(inputTriangulatedFile);
    }

    if (outputTriangulatedFile == NULL) {
      output_tri_file = string(strFileName) + GMTemplate::fileExtension;
    }
    else {
      output_tri_file = string(outputTriangulatedFile);
    }

    if (rePartition && !reTriangulate) {
      infoMsg(IM::Warning,"Warning: rePartition=T option forces -reTriangulate option to be true.\n");
      reTriangulate = true;
    }

    // first check if tri_file exists
    if (rePartition || fsize(input_tri_file.c_str()) == 0) {
      // Then do everything (both partition & triangulation)

      // run partition given options
      triangulator.findPartitions(string(boundaryHeuristic),
				  string(forceLeftRight),
				  string(triangulationHeuristic),
				  findBestBoundary,
				  gm_template);
      if (anyTimeTriangulate == NULL) {
	// just run simple triangulation.
	triangulator.triangulate(string(triangulationHeuristic),
				 jtWeight,
				 gm_template);
      } else {
	// run the anytime algorithm.

	// In this case, here we only run triangulation on the
	// provided new P,C,E partitions until the given amount of time
	// has passed, and we save the best triangulation. This 
	// is seen as a separate step, and would be expected
	// to run for a long while.

	triangulator.anyTimeTriangulate(gm_template,jtWeight);
      }

      backupTriFile(output_tri_file);
      oDataStreamFile os(output_tri_file.c_str());
      fp.writeGMId(os);
      string clStr;
      createCommandLineOptionString(clStr);
      gm_template.writePartitions(os,clStr);
      gm_template.writeMaxCliques(os);
      triangulator.ensurePartitionsAreChordal(gm_template);

    } else if (reTriangulate && !rePartition) {

      // first get the id and partition information.
      {
	iDataStreamFile is(input_tri_file.c_str(),false,false);
	if (!fp.readAndVerifyGMId(is))
	  error("ERROR: triangulation file '%s' does not match graph given in structure file '%s'\n",input_tri_file.c_str(),strFileName);
	gm_template.readPartitions(is);
	gm_template.readMaxCliques(is);
	// read the max cliques but don't triangulate with
	// them here.
      }

      if (!continueTriangulating) {
	// If we are *not* continuing on with old triangulation, then
	// we need to check if the user wants us to save *some* of the
	// previous partition triangulations.
	if (!noReTriP)
	  gm_template.clear_P_Cliques();
	if (!noReTriC)
	  gm_template.clear_C_Cliques();
	if (!noReTriE)
	  gm_template.clear_E_Cliques();
      }
	

      // now using the partition triangulate
      if (anyTimeTriangulate == NULL) {
	// just run simple triangulation.
	triangulator.triangulate(string(triangulationHeuristic),
				 jtWeight,
				 gm_template,
				 !noReTriP,!noReTriC,!noReTriE);
      } else {
	// In this case, here we only run triangulation on the
	// provided new P,C,E partitions until the given amount of time
	// has passed, and we save the best triangulation. This 
	// is seen as a separate step, and would be expected
	// to run for a long while.

	triangulator.anyTimeTriangulate(gm_template,jtWeight,
					!noReTriP,!noReTriC,!noReTriE);
      }

      // write everything out anew
      backupTriFile(output_tri_file);
      oDataStreamFile os(output_tri_file.c_str());

      fp.writeGMId(os);
      string clStr;
      createCommandLineOptionString(clStr);
      gm_template.writePartitions(os,clStr);
      gm_template.writeMaxCliques(os);
      triangulator.ensurePartitionsAreChordal(gm_template);

    } else {

      // 
      // Utilize both the partition information and elimination order
      // information already computed and contained in the file. This
      // enables the program to use external triangulation programs,
      // where this program ensures that the result is triangulated
      // and where it reports the quality of the triangulation.

      iDataStreamFile is(input_tri_file.c_str(),false,false);
      if (!fp.readAndVerifyGMId(is))
	error("ERROR: triangulation file '%s' does not match graph given in structure file '%s'\n",input_tri_file.c_str(),strFileName);

      gm_template.readPartitions(is);
      gm_template.readMaxCliques(is);
      gm_template.triangulatePartitionsByCliqueCompletion();
      triangulator.ensurePartitionsAreChordal(gm_template);
    }

    // At this point, one way or another, we've got the fully triangulated graph in data structures.
    // We now just print this information and the triangulation quality out, and then exit.

    if (printResults) {
	printf("\n--- Printing final clique set and clique weights---\n");

	double p_maxWeight = -1.0;
	double p_totalWeight = -1.0; // starting flag
	printf("  --- Prologue summary, %d cliques\n",gm_template.P.cliques.size());
	for (unsigned i=0;i<gm_template.P.cliques.size();i++) {
	  double curWeight = MaxClique::computeWeight(gm_template.P.cliques[i].nodes);
	  printf("   --- P curWeight = %f\n",curWeight);
	  if (curWeight > p_maxWeight) p_maxWeight = curWeight;
	  if (p_totalWeight == -1.0)
	    p_totalWeight = curWeight;
	  else
	    p_totalWeight = p_totalWeight + log10(1+pow(10,curWeight-p_totalWeight));
	}
	printf("  --- Prologue max clique weight = %f, total weight = %f, jt_weight = %f\n",
	       p_maxWeight,p_totalWeight,
	       JunctionTree::junctionTreeWeight(gm_template.P.cliques,gm_template.PCInterface_in_P));

	double c_maxWeight = -1.0;
	double c_totalWeight = -1.0; // starting flag
	printf("  --- Chunk summary, %d cliques\n",gm_template.C.cliques.size());
	for (unsigned i=0;i<gm_template.C.cliques.size();i++) {
	  double curWeight = MaxClique::computeWeight(gm_template.C.cliques[i].nodes);
	  printf("   --- C curWeight = %f\n",curWeight);
	  if (curWeight > c_maxWeight) c_maxWeight = curWeight;
	  if (c_totalWeight == -1.0)
	    c_totalWeight = curWeight;
	  else
	    c_totalWeight = c_totalWeight + log10(1+pow(10,curWeight-c_totalWeight));
	}
	printf("  --- Chunk max clique weight = %f, total Cx%d weight = %f, per-chunk total C weight = %f, jt_weight = %f\n",
	       c_maxWeight,
	       chunkSkip,
	       c_totalWeight,
	       c_totalWeight - log10((double)chunkSkip),
	       JunctionTree::junctionTreeWeight(gm_template.C.cliques,gm_template.CEInterface_in_C));

	double e_maxWeight = -1.0;
	double e_totalWeight = -1.0; // starting flag
	printf("  --- Epilogue summary, %d cliques\n",gm_template.E.cliques.size());
	for (unsigned i=0;i<gm_template.E.cliques.size();i++) {
	  double curWeight = MaxClique::computeWeight(gm_template.E.cliques[i].nodes);
	  printf("   --- E curWeight = %f\n",curWeight);
	  if (curWeight > e_maxWeight) e_maxWeight = curWeight;
	  if (e_totalWeight == -1.0)
	    e_totalWeight = curWeight;
	  else
	    e_totalWeight = e_totalWeight + log10(1+pow(10,curWeight-e_totalWeight));
	}
	const set <RandomVariable*> emptySet;
	printf("  --- Epilogue max clique weight = %f, total weight = %f, jt_weight = %f\n",
	       e_maxWeight,e_totalWeight,
	       JunctionTree::junctionTreeWeight(gm_template.E.cliques,emptySet));

	double maxWeight
	  = (p_maxWeight>c_maxWeight?p_maxWeight:c_maxWeight);
	maxWeight =
	  (maxWeight>e_maxWeight?maxWeight:e_maxWeight);
	double totalWeight = p_totalWeight;
	// log version of: totalWeight += c_totalWeight
	totalWeight += log10(1+pow(10,c_totalWeight-totalWeight));
	// log version of: totalWeight += e_totalWeight
	totalWeight += log10(1+pow(10,e_totalWeight-totalWeight));

	printf("--- Final set (P,Cx%d,E) has max clique weight = %f, total state space = %f ---\n",
	       chunkSkip,
	       maxWeight,
	       totalWeight);

	// print out a couple of total state spaces for various unrollings
	printf("--- Total weight when unrolling %dx = %f ---\n",maxNumChunksInBoundary+chunkSkip-1,totalWeight);

	totalWeight += log10(1+pow(10,c_totalWeight-totalWeight));	
	printf("--- Total weight when unrolling %dx = %f ---\n",maxNumChunksInBoundary+2*chunkSkip-1,totalWeight);

	totalWeight += log10(1+pow(10,log10(3.0) + c_totalWeight-totalWeight));
	printf("--- Total weight when unrolling %dx = %f ---\n",maxNumChunksInBoundary+5*chunkSkip-1,totalWeight);

	totalWeight += log10(1+pow(10,log10(5.0) + c_totalWeight-totalWeight));
	printf("--- Total weight when unrolling %dx = %f ---\n",maxNumChunksInBoundary+10*chunkSkip-1,totalWeight);

	totalWeight += log10(1+pow(10,log10(10.0) + c_totalWeight-totalWeight));
	printf("--- Total weight when unrolling %dx = %f ---\n",maxNumChunksInBoundary+20*chunkSkip-1,totalWeight);

	totalWeight += log10(1+pow(10,log10(30.0) + c_totalWeight-totalWeight));
	printf("--- Total weight when unrolling %dx = %f ---\n",maxNumChunksInBoundary+50*chunkSkip-1,totalWeight);

	totalWeight += log10(1+pow(10,log10(50.0) + c_totalWeight-totalWeight));
	printf("--- Total weight when unrolling %dx = %f ---\n",maxNumChunksInBoundary+100*chunkSkip-1,totalWeight);

	totalWeight += log10(1+pow(10,log10(400.0) + c_totalWeight-totalWeight));
	printf("--- Total weight when unrolling %dx = %f ---\n",maxNumChunksInBoundary+500*chunkSkip-1,totalWeight);

	totalWeight += log10(1+pow(10,log10(500.0) + c_totalWeight-totalWeight));
	printf("--- Total weight when unrolling %dx = %f ---\n",maxNumChunksInBoundary+1000*chunkSkip-1,totalWeight);

      }
    delete timer;
  }

  exit_program_with_status(0);
}
