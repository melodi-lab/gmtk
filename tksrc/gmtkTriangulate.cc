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
static bool seedme = false;
static float beam=-LZERO;
static char *strFileName=NULL;
static double varFloor = GMTK_DEFAULT_VARIANCE_FLOOR;
static float lldp = 0.001;
static float mnlldp = 0.01;
static int bct=GMTK_DEFAULT_BASECASETHRESHOLD;
static int ns=GMTK_DEFAULT_NUM_SPLITS;
static int startSkip = 0;
static int endSkip = 0;
static int showFrCliques = 0;
static int jut = -1;
static char* anyTimeTriangulate = NULL;
static bool reTriangulate = false;
static bool rePartition = false;
static unsigned M = 1; 
static int allocateDenseCpts=-1;
static char *cppCommandOptions = NULL;
static char* triangulationHeuristic="WFS";
static char* faceHeuristic="SFW";
static bool findBestFace = true;
static char* forceLeftRight="";
// uncomment when reading in for sparse CPTs
// static char *inputMasterFile=NULL;
// static char *inputTrainableParameters=NULL;
// static bool binInputTrainableParameters=false;

Arg Arg::Args[] = {

  /////////////////////////////////////////////////////////////
  // input parameter/structure file handling

  Arg("cppCommandOptions",Arg::Opt,cppCommandOptions,"Additional CPP command line"),
  Arg("strFile",Arg::Req,strFileName,"Graphical Model Structure File"),

  /////////////////////////////////////////////////////////////
  // Triangulation Options
  Arg("triangulationHeuristic",Arg::Opt,triangulationHeuristic,"Elim heuristic, >1 of S=size,T=time,F=fill,W=wght,E=entr,P=pos,H=hint,R=rnd,N=wght-w/o-det"),

  Arg("findBestFace",Arg::Opt,findBestFace,"Run find-best-face (exponential time) algorithm or not."),
  Arg("forceLeftRight",Arg::Opt,forceLeftRight,"Use only either left (L) or right (R) face heuristic, rather than best of both"),
  Arg("faceHeuristic",Arg::Opt,faceHeuristic,"Face heuristic, >1 of S=size,F=fill,W=wght,N=wght-w/o-det,E=entr,M=max-clique,C=max-C-clique,A=st-spc,Q=C-st-spc"),
  Arg("M",Arg::Opt,M,"Mumber of chunks in which to find interface boundary"),

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


  if (jut >= 0) {
    gm_template.unrollAndTriangulate(string(triangulationHeuristic),
				     jut);
  } else {
    GMInfo gm_info;

    gm_info.M = M;

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
      if (gm_info.M != M) {
	error("ERROR: M given in partition file = %d, but M given on command line = %d\n",
		gm_info.M,M);
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

      if (gm_info.M != M) {
	error("ERROR: M given in partition file = %d, but M given on command line = %d\n",
		gm_info.M,M);
      }
      gm_template.triangulatePartitions(is,
					gm_info);
    }

      {
	printf("\n--- Printing final clique set and clique weights---\n");

	float p_maxWeight = -1.0;
	float p_totalWeight = -1.0; // starting flag
	for (unsigned i=0;i<gm_info.Pcliques.size();i++) {
	  float curWeight = gm_template.computeWeight(gm_info.Pcliques[i].nodes);
	  printf("   --- P curWeight = %f\n",curWeight);
	  if (curWeight > p_maxWeight) p_maxWeight = curWeight;
	  if (p_totalWeight == -1.0)
	    p_totalWeight = curWeight;
	  else
	    p_totalWeight = p_totalWeight + log10(1+pow(10,curWeight-p_totalWeight));
	}
	printf("  --- Prologue max clique weight = %f, total weight = %f\n",
	       p_maxWeight,p_totalWeight);

	float c_maxWeight = -1.0;
	float c_totalWeight = -1.0; // starting flag
	for (unsigned i=0;i<gm_info.Ccliques.size();i++) {
	  float curWeight = gm_template.computeWeight(gm_info.Ccliques[i].nodes);
	  printf("   --- C curWeight = %f\n",curWeight);
	  if (curWeight > c_maxWeight) c_maxWeight = curWeight;
	  if (c_totalWeight == -1.0)
	    c_totalWeight = curWeight;
	  else
	    c_totalWeight = c_totalWeight + log10(1+pow(10,curWeight-c_totalWeight));
	}
	printf("  --- Chunk max clique weight = %f, total Cx%d weight = %f, per-chunk total C weight = %f\n",
	       c_maxWeight,
	       M,
	       c_totalWeight,
	       c_totalWeight - log10((double)M));


	float e_maxWeight = -1.0;
	float e_totalWeight = -1.0; // starting flag
	for (unsigned i=0;i<gm_info.Ecliques.size();i++) {
	  float curWeight = gm_template.computeWeight(gm_info.Ecliques[i].nodes);
	  printf("   --- E curWeight = %f\n",curWeight);
	  if (curWeight > e_maxWeight) e_maxWeight = curWeight;
	  if (e_totalWeight == -1.0)
	    e_totalWeight = curWeight;
	  else
	    e_totalWeight = e_totalWeight + log10(1+pow(10,curWeight-e_totalWeight));
	}
	printf("  --- Epilogue max clique weight = %f, total weight = %f\n",
	       e_maxWeight,e_totalWeight);

	float maxWeight
	  = (p_maxWeight>c_maxWeight?p_maxWeight:c_maxWeight);
	maxWeight =
	  (maxWeight>e_maxWeight?maxWeight:e_maxWeight);
	float totalWeight = p_totalWeight;
	totalWeight += log10(1+pow(10,c_totalWeight-totalWeight));
	totalWeight += log10(1+pow(10,e_totalWeight-totalWeight));

	printf("--- Final set (P,Cx%d,E) has max clique weight = %f, total state space = %f ---\n",
	       M,
	       maxWeight,
	       totalWeight);

	// print out a couple of total state spaces for various unrollings
	printf("--- Total weight when unrolling %dx = %f ---\n",2*M-1,totalWeight);

	totalWeight += log10(1+pow(10,c_totalWeight-totalWeight));	
	printf("--- Total weight when unrolling %dx = %f ---\n",3*M-1,totalWeight);

	totalWeight += log10(1+pow(10,log10(3.0) + c_totalWeight-totalWeight));
	printf("--- Total weight when unrolling %dx = %f ---\n",6*M-1,totalWeight);

	totalWeight += log10(1+pow(10,log10(5.0) + c_totalWeight-totalWeight));
	printf("--- Total weight when unrolling %dx = %f ---\n",11*M-1,totalWeight);

	totalWeight += log10(1+pow(10,log10(10.0) + c_totalWeight-totalWeight));
	printf("--- Total weight when unrolling %dx = %f ---\n",21*M-1,totalWeight);

	totalWeight += log10(1+pow(10,log10(30.0) + c_totalWeight-totalWeight));
	printf("--- Total weight when unrolling %dx = %f ---\n",51*M-1,totalWeight);

	totalWeight += log10(1+pow(10,log10(50.0) + c_totalWeight-totalWeight));
	printf("--- Total weight when unrolling %dx = %f ---\n",101*M-1,totalWeight);

	totalWeight += log10(1+pow(10,log10(400.0) + c_totalWeight-totalWeight));
	printf("--- Total weight when unrolling %dx = %f ---\n",501*M-1,totalWeight);

	totalWeight += log10(1+pow(10,log10(500.0) + c_totalWeight-totalWeight));
	printf("--- Total weight when unrolling %dx = %f ---\n",1001*M-1,totalWeight);

      }
  }

  exit_program_with_status(0);
}
