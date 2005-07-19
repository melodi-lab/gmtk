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
//#include "spi.h"
#include "version.h"

VCID("$Header$")

#include "GMTK_FileParser.h"
#include "GMTK_RV.h"
#include "GMTK_DiscRV.h"
#include "GMTK_ContRV.h"
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


/*************************   INPUT TRAINABLE PARAMETER FILE HANDLING  *******************************************/
#define GMTK_ARG_CPP_CMD_OPTS

/*************************   INPUT STRUCTURE PARAMETER FILE HANDLING  *******************************************/
#define GMTK_ARG_STR_FILE
#define GMTK_ARG_INPUT_MASTER_FILE_OPT_ARG
#define GMTK_ARG_INPUT_TRAINABLE_PARAMS
#define GMTK_ARG_INPUT_TRI_FILE
#define GMTK_ARG_OUTPUT_TRI_FILE

/************************            TRIANGULATION OPTIONS             ******************************************/
#define GMTK_ARG_TRIANGULATION_OPTIONS
#define GMTK_ARG_LOAD_PARAMETERS
#define GMTK_ARG_NUM_BACKUP_FILES
#define GMTK_ARG_JTW_UB
#define GMTK_ARG_JT_OPTIONS
#define GMTK_ARG_CROSSOVER_OPTIONS

/*************************   INPUT TRAINABLE PARAMETER FILE HANDLING  *******************************************/
#define GMTK_ARG_ALLOC_DENSE_CPTS

/****************************         GENERAL OPTIONS             ***********************************************/
#define GMTK_ARG_SEED
#define GMTK_ARG_VERB


#define GMTK_ARGUMENTS_DEFINITION
#include "GMTK_Arguments.h"
#undef GMTK_ARGUMENTS_DEFINITION


Arg Arg::Args[] = {

#define GMTK_ARGUMENTS_DOCUMENTATION
#include "GMTK_Arguments.h"
#undef GMTK_ARGUMENTS_DOCUMENTATION

  // final one to signal the end of the list
  Arg()

};


void triangulateCrossover(
 BoundaryTriangulate& triangulator, 
 GMTemplate&          gm_template, 
 FileParser           fp,
 string               input_crossover_tri_file,
 string               output_crossover_tri_file,
 vector<MaxClique>&   input_P_triangulation,
 vector<MaxClique>&   input_C_triangulation,
 vector<MaxClique>&   input_E_triangulation
 );

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

  sprintf(buff,"jtwUB: %s, ",MYBS(JunctionTree::jtWeightUpperBound));
  res += buff;

  sprintf(buff,"jtwPUI: %f, ",JunctionTree::jtWeightPenalizeUnassignedIterated);
  res += buff;

  sprintf(buff,"jtwMC: %s, ",MYBS(JunctionTree::jtWeightMoreConservative));
  res += buff;

  sprintf(buff,"jtwSNSC: %f, ",JunctionTree::jtWeightSparseNodeSepScale);
  res += buff;

  sprintf(buff,"jtwDNSC: %f, ",JunctionTree::jtWeightDenseNodeSepScale);
  res += buff;

  sprintf(buff,"pfCobWeight: %f, ",MaxClique::continuousObservationPerFeaturePenalty);
  res += buff;
  
  sprintf(buff,"findBestBoundary: %s, ",MYBS(findBestBoundary));
  res += buff;
  
  sprintf(buff,"traverseFraction: %f, ",traverseFraction);
  res += buff;

  sprintf(buff,"noBoundaryMemoize: %s, ",MYBS(noBoundaryMemoize));
  res += buff;

  sprintf(buff,"forceLeftRight: %s, ",forceLeftRight);  
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
  string input_tri_file, output_tri_file;
  string input_crossover_tri_file, output_crossover_tri_file;

  ////////////////////////////////////////////
  // set things up so that if an FP exception
  // occurs such as an "invalid" (NaN), overflow
  // or divide by zero, we actually get a FPE
  ieeeFPsetup();
  set_new_handler(memory_error);
  InstallSignalHandlers();

  ////////////////////////////////////////////
  // parse arguments
  bool parse_was_ok = Arg::parse(argc,(char**)argv);
  if(!parse_was_ok) {
    Arg::usage(); exit(-1);
  }


#define GMTK_ARGUMENTS_CHECK_ARGS
#include "GMTK_Arguments.h"
#undef GMTK_ARGUMENTS_CHECK_ARGS

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
  // call with argument 'true' to do thorough graph check.
  fp.ensureValidTemplate(longStrCheck);

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

  BoundaryTriangulate triangulator(fp,maxNumChunksInBoundary,chunkSkip,traverseFraction);

  if (noBoundaryMemoize)
    triangulator.dontMemoizeBoundary();

  //////////////////////////////////////////////////////////////////////
  // Give warnings if crossover paramters are set for non-crossover
  // triangulation methods 
  //////////////////////////////////////////////////////////////////////
  if (string(triangulationHeuristic) != "crossover") {
    if (inputCrossoverTriangulatedFile != NULL) {
      warning("WARNING: inputCrossoverTriangulatedFile only used for triangulationHeuristic crossover");
    }
    if (outputCrossoverTriangulatedFile != NULL) {
      warning("WARNING: outputCrossoverTriangulatedFile is only used for triangulationHeuristic crossover");
    }
  }

  //////////////////////////////////////////////////////////////////////
  // Get name of input triangulation
  //////////////////////////////////////////////////////////////////////
  if (inputTriangulatedFile == NULL) {
    input_tri_file = string(strFileName) + GMTemplate::fileExtension;
  }
  else {
    input_tri_file = string(inputTriangulatedFile);
    if (fsize(input_tri_file.c_str()) == 0) {
      error("ERROR: triangulation file '%s' does not exist or is empty\n",
        input_tri_file.c_str() );
    }
  }

  if (inputCrossoverTriangulatedFile != NULL) {
    input_crossover_tri_file = string(inputCrossoverTriangulatedFile); 
    if (fsize(input_crossover_tri_file.c_str()) == 0) {
      error("ERROR: crossover triangulation file '%s' does not exist or is empty\n", input_crossover_tri_file.c_str() );
    }
  }

  //////////////////////////////////////////////////////////////////////
  // Name the output trifiles
  //////////////////////////////////////////////////////////////////////
  if (string(triangulationHeuristic) == "crossover") {

    if (outputTriangulatedFile == NULL) {
      output_tri_file = string(strFileName) + ".1." + GMTemplate::fileExtension;
    }
    else {
      output_tri_file = string(outputTriangulatedFile);
    }

    if (outputCrossoverTriangulatedFile == NULL) {
      output_crossover_tri_file = 
        string(strFileName) + ".2." + GMTemplate::fileExtension;
    }
    else {
      output_crossover_tri_file = string(outputCrossoverTriangulatedFile);
    }
  }
  else {
    if (outputTriangulatedFile == NULL) {
      output_tri_file = string(strFileName) + GMTemplate::fileExtension;
    }
    else {
      output_tri_file = string(outputTriangulatedFile);
    }
  }

  //////////////////////////////////////////////////////////////////////
  // Give error if one-edge is requested but other conflicting 
  // parameters are given.
  //////////////////////////////////////////////////////////////////////
  if (string(triangulationHeuristic) == "one-edge") {
    if (rePartition) {
      error("ERROR: Can not repartition graph when doing one-edge"); 
    }
    if (fsize(input_tri_file.c_str()) == 0) {
      error("ERROR: An inputTriangulatedFile is required when doing one-edge"); 
    }
  }

  //////////////////////////////////////////////////////////////////////
  // Give errors if crossover is requested but other conflicting 
  // parameters are given.
  //////////////////////////////////////////////////////////////////////
  if (string(triangulationHeuristic) == "crossover") {
    if (rePartition) {
      error("ERROR: Can not repartition graph when doing a crossover"); 
    }

printf("itf:%d  ictf:%d\n", fsize(input_tri_file.c_str()),
         fsize(input_crossover_tri_file.c_str()) );

    if ((inputCrossoverTriangulatedFile == NULL) ||
        (fsize(input_tri_file.c_str()) == 0)     ||
        (fsize(input_crossover_tri_file.c_str()) == 0)) { 
      error("ERROR: An inputTriangulatedFile and an inputCrossoverTriangulatedFile are required when doing a crossover"); 
    }
  }

  //////////////////////////////////////////////////////////////////////
  // Triangulate the graph 
  //////////////////////////////////////////////////////////////////////
  TimerClass* timer = NULL;
  timer = new TimerClass;
  // Initialize the timer if anyTimeTriangulate is selected
  if (anyTimeTriangulate != NULL) {
    time_t given_time;
    if (timeLimit != NULL) {
      time_t t1,t2;
      t1 = timer->parseTimeString( string(anyTimeTriangulate) );      
      t2 = timer->parseTimeString( string(timeLimit) );      
      given_time = min(t1,t2);
    } else
      given_time = timer->parseTimeString( string(anyTimeTriangulate) );
    if (given_time == 0) {
      error("ERROR: Must specify a non-zero amount of time for -anyTimeTriangulate"); 
    }
    infoMsg(IM::Low, "Triangulating for %d seconds\n", (int)given_time);
    timer->Reset(given_time);
  } else { 
    if (timeLimit != NULL) {
      time_t t1;
      t1 = timer->parseTimeString( string(timeLimit) );      
      if (t1 == 0) {
	error("ERROR: -timeLimit option must specify a non-zero amount of time"); 
      }
      infoMsg(IM::Low, "Running for no more than %d seconds\n", (int)t1);
      timer->Reset(t1);
    } else
      timer->DisableTimer();
  }
  triangulator.useTimer(timer);

  if (jut >= 0) {
    // then Just Unroll, Triangulate, and report on quality of triangulation.
    triangulator.unrollAndTriangulate(string(triangulationHeuristic),
				     jut);
  } else {

    GMTemplate gm_template(fp,maxNumChunksInBoundary,chunkSkip);

    BoundaryTriangulate::SavedGraph orgnl_P_graph;
    BoundaryTriangulate::SavedGraph orgnl_C_graph;
    BoundaryTriangulate::SavedGraph orgnl_E_graph;

    if (rePartition && !reTriangulate) {
      infoMsg(IM::Info,"NOTE: rePartition=T option forces -reTriangulate option to be true.\n");
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

      triangulator.saveCurrentNeighbors( gm_template.P.nodes, orgnl_P_graph );
      triangulator.saveCurrentNeighbors( gm_template.C.nodes, orgnl_C_graph );
      triangulator.saveCurrentNeighbors( gm_template.E.nodes, orgnl_E_graph );

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

      vector<MaxClique> input_P_triangulation;
      vector<MaxClique> input_C_triangulation;
      vector<MaxClique> input_E_triangulation;

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

      triangulator.saveCurrentNeighbors( gm_template.P.nodes, orgnl_P_graph );
      triangulator.saveCurrentNeighbors( gm_template.C.nodes, orgnl_C_graph );
      triangulator.saveCurrentNeighbors( gm_template.E.nodes, orgnl_E_graph );

      input_P_triangulation = gm_template.P.cliques;
      input_C_triangulation = gm_template.C.cliques;
      input_E_triangulation = gm_template.E.cliques;

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

        //////////////////////////////////////////////////////////////////////
        // Call the appropriate triangulation interface  
        //////////////////////////////////////////////////////////////////////
        if (string(triangulationHeuristic) == "crossover") {
          triangulateCrossover( triangulator, gm_template, fp, 
            input_crossover_tri_file, output_crossover_tri_file, 
            input_P_triangulation, 
            input_C_triangulation, input_E_triangulation);
        }
        else { 
          triangulator.triangulate(
            string(triangulationHeuristic),
            jtWeight,
            gm_template,
            input_P_triangulation,
            input_C_triangulation,
            input_E_triangulation,
            !noReTriP,!noReTriC,!noReTriE);
        }

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
      // Utilize both the partition information and clique information
      // already computed and contained in the file. This enables the
      // program to use external triangulation programs, where this
      // program ensures that the result is triangulated and where it
      // reports the quality of the triangulation.


      iDataStreamFile is(input_tri_file.c_str(),false,false);
      if (!fp.readAndVerifyGMId(is))
	error("ERROR: triangulation file '%s' does not match graph given in structure file '%s'\n",input_tri_file.c_str(),strFileName);

      gm_template.readPartitions(is);

      triangulator.saveCurrentNeighbors( gm_template.P.nodes, orgnl_P_graph );
      triangulator.saveCurrentNeighbors( gm_template.C.nodes, orgnl_C_graph );
      triangulator.saveCurrentNeighbors( gm_template.E.nodes, orgnl_E_graph );

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
	    p_totalWeight = log10add(curWeight,p_totalWeight);
	}
	printf("   --- Prologue max clique weight = %f, total weight = %f, jt_weight = %f\n",
	       p_maxWeight,p_totalWeight,
	       JunctionTree::junctionTreeWeight(gm_template.P.cliques,
						gm_template.PCInterface_in_P,
						NULL,&gm_template.PCInterface_in_P));
        printf("  --- Prologue triangulation is ");
        if (!triangulator.isEliminationGraph(orgnl_P_graph, 
             gm_template.P.nodes)) {
          printf("not ");
        }
        printf("an elimination graph\n"); 

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
	    c_totalWeight = log10add(curWeight,c_totalWeight);
	}
	printf("   --- Chunk max clique weight = %f, total Cx%d weight = %f, per-chunk total C weight = %f, jt_weight = %f\n",
	       c_maxWeight,
	       chunkSkip,
	       c_totalWeight,
	       c_totalWeight - log10((double)chunkSkip),
	       JunctionTree::junctionTreeWeight(gm_template.C.cliques,
						gm_template.CEInterface_in_C,
						&gm_template.PCInterface_in_C,
						&gm_template.CEInterface_in_C));
        printf("  --- Chunk triangulation is ");
        if (!triangulator.isEliminationGraph(orgnl_C_graph, 
             gm_template.C.nodes)) {
          printf("not ");
        }
        printf("an elimination graph\n"); 

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
	    e_totalWeight = log10add(curWeight,e_totalWeight);
	}
	const set <RV*> emptySet;
	printf("   --- Epilogue max clique weight = %f, total weight = %f, jt_weight = %f\n",
	       e_maxWeight,e_totalWeight,
	       JunctionTree::junctionTreeWeight(gm_template.E.cliques,
						emptySet,
						&gm_template.CEInterface_in_E,NULL));
        printf("  --- Epilogue triangulation is ");
        if (!triangulator.isEliminationGraph(orgnl_E_graph, 
             gm_template.E.nodes)) {
          printf("not ");
        }
        printf("an elimination graph\n"); 

	double maxWeight
	  = (p_maxWeight>c_maxWeight?p_maxWeight:c_maxWeight);
	maxWeight =
	  (maxWeight>e_maxWeight?maxWeight:e_maxWeight);
	double totalWeight = p_totalWeight;
	// log version of: totalWeight += c_totalWeight
	totalWeight = log10add(c_totalWeight,totalWeight);
	// log version of: totalWeight += e_totalWeight
	totalWeight = log10add(e_totalWeight,totalWeight);


	printf("--- Final set (P,Cx%d,E) has max clique weight = %f, total state space = %f ---\n",
	       chunkSkip,
	       maxWeight,
	       totalWeight);

	// print out a couple of total state spaces for various unrollings
	printf("--- Total weight when unrolling %dx = %f ---\n",maxNumChunksInBoundary+chunkSkip-1,totalWeight);

	totalWeight = log10add(c_totalWeight,totalWeight);	
	printf("--- Total weight when unrolling %dx = %f ---\n",maxNumChunksInBoundary+2*chunkSkip-1,totalWeight);

	totalWeight = log10add(log10(3.0) + c_totalWeight,totalWeight);
	printf("--- Total weight when unrolling %dx = %f ---\n",maxNumChunksInBoundary+5*chunkSkip-1,totalWeight);

	totalWeight = log10add(log10(5.0) + c_totalWeight,totalWeight);
	printf("--- Total weight when unrolling %dx = %f ---\n",maxNumChunksInBoundary+10*chunkSkip-1,totalWeight);

	totalWeight = log10add(log10(10.0) + c_totalWeight,totalWeight);
	printf("--- Total weight when unrolling %dx = %f ---\n",maxNumChunksInBoundary+20*chunkSkip-1,totalWeight);

	totalWeight = log10add(log10(30.0) + c_totalWeight,totalWeight);
	printf("--- Total weight when unrolling %dx = %f ---\n",maxNumChunksInBoundary+50*chunkSkip-1,totalWeight);

	totalWeight = log10add(log10(50.0) + c_totalWeight,totalWeight);
	printf("--- Total weight when unrolling %dx = %f ---\n",maxNumChunksInBoundary+100*chunkSkip-1,totalWeight);

	totalWeight = log10add(log10(400.0) + c_totalWeight,totalWeight);
	printf("--- Total weight when unrolling %dx = %f ---\n",maxNumChunksInBoundary+500*chunkSkip-1,totalWeight);

	totalWeight = log10add(log10(500.0) + c_totalWeight,totalWeight);
	printf("--- Total weight when unrolling %dx = %f ---\n",maxNumChunksInBoundary+1000*chunkSkip-1,totalWeight);

      }
    delete timer;
  }

  exit_program_with_status(0);
}



void triangulateCrossover(
 BoundaryTriangulate& triangulator, 
 GMTemplate&          gm_template, 
 FileParser           fp,
 string               input_crossover_tri_file,
 string               output_crossover_tri_file,
 vector<MaxClique>&   input_P_triangulation,
 vector<MaxClique>&   input_C_triangulation,
 vector<MaxClique>&   input_E_triangulation
 )
{
  vector<MaxClique> crossover_P_tri;
  vector<MaxClique> crossover_C_tri;
  vector<MaxClique> crossover_E_tri;

  GMTemplate crossover_gm_template(fp,maxNumChunksInBoundary,chunkSkip);
  iDataStreamFile cis(input_crossover_tri_file.c_str(), false, false);

  if (!fp.readAndVerifyGMId(cis)) {
    error("ERROR: crossover triangulation file '%s' does not match graph given in structure file '%s'\n", input_crossover_tri_file.c_str(), strFileName);
  }

  crossover_gm_template.readPartitions(cis);
  crossover_gm_template.readMaxCliques(cis);
  crossover_gm_template.triangulatePartitionsByCliqueCompletion();

  crossover_P_tri = crossover_gm_template.P.cliques;
  crossover_C_tri = crossover_gm_template.C.cliques;
  crossover_E_tri = crossover_gm_template.E.cliques;

  triangulator.triangulateCrossover(
    gm_template, 
    input_P_triangulation, input_C_triangulation, input_E_triangulation,
    crossover_gm_template,
    crossover_P_tri, crossover_C_tri, crossover_E_tri, 
    crossoverProbability, mutateProbability, 
    !noReTriP, !noReTriC, !noReTriE);

  crossover_gm_template.triangulatePartitionsByCliqueCompletion();

  backupTriFile(output_crossover_tri_file);
  oDataStreamFile cos(output_crossover_tri_file.c_str());
  fp.writeGMId(cos);
  string clStr;
  createCommandLineOptionString(clStr);
  crossover_gm_template.writePartitions(cos, clStr);
  crossover_gm_template.writeMaxCliques(cos);

  triangulator.ensurePartitionsAreChordal(crossover_gm_template);
}


