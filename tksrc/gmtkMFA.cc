/*
 * gmtkMFA.cc
 * triangulate a graph wtih mean field approximation interface factorization
 *
 * Written by Richard Rogers <rprogers@uw.edu>
 *
 * Copyright (C) 2016 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 *
 */

#if HAVE_CONFIG_H
#include <config.h>
#endif
#if HAVE_HG_H
#include "hgstamp.h"
#endif

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

VCID(HGID)


#include "GMTK_FileParser.h"
#include "GMTK_RV.h"
#include "GMTK_DiscRV.h"
#include "GMTK_ContRV.h"
#include "GMTK_GMTemplate.h"
#include "GMTK_GMParms.h"

#include "GMTK_ObservationSource.h"
#include "GMTK_FileSource.h"
#include "GMTK_ASCIIFile.h"
#include "GMTK_Stream.h"

#include "GMTK_MixtureCommon.h"
#include "GMTK_GaussianComponent.h"
#include "GMTK_MeanVector.h"
#include "GMTK_DiagCovarVector.h"
#include "GMTK_DlinkMatrix.h"
#include "GMTK_ProgramDefaultParms.h"
#include "GMTK_Signals.h"
#include "GMTK_BoundaryTriangulate.h"

#include "GMTK_SectionScheduler.h"
#include "GMTK_LinearSectionScheduler.h"

/*************************   INPUT TRAINABLE PARAMETER FILE HANDLING  *******************************************/
#define GMTK_ARG_CPP_CMD_OPTS

/*************************   INPUT STRUCTURE PARAMETER FILE HANDLING  *******************************************/
#define GMTK_ARG_STR_FILE
#define GMTK_ARG_INPUT_MASTER_FILE_OPT_ARG
#define GMTK_ARG_DLOPEN_MAPPERS
#define GMTK_ARG_INPUT_TRAINABLE_PARAMS
#define GMTK_ARG_INPUT_TRI_FILE
#define GMTK_ARG_CHECK_TRI_FILE_CARD
#define GMTK_ARG_OUTPUT_TRI_FILE

/****************************         INFERENCE OPTIONS           ***********************************************/
#define GMTK_ARG_INFERENCE_OPTIONS
#define GMTK_ARG_CLIQUE_VAR_ITER_ORDERS
#define GMTK_ARG_JT_OPTIONS

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
#define GMTK_ARG_SKIP_STARTUP_CHECKS
#define GMTK_ARG_VERSION
#define GMTK_ARG_VERB
#define GMTK_ARG_HELP

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

  sprintf(buff,"jtwUB: %s, ",MYBS(SectionScheduler::jtWeightUpperBound));
  res += buff;

  sprintf(buff,"jtwPUI: %f, ",SectionScheduler::jtWeightPenalizeUnassignedIterated);
  res += buff;

  sprintf(buff,"jtwMC: %s, ",MYBS(SectionScheduler::jtWeightMoreConservative));
  res += buff;

  sprintf(buff,"jtwSNSC: %f, ",SectionScheduler::jtWeightSparseNodeSepScale);
  res += buff;

  sprintf(buff,"jtwDNSC: %f, ",SectionScheduler::jtWeightDenseNodeSepScale);
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
#if 0
ObservationMatrix globalObservationMatrix;
#else
FileSource fileSource;
FileSource *gomFS = &fileSource;
ObservationSource *globalObservationMatrix = &fileSource;
#endif

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
  InstallSignalHandlersTime();

  ////////////////////////////////////////////
  // parse arguments
  bool parse_was_ok = Arg::parse(argc,(char**)argv,
"\nThis program analyzes the graphical structure of a model to determine\n"
"an efficient way to perform inference on it (Mean Field Approximation version).\n");
  if(!parse_was_ok) {
    Arg::usage(); exit(-1);
  }


#define GMTK_ARGUMENTS_CHECK_ARGS
#include "GMTK_Arguments.h"
#undef GMTK_ARGUMENTS_CHECK_ARGS

  /////////////////////////////////////////////
  if (loadParameters) {
    // read in all the parameters
    dlopenDeterministicMaps(dlopenFilenames, MAX_NUM_DLOPENED_FILES);
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
      error("ERROR: inference architecture file '%s' does not exist or is empty\n",
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

    printf("itf:%ld  ictf:%ld\n", (long)fsize(input_tri_file.c_str()),
	   (long)fsize(input_crossover_tri_file.c_str()) );

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







  GMTemplate gm_template(fp,maxNumChunksInBoundary,chunkSkip);

  BoundaryTriangulate::SavedGraph orgnl_P_graph;
  BoundaryTriangulate::SavedGraph orgnl_C_graph;
  BoundaryTriangulate::SavedGraph orgnl_E_graph;

#if 0  
    if (rePartition && !reTriangulate) {
      infoMsg(IM::Info,"NOTE: rePartition=T option forces -reTriangulate option to be true.\n");
      reTriangulate = true;
    }

    // first check if tri_file exists
    if (rePartition || fsize(input_tri_file.c_str()) == 0) {
      // Then do everything (both partition & triangulation)
#endif

  // run partition given options
  triangulator.findPartitions(string(boundaryHeuristic),
			      string(forceLeftRight),
			      string(triangulationHeuristic),
			      findBestBoundary,
			      gm_template);

  triangulator.saveCurrentNeighbors( gm_template.P.nodes, orgnl_P_graph );
  triangulator.saveCurrentNeighbors( gm_template.C.nodes, orgnl_C_graph );
  triangulator.saveCurrentNeighbors( gm_template.E.nodes, orgnl_E_graph );
  
  triangulator.triangulate(string(triangulationHeuristic),
			   jtWeight,
			   gm_template);

  SectionScheduler *myjt = new LinearSectionScheduler(gm_template, fp, gomFS);

  myjt->setUpDataStructures(varSectionAssignmentPrior,varCliqueAssignmentPrior);

  printf("writing IA to %s\n", output_tri_file.c_str());
  backupTriFile(output_tri_file);
  myjt->printAllIAInfo(output_tri_file.c_str(), writeComments);

  exit_program_with_status(0);
}


#if 0
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

  if (!fp.readAndVerifyGMId(cis,checkTriFileCards)) {
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
  cos.setWriteCommentsStatus(writeComments);
  fp.writeGMId(cos);
  string clStr;
  createCommandLineOptionString(clStr);
  crossover_gm_template.writePartitions(cos, clStr);
  crossover_gm_template.writeMaxCliques(cos);

  triangulator.ensurePartitionsAreChordal(crossover_gm_template);
}
#endif

