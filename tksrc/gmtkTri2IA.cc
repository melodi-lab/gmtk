/*
 * gmtkJT.cc
 * produce a junction tree
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (C) 2001 Jeff Bilmes
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
#include <time.h>

#include "general.h"
#include "error.h"
#include "debug.h"
#include "rand.h"
#include "arguments.h"
#include "ieeeFPsetup.h"

#include "GMTK_WordOrganization.h"

#include "GMTK_FileParser.h"
#include "GMTK_RV.h"
#include "GMTK_DiscRV.h"
#include "GMTK_ContRV.h"
#include "GMTK_GMTemplate.h"
#include "GMTK_GMParms.h"
#include "GMTK_MixtureCommon.h"
#include "GMTK_GaussianComponent.h"
#include "GMTK_MeanVector.h"
#include "GMTK_DiagCovarVector.h"
#include "GMTK_Dlinks.h"
#include "GMTK_DlinkMatrix.h"
#include "GMTK_ProgramDefaultParms.h"
#include "GMTK_BoundaryTriangulate.h"
#include "GMTK_JunctionTree.h"
#include "GMTK_MaxClique.h"

#include "GMTK_ObservationSource.h"

VCID(HGID)




/*************************   INPUT TRAINABLE PARAMETER FILE HANDLING  *******************************************/
#define GMTK_ARG_INPUT_TRAINABLE_FILE_HANDLING
#define GMTK_ARG_CPP_CMD_OPTS
#define GMTK_ARG_INPUT_MASTER_FILE
#define GMTK_ARG_DLOPEN_MAPPERS
#define GMTK_ARG_INPUT_TRAINABLE_PARAMS
#define GMTK_ARG_ALLOC_DENSE_CPTS
#define GMTK_ARG_CPT_NORM_THRES

/*************************   INPUT STRUCTURE PARAMETER FILE HANDLING  *******************************************/
#define GMTK_ARG_INPUT_MODEL_FILE_HANDLING
#define GMTK_ARG_STR_FILE
#define GMTK_ARG_TRI_FILE
#define GMTK_ARG_IA_FILE
#define GMTK_ARG_CHECK_TRI_FILE_CARD
#define GMTK_ARG_JT_INFO_FILE
#define GMTK_ARG_JTW_UB
#define GMTK_ARG_LATTICE_PARAMS

/****************************         INFERENCE OPTIONS           ***********************************************/
#define GMTK_ARG_CLIQUE_VAR_ITER_ORDERS

/*************************   CONTINUOUS RANDOM VARIABLE OPTIONS       *******************************************/
#define GMTK_ARG_CONTINUOUS_RANDOM_VAR_OPTIONS
#define GMTK_ARG_VAR_FLOOR
#define GMTK_ARG_VAR_FLOOR_ON_READ

/*************************          MEMORY MANAGEMENT OPTIONS         *******************************************/
#define GMTK_ARG_MEMORY_MANAGEMENT_OPTIONS
#define GMTK_ARG_MEM_GROWTH


/****************************         GENERAL OPTIONS             ***********************************************/
#define GMTK_ARG_GENERAL_OPTIONS
#define GMTK_ARG_SEED
#define GMTK_ARG_SKIP_STARTUP_CHECKS
#define GMTK_ARG_VERB
#define GMTK_ARG_HELP
#define GMTK_ARG_VERSION
#define GMTK_ARG_CLIQUE_PRINT




#define GMTK_ARGUMENTS_DEFINITION
#include "GMTK_Arguments.h"
#undef GMTK_ARGUMENTS_DEFINITION


Arg Arg::Args[] = {

#define GMTK_ARGUMENTS_DOCUMENTATION
#include "GMTK_Arguments.h"
#undef GMTK_ARGUMENTS_DOCUMENTATION

#if 0
  Arg("boostVerbosity",Arg::Opt,boostVerbosity,"Verbosity (0 <= v <= 100) during boost verb partitions"),
  Arg("boostRng",Arg::Opt,boostVerbosityRng,"Range to boost verbosity"),
#endif

  // final one to signal the end of the list
  Arg()

};


/*
 *
 * definition of needed global arguments
 *
 */
RAND rnd(seedme);
GMParms GM_Parms;

ObservationSource *globalObservationMatrix;

int
main(int argc,char*argv[]) {
  try { // for catching std::bad_alloc(), indicating memory exhaustion

    { // use double so that we can destruct objects at end.

      ////////////////////////////////////////////
      // set things up so that if an FP exception
      // occurs such as an "invalid" (NaN), overflow
      // or divide by zero, we actually get a FPE
      ieeeFPsetup();
      set_new_handler(memory_error);
#if 0
      CODE_TO_COMPUTE_ENDIAN;
#endif
      
      ////////////////////////////////////////////
      // parse arguments
      bool parse_was_ok = Arg::parse(argc,(char**)argv,
				     "\nThis program computes the probability of evidence, and can also\n"
				     "compute the posterior distributions of the hidden variables\n");
      if(!parse_was_ok) {
	// Arg::usage(); 
	exit(-1);
      }
      
      infoMsg(IM::Max,"Finished parsing arguments\n");
      
#define GMTK_ARGUMENTS_CHECK_ARGS
#include "GMTK_Arguments.h"
#undef GMTK_ARGUMENTS_CHECK_ARGS
      
      
      
      
      /////////////////////////////////////////////
      // read in all the parameters
      
      dlopenDeterministicMaps(dlopenFilenames, MAX_NUM_DLOPENED_FILES);
      if (inputMasterFile) {
	// flat, where everything is contained in one file, always ASCII
	infoMsg(IM::Max,"Reading master file...\n");
	iDataStreamFile pf(inputMasterFile,false,true,cppCommandOptions);
	GM_Parms.read(pf);
	infoMsg(IM::Max,"Finished reading master file.\n");
      }
      if (inputTrainableParameters) {
	// flat, where everything is contained in one file
	infoMsg(IM::Max,"Reading trainable file...\n");
	iDataStreamFile pf(inputTrainableParameters,binInputTrainableParameters,true,cppCommandOptions);
	GM_Parms.readTrainable(pf);
	infoMsg(IM::Max,"Finished reading trainable file.\n");
      }
      GM_Parms.finalizeParameters();

      
      /////////////////////////////
      // read in the structure of the GM, this will
      // die if the file does not exist.
      infoMsg(IM::Max,"Reading structure file...\n");
      FileParser fp(strFileName,cppCommandOptions);
      infoMsg(IM::Tiny,"Finished reading in all parameters and structures\n");
      
      // parse the file
      infoMsg(IM::Max,"Parsing structure file...\n");
      fp.parseGraphicalModel();
      
      // create the rv variable objects
      infoMsg(IM::Max,"Creating rv objects...\n");
      fp.createRandomVariableGraph();
      

      // Make sure that there are no directed loops in the graph.
      infoMsg(IM::Max,"Checking template...\n");
      fp.ensureValidTemplate();
      
      // link the RVs with the parameters.
      infoMsg(IM::Max,"Allocating cpts...\n");
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


      // Utilize both the partition information and elimination order
      // information already computed and contained in the file. This
      // enables the program to use external triangulation programs,
      // where this program ensures that the result is triangulated
      // and where it reports the quality of the triangulation.
      
      string tri_file;
      if (triFileName == NULL) 
	tri_file = string(strFileName) + GMTemplate::fileExtension;
      else 
	tri_file = string(triFileName);
      
      string ia_file; 
      if (iaFileName == NULL)
	ia_file = string(strFileName) + string(".ia");
      else
	ia_file = string(iaFileName);
      
      infoMsg(IM::Max,"Creating template...\n");
      GMTemplate gm_template(fp);
      {
	infoMsg(IM::Max,"Reading triangulation file...\n");
	
	// do this in scope so that is gets deleted now rather than later.
	iDataStreamFile is(tri_file.c_str());
	if (!fp.readAndVerifyGMId(is,checkTriFileCards))
	  error("ERROR: triangulation file '%s' does not match graph given in structure file '%s'\n",tri_file.c_str(),strFileName);
	
	gm_template.readPartitions(is);
	gm_template.readMaxCliques(is);
	
      }

      infoMsg(IM::Max,"Triangulating graph...\n");
      gm_template.triangulatePartitionsByCliqueCompletion();
      if (1) { 
	// check that graph is indeed triangulated.
	// TODO: perhaps take this check out so that inference code does
	// not need to link to the triangulation code (either that, or put
	// the triangulation check in a different file, so that we only
	// link to tri check code).
	BoundaryTriangulate triangulator(fp,
					 gm_template.maxNumChunksInBoundary(),
					 gm_template.chunkSkip(),1.0);
	triangulator.ensurePartitionsAreChordal(gm_template);
      }
      

      ////////////////////////////////////////////////////////////////////
      // CREATE JUNCTION TREE DATA STRUCTURES
      infoMsg(IM::Default,"Creating Junction Tree\n"); fflush(stdout);
      JunctionTree myjt(gm_template);
      
      myjt.setUpDataStructures(varPartitionAssignmentPrior,varCliqueAssignmentPrior);
      
      myjt.prepareForUnrolling();
      
      myjt.printAllIAInfo(ia_file.c_str(), writeComments);
      
      infoMsg(IM::Default,"DONE creating Junction Tree\n"); fflush(stdout);
      ////////////////////////////////////////////////////////////////////
      
    } // close brace to cause a destruct on valid end of program.
    exit_program_with_status(0); 
  } catch (std::bad_alloc const &e) {
    memory_error();
  }
}

