/*
 * gmtkKernel.cc
 *  Use a DGM as a Kernel (e.g., Fisher Kernel, or something else).
 *  In other words, for each observation file, we will write out 
 *  a vector that is on the order of the number of current system parameters.
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
#include "version.h"

#include "GMTK_WordOrganization.h"

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
#include "GMTK_BoundaryTriangulate.h"
#include "GMTK_JunctionTree.h"
#include "GMTK_MaxClique.h"


/*****************************   OBSERVATION INPUT FILE HANDLING   **********************************************/
#define GMTK_ARG_OBS_FILES

/*************************   INPUT TRAINABLE PARAMETER FILE HANDLING  *******************************************/
#define GMTK_ARG_CPP_CMD_OPTS
#define GMTK_ARG_INPUT_MASTER_FILE
#define GMTK_ARG_OUTPUT_MASTER_FILE
#define GMTK_ARG_INPUT_TRAINABLE_PARAMS
#define GMTK_ARG_OUTPUT_TRAINABLE_PARAMS
#define GMTK_ARG_WPAEEI
#define GMTK_ARG_ALLOC_DENSE_CPTS
#define GMTK_ARG_CPT_NORM_THRES

/*************************   INPUT STRUCTURE PARAMETER FILE HANDLING  *******************************************/
#define GMTK_ARG_STR_FILE
#define GMTK_ARG_TRI_FILE
#define GMTK_ARG_CHECK_TRI_FILE_CARD
#define GMTK_ARG_JT_INFO_FILE
#define GMTK_ARG_JTW_UB


/*************************   CONTINUOUS RANDOM VARIABLE OPTIONS       *******************************************/
#define GMTK_ARG_VAR_FLOOR
#define GMTK_ARG_VAR_FLOOR_ON_READ


/*************************          BEAM PRUNING OPTIONS              *******************************************/
#define GMTK_ARG_CBEAM
#define GMTK_ARG_CPBEAM
#define GMTK_ARG_CKBEAM
#define GMTK_ARG_CCBEAM
#define GMTK_ARG_CRBEAM
#define GMTK_ARG_CMBEAM
#define GMTK_ARG_SBEAM
#define GMTK_ARG_EBEAM

/*************************          MEMORY MANAGEMENT OPTIONS         *******************************************/
#define GMTK_ARG_HASH_LOAD_FACTOR
#define GMTK_ARG_STORE_DETERMINISTIC_CHILDREN
#define GMTK_ARG_CLEAR_CLIQUE_VAL_MEM


/****************************      FILE RANGE OPTIONS             ***********************************************/
#define GMTK_ARG_DCDRNG
#define GMTK_ARG_START_END_SKIP

/****************************         GENERAL OPTIONS             ***********************************************/
#define GMTK_ARG_SEED
#define GMTK_ARG_VERB
#define GMTK_ARG_HELP
#define GMTK_ARG_VERSION

/****************************         INFERENCE OPTIONS           ***********************************************/
#define GMTK_ARG_ISLAND
#define GMTK_ARG_CLIQUE_TABLE_NORMALIZE
#define GMTK_ARG_CE_SEP_DRIVEN
#define GMTK_ARG_COMPONENT_CACHE
#define GMTK_ARG_CLIQUE_VAR_ITER_ORDERS
#define GMTK_ARG_JT_OPTIONS
#define GMTK_ARG_VE_SEPS

/****************************         EM TRAINING OPTIONS         ***********************************************/
#define GMTK_ARG_KERNEL_PARAMS

/************************  OBSERVATION MATRIX TRANSFORMATION OPTIONS   ******************************************/
#define GMTK_ARG_OBS_MATRIX_XFORMATION


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


/*
 *
 * definition of needed global arguments
 *
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

  CODE_TO_COMPUTE_ENDIAN;

  ////////////////////////////////////////////
  // parse arguments
  bool parse_was_ok = Arg::parse(argc,(char**)argv);
  if(!parse_was_ok) {
    Arg::usage(); exit(-1);
  }

#define GMTK_ARGUMENTS_CHECK_ARGS
#include "GMTK_Arguments.h"
#undef GMTK_ARGUMENTS_CHECK_ARGS


  // Check if we want to do the Fisher kernel rather than the
  // accumulator kernel.
  EMable::fisherKernelMode = fisherKernelP;

  globalObservationMatrix.openFiles(nfiles,
				    (const char**)&ofs,
				    (const char**)&frs,
				    (const char**)&irs,
				    (unsigned*)&nfs,
				    (unsigned*)&nis,
				    (unsigned*)&ifmts,
				    (bool*)&iswp,
				    startSkip,
				    endSkip,
				    Cpp_If_Ascii,
				    cppCommandOptions,
				    (const char**)&postpr,  //Frame_Range_Str,
				    Action_If_Diff_Num_Frames,
				    Action_If_Diff_Num_Sents,
				    Per_Stream_Transforms,
				    Post_Transforms,
				    Ftr_Combo,
				    (const char**)&sr,
				    (const char**)&prepr,
				    gpr_str
				    );


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
  // comment for now Sun Jan 11 09:47:23 2004
  GM_Parms.finalizeParameters();
  GM_Parms.markObjectsToNotTrain(objsToNotUtilizeFile,cppCommandOptions);

  /////////////////////////////
  // read in the structure of the GM, this will
  // die if the file does not exist.
  FileParser fp(strFileName,cppCommandOptions);
  infoMsg(IM::Tiny,"Finished reading in all parameters and structures\n");

  // parse the file
  fp.parseGraphicalModel();
  // create the rv variable objects
  fp.createRandomVariableGraph();
  // Make sure that there are no directed loops in the graph.
  fp.ensureValidTemplate();

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
  fp.checkConsistentWithGlobalObservationStream();
  GM_Parms.checkConsistentWithGlobalObservationStream();

  GM_Parms.setStride(globalObservationMatrix.stride());

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
  GMTemplate gm_template(fp);
  {
    // do this in scope so that is gets deleted now rather than later.
    iDataStreamFile is(tri_file.c_str());
    if (!fp.readAndVerifyGMId(is,checkTriFileCards))
      error("ERROR: triangulation file '%s' does not match graph given in structure file '%s'\n",tri_file.c_str(),strFileName);
    gm_template.readPartitions(is);
    gm_template.readMaxCliques(is);
  }
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
  if (jtFileName != NULL)
    myjt.printAllJTInfo(jtFileName);
  infoMsg(IM::Default,"DONE creating Junction Tree\n"); fflush(stdout);
  ////////////////////////////////////////////////////////////////////

  if (randomizeParams) {
    infoMsg(IM::Default,"WARNING: GMTK is randomizing all trainable parameters and writing them to file 'random.gmp'\n");
    GM_Parms.makeRandom();
    oDataStreamFile of("random.gmp");
    GM_Parms.writeTrainable(of);
  }

  if (globalObservationMatrix.numSegments()==0) {
    infoMsg(IM::Default,"ERROR: no segments are available in observation file. Exiting...");
    exit_program_with_status(0);
  }

  Range* dcdrng = new Range(dcdrng_str,0,globalObservationMatrix.numSegments());
  if (dcdrng->length() <= 0) {
    infoMsg(IM::Default,"Training range '%s' specifies empty set. Exiting...\n",
	  dcdrng_str);
    exit_program_with_status(0);
  }



  struct rusage rus; /* starting time */
  struct rusage rue; /* ending time */
  getrusage(RUSAGE_SELF,&rus);

  // Now, do CE/DE iterations, writing out the parameter "increments" one per segment.

  Range::iterator* dcdrng_it = new Range::iterator(dcdrng->begin());
  {
    oDataStreamFile outf(storeFeatureFile,transFileIsBinary);
    bool firstTime = true;
    while (!dcdrng_it->at_end()) {
      const unsigned segment = (unsigned)(*(*dcdrng_it));
      if (globalObservationMatrix.numSegments() < (segment+1)) 
	error("ERROR: only %d segments in file, segment must be in range [%d,%d]\n",
	      globalObservationMatrix.numSegments(),
	      0,globalObservationMatrix.numSegments()-1);

      const unsigned numFrames = GM_Parms.setSegment(segment);

      logpr data_prob = 1.0;
      GM_Parms.emInitAccumulators(firstTime);

      unsigned numUsableFrames;
      if (island) {
	myjt.collectDistributeIsland(numFrames,
				     numUsableFrames,
				     base,
				     lst,
				     true, // run EM algorithm,
				     false, // run Viterbi algorithm
				     localCliqueNormalization);
	printf("Segment %d, after Island, log(prob(evidence)) = %f, per frame =%f, per numUFrams = %f, ",
	       segment,
	       myjt.curProbEvidenceIsland().val(),
	       myjt.curProbEvidenceIsland().val()/numFrames,
	       myjt.curProbEvidenceIsland().val()/numUsableFrames);
	data_prob = myjt.curProbEvidenceIsland();
      } else {
	numUsableFrames = myjt.unroll(numFrames);
	infoMsg(IM::Low,"Collecting Evidence\n");
	myjt.collectEvidence();
	infoMsg(IM::Low,"Done Collecting Evidence\n");
	data_prob = myjt.probEvidence();

	infoMsg(IM::Low,"Distributing Evidence\n");
	myjt.distributeEvidence();
	infoMsg(IM::Low,"Done Distributing Evidence\n");

	printf("Segment %d, after CE/DE, log(prob(evidence)) = %f, per frame =%f, per numUFrams = %f, ",
	       segment,
	       data_prob.val(),
	       data_prob.val()/numFrames,
	       data_prob.val()/numUsableFrames);
      
	myjt.emIncrement(data_prob,localCliqueNormalization);

      }


      printf("writing %s-kernel feature space vector ...\n",(fisherKernelP?"Fisher":"accumulator"));
      if (annotateTransformationOutput) {
	char buff[1024];
	sprintf(buff,"Segment %d : %d frames, %d usable frames, log(PE) = %f",segment,numFrames,numUsableFrames,data_prob.val());
	outf.write(buff);
	outf.nl();
      };
      outf.write(data_prob.val());
      GM_Parms.emWriteUnencodedAccumulators(outf,writeLogVals);
      outf.nl();
      
      (*dcdrng_it)++;
      firstTime = false;
    }
  }

  getrusage(RUSAGE_SELF,&rue);
  if (IM::messageGlb(IM::Default)) { 
    infoMsg(IM::Default,"### Final time (seconds) just for kernel computing stage: ");
    double userTime,sysTime;
    reportTiming(rus,rue,userTime,sysTime,stdout);
  }

  exit_program_with_status(0);
}
