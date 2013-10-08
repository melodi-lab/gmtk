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

#include "general.h"
#include "error.h"
#include "debug.h"
#include "rand.h"
#include "arguments.h"
#include "ieeeFPsetup.h"
#include "version.h"

#include "GMTK_WordOrganization.h"

VCID(HGID)


#include "GMTK_FileParser.h"
#include "GMTK_RV.h"
#include "GMTK_DiscRV.h"
#include "GMTK_ContRV.h"
#include "GMTK_GMTemplate.h"
#include "GMTK_GMParms.h"
#if 0
#  include "GMTK_ObservationMatrix.h"
#else
#  include "GMTK_ObservationSource.h"
#  include "GMTK_FileSource.h"
#  include "GMTK_CreateFileSource.h"
#  include "GMTK_ASCIIFile.h"
#  include "GMTK_FlatASCIIFile.h"
#  include "GMTK_PFileFile.h"
#  include "GMTK_HTKFile.h"
#  include "GMTK_HDF5File.h"
#  include "GMTK_BinaryFile.h"
#  include "GMTK_Filter.h"
#  include "GMTK_Stream.h"
#endif
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
#define GMTK_ARG_INPUT_TRAINABLE_FILE_HANDLING
#define GMTK_ARG_CPP_CMD_OPTS
#define GMTK_ARG_INPUT_MASTER_FILE
#define GMTK_ARG_OUTPUT_MASTER_FILE
#define GMTK_ARG_DLOPEN_MAPPERS
#define GMTK_ARG_INPUT_TRAINABLE_PARAMS
#define GMTK_ARG_OUTPUT_TRAINABLE_PARAMS
#define GMTK_ARG_WPAEEI
#define GMTK_ARG_ALLOC_DENSE_CPTS
#define GMTK_ARG_CPT_NORM_THRES

/*************************   INPUT STRUCTURE PARAMETER FILE HANDLING  *******************************************/
#define GMTK_ARG_INPUT_MODEL_FILE_HANDLING
#define GMTK_ARG_STR_FILE
#define GMTK_ARG_TRI_FILE
#define GMTK_ARG_CHECK_TRI_FILE_CARD
#define GMTK_ARG_JT_INFO_FILE
#define GMTK_ARG_JTW_UB
#define GMTK_ARG_LATTICE_PARAMS


/*************************   CONTINUOUS RANDOM VARIABLE OPTIONS       *******************************************/
#define GMTK_ARG_CONTINUOUS_RANDOM_VAR_OPTIONS
#define GMTK_ARG_VAR_FLOOR
#define GMTK_ARG_VAR_FLOOR_ON_READ


/*************************          BEAM PRUNING OPTIONS              *******************************************/
#define GMTK_ARG_BEAM_PRUNING_OPTIONS
#define GMTK_ARG_CBEAM
#define GMTK_ARG_CPBEAM
#define GMTK_ARG_CKBEAM
#define GMTK_ARG_CCBEAM
#define GMTK_ARG_CRBEAM
#define GMTK_ARG_CMBEAM
#define GMTK_ARG_SBEAM
#define GMTK_ARG_EBEAM

/*************************          MEMORY MANAGEMENT OPTIONS         *******************************************/
#define GMTK_ARG_MEMORY_MANAGEMENT_OPTIONS
#define GMTK_ARG_HASH_LOAD_FACTOR
#define GMTK_ARG_STORE_DETERMINISTIC_CHILDREN
#define GMTK_ARG_CLEAR_CLIQUE_VAL_MEM
#define GMTK_ARG_MEM_GROWTH
#define GMTK_ARG_USE_MMAP

/****************************      FILE RANGE OPTIONS             ***********************************************/
#define GMTK_ARG_FILE_RANGE_OPTIONS
#define GMTK_ARG_TRRNG
#define GMTK_ARG_START_END_SKIP

/****************************         GENERAL OPTIONS             ***********************************************/
#define GMTK_ARG_GENERAL_OPTIONS
#define GMTK_ARG_SEED
#define GMTK_ARG_SKIP_STARTUP_CHECKS
#define GMTK_ARG_VERB
#define GMTK_ARG_HELP
#define GMTK_ARG_VERSION

/****************************         INFERENCE OPTIONS           ***********************************************/
#define GMTK_ARG_INFERENCE_OPTIONS
#define GMTK_ARG_ISLAND
#define GMTK_ARG_DEBUG_PART_RNG
#define GMTK_ARG_DEBUG_INCREMENT
#define GMTK_ARG_CLIQUE_TABLE_NORMALIZE
#define GMTK_ARG_CE_SEP_DRIVEN
#define GMTK_ARG_COMPONENT_CACHE
#define GMTK_ARG_CLIQUE_VAR_ITER_ORDERS
#define GMTK_ARG_JT_OPTIONS
#define GMTK_ARG_VE_SEPS
#define GMTK_ARG_FAIL_ON_ZERO_CLIQUE

/****************************         EM TRAINING OPTIONS         ***********************************************/
#define GMTK_ARG_EM_TRAINING_OPTIONS
#define GMTK_ARG_EM_TRAINING_PARAMS

/************************  OBSERVATION MATRIX TRANSFORMATION OPTIONS   ******************************************/
#define GMTK_ARG_OBS_MATRIX_OPTIONS
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
#if 0
ObservationMatrix globalObservationMatrix;
#endif

FileSource *gomFS;
ObservationSource *globalObservationMatrix;

int
main(int argc,char*argv[])
{{ // use double so that we can destruct objects at end.

  // printf("EMable::minIncrementProbabilty.v = %f,lfm=%f\n",EMable::minIncrementProbabilty.val(),log_FLT_MIN);

  ////////////////////////////////////////////
  // set things up so that if an FP exception
  // occurs such as an "invalid" (NaN), overflow
  // or divide by zero, we actually get a FPE
  ieeeFPsetup();
  set_new_handler(memory_error);

  CODE_TO_COMPUTE_ENDIAN;

  ////////////////////////////////////////////
  // parse arguments
  bool parse_was_ok = Arg::parse(argc,(char**)argv,
"\nThis program learns the optimal model parameters from training data\n");
  if(!parse_was_ok) {
    Arg::usage(); exit(-1);
  }

#define GMTK_ARGUMENTS_CHECK_ARGS
#include "GMTK_Arguments.h"
#undef GMTK_ARGUMENTS_CHECK_ARGS

  gomFS = instantiateFileSource();
  globalObservationMatrix = gomFS;

  /////////////////////////////////////////////
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

  // comment for now Sun Jan 11 09:47:23 2004
  GM_Parms.finalizeParameters();
  GM_Parms.markObjectsToNotTrain(objsToNotTrainFile,cppCommandOptions);

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

  GM_Parms.setStride(gomFS->stride());

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

  //  printf("Dlinks: min lag %d    max lag %d\n", Dlinks::globalMinLag(), Dlinks::globalMaxLag());
  // FIXME - min past = min(dlinkPast, VECPTPast), likewise for future
  int dlinkPast = Dlinks::globalMinLag();
  dlinkPast = (dlinkPast < 0) ? -dlinkPast : 0;
  gomFS->setMinPastFrames( dlinkPast );
  
  int dlinkFuture = Dlinks::globalMaxLag();
  dlinkFuture = (dlinkFuture > 0) ? dlinkFuture : 0;
  gomFS->setMinFutureFrames( dlinkFuture );


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

  Range* pdbrng = new Range(pdbrng_str,0,0x7FFFFFFF);
  myjt.setPartitionDebugRange(*pdbrng);

  if (randomizeParams) {
    infoMsg(IM::Default,"WARNING: GMTK is randomizing all trainable parameters and writing them to file 'random.gmp'\n");
    GM_Parms.makeRandom();
    oDataStreamFile of("random.gmp");
    GM_Parms.writeTrainable(of);
  }

  if (gomFS->numSegments()==0) {
    infoMsg(IM::Default,"ERROR: no segments are available in observation file. Exiting...");
    exit_program_with_status(0);
  }

  Range* trrng = new Range(trrng_str,0,gomFS->numSegments());
#if 0
  if (trrng->length() <= 0) {
    infoMsg(IM::Default,"Training range '%s' specifies empty set. Exiting...\n",
	  trrng_str);
    exit_program_with_status(0);
  }
#endif

  if (trrng->length() == 0 && loadAccFile == NULL) {
    error("ERROR: with EM training. Either must specify segments to train or must load accumulatores (or both).");
  }

  if (island) {
    // island needs to know about this internally.
    myjt.setCurEMTrainingBeam(emTrainingBeam);
  }


  logpr total_data_prob = 1.0;

  struct rusage rus; /* starting time */
  struct rusage rue; /* ending time */
  getrusage(RUSAGE_SELF,&rus);

  /////////////////////////////////////////////////////////
  // first load any and all accumulators
  if (loadAccFile != NULL) {
    if (loadAccRange == NULL) {
      infoMsg(IM::Default,"Loading accumulators from '%s'\n",loadAccFile);
      iDataStreamFile inf(loadAccFile,accFileIsBinary);
      inf.read(total_data_prob.valref());
      GM_Parms.emLoadAccumulators(inf);
    } else {
      Range lfrng(loadAccRange,0,10000);
      for (Range::iterator lfit=lfrng.begin();
	   !lfit.at_end();
	   lfit++) {
	const int bufsize = 2048;
	char buff[bufsize];
	copyStringWithTag(buff,loadAccFile,(*lfit),bufsize);
	iDataStreamFile inf(buff,accFileIsBinary);
	if (lfit == lfrng.begin()) {
	  infoMsg(IM::Default,"Loading accumulators from '%s'\n",buff);
	  inf.read(total_data_prob.valref());
	  GM_Parms.emLoadAccumulators(inf);
	} else {
	  infoMsg(IM::Default,"Accumulating accumulators from '%s'\n",buff);
	  logpr tmp;
	  inf.read(tmp.valref());
	  total_data_prob *= tmp;
	  GM_Parms.emAccumulateAccumulators(inf);
	}
      }
    }
  }

  // Now, do EM training iterations
  logpr previous_dp;
  previous_dp.set_to_almost_zero();
  if (fsize(llStoreFile) == sizeof(logpr)) {
    // get the previous log likelyhood if it exists in a file
    // so we can keep track of how much the likelihood is changing
    // when we do externally driven EM iteration (say during
    // parallel training).
    iDataStreamFile inf(llStoreFile,false);
    inf.read(previous_dp.valref());
  }

  double llDiffPerc = 100.0;
  
  unsigned i=0;
  for (; i<maxEMIterations; i++)  {

    unsigned total_num_frames = 0;

    if (trrng->length() > 0) {
      total_data_prob = 1.0;
      Range::iterator* trrng_it = new Range::iterator(trrng->begin());
      while (!trrng_it->at_end()) {
	const unsigned segment = (unsigned)(*(*trrng_it));
	try {
	  if (gomFS->numSegments() < (segment+1)) 
	    error("ERROR: only %d segments in file, training range must be in range [%d,%d] inclusive\n",
		  gomFS->numSegments(),
		  0,gomFS->numSegments()-1);

	  const unsigned numFrames = GM_Parms.setSegment(segment);
#if 0
	  if (gomFS->active()) {
	    gomFS->printSegmentInfo();
	    ::fflush(stdout);
	  }
#endif


	  if (island) {
	    unsigned numUsableFrames;
	    myjt.collectDistributeIsland(numFrames,
					 numUsableFrames,
					 base,
					 lst,
					 rootBase, islandRootPower,
					 true, // run EM algorithm
					 false, // run Viterbi algorithm
					 localCliqueNormalization);
	    total_num_frames += numUsableFrames;
	    printf("Segment %d, after Island, log(prob(evidence)) = %f, per frame =%f, per numUFrams = %f\n",
		   segment,
		   myjt.curProbEvidenceIsland().val(),
		   myjt.curProbEvidenceIsland().val()/numFrames,
		   myjt.curProbEvidenceIsland().val()/numUsableFrames);
	    if (myjt.curProbEvidenceIsland().not_essentially_zero()) {
	      total_data_prob *= myjt.curProbEvidenceIsland();
	    }
	  } else {
	    unsigned numUsableFrames = myjt.unroll(numFrames);
	    gomFS->justifySegment(numUsableFrames);
	    total_num_frames += numUsableFrames;
	    infoMsg(IM::Low,"Collecting Evidence\n");
	    myjt.collectEvidence();
	    infoMsg(IM::Low,"Done Collecting Evidence\n");
	    logpr probe = myjt.probEvidence();
	    printf("Segment %d, after CE, log(prob(evidence)) = %f, per frame =%f, per numUFrams = %f\n",
		   segment,
		   probe.val(),
		   probe.val()/numFrames,
		   probe.val()/numUsableFrames);
	    if (probe.essentially_zero()) {
	      infoMsg(IM::Default,"Not training segment since probability is essentially zero\n");
	    } else {
	      total_data_prob *= probe;
	      infoMsg(IM::Low,"Distributing Evidence\n");
	      myjt.distributeEvidence();
	      infoMsg(IM::Low,"Done Distributing Evidence\n");
	    
	      if (IM::messageGlb(IM::Huge)) {
		// print out all the clique probabilities. In the ideal
		// case, they should be the same.
		myjt.printProbEvidenceAccordingToAllCliques();
	      }
	      // And actually train with EM.
	      infoMsg(IM::Low,"Incrementing EM Accumulators\n");
	      myjt.emIncrement(probe,localCliqueNormalization,emTrainingBeam);
	    }
	  }
	} catch (ZeroCliqueException &e) {
	  warning("Segment %d aborted due to zero clique\n", segment);
	}
	(*trrng_it)++;
      }
      infoMsg(IM::Default,"EMIter%d: Total data log prob from %d frames processed is: %1.9e\n",
	      i,
	      total_num_frames,total_data_prob.val());
    }


    if (storeAccFile != NULL) {
      // just store the accumulators and exit.
      warning("NOTE: storing current accumulators (from training %d segments) to file '%s' and exiting.",
	      trrng->length(),storeAccFile);
      oDataStreamFile outf(storeAccFile,accFileIsBinary);
      outf.write(total_data_prob.val());
      GM_Parms.emStoreAccumulators(outf);
      exit_program_with_status(0);
    }

    // at this point, either we should have
    // done some training or we should have
    // accumulated something.

    GM_Parms.emEndIteration();
    // if (total_data_prob > previous_dp)
    GM_Parms.emSwapCurAndNew();

    /////////////////////////////////////////////////////////
    // the basic parameters after each iteration 
    if (writeParametersAfterEachEMIteration 
	&& outputTrainableParameters != NULL) {
      char buff[2048];
      copyStringWithTag(buff,outputTrainableParameters,
			i,2048);
      oDataStreamFile of(buff,binOutputTrainableParameters);
      GM_Parms.writeTrainable(of);
    }
    // also write according to output master
    GM_Parms.write(outputMasterFile,cppCommandOptions,i);  

    // store the current total data probability to a file.
    if (llStoreFile != NULL) {
      oDataStreamFile of(llStoreFile,false);
      of.write(total_data_prob.val());
    }

    // compute the log likelihood difference percentage
    if (previous_dp.val() == 0) {
      // this means that the data has probability 1, since log(1) = 0
      if (total_data_prob.val() == previous_dp.val())
	llDiffPerc = 0.0;
      else {
	previous_dp.valref() = ::exp((double)-20.0);
	llDiffPerc = 
	  100.0*fabs((total_data_prob.val() - previous_dp.val())/previous_dp.val());
      }
    } else {
      llDiffPerc = 
	100.0*fabs((total_data_prob.val() - previous_dp.val())/previous_dp.val());
    }
    previous_dp = total_data_prob;

    if (llDiffPerc < lldp) {
      printf("Log likelihood difference percentage (%e) fell below threshold (%e). Ending EM training.\n",llDiffPerc,lldp);
      break;
    }
  }
  if (i == maxEMIterations) {
    printf("Maximum number of EM iterations reached (%d). Stopping\n",maxEMIterations);
  }

  /////////////////////////////////////////////////////////
  // finally, write out the final basic parameters
  // w/o any numeric tag.
  if (outputTrainableParameters != NULL) {
    char buff[2048];
    copyStringWithTag(buff,outputTrainableParameters,
		      CSWT_EMPTY_TAG,2048);
    oDataStreamFile of(buff,binOutputTrainableParameters);
    GM_Parms.writeTrainable(of);
  }
  // also write according to output master
  GM_Parms.write(outputMasterFile,cppCommandOptions);  

  getrusage(RUSAGE_SELF,&rue);
  if (IM::messageGlb(IM::Default)) { 
    infoMsg(IM::Default,"### Final time (seconds) just for EM stage: ");
    double userTime,sysTime;
    reportTiming(rus,rue,userTime,sysTime,stdout);
  }


} // close brace to cause a destruct on valid end of program before print message.
 exit_program_with_status(0); 
}


