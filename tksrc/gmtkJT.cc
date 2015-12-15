/*
 * gmtkJT.cc
 *   compute probability of evidence or clique posteriors
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *            Richard Rogers <rprogers@uw.edu>
 *
 * Copyright (C) 2001, 2015 Jeff Bilmes
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

#include "GMTK_ObservationSource.h"
#include "GMTK_FileSource.h"
#include "GMTK_CreateFileSource.h"
#include "GMTK_ASCIIFile.h"
#include "GMTK_FlatASCIIFile.h"
#include "GMTK_PFileFile.h"
#include "GMTK_HTKFile.h"
#include "GMTK_HDF5File.h"
#include "GMTK_BinaryFile.h"
#include "GMTK_Filter.h"
#include "GMTK_Stream.h"
#include "GMTK_ZeroCliqueException.h"

// Supported inference tasks
#include "GMTK_ProbEvidenceTask.h"
#include "GMTK_ForwardBackwardTask.h"

// Section scheduling
#include "GMTK_LinearSectionScheduler.h"
#include "GMTK_IslandSectionScheduler.h"
#include "GMTK_ArchipelagosSectionScheduler.h"

// Supported within-sectin inference algorithms
#include "GMTK_SparseJoinInference.h"
#include "GMTK_PedagogicalInference.h"

#include "GMTK_MixtureCommon.h"
#include "GMTK_GaussianComponent.h"
#include "GMTK_MeanVector.h"
#include "GMTK_DiagCovarVector.h"
#include "GMTK_Dlinks.h"
#include "GMTK_DlinkMatrix.h"
#include "GMTK_ProgramDefaultParms.h"
#include "GMTK_BoundaryTriangulate.h"
#include "GMTK_MaxClique.h"


VCID(HGID)



/*****************************   OBSERVATION INPUT FILE HANDLING   **********************************************/
#define GMTK_ARG_OBS_FILES

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

/*************************          MEMORY MANAGEMENT OPTIONS         *******************************************/
#define GMTK_ARG_MEMORY_MANAGEMENT_OPTIONS
#define GMTK_ARG_HASH_LOAD_FACTOR
#define GMTK_ARG_STORE_DETERMINISTIC_CHILDREN
#define GMTK_ARG_CLEAR_CLIQUE_VAL_MEM
#define GMTK_ARG_MEM_GROWTH
#define GMTK_ARG_USE_MMAP

/****************************      FILE RANGE OPTIONS             ***********************************************/
#define GMTK_ARG_FILE_RANGE_OPTIONS
#define GMTK_ARG_DCDRNG
#define GMTK_ARG_START_END_SKIP

/****************************         GENERAL OPTIONS             ***********************************************/
#define GMTK_ARG_GENERAL_OPTIONS
#define GMTK_ARG_SEED
#define GMTK_ARG_SKIP_STARTUP_CHECKS
#define GMTK_ARG_VERB
#define GMTK_ARG_HELP
#define GMTK_ARG_VERSION
#define GMTK_ARG_CLIQUE_PRINT


/****************************         INFERENCE OPTIONS           ***********************************************/
#define GMTK_ARG_INFERENCE_OPTIONS
#define GMTK_ARG_DO_DIST_EVIDENCE
#define GMTK_ARG_PROB_EVIDENCE
#define GMTK_ARG_ONLY_KEEP_SEPS
#define GMTK_ARG_ISLAND
#define GMTK_ARG_DEBUG_PART_RNG
#define GMTK_ARG_DEBUG_INCREMENT
#define GMTK_ARG_CLIQUE_TABLE_NORMALIZE
#define GMTK_ARG_CE_SEP_DRIVEN
#define GMTK_ARG_MIXTURE_CACHE
#define GMTK_ARG_VITERBI_SCORE
#define GMTK_ARG_CLIQUE_VAR_ITER_ORDERS
#define GMTK_ARG_JT_OPTIONS
#define GMTK_ARG_VE_SEPS
#define GMTK_ARG_FAIL_ON_ZERO_CLIQUE

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
FileSource *gomFS;
ObservationSource *globalObservationMatrix;


int
main(int argc,char*argv[])
{{ // use double so that we can destruct objects at end.

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

  infoMsg(IM::Max,"Opening Files ...\n");
  gomFS = instantiateFileSource();
  globalObservationMatrix = gomFS;

  infoMsg(IM::Max,"Finished opening files.\n");


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

  // make sure that all observation variables work with the global observation stream.
  infoMsg(IM::Max,"Checking consistency between cpts and observations...\n");
  fp.checkConsistentWithGlobalObservationStream();
  GM_Parms.checkConsistentWithGlobalObservationStream();

  GM_Parms.setStride(gomFS->stride());


  /////
  // TODO: check that beam is a valid value.
  // logpr pruneRatio;
  // pruneRatio.valref() = -beam;


  //  printf("Dlinks: min lag %d    max lag %d\n", Dlinks::globalMinLag(), Dlinks::globalMaxLag());
  // FIXME - min past = min(dlinkPast, VECPTPast), likewise for future
  int dlinkPast = Dlinks::globalMinLag();
  dlinkPast = (dlinkPast < 0) ? -dlinkPast : 0;
  gomFS->setMinPastFrames( dlinkPast );
  
  int dlinkFuture = Dlinks::globalMaxLag();
  dlinkFuture = (dlinkFuture > 0) ? dlinkFuture : 0;
  gomFS->setMinFutureFrames( dlinkFuture );

  if (gomFS->numSegments()==0)
    error("ERROR: no segments are available in observation file");

  Range* dcdrng = new Range(dcdrng_str,0,gomFS->numSegments());
  if (dcdrng->length() <= 0) {
    infoMsg(IM::Default,"Training range '%s' specifies empty set. Exiting...\n",
	  dcdrng_str);
    exit_program_with_status(0);
  }


  infoMsg(IM::Max,"Creating template...\n");
  GMTemplate gm_template(fp);


  // Instantiate the requested inference algorithm for the time series as a whole

  SectionScheduler *section_scheduler = NULL;

  if (false) {
#ifdef GMTK_ISLANDSECTIONSCHEDULER_H
  } else if (island) {
    section_scheduler = new IslandSectionScheduler(gm_template, fp, gomFS);
#endif
#ifdef GMTK_ARCHIPELAGOSSECTIONSCHEDULER_H
  } else if (false /*archipelagos*/) {
    section_scheduler = new ArchipelagosSectionScheduler(gm_template, fp, gomFS);
#endif
#ifdef GMTK_LINEARSECTIONSCHEDULER_H
  } else {
    section_scheduler = new LinearSectionScheduler(gm_template, fp, gomFS);
#endif
  }


  // Setup the within-section inference implementation

  // FIXME - move these to arguments
  bool pedagogical = false, sparse_join = true;

  SectionInferenceAlgorithm *section_inference_alg = NULL;
  if (false) {
#ifdef GMTK_PEDAGOGICALINFERENCE_H
  } else if (pedagogical) {
    section_inference_alg = new PedagogicalInference(section_scheduler);
#endif
#ifdef GMTK_LOOPYBELIEFINFERENCE_H
  } else if (loopy) {
    section_inference_alg = new LoopyBeliefInference(section_scheduler);
#endif
#ifdef GMTK_SPARSEJOININFERENCE_H
  } else if (sparse_join) {
    section_inference_alg = new SparseJoinInference(section_scheduler); // current "standard" algorithm
#endif
  }
  assert(section_inference_alg);


  string tri_file;
  if (triFileName == NULL) 
    tri_file = string(strFileName) + GMTemplate::fileExtension;
  else 
    tri_file = string(triFileName);

  {
    // do this in scope so that is gets deleted now rather than later.
    iDataStreamFile is(tri_file.c_str());
    section_scheduler->setUpDataStructures(is, varPartitionAssignmentPrior,varCliqueAssignmentPrior, checkTriFileCards);
  }

  ForwardBackwardTask *fwd_bkwd_alg = NULL;
  ProbEvidenceTask    *probE_alg    = NULL;

  if (!probE || doDistributeEvidence) {
    // doing the forward/backward task
    fwd_bkwd_alg = dynamic_cast<ForwardBackwardTask *>(section_scheduler);
    assert(fwd_bkwd_alg); // The selected section inference algorithm must implement the ForwardBackwardTask API.
                          // The argument checking logic must prevent illegal combinations of inference task & algorithm,
                          // e.g., OnlineInference can't do the ForwardBackwardTask because it can't seek backwards
                          // arbitrarily far in the observation stream.
  } else {
    // doing the probability of evidence task (forward pass only)
    probE_alg = dynamic_cast<ProbEvidenceTask *>(section_scheduler);
    assert(probE_alg); // The selected section inference algorithm must implement the ProbEvidenceTask API.
                       // The argument checking logic must prevent illegal combinations of inference task & algorithm,
                       // e.g., island can't do ProbEvidenceTask because it by definition does a backward pass.
  }


  if (jtFileName != NULL)
    section_scheduler->printInferencePlanSummary(jtFileName);


  if (IM::messageGlb(IM::Giga)) { 
    section_scheduler->reportScoreStats();
  }

  section_scheduler->setCliquePrintRanges(pSectionCliquePrintRange, cSectionCliquePrintRange, eSectionCliquePrintRange);

  // setup enhanced verbosity for selected sections
  Range* pdbrng = new Range(pdbrng_str,0,0x7FFFFFFF);
  section_scheduler->setSectionDebugRange(*pdbrng);

  // Output file in one of the supported observation file formats to hold the clique posteriors
  ObservationFile *clique_posterior_file = NULL; 
  // Only use an observation file format if the user specified a file name and selected some cliques.
  // Otherwise, output (if any) goes to stdout in ASCII format
  if (cliqueOutputName && (pSectionCliquePrintRange || cSectionCliquePrintRange || eSectionCliquePrintRange) ) {

      unsigned p_size, c_size, e_size;
      section_scheduler->getCliquePosteriorSize(p_size, c_size, e_size);
      unsigned clique_size = (p_size > c_size) ? p_size : c_size;
      clique_size = (clique_size > e_size) ? clique_size : e_size;
      
      // For the sections (P', C', E') that have selected cliques for output, the selected
      // cliques in each section must be the same size.
      if (pSectionCliquePrintRange && p_size != clique_size) {
	error("ERROR: incompatible prologue cliques selected for file output: selected P cliques are size %u, other clique size %u\n", p_size, clique_size);
      }
      if (cSectionCliquePrintRange && c_size != clique_size) {
	error("ERROR: incompatible chunk cliques selected for file output: selected C cliques are size %u, other clique size %u\n", c_size, clique_size);
      }
      if (eSectionCliquePrintRange && e_size != clique_size) {
	error("ERROR: incompatible epilogue cliques selected for file output: selected E cliques are size %u, other clique %u\n", e_size, clique_size);
      }
      section_scheduler->printCliqueOrders(stdout);
      clique_posterior_file = instantiateWriteFile(cliqueListName, cliqueOutputName, cliquePrintSeparator,
						   cliquePrintFormat, clique_size, 0, cliquePrintSwap);

  }

  // trac inference time
  struct rusage rus; /* starting time */
  struct rusage rue; /* ending time */
  getrusage(RUSAGE_SELF,&rus);

  // iterate over selected portions of observation data performing requested inference
  Range::iterator* dcdrng_it = new Range::iterator(dcdrng->begin());
  while (!dcdrng_it->at_end()) {
    const unsigned segment = (unsigned)(*(*dcdrng_it));
    if (gomFS->numSegments() < (segment+1)) 
      error("ERROR: only %d segments in file, segment must be in range [%d,%d]\n",
	    gomFS->numSegments(),
	    0,gomFS->numSegments()-1);

    infoMsg(IM::Max,"Loading segment %d ...\n",segment);
    const unsigned numFrames = GM_Parms.setSegment(segment);
    infoMsg(IM::Max,"Finished loading segment %d with %d frames.\n",segment,numFrames);
    
    try {
      unsigned numUsableFrames;
      logpr probe;

      // TODO:  maybe  if (fwd_bkwd_alg) { ... }  if (probE_alg) { ... }

      if (!probE || doDistributeEvidence) { // doing the forward/backward task
	assert(fwd_bkwd_alg);

	probe = fwd_bkwd_alg->forwardBackward(section_inference_alg,
					      &numUsableFrames,
					      cliquePosteriorNormalize, 
					      cliquePosteriorUnlog,
					      clique_posterior_file,
					      doDistributeEvidence);
	
      } else {  // doing the probability of evidence task (forward pass only)

	infoMsg(IM::Max,"Beginning call to probability of evidence.\n"); // TODO: move to probEvidence() impl
	assert(probE_alg);

	probe = probE_alg->probEvidence(section_inference_alg,
					&numUsableFrames,
					NULL, // returns # of modified sections used for the current segment
					false,  // impose a time limit
					false,  // skip inference on E'
					cliquePosteriorNormalize,
					cliquePosteriorUnlog,
					clique_posterior_file);
	// FIXME - move this to probEvidence() impl ?
	printf("Segment %d, after Prob E: log(prob(evidence)) = %f, per frame =%f, per numUFrams = %f\n",
	       segment,
	       probe.val(),
	       probe.val()/numFrames,
	       probe.val()/numUsableFrames);

      }
    } catch (ZeroCliqueException &e) {
      warning("Segment %d aborted due to zero clique\n", segment);
    }

    // if we're writing clique posteriors to a file, terminate the current output segment
    if (clique_posterior_file) {
      clique_posterior_file->endOfSegment();
    }

    (*dcdrng_it)++; // increment current segment
  }
  
  delete clique_posterior_file; // close the clique posterior output file

  getrusage(RUSAGE_SELF,&rue);
  if (IM::messageGlb(IM::Default)) { 
    infoMsg(IM::Default,"### Final time (seconds) just for inference: ");
    double userTime,sysTime;
    reportTiming(rus,rue,userTime,sysTime,stdout);
  }

} // close brace to cause a destruct on valid end of program.
 exit_program_with_status(0); 
}

