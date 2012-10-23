/*
 * gmtkJT.cc
 * produce a junction tree
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
#include "version.h"

#include "GMTK_WordOrganization.h"

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
#include "GMTK_Dlinks.h"
#include "GMTK_DlinkMatrix.h"
#include "GMTK_ProgramDefaultParms.h"
#include "GMTK_BoundaryTriangulate.h"
#include "GMTK_JunctionTree.h"
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

#if 0
static unsigned boostVerbosity=0;
const static char *boostVerbosityRng=NULL;
#endif

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
#if 0
ObservationMatrix globalObservationMatrix;
#endif
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



  // make sure that all observation variables work
  // with the global observation stream.
  infoMsg(IM::Max,"Checking consistency between cpts and observations...\n");
  fp.checkConsistentWithGlobalObservationStream();
  GM_Parms.checkConsistentWithGlobalObservationStream();

  GM_Parms.setStride(gomFS->stride());


  /////
  // TODO: check that beam is a valid value.
  // logpr pruneRatio;
  // pruneRatio.valref() = -beam;

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

  myjt.setCliquePrintRanges(pPartCliquePrintRange,cPartCliquePrintRange,ePartCliquePrintRange);
  infoMsg(IM::Default,"DONE creating Junction Tree\n"); fflush(stdout);
  ////////////////////////////////////////////////////////////////////

  if (gomFS->numSegments()==0)
    error("ERROR: no segments are available in observation file");

  if (IM::messageGlb(IM::Giga)) { 
    gm_template.reportScoreStats();
  }

  Range* dcdrng = new Range(dcdrng_str,0,gomFS->numSegments());
  if (dcdrng->length() <= 0) {
    infoMsg(IM::Default,"Training range '%s' specifies empty set. Exiting...\n",
	  dcdrng_str);
    exit_program_with_status(0);
  }

  
  Range* pdbrng = new Range(pdbrng_str,0,0x7FFFFFFF);
  myjt.setPartitionDebugRange(*pdbrng);


  struct rusage rus; /* starting time */
  struct rusage rue; /* ending time */
  getrusage(RUSAGE_SELF,&rus);


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
      if (probE) {
	unsigned numUsableFrames;
      
	// Range* bvrng = NULL;
	// if (boostVerbosityRng != NULL)
	// bvrng = new Range(bvrng,0,gomFS->numSegments());

	// logpr probe = myjt.probEvidence(numFrames,numUsableFrames,bvrng,boostVerbosity);

	infoMsg(IM::Max,"Beginning call to probability of evidence.\n");
	logpr probe = myjt.probEvidenceFixedUnroll(numFrames,&numUsableFrames);
	printf("Segment %d, after Prob E: log(prob(evidence)) = %f, per frame =%f, per numUFrams = %f\n",
	       segment,
	       probe.val(),
	       probe.val()/numFrames,
	       probe.val()/numUsableFrames);
      } else if (island) {
	unsigned numUsableFrames;

	infoMsg(IM::Max,"Beginning call to island collect/distribute evidence.\n");
	logpr probe = myjt.collectDistributeIsland(numFrames,
						   numUsableFrames,
						   base,
						   lst,
						   rootBase, islandRootPower);

	printf("Segment %d, after island Prob E: log(prob(evidence)) = %f, per frame =%f, per numUFrams = %f\n",
	       segment,
	       probe.val(),
	       probe.val()/numFrames,
	       probe.val()/numUsableFrames);

      } else {

	infoMsg(IM::Max,"Beginning call to unroll\n");
	unsigned numUsableFrames = myjt.unroll(numFrames);
	gomFS->justifySegment(numUsableFrames);

	infoMsg(IM::Low,"Collecting Evidence\n");
	myjt.collectEvidence();
	infoMsg(IM::Low,"Done Collecting Evidence\n");
	logpr probe = myjt.probEvidence();
	printf("Segment %d, after CE, log(prob(evidence)) = %f, per frame =%f, per numUFrams = %f\n",
	       segment,
	       probe.val(),
	       probe.val()/numFrames,
	       probe.val()/numUsableFrames);

	if (doDistributeEvidence) {
	  infoMsg(IM::Low,"Distributing Evidence\n");
	  myjt.distributeEvidence();
	  infoMsg(IM::Low,"Done Distributing Evidence\n");

	  if (JunctionTree::viterbiScore)
	    infoMsg(IM::SoftWarning,"NOTE: Clique sums will be different since viteri option is active\n");
	  if (IM::messageGlb(IM::Low)) {
	    myjt.printProbEvidenceAccordingToAllCliques();
	    probe = myjt.probEvidence();
	    printf("Segment %d, after DE, log(prob(evidence)) = %f, per frame =%f, per numUFrams = %f\n",
		   segment,
		   probe.val(),
		   probe.val()/numFrames,
		   probe.val()/numUsableFrames);
	  }
	}
	if (pPartCliquePrintRange || cPartCliquePrintRange || ePartCliquePrintRange)
	  myjt.printAllCliques(stdout,true,cliquePrintOnlyEntropy);
      }
    } catch (ZeroCliqueException &e) {
      warning("Segment %d aborted due to zero clique\n", segment);
    }
    (*dcdrng_it)++;
  }
  
  getrusage(RUSAGE_SELF,&rue);
  if (IM::messageGlb(IM::Default)) { 
    infoMsg(IM::Default,"### Final time (seconds) just for inference: ");
    double userTime,sysTime;
    reportTiming(rus,rue,userTime,sysTime,stdout);
  }

} // close brace to cause a destruct on valid end of program.
 exit_program_with_status(0); 
}

