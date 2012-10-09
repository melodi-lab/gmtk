
/*
 * gmtkOnline
 * dynamic graphical model filtering and smoothing
 *
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2012, < fill in later >
 * 
 * < License reference >
 * < Disclaimer >
 *
 */

#if HAVE_CONFIG_H
#include <config.h>
#endif
#if HAVE_HG_H
#include "hgstamp.h"
#endif


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
#include "GMTK_Dlinks.h"

#include "GMTK_Filter.h"

#include "GMTK_ObservationSource.h"
#include "GMTK_StreamSource.h"
#include "GMTK_ASCIIStream.h"
#include "GMTK_BinStream.h"
#if 0
#include "GMTK_FilterStream.h"
#endif

#include "GMTK_Stream.h"

#include "GMTK_MixtureCommon.h"
#include "GMTK_GaussianComponent.h"
#include "GMTK_MeanVector.h"
#include "GMTK_DiagCovarVector.h"
#include "GMTK_DlinkMatrix.h"
#include "GMTK_ProgramDefaultParms.h"
#include "GMTK_BoundaryTriangulate.h"
#include "GMTK_JunctionTree.h"
#include "GMTK_MaxClique.h"

VCID(HGID)

/*************************   INPUT OBSERVATION STREAM HANDLING  *******************************************/
#define GMTK_ARG_STREAM_INPUT
/*************************   INPUT TRAINABLE PARAMETER FILE HANDLING  *******************************************/
#define GMTK_ARG_INPUT_TRAINABLE_FILE_HANDLING
#define GMTK_ARG_CPP_CMD_OPTS
#define GMTK_ARG_INPUT_MASTER_FILE
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

#define GMTK_ARG_START_END_SKIP
#if 0
/****************************      FILE RANGE OPTIONS             ***********************************************/
#define GMTK_ARG_FILE_RANGE_OPTIONS
#define GMTK_ARG_DCDRNG
#endif


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
#define GMTK_ARG_DEBUG_PART_RNG
#define GMTK_ARG_CLIQUE_TABLE_NORMALIZE
#define GMTK_ARG_CE_SEP_DRIVEN
#define GMTK_ARG_MIXTURE_CACHE
#define GMTK_ARG_VITERBI_SCORE
#define GMTK_ARG_CLIQUE_VAR_ITER_ORDERS
#define GMTK_ARG_JT_OPTIONS
#define GMTK_ARG_VE_SEPS
#if 0
#define GMTK_ARG_ISLAND
#endif

#if 0
/************************  OBSERVATION MATRIX TRANSFORMATION OPTIONS   ******************************************/
#define GMTK_ARG_OBS_MATRIX_OPTIONS
#define GMTK_ARG_OBS_MATRIX_XFORMATION
#endif

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

#if 1
ObservationStream *stream = NULL;
StreamSource *gomSS = NULL;
ObservationSource *globalObservationMatrix = NULL;
#endif

int
main(int argc,char*argv[])
{{ // use double so that we can destruct objects at end.

  ////////////////////////////////////////////
  // set things up so that if an FP exception
  // occurs such as an "invalid" (NaN), overflow
  // or divide by zero, we actually get a FPE
  ieeeFPsetup();
  set_new_handler(memory_error);

  //  CODE_TO_COMPUTE_ENDIAN;


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

    FILE *inFile;

    if (strcmp("-", oss[0])) {
      inFile = fopen(oss[0], (ifmts[0] == RAWBIN) ? "rb" : "r");
    } else {
      inFile = stdin;
    }

    if (!inFile) {
      error("ERROR: '%s' %s", oss[0], strerror(errno));
    }
    
    // FIXME - use StreamSource & support posttrans ?

    // FIXME - add -sfr and -sir for feature range selection ?

    if (ifmts[0] == RAWBIN) {
      stream = new BinaryStream(inFile, nfs[0], nis[0], NULL /* sfr */ , NULL /* sir */, inputNetByteOrder[0]);
    } else if (ifmts[0] == RAWASC) {
      stream = new  ASCIIStream(inFile, nfs[0], nis[0], NULL /* sfr */ , NULL /* sir */);
    } else {
      error("ERROR: -fmt1 must be 'binary' or 'ascii', got '%s'", fmts[0]);
    }
    assert(stream);
    gomSS = new StreamSource(1, &stream, streamBufferSize);
    globalObservationMatrix = gomSS;


  /////////////////////////////////////////////
  // read in all the parameters

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

#if 0

  // make sure that all observation variables work
  // with the global observation stream.
  infoMsg(IM::Max,"Checking consistency between cpts and observations...\n");
  fp.checkConsistentWithGlobalObservationStream();
  GM_Parms.checkConsistentWithGlobalObservationStream();

  GM_Parms.setStride(gomSS->stride());

#endif

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

  if (IM::messageGlb(IM::Giga)) { 
    gm_template.reportScoreStats();
  }

#if 0
  Range* dcdrng = new Range(dcdrng_str,0,gomSS->numSegments());
  if (dcdrng->length() <= 0) {
    infoMsg(IM::Default,"Training range '%s' specifies empty set. Exiting...\n",
	  dcdrng_str);
    exit_program_with_status(0);
  }
#endif
  
  Range* pdbrng = new Range(pdbrng_str,0,0x7FFFFFFF);
  myjt.setPartitionDebugRange(*pdbrng);


  struct rusage rus; /* starting time */
  struct rusage rue; /* ending time */
  getrusage(RUSAGE_SELF,&rus);

  printf("Dlink lags: min %d   max %d\n", Dlinks::globalMinLag(), Dlinks::globalMaxLag());

  unsigned segNum = 0;
  unsigned frameNum;
  // compute min # of frames
  // startSkip? + max(frames(P'),frames(C'),pastFrames(dlinks)) + 
  //              \tau * frames(C') + 
  // endSkip?   + max(frames(C'),frames(E'),future(dlinks))
  // do P'
  for (; !gomSS->EOS(); ) {
    gomSS->preloadFrames(1); // FIXME - min # frames
    for (frameNum = 0; gomSS->numFrames() == 0 || 
	               frameNum < gomSS->numFrames(); 
	 // FIXME - frameNum < seg length - frames(E')
	 frameNum += 1)
    {
      gomSS->enqueueFrames(1); // frames(C')
    }
    // do E'
    printf("Seg %u  %u frames  %u\n", segNum++, frameNum, gomSS->numFrames());
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

