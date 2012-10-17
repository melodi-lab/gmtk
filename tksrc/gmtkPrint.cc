/*
 * gmtkPrint.cc
 *
 * A GMTK program to print out Viterbi values from files written
 * by gmtkViterbi with the -binaryViterbiFile option
 *
 *
 * Copyright (c) 2012, < fill in later >
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
#include <regex.h>

// TODO: remove next 2 eventually 
#include <iostream>
#include <fstream>

#include "general.h"
#include "error.h"
#include "debug.h"
#include "rand.h"
#include "arguments.h"
#include "ieeeFPsetup.h"
#include "version.h"
#include "file_utils.h"

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
#include "GMTK_Signals.h"

#define GMTK_ARGUMENTS_REQUIRE_BINARY_VIT_FILE

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


/*************************          MEMORY MANAGEMENT OPTIONS         *******************************************/
#define GMTK_ARG_MEMORY_MANAGEMENT_OPTIONS
#define GMTK_ARG_HASH_LOAD_FACTOR
#define GMTK_ARG_STORE_DETERMINISTIC_CHILDREN
#define GMTK_ARG_CLEAR_CLIQUE_VAL_MEM
#define GMTK_ARG_MEM_GROWTH
#define GMTK_ARG_USE_MMAP

/*************************          INFERENCE OPTIONS                 *******************************************/
#define GMTK_ARG_CLIQUE_VAR_ITER_ORDERS
#define GMTK_ARG_DEBUG_PART_RNG

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

/************************  OBSERVATION MATRIX TRANSFORMATION OPTIONS   ******************************************/
#define GMTK_ARG_OBS_MATRIX_OPTIONS
#define GMTK_ARG_OBS_MATRIX_XFORMATION

/************************            DECODING OPTIONS                  ******************************************/
#define GMTK_ARG_NEW_DECODING_OPTIONS

// should be made conditional on having setrlimit available
#define GMTK_ARG_RESOURCE_OPTIONS
#define GMTK_ARG_RLIMIT_PARAMS

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
{

  ////////////////////////////////////////////
  // set things up so that if an FP exception
  // occurs such as an "invalid" (NaN), overflow
  // or divide by zero, we actually get a FPE
  ieeeFPsetup();
  set_new_handler(memory_error);
  InstallSignalHandlers();

  CODE_TO_COMPUTE_ENDIAN;
   
  ////////////////////////////////////////////
  // parse arguments
  bool parse_was_ok = Arg::parse(argc,(char**)argv,
"\nThis program determines the most likely values of the hidden variabls\n");
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

  myjt.setCliquePrintRanges(pPartCliquePrintRange,cPartCliquePrintRange,ePartCliquePrintRange);
  infoMsg(IM::Default,"DONE creating Junction Tree\n"); fflush(stdout);
  ////////////////////////////////////////////////////////////////////


  if (gomFS->numSegments()==0) {
    infoMsg(IM::Default,"ERROR: no segments are available in observation file. Exiting...");
    exit_program_with_status(0);
  }

  Range* dcdrng = new Range(dcdrng_str,0,gomFS->numSegments());
  if (dcdrng->length() == 0) { 
    error("Decoding range must specify a non-zero length range. Range given is %s\n",
	  dcdrng_str);
  }

  logpr total_data_prob = 1.0;

  if (!JunctionTree::binaryViterbiFile) {
    error("Argument Error: Missing REQUIRED argument: -binaryViterbiFile <str>\n");
  }

  char cookie[GMTK_VITERBI_COOKIE_LENGTH+1];

  if (fgets(cookie, GMTK_VITERBI_COOKIE_LENGTH+1, JunctionTree::binaryViterbiFile) != cookie) {
    error("ERROR: Binary viterbi file '%s' did not begin with %s", 
	  JunctionTree::binaryViterbiFilename, GMTK_VITERBI_COOKIE);
  }
  if (strcmp(cookie, GMTK_VITERBI_COOKIE)) {
    error("ERROR: Binary viterbi file '%s' did not begin with %s", 
	  JunctionTree::binaryViterbiFilename, GMTK_VITERBI_COOKIE);
  }
  
  unsigned num_segments_in_file;
  if (fread(&num_segments_in_file, sizeof(num_segments_in_file), 1, JunctionTree::binaryViterbiFile) != 1) {
    error("ERROR: failed to read # of segments from '%s'\n", JunctionTree::binaryViterbiFilename);
  }
  if (num_segments_in_file != gomFS->numSegments()) {
    error("ERROR: '%s' contains %u segments, but the current observation files contain %u\n",
	  JunctionTree::binaryViterbiFilename, num_segments_in_file, gomFS->numSegments());
  }

  unsigned N_best;
  if (fread(&N_best, sizeof(N_best), 1, JunctionTree::binaryViterbiFile) != 1) {
    error("ERROR: failed to read N (for N-best) from '%s'\n", JunctionTree::binaryViterbiFilename);
  }


  // What is the difference between the pVit* files and the vit*
  // files?  The pVit* files are similar to the vit* files, but the
  // pVit print via partitions, while the vit* does a bunch more
  // processing to give the illusion that there are only frames/slices
  // rather than GMTK partitions (which could span any number of
  // frames, even partial frames). Also, the vit files are more meant
  // for users, while the pVit files are more meant for algorithm
  // writters/debugging which is the reason that we have both of them
  // here.
  FILE* pVitValsFile = NULL;
  if (pVitValsFileName) {
    if (strcmp("-",pVitValsFileName) == 0)
      pVitValsFile = stdout;
    else {
      if ((pVitValsFile = fopen(pVitValsFileName, "w")) == NULL)
	error("Can't open file '%s' for writing\n",pVitValsFileName);
    }
  }
#if 1
  FILE* vitValsFile = NULL;
  if (vitValsFileName) {
    if (vitValsFileName && strcmp("-",vitValsFileName) == 0)
      vitValsFile = stdout;
    else {
      if ((vitValsFile = fopen(vitValsFileName, "w")) == NULL)
	error("Can't open file '%s' for writing\n",vitValsFileName);
    }
  }
  if (!pVitValsFile && !vitValsFile) {
    error("Argument Error: Missing REQUIRED argument: -pVitValsFile <str>  OR  -vitValsFile <str>\n");
  }
#endif

  
  Range* pdbrng = new Range(pdbrng_str,0,0x7FFFFFFF);
  myjt.setPartitionDebugRange(*pdbrng);

  // We always do viterbi scoring/option in this program.
  JunctionTree::viterbiScore = true;

  struct rusage rus; /* starting time */
  struct rusage rue; /* ending time */
  getrusage(RUSAGE_SELF,&rus);

  Range::iterator* dcdrng_it = new Range::iterator(dcdrng->begin());
    
  regex_t *pVitPreg = NULL;
  if (pVitRegexFilter != NULL) {
    pVitPreg = (regex_t*) malloc(sizeof(regex_t));
    const unsigned case_ignore = (pVitCaseSensitiveRegexFilter? 0 : REG_ICASE);
    if (regcomp(pVitPreg,pVitRegexFilter,
		REG_EXTENDED
		| case_ignore
		| REG_NOSUB
		)) {
      error("ERROR: problem with regular expression filter string '%s'\n",pVitRegexFilter);
    }
  }

  regex_t *vitPreg = NULL;
  if (vitRegexFilter != NULL) {
    vitPreg = (regex_t*) malloc(sizeof(regex_t));
    const unsigned case_ignore = (vitCaseSensitiveRegexFilter? 0 : REG_ICASE);
    if (regcomp(vitPreg,vitRegexFilter,
		REG_EXTENDED
		| case_ignore
		| REG_NOSUB
		)) {
      error("ERROR: problem with regular expression filter string '%s'\n",vitRegexFilter);
    }
  }

  while (!dcdrng_it->at_end()) {
    const unsigned segment = (unsigned)(*(*dcdrng_it));

    gmtk_off_t indexOff;
    gmtk_off_t off;
    float score;

    indexOff = (gmtk_off_t) ( GMTK_VITERBI_HEADER_SIZE + segment * (sizeof(gmtk_off_t) + sizeof(float)) );
    if (gmtk_fseek(JunctionTree::binaryViterbiFile, indexOff, SEEK_SET)) {
      char *err = strerror(errno);
      error("Error seeking in '%s': %s\n", JunctionTree::binaryViterbiFilename, err);
    }
    if (fread(&off, sizeof(off), 1, JunctionTree::binaryViterbiFile) != 1) {
      char *err = strerror(errno);
      error("Error reading from '%s': %s\n", JunctionTree::binaryViterbiFilename, err);
    }
    JunctionTree::binaryViterbiOffset = off;
    if (fread(&score, sizeof(score), 1, JunctionTree::binaryViterbiFile) != 1) {
      char *err = strerror(errno);
      error("Error reading from '%s': %s\n", JunctionTree::binaryViterbiFilename, err);
    }
    if (gmtk_fseek(JunctionTree::binaryViterbiFile, off, SEEK_SET)) {
      char *err = strerror(errno);
      error("Error seeking in '%s': %s\n", JunctionTree::binaryViterbiFilename, err);
    }
//printf("idx seg %03x -> %04llx @ %04llx    %llx\n", segment, off, indexOff, ftello(JunctionTree::binaryViterbiFile));

    if (gomFS->numSegments() < (segment+1)) 
      error("ERROR: only %d segments in file, decode range must be in range [%d,%d] inclusive\n",
	    gomFS->numSegments(),
	    0,gomFS->numSegments()-1);

    const unsigned numFrames = GM_Parms.setSegment(segment);

    logpr probe(NULL, score);

    total_data_prob *= probe;

    if (pVitValsFile) {
      fprintf(pVitValsFile,"========\nSegment %d, number of frames = %d, viterbi-score = %f\n",
	      segment, numFrames, probe.val());
      myjt.printSavedPartitionViterbiValues(numFrames,
					    JunctionTree::binaryViterbiFile,
					    pVitValsFile,
					    pVitAlsoPrintObservedVariables,
					    pVitPreg,
					    pVitPartRangeFilter);
    }

    if (pVitValsFile || pPartCliquePrintRange || cPartCliquePrintRange || ePartCliquePrintRange)
      myjt.resetViterbiPrinting();
    if (vitValsFile) {
      fprintf(vitValsFile,"========\nSegment %d, number of frames = %d, viterbi-score = %f\n",
	      segment, numFrames, score);
      if (!vitFrameRangeFilter) {
	myjt.printSavedViterbiValues(numFrames, 
				     vitValsFile,
				     JunctionTree::binaryViterbiFile,
				     vitAlsoPrintObservedVariables,
				     vitPreg,
				     vitPartRangeFilter);
      } else {
	myjt.printSavedViterbiFrames(numFrames, vitValsFile, NULL,
				     vitAlsoPrintObservedVariables,
				     vitPreg,
				     vitFrameRangeFilter);
      }
    }
    (*dcdrng_it)++;
  }

  if (pVitPreg != NULL) {
    regfree(pVitPreg);
    free(pVitPreg);
  }

  infoMsg(IM::Default,"Total data log prob for all segments is: %1.9e\n",
	  total_data_prob.val());

  getrusage(RUSAGE_SELF,&rue);
  if (IM::messageGlb(IM::Default)) { 
    infoMsg(IM::Default,"### Final time (seconds) just for decoding stage: ");
    double userTime,sysTime;
    reportTiming(rus,rue,userTime,sysTime,stdout);
  }

  if (pVitValsFile && pVitValsFile != stdout)
    fclose(pVitValsFile);
#if 1
  if (vitValsFile && vitValsFile != stdout)
    fclose(vitValsFile);
#endif
  if (JunctionTree::binaryViterbiFile)
    fclose(JunctionTree::binaryViterbiFile);

  exit_program_with_status(0);
}
