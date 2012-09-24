/*
 * gmtkViterbi.cc
 *
 * A GMTK program for various Viterbi decoding options. I.e., this program can:
 *
 *   1) Write out all the k-best random variable values corresponding
 *      to all (or a subset of) the hidden variables in a graph, and
 *      do this for a number of segments.
 *
 *   2) Do a compressed viterbi output. I.e., write out all variable
 *      values for all times of a child random variable (the query
 *      variable, or the word variable on speech) whenever a parent
 *      variable (a transition variable, of cardinality 2) has a value of 1.
 *
      // TODO: 1) DONE: have a regex to filter which random variables are printed
      //       2) DONE: have a range spec to filter which frames are printed
      //       3) trigger: have a trigger condition (i.e., variable(offset) = value)
      //          to filter which frames are printed, where 'offset'
      //          can be either -1, 0, or 1. (-triggerVariable1 -triggerValue1 -triggerVariable2 -triggerValue2 ...), the values should be int sets.
      //       4) compression: have the option to print only when one of the values
      //          has changed from last time (so, print first C, and last E always).
      //          Note this is not a substitute for 3 as the trigger might be true multiple
      //          frames in a row, while the value stays the same (e.g., transition to the same word)
      //    allow 'ands' and 'ors' of the above. 
      //      5) DONE: the option to always print ints, even if symbol table exists.
 *      
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


#include "GMTK_WordOrganization.h"

VCID(HGID)


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
#include "GMTK_Signals.h"


/*****************************   OBSERVATION INPUT FILE HANDLING   **********************************************/
#define GMTK_ARG_OBS_FILES

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
#define GMTK_ARG_ISLAND
#define GMTK_ARG_DEBUG_PART_RNG
#define GMTK_ARG_DEBUG_INCREMENT
#define GMTK_ARG_CLIQUE_TABLE_NORMALIZE
#define GMTK_ARG_CE_SEP_DRIVEN
#define GMTK_ARG_MIXTURE_CACHE
#define GMTK_ARG_CLIQUE_VAR_ITER_ORDERS
#define GMTK_ARG_JT_OPTIONS
#define GMTK_ARG_VE_SEPS

/************************  OBSERVATION MATRIX TRANSFORMATION OPTIONS   ******************************************/
#define GMTK_ARG_OBS_MATRIX_OPTIONS
#define GMTK_ARG_OBS_MATRIX_XFORMATION

/************************            DECODING OPTIONS                  ******************************************/
#define GMTK_VITERBI_FILE_WRITE
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
  InstallSignalHandlers();

  CODE_TO_COMPUTE_ENDIAN;
   
  ////////////////////////////////////////////
  // parse arguments
  bool parse_was_ok = Arg::parse(argc,(char**)argv,
"\nThis program determines the most likely values of the hidden variabls\n");
  if(!parse_was_ok) {
    Arg::usage(); exit(-1);
  }
  
  // if (island == false) {
  // error("Program must currently be run with -island T\n");
  // }

#define GMTK_ARGUMENTS_CHECK_ARGS
#include "GMTK_Arguments.h"
#undef GMTK_ARGUMENTS_CHECK_ARGS

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

  myjt.setCliquePrintRanges(pPartCliquePrintRange,cPartCliquePrintRange,ePartCliquePrintRange);
  infoMsg(IM::Default,"DONE creating Junction Tree\n"); fflush(stdout);
  ////////////////////////////////////////////////////////////////////


  if (globalObservationMatrix.numSegments()==0) {
    infoMsg(IM::Default,"ERROR: no segments are available in observation file. Exiting...");
    exit_program_with_status(0);
  }

  Range* dcdrng = new Range(dcdrng_str,0,globalObservationMatrix.numSegments());
  if (dcdrng->length() == 0) { 
    error("Decoding range must specify a non-zero length range. Range given is %s\n",
	  dcdrng_str);
  }

  logpr total_data_prob = 1.0;

  if (JunctionTree::binaryViterbiFile && (pVitValsFileName || vitValsFileName)) {
    error("Can't do both binary Viterbi output and -pVitValsFileName or -vitValsFileName");
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
  if (!pVitValsFile && !vitValsFile && !JunctionTree::binaryViterbiFile) {
    error("Argument Error: Missing REQUIRED argument: -pVitValsFile <str>  OR  -vitValsFile <str> OR -binaryViterbiFile <str>\n");
  }
#endif

  
  Range* pdbrng = new Range(pdbrng_str,0,0x7FFFFFFF);
  myjt.setPartitionDebugRange(*pdbrng);

  // We always do viterbi scoring/option in this program.
  JunctionTree::viterbiScore = true;

  struct rusage rus; /* starting time */
  struct rusage rue; /* ending time */
  getrusage(RUSAGE_SELF,&rus);

  total_data_prob = 1.0;
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

  if (JunctionTree::binaryViterbiFile) {
//printf("Writing %x segments\n", numSegments);
    if (fputs(GMTK_VITERBI_COOKIE, JunctionTree::binaryViterbiFile) == EOF) {
      char *err = strerror(errno);
      error("Error writing to '%s': %s\n", JunctionTree::binaryViterbiFilename, err);
    }
    unsigned numSegments = globalObservationMatrix.numSegments();
    if (fwrite(&numSegments, sizeof(numSegments), 1, JunctionTree::binaryViterbiFile) != 1) {
      char *err = strerror(errno);
      error("Error writing to '%s': %s\n", JunctionTree::binaryViterbiFilename, err);
    }
    unsigned N_best = 1;
    if (fwrite(&N_best, sizeof(N_best), 1, JunctionTree::binaryViterbiFile) != 1) {
      char *err = strerror(errno);
      error("Error writing to '%s': %s\n", JunctionTree::binaryViterbiFilename, err);
    }
    off_t off = (off_t) 0;
    float score = 0.0/0.0; // initially nan so I can check that the index is written correctly (non-nan)
    for (unsigned i=0; i < numSegments; i+=1) {
      if (fwrite(&off, sizeof(off), 1, JunctionTree::binaryViterbiFile) != 1) {
	char *err = strerror(errno);
	error("Error writing to '%s': %s\n", JunctionTree::binaryViterbiFilename, err);
      }
      if (fwrite(&score, sizeof(score), 1, JunctionTree::binaryViterbiFile) != 1) {
	char *err = strerror(errno);
	error("Error writing to '%s': %s\n", JunctionTree::binaryViterbiFilename, err);
      }
    }
    JunctionTree::nextViterbiOffset = ftello(JunctionTree::binaryViterbiFile);
  }
//printf("seg 0 starts @ %llx\n", ftello(JunctionTree::binaryViterbiFile));

  while (!dcdrng_it->at_end()) {
    const unsigned segment = (unsigned)(*(*dcdrng_it));

    off_t indexOff;
    off_t off;
    float score;
    if (JunctionTree::binaryViterbiFile) {
      off = JunctionTree::nextViterbiOffset;
      indexOff = (off_t) ( GMTK_VITERBI_HEADER_SIZE + segment * (sizeof(off_t) + sizeof(float)) );
      if (fseeko(JunctionTree::binaryViterbiFile, indexOff, SEEK_SET)) {
	char *err = strerror(errno);
	error("Error seeking in '%s': %s\n", JunctionTree::binaryViterbiFilename, err);
      }
//printf("idx seg %03x -> %04llx @ %04llx = %llx\n", segment, off, indexOff, ftello(JunctionTree::binaryViterbiFile));
      if (fwrite(&off, sizeof(off), 1, JunctionTree::binaryViterbiFile) != 1) {
	char *err = strerror(errno);
	error("Error writing to '%s': %s\n", JunctionTree::binaryViterbiFilename, err);
      }
      if (fseeko(JunctionTree::binaryViterbiFile, off, SEEK_SET)) {
	char *err = strerror(errno);
	error("Error seeking in '%s': %s\n", JunctionTree::binaryViterbiFilename, err);
      }
      JunctionTree::binaryViterbiOffset = off;
    }

    if (globalObservationMatrix.numSegments() < (segment+1)) 
      error("ERROR: only %d segments in file, decode range must be in range [%d,%d] inclusive\n",
	    globalObservationMatrix.numSegments(),
	    0,globalObservationMatrix.numSegments()-1);

    const unsigned numFrames = GM_Parms.setSegment(segment);

    logpr probe;
    if (island) {
      unsigned numUsableFrames;
      myjt.collectDistributeIsland(numFrames,
				   numUsableFrames,
				   base,
				   lst,
				   sqrtBase,
				   false, // run EM algorithm
				   true,  // run viterbi algorithm
				   false  // localCliqueNormalization, unused here.
				   );
      probe = myjt.curProbEvidenceIsland();
      printf("Segment %d, after Island, viterbi log(prob(evidence)) = %f, per frame =%f, per numUFrams = %f\n",
	     segment,
	     probe.val(),
	     probe.val()/numFrames,
	     probe.val()/numUsableFrames);
      if (probe.not_essentially_zero()) {
	total_data_prob *= probe;
      }
    } else {
      // linear space inference
      unsigned numUsableFrames = myjt.unroll(numFrames);
      infoMsg(IM::Inference, IM::Med,"Collecting Evidence\n");
      myjt.collectEvidence();
      infoMsg(IM::Inference, IM::Med,"Done Collecting Evidence\n");
      probe = myjt.probEvidence();
      infoMsg(IM::Default,"Segment %d, after CE, viterbi log(prob(evidence)) = %f, per frame =%f, per numUFrams = %f\n",
	     segment,
	     probe.val(),
	     probe.val()/numFrames,
	     probe.val()/numUsableFrames);
      if (probe.essentially_zero()) {
	infoMsg(IM::Default,"Skipping segment %d since probability is essentially zero\n",
		segment);
      } else {
	myjt.setRootToMaxCliqueValue();
	total_data_prob *= probe;
	infoMsg(IM::Inference, IM::Low,"Distributing Evidence\n");
	myjt.distributeEvidence();
	infoMsg(IM::Inference, IM::Low,"Done Distributing Evidence\n");
      }
    }

    if (JunctionTree::binaryViterbiFile) {
      indexOff = (off_t) ( GMTK_VITERBI_HEADER_SIZE + sizeof(off_t) + segment * (sizeof(off_t) + sizeof(float)) );
      if (fseeko(JunctionTree::binaryViterbiFile, indexOff, SEEK_SET)) {
	char *err = strerror(errno);
	error("Error seeking in '%s': %s\n", JunctionTree::binaryViterbiFilename, err);
      }
      score = probe.val();
      if (fwrite(&score, sizeof(score), 1, JunctionTree::binaryViterbiFile) != 1) {
	char *err = strerror(errno);
	error("Error writing to '%s': %s\n", JunctionTree::binaryViterbiFilename, err);
      }
    }
    
    if (probe.essentially_zero())
      warning("Segment %d: Not printing Viterbi values since segment has zero probability\n",
	      segment);
    else {
      if (pPartCliquePrintRange || cPartCliquePrintRange || ePartCliquePrintRange)
	myjt.printAllCliques(stdout,true,cliquePrintOnlyEntropy);

      if (pVitValsFile) {
	fprintf(pVitValsFile,"========\nSegment %d, number of frames = %d, viterbi-score = %f\n",
		segment,numFrames,probe.val());
	myjt.printSavedPartitionViterbiValues(pVitValsFile,
					      pVitAlsoPrintObservedVariables,
					      pVitPreg,
					      pVitPartRangeFilter);
      }

#if 1     
      if (pVitValsFile || pPartCliquePrintRange || cPartCliquePrintRange || ePartCliquePrintRange)
	myjt.resetViterbiPrinting();
      if (vitValsFile) {
	fprintf(vitValsFile,"========\nSegment %d, number of frames = %d, viterbi-score = %f\n",
		segment,numFrames,probe.val());
	myjt.printSavedViterbiValues(numFrames, vitValsFile, NULL,
				     vitAlsoPrintObservedVariables,
				     vitPreg,
				     vitPartRangeFilter);
      }
#endif

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
