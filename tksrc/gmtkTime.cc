/*
 * gmtkTime.cc
 *   run inference for a given fixed amount of absolute time and report back the amount of work that was done in that time.
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
#include <signal.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/wait.h>
#include <unistd.h>


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
#define GMTK_ARG_INPUT_TRAINABLE_PARAMS
#define GMTK_ARG_ALLOC_DENSE_CPTS
#define GMTK_ARG_CPT_NORM_THRES

/*************************   INPUT STRUCTURE PARAMETER FILE HANDLING  *******************************************/
#define GMTK_ARG_STR_FILE
#define GMTK_ARG_TRI_FILE
#define GMTK_ARG_CHECK_TRI_FILE_CARD
#define GMTK_ARG_JT_INFO_FILE
#define GMTK_ARG_JT_INFO_FILE_DEF_VAL NULL
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
#define GMTK_ARG_VERB_DEF_VAL (IM::Default-1)
#define GMTK_ARG_HELP
#define GMTK_ARG_VERSION

/****************************         INFERENCE OPTIONS           ***********************************************/
#define GMTK_ARG_DO_DIST_EVIDENCE
#define GMTK_ARG_PROB_EVIDENCE
#define GMTK_ARG_ISLAND
#define GMTK_ARG_CLIQUE_TABLE_NORMALIZE
#define GMTK_ARG_CE_SEP_DRIVEN
#define GMTK_ARG_MIXTURE_CACHE
#define GMTK_ARG_CLIQUE_VAR_ITER_ORDERS
#define GMTK_ARG_JT_OPTIONS
#define GMTK_ARG_VE_SEPS

/************************  OBSERVATION MATRIX TRANSFORMATION OPTIONS   ******************************************/
#define GMTK_ARG_OBS_MATRIX_XFORMATION

/************************            TIMING OPTIONS                    ******************************************/
#define GMTK_ARG_TIMING

// should be made conditional on having setrlimit available
#define GMTK_ARG_RLIMIT_PARAMS

/////////////////////////////////////////////////////////////
// General Options



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

/*
 *  Signal handler to set JunctionTree's probEvidenceTime expired timer.
 */
void jtExpiredSigHandler(int arg) {
  JunctionTree::probEvidenceTimeExpired = true;
}

#define MAX(x,y) ((x)>(y)?(x):(y))

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
    Arg::usage(); 
    exit(EXIT_FAILURE);
  }

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
  GM_Parms.finalizeParameters();

  /////////////////////////////
  // read in the structure of the GM, this will
  // die if the file does not exist.
  FileParser fp(strFileName,cppCommandOptions);
  infoMsg(IM::Tiny,"Finished reading in all parameters and structures\n");
  infoMsg(IM::Tiny,"Finished reading in structure\n");

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

  /////
  // TODO: check that beam is a valid value.
  // logpr pruneRatio;
  // pruneRatio.valref() = -beam;
  if (!multiTest) {

    // we're not in multitest mode.
    string tri_file;
    if (triFileName == NULL) 
      tri_file = string(strFileName) + GMTemplate::fileExtension;
    else 
      tri_file = string(triFileName);

    GMTemplate gm_template(fp);
    {
      // do this in scope so that is gets deleted now rather than later.
      iDataStreamFile is(tri_file.c_str(),false,false);
      if (!fp.readAndVerifyGMId(is,checkTriFileCards))
	error("ERROR: triangulation file '%s' does not match graph given in structure file '%s'\n",tri_file.c_str(),strFileName);
    
      gm_template.readPartitions(is);
      gm_template.readMaxCliques(is);
    }
    gm_template.triangulatePartitionsByCliqueCompletion();
    { 
      // check that graph is indeed triangulated.
      BoundaryTriangulate triangulator(fp,
				       gm_template.maxNumChunksInBoundary(),
				       gm_template.chunkSkip(),1.0);
      if (!triangulator.ensurePartitionsAreChordal(gm_template)) {
	error("ERROR: triangulation file '%s' is not chordal",
	      tri_file.c_str());
      }
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
    
    if (globalObservationMatrix.numSegments()==0)
      error("ERROR: no segments are available in observation file");
    
    Range* dcdrng = new Range(dcdrng_str,0,globalObservationMatrix.numSegments());
    if (dcdrng->length() <= 0) {
      infoMsg(IM::Default,"Training range '%s' specifies empty set. Exiting...\n",
	      dcdrng_str);
      exit_program_with_status(0);
    }


    JunctionTree::probEvidenceTimeExpired = false;
    signal(SIGALRM,jtExpiredSigHandler);

    unsigned curTime;
    for (curTime = 0; curTime < numTimes; curTime ++) {

      JunctionTree::probEvidenceTimeExpired = false;
      printf("Run %d: Running program for approximately %d seconds\n",curTime,seconds);
      fflush(stdout);

      alarm(seconds);
      struct rusage rus; /* starting time */
      struct rusage rue; /* ending time */
      getrusage(RUSAGE_SELF,&rus);

      unsigned totalNumberPartitionsDone = 0;
      unsigned totalNumberSegmentsDone = 0;
      unsigned numCurPartitionsDone = 0;
      while (1) {
	Range::iterator* dcdrng_it = new Range::iterator(dcdrng->begin());
	while (!dcdrng_it->at_end()) {
	  const unsigned segment = (unsigned)(*(*dcdrng_it));
	  if (globalObservationMatrix.numSegments() < (segment+1)) 
	    error("ERROR: only %d segments in file, segment must be in range [%d,%d]\n",
		  globalObservationMatrix.numSegments(),
		  0,globalObservationMatrix.numSegments()-1);

	  const unsigned numFrames = GM_Parms.setSegment(segment);

	  unsigned numUsableFrames;
	  numCurPartitionsDone = 0;
	  if (probE && !island) {
	    logpr probe = myjt.probEvidenceTime(numFrames,numUsableFrames,numCurPartitionsDone,noEPartition);
	    totalNumberPartitionsDone += numCurPartitionsDone;
	    infoMsg(IM::Info,"Segment %d, after Prob E: log(prob(evidence)) = %f, per frame =%f, per numUFrams = %f\n",
		    segment,
		    probe.val(),
		    probe.val()/numFrames,
		    probe.val()/numUsableFrames);
	  } else if (island) {
	    myjt.collectDistributeIsland(numFrames,
					 numUsableFrames,
					 base,
					 lst);
	    // TODO: note that frames not always equal to partitions but
	    // do this for now. Ultimately fix this.
	    totalNumberPartitionsDone += numUsableFrames;
	  } else {
	    error("gmtkTime doesn't currently support linear full-mem collect/distribute evidence. Use either '-probE' option, or as a simulation, '-island -lst HUGE_INT'\n");
	  }
      
	  if (JunctionTree::probEvidenceTimeExpired)
	    break;

	  (*dcdrng_it)++;
	  totalNumberSegmentsDone ++;
	}
	delete dcdrng_it;
	if (JunctionTree::probEvidenceTimeExpired)
	  break;
      }
  
      getrusage(RUSAGE_SELF,&rue);
      alarm(0);
      double userTime,sysTime;
      reportTiming(rus,rue,userTime,sysTime,stdout);
      printf("Inference stats: %0.2f seconds, %d segments + %d residual partitions, %d total partitions, %0.3e partitions/sec\n",
	     userTime,
	     totalNumberSegmentsDone,
	     numCurPartitionsDone,
	     totalNumberPartitionsDone,
	     (double)totalNumberPartitionsDone/userTime);
    }

  } else {
    // we're in multitest mode.

    printf("Running in multi-test mode\n");
    fflush(stdout);

    unsigned iteration = 0;
    bool first = true;

    // TODO: these ultimately should come from the same place.
    string orig_vpap_str = varPartitionAssignmentPrior;
    string orig_vcap_str = varCliqueAssignmentPrior;
    string orig_jcap_str = JunctionTree::junctionTreeMSTpriorityStr;
    string orig_icap_str = JunctionTree::interfaceCliquePriorityStr;
    string orig_tri_file;
    unsigned orig_seconds = seconds;
    bool orig_comp_cache = MixtureCommon::cacheMixtureProbabilities;

    if (triFileName == NULL) 
      orig_tri_file = string(strFileName) + GMTemplate::fileExtension;
    else 
      orig_tri_file = string(triFileName);

    double bestRate = 0.0;
    string best_tri_file;
    string best_vpap_str;
    string best_vcap_str;
    string best_jcap_str;
    string best_icap_str;
    bool best_cc = MixtureCommon::cacheMixtureProbabilities;


    while (1) {

      // Utilize both the partition information and elimination order
      // information already computed and contained in the file. This
      // enables the program to use external triangulation programs,
      // where this program ensures that the result is triangulated
      // and where it reports the quality of the triangulation.

      // restore original options to overrided only if need be.
      string tri_file = orig_tri_file;
      string vpap_str  = orig_vpap_str ;
      string vcap_str =  orig_vcap_str;
      string jcap_str =  orig_jcap_str;
      string icap_str =  orig_icap_str;
      seconds = orig_seconds;
      MixtureCommon::cacheMixtureProbabilities = orig_comp_cache;

      // get name of triangulation file (and other options) from the command line.
      //   If trifile is named 'end' then, stop processing.
      //   If an empty line occurs, then use triFileName approach, assuming trifile is re-written.

      char buff[16384];
      fflush(stdout);
      if (!fgets(buff,sizeof(buff),stdin)) {
	// we're done.
	break;
      }
      unsigned n = strlen(buff);
      if (buff[n-1] == '\n')
	buff[n-1] = '\0';
      if (strlen(buff) == 0) {
	// do nothing since defaults are already set.
      } else if (strcmp(buff,"END") == 0 || 		 
		 strcmp(buff,"end") == 0 ||
		 strcmp(buff,"quit") == 0 ||
		 strcmp(buff,".") == 0) {
	multiTest = false;
	break; // out of enclosing do loop 
      } else {
	// TODO: support vcap, etc. options. vcap=DBD, etc.
	// parse the buffer string using the following format.
	// option1="value" option2="value2" option3="value3" etc.
	// double quotes can be excaped (part of the value) by using \" character.
	// I.e., we can do trifile="foobar\"baz" will ask for a file named foobar"baz
	// Options currently supported:
	//   trifile="file"
	//   vcap="vcap options"
	//   vpap="vpap options"
	//   jcap="jpap options"
	//   icap="icap options"

	char* buffp=buff;
	while (*buffp) {
	  // get next option.

	  // skip white space
	  while ( *buffp && isspace(*buffp) )
	    buffp++;

	  // end if this is the end of the string.
	  if (!*buffp)
	    break;

	  // TODO: do this in a better way.
	  if (!strncmp("vcap=\"",buffp,6)) {
	    buffp += 6; // skip to option.
	    char* buffpp = buffp+1;
	    while (*buffpp && ((*buffpp != '"') || (buffpp[-1] == '\\'))) {
	      buffpp++;
	    }
	    if (!*buffpp)
	      break; // missing end quote, so no option.
	    
	    *buffpp = '\0';
	    vcap_str = buffp;
	    // printf("vcap_str=(%s)\n",vcap_str.c_str());
	    *buffpp = '"';
	    buffp = buffpp+1;
	  } else if (!strncmp("vpap=\"",buffp,6)) {
	    buffp += 6; // skip to option.
	    char* buffpp = buffp+1;
	    while (*buffpp && ((*buffpp != '"') || (buffpp[-1] == '\\'))) {
	      buffpp++;
	    }
	    if (!*buffpp)
	      break; // missing end quote, so no option.
	    
	    *buffpp = '\0';
	    vpap_str = buffp;
	    *buffpp = '"';
	    buffp = buffpp+1;
	  } else if (!strncmp("jcap=\"",buffp,6)) {
	    buffp += 6; // skip to option.
	    char* buffpp = buffp+1;
	    while (*buffpp && ((*buffpp != '"') || (buffpp[-1] == '\\'))) {
	      buffpp++;
	    }
	    if (!*buffpp)
	      break; // missing end quote, so no option.
	    
	    *buffpp = '\0';
	    jcap_str = buffp;
	    *buffpp = '"';
	    buffp = buffpp+1;
	  } else if (!strncmp("icap=\"",buffp,6)) {
	    buffp += 6; // skip to option.
	    char* buffpp = buffp+1;
	    while (*buffpp && ((*buffpp != '"') || (buffpp[-1] == '\\'))) {
	      buffpp++;
	    }
	    if (!*buffpp)
	      break; // missing end quote, so no option.
	    
	    *buffpp = '\0';
	    icap_str = buffp;
	    *buffpp = '"';
	    buffp = buffpp+1;
	  } else if (!strncmp("trifile=\"",buffp,9)) {
	    buffp += 9; // skip to option.
	    char* buffpp = buffp+1;
	    while (*buffpp && ((*buffpp != '"') || (buffpp[-1] == '\\'))) {
	      buffpp++;
	    }
	    if (!*buffpp)
	      break; // missing end quote, so no option.
	    
	    *buffpp = '\0';
	    tri_file = buffp;
	    // printf("trifile=(%s)\n",tri_file.c_str());
	    *buffpp = '"';
	    buffp = buffpp+1;
	  } else if (!strncmp("seconds=\"",buffp,9)) {
	    buffp += 9; // skip to option.
	    char* buffpp = buffp+1;
	    while (*buffpp && ((*buffpp != '"') || (buffpp[-1] == '\\'))) {
	      buffpp++;
	    }
	    if (!*buffpp)
	      break; // missing end quote, so no option.
	    
	    *buffpp = '\0';

	    char* endptr;
	    unsigned val = (unsigned) strtol(buffp, &endptr, 0);
	    if ( endptr == buffp ) {
	      // fail
	      fprintf(stderr,"WARNING: bad option to seconds=\"%s\"\n",buffp);
	      buffp++;
	    } else {
	      seconds = val;
	      buffp = buffpp+1;
	    }
	    *buffpp = '"';
	  } else if (!strncmp("componentCache=\"",buffp,16)) {
	    buffp += 16; // skip to option.
	    char* buffpp = buffp+1;
	    while (*buffpp && ((*buffpp != '"') || (buffpp[-1] == '\\'))) {
	      buffpp++;
	    }
	    if (!*buffpp)
	      break; // missing end quote, so no option.
	    
	    *buffpp = '\0';
	    if (*buffp == 'T' || *buffp == 't') {
	      MixtureCommon::cacheMixtureProbabilities = true;
	      buffp+=2;
	    } else if (*buffp == 'F' || *buffp == 'f') {
	      MixtureCommon::cacheMixtureProbabilities = false;
	      buffp+=2;
	    } else {
	      // fail
	      fprintf(stderr,"WARNING: bad option to componentCache=\"%c\"",*buffp);
	      buffp++;
	    }
	    *buffpp = '"';
	  } else {
	    fprintf(stderr,"WARNING: unrecognized string option (%s)\n",buffp);
	    // skip
	    buffp++;
	  }
	}
      }

      // TODO: clean this up. code below can't retain these ptrs.
      JunctionTree::junctionTreeMSTpriorityStr=(char*)jcap_str.c_str();
      JunctionTree::interfaceCliquePriorityStr=(char*)icap_str.c_str();

      if (first) {
	best_tri_file = tri_file;
	best_vpap_str = vpap_str;
	best_vcap_str =	vcap_str;
	best_jcap_str = jcap_str;
	best_icap_str = icap_str;
	best_cc = MixtureCommon::cacheMixtureProbabilities;
      }


      GMTemplate gm_template(fp);
      {
	// do this in scope so that is gets deleted now rather than later.
	iDataStreamFile is(tri_file.c_str(),false,false);
	if (!fp.readAndVerifyGMId(is,checkTriFileCards))
	  error("ERROR: triangulation file '%s' does not match graph given in structure file '%s'\n",tri_file.c_str(),strFileName);
    
	gm_template.readPartitions(is);
	gm_template.readMaxCliques(is);
      }
      gm_template.triangulatePartitionsByCliqueCompletion();
      { 
	// check that graph is indeed triangulated.
	BoundaryTriangulate triangulator(fp,
					 gm_template.maxNumChunksInBoundary(),
					 gm_template.chunkSkip(),1.0);
	if (!triangulator.ensurePartitionsAreChordal(gm_template)) {
	  error("ERROR: triangulation file '%s' is not chordal",
		tri_file.c_str());
	}
      }


      ////////////////////////////////////////////////////////////////////
      // CREATE JUNCTION TREE DATA STRUCTURES
      infoMsg(IM::Default,"Creating Junction Tree\n"); fflush(stdout);
      JunctionTree myjt(gm_template);
      myjt.setUpDataStructures(vpap_str.c_str(),vcap_str.c_str());
      myjt.prepareForUnrolling();
      if (jtFileName != NULL)
	myjt.printAllJTInfo(jtFileName);
      infoMsg(IM::Default,"DONE creating Junction Tree\n"); fflush(stdout);
      ////////////////////////////////////////////////////////////////////

      if (globalObservationMatrix.numSegments()==0)
	error("ERROR: no segments are available in observation file");

      Range* dcdrng = new Range(dcdrng_str,0,globalObservationMatrix.numSegments());
      if (dcdrng->length() <= 0) {
	infoMsg(IM::Default,"Training range '%s' specifies empty set. Exiting...\n",
		dcdrng_str);
	exit_program_with_status(0);
      }


      JunctionTree::probEvidenceTimeExpired = false;
      signal(SIGALRM,jtExpiredSigHandler);


      struct ipc_struct {
	unsigned totalNumberPartitionsDone;
	unsigned totalNumberSegmentsDone;
	unsigned numCurPartitionsDone ;
	ipc_struct() {
	  totalNumberPartitionsDone = 0;
	  totalNumberSegmentsDone = 0;
	  numCurPartitionsDone  = 0;
	};
      };

      // create a simple pipe for the child to communicate some stuff
      // to the parent.
      int filedes[2];
      if (pipe(&filedes[0])) {
	error("ERROR: can't create pipe. errno = %d, %s\n",errno,strerror(errno));
      }

      int& read_fd = filedes[0];
      int& write_fd = filedes[1];

      const int limitTime = MAX(seconds+rlimitSlop,1);

      printf("--------\n%d: Operating on trifile '%s'\n",iteration,tri_file.c_str());
      printf("%d: Other options: vpap=%s,vcap=%s,jcap=%s,icap=%s,cc=%d\n",
	     iteration,
	     vpap_str.c_str(),vcap_str.c_str(),
	     jcap_str.c_str(),icap_str.c_str(),MixtureCommon::cacheMixtureProbabilities);
      printf("%d: ",iteration); 
      printf("Running program for approximately %d seconds, not to exceed %d CPU seconds.\n",seconds,limitTime);
      fflush(stdout);


      pid_t pid = fork();
      if (pid != 0) {

	// this is the parent process.
	// printf("%d: child process = %d\n",iteration,pid);

	// close, so that we don't keep accumulating open fds. 
	close(write_fd);

	int status;
	struct rusage rus; /* starting time */
	struct rusage rue; /* ending time */
	int rc;
	if ((rc=getrusage(RUSAGE_CHILDREN,&rus)))
	  error("ERROR: parent process can't call getruage start, returned %d\n",rc);

	fflush(stdout);
	// wait for the child.
	waitpid(pid,&status,0);

	if (WEXITSTATUS(status)== EXIT_SUCCESS && !WIFSIGNALED(status)) {

	  // Assume process ended normally.
	  if ((rc=getrusage(RUSAGE_CHILDREN,&rue)))
	    error("ERROR: parent process can't call getrusage end, returned %d\n",rc);

	  double userTime,sysTime;
	  printf("%d: Actual running time: ",iteration);
	  reportTiming(rus,rue,userTime,sysTime,stdout);
	  fflush(stdout);

	  // get parameters from file that child must have written.
	  ipc_struct child_info;
	  if ((rc=read(read_fd,(void*)&child_info,sizeof(child_info))) != sizeof(child_info)) {
	    error("ERROR: can't read from child pipe, errno = %d, %s",errno,strerror(errno));
	  }
	  printf("%d: ",iteration);

	  double curRate;
	  if (userTime > 0.0) 
	    curRate = (double)child_info.totalNumberPartitionsDone/userTime;
	  else
	    curRate = 0.0;
	  printf("Inference stats: %0.2f seconds, %d segments + %d residual partitions, %d total partitions, %0.3e partitions/sec\n",
		 userTime,
		 child_info.totalNumberSegmentsDone,
		 child_info.numCurPartitionsDone,
		 child_info.totalNumberPartitionsDone,
		 curRate);
	  fflush(stdout);
	  
	  if (curRate > bestRate) {
	    best_tri_file = tri_file;
	    best_vpap_str = vpap_str;
	    best_vcap_str = vcap_str;
	    best_jcap_str = jcap_str;
	    best_icap_str = icap_str;
	    bestRate = curRate;
	  }

	} else {
	  // child exited abnormally, probably ran out of time.
	  printf("%d: NOTICE: Triangulation failed to complete in alloted time of %d seconds, or process failed (status = 0x%X): ",
		 iteration,limitTime,status);
	  if ((rc=getrusage(RUSAGE_CHILDREN,&rue)))
	    error("ERROR: parent process can't call getrusage end, returned %d\n",rc);
	  double userTime,sysTime;
	  reportTiming(rus,rue,userTime,sysTime,stdout);
	  fflush(stdout);
	}
	// close down the pipe in any case.
	close(read_fd);

      } else {
	// this is the child process.

	close(read_fd);

	// limit the amount of time we run.
	struct rlimit rlim;
	rlim.rlim_cur = limitTime;
	rlim.rlim_max = limitTime;
	int rc;

	if ((rc = setrlimit(RLIMIT_CPU,&rlim)))
	  warning("WARNING: child process can't set limit to %d seconds, setrlimit() returned %d. No hard limit on process time!!!\n",
		  rlim.rlim_cur,rc);

	alarm(seconds);

	// struct rusage rus; /* starting time */
	//struct rusage rue; /* ending time */
	// getrusage(RUSAGE_SELF,&rus);

	ipc_struct child_info;
	while (1) {
	  Range::iterator* dcdrng_it = new Range::iterator(dcdrng->begin());
	  while (!dcdrng_it->at_end()) {
	    const unsigned segment = (unsigned)(*(*dcdrng_it));
	    if (globalObservationMatrix.numSegments() < (segment+1)) 
	      error("ERROR: only %d segments in file, segment must be in range [%d,%d]\n",
		    globalObservationMatrix.numSegments(),
		    0,globalObservationMatrix.numSegments()-1);

	    const unsigned numFrames = GM_Parms.setSegment(segment);

	    unsigned numUsableFrames;
	    child_info.numCurPartitionsDone = 0;
	    if (probE && !island) {
	      logpr probe = myjt.probEvidenceTime(numFrames,numUsableFrames,child_info.numCurPartitionsDone,noEPartition);
	      child_info.totalNumberPartitionsDone += child_info.numCurPartitionsDone;
	      infoMsg(IM::Info,"Segment %d, after Prob E: log(prob(evidence)) = %f, per frame =%f, per numUFrams = %f\n",
		      segment,
		      probe.val(),
		      probe.val()/numFrames,
		      probe.val()/numUsableFrames);
	    } else if (island) {
	      myjt.collectDistributeIsland(numFrames,
					   numUsableFrames,
					   base,
					   lst);
	      // TODO: note that frames not always equal to partitions but
	      // do this for now. Ultimately fix this.
	      child_info.totalNumberPartitionsDone += numUsableFrames;
	    } else {
	      error("gmtkTime doesn't currently support linear full-mem collect/distribute evidence\n");
	    }
      
	    if (JunctionTree::probEvidenceTimeExpired)
	      break;

	    (*dcdrng_it)++;
	    child_info.totalNumberSegmentsDone ++;
	  }
	  delete dcdrng_it;
	  if (JunctionTree::probEvidenceTimeExpired)
	    break;
	}
	// turn off the signal.
	alarm(0);

	// write stuff to pipe.
	if ((rc=write(write_fd,(void*)&child_info,sizeof(child_info))) != sizeof(child_info)) {
	  error("ERROR: can't write to parent pipe, errno = %d, %s",errno,strerror(errno));
	}

	// gracefully close
	close(write_fd);

	// exit normally, so parent realizes this.
	exit(EXIT_SUCCESS);
	// END OF CHILD PROCESS
      }
      iteration++;
      first = false;
    }

    printf("--------\n");
    printf("Best trifile found at %0.3e partitions/sec is '%s'\n",bestRate,best_tri_file.c_str());
    printf("Best options: vpap=%s, vcap=%s, jcap=%s, icap=%s, cc=%d\n",best_vpap_str.c_str(),best_vcap_str.c_str(),
	   best_jcap_str.c_str(),best_icap_str.c_str(),best_cc);
    printf("--------\n");

  } // end of multi-test section.

  exit_program_with_status(EXIT_SUCCESS);
}
