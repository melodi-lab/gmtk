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
#include "spi.h"
#include "version.h"

VCID("$Header$");

#include "GMTK_FileParser.h"
#include "GMTK_RandomVariable.h"
#include "GMTK_DiscreteRandomVariable.h"
#include "GMTK_ContinuousRandomVariable.h"
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

/*
 * command line arguments
 */
static bool seedme = false;
static char *strFileName=NULL;
static char *triFileName=NULL;
static double varFloor = GMTK_DEFAULT_VARIANCE_FLOOR;

static int allocateDenseCpts=0;
static char *cppCommandOptions = NULL;

static unsigned verbosity = IM::Default;

static unsigned unroll_k = 0;
int startSkip = 0;
int endSkip = 0;
float lldp = 0.001;
float mnlldp = 0.01;
float beam=-LZERO;

bool doDistributeEvidence=false;

// char *outputTrainableParameters="outParms%d.gmp";
char *outputTrainableParameters=NULL;
bool binOutputTrainableParameters=false;
bool writeParametersAfterEachEMIteration=true;


char *inputMasterFile=NULL;
char *outputMasterFile=NULL;
char *inputTrainableParameters=NULL;
bool binInputTrainableParameters=false;
char *objsToNotTrainFile=NULL;


#define MAX_NUM_OBS_FILES (3)
char *ofs[MAX_NUM_OBS_FILES] = { NULL, NULL, NULL }; 
unsigned nfs[MAX_NUM_OBS_FILES] = { 0, 0, 0 };
unsigned nis[MAX_NUM_OBS_FILES] = { 0, 0, 0 };
char *frs[MAX_NUM_OBS_FILES] = { "all", "all", "all" };
char *irs[MAX_NUM_OBS_FILES] = { "all", "all", "all" };
char *fmts[MAX_NUM_OBS_FILES] = { "pfile", "pfile", "pfile" };
bool iswps[MAX_NUM_OBS_FILES] = { false, false, false };

static unsigned segment=0;

bool print_version_and_exit = false;
Arg Arg::Args[] = {

  /////////////////////////////////////////////////////////////
  // observation input file handling

  Arg("of1",Arg::Req,ofs[0],"Observation File 1"),
  Arg("nf1",Arg::Opt,nfs[0],"Number of floats in observation file 1"),
  Arg("ni1",Arg::Opt,nis[0],"Number of ints in observation file 1"),
  Arg("fr1",Arg::Opt,frs[0],"Float range for observation file 1"),
  Arg("ir1",Arg::Opt,irs[0],"Int range for observation file 1"),
  Arg("fmt1",Arg::Opt,fmts[0],"Format (htk,bin,asc,pfile) for observation file 1"),
  Arg("iswp1",Arg::Opt,iswps[0],"Endian swap condition for observation file 1"),


  Arg("of2",Arg::Opt,ofs[1],"Observation File 1"),
  Arg("nf2",Arg::Opt,nfs[1],"Number of floats in observation file 1"),
  Arg("ni2",Arg::Opt,nis[1],"Number of ints in observation file 1"),
  Arg("fr2",Arg::Opt,frs[1],"Float range for observation file 1"),
  Arg("ir2",Arg::Opt,irs[1],"Int range for observation file 1"),
  Arg("fmt2",Arg::Opt,fmts[1],"Format (htk,bin,asc,pfile) for observation file 1"),
  Arg("iswp2",Arg::Opt,iswps[1],"Endian swap condition for observation file 1"),

  Arg("doDistributeEvidence",Arg::Opt,doDistributeEvidence,"Do distribute evidence also"),

  Arg("segment",Arg::Opt,segment,"Which segment to do"),

  /////////////////////////////////////////////////////////////
  // input parameter/structure file handling

  Arg("cppCommandOptions",Arg::Opt,cppCommandOptions,"Additional CPP command line"),

  Arg("inputMasterFile",Arg::Req,inputMasterFile,"Input file of multi-level master CPP processed GM input parameters"),
  Arg("outputMasterFile",Arg::Opt,outputMasterFile,"Output file to place master CPP processed GM output parameters"),

  Arg("inputTrainableParameters",Arg::Opt,inputTrainableParameters,"File of only and all trainable parameters"),
  Arg("binInputTrainableParameters",Arg::Opt,binInputTrainableParameters,"Binary condition of trainable parameters file"),

  Arg("outputTrainableParameters",Arg::Opt,outputTrainableParameters,"File to place only and all trainable output parametes"),
  Arg("binOutputTrainableParameters",Arg::Opt,binOutputTrainableParameters,"Binary condition of output trainable parameters?"),


  Arg("strFile",Arg::Req,strFileName,"Graphical Model Structure File"),
  Arg("triFile",Arg::Opt,triFileName,"Triangulation file for strFile"),

  /////////////////////////////////////////////////////////////
  // General Options

  Arg("allocateDenseCpts",Arg::Opt,allocateDenseCpts,"Automatically allocate any undefined CPTs. arg = -1, no read params, arg = 0 noallocate, arg = 1 means use random initial CPT values. arg = 2, use uniform values"),

  Arg("seed",Arg::Opt,seedme,"Seed the random number generator"),
  Arg("verbosity",Arg::Opt,verbosity,"Verbosity (0 <= v <= 100) of informational/debugging msgs"),

  // temp arg
  Arg("k",Arg::Opt,unroll_k,"amount to unroll"),

  Arg("version",Arg::Opt,print_version_and_exit,"Print GMTK version number and exit."),

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

  ////////////////////////////////////////////
  // parse arguments
  Arg::parse(argc,argv);
  (void) IM::setGlbMsgLevel(verbosity);

  if (print_version_and_exit) {
    printf("%s\n",gmtk_version_id);
    exit(0);
  }


  ////////////////////////////////////////////
  // check for valid argument values.
  int nfiles = 0;
  unsigned ifmts[MAX_NUM_OBS_FILES];
  for (int i=0;i<MAX_NUM_OBS_FILES;i++) {
    if (ofs[i] != NULL && nfs[i] == 0 && nis[i] == 0)
      error("ERROR: command line parameters must specify one of nf%d and ni%d as not zero",
	    i+1,i+1);
    nfiles += (ofs[i] != NULL);
    if (strcmp(fmts[i],"htk") == 0)
      ifmts[i] = HTK;
    else if (strcmp(fmts[i],"binary") == 0)
      ifmts[i] = RAWBIN;
    else if (strcmp(fmts[i],"ascii") == 0)
      ifmts[i] = RAWASC;
    else if (strcmp(fmts[i],"pfile") == 0)
      ifmts[i] = PFILE;
    else
      error("ERROR: Unknown observation file format type: '%s'\n",fmts[i]);
  }

  globalObservationMatrix.openFiles(nfiles,
				    (const char**)&ofs,
				    (const char**)&frs,
				    (const char**)&irs,
				    (unsigned*)&nfs,
				    (unsigned*)&nis,
				    (unsigned*)&ifmts,
				    (bool*)&iswps,
				    startSkip,
				    endSkip);



  MixtureCommon::checkForValidRatioValues();
  MeanVector::checkForValidValues();
  DiagCovarVector::checkForValidValues();
  DlinkMatrix::checkForValidValues();
  if (lldp < 0.0 || mnlldp < 0.0)
    error("lldp & mnlldp must be >= 0");
  if (beam < 0.0)
    error("beam must be >= 0");
  if (startSkip < 0 || endSkip < 0)
    error("startSkip/endSkip must be >= 0");


  ////////////////////////////////////////////
  // set global variables/change global state from args
  GaussianComponent::setVarianceFloor(varFloor);
  if (seedme)
    rnd.seed();

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
  GM_Parms.loadGlobal();
  // comment for now Sun Jan 11 09:47:23 2004
  // GM_Parms.markObjectsToNotTrain(objsToNotTrainFile,cppCommandOptions);

  /////////////////////////////
  // read in the structure of the GM, this will
  // die if the file does not exist.
  FileParser fp(strFileName,cppCommandOptions);
  printf("Finished reading in all parameters and structures\n");
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
  logpr pruneRatio;
  pruneRatio.valref() = -beam;

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
  iDataStreamFile is(tri_file.c_str());
  if (!fp.readAndVerifyGMId(is))
    error("ERROR: triangulation file '%s' does not match graph given in structure file '%s'\n",tri_file.c_str(),strFileName);
  
  gm_template.readPartitions(is);
  gm_template.readMaxCliques(is);
  gm_template.triangulatePartitionsByCliqueCompletion();
  if (1) { 
    // check that it is triangulated for now, 
    // 
    // TODO: ultimately take this check out so that inference code
    // does not need to link to the triangulation code (either that,
    // or put the triangulation check in a different file, so that we
    // only link to tri check code).
    BoundaryTriangulate triangulator(fp,
				     gm_template.maxNumChunksInBoundary(),
				     gm_template.chunkSkip(),1.0);
    triangulator.ensurePartitionsAreChordal(gm_template);
  }

  printf("Creating Junction Tree\n"); fflush(stdout);
  JunctionTree myjt(gm_template);
  myjt.setUpDataStructures();

  // this will cause problems when printing number of bits.
  // myjt.printAllJTInfo("jt_info_before_unrolling.txt");
  myjt.prepareForUnrolling();
  // TODO: allow user input file name
  myjt.printAllJTInfo("jt_info.txt");
  printf("DONE creating Junction Tree\n"); fflush(stdout);

  // fprintf(stderr,"starting unrolling\n");
  // myjt.unroll(unroll_k);
  // fprintf(stderr,"ending unrolling\n");

  if (globalObservationMatrix.numSegments()==0)
    error("ERROR: no segments are available in observation file");

  if (globalObservationMatrix.numSegments() < (segment+1)) 
    error("ERROR: only %d segments in file, segment must be in range [%d,%d]\n",
	  globalObservationMatrix.numSegments(),
	  0,globalObservationMatrix.numSegments()-1);

  GM_Parms.clampFirstExample();
  for (unsigned i=0;i<segment;i++) {
    GM_Parms.clampNextExample();
  }
  
  globalObservationMatrix.loadSegment(segment);

  int frames = globalObservationMatrix.numFrames();

  // myjt.unroll(partition_unroll_amount);
  unsigned numUsableFrames = myjt.unroll(frames);
  printf("Collecting Evidence\n"); fflush(stdout);
  myjt.collectEvidence();
  printf("Done Collecting Evidence\n"); fflush(stdout);
  // myjt.printAllCliquesProbEvidence();

  logpr probe = myjt.probEvidence();
  printf("after CE, log(prob(evidence)) = %f, per frame =%f, per numUFrams = %f\n",
	 probe.val(),
	 probe.val()/frames,
	 probe.val()/numUsableFrames);

  if (doDistributeEvidence) {
    printf("Distributing Evidence\n");
    myjt.distributeEvidence();
    printf("DONE Distributing Evidence\n");
    myjt.printAllCliquesProbEvidence();
  
    probe = myjt.probEvidence();
    printf("after DE, log(prob(evidence)) = %f, per frame =%f, per numUFrams = %f\n",
	   probe.val(),
	   probe.val()/frames,
	   probe.val()/numUsableFrames);
  }

  
#if 0
  fprintf(stderr,"Hit return when ready:");
  {
    scanf("\n");
  }
#endif

  exit_program_with_status(0);
}
