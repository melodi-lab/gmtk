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
//#include "spi.h"
#include "version.h"

#include "GMTK_WordOrganization.h"

VCID("$Header$");

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

/*
 * command line arguments
 */

/////////////////////////////////////////////////////////////
// observation input file handling
#define MAX_NUM_OBS_FILES (5)
char *ofs[MAX_NUM_OBS_FILES] = { NULL, NULL, NULL, NULL,NULL }; 
unsigned nfs[MAX_NUM_OBS_FILES] = { 0, 0, 0,0,0 };
unsigned nis[MAX_NUM_OBS_FILES] = { 0, 0, 0,0,0 };
char *frs[MAX_NUM_OBS_FILES] = { "all", "all", "all","all","all" };
char *irs[MAX_NUM_OBS_FILES] = { "all", "all", "all","all","all" };
char *fmts[MAX_NUM_OBS_FILES] = { "pfile", "pfile", "pfile","pfile","pfile" };
bool iswps[MAX_NUM_OBS_FILES] = { false, false, false,false,false };


/////////////////////////////////////////////////////////////
// input parameter/structure file handling
static char *cppCommandOptions = NULL;
static char *inputMasterFile=NULL;
static char *outputMasterFile=NULL;
static char *inputTrainableParameters=NULL;
static bool binInputTrainableParameters=false;
// static char *outputTrainableParameters="outParms%d.gmp";
static char *outputTrainableParameters=NULL;
static bool binOutputTrainableParameters=false;
static bool writeParametersAfterEachEMIteration=true;
static int allocateDenseCpts=0;


/////////////////////////////////////////////////////////////
// Structure file, Triangulation File, and Junction Tree Options.
static char *strFileName=NULL;
static char *triFileName=NULL;
static char *jtFileName="jt_info.txt";

/////////////////////////////////////////////////////////////
// Continuous RV Options
static double varFloor = GMTK_DEFAULT_VARIANCE_FLOOR;

/////////////////////////////////////////////////////////////
// Beam Options
static double emTrainingBeam=-LZERO;

/////////////////////////////////////////////////////////////
// File Range Options
static char *trrng_str="all";
static int startSkip = 0;
static int endSkip = 0;

/////////////////////////////////////////////////////////////
// General Options
static bool seedme = false;
static unsigned verbosity = IM::Default;
static bool help =false;
static bool print_version_and_exit = false;

/////////////////////////////////////////////////////////////
// Inference Options
static bool island=false;
static unsigned base=3;
static unsigned lst=100;
static char* varPartitionAssignmentPrior = "COI";
static char* varCliqueAssignmentPrior = "COI";



/////////////////////////////////////////////////////////////
// EM Training Options
static unsigned maxEMIterations;
static bool randomizeParams = false;
static float lldp = 0.001;
static float mnlldp = 0.01;
static char *loadAccFile = NULL;
static char *loadAccRange = NULL;
static char *storeAccFile = NULL;
static bool accFileIsBinary = true;
static char *llStoreFile = NULL;
static char *objsToNotTrainFile=NULL;
static bool localCliqueNormalization = false;

//////////////////////////////////////////////////////////////
// Observation matrix options

bool     Cpp_If_Ascii        = false;
char     *pr_str[MAX_NUM_OBS_FILES] = {NULL,NULL,NULL,NULL,NULL};   // per stream per sentence range string 

char*    Action_If_Diff_Num_Frames_Str[MAX_NUM_OBS_FILES]={"er","er","er","er","er"};   // 
unsigned Action_If_Diff_Num_Frames[MAX_NUM_OBS_FILES]={ERROR,ERROR,ERROR,ERROR,ERROR};   // 
char*    Action_If_Diff_Num_Sents_Str[MAX_NUM_OBS_FILES]={"te","te","te","te","te"}; 
unsigned Action_If_Diff_Num_Sents[MAX_NUM_OBS_FILES]={TRUNCATE_FROM_END,TRUNCATE_FROM_END,TRUNCATE_FROM_END,TRUNCATE_FROM_END,TRUNCATE_FROM_END};   // 

char    *Per_Stream_Transforms[MAX_NUM_OBS_FILES]={NULL,NULL,NULL,NULL,NULL};   // 
char    *Post_Transforms=NULL;

char    *Ftr_Combo_Str="none";
unsigned Ftr_Combo=FTROP_NONE;
 
#ifdef INTV_WORDS_BIGENDIAN
bool iswp[MAX_NUM_OBS_FILES] = {true,true,true,true,true};
#else
bool iswp[MAX_NUM_OBS_FILES] = {false,false,false,false,false};
#endif 


Arg Arg::Args[] = {

  /////////////////////////////////////////////////////////////
  // observation input file handling
  Arg("of",Arg::Req,ofs,"Observation File.  Replace X with the file number.",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("nf",Arg::Opt,nfs,"Number of floats in observation file X",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("ni",Arg::Opt,nis,"Number of ints in observation file X",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("fr",Arg::Opt,frs,"Float range for observation file X",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("ir",Arg::Opt,irs,"Int range for observation file X",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("fmt",Arg::Opt,fmts,"Format (htk,binary,ascii,pfile) for observation file X",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("iswp",Arg::Opt,iswps,"Endian swap condition for observation file X",Arg::ARRAY,MAX_NUM_OBS_FILES),

  /////////////////////////////////////////////////////////////
  // input parameter/structure file handling
  Arg("cppCommandOptions",Arg::Opt,cppCommandOptions,"Additional CPP command line"),
  Arg("inputMasterFile",Arg::Req,inputMasterFile,"Input file of multi-level master CPP processed GM input parameters"),
  Arg("outputMasterFile",Arg::Opt,outputMasterFile,"Output file to place master CPP processed GM output parameters"),
  Arg("inputTrainableParameters",Arg::Opt,inputTrainableParameters,"File of only and all trainable parameters"),
  Arg("binInputTrainableParameters",Arg::Opt,binInputTrainableParameters,"Binary condition of trainable parameters file"),
  Arg("outputTrainableParameters",Arg::Opt,outputTrainableParameters,"File to place only and all trainable output parametes"),
  Arg("binOutputTrainableParameters",Arg::Opt,binOutputTrainableParameters,"Binary condition of output trainable parameters?"),
  Arg("wpaeei",Arg::Opt,writeParametersAfterEachEMIteration,"Write Parameters After Each EM Iteration Completes"),
  Arg("allocateDenseCpts",Arg::Opt,allocateDenseCpts,"Auto allocate undef CPTs. arg = -1, no read params, arg = 0 noallocate, arg = 1 use random initial CPT values. arg = 2, use uniform values"),
  Arg("cptNormThreshold",Arg::Opt,CPT::normalizationThreshold,"Read error if |Sum-1.0|/card > norm_threshold"),

  /////////////////////////////////////////////////////////////
  // Structure file, Triangulation File, and Junction Tree Options.
  Arg("strFile",Arg::Req,strFileName,"Graphical Model Structure File"),
  Arg("triFile",Arg::Opt,triFileName,"Triangulation file for strFile"),
  Arg("jtFile",Arg::Opt,jtFileName,"Name of file to write junction tree information"),
  // this is here only to affect printing of jt_info.txt
  Arg("jtwUB",
      Arg::Opt,JunctionTree::jtWeightUpperBound,
      "True means jtWeight is allways an upper bound on true JT weight, false means jtWeight is estimate"),

  /////////////////////////////////////////////////////////////
  // Continuous RV Options
  Arg("varFloor",Arg::Opt,varFloor,"Variance Floor"),
  Arg("floorVarOnRead",Arg::Opt,DiagCovarVector::floorVariancesWhenReadIn,
       "Floor the variances to varFloor when they are read in"),

  /////////////////////////////////////////////////////////////
  // Beam Options
  Arg("cbeam",Arg::Opt,MaxClique::cliqueBeam,"Clique beam width pruning log value"),
  Arg("ckbeam",Arg::Opt,MaxClique::cliqueBeamMaxNumStates,"Prune to this clique max state space (0 = no pruning)"),
  Arg("crbeam",Arg::Opt,MaxClique::cliqueBeamRetainFraction,"Fraction of clique state space to retain. Range: 0 < v <= 1. v = 1 means no pruning"),
  Arg("sbeam",Arg::Opt,SeparatorClique::separatorBeam,"Separator beam width pruning log value"),
  Arg("ebeam",Arg::Opt,emTrainingBeam,"EM training beam width"),

  /////////////////////////////////////////////////////////////
  // Memory management options
  Arg("clearCliqueValMem",Arg::Opt,MaxClique::perSegmentClearCliqueValueCache,"Free clique/separator value cache for each segment"),


  /////////////////////////////////////////////////////////////
  // File Range Options
  Arg("trrng",Arg::Opt,trrng_str,"Range to decode over segment file"),
  Arg("startSkip",Arg::Opt,startSkip,"Frames to skip at beginning (i.e., first frame is buff[startSkip])"),
  Arg("endSkip",Arg::Opt,endSkip,"Frames to skip at end (i.e., last frame is buff[len-1-endSkip])"),

  /////////////////////////////////////////////////////////////
  // General Options
  Arg("seed",Arg::Opt,seedme,"Seed the random number generator"),
  Arg("verbosity",Arg::Opt,verbosity,"Verbosity (0 <= v <= 100) of informational/debugging msgs"),
  Arg("help",   Arg::Tog, help,  "Print this message"),
  Arg("version",Arg::Opt,print_version_and_exit,"Print GMTK version number and exit."),

  /////////////////////////////////////////////////////////////
  // Inference Options
  Arg("island",Arg::Opt,island,"Run island algorithm"),
  Arg("base",Arg::Opt,base,"Island algorithm logarithm base"),
  Arg("lst",Arg::Opt,lst,"Island algorithm linear segment threshold"),
  Arg("ceSepDriven",Arg::Opt,MaxClique::ceSeparatorDrivenInference,"Do separator driven inference (=true) or clique driven (=false)"),
  Arg("componentCache",Arg::Opt,MixtureCommon::cacheMixtureProbabilities,"Cache mixture and component probabilities, faster but uses more memory."),
  Arg("vpap",Arg::Opt,varPartitionAssignmentPrior,"Variable partition assignment priority. Sequence of chars in set [C,D,O,B,I]"),  
  Arg("vcap",Arg::Opt,varCliqueAssignmentPrior,"Variable clique sorting priority. Sequence of chars in set [C,D,O,B,I]"),


  /////////////////////////////////////////////////////////////
  // EM Training Options
  Arg("maxEmIters",Arg::Opt,maxEMIterations,"Max number of EM iterations to do"),
  Arg("random",Arg::Opt,randomizeParams,"Randomize the parameters"),
  // support for vanishing
  Arg("mcvr",Arg::Opt,MixtureCommon::mixCoeffVanishRatio,"Mixture Coefficient Vanishing Ratio"),
  Arg("botForceVanish",Arg::Opt,MixtureCommon::numBottomToForceVanish,"Number of bottom mixture components to force vanish"),
  // support for splitting
  Arg("mcsr",Arg::Opt,MixtureCommon::mixCoeffSplitRatio,"Mixture Coefficient Splitting Ratio"),
  Arg("topForceSplit",Arg::Opt,MixtureCommon::numTopToForceSplit,"Number of top mixture components to force split"),
  Arg("meanCloneSTDfrac",Arg::Opt,MeanVector::cloneSTDfrac,"Fraction of mean to use for STD in mean clone"),
  Arg("covarCloneSTDfrac",Arg::Opt,DiagCovarVector::cloneSTDfrac,"Fraction of var to use for STD in covar clone"),
  Arg("dlinkCloneSTDfrac",Arg::Opt,DlinkMatrix::cloneSTDfrac,"Fraction of var to use for STD in covar clone"),
  Arg("cloneShareMeans",Arg::Opt,GaussianComponent::cloneShareMeans,"Gaussian component clone shares parent mean"),
  Arg("cloneShareCovars",Arg::Opt,GaussianComponent::cloneShareCovars,"Gaussian component clone shares parent covars"),
  Arg("cloneShareDlinks",Arg::Opt,GaussianComponent::cloneShareDlinks,"Gaussian component clone shares parent dlinks"),
  // likelihood difference thresholds
  Arg("lldp",Arg::Opt,lldp,"Log Likelihood difference percentage for termination"),
  Arg("mnlldp",Arg::Opt,mnlldp,"Absolute value of max negative Log Likelihood difference percentage for termination"),
  // EM accumulator file support
  Arg("storeAccFile",Arg::Opt,storeAccFile,"Store accumulators file"),
  Arg("loadAccFile",Arg::Opt,loadAccFile,"Load accumulators file"), 
  Arg("loadAccRange",Arg::Opt,loadAccRange,"Load accumulators file range"), 
  Arg("accFileIsBinary",Arg::Opt,accFileIsBinary,"Binary accumulator files"), 
  // log likelihood store file
  Arg("llStoreFile",Arg::Opt,llStoreFile,"File to store previous sum LL's"), 
  Arg("objsNotToTrain",Arg::Opt,objsToNotTrainFile,"File listing trainable parameter objects to not train."),
  Arg("localCliqueNorm",Arg::Opt,localCliqueNormalization,"Use local clique sum for posterior normalization."),

  // Observation Matrix options
  Arg("pr",  Arg::Opt, pr_str,"per stream per-sentence range",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("fdiffact", Arg::Opt, Action_If_Diff_Num_Frames_Str ,"Action if different number of frames in streams: error (er), repeat last frame (rl), first frame (rf), segmentally expand (se), truncate from start (ts), truncate from end (te)",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("sdiffact", Arg::Opt, Action_If_Diff_Num_Sents_Str ,"Action if different number of sentences in streams: error (er), truncate from end (te), repeat last sent (rl), and wrap around (wa).",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("cppifascii",        Arg::Tog, Cpp_If_Ascii,"Pre-process ASCII files using CPP"),
  Arg("trans",  Arg::Opt,Per_Stream_Transforms ,"per stream transformations string",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("posttrans",  Arg::Opt,Post_Transforms ,"Final global transformations string"),
  Arg("comb",      Arg::Opt, Ftr_Combo_Str,"Combine float features (none: no combination, add, sub, mul,div"),
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
  // Figure out the Endian of the machine this is running on and set the swap defaults accordingly
  bool doWeSwap;

  ByteEndian byteEndian = getWordOrganization();
  switch(byteEndian) {
  case BYTE_BIG_ENDIAN:
    doWeSwap=false;
    break;
  case BYTE_LITTLE_ENDIAN:
    doWeSwap=true;
    break;
  default:
    // We weren't able to figure the Endian out.  Leave the swap defaults as they are.
#ifdef INTV_WORDS_BIGENDIAN
    doWeSwap=true;
#else
    doWeSwap=false;
#endif
  }
  
  for(int i=0; i<MAX_NUM_OBS_FILES; ++i) {
    iswp[i]=doWeSwap;
  }

  ////////////////////////////////////////////
  // parse arguments
  bool parse_was_ok = Arg::parse(argc,(char**)argv);
  if(help) {
    Arg::usage();
    exit(0);
  }
  
  if (print_version_and_exit) {
    printf("%s\n",gmtk_version_id);
    exit(0);
  }

  if(!parse_was_ok) {
    Arg::usage(); exit(-1);
  }

  (void) IM::setGlbMsgLevel(verbosity);


  if (MixtureCommon::cacheMixtureProbabilities)
    MixtureCommon::cacheComponentsInEmTraining = true;
  else 
    MixtureCommon::cacheComponentsInEmTraining = false;

  ////////////////////////////////////////////
  // check for valid argument values.
  int nfiles = 0;
  unsigned ifmts[MAX_NUM_OBS_FILES];
  for (int i=0;i<MAX_NUM_OBS_FILES;i++) {

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

    if (ofs[i] != NULL && ifmts[i]!=PFILE && nfs[i] == 0 && nis[i] == 0)
      error("ERROR: command line parameters must specify one of nf%d and ni%d as not zero",
	    i+1,i+1);
    
    if(ofs[i] != NULL && ifmts[i]==PFILE) {
      FILE *in_fp = fopen(ofs[i], "r");
      if (in_fp==NULL) error("Couldn't open input pfile for reading.");
      bool debug_level=0;
      InFtrLabStream_PFile* in_streamp = new InFtrLabStream_PFile(debug_level,"",in_fp,1,iswp[i]);
      unsigned num_labs=in_streamp->num_labs();
      unsigned num_ftrs=in_streamp->num_ftrs();
      if(nis[i] != 0 && nis[i] != num_labs) error("ERROR: command line parameter ni%d (%d) is different from the one found in the pfile (%d)",i+1,nis[i],num_labs); 
      if(nfs[i] != 0 && nfs[i] != num_ftrs) error("ERROR: command line parameter nf%d (%d) is different from the one found in the pfile (%d)",i+1,nfs[i],num_ftrs); 
      nis[i]=num_labs;
      nfs[i]=num_ftrs;

      if (fclose(in_fp)) error("Couldn't close input pfile.");
      delete in_streamp;
    }
    
    nfiles += (ofs[i] != NULL);
  }

  if (strcmp(Ftr_Combo_Str,"none") == 0)     Ftr_Combo = FTROP_NONE;
  else if (strcmp(Ftr_Combo_Str,"add") == 0) Ftr_Combo = FTROP_ADD;
  else if (strcmp(Ftr_Combo_Str,"sub") == 0) Ftr_Combo = FTROP_SUB;
  else if (strcmp(Ftr_Combo_Str,"mul") == 0) Ftr_Combo = FTROP_MUL;
  else if (strcmp(Ftr_Combo_Str,"div") == 0) Ftr_Combo = FTROP_DIV;
  else error("ERROR: Unknown feature combination type: '%s'\n",Ftr_Combo_Str);
  
  for(int i=0; i < MAX_NUM_OBS_FILES; ++i) {
    if(ofs[i]!=NULL) {
      if (strcmp(Action_If_Diff_Num_Frames_Str[i],"er") == 0)      Action_If_Diff_Num_Frames[i] = ERROR;
      else if (strcmp(Action_If_Diff_Num_Frames_Str[i],"rl") == 0) Action_If_Diff_Num_Frames[i] = REPEAT_LAST;
      else if (strcmp(Action_If_Diff_Num_Frames_Str[i],"rf") == 0) Action_If_Diff_Num_Frames[i] = REPEAT_FIRST;
      else if (strcmp(Action_If_Diff_Num_Frames_Str[i],"se") == 0) Action_If_Diff_Num_Frames[i] = EXPAND_SEGMENTALLY;
      else if (strcmp(Action_If_Diff_Num_Frames_Str[i],"ts") == 0) Action_If_Diff_Num_Frames[i] = TRUNCATE_FROM_START;
      else if (strcmp(Action_If_Diff_Num_Frames_Str[i],"te") == 0) Action_If_Diff_Num_Frames[i] = TRUNCATE_FROM_END;
      else error("ERROR: Unknown action when diff num of frames: '%s'\n",Action_If_Diff_Num_Frames_Str[i]);
    }
  }
  
  for(int i=0; i < MAX_NUM_OBS_FILES; ++i) {
    if(ofs[i]!=NULL) {
      if (strcmp(Action_If_Diff_Num_Sents_Str[i],"er") == 0)      Action_If_Diff_Num_Sents[i] = ERROR;
      else if (strcmp(Action_If_Diff_Num_Sents_Str[i],"rl") == 0) Action_If_Diff_Num_Sents[i] = REPEAT_LAST;
      else if (strcmp(Action_If_Diff_Num_Sents_Str[i],"wa") == 0) Action_If_Diff_Num_Sents[i] = WRAP_AROUND;
      else if (strcmp(Action_If_Diff_Num_Sents_Str[i],"te") == 0) Action_If_Diff_Num_Sents[i] = TRUNCATE_FROM_END;
      else error("ERROR: Unknown action when diff num of sentences: '%s'\n",Action_If_Diff_Num_Sents_Str[i]);
    }
  }

  if (startSkip < 0 || endSkip < 0)
    error("startSkip/endSkip must be >= 0");

  globalObservationMatrix.openFiles(nfiles,
				    (const char**)&ofs,
				    (const char**)&frs,
				    (const char**)&irs,
				    (unsigned*)&nfs,
				    (unsigned*)&nis,
				    (unsigned*)&ifmts,
				    (bool*)&iswps,
				    startSkip,
				    endSkip,
				    Cpp_If_Ascii,
				    cppCommandOptions,
				    (const char**)&pr_str,  //Frame_Range_Str,
				    Action_If_Diff_Num_Frames,
				    Action_If_Diff_Num_Sents,
				    Per_Stream_Transforms,
				    Post_Transforms,
				    Ftr_Combo);


  MixtureCommon::checkForValidRatioValues();
  MeanVector::checkForValidValues();
  DiagCovarVector::checkForValidValues();
  DlinkMatrix::checkForValidValues();
  if (MaxClique::cliqueBeam < 0.0)
    error("cliqueBeam argument must be >= 0");
  if (MaxClique::cliqueBeamRetainFraction <= 0.0 || MaxClique::cliqueBeamRetainFraction > 1.0)
    error("crbeam argument must be: 0.0 < v <= 1.0");
  if (SeparatorClique::separatorBeam < 0.0)
    error("separatorBeam must be >= 0");


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
    if (!fp.readAndVerifyGMId(is))
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
    infoMsg(IM::Default,"WARNING: GMTK is randomizing all trainable parameters and writing them to file 'random.gmp'");
    GM_Parms.makeRandom();
    oDataStreamFile of("random.gmp");
    GM_Parms.writeTrainable(of);
  }

  if (globalObservationMatrix.numSegments()==0) {
    infoMsg(IM::Default,"ERROR: no segments are available in observation file. Exiting...");
    exit_program_with_status(0);
  }

  BP_Range* trrng = new BP_Range(trrng_str,0,globalObservationMatrix.numSegments());
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
      BP_Range lfrng(loadAccRange,0,1000);
      for (BP_Range::iterator lfit=lfrng.begin();
	   lfit<=lfrng.max();
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

  for (unsigned i=0; i<maxEMIterations; i++)  {

    unsigned total_num_frames = 0;

    if (trrng->length() > 0) {
      total_data_prob = 1.0;
      BP_Range::iterator* trrng_it = new BP_Range::iterator(trrng->begin());
      while ((*trrng_it) <= trrng->max()) {
	const unsigned segment = (unsigned)(*(*trrng_it));
	if (globalObservationMatrix.numSegments() < (segment+1)) 
	  error("ERROR: only %d segments in file, training range must be in range [%d,%d] inclusive\n",
		globalObservationMatrix.numSegments(),
		0,globalObservationMatrix.numSegments()-1);

	const unsigned numFrames = GM_Parms.setSegment(segment);
#if 0
	if (globalObservationMatrix.active()) {
	  globalObservationMatrix.printSegmentInfo();
	  ::fflush(stdout);
	}
#endif


	if (island) {
	  if (MixtureCommon::cacheMixtureProbabilities == true) {
	    infoMsg(IM::Default,"NOTE: with island algorithm, might want to also try turning off Gaussian component caching with '-componentCache F'\n"); 
	    fflush(stdout);
	  }
	  unsigned numUsableFrames;
	  myjt.collectDistributeIsland(numFrames,
				       numUsableFrames,
				       base,
				       lst,
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
	      myjt.printAllCliquesProbEvidence();
	    }
	    // And actually train with EM.
	    infoMsg(IM::Low,"Incrementing EM Accumulators\n");
	    myjt.emIncrement(probe,localCliqueNormalization,emTrainingBeam);
	  }
	}
	(*trrng_it)++;
      }
      infoMsg(IM::Default,"Total data log prob from %d frames processed is: %1.9e\n",
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
    GM_Parms.write(outputMasterFile,i);  

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
	previous_dp.valref() = std::exp((double)-20.0);
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
  GM_Parms.write(outputMasterFile);  

  getrusage(RUSAGE_SELF,&rue);
  if (IM::messageGlb(IM::Default)) { 
    infoMsg(IM::Default,"### Final time (seconds) just for EM stage: ");
    double userTime,sysTime;
    reportTiming(rus,rue,userTime,sysTime,stdout);
  }

  exit_program_with_status(0);
}
