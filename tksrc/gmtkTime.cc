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

/*
 * command line arguments
 */

/////////////////////////////////////////////////////////////
// observation input file handling
#define MAX_NUM_OBS_FILES (5)
char    *ofs[MAX_NUM_OBS_FILES] = { NULL, NULL, NULL, NULL,NULL }; 
unsigned nfs[MAX_NUM_OBS_FILES] = { 0, 0, 0,0,0 };
unsigned nis[MAX_NUM_OBS_FILES] = { 0, 0, 0,0,0 };
char   *fmts[MAX_NUM_OBS_FILES] = { "pfile", "pfile", "pfile","pfile","pfile" };
char    *frs[MAX_NUM_OBS_FILES] = { "all", "all", "all","all","all" };
char    *irs[MAX_NUM_OBS_FILES] = { "all", "all", "all","all","all" };
char     *sr[MAX_NUM_OBS_FILES] = { "all", "all", "all","all","all" };
// per stream frame range string before any tranformations are applied
char  *prepr[MAX_NUM_OBS_FILES] = {NULL,NULL,NULL,NULL,NULL};   
// per stream frame range string after per-stream transformations are applied
char *postpr[MAX_NUM_OBS_FILES] = {NULL,NULL,NULL,NULL,NULL};   
char *gpr_str                   = NULL;   // global final frame range string

/////////////////////////////////////////////////////////////
// input parameter/structure file handling
static char *cppCommandOptions = NULL;
static char *inputMasterFile=NULL;
static char *outputMasterFile=NULL;
static char *inputTrainableParameters=NULL;
static bool binInputTrainableParameters=false;
// static char *outputTrainableParameters="outParms%d.gmp";
// static char *outputTrainableParameters=NULL;
// static bool binOutputTrainableParameters=false;
// static bool writeParametersAfterEachEMIteration=true;
// static char *objsToNotTrainFile=NULL;
static int allocateDenseCpts=0;


/////////////////////////////////////////////////////////////
// Structure file, Triangulation File, and Junction Tree Options.
static char *strFileName=NULL;
static char *triFileName=NULL;
static char *jtFileName=NULL;

/////////////////////////////////////////////////////////////
// Continuous RV Options
static double varFloor = GMTK_DEFAULT_VARIANCE_FLOOR;

/////////////////////////////////////////////////////////////
// Beam Options
// static double cliqueBeam=-LZERO;
// static double separatorBeam=-LZERO;

/////////////////////////////////////////////////////////////
// File Range Options
static char *dcdrng_str="all";
static int startSkip = 0;
static int endSkip = 0;

/////////////////////////////////////////////////////////////
// General Options
#define DEF_VERBOSITY (IM::Default-1)
static bool seedme = false;
static unsigned verbosity = DEF_VERBOSITY;
static bool print_version_and_exit = false;
static unsigned seconds = 10;
static bool noEPartition = false;
static unsigned numTimes = 1;
static bool multiTest = false;
static int rlimitSlop = 2;
// static bool limitBest = true;

static unsigned help = 0;  // 0: no help; HIGHEST_PRIORITY (1) ... LOWEST_PRIORITY (5) : increasing levels of help.  The priority levels are defined in arguments.h 

/////////////////////////////////////////////////////////////
// Inference Options
static bool doDistributeEvidence=false;
static bool probE=true;
static bool island=false;
static unsigned base=3;
static unsigned lst=100;
static char* varPartitionAssignmentPrior = "COI";
static char* varCliqueAssignmentPrior = "COT";

////////////////////////////////////
// Observation matrix options
bool     Cpp_If_Ascii        = false;
char*    Action_If_Diff_Num_Frames_Str[MAX_NUM_OBS_FILES]={"er","er","er","er","er"};
unsigned Action_If_Diff_Num_Frames[MAX_NUM_OBS_FILES]={FRAMEMATCH_ERROR,FRAMEMATCH_ERROR,FRAMEMATCH_ERROR,FRAMEMATCH_ERROR,FRAMEMATCH_ERROR};
char*    Action_If_Diff_Num_Sents_Str[MAX_NUM_OBS_FILES]={"te","te","te","te","te"};
unsigned Action_If_Diff_Num_Sents[MAX_NUM_OBS_FILES]={SEGMATCH_TRUNCATE_FROM_END,SEGMATCH_TRUNCATE_FROM_END,SEGMATCH_TRUNCATE_FROM_END,SEGMATCH_TRUNCATE_FROM_END,SEGMATCH_TRUNCATE_FROM_END};

char    *Per_Stream_Transforms[MAX_NUM_OBS_FILES]={NULL,NULL,NULL,NULL,NULL};
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
  Arg("of",  Arg::Req,ofs,"Observation File.  Replace X with the file number",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("nf",  Arg::Opt,nfs,"Number of floats in observation file X",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("ni",  Arg::Opt,nis,"Number of ints in observation file X",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("fmt", Arg::Opt,fmts,"Format (htk,binary,ascii,pfile) for observation file X",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("iswp",Arg::Opt,iswp,"Endian swap condition for observation file X",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("fr",  Arg::Opt,frs,"Float range for observation file X",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("ir",  Arg::Opt,irs,"Int range for observation file X",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("sr",  Arg::Opt,sr,"Sentence range for observation file X",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("prepr", Arg::Opt, prepr,"Frame range for obs file X before any transforms are applied",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("postpr",Arg::Opt, postpr,"Frame range for obs file X after per-stream transforms are applied",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("gpr",    Arg::Opt, gpr_str,"Global final frame range"),

  /////////////////////////////////////////////////////////////
  // input parameter/structure file handling
  Arg("cppCommandOptions",Arg::Opt,cppCommandOptions,"Additional CPP command line"),
  Arg("inputMasterFile",Arg::Req,inputMasterFile,"Input file of multi-level master CPP processed GM input parameters"),
  Arg("outputMasterFile",Arg::Opt,outputMasterFile,"Output file to place master CPP processed GM output parameters"),
  Arg("inputTrainableParameters",Arg::Opt,inputTrainableParameters,"File of only and all trainable parameters"),
  Arg("binInputTrainableParameters",Arg::Opt,binInputTrainableParameters,"Binary condition of trainable parameters file"),
  // Arg("outputTrainableParameters",Arg::Opt,outputTrainableParameters,"File to place only and all trainable output parametes"),
  // Arg("binOutputTrainableParameters",Arg::Opt,binOutputTrainableParameters,"Binary condition of output trainable parameters?"),
  Arg("allocateDenseCpts",Arg::Opt,allocateDenseCpts,"Automatically allocate any undefined CPTs. arg = -1, no read params, arg = 0 noallocate, arg = 1 means use random initial CPT values. arg = 2, use uniform values"),
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
  Arg("cbeam",Arg::Opt,MaxClique::cliqueBeam,"Clique beam width to prune clique (log value)"),
  Arg("cpbeam",Arg::Opt,MaxClique::cliqueBeamBuildBeam,"Clique beam width while building cliques (log value)"),
  Arg("ckbeam",Arg::Opt,MaxClique::cliqueBeamMaxNumStates,"Prune to this clique max state space (0 = no pruning)"),
  Arg("crbeam",Arg::Opt,MaxClique::cliqueBeamRetainFraction,"Fraction of clique state space to retain. Range: 0 < v <= 1. v = 1 means no pruning"),
  Arg("sbeam",Arg::Opt,SeparatorClique::separatorBeam,"Separator beam width pruning log value"),

  /////////////////////////////////////////////////////////////
  // Memory management options
  Arg("clearCliqueValMem",Arg::Opt,MaxClique::perSegmentClearCliqueValueCache,"Free clique/separator value cache for each segment"),


  /////////////////////////////////////////////////////////////
  // File Range Options
  Arg("dcdrng",Arg::Opt,dcdrng_str,"Range to decode over segment file"),
  Arg("startSkip",Arg::Opt,startSkip,"Frames to skip at beginning (i.e., first frame is buff[startSkip])"),
  Arg("endSkip",Arg::Opt,endSkip,"Frames to skip at end (i.e., last frame is buff[len-1-endSkip])"),

  /////////////////////////////////////////////////////////////
  // General Options
  Arg("seed",Arg::Opt,seedme,"Seed the random number generator"),
  Arg("verbosity",Arg::Opt,verbosity,"Verbosity (0 <= v <= 100) of informational/debugging msgs"),
  Arg("version",Arg::Opt,print_version_and_exit,"Print GMTK version number and exit."),
  Arg("seconds",Arg::Opt,seconds,"Number of seconds to run and then exit."),
  Arg("times",Arg::Opt,numTimes,"Number of times to run program seconds seconds long (not multitest mode)."),
  Arg("multiTest",Arg::Opt,multiTest,"Run gmtkTime in multi-test mode, taking triangulation file names from command line."),
  Arg("slop",Arg::Opt,rlimitSlop,"In multiTest mode, number of additional seconds before fail-terminate is forced."),
  Arg("noEPartition",Arg::Opt,noEPartition,"If true, do not run E partition (only [P C C ... C] skipping E)"),


  // Arg("limitBest",Arg::Opt,limitBest,"Limit running time to be approximately best seen so far.."),

  Arg("hashLoadFactor",Arg::Opt,hash_abstract::loadFactor,"Hash table load factor, in range 0.5 <= lf <= 0.95"),
  Arg("help",  Arg::Help, help,  "Print this message. Add an argument from 1 to 5 for increasing help info."),

  /////////////////////////////////////////////////////////////
  // Inference Options
  Arg("doDistributeEvidence",Arg::Opt,doDistributeEvidence,"Do distribute evidence also"),
  Arg("probE",Arg::Opt,probE,"Run the const memory probE function"),
  Arg("island",Arg::Opt,island,"Run island algorithm"),
  Arg("base",Arg::Opt,base,"Island algorithm logarithm base"),
  Arg("lst",Arg::Opt,lst,"Island algorithm linear segment threshold"),
  Arg("ceSepDriven",Arg::Opt,MaxClique::ceSeparatorDrivenInference,"Do separator driven inference (=true) or clique driven (=false)"),
  Arg("componentCache",Arg::Opt,MixtureCommon::cacheMixtureProbabilities,"Cache mixture probabilities, faster but uses more memory."),
  Arg("vpap",Arg::Opt,varPartitionAssignmentPrior,"Variable partition assignment priority. Sequence of chars in set [C,D,O,B,S,I,A,F,N]"),  
  Arg("vcap",Arg::Opt,varCliqueAssignmentPrior,"Variable clique sorting priority. Sequence of chars in set [C,D,O,B,S,I,A,F,N,T,M,+,.]"),
  Arg("jcap",Arg::Opt,JunctionTree::junctionTreeMSTpriorityStr,"Junction Tree Clique MST Sorting Priority. From Set: [D,E,S,U,V,W,H,O,L,Q]"),
  Arg("icap",Arg::Opt,JunctionTree::interfaceCliquePriorityStr,"Interface Clique Priority Determiner Priority. From Set: [W,D,H,O,I]"),
  Arg("useVESeparators",
      Arg::Opt,JunctionTree::useVESeparators,
      "Use Virtual Evidence (VE) Separators (if any are available) during inference"),
  Arg("veSepWhere",
      Arg::Opt,JunctionTree::veSeparatorWhere,
      "Where to use VE seps. Bitwise or of 0x1 (P), 0x2 (C), 0x4 (E)"),
  Arg("veSepFileName",
      Arg::Opt,SeparatorClique::veSeparatorFileName,
      "Name of VE separators file to store VE sep/read previous VE sep info"),
  Arg("veSepRecompute",
      Arg::Opt,SeparatorClique::recomputeVESeparatorTables,
      "Force a re-compute of VE separator information"),  
  Arg("veSepLogProdCardLimit",
      Arg::Opt,SeparatorClique::veSeparatorLogProdCardLimit,
      "The log (base 10) upper limit on a VE sep variable cardinality product"),

  // Observation Matrix transformation options
  Arg("fdiffact",  Arg::Opt, Action_If_Diff_Num_Frames_Str ,"Action if different number of frames in streams: error (er), repeat last frame (rl), first frame (rf), segmentally expand (se), truncate from start (ts), truncate from end (te)",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("sdiffact",  Arg::Opt, Action_If_Diff_Num_Sents_Str ,"Action if different number of sentences in streams: error (er), truncate from end (te), repeat last sent (rl), and wrap around (wa).",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("cppifascii",Arg::Tog, Cpp_If_Ascii,"Pre-process ASCII files using CPP"),
  Arg("trans",     Arg::Opt,Per_Stream_Transforms ,"per stream transformations string",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("posttrans", Arg::Opt,Post_Transforms ,"Final global transformations string"),
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
  
  if (print_version_and_exit) {
    printf("%s\n",gmtk_version_id);
    exit(EXIT_SUCCESS);
  }
  
  if(!parse_was_ok) {
    Arg::usage(); 
    exit(EXIT_FAILURE);
  }

  (void) IM::setGlbMsgLevel(verbosity);
  GM_Parms.setMsgLevel(verbosity);

  // Make sure not to cache the mixture component probabilities as it
  // is only needed in EM training.
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

    if (ofs[i] != NULL && ifmts[i] != PFILE && nfs[i] == 0 && nis[i] == 0)
      error("ERROR: command line parameters must specify one of nf%d and ni%d as not zero",
	    i+1,i+1);

    if(ofs[i] != NULL && ifmts[i]==PFILE) {
      FILE *in_fp = fopen(ofs[i], "r");
      if (in_fp==NULL) error("Couldn't open input pfile for reading.");
      bool debug_level=0;
      InFtrLabStream_PFile* in_streamp = new InFtrLabStream_PFile(debug_level,"",in_fp,1,iswp[i]);
      unsigned num_labs=in_streamp->num_labs();
      unsigned num_ftrs=in_streamp->num_ftrs();

      ////////////////////////////////////////////////////////////
      // Check consistency between pfile and supplied arguments //
      char search_str[]="nXXXXX";
      sprintf(search_str,"-ni%d",i+1);
      bool found=false;
      for(int j=1; j < argc; ++j) {
	if(strcmp(argv[j],search_str)==0) found=true;
      }
      if(found && nis[i] != num_labs) error("ERROR: command line parameter ni%d (%d) is different from the one found in the pfile (%d)",i+1,nis[i],num_labs); 
      sprintf(search_str,"-nf%d",i+1);
      found=false;
      for(int j=1; j < argc; ++j) {
	if(strcmp(argv[j],search_str)==0) found=true;
      }
      if(found && nfs[i] != num_ftrs) error("ERROR: command line parameter nf%d (%d) is different from the one found in the pfile (%d)",i+1,nfs[i],num_ftrs); 
      ////////////////////////////////////////////////////////////
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
       if (strcmp(Action_If_Diff_Num_Frames_Str[i],"er") == 0)      Action_If_Diff_Num_Frames[i] = FRAMEMATCH_ERROR;
       else if (strcmp(Action_If_Diff_Num_Frames_Str[i],"rl") == 0) Action_If_Diff_Num_Frames[i] = FRAMEMATCH_REPEAT_LAST;
       else if (strcmp(Action_If_Diff_Num_Frames_Str[i],"rf") == 0) Action_If_Diff_Num_Frames[i] = FRAMEMATCH_REPEAT_FIRST;
       else if (strcmp(Action_If_Diff_Num_Frames_Str[i],"se") == 0) Action_If_Diff_Num_Frames[i] = FRAMEMATCH_EXPAND_SEGMENTALLY;
       else if (strcmp(Action_If_Diff_Num_Frames_Str[i],"ts") == 0) Action_If_Diff_Num_Frames[i] = FRAMEMATCH_TRUNCATE_FROM_START;
       else if (strcmp(Action_If_Diff_Num_Frames_Str[i],"te") == 0) Action_If_Diff_Num_Frames[i] = FRAMEMATCH_TRUNCATE_FROM_END;
       else error("ERROR: Unknown action when diff num of frames: '%s'\n",Action_If_Diff_Num_Frames_Str[i]);
     }
   }

   for(int i=0; i < MAX_NUM_OBS_FILES; ++i) {
     if(ofs[i]!=NULL) {
       if (strcmp(Action_If_Diff_Num_Sents_Str[i],"er") == 0)      Action_If_Diff_Num_Sents[i] = SEGMATCH_ERROR;
       else if (strcmp(Action_If_Diff_Num_Sents_Str[i],"rl") == 0) Action_If_Diff_Num_Sents[i] = SEGMATCH_REPEAT_LAST;
       else if (strcmp(Action_If_Diff_Num_Sents_Str[i],"wa") == 0) Action_If_Diff_Num_Sents[i] = SEGMATCH_WRAP_AROUND;
       else if (strcmp(Action_If_Diff_Num_Sents_Str[i],"te") == 0) Action_If_Diff_Num_Sents[i] = SEGMATCH_TRUNCATE_FROM_END;
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


  MixtureCommon::checkForValidRatioValues();
  MeanVector::checkForValidValues();
  DiagCovarVector::checkForValidValues();
  DlinkMatrix::checkForValidValues();
  if (MaxClique::cliqueBeam < 0.0)
    error("cliqueBeam argument must be >= 0");
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
      if (!fp.readAndVerifyGMId(is))
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
	    error("gmtkTime doesn't currently support linear full-mem collect/distribute evidence\n");
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
	if (!fp.readAndVerifyGMId(is))
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
