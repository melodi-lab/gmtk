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

// TODO: remove next 2 eventually 
#include <iostream.h>
#include <fstream.h>

#include "general.h"
#include "error.h"
#include "debug.h"
#include "rand.h"
#include "arguments.h"
#include "ieeeFPsetup.h"
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
char    *ofs[MAX_NUM_OBS_FILES] = { NULL, NULL, NULL, NULL,NULL }; 
unsigned nfs[MAX_NUM_OBS_FILES] = { 0, 0, 0,0,0 };
unsigned nis[MAX_NUM_OBS_FILES] = { 0, 0, 0,0,0 };
char   *fmts[MAX_NUM_OBS_FILES] = { "pfile", "pfile", "pfile","pfile","pfile" };
bool   iswps[MAX_NUM_OBS_FILES] = { false, false, false,false,false };
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
static char *inputTrainableParameters=NULL;
static bool binInputTrainableParameters=false;
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
// static double cliqueBeam=-LZERO;
// static double separatorBeam=-LZERO;

/////////////////////////////////////////////////////////////
// File Range Options
static char *dcdrng_str="all";
static int startSkip = 0;
static int endSkip = 0;

/////////////////////////////////////////////////////////////
// General Options
static bool seedme = false;
static unsigned verbosity = IM::Default;
static bool help = false;
static bool print_version_and_exit = false;

/////////////////////////////////////////////////////////////
// Inference Options
static bool island=false;
static unsigned base=3;
static unsigned lst=100;
static char* varPartitionAssignmentPrior = "COI";
static char* varCliqueAssignmentPrior = "COI";


/////////////////////////////////////////////////////////////
// Decoding Options
static char *dumpNames = NULL;
static char *ofilelist = NULL;
static char *wordVar=NULL;
static char *varMapFile=NULL;
static char *transitionLabel=NULL;
static char* showVitVals = NULL;

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
  /////////////////////////////////////////////////////////////
  // observation input file handling
  Arg("of",  Arg::Req,ofs,"Observation File.  Replace X with the file number",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("nf",  Arg::Opt,nfs,"Number of floats in observation file X",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("ni",  Arg::Opt,nis,"Number of ints in observation file X",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("fmt", Arg::Opt,fmts,"Format (htk,binary,ascii,pfile) for observation file X",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("iswp",Arg::Opt,iswps,"Endian swap condition for observation file X",Arg::ARRAY,MAX_NUM_OBS_FILES),
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
  Arg("inputTrainableParameters",Arg::Opt,inputTrainableParameters,"File of only and all trainable parameters"),
  Arg("binInputTrainableParameters",Arg::Opt,binInputTrainableParameters,"Binary condition of trainable parameters file"),
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
  // Decoding Options
  Arg("dumpNames",Arg::Opt,dumpNames,"File containing the names of the variables to save to a file"),
  Arg("ofilelist",Arg::Opt,ofilelist,"List of filenames to dump the hidden variable values to"),
  // These 3 must be used together or not at all
  Arg("printWordVar",Arg::Opt,wordVar,"Print the word var - which has this label"),
  Arg("varMap",Arg::Opt,varMapFile,"Use this file to map from word-index to string"),
  Arg("transitionLabel",Arg::Opt,transitionLabel,"The label of the word transition variable"),
  Arg("showVitVals",Arg::Opt,showVitVals,"File to print viterbi values, '-' for stdout"),
   
  // Observation Matrix transformation options
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

  if (dumpNames)
    if (ofilelist==NULL) error("Must also specify output files for binary writing");

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
				    (bool*)&iswps,
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


  if (globalObservationMatrix.numSegments()==0) {
    infoMsg(IM::Default,"ERROR: no segments are available in observation file. Exiting...");
    exit_program_with_status(0);
  }

  Range* dcdrng = new Range(dcdrng_str,0,globalObservationMatrix.numSegments());
#if 0
  if (dcdrng->length() <= 0) {
    infoMsg(IM::Default,"Training range '%s' specifies empty set. Exiting...\n",
	  dcdrng_str);
    exit_program_with_status(0);
  }
#endif

  if (dcdrng->length() == 0) { 
    error("Decoding range must specify a non-zero length range. Range given is %s\n",
	  dcdrng_str);
  }

  logpr total_data_prob = 1.0;

  FILE* vitValsFile;
  if (showVitVals) {
    if (strcmp(showVitVals,"-") == 0)
      vitValsFile = stdout;
    else {
      if ((vitValsFile = fopen(showVitVals, "w")) == NULL)
	error("Can't open file '%s' for writing\n",vitValsFile);
    }
  } else {
    vitValsFile = NULL;
  }

  map<int, string> word_map;
  if (wordVar != NULL)
  {
    if (varMapFile==NULL)
      error("File to map from word index to string not specified.");
    if(transitionLabel == NULL)
      error("The label of the transition variable was not specified.");
    ifstream in(varMapFile);
    if (!in) { cout << "Unable to open " << varMapFile << endl; exit(1); }
    string name;
    int val;
    while (!in.eof())
    {
      in >> val >> name >> ws;
      // printf("reading %d word = (%s)\n",val,name.c_str());
      word_map[val] = name;
    }
    in.close();
  }

  set<string> dumpVars;
  map<string,int> posFor;
  int ndv=0;
  if (dumpNames)
  {
      ifstream din(dumpNames);
      if (!din) {cout << "Unable to open " << dumpNames << endl; exit(1);}
      string s;
      while (!din.eof())
      {
          din >> s >> ws;
          dumpVars.insert(s);
          posFor[s] = ndv++;  // which position within a frame to put it
      }
      din.close();
  }

  vector<string> ofiles;
  if (ofilelist)
  {
      ifstream oin(ofilelist);
      if (!oin) {cout << "Unable to open " << ofilelist << endl; exit(1);} 
      string s;
      while (!oin.eof())
      {
          oin >> s >> ws;
          ofiles.push_back(s);
      }
      oin.close();
  }
  // We always do viterbi scoring/option in this program.
  JunctionTree::viterbiScore = true;

  struct rusage rus; /* starting time */
  struct rusage rue; /* ending time */
  getrusage(RUSAGE_SELF,&rus);

  total_data_prob = 1.0;
  Range::iterator* dcdrng_it = new Range::iterator(dcdrng->begin());

  if (island) {
    if (MixtureCommon::cacheMixtureProbabilities == true) {
      infoMsg(IM::Default,"NOTE: with island algorithm, might want to also try turning off Gaussian component caching with '-componentCache F'\n"); 
      fflush(stdout);
    }
  }
    
  while (!dcdrng_it->at_end()) {
    const unsigned segment = (unsigned)(*(*dcdrng_it));
    if (globalObservationMatrix.numSegments() < (segment+1)) 
      error("ERROR: only %d segments in file, decode range must be in range [%d,%d] inclusive\n",
	    globalObservationMatrix.numSegments(),
	    0,globalObservationMatrix.numSegments()-1);

    const unsigned numFrames = GM_Parms.setSegment(segment);

    logpr probe;
    if (island) {
      // warning("WARNING:: Island algorithm for decoding not yet debugged (but almost). Use at your own risk!!\n");
      unsigned numUsableFrames;
      myjt.collectDistributeIsland(numFrames,
				   numUsableFrames,
				   base,
				   lst,
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
      infoMsg(IM::Med,"Collecting Evidence\n");
      myjt.collectEvidence();
      infoMsg(IM::Med,"Done Collecting Evidence\n");
      probe = myjt.probEvidence();
      infoMsg(IM::Default,"Segment %d, after CE, viterbi log(prob(evidence)) = %f, per frame =%f, per numUFrams = %f\n",
	     segment,
	     probe.val(),
	     probe.val()/numFrames,
	     probe.val()/numUsableFrames);
      if (probe.essentially_zero()) {
	infoMsg(IM::Default,"Skipping segment since probability is essentially zero\n");
      } else {
	myjt.setRootToMaxCliqueValue();
	total_data_prob *= probe;
	infoMsg(IM::Low,"Distributing Evidence\n");
	myjt.distributeEvidence();
	infoMsg(IM::Low,"Done Distributing Evidence\n");
      }
    }
    if (vitValsFile)
      myjt.printCurrentRVValues(vitValsFile);

    if (dumpNames)
    {
        if (segment >= ofiles.size()) 
	  error("More utterances than output files");
        FILE *fp = fopen(ofiles[segment].c_str(), "wb");
        if (!fp) {cout << "Unable to open " << ofiles[segment] << endl; exit(1);}
        int *vals = new int[myjt.curNodes().size()];
        // just in case a variable isn't present in a slice, write a -1
        for (int i=0; i<int(myjt.curNodes().size()); i++) vals[i] = -1;
        int nv=0;
        for (int i=0; i<int(myjt.curNodes().size()); i++)
            if (dumpVars.count(myjt.curNodes()[i]->name())) 
            {
                if (!myjt.curNodes()[i]->discrete() || !myjt.curNodes()[i]->hidden()) 
                    error("variables to dump must be discrete and hidden");
                int p = myjt.curNodes()[i]->frame()*dumpVars.size() 
                        + posFor[myjt.curNodes()[i]->name()];
                assert(p<int(myjt.curNodes().size()));
                vals[p] = RV2DRV(myjt.curNodes()[i])->val; 
                nv++;
            }
        fwrite((void *) vals, sizeof(int), nv, fp); 
        delete [] vals;
        fclose(fp);
    }
    if (wordVar && probe.not_essentially_zero())
    {
      // print the sequence of values for this variable
      // compress consecutive values into a single instance
      // the times are right if a word transition at time t means there is
      // a new word at t+1
      string pvn = string(wordVar);
      string tl = string(transitionLabel);
      for (int i=0, lv=-1, lf=0; i<int(myjt.curNodes().size()); i++)
      {
        if (myjt.curNodes()[i]->name() == pvn)
        {
          if (!myjt.curNodes()[i]->discrete()) 
            error("Can only print Viterbi values for discrete variables");
          if (RV2DRV(myjt.curNodes()[i])->cardinality != word_map.size())
            error("Word-val to string map file size %d does not match the number of words %d.",int(word_map.size()),RV2DRV(myjt.curNodes()[i])->cardinality);
          lv = RV2DRV(myjt.curNodes()[i])->val;
        }
        else if (myjt.curNodes()[i]->name()==tl)
        {
          if (!myjt.curNodes()[i]->discrete()) 
            error("Word transition variable should be discrete");
          if (RV2DRV(myjt.curNodes()[i])->cardinality != 2) 
            error("Word transition variable should have two values");
          if (RV2DRV(myjt.curNodes()[i])->val==1)  // a word transition
          {
            cout << word_map[lv] << " (" << lf << "-" 
                 << myjt.curNodes()[i]->frame() << ")\n" << flush;  
            lf = myjt.curNodes()[i]->frame()+1;
          }
        }
      }
    }
    (*dcdrng_it)++;
  }

  infoMsg(IM::Default,"Total data log prob is: %1.9e\n",
	  total_data_prob.val());

  if (vitValsFile && vitValsFile != stdout)
    fclose(vitValsFile);

  getrusage(RUSAGE_SELF,&rue);
  if (IM::messageGlb(IM::Default)) { 
    infoMsg(IM::Default,"### Final time (seconds) just for decoding stage: ");
    double userTime,sysTime;
    reportTiming(rus,rue,userTime,sysTime,stdout);
  }

  exit_program_with_status(0);
}
