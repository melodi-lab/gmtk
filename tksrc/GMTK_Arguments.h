/*
 * GMTK_Arguments.h
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *  $Header$
 *
 *
 *   Defines all the arguments for all gmtk programs in one place. Defines all three of:
 *       1) the (static) variable names, if not defined elsewhere in GMTK.
 *       2) the argument names, definition, and documentation .
 *       3) code to check that the values of the arguments as specified on the command line is correct.
 * 
 *   The user of this file includes it three times, each time with one of the below defined.
 *       GMTK_ARGUMENTS_DEFINITION     // defines the argument as a static variable if needed
 *       GMTK_ARGUMENTS_DOCUMENTATION  // specifies documentation to arguments.h
 *       GMTK_ARGUMENTS_CHECK_ARGS     // emits code to check that the argument was set ok according to command line, and dies with error if not.
 * 
 *  Each argument is obtained by defining an appropriate #define with a number (0,1,2) with the arguments priority for arguments.h
 *
 */


/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)


/*==============================================================================================================*/
/****************************************************************************************************************/
/****************************************************************************************************************/
/*****************************                                     **********************************************/
/*****************************   OBSERVATION INPUT FILE HANDLING   **********************************************/
/*****************************                                     **********************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_OBS_FILES)

// This next code is used by a number of routines to compute and set
// the default endian swapping condition associated with the
// arguments. We figure out the Endian of the machine this is running
// on and set the swap defaults accordingly.

#define DEF_CODE_TO_COMPUTE_ENDIAN(DEFAULT_SWAP_VALUE)   \
  bool doWeSwap; \
  ByteEndian byteEndian = getWordOrganization(); \
  switch(byteEndian) { \
  case BYTE_BIG_ENDIAN: \
    doWeSwap=false; \
    break; \
  case BYTE_LITTLE_ENDIAN: \
     doWeSwap=true; \
     break; \
  default: \
    /* We weren't able to figure the Endian out.  Leave the swap defaults as they are. */ \
    doWeSwap=DEFAULT_SWAP_VALUE; \
  } \
  \
   for(int i=0; i<MAX_NUM_OBS_FILES; ++i) { \
     iswp[i]=doWeSwap; \
  }

#ifdef INTV_WORDS_BIGENDIAN
#define CODE_TO_COMPUTE_ENDIAN DEF_CODE_TO_COMPUTE_ENDIAN(true) 
#else
#define CODE_TO_COMPUTE_ENDIAN DEF_CODE_TO_COMPUTE_ENDIAN(false) 
#endif


#if defined(GMTK_ARGUMENTS_DEFINITION)

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

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

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
  Arg("gpr",   Arg::Opt, gpr_str,"Global final frame range"),


#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

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

#else
#endif
#endif // defined(GMTK_ARG_OBS_FILES)


/*==============================================================================================================*/
/****************************************************************************************************************/
/****************************************************************************************************************/
/*************************                                            *******************************************/
/*************************   INPUT TRAINABLE PARAMETER FILE HANDLING  *******************************************/
/*************************                                            *******************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/

#if defined(GMTK_ARG_CPP_CMD_OPTS)
#if defined(GMTK_ARGUMENTS_DEFINITION)

  static char *cppCommandOptions = NULL;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("cppCommandOptions",Arg::Opt,cppCommandOptions,"Additional CPP command line"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#else
#endif
#endif // defined(GMTK_ARG_CPP_CMD_OPTS)


/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_INPUT_MASTER_FILE)
#if defined(GMTK_ARGUMENTS_DEFINITION)

   static char *inputMasterFile=NULL;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

   Arg("inputMasterFile",Arg::Req,inputMasterFile,"Input file of multi-level master CPP processed GM input parameters"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGSo)

#else
#endif
#endif // defined(GMTK_ARG_INPUT_MASTER_FILE)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_OUTPUT_MASTER_FILE)
#if defined(GMTK_ARGUMENTS_DEFINITION)

  static char *outputMasterFile=NULL;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("outputMasterFile",Arg::Opt,outputMasterFile,"Output file to place master CPP processed GM output parameters"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#else
#endif
#endif // defined(GMTK_ARG_OUTPUT_MASTER_FILE)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_INPUT_TRAINABLE_PARAMS)
#if defined(GMTK_ARGUMENTS_DEFINITION)

   static char *inputTrainableParameters=NULL;
   static bool binInputTrainableParameters=false;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("inputTrainableParameters",Arg::Opt,inputTrainableParameters,"File of only and all trainable parameters"),
  Arg("binInputTrainableParameters",Arg::Opt,binInputTrainableParameters,"Binary condition of trainable parameters file"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
 
#else
#endif
#endif // defined(GMTK_ARG_INPUT_TRAINABLE_PARAMS

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/

#if defined(GMTK_ARG_OUTPUT_TRAINABLE_PARAMS)
#if defined(GMTK_ARGUMENTS_DEFINITION)

  static char *outputTrainableParameters=NULL;
  static bool binOutputTrainableParameters=false;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("outputTrainableParameters",Arg::Opt,outputTrainableParameters,"File to place only and all trainable output parametes"),
  Arg("binOutputTrainableParameters",Arg::Opt,binOutputTrainableParameters,"Binary condition of output trainable parameters?"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
 
#else
#endif
#endif // defined(GMTK_ARG_OUTPUT_TRAINABLE_PARAMS)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/

#if defined(GMTK_ARG_WPAEEI)
#if defined(GMTK_ARGUMENTS_DEFINITION)

   static bool writeParametersAfterEachEMIteration=true;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

   Arg("wpaeei",Arg::Opt,writeParametersAfterEachEMIteration,"Write Parameters After Each EM Iteration Completes"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
 
#else
#endif
#endif // defined(GMTK_ARG_WPAEEI)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/

#if defined(GMTK_ARG_ALLOC_DENSE_CPTS)
#if defined(GMTK_ARGUMENTS_DEFINITION)

   static int allocateDenseCpts=0;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("allocateDenseCpts",Arg::Opt,allocateDenseCpts,"Automatically allocate undefined CPTs. (-1) = don't read params, (0) = noallocate, (1) = use random initial CPT values, (2) = use uniform values"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
     
  if (allocateDenseCpts != -1 && allocateDenseCpts != 0 && allocateDenseCpts != 1 && allocateDenseCpts != 2)
    error("-allocateDenseCpts argument must be in {-1,0,1,2}\n") ;

#else
#endif
#endif // defined(GMTK_ARG_ALLOC_DENSE_CPTS)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/

#if defined(GMTK_ARG_CPT_NORM_THRES)
#if defined(GMTK_ARGUMENTS_DEFINITION)

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("cptNormThreshold",Arg::Opt,CPT::normalizationThreshold,"Read error if |Sum-1.0|/card > norm_threshold"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#else
#endif
#endif // defined(GMTK_ARG_CPT_NORM_THRES)

/*==============================================================================================================*/
/****************************************************************************************************************/
/****************************************************************************************************************/
/*************************                                            *******************************************/
/*************************   INPUT STRUCTURE PARAMETER FILE HANDLING  *******************************************/
/*************************                                            *******************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/


/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_STR_FILE)
#if defined(GMTK_ARGUMENTS_DEFINITION)

  static char *strFileName=NULL;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("strFile",Arg::Req,strFileName,"Graphical Model Structure File"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#else
#endif
#endif // defined(GMTK_ARG_STR_FILE)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_TRI_FILE)
#if defined(GMTK_ARGUMENTS_DEFINITION)

  static char *triFileName=NULL;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("triFile",Arg::Opt,triFileName,"Triangulation file for strFile"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#else
#endif
#endif // defined(GMTK_ARG_TRI_FILE)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_JT_INFO_FILE)
#if defined(GMTK_ARGUMENTS_DEFINITION)

   static char *jtFileName="jt_info.txt";

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("jtFile",Arg::Opt,jtFileName,"Name of file to write junction tree information"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#else
#endif
#endif // defined(GMTK_ARG_JT_INFO_FILE)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_JTW_UB)
#if defined(GMTK_ARGUMENTS_DEFINITION)


#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("jtwUB",
      Arg::Opt,JunctionTree::jtWeightUpperBound,
      "True means jtWeight is allways an upper bound on true JT weight, false means jtWeight is estimate"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#else
#endif
#endif // defined(GMTK_ARG_JTW_UB)


/*==============================================================================================================*/
/****************************************************************************************************************/
/****************************************************************************************************************/
/*************************                                            *******************************************/
/*************************   CONTINUOUS RANDOM VARIABLE OPTIONS       *******************************************/
/*************************                                            *******************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/


/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_VAR_FLOOR)
#if defined(GMTK_ARGUMENTS_DEFINITION)

  static double varFloor = GMTK_DEFAULT_VARIANCE_FLOOR;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("varFloor",Arg::Opt,varFloor,"Variance Floor (variances can't fall below this value)"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

  // set global variables/change global state from args
  GaussianComponent::setVarianceFloor(varFloor);

#else
#endif
#endif // defined(GMTK_ARG_VAR_FLOOR)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_VAR_FLOOR_ON_READ)
#if defined(GMTK_ARGUMENTS_DEFINITION)

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("floorVarOnRead",Arg::Opt,DiagCovarVector::floorVariancesWhenReadIn,
       "Floor the variances to varFloor when they are read in"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#else
#endif
#endif // defined(GMTK_ARG_VAR_FLOOR_ON_READ)


/*==============================================================================================================*/
/****************************************************************************************************************/
/****************************************************************************************************************/
/*************************                                            *******************************************/
/*************************          BEAM PRUNING OPTIONS              *******************************************/
/*************************                                            *******************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_CBEAM)
#if defined(GMTK_ARGUMENTS_DEFINITION)

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("cbeam",Arg::Opt,MaxClique::cliqueBeam,"Clique beam width to prune clique (log value)"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

  if (MaxClique::cliqueBeam < 0.0)
    error("ERROR: argument cliqueBeam=%f argument must be >= 0",MaxClique::cliqueBeam);

#else
#endif
#endif // defined(GMTK_ARG_CBEAM)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_CPBEAM)
#if defined(GMTK_ARGUMENTS_DEFINITION)

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("cpbeam",Arg::Opt,MaxClique::cliqueBeamBuildBeam,"Clique beam width while building cliques (log value)"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#else
#endif
#endif // defined(GMTK_ARG_CPBEAM)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_CKBEAM)
#if defined(GMTK_ARGUMENTS_DEFINITION)

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("ckbeam",Arg::Opt,MaxClique::cliqueBeamMaxNumStates,"Prune to this clique max state space (0 = no pruning)"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

  if (MaxClique::cliqueBeamMaxNumStates <= 0)
    error("ERROR: ckbeam option must be a positive integer\n");

#else
#endif
#endif // defined(GMTK_ARG_CKBEAM)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_CRBEAM)
#if defined(GMTK_ARGUMENTS_DEFINITION)

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("crbeam",Arg::Opt,MaxClique::cliqueBeamRetainFraction,"Fraction of clique state space to retain. Range: 0 < v <= 1. v = 1 means no pruning"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

  if (MaxClique::cliqueBeamRetainFraction <= 0.0 || MaxClique::cliqueBeamRetainFraction > 1.0)
    error("crbeam argument must be: 0.0 < v <= 1.0");

#else
#endif
#endif // defined(GMTK_ARG_CRBEAM)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/

#if defined(GMTK_ARG_SBEAM)
#if defined(GMTK_ARGUMENTS_DEFINITION)

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("sbeam",Arg::Opt,SeparatorClique::separatorBeam,"Separator beam width pruning log value"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

  if (SeparatorClique::separatorBeam < 0.0)
    error("separatorBeam must be >= 0");

#else
#endif
#endif // defined(GMTK_ARG_SBEAM)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_EBEAM)
#if defined(GMTK_ARGUMENTS_DEFINITION)

  static double emTrainingBeam=-LZERO;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("ebeam",Arg::Opt,emTrainingBeam,"EM training beam width"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#else
#endif
#endif // defined(GMTK_ARG_EBEAM)

/*==============================================================================================================*/
/****************************************************************************************************************/
/****************************************************************************************************************/
/*************************                                            *******************************************/
/*************************          MEMORY MANAGEMENT OPTIONS         *******************************************/
/*************************                                            *******************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_HASH_LOAD_FACTOR)
#if defined(GMTK_ARGUMENTS_DEFINITION)

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("hashLoadFactor",Arg::Opt,hash_abstract::loadFactor,"Hash table load factor, in range 0.05 <= lf <= 1.0"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

    if (hash_abstract::loadFactor < 0.05 || hash_abstract::loadFactor >= 1.0) 
      error("ARG ERROR: hashLoadFactor must be between 0.05 and 1.0 non-inclusive");

#else
#endif
#endif // defined(GMTK_ARG_HASH_LOAD_FACTOR)


/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/

#if defined(GMTK_ARG_CLEAR_CLIQUE_VAL_MEM)
#if defined(GMTK_ARGUMENTS_DEFINITION)

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("clearCliqueValMem",Arg::Opt,MaxClique::perSegmentClearCliqueValueCache,"Free clique/separator value cache for each segment"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#else
#endif
#endif // defined(GMTK_ARG_CLEAR_CLIQUE_VAL_MEM)

/*==============================================================================================================*/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************                                     ***********************************************/
/****************************      FILE RANGE OPTIONS             ***********************************************/
/****************************                                     ***********************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/


/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_DCDRNG)
#if defined(GMTK_ARGUMENTS_DEFINITION)

  static char *dcdrng_str="all";

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("dcdrng",Arg::Opt,dcdrng_str,"Range to decode over segment file"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#else
#endif
#endif // defined(GMTK_ARG_DCDRNG)


/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_TRRNG)
#if defined(GMTK_ARGUMENTS_DEFINITION)

  static char *trrng_str="all";

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("trrng",Arg::Opt,trrng_str,"Range to decode over segment file"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#else
#endif
#endif // defined(GMTK_ARG_TRRNG)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_START_END_SKIP)
#if defined(GMTK_ARGUMENTS_DEFINITION)

  static int startSkip = 0;
  static int endSkip = 0;


#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("startSkip",Arg::Opt,startSkip,"Frames to skip at beginning (i.e., first frame is buff[startSkip])"),
  Arg("endSkip",Arg::Opt,endSkip,"Frames to skip at end (i.e., last frame is buff[len-1-endSkip])"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

  if (startSkip < 0 || endSkip < 0)
    error("ERROR: arguments startSkip=%d/endSkip=%d must both be >= 0",startSkip,endSkip);

#else
#endif
#endif // defined(GMTK_ARG_START_END_SKIP)

/*==============================================================================================================*/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************                                     ***********************************************/
/****************************         GENERAL OPTIONS             ***********************************************/
/****************************                                     ***********************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/


/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_SEED)
#if defined(GMTK_ARGUMENTS_DEFINITION)

  static bool seedme = false;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("seed",Arg::Opt,seedme,"Seed the random number generator"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

  if (seedme)
    rnd.seed();

#else
#endif
#endif // defined(GMTK_ARG_SEED)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_VERB)
#if defined(GMTK_ARGUMENTS_DEFINITION)

  static unsigned verbosity = IM::Default;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("verbosity",Arg::Opt,verbosity,"Verbosity (0 <= v <= 100) of informational/debugging msgs"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

  (void) IM::setGlbMsgLevel(verbosity);
  GM_Parms.setMsgLevel(verbosity);

#else
#endif
#endif // defined(GMTK_ARG_VERB)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_HELP)
#if defined(GMTK_ARGUMENTS_DEFINITION)

   // 0: no help; HIGHEST_PRIORITY (1) ... LOWEST_PRIORITY (5) : increasing levels of help.  The priority levels are defined in arguments.h 
   static unsigned help = 0;  

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("help",  Arg::Help, help,  "Print this message. Add an argument from 1 to 5 for increasing help info."),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

  if(help) {
    Arg::usage();
    exit(0);
  }

#else
#endif
#endif // defined(GMTK_ARG_HELP)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_VERSION)
#if defined(GMTK_ARGUMENTS_DEFINITION)

  static bool print_version_and_exit = false;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("version",Arg::Opt,print_version_and_exit,"Print GMTK version number and exit."),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

  if (print_version_and_exit) {
    printf("%s\n",gmtk_version_id);
    exit(0);
  }


#else
#endif
#endif // defined(GMTK_ARG_VERSION)



/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_CLIQUE_PRINT)
#if defined(GMTK_ARGUMENTS_DEFINITION)

static char* pPartCliquePrintRange = NULL;
static char* cPartCliquePrintRange = NULL;
static char* ePartCliquePrintRange = NULL;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("pCliquePrintRange",Arg::Opt,pPartCliquePrintRange,"With CE/DE, print range cliques from P partition."),
  Arg("cCliquePrintRange",Arg::Opt,cPartCliquePrintRange,"With CE/DE, print range cliques from C partition."),
  Arg("eCliquePrintRange",Arg::Opt,ePartCliquePrintRange,"With CE/DE, print range cliques from E partition."),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)



#else
#endif
#endif // defined(GMTK_ARG_CLIQUE_PRINT)


/*==============================================================================================================*/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************                                     ***********************************************/
/****************************         INFERENCE OPTIONS           ***********************************************/
/****************************                                     ***********************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_ISLAND)
#if defined(GMTK_ARGUMENTS_DEFINITION)

  static bool island=false;
  static unsigned base=3;
  static unsigned lst=100;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("island",Arg::Opt,island,"Run island algorithm"),
  Arg("base",Arg::Opt,base,"Island algorithm logarithm base"),
  Arg("lst",Arg::Opt,lst,"Island algorithm linear segment threshold"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

  // these options are checked by the island algorithm code.

#else
#endif
#endif // defined(GMTK_ARG_ISLAND)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_CE_SEP_DRIVEN)
#if defined(GMTK_ARGUMENTS_DEFINITION)

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("ceSepDriven",Arg::Opt,MaxClique::ceSeparatorDrivenInference,"Do separator driven inference (=true) or clique driven (=false)"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#else
#endif
#endif // defined(GMTK_ARG_CE_SEP_DRIVEN)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_COMPONENT_CACHE)
#if defined(GMTK_ARGUMENTS_DEFINITION)

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("componentCache",Arg::Opt,MixtureCommon::cacheMixtureProbabilities,"Cache mixture and component probabilities, faster but uses more memory."),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

  if (MixtureCommon::cacheMixtureProbabilities)
    MixtureCommon::cacheComponentsInEmTraining = true;
  else 
    MixtureCommon::cacheComponentsInEmTraining = false;

#else
#endif
#endif // defined(GMTK_ARG_COMPONENT_CACHE)



/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_MIXTURE_CACHE)
#if defined(GMTK_ARGUMENTS_DEFINITION)

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("componentCache",Arg::Opt,MixtureCommon::cacheMixtureProbabilities,"Cache mixture probabilities, faster but uses more memory."),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

  // Make sure not to cache the mixture component probabilities as it
  // is only needed in EM training.
  MixtureCommon::cacheComponentsInEmTraining = false;

#else
#endif
#endif // defined(GMTK_ARG_MIXTURE_CACHE)


/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_CLIQUE_VAR_ITER_ORDERS)
#if defined(GMTK_ARGUMENTS_DEFINITION)

static char* varPartitionAssignmentPrior = "COI";
static char* varCliqueAssignmentPrior = "COT";

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("vpap",Arg::Opt,varPartitionAssignmentPrior,"Variable partition assignment priority. Sequence of chars in set [C,D,O,B,S,I,A,F,N]"),  
  Arg("vcap",Arg::Opt,varCliqueAssignmentPrior,"Variable clique sorting priority. Sequence of chars in set [C,D,O,B,S,I,A,F,N,T,M,+,.]"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#else
#endif
#endif // defined(GMTK_ARG_CLIQUE_VAR_ITER_ORDERS)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_JT_OPTIONS)
#if defined(GMTK_ARGUMENTS_DEFINITION)

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("jcap",Arg::Opt,JunctionTree::junctionTreeMSTpriorityStr,"Junction Tree Clique MST Sorting Priority. From Set: [D,E,S,U,V,W,H,O,L,Q]"),
  Arg("icap",Arg::Opt,JunctionTree::interfaceCliquePriorityStr,"Interface Clique Priority Determiner Priority. From Set: [W,D,H,O,I]"),


#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#else
#endif
#endif // defined(GMTK_ARG_JT_OPTIONS)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/

#if defined(GMTK_ARG_VE_SEPS)
#if defined(GMTK_ARGUMENTS_DEFINITION)

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

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
  
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#else
#endif
#endif // defined(GMTK_ARG_VE_SEPS)



/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_VITERBI_SCORE)
#if defined(GMTK_ARGUMENTS_DEFINITION)

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("viterbiScore",Arg::Opt,JunctionTree::viterbiScore,"Compute p(o,h_max) (rather than sum_h p(o,h))"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#else
#endif
#endif // defined(GMTK_ARG_VITERBI_SCORE)



/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_DO_DIST_EVIDENCE)
#if defined(GMTK_ARGUMENTS_DEFINITION)

  static bool doDistributeEvidence=false;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("doDistributeEvidence",Arg::Opt,doDistributeEvidence,"Also run distribute-evidence"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#else
#endif
#endif // defined(GMTK_ARG_DO_DIST_EVIDENCE)



/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_PROB_EVIDENCE)
#if defined(GMTK_ARGUMENTS_DEFINITION)

  static bool probE=false;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("probE",Arg::Opt,probE,"Run the constant memory prob(evidence) function"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#else
#endif
#endif // defined(GMTK_ARG_PROB_EVIDENCE)


/*==============================================================================================================*/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************                                     ***********************************************/
/****************************         EM TRAINING OPTIONS         ***********************************************/
/****************************                                     ***********************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_EM_TRAINING_PARAMS)
#if defined(GMTK_ARGUMENTS_DEFINITION)

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

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

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


#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

  MixtureCommon::checkForValidRatioValues();
  MeanVector::checkForValidValues();
  DiagCovarVector::checkForValidValues();
  DlinkMatrix::checkForValidValues();



#else
#endif
#endif // defined(GMTK_ARG_EM_TRAINING_PARAMS)

/*==============================================================================================================*/
/****************************************************************************************************************/
/****************************************************************************************************************/
/************************                                              ******************************************/
/************************  OBSERVATION MATRIX TRANSFORMATION OPTIONS   ******************************************/
/************************                                              ******************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_OBS_MATRIX_XFORMATION)
#if defined(GMTK_ARGUMENTS_DEFINITION)

bool     Cpp_If_Ascii        = false;

char*    Action_If_Diff_Num_Frames_Str[MAX_NUM_OBS_FILES]={"er","er","er","er","er"};   // 
unsigned Action_If_Diff_Num_Frames[MAX_NUM_OBS_FILES]={FRAMEMATCH_ERROR,FRAMEMATCH_ERROR,FRAMEMATCH_ERROR,FRAMEMATCH_ERROR,FRAMEMATCH_ERROR};   // 
char*    Action_If_Diff_Num_Sents_Str[MAX_NUM_OBS_FILES]={"te","te","te","te","te"}; 
unsigned Action_If_Diff_Num_Sents[MAX_NUM_OBS_FILES]={SEGMATCH_TRUNCATE_FROM_END,SEGMATCH_TRUNCATE_FROM_END,SEGMATCH_TRUNCATE_FROM_END,SEGMATCH_TRUNCATE_FROM_END,SEGMATCH_TRUNCATE_FROM_END};   // 

char    *Per_Stream_Transforms[MAX_NUM_OBS_FILES]={NULL,NULL,NULL,NULL,NULL};   // 
char    *Post_Transforms=NULL;

char    *Ftr_Combo_Str="none";
unsigned Ftr_Combo=FTROP_NONE;
 
#ifdef INTV_WORDS_BIGENDIAN
bool iswp[MAX_NUM_OBS_FILES] = {true,true,true,true,true};
#else
bool iswp[MAX_NUM_OBS_FILES] = {false,false,false,false,false};
#endif 

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("fdiffact",  Arg::Opt, Action_If_Diff_Num_Frames_Str ,"Action if different number of frames in streams: error (er), repeat last frame (rl), first frame (rf), segmentally expand (se), truncate from start (ts), truncate from end (te)",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("sdiffact",  Arg::Opt, Action_If_Diff_Num_Sents_Str ,"Action if different number of sentences in streams: error (er), truncate from end (te), repeat last sent (rl), and wrap around (wa).",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("cppifascii",Arg::Tog, Cpp_If_Ascii,"Pre-process ASCII files using CPP"),
  Arg("trans",     Arg::Opt,Per_Stream_Transforms ,"per stream transformations string",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("posttrans", Arg::Opt,Post_Transforms ,"Final global transformations string"),
  Arg("comb",      Arg::Opt, Ftr_Combo_Str,"Combine float features (none: no combination, add, sub, mul,div"),


#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

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


#else
#endif
#endif // defined(GMTK_ARG_OBS_MATRIX_XFORMATION)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/




#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_XX_XX)
#if defined(GMTK_ARGUMENTS_DEFINITION)
--
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
--
#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
--
#else
#endif
#endif // defined(GMTK_ARG_XX_XX)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/

