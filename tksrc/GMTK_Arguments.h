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

#ifdef HAVE_CONFIG_H
#include <config.h>

#ifdef HAVE_HG_H
#include "hgstamp.h"
#endif
#endif

#include "GMTK_ProgramDefaultParms.h"

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/

/* initial definitions commonto all arguments */

#if defined(GMTK_ARGUMENTS_DEFINITION)

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

const char*const argerr = "ARG ERROR";
// include nop statement to avoid warning message.
(void)argerr;


#else
#endif



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
#if defined(GMTK_ARGUMENTS_DEFINITION)


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


   // observation input file handling
#define MAX_NUM_OBS_FILES (5)
   char    *ofs[MAX_NUM_OBS_FILES] = { NULL, NULL, NULL, NULL,NULL }; 
   unsigned nfs[MAX_NUM_OBS_FILES] = { 0, 0, 0,0,0 };
   unsigned nis[MAX_NUM_OBS_FILES] = { 0, 0, 0,0,0 };
   const char   *fmts[MAX_NUM_OBS_FILES] = { "pfile", "pfile", "pfile","pfile","pfile" };
   const char    *frs[MAX_NUM_OBS_FILES] = { "all", "all", "all","all","all" };
   const char    *irs[MAX_NUM_OBS_FILES] = { "all", "all", "all","all","all" };
   const char     *sr[MAX_NUM_OBS_FILES] = { "all", "all", "all","all","all" };
   // per stream frame range string before any tranformations are applied
   char  *prepr[MAX_NUM_OBS_FILES] = {NULL,NULL,NULL,NULL,NULL};   
   // per stream frame range string after per-stream transformations are applied
   char *postpr[MAX_NUM_OBS_FILES] = {NULL,NULL,NULL,NULL,NULL};   
   char *gpr_str                   = NULL;   // global final frame range string

extern bool ObservationsAllowNan;

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
  Arg("prepr", Arg::Opt, prepr,"Pre Per-segment frame Range for obs file X before any transforms are applied",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("postpr",Arg::Opt, postpr,"Post Per-segment frame Range for obs file X after per-stream transforms are applied",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("gpr",   Arg::Opt, gpr_str," Global Per-segment final frame Range"),
  Arg("obsNAN",   Arg::Opt, ObservationsAllowNan," True if observation files allow FP NAN values"),


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
      error("%s: Unknown observation file format type: '%s'\n",argerr,fmts[i]);

    if (ofs[i] != NULL && ifmts[i]!=PFILE && nfs[i] == 0 && nis[i] == 0)
      error("%s: command line parameters must specify one of nf%d and ni%d as not zero",argerr,
	    i+1,i+1);
    
    if(ofs[i] != NULL && ifmts[i]==PFILE) {
      FILE *in_fp = fopen(ofs[i], "r");
      if (in_fp==NULL) 
	error("Couldn't open input pfile for reading.");
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
      if(found && nis[i] != num_labs) 
	error("%s: command line parameter ni%d (%d) is different from the one found in the pfile (%d)",argerr,
	      i+1,nis[i],num_labs); 
      sprintf(search_str,"-nf%d",i+1);
      found=false;
      for(int j=1; j < argc; ++j) {
	if(strcmp(argv[j],search_str)==0) found=true;
      }
      if(found && nfs[i] != num_ftrs) 
	error("%s: command line parameter nf%d (%d) is different from the one found in the pfile (%d)",
	      argerr,i+1,nfs[i],num_ftrs); 
      ////////////////////////////////////////////////////////////
      nis[i]=num_labs;
      nfs[i]=num_ftrs;

      if (fclose(in_fp)) 
	error("Couldn't close input pfile.");
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

  Arg("cppCommandOptions",Arg::Opt,cppCommandOptions,"Command line options to give to 'cpp'"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#else
#endif
#endif // defined(GMTK_ARG_CPP_CMD_OPTS)


/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_INPUT_MASTER_FILE) || defined(GMTK_ARG_INPUT_MASTER_FILE_OPT_ARG)
#if defined(GMTK_ARGUMENTS_DEFINITION)

   static char *inputMasterFile=NULL;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

#ifdef GMTK_ARG_INPUT_MASTER_FILE_OPT_ARG
   Arg("inputMasterFile",Arg::Opt,inputMasterFile,"Input file of multi-level master CPP processed GM input parameters"),
#else
   Arg("inputMasterFile",Arg::Req,inputMasterFile,"Input file of multi-level master CPP processed GM input parameters"),
#endif

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

  Arg("allocateDenseCpts",Arg::Opt,allocateDenseCpts,"Automatically allocate undefined CPTs. (-1) = don't read params, (0) = don't allocate, (1) = use random initial CPT values, (2) = use uniform values"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
     
  if (allocateDenseCpts != -1 && allocateDenseCpts != 0 && allocateDenseCpts != 1 && allocateDenseCpts != 2)
    error("%s: -allocateDenseCpts argument must be in {-1,0,1,2}\n",argerr) ;

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
/*************************   INPUT STRUCTURE/TRI/JT FILE HANDLING     *******************************************/
/*************************                                            *******************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/


/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_STR_FILE) || defined(GMTK_ARG_STR_FILE_OPT_ARG)
#if defined(GMTK_ARGUMENTS_DEFINITION)

  static char *strFileName=NULL;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

#ifdef GMTK_ARG_STR_FILE_OPT_ARG
  Arg("strFile",Arg::Opt,strFileName,"Graphical Model Structure File"),
#else
  Arg("strFile",Arg::Req,strFileName,"Graphical Model Structure File"),
#endif

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#else
#endif
#endif // defined(GMTK_ARG_STR_FILE)



/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/

#if defined(GMTK_ARG_CHECK_TRI_FILE_CARD)
#if defined(GMTK_ARGUMENTS_DEFINITION)

  static bool checkTriFileCards=true;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("checkTriFileCards",Arg::Opt,checkTriFileCards,"Verify rv cardinalities in triangulation file match .str file"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#else
#endif
#endif // defined(GMTK_ARG_CHECK_TRI_FILE_CARD)


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

// 
// This option is meant just for gmtkTriangulate.cc, a program that
// normally outputs a triangulation file, but sometimes (e.g., when
// called by scripts, when wanting to triangulate just one partition
// at a time, or when wanting just to check various stats about a
// trifile such as clique weight), we need to read in a pre-existing
// tri-file. These options are for that purpose.  Note that we
// shouldn't have both *INPUT_TRI_FILE and *TRI_FILE defined in the
// same file.

#if defined(GMTK_ARG_INPUT_TRI_FILE)
#if defined(GMTK_ARGUMENTS_DEFINITION)

  static char *inputTriangulatedFile=NULL;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("inputTriangulatedFile",Arg::Opt,inputTriangulatedFile,"Non-default previous triangulated file to start with"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#else
#endif
#endif // defined(GMTK_ARG_INPUT_TRI_FILE)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_OUTPUT_TRI_FILE)
#if defined(GMTK_ARGUMENTS_DEFINITION)

  static char *outputTriangulatedFile=NULL;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("outputTriangulatedFile",Arg::Opt,outputTriangulatedFile,"File name to write resulting triangulation to"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#else
#endif
#endif // defined(GMTK_ARG_OUTPUT_TRI_FILE)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_CROSSOVER_OPTIONS)
#if defined(GMTK_ARGUMENTS_DEFINITION)

  static char *inputCrossoverTriangulatedFile=NULL;
  static char *outputCrossoverTriangulatedFile=NULL;

  static float crossoverProbability = 0.2;
  static float mutateProbability = 0.7;
                                                                                
#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("inputCrossoverTriangulatedFile",Arg::Opt,inputCrossoverTriangulatedFile,
    "Non-default previous triangulated file to start with"),
  Arg("outputCrossoverTriangulatedFile",Arg::Opt,outputCrossoverTriangulatedFile, 
    "File name to write second resulting triangulation to"),
  Arg("crossoverProbability", Arg::Opt, crossoverProbability, 
    "Probability of an edge swap when using crossover"),
  Arg("mutateProbability", Arg::Opt, mutateProbability, 
    "Probability of an edge mutation when using crossover"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#else
#endif
#endif 

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_JT_INFO_FILE)
#if defined(GMTK_ARGUMENTS_DEFINITION)

#ifdef GMTK_ARG_JT_INFO_FILE_DEF_VAL
    static char *jtFileName = GMTK_ARG_JT_INFO_FILE_DEF_VAL;
#else
    const static char *jtFileName = "jt_info.txt";
#endif

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
    error("%s: argument cliqueBeam=%f argument must be >= 0",argerr,MaxClique::cliqueBeam);

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
  Arg("cpfilter",Arg::Opt,MaxClique::cliqueBeamBuildFilter,"Adaptive filter to use for clique bild pruning"),
    Arg("cpch",Arg::Opt,MaxClique::cliqueBeamContinuationHeuristic,"For clique beam build pruning, use a continuation heurisic within the rest of the clique"),

    Arg("cpef",Arg::Opt,MaxClique::cliqueBeamBuildExpansionFactor,"For clique beam build pruning, the amount that the beam expands each time we find a zero clique"),
    Arg("cpme",Arg::Opt,MaxClique::cliqueBeamBuildMaxExpansions,"For clique beam build pruning, the maximum number of expansions before we fail and exit"),


#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

    if (MaxClique::cliqueBeamBuildBeam == (-LZERO) ) {
      // then pruning is turned off, so set remaining params to sensibl values.
      MaxClique::cliqueBeamBuildMaxExpansions = 1;      
      MaxClique::cliqueBeamBuildExpansionFactor = 1.0;
    }

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
  Arg("cusample",Arg::Opt,MaxClique::cliqueBeamUniformSampleAmount,"Uniformly sample pruned clique (0<v=<=1 fraction, > 1 number)"),


#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

    if (MaxClique::cliqueBeamUniformSampleAmount < 0.0) {
      error("ERROR: -cusample option must be non-negative");
    }

#else
#endif
#endif // defined(GMTK_ARG_CKBEAM)


/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_CCBEAM)
#if defined(GMTK_ARGUMENTS_DEFINITION)

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("ccclusters",Arg::Opt,MaxClique::cliqueBeamClusterPruningNumClusters,"Number of clusters to use in cluster pruning"),
  Arg("ccbeam",Arg::Opt,MaxClique::cliqueBeamClusterBeam,"Clique cluster beam width to prune clique clusters (log value)"),

  Arg("cckbeam",Arg::Opt,MaxClique::cliqueBeamClusterMaxNumStates,"Max number of states in each cluster in cluster pruning"),

  Arg("ccrbeam",Arg::Opt,MaxClique::cliqueBeamClusterRetainFraction,"Fraction of in-cluster clique state space to retain. Range: 0 < v <= 1."),

  Arg("ccmbeam",Arg::Opt,MaxClique::cliqueBeamClusterMassRetainFraction,"Percentage of clique cluster mass to retain. Range: 0 < v <= 1. v = 1.0 means no pruning"),
  Arg("ccmexp",Arg::Opt,MaxClique::cliqueBeamClusterMassExponentiate,"Exponent to apply to clique cluster scores when doing mass pruning. Must be non-negative."),
  Arg("ccmmin",Arg::Opt,MaxClique::cliqueBeamClusterMassMinSize,"When using -cmbeam, min possible resulting clique cluster state size (>= 1)"),
  Arg("ccmfurther",Arg::Opt,MaxClique::cliqueBeamClusterMassFurtherBeam,"When using -ccmbeam, additional beam to use after mass has been acounted for (>= 0)"),



#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)


  if (MaxClique::cliqueBeamClusterMassRetainFraction <= 0.0 || MaxClique::cliqueBeamClusterMassRetainFraction > 1.0)
    error("%s: ccmbeam argument must be: 0.0 < v <= 1.0",argerr);
  if (MaxClique::cliqueBeamClusterMassMinSize <= 0)
    error("%s: -ccmmin option must be at least unity.",argerr);
  if (MaxClique::cliqueBeamClusterMassFurtherBeam < 0)
    error("%s: -ccmfurther option must be >= 0.",argerr);
  if (MaxClique::cliqueBeamClusterMassExponentiate < 0.0) 
    error("%s: -ccmexp option must be >= 0.",argerr);



#else
#endif
#endif // defined(GMTK_ARG_CCBEAM)



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
    error("%s: crbeam argument must be: 0.0 < v <= 1.0",argerr);

#else
#endif
#endif // defined(GMTK_ARG_CRBEAM)


/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/

#if defined(GMTK_ARG_CMBEAM)
#if defined(GMTK_ARGUMENTS_DEFINITION)

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("cmbeam",Arg::Opt,MaxClique::cliqueBeamMassRetainFraction,"Percentage of clique mass to relinquish. Range: 0 < v <= 1. v = 1.0 means no pruning"),
  Arg("cmexp",Arg::Opt,MaxClique::cliqueBeamMassExponentiate,"Exponent to apply to clique scores when doing mass pruning. Must be non-negative."),
  Arg("cmmin",Arg::Opt,MaxClique::cliqueBeamMassMinSize,"When using -cmbeam, min possible resulting clique state size (>= 1)"),
  Arg("cmfurther",Arg::Opt,MaxClique::cliqueBeamMassFurtherBeam,"When using -cmbeam, additional beam to use after mass has been acounted for (>= 0)"),


#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

  if (MaxClique::cliqueBeamMassRetainFraction <= 0.0 || MaxClique::cliqueBeamMassRetainFraction > 1.0)
    error("%s: cmbeam argument must be: 0.0 < v <= 1.0",argerr);
  if (MaxClique::cliqueBeamMassMinSize <= 0)
    error("%s: -cmmin option must be at least unity.",argerr);
  if (MaxClique::cliqueBeamMassFurtherBeam < 0)
    error("%s: -cmfurther option must be >= 0.",argerr);
  if (MaxClique::cliqueBeamMassExponentiate < 0.0) 
    error("%s: -cmexp option must be >= 0.",argerr);


#else
#endif
#endif // defined(GMTK_ARG_CMBEAM)



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
    error("%s: separatorBeam must be >= 0",argerr);

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

  // We could make this avalable to the command line, but this is really meant for an
  // internal variable. If desired, you can enable this by uncommenting.
  // Arg("minEMIncrementProb",Arg::Opt,EMable::minIncrementProbabilty.v,"Natural log of minumum EM increment posterior prob"),

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
      error("%s: hashLoadFactor must be between 0.05 and 1.0 non-inclusive",argerr);

#else
#endif
#endif // defined(GMTK_ARG_HASH_LOAD_FACTOR)



#if defined(GMTK_ARG_STORE_DETERMINISTIC_CHILDREN)
#if defined(GMTK_ARGUMENTS_DEFINITION)

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("deterministicChildrenStore",Arg::Opt,MaxClique::storeDeterministicChildrenInClique,"Store deterministic children in clique memory"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)


#else
#endif
#endif // defined(GMTK_ARG_STORE_DETERMINISTIC_CHILDREN)




/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/

#if defined(GMTK_ARG_CLEAR_CLIQUE_VAL_MEM)
#if defined(GMTK_ARGUMENTS_DEFINITION)

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("clearCliqueValMem",Arg::Opt,JunctionTree::perSegmentClearCliqueValueCache,"Free clique/separator value cache for each segment"),

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

  const static char *dcdrng_str="all";

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

  static const char *trrng_str="all";

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
    error("%s: arguments startSkip=%d/endSkip=%d must both be >= 0",argerr,startSkip,endSkip);

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

#ifdef GMTK_ARG_VERB_DEF_VAL
  static unsigned verbosity = GMTK_ARG_VERB_DEF_VAL;
#else
  static unsigned verbosity = IM::Default;
#endif

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("verbosity",Arg::Opt,verbosity,"Verbosity (0 <= v <= 100) of informational/debugging msgs"),
  Arg("printIntValues",Arg::Opt,RV::alwaysPrintIntegerRVValues,"always print rv values as integer rather than symbols"),

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
#ifdef HAVE_CONFIG_H
    printf("%s (Mercurial id: %s)\n",gmtk_version_id,HGID);
#else
    printf("%s\n", gmtk_version_id);
#endif
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
static bool  cliquePrintOnlyEntropy = false;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("pCliquePrintRange",Arg::Opt,pPartCliquePrintRange,"With CE/DE, print range cliques from P partition."),
  Arg("cCliquePrintRange",Arg::Opt,cPartCliquePrintRange,"With CE/DE, print range cliques from C partition."),
  Arg("eCliquePrintRange",Arg::Opt,ePartCliquePrintRange,"With CE/DE, print range cliques from E partition."),
  Arg("cliquePrintOnlyEntropy",Arg::Opt,cliquePrintOnlyEntropy,"With CE/DE, print only clique entropies."),

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
  if (island) {
    infoMsg(IM::Default,"NOTE: running island algorithm, turning off component caching '-componentCache F', setting hash load factor to at least 0.98 '-hashLoadFactor 0.98', and not storing deterministic children '-deterministicChildrenStore F'\n"); 
    fflush(stdout);
    MixtureCommon::cacheMixtureProbabilities = false;
    // make sure to use other low memory options.
    if (hash_abstract::loadFactor < 0.98)
      hash_abstract::loadFactor = 0.98;
    MaxClique::storeDeterministicChildrenInClique = false;
  }


#else
#endif
#endif // defined(GMTK_ARG_ISLAND)

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

static const char* varPartitionAssignmentPrior = "COI";
static const char* varCliqueAssignmentPrior = "COT";

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
      "Use Virtual Evidence (VE) Separators (if any are available) during inference (Bitwise or of 0x1 (PC) or PCG (0x2)"),
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
  Arg("localCliqueNorm",Arg::Opt,localCliqueNormalization,"Use local clique sum for EM posterior normalization."),
  Arg("dirichletPriors",Arg::Opt,EMable::useDirichletPriors,"Enable the use of Dirichlet priors for this process."),

  Arg("gmarCoeffL2",Arg::Opt,GaussianComponent::gmarCoeffL2,"Gaussian mean l2 accuracy-regularization tradeoff coeff (ie, prior concentration)"),
  Arg("gdarCoeffL2",Arg::Opt,GaussianComponent::gdarCoeffL2,"Gaussian dlink l2 accuracy-regularziation tradeoff coeff (ie, prior concentration)"),


#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

  MixtureCommon::checkForValidRatioValues();
  MeanVector::checkForValidValues();
  DiagCovarVector::checkForValidValues();
  DlinkMatrix::checkForValidValues();


#else
#endif
#endif // defined(GMTK_ARG_EM_TRAINING_PARAMS)



/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/

/* this next argument is applicable to all inference, but we add it 
 * here as an error check needs to be done only in the EM and/or training
 * case.
 */

#if defined(GMTK_ARG_CLIQUE_TABLE_NORMALIZE)
#if defined(GMTK_ARGUMENTS_DEFINITION)

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("cliqueTableNormalize",Arg::Opt,MaxClique::normalizeScoreEachClique,"Normalize scores of each clique right after its creation (increases dynamic range)."),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

    if (MaxClique::normalizeScoreEachClique < 0.0) {
      error("ERROR: -cliqueTableNormalize option must be non-negative\n");
    }


#if defined(GMTK_ARG_EM_TRAINING_PARAMS)
    if (MaxClique::normalizeScoreEachClique != 0.0 && localCliqueNormalization == false) {
      // EM training won't work in this case unless it does local clique normalization as well.
      localCliqueNormalization = true;
      infoMsg(IM::SoftWarning,"Turning on EM local clique normalization since clique table score normalization is on.\n");
    }
#endif


#else
#endif
#endif // defined(GMTK_ARG_CLIQUE_TABLE_NORMALIZE)




/*==============================================================================================================*/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************                                     ***********************************************/
/****************************       GMTK KERNEL OPTIONS           ***********************************************/
/****************************                                     ***********************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_KERNEL_PARAMS)
#if defined(GMTK_ARGUMENTS_DEFINITION)

static bool randomizeParams = false;
static bool transFileIsBinary = false;
static char *objsToNotUtilizeFile=NULL;
static bool localCliqueNormalization = false;

// WARNING: we set the default behaivor of any program that includes the kernel code to true (to default to fisher kernel)
// but this means that we should not include this in any program that wants to do EM training by default.
static bool fisherKernelP = true;
static char *storeFeatureFile = NULL;
static bool annotateTransformationOutput = true;
static bool writeLogVals = false;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("random",Arg::Opt,randomizeParams,"Randomize the parameters"),
  // kernel accumulator file support
  Arg("transFileIsBinary",Arg::Opt,transFileIsBinary,"Use binary to write the parameters in transformed space"), 
  Arg("objsNotToUtilize",Arg::Opt,objsToNotUtilizeFile,"File listing trainable parameter objects to not utilize in transformed space."),
  Arg("storeFeatureFile",Arg::Req,storeFeatureFile,"File to store feature space values in"),
  Arg("localCliqueNorm",Arg::Opt,localCliqueNormalization,"Use local clique sum for posterior normalization."),
  Arg("fisherKernel",Arg::Opt,fisherKernelP,"Compute the fisher kernel transformation."),
  Arg("annotateTransformationOutput",Arg::Opt,annotateTransformationOutput,"Annotation the output tranformation matrix (one per line)."),
  Arg("writeLogVals",Arg::Opt,writeLogVals,"Write log(p) rather than just p for prob values."),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)



#else
#endif
#endif // defined(GMTK_ARG_KERNEL_PARAMS)




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

const char*    Action_If_Diff_Num_Frames_Str[MAX_NUM_OBS_FILES]={"er","er","er","er","er"};   // 
unsigned Action_If_Diff_Num_Frames[MAX_NUM_OBS_FILES]={FRAMEMATCH_ERROR,FRAMEMATCH_ERROR,FRAMEMATCH_ERROR,FRAMEMATCH_ERROR,FRAMEMATCH_ERROR};   // 
const char*    Action_If_Diff_Num_Sents_Str[MAX_NUM_OBS_FILES]={"te","te","te","te","te"}; 
unsigned Action_If_Diff_Num_Sents[MAX_NUM_OBS_FILES]={SEGMATCH_TRUNCATE_FROM_END,SEGMATCH_TRUNCATE_FROM_END,SEGMATCH_TRUNCATE_FROM_END,SEGMATCH_TRUNCATE_FROM_END,SEGMATCH_TRUNCATE_FROM_END};   // 

char    *Per_Stream_Transforms[MAX_NUM_OBS_FILES]={NULL,NULL,NULL,NULL,NULL};   // 
char    *Post_Transforms=NULL;

const char    *Ftr_Combo_Str="none";
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
  else error("%s: Unknown feature combination type: '%s'\n",argerr,Ftr_Combo_Str);
  
  for(int i=0; i < MAX_NUM_OBS_FILES; ++i) {
    if(ofs[i]!=NULL) {
      if (strcmp(Action_If_Diff_Num_Frames_Str[i],"er") == 0)      Action_If_Diff_Num_Frames[i] = FRAMEMATCH_ERROR;
      else if (strcmp(Action_If_Diff_Num_Frames_Str[i],"rl") == 0) Action_If_Diff_Num_Frames[i] = FRAMEMATCH_REPEAT_LAST;
      else if (strcmp(Action_If_Diff_Num_Frames_Str[i],"rf") == 0) Action_If_Diff_Num_Frames[i] = FRAMEMATCH_REPEAT_FIRST;
      else if (strcmp(Action_If_Diff_Num_Frames_Str[i],"se") == 0) Action_If_Diff_Num_Frames[i] = FRAMEMATCH_EXPAND_SEGMENTALLY;
      else if (strcmp(Action_If_Diff_Num_Frames_Str[i],"ts") == 0) Action_If_Diff_Num_Frames[i] = FRAMEMATCH_TRUNCATE_FROM_START;
      else if (strcmp(Action_If_Diff_Num_Frames_Str[i],"te") == 0) Action_If_Diff_Num_Frames[i] = FRAMEMATCH_TRUNCATE_FROM_END;
      else error("%s: Unknown action when diff num of frames: '%s'\n",argerr,Action_If_Diff_Num_Frames_Str[i]);
    }
  }
  
  for(int i=0; i < MAX_NUM_OBS_FILES; ++i) {
    if(ofs[i]!=NULL) {
      if (strcmp(Action_If_Diff_Num_Sents_Str[i],"er") == 0)      Action_If_Diff_Num_Sents[i] = SEGMATCH_ERROR;
      else if (strcmp(Action_If_Diff_Num_Sents_Str[i],"rl") == 0) Action_If_Diff_Num_Sents[i] = SEGMATCH_REPEAT_LAST;
      else if (strcmp(Action_If_Diff_Num_Sents_Str[i],"wa") == 0) Action_If_Diff_Num_Sents[i] = SEGMATCH_WRAP_AROUND;
      else if (strcmp(Action_If_Diff_Num_Sents_Str[i],"te") == 0) Action_If_Diff_Num_Sents[i] = SEGMATCH_TRUNCATE_FROM_END;
      else error("%s: Unknown action when diff num of sentences: '%s'\n",argerr,
		 Action_If_Diff_Num_Sents_Str[i]);
    }
  }


#else
#endif
#endif // defined(GMTK_ARG_OBS_MATRIX_XFORMATION)


/*==============================================================================================================*/
/****************************************************************************************************************/
/****************************************************************************************************************/
/************************                                              ******************************************/
/************************            DECODING OPTIONS                  ******************************************/
/************************                                              ******************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/


/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_DECODING_OPTIONS)
#if defined(GMTK_ARGUMENTS_DEFINITION)

  static char *dumpNames = NULL;
  static char *ofilelist = NULL;
  static char *wordVar=NULL;
  static char *varMapFile=NULL;
  static char *transitionLabel=NULL;
  static char* showVitVals = NULL;


#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("dumpNames",Arg::Opt,dumpNames,"File containing the names of the variables to save to a file"),
  Arg("ofilelist",Arg::Opt,ofilelist,"List of filenames to dump the hidden variable values to"),
  // These 3 must be used together or not at all
  Arg("printWordVar",Arg::Opt,wordVar,"Print the word var - which has this label"),
  Arg("varMap",Arg::Opt,varMapFile,"Use this file to map from word-index to string"),
  Arg("transitionLabel",Arg::Opt,transitionLabel,"The label of the word transition variable"),
  Arg("showVitVals",Arg::Opt,showVitVals,"File to print viterbi values, '-' for stdout"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

  if (dumpNames)
    if (ofilelist==NULL) 
      error("%s: Must also specify output files for binary writing",argerr);

#else
#endif
#endif // defined(GMTK_ARG_DECODING_OPTIONS)

// Options for the transition to the new gmtkViterbi front end.

#if defined(GMTK_ARG_NEW_DECODING_OPTIONS)
#if defined(GMTK_ARGUMENTS_DEFINITION)


  // arguments for partition based Viterbi printing.
  static char* pVitValsFileName = NULL;
  // TODO: get binary printing working
  // static bool pVitValsFileBinp = false;
  static char* pVitRegexFilter = NULL;
  static bool pVitCaseSensitiveRegexFilter = false;
  static char* pVitPartRangeFilter = NULL;
  static bool pVitAlsoPrintObservedVariables = false;

  // arguments for frame-based Viterbi printing.
//  static char* vitValsFileName = NULL;
  // TODO: get binary printing working
  // static bool vitValsFileBinp = false;
  static char* vitRegexFilter = NULL;
  static bool vitCaseSensitiveRegexFilter = false;
//  static bool vitAlsoPrintObservedVariables = false;
//  static bool vitReverseOrder = false;
#define MAX_VITERBI_TRIGGERS 3
  const char   *vitTriggerVariables[MAX_VITERBI_TRIGGERS] = { NULL, NULL, NULL };
  const char   *vitTriggerSets[MAX_VITERBI_TRIGGERS] = { NULL, NULL, NULL };



#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  // partition based 
  Arg("pVitValsFile",Arg::Req,pVitValsFileName,"Partition Vit: file to print viterbi values, '-' for stdout"),
  // TODO: not currently used, but should add.
  // Arg("pVitBinVitValsFile",Arg::Opt,pVitValsFileBinp,"Partition Vit: Should file to print viterbi values be binary? (T/F) "),

  Arg("pVitRegexFilter",Arg::Opt,pVitRegexFilter,"Partition Vit: Regular expression to filter variable names."),
  Arg("pVitCaseSensitiveRegexFilter",Arg::Opt,pVitCaseSensitiveRegexFilter,"Partition Vit: Case sensitivity of the rv regular expression filter."),

  Arg("pVitPrintRange",Arg::Opt,pVitPartRangeFilter,"Partition Vit: value printing, integer range filter for partitions (e.g., frames, slices) to print."),

  Arg("pVitPrintObservedVariables",Arg::Opt,pVitAlsoPrintObservedVariables,"Partition Vit: also print observed random variables in addtion to hidden"),

#if 0
  // this is not implemented yet.
  // frame based
  Arg("vitValsFile",Arg::Req,vitValsFileName,"Vit: file to print viterbi values, '-' for stdout"),
  // TODO: not currently used, but should add.
  // Arg("vitBinVitValsFile",Arg::Opt,vitValsFileBinp,"Vit: Should file to print viterbi values be binary? (T/F) "),

  Arg("vitRegexFilter",Arg::Opt,vitRegexFilter,"Vit: Regular expression to filter variable names."),
  Arg("vitCaseSensitiveRegexFilter",Arg::Opt,vitCaseSensitiveRegexFilter,"Vit: Case sensitivity of the rv regular expression filter."),

  Arg("vitPrintObservedVariables",Arg::Opt,vitAlsoPrintObservedVariables,"Vit: also print observed random variables in addtion to hidden"),

  Arg("vitReverseOrder",Arg::Opt,vitReverseOrder,"Vit: print values in reverse order."),

  Arg("vitTriggerVar", Arg::Opt,vitTriggerVariables,"Viterbi: Trigger variable. Replace X with trigger variable number",Arg::ARRAY,MAX_VITERBI_TRIGGERS),
  Arg("vitTriggerSet", Arg::Opt,vitTriggerSets,"Viterbi: Trigger value set. Replace X with trigger set number",Arg::ARRAY,MAX_VITERBI_TRIGGERS),
#endif

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)


#else
#endif
#endif // defined(GMTK_ARG_NEW_DECODING_OPTIONS)


/*==============================================================================================================*/
/****************************************************************************************************************/
/****************************************************************************************************************/
/************************                                              ******************************************/
/************************            TIMING OPTIONS                    ******************************************/
/************************                                              ******************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/



/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/




#if defined(GMTK_ARG_TIMING)
#if defined(GMTK_ARGUMENTS_DEFINITION)

  static unsigned seconds = 10;
  static bool noEPartition = false;
  static unsigned numTimes = 1;
  static bool multiTest = false;
  static int rlimitSlop = 2;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("seconds",Arg::Opt,seconds,"Number of seconds to run and then exit."),
  Arg("times",Arg::Opt,numTimes,"Number of times to run program seconds seconds long (not multitest mode)."),
  Arg("multiTest",Arg::Opt,multiTest,"Run gmtkTime in multi-test mode, taking triangulation file names from command line."),
  Arg("slop",Arg::Opt,rlimitSlop,"In multiTest mode, number of additional seconds before fail-terminate is forced."),
  Arg("noEPartition",Arg::Opt,noEPartition,"If true, do not run E partition (only [P C C ... C] skipping E)"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#else
#endif
#endif // defined(GMTK_ARG_TIMING)


/*==============================================================================================================*/
/****************************************************************************************************************/
/****************************************************************************************************************/
/************************                                              ******************************************/
/************************            TRIANGULATION OPTIONS             ******************************************/
/************************                                              ******************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/


/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/



#if defined(GMTK_ARG_TRIANGULATION_OPTIONS)
#if defined(GMTK_ARGUMENTS_DEFINITION)

  const static char* triangulationHeuristic="completed";
  static bool jtWeight = true;
  static double traverseFraction = 1.0;
  static bool noBoundaryMemoize = false;
  const static char* forceLeftRight="";
  const static char* boundaryHeuristic="S";
  static unsigned maxNumChunksInBoundary = 1; 
  static unsigned chunkSkip = 1; 
  static int jut = -1;
  static char* anyTimeTriangulate = NULL;
  static char* timeLimit = NULL;
  static bool rePartition = false;
  static bool reTriangulate = false;
  static bool continueTriangulating = false;
  static bool noReTriP = false;
  static bool noReTriC = false;
  static bool noReTriE = false;
  static bool printResults = false;

  static bool longStrCheck = false;
  static bool findBestBoundary = true;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("triangulationHeuristic",
      Arg::Opt,triangulationHeuristic,
      "Triang. heuristic, >1 of S=size,T=time,F=fill,W=wght,X=rev-time,P=pos,H=hint,R=rnd,N=wght-w/o-det"),

  Arg("jtWeight",
      Arg::Opt,jtWeight,
      "True means use an estimate of the JT weight to score triangulation rather than sum of weight"),


  Arg("jtwPUI",
      Arg::Opt,JunctionTree::jtWeightPenalizeUnassignedIterated,
      "Amount jtWeight should penalize cliques with unassigned iterated nodes (0.0 means no penalty)"),

  Arg("jtwMC",
      Arg::Opt,JunctionTree::jtWeightMoreConservative,
      "True means jtWeight should be more conservative (more upper bound like) regarding charges to some nodes"),

  Arg("jtwSNSC",
      Arg::Opt,JunctionTree::jtWeightSparseNodeSepScale,
      "Amount to scale charge of a sparse node in a clique's incomming separator"),

  Arg("jtwDNSC",
      Arg::Opt,JunctionTree::jtWeightDenseNodeSepScale,
      "Amount to scale charge of a dense node in a clique's incomming separator"),

  Arg("pfCobWeight",
      Arg::Opt,MaxClique::continuousObservationPerFeaturePenalty,
      "Per-Feature Dimension Continuous Observation Log penalty to use in clique weight calc"),

  Arg("findBestBoundary",
      Arg::Opt,findBestBoundary,
      "Run the (exponential time) boundary algorithm or not."),

  Arg("traverseFraction",
      Arg::Opt,traverseFraction,
      "Fraction of current interface to traverse in boundary recursion."),

  Arg("noBoundaryMemoize",
      Arg::Opt,noBoundaryMemoize,
      "Do not memoize boundaries (less memory but runs slower)"),

  Arg("forceLeftRight",
      Arg::Opt,forceLeftRight,
      "Run boundary algorithm only for either left (L) or right (R) interface, rather than both"),

  Arg("boundaryHeuristic",
      Arg::Opt,boundaryHeuristic,
      "Boundary heuristic, >1 of S=size,F=fill,W=wght,N=wght-w/o-det,M=max-clique,C=max-C-clique,A=st-spc,Q=C-st-spc"),

  Arg("M",
      Arg::Opt,maxNumChunksInBoundary,
      "Max number simultaneous chunks in which boundary may simultaneously exist"),

  Arg("S",
      Arg::Opt,chunkSkip,
      "Number of chunks that should exist between boundaries"),

  Arg("disconnectFromObservedParent",
      Arg::Opt,RV::disconnectChildrenOfObservedParents,
      "In going to UGM, disconnect children from observed parents when possible"),


  Arg("unroll",
      Arg::Opt,jut,
      "Unroll graph & triangulate using heuristics. DON'T use P,C,E constrained triangulation."),

  Arg("anyTimeTriangulate",
      Arg::Opt,anyTimeTriangulate,
      "Run the any-time triangulation algorithm for given duration."),

  Arg("timeLimit",
      Arg::Opt,timeLimit,
      "Do not run for longer than the given amount of time."),

  Arg("rePartition",
      Arg::Opt,rePartition,
      "Re-Run the boundary algorithm even if .str.trifile exists to produce new partition and new triangulation."),

  Arg("reTriangulate",
      Arg::Opt,reTriangulate,
      "Re-Run only triangluation using existing partition given in .trifile."),

  Arg("continueTriangulating",
      Arg::Opt,continueTriangulating,
      "When re-triangulating existing .tri file, continue besting existing triangulations"),

  Arg("noReTriP",
      Arg::Opt,noReTriP,
      "When re-triangulating existing .tri file, don't re-triangulate P, keep old"),
  Arg("noReTriC",
      Arg::Opt,noReTriC,
      "When re-triangulating existing .tri file, don't re-triangulate C, keep old"),
  Arg("noReTriE",
      Arg::Opt,noReTriE,
      "When re-triangulating existing .tri file, don't re-triangulate E, keep old"),

  Arg("printResults",Arg::Opt,printResults,"Print information about result of triangulation."),

  Arg("longStrCheck",Arg::Opt,longStrCheck,"Set to true to do the long check for structure file validity"),


#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

  if (chunkSkip < 1)
    error("%s: chunk skip parameter S must be >= 1\n",argerr);
  if (maxNumChunksInBoundary < 1)
    error("%s: max number chunks in boundary parameter M must be >= 1\n",argerr);
  if (fabs(MaxClique::continuousObservationPerFeaturePenalty) > 1.0) {
    infoMsg(IM::Warning,"###\n### !!!DANGER WILL ROBINSON!! LARGE -pfCobWeight VALUE %f MIGHT CAUSE FLOATING POINT EXCEPTION. SUGGEST REDUCE IT IF FPE OCCURS!! ###\n###\n",MaxClique::continuousObservationPerFeaturePenalty);
  }


#else
#endif
#endif // defined(GMTK_ARG_TRIANGULATION_OPTIONS)



/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_LOAD_PARAMETERS)
#if defined(GMTK_ARGUMENTS_DEFINITION)

  static bool loadParameters = false;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("loadParameters",Arg::Opt,loadParameters,"Also load in all trainable parameters."),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#else
#endif
#endif // defined(GMTK_ARG_LOAD_PARAMETERS)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_NUM_BACKUP_FILES)
#if defined(GMTK_ARGUMENTS_DEFINITION)

  static unsigned numBackupFiles = 10;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("numBackupFiles",Arg::Opt,numBackupFiles,"Number of backup .trifiles (_bak0,_bak1,etc.) to keep."),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)


#else
#endif
#endif // defined(GMTK_ARG_NUM_BACKUP_FILES)



/*==============================================================================================================*/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************                                     ***********************************************/
/****************************         TYING OPTIONS               ***********************************************/
/****************************                                     ***********************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_TYING_PARAMS)
#if defined(GMTK_ARGUMENTS_DEFINITION)

static char *loadAccFile = NULL;
static char *loadAccRange = NULL;
static bool accFileIsBinary = true;
static char *loadCmdFile = NULL;


#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  // EM accumulator file support - loading only
  Arg("loadAccFile",Arg::Opt,loadAccFile,"Load accumulators file"), 
  Arg("loadAccRange",Arg::Opt,loadAccRange,"Load accumulators file range"), 
  Arg("accFileIsBinary",Arg::Opt,accFileIsBinary,"Binary accumulator files"), 

  Arg("loadCmdFile",Arg::Opt,loadCmdFile,"Load tying command file"), 




#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

  MixtureCommon::checkForValidRatioValues();
  MeanVector::checkForValidValues();
  DiagCovarVector::checkForValidValues();
  DlinkMatrix::checkForValidValues();



#else
#endif
#endif // defined(GMTK_ARG_TYING_PARAMS)



/*==============================================================================================================*/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************                                     ***********************************************/
/****************************         LATTICE OPTIONS             ***********************************************/
/****************************                                     ***********************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_LATTICE_PARAMS)
#include "GMTK_LatticeADT.h"

#if defined(GMTK_ARGUMENTS_DEFINITION)


#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("latticeDefaultFrameRate",Arg::Opt,LatticeADT::_defaultFrameRate,"Lattice, default frame rate (if negative, compute from file)"),
  Arg("latticeUseMaxScore",Arg::Opt,LatticeADT::_latticeNodeUseMaxScore,"Lattice, use max edge score for node CPT"),
  Arg("latticeIgnoreNodeTimeMarks",Arg::Opt,LatticeADT::_ignoreLatticeNodeTimeMarks,"Lattice, ignore lattice node time marks in all lattices"),


#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)


#else
#endif
#endif // defined(GMTK_ARG_LATTICE_PARAMS)




/*==============================================================================================================*/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************                                     ***********************************************/
/****************************     RESOURCE LIMITNG OPTIONS        ***********************************************/
/****************************                                     ***********************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/


/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_RLIMIT_PARAMS)
#if defined(GMTK_ARGUMENTS_DEFINITION)
  static int rlimitMaxMem = 0;
  static int rlimitMaxTime = 0;
  static int rlimitMaxCore = -1;


#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)
  Arg("maxMem",Arg::Opt,rlimitMaxMem,  "Maximum virtual memory  (Bytes); 0 means 'unlimited'"), 
  Arg("maxTime",Arg::Opt,rlimitMaxTime,"Maximum CPU time      (seconds); 0 means 'unlimited'"), 
  Arg("maxCore",Arg::Opt,rlimitMaxCore,"Maximum core file size  (Bytes);-1 means 'unlimited'"), 


#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
    if (rlimitMaxMem > 0) {
      infoMsg(IM::Tiny,"Setting max allowed virtual memory\n");
      struct rlimit rl;
      rl.rlim_cur = rlimitMaxMem;
      rl.rlim_max = rlimitMaxMem;
      if (setrlimit(RLIMIT_AS,&rl) != 0)
	error("Error: tried to set max memory beyond the allowed maximum");
    }

    if (rlimitMaxTime > 0) {
      infoMsg(IM::Tiny,"Setting max allowed CPU time\n");
      struct rlimit rl;

      // soft limit
      rl.rlim_cur = rlimitMaxTime;

      // hard limit must be higher than soft limit, otherwise an
      // immediate SIGKILL is sent which doesn't allow an orderly
      // exit; so, simply add 10 seconds to allow time for the SIGXCPU
      // handler to be called
      rl.rlim_max = rlimitMaxTime + 10; 
      if (setrlimit(RLIMIT_CPU,&rl) != 0)
	error("Error: tried to set max CPU time beyond the allowed maximum");
    };

    if (rlimitMaxCore >= 0) {
      infoMsg(IM::Tiny,"Setting max allowed core file size\n");
      struct rlimit rl;
      rl.rlim_cur = rlimitMaxCore;
      rl.rlim_max = rlimitMaxCore;
      if (setrlimit(RLIMIT_CORE,&rl) != 0)
	error("Error: tried to set max core file size beyond the allowed maximum");
    }

#else
#endif
#endif // defined(GMTK_ARG_RLIMIT_PARAMS)

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

