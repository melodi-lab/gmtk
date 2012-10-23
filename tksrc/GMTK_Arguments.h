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

#include "GMTK_Filter.h"
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


#include "GMTK_ObservationArguments.h"



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

#if defined(GMTK_ARG_INPUT_TRAINABLE_FILE_HANDLING)
#if defined(GMTK_ARGUMENTS_DOCUMENTATION)
  Arg("\n*** Input trainable parameter file handling ***\n"),
#endif
#endif


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

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#else
#endif
#endif // defined(GMTK_ARG_INPUT_MASTER_FILE)


/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_DLOPEN_MAPPERS)
#if defined(GMTK_ARGUMENTS_DEFINITION)

#include <vector>
#include <dlfcn.h>

#define MAX_NUM_DLOPENED_FILES (5)
   char *dlopenFilenames[MAX_NUM_DLOPENED_FILES] = {NULL,NULL,NULL,NULL,NULL};

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

Arg("map",Arg::Opt,dlopenFilenames,"Deterministic mapping dynamic library file. Replace X with the file number",Arg::ARRAY,MAX_NUM_DLOPENED_FILES),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#else
#endif
#endif // defined(GMTK_ARG_DLOPEN_MAPPERS)


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

#if defined(GMTK_ARG_INPUT_MODEL_FILE_HANDLING)
#if defined(GMTK_ARGUMENTS_DOCUMENTATION)
  Arg("\n*** Input model file handling ***\n"),
#endif
#endif


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

#if defined(GMTK_ARG_CONTINUOUS_RANDOM_VAR_OPTIONS)
#if defined(GMTK_ARGUMENTS_DOCUMENTATION)
  Arg("\n*** Continuous random variable options ***\n"),
#endif
#endif


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

#if defined(GMTK_ARG_BEAM_PRUNING_OPTIONS)
#if defined(GMTK_ARGUMENTS_DOCUMENTATION)
  Arg("\n*** Beam pruning options ***\n"),
#endif
#endif

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

#if defined(GMTK_ARG_MEMORY_MANAGEMENT_OPTIONS)
#if defined(GMTK_ARGUMENTS_DOCUMENTATION)
  Arg("\n*** Memory management options ***\n"),
#endif
#endif

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_MEM_GROWTH)
#if defined(GMTK_ARGUMENTS_DEFINITION)

  static char const *memGrowthOption = "default";

#define GMTK_MEM_GROWTH_CONSERVATIVE 0
#define GMTK_MEM_GROWTH_DEFAULT      1
#define GMTK_MEM_GROWTH_AGGRESSIVE   2

  static unsigned memGrowthStrategy = GMTK_MEM_GROWTH_DEFAULT;

#define GMTK_MEM_CONSERVATIVE_START_SIZE  1
#define GMTK_MEM_CONSERVATIVE_GROWTH_RATE 1.05
#define GMTK_MEM_CONSERVATIVE_DECAY_RATE 0.0

#define GMTK_MEM_DEFAULT_START_SIZE  23
#define GMTK_MEM_DEFAULT_GROWTH_RATE 1.25
#define GMTK_MEM_DEFAULT_DECAY_RATE 0.0

#define GMTK_MEM_AGGRESSIVE_START_SIZE  23
#define GMTK_MEM_AGGRESSIVE_GROWTH_RATE 2.0
#define GMTK_MEM_AGGRESSIVE_DECAY_RATE 0.0


#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("memoryGrowth",Arg::Opt,memGrowthOption,"Rate to grow data structures (conservative, default, aggressive)"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

  if (strncasecmp(memGrowthOption, "conservative", 13) == 0) {
    infoMsg(IM::Default,"NOTE: using conservative memory strategy - turning off component caching '-componentCache F', setting hash load factor to at least 0.98 '-hashLoadFactor 0.98', and not storing deterministic children '-deterministicChildrenStore F'\n"); 
    fflush(stdout);
    MixtureCommon::cacheMixtureProbabilities = false;
    // make sure to use other low memory options.
    if (hash_abstract::loadFactor < 0.98)
      hash_abstract::loadFactor = 0.98;
    MaxClique::storeDeterministicChildrenInClique = false;

    memGrowthStrategy = GMTK_MEM_GROWTH_CONSERVATIVE;
    
    CliqueValueHolder::defaultAllocationUnitChunkSize = GMTK_MEM_CONSERVATIVE_START_SIZE;
    CliqueValueHolder::defaultGrowthFactor            = GMTK_MEM_CONSERVATIVE_GROWTH_RATE;
    
    SeparatorClique::aiStartingSize                   = GMTK_MEM_CONSERVATIVE_START_SIZE;
    SeparatorClique::aiGrowthFactor                   = GMTK_MEM_CONSERVATIVE_GROWTH_RATE;
    
    SeparatorClique::remStartingSize                  = GMTK_MEM_CONSERVATIVE_START_SIZE;
    SeparatorClique::remGrowthFactor                  = GMTK_MEM_CONSERVATIVE_GROWTH_RATE;
    
    SeparatorClique::sepSpaceMgrStartingSize          = GMTK_MEM_CONSERVATIVE_START_SIZE;
    SeparatorClique::sepSpaceMgrGrowthRate            = GMTK_MEM_CONSERVATIVE_GROWTH_RATE;
    SeparatorClique::sepSpaceMgrDecayRate             = GMTK_MEM_CONSERVATIVE_DECAY_RATE;
    
    SeparatorClique::remSpaceMgrStartingSize          = GMTK_MEM_CONSERVATIVE_START_SIZE;
    SeparatorClique::remSpaceMgrGrowthRate            = GMTK_MEM_CONSERVATIVE_GROWTH_RATE;
    SeparatorClique::remSpaceMgrDecayRate             = GMTK_MEM_CONSERVATIVE_DECAY_RATE;
    
    ConditionalSeparatorTable::remHashMapStartingSize = GMTK_MEM_CONSERVATIVE_START_SIZE;
    
    MaxClique::spaceMgrStartingSize                   = GMTK_MEM_CONSERVATIVE_START_SIZE;
    MaxClique::spaceMgrGrowthRate                     = GMTK_MEM_CONSERVATIVE_GROWTH_RATE;
    MaxClique::spaceMgrDecayRate                      = GMTK_MEM_CONSERVATIVE_DECAY_RATE;
  
    MaxCliqueTable::valuePoolGrowthRate               = GMTK_MEM_CONSERVATIVE_GROWTH_RATE;

  } else if (strncasecmp(memGrowthOption, "default", 8) == 0) {
    memGrowthStrategy = GMTK_MEM_GROWTH_DEFAULT;

    CliqueValueHolder::defaultAllocationUnitChunkSize = GMTK_MEM_DEFAULT_START_SIZE;
    CliqueValueHolder::defaultGrowthFactor            = GMTK_MEM_DEFAULT_GROWTH_RATE;
    
    SeparatorClique::aiStartingSize                   = GMTK_MEM_DEFAULT_START_SIZE;
    SeparatorClique::aiGrowthFactor                   = GMTK_MEM_DEFAULT_GROWTH_RATE;
    
    SeparatorClique::remStartingSize                  = GMTK_MEM_DEFAULT_START_SIZE;
    SeparatorClique::remGrowthFactor                  = GMTK_MEM_DEFAULT_GROWTH_RATE;
    
    SeparatorClique::sepSpaceMgrStartingSize          = GMTK_MEM_DEFAULT_START_SIZE;
    SeparatorClique::sepSpaceMgrGrowthRate            = GMTK_MEM_DEFAULT_GROWTH_RATE;
    SeparatorClique::sepSpaceMgrDecayRate             = GMTK_MEM_DEFAULT_DECAY_RATE;
    
    SeparatorClique::remSpaceMgrStartingSize          = GMTK_MEM_DEFAULT_START_SIZE;
    SeparatorClique::remSpaceMgrGrowthRate            = GMTK_MEM_DEFAULT_GROWTH_RATE;
    SeparatorClique::remSpaceMgrDecayRate             = GMTK_MEM_DEFAULT_DECAY_RATE;
    
    ConditionalSeparatorTable::remHashMapStartingSize = GMTK_MEM_DEFAULT_START_SIZE;
    
    MaxClique::spaceMgrStartingSize                   = GMTK_MEM_DEFAULT_START_SIZE;
    MaxClique::spaceMgrGrowthRate                     = GMTK_MEM_DEFAULT_GROWTH_RATE;
    MaxClique::spaceMgrDecayRate                      = GMTK_MEM_DEFAULT_DECAY_RATE;
  
    MaxCliqueTable::valuePoolGrowthRate               = GMTK_MEM_DEFAULT_GROWTH_RATE;

  } else if (strncasecmp(memGrowthOption, "aggressive", 11) == 0) {
    memGrowthStrategy = GMTK_MEM_GROWTH_AGGRESSIVE;

    CliqueValueHolder::defaultAllocationUnitChunkSize = GMTK_MEM_AGGRESSIVE_START_SIZE;
    CliqueValueHolder::defaultGrowthFactor            = GMTK_MEM_AGGRESSIVE_GROWTH_RATE;
    
    SeparatorClique::aiStartingSize                   = GMTK_MEM_AGGRESSIVE_START_SIZE;
    SeparatorClique::aiGrowthFactor                   = GMTK_MEM_AGGRESSIVE_GROWTH_RATE;
    
    SeparatorClique::remStartingSize                  = GMTK_MEM_AGGRESSIVE_START_SIZE;
    SeparatorClique::remGrowthFactor                  = GMTK_MEM_AGGRESSIVE_GROWTH_RATE;
    
    SeparatorClique::sepSpaceMgrStartingSize          = GMTK_MEM_AGGRESSIVE_START_SIZE;
    SeparatorClique::sepSpaceMgrGrowthRate            = GMTK_MEM_AGGRESSIVE_GROWTH_RATE;
    SeparatorClique::sepSpaceMgrDecayRate             = GMTK_MEM_AGGRESSIVE_DECAY_RATE;
    
    SeparatorClique::remSpaceMgrStartingSize          = GMTK_MEM_AGGRESSIVE_START_SIZE;
    SeparatorClique::remSpaceMgrGrowthRate            = GMTK_MEM_AGGRESSIVE_GROWTH_RATE;
    SeparatorClique::remSpaceMgrDecayRate             = GMTK_MEM_AGGRESSIVE_DECAY_RATE;
    
    ConditionalSeparatorTable::remHashMapStartingSize = GMTK_MEM_AGGRESSIVE_START_SIZE;
    
    MaxClique::spaceMgrStartingSize                   = GMTK_MEM_AGGRESSIVE_START_SIZE;
    MaxClique::spaceMgrGrowthRate                     = GMTK_MEM_AGGRESSIVE_GROWTH_RATE;
    MaxClique::spaceMgrDecayRate                      = GMTK_MEM_AGGRESSIVE_DECAY_RATE;
  
    MaxCliqueTable::valuePoolGrowthRate               = GMTK_MEM_AGGRESSIVE_GROWTH_RATE;

  } else {
    error("%s: Unknown -memoryGrowth option '%s', must be 'conservative', 'default', or 'aggressive'", argerr, memGrowthOption);
  }

#else
#endif
#endif // defined(GMTK_ARG_MEM_GROWTH)


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



/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/

#if defined(GMTK_ARG_USE_MMAP)
#if defined(GMTK_ARGUMENTS_DEFINITION)

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("mmapViterbiValues",Arg::Opt,JunctionTree::mmapViterbi,"Use mmap() to get memory to hold Viterbi values"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#else
#endif
#endif // defined(GMTK_ARG_USE_MMAP)



/*==============================================================================================================*/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************                                     ***********************************************/
/****************************         GENERAL OPTIONS             ***********************************************/
/****************************                                     ***********************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/

#if defined(GMTK_ARG_GENERAL_OPTIONS)
#if defined(GMTK_ARGUMENTS_DOCUMENTATION)
  Arg("\n*** General options ***\n"),
#endif
#endif

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/

#if defined(GMTK_ARG_INFOSEPARATOR)
#if defined(GMTK_ARGUMENTS_DEFINITION)

  static const char *fieldSeparator = "\n";

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("fieldSeparator", Arg::Opt,fieldSeparator,"String that separates fields"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
#else
#endif
#endif // defined(GMTK_ARG_INFOSEPARATOR)


/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/

#if defined(GMTK_ARG_INFOFIELDFILE)
#if defined(GMTK_ARGUMENTS_DEFINITION)

  static char *fieldFile = NULL;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("fieldFile", Arg::Opt,fieldFile,"File listing model info field order"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
#else
#endif
#endif // defined(GMTK_ARG_INFOSEPARATOR)


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

#if defined(GMTK_ARG_SKIP_STARTUP_CHECKS)
#if defined(GMTK_ARGUMENTS_DEFINITION)

  static bool skipStartupChecks = false;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("skipStartupChecks",Arg::Opt,skipStartupChecks, "Skip expensive model validity checks performed at GMTK startup"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
#else
#endif
#endif // defined(GMTK_ARG_SKIP_STARTUP_CHECKS)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_VERB)
#if defined(GMTK_ARGUMENTS_DEFINITION)

#include "debug.h"

#ifdef GMTK_ARG_VERB_DEF_VAL
  static unsigned verbosity = GMTK_ARG_VERB_DEF_VAL;
#else
  static unsigned verbosity = IM::Default;
#endif

  static const char *modularVerbosity = NULL;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("verbosity",Arg::Opt,modularVerbosity,"Verbosity - coma separated list of m=v, where m is all, " moduleHelpString "; 0 <= v <= 100"),
  Arg("printIntValues",Arg::Opt,RV::alwaysPrintIntegerRVValues,"Always print rv values as integer rather than symbols"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

  (void) IM::setGlbMsgLevel(verbosity);
  GM_Parms.setMsgLevel(verbosity);
  for (unsigned i= 0; i < IM::ModuleCount; i+=1) {
    IM::setGlbMsgLevel((IM::ModuleName)i, verbosity);
    GM_Parms.setMsgLevel((IM::ModuleName)i, verbosity);
  }
  
  if (modularVerbosity) {
    char *token, *copy;
    const char delimiters[] = ",";
    copy = strdup(modularVerbosity);
    for (token=strtok(copy, delimiters); token; token=strtok(NULL, delimiters)) {
      (void) IM::setGlbMsgLevel(token);
//      GM_Parms.setMsgLevel(token);
    }
  }

#else
#endif
#endif // defined(GMTK_ARG_VERB)



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




/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/

#if defined(GMTK_ARG_CLIQUE_PRINT_NORMALIZE)
#if defined(GMTK_ARGUMENTS_DEFINITION)

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("cliquePrintNormalize",Arg::Opt,JunctionTree::normalizePrintedCliques,"Normalize scores of each printed clique to probabilities."),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#else
#endif
#endif // defined(GMTK_ARG_CLIQUE_PRINT_NORMALIZE)


/*==============================================================================================================*/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************                                     ***********************************************/
/****************************         INFERENCE OPTIONS           ***********************************************/
/****************************                                     ***********************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/

#if defined(GMTK_ARG_INFERENCE_OPTIONS)
#if defined(GMTK_ARGUMENTS_DOCUMENTATION)
  Arg("\n*** Inference options ***\n"),
#endif
#endif

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_ISLAND)
#if defined(GMTK_ARGUMENTS_DEFINITION)

  static bool island=false;
  static unsigned base=3;
  const static char* baseString = "3";
  static bool rootBase=false; // true iff we should use \sqrt T as the logarithm base, otherwise it's constant
  static float islandRootPower=0.5; // allow arbitrary root of T as log base, default is square root
  static unsigned lst=100;

#define GMTK_SQRT_BASE_STRING "root"

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("island",Arg::Opt,island,"Run island algorithm"),
  Arg("base",Arg::Opt,baseString,"Island algorithm logarithm base (integer or 'root')"),
  Arg("root",Arg::Opt,islandRootPower,"use T^r as the island logarithm base, where T is the number of frames"),
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

    if (strncasecmp(baseString, GMTK_SQRT_BASE_STRING, strlen(GMTK_SQRT_BASE_STRING)+1 ) == 0) {
      rootBase = true;
      if (islandRootPower < 0.0 || 1.0 < islandRootPower) {
	error("%s: -root %f must be between 0 and 1", argerr, islandRootPower);
      }
    } else {
      int tmp = atoi(baseString);
      if (tmp < 2) {
	error("%s: -base %d is too small, it must be >= 2", argerr, tmp);
      }
      base = (unsigned) tmp;
    }
  }


#else
#endif
#endif // defined(GMTK_ARG_ISLAND)

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_DEBUG_PART_RNG)
#if defined(GMTK_ARGUMENTS_DEFINITION)

  const static char *pdbrng_str="all";

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("debugPartitions",Arg::Opt,pdbrng_str,"Partition range to generate debug output"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#else
#endif
#endif // defined(GMTK_ARG_DEBUG_PART_RNG)



/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/


#if defined(GMTK_ARG_DEBUG_INCREMENT)
#if defined(GMTK_ARGUMENTS_DEFINITION)

  extern int debugIncrement;

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("debugIncrement",Arg::Opt,debugIncrement,"Increment to adjust inference verbosity on USR1/2 signals"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#else
#endif
#endif // defined(GMTK_ARG_DEBUG_INCREMENT)


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


#if defined(GMTK_ARG_FAIL_ON_ZERO_CLIQUE)
#if defined(GMTK_ARGUMENTS_DEFINITION)

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("failOnZeroClique",Arg::Opt,MaxClique::failOnZeroClique,"abort GMTK program on zero clique errors"),

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

#else
#endif
#endif // defined(GMTK_ARG_FAIL_ON_ZERO_CLIQUE)



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

#if defined(GMTK_ARG_EM_TRAINING_OPTIONS)
#if defined(GMTK_ARGUMENTS_DOCUMENTATION)
  Arg("\n*** EM training options ***\n"),
#endif
#endif

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

#if defined(GMTK_ARGUMENTS_ONLINE_NORMALIZATION)
static float normalizeScoreEachClique = 0.0;
#else
static float normalizeScoreEachClique = MaxClique::normalizeScoreEachClique;
#endif

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("cliqueTableNormalize",Arg::Opt,normalizeScoreEachClique,"Normalize scores of each clique right after its creation (increases dynamic range)."),


#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)
    MaxClique::normalizeScoreEachClique = normalizeScoreEachClique;

    if (MaxClique::normalizeScoreEachClique < 0.0) {
      error("ERROR: -cliqueTableNormalize option must be non-negative\n");
    }

#if defined(GMTK_ARG_EM_TRAINING_PARAMS)
    if (MaxClique::normalizeScoreEachClique != 1.0 && localCliqueNormalization == false) {
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

#if defined(GMTK_ARG_KERNEL_OPTIONS)
#if defined(GMTK_ARGUMENTS_DOCUMENTATION)
  Arg("\n*** Kernel options ***\n"),
#endif
#endif

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
/************************            DECODING OPTIONS                  ******************************************/
/************************                                              ******************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/

#if defined(GMTK_ARG_DECODING_OPTIONS)
#if defined(GMTK_ARGUMENTS_DOCUMENTATION)
  Arg("\n*** Decoding options ***\n"),
#endif
#endif

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
  static char* vitValsFileName = NULL;
  // TODO: get binary printing working
  // static bool vitValsFileBinp = false;
  static char* vitRegexFilter = NULL;
  static bool vitCaseSensitiveRegexFilter = false;
  static char* vitPartRangeFilter = NULL;
  static char* vitFrameRangeFilter = NULL;
  static bool vitAlsoPrintObservedVariables = false;
//  static bool vitReverseOrder = false;
#define MAX_VITERBI_TRIGGERS 3
  const char   *vitTriggerVariables[MAX_VITERBI_TRIGGERS] = { NULL, NULL, NULL };
  const char   *vitTriggerSets[MAX_VITERBI_TRIGGERS] = { NULL, NULL, NULL };



#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

  Arg("\n*** Decoding options ***\n"),
  // partition based 
  Arg("pVitValsFile",Arg::Opt,pVitValsFileName,"Partition Vit: file to print viterbi values, '-' for stdout"),
  // TODO: not currently used, but should add.
  // Arg("pVitBinVitValsFile",Arg::Opt,pVitValsFileBinp,"Partition Vit: Should file to print viterbi values be binary? (T/F) "),

  Arg("pVitRegexFilter",Arg::Opt,pVitRegexFilter,"Partition Vit: Regular expression to filter variable names."),
  Arg("pVitCaseSensitiveRegexFilter",Arg::Opt,pVitCaseSensitiveRegexFilter,"Partition Vit: Case sensitivity of the rv regular expression filter."),

  Arg("pVitPrintRange",Arg::Opt,pVitPartRangeFilter,"Partition Vit: value printing, integer range filter for modified partitions (e.g., frames, slices) to print."),

  Arg("pVitPrintObservedVariables",Arg::Opt,pVitAlsoPrintObservedVariables,"Partition Vit: also print observed random variables in addtion to hidden"),

  // frame based
  Arg("vitValsFile",Arg::Opt,vitValsFileName,"Vit: file to print viterbi values, '-' for stdout"),
  // TODO: not currently used, but should add.
  // Arg("vitBinVitValsFile",Arg::Opt,vitValsFileBinp,"Vit: Should file to print viterbi values be binary? (T/F) "),

  Arg("vitRegexFilter",Arg::Opt,vitRegexFilter,"Vit: Regular expression to filter variable names."),
  Arg("vitCaseSensitiveRegexFilter",Arg::Opt,vitCaseSensitiveRegexFilter,"Vit: Case sensitivity of the rv regular expression filter."),

  Arg("vitPrintRange",Arg::Opt,vitPartRangeFilter,"Vit: value printing, integer range filter for original partitions (e.g., frames, slices) to print."),
  Arg("vitFrameRange",Arg::Opt,vitFrameRangeFilter,"Vit: value printing, integer range filter for frames to print."),

  Arg("vitPrintObservedVariables",Arg::Opt,vitAlsoPrintObservedVariables,"Vit: also print observed random variables in addtion to hidden"),

#if defined(GMTK_ARGUMENTS_REQUIRE_BINARY_VIT_FILE)
  Arg("binaryVitFile",Arg::Req,JunctionTree::binaryViterbiFilename,"File containing binary Viterbi values for printing"),
#else
  Arg("binaryVitFile",Arg::Opt,JunctionTree::binaryViterbiFilename,"File to write binary Viterbi values for later printing"),
#endif
#if 0
  // this is not implemented yet.
  Arg("vitReverseOrder",Arg::Opt,vitReverseOrder,"Vit: print values in reverse order."),

  Arg("vitTriggerVar", Arg::Opt,vitTriggerVariables,"Viterbi: Trigger variable. Replace X with trigger variable number",Arg::ARRAY,MAX_VITERBI_TRIGGERS),
  Arg("vitTriggerSet", Arg::Opt,vitTriggerSets,"Viterbi: Trigger value set. Replace X with trigger set number",Arg::ARRAY,MAX_VITERBI_TRIGGERS),
#endif

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

  if (JunctionTree::binaryViterbiFilename) {
#if defined(GMTK_VITERBI_FILE_WRITE)
    JunctionTree::binaryViterbiFile = fopen(JunctionTree::binaryViterbiFilename, "w+b");
#else
    JunctionTree::binaryViterbiFile = fopen(JunctionTree::binaryViterbiFilename, "rb");
#endif
    if (!JunctionTree::binaryViterbiFile) {
      char *err = strerror(errno);
      error("ERROR: Failed to open '%s': %s\n", JunctionTree::binaryViterbiFilename, err);
    }
  }

  if (vitPartRangeFilter && vitFrameRangeFilter) {
    error("%s: Can't use both -vitPrintRange and -vitFrameRange\n", argerr);
  }

  if ( (vitPartRangeFilter || vitFrameRangeFilter) && ! vitValsFileName ) {
    error("%s: -vitPrintRange and -vitFrameRange require -vitValsFile to be specified\n", argerr);
  }
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

#if defined(GMTK_ARG_TIMING_OPTIONS)
#if defined(GMTK_ARGUMENTS_DOCUMENTATION)
  Arg("\n*** Timing options ***\n"),
#endif
#endif

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


#if defined(GMTK_ARG_TRIANGULATION_OPTIONS)
#if defined(GMTK_ARGUMENTS_DOCUMENTATION)
  Arg("\n*** Triangulation options ***\n"),
#endif
#endif

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

  static bool writeComments = true;

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
  Arg("writeComments",Arg::Opt,writeComments,"true/false status of if we should write comments in ouptut trifiles"),


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

#if defined(GMTK_ARG_TYING_OPTIONS)
#if defined(GMTK_ARGUMENTS_DOCUMENTATION)
  Arg("\n*** Tying options ***\n"),
#endif
#endif

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

  Arg("\n*** Lattice options ***\n"),
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

#if defined(GMTK_ARG_RESOURCE_OPTIONS)
#if defined(GMTK_ARGUMENTS_DOCUMENTATION)
  Arg("\n*** Resource limiting options ***\n"),
#endif
#endif

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

