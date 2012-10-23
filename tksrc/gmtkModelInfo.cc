/*
 * gmtkModelInfo.cc
 * produce a junction tree
 *
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (c) 2011, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about the suitability of this software
 * for any purpose. It is provided "as is" without express or implied warranty.
 *
 */

/*
 * This program converts from ascii trainable parameters to binary
 * and vice versa.
 *
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
#include "rand.h"
#include "arguments.h"
#include "ieeeFPsetup.h"
#include "debug.h"
//#include "spi.h"
#include "version.h"

#if HAVE_CONFIG_H
#include <config.h>
#endif
#if HAVE_HG_H
#include "hgstamp.h"
#endif
VCID(HGID)


#include "GMTK_FileParser.h"
#include "GMTK_RV.h"
#include "GMTK_DiscRV.h"
#include "GMTK_ContRV.h"
#include "GMTK_GMParms.h"
#include "GMTK_GMTemplate.h"
#include "GMTK_Partition.h"
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
#include "GMTK_LinMeanCondDiagGaussian.h"
#include "GMTK_MeanVector.h"
#include "GMTK_DiagCovarVector.h"
#include "GMTK_DlinkMatrix.h"

#include "GMTK_WordOrganization.h"


#define GMTK_ARG_OBS_FILES
/****************************      FILE RANGE OPTIONS             ***********************************************/
#define GMTK_ARG_FILE_RANGE_OPTIONS
#define GMTK_ARG_DCDRNG
#define GMTK_ARG_START_END_SKIP
/************************  OBSERVATION MATRIX TRANSFORMATION OPTIONS   ******************************************/
#define GMTK_ARG_OBS_MATRIX_OPTIONS
#define GMTK_ARG_OBS_MATRIX_XFORMATION


/*************************   INPUT TRAINABLE PARAMETER FILE HANDLING  *******************************************/
#define GMTK_ARG_INPUT_TRAINABLE_FILE_HANDLING
#define GMTK_ARG_INPUT_MASTER_FILE_OPT_ARG
#define GMTK_ARG_DLOPEN_MAPPERS
#define GMTK_ARG_OUTPUT_MASTER_FILE
#define GMTK_ARG_OUTPUT_TRAINABLE_PARAMS
#define GMTK_ARG_INPUT_TRAINABLE_PARAMS
#define GMTK_ARG_CPP_CMD_OPTS
#define GMTK_ARG_ALLOC_DENSE_CPTS
#define GMTK_ARG_CPT_NORM_THRES

/*************************   INPUT STRUCTURE PARAMETER FILE HANDLING  *******************************************/
#define GMTK_ARG_INPUT_MODEL_FILE_HANDLING
#define GMTK_ARG_STR_FILE_OPT_ARG

/*************************   CONTINUOUS RANDOM VARIABLE OPTIONS       *******************************************/
#define GMTK_ARG_CONTINUOUS_RANDOM_VAR_OPTIONS
#define GMTK_ARG_VAR_FLOOR
#define GMTK_ARG_VAR_FLOOR_ON_READ


#define GMTK_ARG_GENERAL_OPTIONS
#define GMTK_ARG_SEED
#define GMTK_ARG_VERB
#define GMTK_ARG_VERSION
#define GMTK_ARG_HELP

#define GMTK_ARG_INFOSEPARATOR
#define GMTK_ARG_INFOFIELDFILE

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
 * definition of needed global arguments
 */
RAND rnd(false);
GMParms GM_Parms;
#if 0
ObservationMatrix globalObservationMatrix;
#endif

FileSource *gomFS;
ObservationSource *globalObservationMatrix;

  // these are the available fields
  enum fields {
    nRVsInP_F,
    nRVsInC_F,
    nRVsInE_F,
    nDiscRVs_F,
    nContRVs_F,
    nObsRVs_F,
    nHidRVs_F,
    maxRVcard_F,
    minRVcard_F,
    maxRVdim_F,
    minRVdim_F,
    graphComp_F,
    nRVs_F,
    NP_F,
    NC_F,
    NE_F,
    minFrames_F,
    maxFrames_F,
    swPar_F,
    swWght_F,
    symTab_F,
    denseCPT_F,
    sparseCPT_F,
    determCPT_F,
    diagGauss_F,
    decTree_F,
    itDT_F,
    internalDT_F,
    Fngram_F,
    ngram_F,
    lattice_F,
    sparseGauss_F,
    l1Reg_F,
    l2Reg_F,
    expDist_F,
    gammaDist_F,
    betaDist_F,
    veCPT_F,
    veSep_F,
    format_F,
    unknown_F, // use this for fields you don't know how to compute
    done}; // done must always be last!

fields
lookup(const char*list[], char *target) {
  unsigned i;
  for (i=0; strcmp(target, list[i]) && strcmp("done", list[i]); i+=1)
    ;
  if (strcmp(target, list[i]) == 0)
    return (fields)i;
  else
    return done;
}


int
main(int argc,char*argv[])
{
  ////////////////////////////////////////////
  // set things up so that if an FP exception
  // occurs such as an "invalid" (NaN), overflow
  // or divide by zero, we actually get a FPE
  ieeeFPsetup();

  CODE_TO_COMPUTE_ENDIAN;

  ////////////////////////////////////////////
  // parse arguments
  bool parse_was_ok = Arg::parse(argc,(char**)argv,
"\nThis program prints out some information about the number of variables\n"
"in a model and which GMTK features the model uses\n");
  if(!parse_was_ok) {
    Arg::usage(); exit(-1);
  }


#define GMTK_ARGUMENTS_CHECK_ARGS
#include "GMTK_Arguments.h"
#undef GMTK_ARGUMENTS_CHECK_ARGS


  /////////////////////////////////////////////

  if (strFileName == NULL || allocateDenseCpts == 0) {
    // check this only if the structure file is not given or
    // if we are not allocating Dense CPTs, since
    // in that case there is no way the user could specify
    // automatic allocation of such CPTs.
    if ((inputMasterFile == NULL) && (inputTrainableParameters == NULL)) {
      warning("ERROR: need to specify command line parameters inputMasterFile or inputTrainableParameters (or both) when no structure file is given");
      Arg::usage();
      error("");
    }
  }

  infoMsg(IM::Max,"Opening Files ...\n");
  gomFS = instantiateFileSource();
  globalObservationMatrix = gomFS;
  infoMsg(IM::Max,"Finished opening files.\n");


  ////////////////////////////////////////////
  dlopenDeterministicMaps(dlopenFilenames, MAX_NUM_DLOPENED_FILES);
  if (inputMasterFile != NULL) {
    iDataStreamFile pf(inputMasterFile,false,true,cppCommandOptions);
    GM_Parms.read(pf);
  }
  if (inputTrainableParameters != NULL) {
    // flat, where everything is contained in one file
    iDataStreamFile pf(inputTrainableParameters,binInputTrainableParameters,true,cppCommandOptions);
    GM_Parms.readTrainable(pf);
  }
  GM_Parms.finalizeParameters();  


  // load up the structure file as we might want
  // it to allocate some Dense CPTs.
  FileParser fp(strFileName,cppCommandOptions);
  infoMsg(IM::Tiny,"Finished reading in all parameters and structures\n");
    
  // parse the file
  infoMsg(IM::Max,"Parsing structure file...\n");
  fp.parseGraphicalModel();
  // create the rv variable objects
  infoMsg(IM::Max,"Creating rv objects...\n");
  fp.createRandomVariableGraph();

  // Make sure that there are no directed loops in the graph.
  infoMsg(IM::Max,"Checking template...\n");
  fp.ensureValidTemplate();

  // link the RVs with the parameters that are contained in
  // the bn1_gm.dt file.
  if (allocateDenseCpts == 0)
    fp.associateWithDataParams(FileParser::noAllocate);
  else if (allocateDenseCpts == 1)
    fp.associateWithDataParams(FileParser::allocateRandom);
  else if (allocateDenseCpts == 2)
    fp.associateWithDataParams(FileParser::allocateUniform);
  else
    error("Error: command line argument '-allocateDenseCpts d', must have d = {0,1,2}\n");
  
  
  printf("Finished reading in all parameters and structures\n");
  printf("Total number of trainable parameters in input files = %u\n",
	 GM_Parms.totalNumberParameters());
  GM_Parms.markUsedMixtureComponents();
  
  
  infoMsg(IM::Max,"Creating template...\n");
  //  GMTemplate gm_template(fp);
  
  vector<RV*> unrolled_rvs;
  map<RVInfo::rvParent, unsigned> unrolled_map;
  fp.unroll(0, unrolled_rvs, unrolled_map);

  
  unsigned nRVsInP = 0;
  unsigned nRVsInC = 0;
  unsigned nRVsInE = 0;
  unsigned nDiscRVs = 0;
  unsigned nContRVs = 0;
  unsigned nObsRVs = 0;
  unsigned nHidRVs = 0;
  unsigned maxRVcard = 0;
  unsigned minRVcard = (unsigned)(-1);
  unsigned maxRVdim = 0;
  unsigned minRVdim = (unsigned)(-1);
  unsigned nRVs = unrolled_rvs.size();
  unsigned NP = fp.numFramesInP();
  unsigned NC = fp.numFramesInC();
  unsigned NE = fp.numFramesInE();
  bool swPar = false;
  bool scale = false;
  bool penalty = false;
  bool shift = false;
  bool symTab = false;
  bool denseCPT  = GM_Parms.mdCpts.size() > 0;
  bool sparseCPT = GM_Parms.msCpts.size() > 0;
  bool determCPT = GM_Parms.mtCpts.size() > 0;
  bool diagGauss = GM_Parms.covars.size() > 0;
  unsigned externalDTcount = 0;
  for (GMParms::ObjectMapType::iterator it = GM_Parms.dtsMap.begin(); 
       it != GM_Parms.dtsMap.end(); 
       it++)
    {
      if (it->first.find("internal:") == string::npos) {
	externalDTcount += 1;
      }
      //fprintf(stderr, "dt %3u:   %u %s\n", it->second, it->first.find("internal:"), it->first.c_str());
    }
#if 0
  bool decTree = GM_Parms.dts.size() > 0;
#else
  bool decTree = externalDTcount > 0;
#endif
  //fprintf(stderr, "%u dts, %u mapsize, ext %u  DT "BF"\n", GM_Parms.dts.size(), GM_Parms.dtsMap.size(), externalDTcount, PB(decTree));
  bool itDT = GM_Parms.iterableDts.size() > 0;
  bool internalDT = false;
  bool Fngram = GM_Parms.fngramCpts.size() > 0;
  bool ngram  = GM_Parms.ngramCpts.size() > 0;
  bool lattice = GM_Parms.latticeAdts.size() > 0 ||
                 GM_Parms.iterableLatticeAdts.size() > 0 ||
                 GM_Parms.latticeNodeCpts.size() > 0 ||
                 GM_Parms.latticeEdgeCpts.size() > 0;
  bool sparseGauss = false; // GM_Parms.sPmfs.size() > 0;
  bool l1Reg = false; // L1 regularization not currently implemented
  bool l2Reg = false;
  for (vector< Component* >::iterator it = GM_Parms.components.begin();
       it != GM_Parms.components.end(); 
       it++)
  {
    LinMeanCondDiagGaussian *comp = dynamic_cast<LinMeanCondDiagGaussian *>(*it);
    if (comp != NULL) {
      l2Reg = l2Reg || comp->regularized();
    }
  }
  bool expDist = false;
  bool gammaDist = false;
  bool betaDist = false;
  bool veCPT = GM_Parms.veCpts.size() > 0;
  bool veSep = false;

  for (vector<RV*>::iterator it = unrolled_rvs.begin(); it != unrolled_rvs.end(); ++it) {
    if ((*it)->frame() < fp.firstChunkFrame()) {
      nRVsInP += 1;
    } else if ((*it)->frame() <= fp.lastChunkFrame()) {
      nRVsInC += 1;
    } else {
      nRVsInE += 1;
    }
    if ((*it)->discrete()) {
      nDiscRVs += 1;
      unsigned card = ((DiscRV *)(*it))->cardinality;
      if (card > maxRVcard) maxRVcard = card;
      if (card < minRVcard) minRVcard = card;
      DiscRV *drv = (DiscRV *) (*it);
      symTab = symTab || drv->symbolTable();
    } else if ((*it)->continuous()) {
      nContRVs += 1;
      unsigned dim = (unsigned) ((ContRV*)(*it))->dimensionality();
      if (dim > maxRVdim) maxRVdim = dim;
      if (dim < minRVdim) minRVdim = dim;
    }
    if ((*it)->hidden()) {
      nHidRVs += 1;
    } else if ((*it)->observed()) {
      nObsRVs += 1;
    }
    if ((*it)->switching()) {
      swPar = true;
    }
    if ((*it)->scale()) {
      scale = true;
    }
    if ((*it)->penalty()) {
      penalty = true;
    }
    if ((*it)->shift()) {
      shift = true;
    }
      
  }
  if (nContRVs == 0) {
    minRVdim = 0;
  }

  char swWght[4] = {0,0,0,0};
  int idx = 0;
  if (scale)   swWght[idx++] = 's';
  if (penalty) swWght[idx++] = 'p';
  if (shift)   swWght[idx++] = 'h';
//#rvs  #rvs in P C E  #Drvs  #Crvs  min card  max card  #hid  #obs  min Crv dim  max Crv dim  gr comp.  #frames P C E  min fr  max fr
//sw par  sw wght  sym  dense  sparse  det  diagGauss DT  ItDT  IDT  Fng  Ng  Lat  sparseG  l1  l2  exp  gamma  beta  VECPT  VESEP format

#if 0
#define BF "%u"
#define PB(b) ((b) ? 1 : 0)
  printf("%u %u %u %u %u %u %u %u %u %u %u %u XXX %u %u %u XXX XXX "BF" '%s' "BF" "BF" "BF" "BF" "BF" "BF" "BF" X "BF" "BF" "BF" X "BF" "BF" X X X "BF" X %s\n", 
	 nRVs, nRVsInP, nRVsInC, nRVsInE, nDiscRVs, nContRVs, minRVcard, maxRVcard, nHidRVs, nObsRVs, minRVdim, maxRVdim, NP, NC, NE, 
	 PB(swPar), swWght, PB(symTab), PB(denseCPT), PB(sparseCPT), PB(determCPT), PB(diagGauss), PB(decTree), PB(itDT), PB(Fngram), PB(ngram), PB(lattice), /*PB(sparseGauss),*/ l1Reg, l2Reg, PB(veCPT), fmts[0]);
#else

  // short name for each field
  const char* fieldNames[done+1] = {
    "nRVsInP",
    "nRVsInC",
    "nRVsInE",
    "nDiscRVs",
    "nContRVs",
    "nObsRVs",
    "nHidRVs",
    "maxRVcard",
    "minRVcard",
    "maxRVdim",
    "minRVdim",
    "graphComp",
    "nRVs",
    "NP",
    "NC",
    "NE",
    "minFrames",
    "maxFrames",
    "swPar",
    "swWght",
    "symTab",
    "denseCPT",
    "sparseCPT",
    "determCPT",
    "diagGauss",
    "decTree",
    "itDT",
    "internalDT",
    "Fngram",
    "ngram",
    "lattice",
    "sparseGauss",
    "l1Reg",
    "l2Reg",
    "expDist",
    "gammaDist",
    "betaDist",
    "veCPT",
    "veSep",
    "format",
    "unknown",
    "done"};

  // verbose name for each field
  const char* longNames[done+1] = {
    "# RVs in prolog",
    "# RVs in chunk",
    "# RVs in epilog",
    "# discrete RVs",
    "# continuous RVs",
    "# observed RVs",
    "# hidden RVs",
    "max RV cardinality",
    "min RV cardinality",
    "max RV dimension",
    "min RV dimension",
    "graph complexity",
    "# total RVs",
    "# frames in proglog",
    "# frames in chunk",
    "# frames in epilog",
    "min frames in data",
    "max frames in data",
    "switching parents",
    "switching weights",
    "symbol tables",
    "dense CPTs",
    "sparse CPTs",
    "deterministic CPTs",
    "diagonal Gaussians",
    "decision trees",
    "iterated decision trees",
    "internal decision trees",
    "Fngram CPTs",
    "Ngram CPTs",
    "lattice CPTs",
    "sparse Gaussians",
    "L1 regularization",
    "L2 regularization",
    "exponential distribution",
    "gamma distribution",
    "beta distribution",
    "VE CPTs",
    "VE separators",
    "data format",
    "unknown field",
    "done"};

  // output the fields in this order
  fields outputOrder[done+1] = {
#if 1
    nRVs_F, 
    nRVsInP_F, 
    nRVsInC_F, 
    nRVsInE_F, 
    nDiscRVs_F, 
    nContRVs_F, 
    minRVcard_F, 
    maxRVcard_F, 
    nHidRVs_F, 
    nObsRVs_F, 
    minRVdim_F, 
    maxRVdim_F, 
    NP_F,
    NC_F, 
    NE_F, 
    swPar_F, 
    swWght_F,
    symTab_F,
    denseCPT_F,
    sparseCPT_F,
    determCPT_F, 
    diagGauss_F, 
    decTree_F, 
    itDT_F, 
    Fngram_F, 
    ngram_F, 
    lattice_F, 
    l1Reg_F, 
    l2Reg_F, 
    veCPT_F, 
    format_F,
#else
    nRVs_F, 
    nRVsInP_F, 
    nRVsInC_F, 
    nRVsInE_F, 
    nDiscRVs_F, 
    nContRVs_F, 
    minRVcard_F, 
    maxRVcard_F, 
    nHidRVs_F, 
    nObsRVs_F, 
    minRVdim_F, 
    maxRVdim_F, 
    unknown_F, // graphComp_F,
    NP_F,
    NC_F, 
    NE_F, 
    unknown_F, // minFrames_F,
    unknown_F, // maxFrames_F,
    swPar_F, 
    swWght_F,
    symTab_F,
    denseCPT_F,
    sparseCPT_F,
    determCPT_F, 
    diagGauss_F, 
    decTree_F, 
    itDT_F, 
    unknown_F, // internalDT_F,
    Fngram_F, 
    ngram_F, 
    lattice_F, 
    unknown_F, // sparseGauss_F,
    l1Reg_F, 
    l2Reg_F, 
    unknown_F, // exp dist
    unknown_F, // gamma dist
    unknown_F, // beta dist
    veCPT_F, 
    unknown_F, // ve sep
    format_F,
#endif
    done};

  //  const char *fieldsSeparator = " ";
  enum outputFormat {terse, normal, verbose};
  outputFormat mode = normal;

#define MAX_FIELD_LENGTH 32
  if (fieldFile) {
    FILE *f = fopen(fieldFile, "r");
    if (!f) error("Failed to open file: %s", fieldFile);
    char field[MAX_FIELD_LENGTH];
    unsigned i = 0;
    for (; i < done-1 && fgets(field, MAX_FIELD_LENGTH, f) != NULL; ) {
      field[strlen(field)-1] = 0;
      outputOrder[i++] = lookup(fieldNames, field);
    }
    outputOrder[i] = done;
  }
  if (IM::glbMsgLevel(IM::ModelInfo) < IM::Default) 
    mode = terse; 
  else if (IM::glbMsgLevel(IM::ModelInfo) >= IM::Moderate)
    mode = verbose;

  FILE *outF = stdout;

#define BF "%u"
#define PB(b) \
switch (mode) {\
 case terse: fprintf(outF, "%s"BF, sep, ((b) ? 1 : 0)); break; \
 case normal: fprintf(outF, "%s%s = "BF, sep, names[outputOrder[i]], ((b) ? 1 : 0)); break; \
 case verbose: fprintf(outF, "%s%s = "BF, sep, names[outputOrder[i]], ((b) ? 1 : 0)); \
}

#define PU(x) \
switch (mode) {\
 case terse: fprintf(outF, "%s%u", sep, (x)); break; \
 case normal: fprintf(outF, "%s%s = %u", sep, names[outputOrder[i]], (x)); break; \
 case verbose: fprintf(outF, "%s%s = %u", sep, names[outputOrder[i]], (x)); \
}

#define PS(x) \
switch (mode) {\
 case terse: fprintf(outF, "%s%s", sep, (x)); break; \
 case normal: fprintf(outF, "%s%s = %s", sep, names[outputOrder[i]], (x)); break; \
 case verbose: fprintf(outF, "%s%s = %s", sep, names[outputOrder[i]], (x)); \
}

  bool printSep = false;
  const char *sep;
  const char **names = NULL;
  if (mode == normal) {
    names = fieldNames;
    //    fieldsSeparator = "\n";
  } else if (mode == verbose) {
    names = longNames;
    //    fieldsSeparator = "\n";
  }

  for (unsigned i=0; outputOrder[i] != done; i+=1) {
    //fprintf (stderr, "%u %u %u\n", i, outputOrder[i], done);
    if (printSep)
      sep = fieldSeparator;
    else
      sep = "";
    switch (outputOrder[i]) {
    case nRVsInP_F:     PU(nRVsInP); break;
    case nRVsInC_F:     PU(nRVsInC); break;
    case nRVsInE_F:     PU(nRVsInE); break;
    case nDiscRVs_F:    PU(nDiscRVs); break;
    case nContRVs_F:    PU(nContRVs); break;
    case nObsRVs_F:     PU(nObsRVs); break;
    case nHidRVs_F:     PU(nHidRVs); break;
    case maxRVcard_F:   PU(maxRVcard); break;
    case minRVcard_F:   PU(minRVcard); break;
    case maxRVdim_F:    PU(maxRVdim); break;
    case minRVdim_F:    PU(minRVdim); break;
    case nRVs_F:        PU(nRVs); break;
    case NP_F:          PU(NP); break;
    case NC_F:          PU(NC); break;
    case NE_F:          PU(NE); break;
    case swPar_F:       PB(swPar); break;
    case swWght_F:      if (mode == terse) {
	                  fprintf(outF,"%s'%s'", sep, swWght);
                        } else {
			  fprintf(outF,"%s%s = '%s'", sep, names[outputOrder[i]], swWght); 
                        }
#if 0
	  fprintf(outF,"%s%s'", sep, (mode==terse ? "" : names[outputOrder[i]]" = ")); 
                        if (scale) fprintf(outF, "s");
                        if (penalty) fprintf(outF, "p");
                        if (shift) fprintf(outF, "h");
                        fprintf(outF, "'"); break;
#endif
			break;
    case symTab_F:      PB(symTab); break;
    case denseCPT_F:    PB(denseCPT); break;
    case sparseCPT_F:   PB(sparseCPT); break;
    case determCPT_F:   PB(determCPT); break;
    case diagGauss_F:   PB(diagGauss); break;
    case decTree_F:     PB(decTree); break;
    case itDT_F:        PB(itDT); break;
    case internalDT_F:  PB(internalDT); break;
    case Fngram_F:      PB(Fngram); break;
    case ngram_F:       PB(ngram); break;
    case lattice_F:     PB(lattice); break;
    case sparseGauss_F: PB(sparseGauss); break;
    case l1Reg_F:       PB(l1Reg); break;
    case l2Reg_F:       PB(l2Reg); break;
    case expDist_F:     PB(expDist); break;
    case gammaDist_F:   PB(gammaDist); break;
    case betaDist_F:    PB(betaDist); break;
    case veCPT_F:       PB(veCPT); break;
    case veSep_F:       PB(veSep); break;
    case format_F:      PS(fmts[0]); break;
    case unknown_F:     PS("?"); break;
    case done : break;
    default: error("forget to update the switch statement after adding a field? %u", outputOrder[i]);
    }
    printSep = true;
  }
  fprintf(outF, "\n");
#endif
  assert(nRVs = nRVsInP + nRVsInC + nRVsInE);
  assert(nRVs = nObsRVs + nHidRVs);
  assert(nRVs = nDiscRVs + nContRVs);

  exit_program_with_status(0);
}

