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
#include "GMTK_ObservationMatrix.h"
#include "GMTK_MixtureCommon.h"
#include "GMTK_GaussianComponent.h"
#include "GMTK_MeanVector.h"
#include "GMTK_DiagCovarVector.h"
#include "GMTK_DlinkMatrix.h"

#include "GMTK_WordOrganization.h"


#define GMTK_ARG_OBS_FILES
/****************************      FILE RANGE OPTIONS             ***********************************************/
#define GMTK_ARG_DCDRNG
#define GMTK_ARG_START_END_SKIP
/************************  OBSERVATION MATRIX TRANSFORMATION OPTIONS   ******************************************/
#define GMTK_ARG_OBS_MATRIX_XFORMATION


#define GMTK_ARG_INPUT_MASTER_FILE_OPT_ARG
#define GMTK_ARG_OUTPUT_MASTER_FILE
#define GMTK_ARG_OUTPUT_TRAINABLE_PARAMS
#define GMTK_ARG_INPUT_TRAINABLE_PARAMS
#define GMTK_ARG_CPP_CMD_OPTS
#define GMTK_ARG_VERB
#define GMTK_ARG_ALLOC_DENSE_CPTS
#define GMTK_ARG_SEED
#define GMTK_ARG_VAR_FLOOR
#define GMTK_ARG_STR_FILE_OPT_ARG
#define GMTK_ARG_VAR_FLOOR
#define GMTK_ARG_VERSION
#define GMTK_ARG_VAR_FLOOR_ON_READ
#define GMTK_ARG_CPT_NORM_THRES

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
ObservationMatrix globalObservationMatrix;

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
  bool parse_was_ok = Arg::parse(argc,(char**)argv);
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
				    gpr_str);
  infoMsg(IM::Max,"Finished opening files.\n");


  ////////////////////////////////////////////
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
  bool decTree = GM_Parms.dts.size() > 0;
  bool itDT = GM_Parms.iterableDts.size() > 0;
  bool internalDT = false;
  bool Fngram = GM_Parms.fngramCpts.size() > 0;
  bool ngram  = GM_Parms.ngramCpts.size() > 0;
  bool lattice = GM_Parms.latticeAdts.size() > 0 ||
                 GM_Parms.iterableLatticeAdts.size() > 0 ||
                 GM_Parms.latticeNodeCpts.size() > 0 ||
                 GM_Parms.latticeEdgeCpts.size() > 0;
  bool sparseGauss = false; // GM_Parms.sPmfs.size() > 0;
  bool l1Reg = false;
  bool l2Reg = false;
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

#define BF "%u"
#define PB(b) ((b) ? 1 : 0)

  char swWght[4] = {0,0,0,0};
  int idx = 0;
  if (scale)   swWght[idx++] = 's';
  if (penalty) swWght[idx++] = 'p';
  if (shift)   swWght[idx++] = 'h';
//#rvs  #rvs in P C E  #Drvs  #Crvs  min card  max card  #hid  #obs  min Crv dim  max Crv dim  gr comp.  #frames P C E  min fr  max fr
//sw par  sw wght  sym  dense  sparse  det  DT  ItDT  IDT  Fng  Ng  Lat  sparseG  l1  l2  exp  gamma  beta  VECPT  VESEP format
  printf("%u %u %u %u %u %u %u %u %u %u %u %u XXX %u %u %u XXX XXX "BF" '%s' "BF" "BF" "BF" "BF" "BF" "BF" "BF" "BF" "BF" "BF" X X X X X X "BF" X %s\n", 
	 nRVs, nRVsInP, nRVsInC, nRVsInE, nDiscRVs, nContRVs, minRVcard, maxRVcard, nHidRVs, nObsRVs, minRVdim, maxRVdim, NP, NC, NE, 
	 PB(swPar), swWght, PB(symTab), PB(denseCPT), PB(sparseCPT), PB(determCPT), PB(diagGauss), PB(decTree), PB(itDT), PB(Fngram), PB(ngram), PB(lattice), /*PB(sparseGauss),*/ PB(veCPT), fmts[0]);
  assert(nRVs = nRVsInP + nRVsInC + nRVsInE);
  assert(nRVs = nObsRVs + nHidRVs);
  assert(nRVs = nDiscRVs + nContRVs);

  if (outputMasterFile != NULL) {
    GM_Parms.write(outputMasterFile,cppCommandOptions);
  }
  if (outputTrainableParameters != NULL) {
    oDataStreamFile of(outputTrainableParameters,binOutputTrainableParameters);
    GM_Parms.writeTrainable(of);
  }

  exit_program_with_status(0);
}

