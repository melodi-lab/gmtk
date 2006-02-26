/*
 * gmtkTie.cc
 * perform model parameter tying using one of a variety of methods
 *
 * Written by Simon King <Simon.King@ed.ac.uk>
 *
 * Copyright (c) 2006, < fill in later >
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
 * This program performs model parameter tying
 *
 * 1 Currently, these parameter types can be tied:
 *  1.1 Gaussian mixture distributions
 *  1.2 MeanVectors
 * 2 Supported tying methods are:

 *  - Basic data-driven clustering using these distance metrics:

 *   + for single (diagonal only?) Gaussians, a weighted (by the
 *   variance) Euclidean distance between the means (equation 3.1 from
 *   Odell's thesis)

 *   + for GMMs, the average probability of each component mean in
 *   GMM1, with respect to GMM2, plus the reverse

 *
 *
 * Future things to implement:
 * - other distribution types
 * - decision-tree based tying using these clustering criteria
 *  + likelihood based criterion: equations 3.2-3.14 from Odell's thesis
 * - outlier merging
 * - tied-mixture systems, with tying based on mixture weight values only
 *
 *
 *
 *
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
#include <list.h>

#include "general.h"
#include "error.h"
#include "rand.h"
#include "arguments.h"
#include "ieeeFPsetup.h"
#include "version.h"

#include "GMTK_FileParser.h"
#include "GMTK_RV.h"
#include "GMTK_DiscRV.h"
#include "GMTK_ContRV.h"
#include "GMTK_GMParms.h"
#include "GMTK_ObservationMatrix.h"
#include "GMTK_MixtureCommon.h"
#include "GMTK_GaussianComponent.h"
#include "GMTK_MeanVector.h"
#include "GMTK_DiagCovarVector.h"
#include "GMTK_DlinkMatrix.h"
#include "GMTK_Mixture.h"
#include "GMTK_Tie.h"

#include "tieSupport.h"

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


////////////////////////////////////////////
// command line arguments specific to gmtkTie
#define GMTK_ARG_TYING_PARAMS

#define GMTK_ARGUMENTS_DEFINITION
#include "GMTK_Arguments.h"
#undef GMTK_ARGUMENTS_DEFINITION


Arg Arg::Args[] = {


#define GMTK_ARGUMENTS_DOCUMENTATION
#include "GMTK_Arguments.h"
#undef GMTK_ARGUMENTS_DOCUMENTATION


#define NAMED_COLLECTION_GLOBAL_NAME "global"

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
  // we don't currently need the structure file, but will in the
  // future
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

  /////////////////////////////////////////////
  // only run if there is somewhere to save the result
  if (outputTrainableParameters == NULL)
    error("No output traininable parameters file specified\n");

  ////////////////////////////////////////////
  if (inputMasterFile != NULL) {
    iDataStreamFile pf(inputMasterFile,false,true,cppCommandOptions);
    GM_Parms.read(pf);
    infoMsg(IM::Tiny,"Finished reading master file\n");
  }

  if (inputTrainableParameters != NULL) {
    // flat, where everything is contained in one file
    iDataStreamFile pf(inputTrainableParameters,binInputTrainableParameters,true,cppCommandOptions);
    GM_Parms.readTrainable(pf);
    infoMsg(IM::Tiny,"Finished reading trainable params file\n");
  }
  GM_Parms.finalizeParameters();  

  if (strFileName != NULL) {
    // load up the structure file as we might want
    // it to allocate some Dense CPTs.
    FileParser fp(strFileName,cppCommandOptions);

    // parse the file
    fp.parseGraphicalModel();
    // create the rv variable objects
    fp.createRandomVariableGraph();
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
  }
  infoMsg(IM::Tiny,"Finished reading structure file\n");


  ////////////////////////////////////////////
  // the object responsible for all tying processes
  GMTK_Tie tie(&GM_Parms);


  ////////////////////////////////////////////
  // the file containing the user-supplied commands and their optional
  // parameters
  static char *tieCommandsFile=loadCmdFile;

  ////////////////////////////////////////////
  // load in the tying command file
  if (tieCommandsFile != NULL) {
    iDataStreamFile cf(tieCommandsFile,false,true,cppCommandOptions);

    ////////////////////////////////////////////
    // load and parse the command file
    if (!tie.read_commands(cf)) 
      error("Errors in command file");
    infoMsg(IM::Tiny,"Finished reading tie command file\n");

    ////////////////////////////////////////////
    // do some preliminary error checking
    if (!tie.validate_command_parameters())
      error("Errors in parameters named in command file");
    infoMsg(IM::Mod,"Finished validating tie command parameters\n");

  } 
  else
    ////////////////////////////////////////////
    // it makes no sense to run this program with no commands
    error("Error: must provide a tie command file");


  ////////////////////////////////////////////
  // work out if we need to load accumulator file(s)
  bool needAccFile = tie.need_occupancy_counts();

  ////////////////////////////////////////////
  // if we do, then load them
  logpr total_data_prob = 1.0;
  if (needAccFile){
    if(loadAccFile != NULL) {
      if (loadAccRange == NULL) {
	infoMsg(IM::Default,"Loading accumulators from '%s'\n",loadAccFile);
	iDataStreamFile inf(loadAccFile,accFileIsBinary);
	inf.read(total_data_prob.valref());
	GM_Parms.emLoadAccumulators(inf);
      } else {
	Range lfrng(loadAccRange,0,1000);
	for (Range::iterator lfit=lfrng.begin();
	     !lfit.at_end();
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
      infoMsg(IM::Tiny,"Finished loading all accumulator files\n");
    } else
      error("Accumulator file(s) needed but not specified on command line");
  } else if (loadAccFile != NULL)
      warning("Accumulator file(s) specified on command line are not required - not loading them");

  infoMsg(IM::Tiny,"Finished reading in all files\n");

  ////////////////////////////////////////////
  // remove any unused parameters from GM_Parms
  //
  // should make this optional - there may be some situation where
  // unused parameters should be retained to take part in the tying
  // procedure
  infoMsg(IM::Tiny,"Note: unused parameters will NOT be saved\n");
  GM_Parms.markUsedMixtureComponents();
  infoMsg(IM::Mod,"Finished marking used params\n");

  ////////////////////////////////////////////
  // go through the commands and expand out the list of matching
  // parameters from GM_Parms
  tie.find_matching_command_parameters();
  
  ////////////////////////////////////////////
  // execute the commands in the order given in the tying command file
  infoMsg(IM::Tiny,"Executing a list of %d commands\n",tie.commands.size());
  for (unsigned i=0;i<tie.commands.size();i++)
    if (!tie.execute_command(i))
      ////////////////////////////////////////////
      // do not reduce this to a warning because subsequent commands
      // could reply on the success of earlier ones
      error("Command %d failed",i);
	
  /*
    no changes are currently made to the master file
  if (outputMasterFile != NULL) {
    GM_Parms.write(outputMasterFile,cppCommandOptions);
  }
  */

  oDataStreamFile of(outputTrainableParameters,binOutputTrainableParameters);
  GM_Parms.writeTrainable(of);

  exit_program_with_status(0);
}






