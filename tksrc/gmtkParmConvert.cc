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

/*
 * This program converts from ascii trainable parameters to binary
 * and vice versa.
 *
 *
 */

#if HAVE_CONFIG_H
#include <config.h>
#endif
#if HAVE_HG_H
#include "hgstamp.h"
#endif

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

VCID(HGID)


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
#define GMTK_ARG_SKIP_STARTUP_CHECKS
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


/*
 *  To read the file containing the HMM names. Note that this file also contains 
 *  the size of each HMM (i.e. number of states in the HMM) 
 *
 */
void
read_and_parse_hmm_list_file(char *hmmListFileName,vector<string> &hmm_names,
			     vector<int> &numStates, vector<int> &enable, bool endis){
  printf("Reading HMM Names from %s\n",hmmListFileName);
  fflush(stdout);
  iDataStreamFile iif(hmmListFileName,false,false);
  int length;
  iif.read(length,"Can't read hmm names list length");

  hmm_names.resize(length);
  numStates.resize(length);
  enable.resize(length);
  for (int i=0;i<length;i++) {
    iif.read(hmm_names[i],"Can't read hmm names");
    iif.read(numStates[i],"Can't read number of states");
    if (endis) iif.read(enable[i],"Can't read enable/disable information");
    else enable[i] = 1; // if there is no info, then simply enable. 
    //    printf ("%d\n",enable[i]);
  }
}


/*
 *  To read the file containing the HTK Header information and write this out
 *  as a header to the output parameter file. 
 *
 */
void 
read_header_and_writeout(oDataStreamFile& os,char *htkHeaderFile){

  iDataStreamFile iif(htkHeaderFile,false,false);
  vector<string> header;
  header.resize(100); // max number of lines in header = 100.

  string tmp;
  unsigned int counter = 0;
  do{
    iif.readStringUntil(tmp,'\n',false,
			"Read HTK Header: Unable to read string in file");
    header[counter++] = tmp;
  }while ((iif.prepareNext()));
  
  for (unsigned int ii = 0;ii < counter;ii++) {
    os.write(header[ii]);
    os.nl();
  }
}

/*
 *  To read the file containing the HTK Footer information and write this out
 *  as a footer to the output parameter file. 
 *
 */
void 
read_footer_and_writeout(oDataStreamFile& os,char *htkFooterFile){

  iDataStreamFile iif(htkFooterFile,false,false);
  vector<string> header;
  header.resize(100); // max number of lines in header = 100.

  string tmp;
  unsigned int counter = 0;
  do{
    iif.readStringUntil(tmp,'\n',false,
			"Read HTK Footer: Unable to read string in file");
    header[counter++] = tmp;
  }while ((iif.prepareNext()));
  
  for (unsigned int ii = 0;ii < counter;ii++) {
    os.write(header[ii]);
    os.nl();
  }
}

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


  printf("Finished reading in all parameters and structures\n");
  printf("Total number of trainable parameters in input files = %u\n",
	 GM_Parms.totalNumberParameters());

  GM_Parms.markUsedMixtureComponents();

  if (outputMasterFile != NULL) {
    GM_Parms.write(outputMasterFile,cppCommandOptions);
  }

  if (outputTrainableParameters != NULL) {

    oDataStreamFile of(outputTrainableParameters,binOutputTrainableParameters);

    if (htkoutputTrainableParameters) {
      printf("Writing Parameters in HTK format\n");
      fflush(stdout);
      vector<string> hmm_names;
      vector<int> numStates;
      vector<int> enable;
      if (hmmListFileName != NULL)
	read_and_parse_hmm_list_file(hmmListFileName,hmm_names,
				     numStates,enable,clusterHMMs);
      if (htkHeaderFile != NULL)
	read_header_and_writeout(of,htkHeaderFile);
            
      GM_Parms.writeTrainableHTK(of,transitionMatrixName,triphoneCollectionName,
				 hmm_names,numStates,enable,teeModelforsp);
      if (htkFooterFile != NULL)
	read_footer_and_writeout(of,htkFooterFile);

    }
    else  
      GM_Parms.writeTrainable(of);
  }

  exit_program_with_status(0);
}

