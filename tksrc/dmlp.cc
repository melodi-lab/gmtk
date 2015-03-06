/*
 * dmlp.cc
 *
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (C) 2014 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
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
#include <signal.h>


#include "general.h"
#include "error.h"
#include "rand.h"
#include "arguments.h"
#include "ieeeFPsetup.h"
#include "debug.h"
#if 0
#  include "version.h"
#endif

#include "GMTK_FileParser.h"
#include "GMTK_RV.h"
#include "GMTK_ObsDiscRV.h"
#include "GMTK_ContRV.h"
#include "GMTK_GMParms.h"
#include "GMTK_GMTemplate.h"
#include "GMTK_Partition.h"

#include "GMTK_ObservationSource.h"
#include "GMTK_FileSource.h"
#include "GMTK_CreateFileSource.h"
#include "GMTK_ASCIIFile.h"
#include "GMTK_FlatASCIIFile.h"
#include "GMTK_PFileFile.h"
#include "GMTK_HTKFile.h"
#include "GMTK_HDF5File.h"
#include "GMTK_BinaryFile.h"
#include "GMTK_Filter.h"
#include "GMTK_Stream.h"

#include "GMTK_MixtureCommon.h"
#include "GMTK_GaussianComponent.h"
#include "GMTK_LinMeanCondDiagGaussian.h"
#include "GMTK_MeanVector.h"
#include "GMTK_DiagCovarVector.h"
#include "GMTK_DlinkMatrix.h"

#include "GMTK_WordOrganization.h"

#include "GMTK_DeepVECPT.h"

VCID(HGID)
static char const *dVeCptName           = NULL;

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

#define GMTK_ARG_OUTPUT_MASTER_FILE
#define GMTK_ARG_OUTPUT_TRAINABLE_PARAMS
#define GMTK_ARG_INPUT_TRAINABLE_PARAMS
#define GMTK_ARG_CPP_CMD_OPTS
#define GMTK_ARG_ALLOC_DENSE_CPTS
#define GMTK_ARG_CPT_NORM_THRES


/*************************   CONTINUOUS RANDOM VARIABLE OPTIONS       *******************************************/
#define GMTK_ARG_CONTINUOUS_RANDOM_VAR_OPTIONS
#define GMTK_ARG_VAR_FLOOR
#define GMTK_ARG_VAR_FLOOR_ON_READ


/*************************   GENERAL OPTIONS                          *******************************************/

#define GMTK_ARG_GENERAL_OPTIONS
#define GMTK_ARG_SEED
#define GMTK_ARG_VERB
#define GMTK_ARG_VERSION
#define GMTK_ARG_HELP


#define GMTK_ARGUMENTS_DEFINITION
#include "GMTK_Arguments.h"
#undef GMTK_ARGUMENTS_DEFINITION


Arg Arg::Args[] = {


#define GMTK_ARGUMENTS_DOCUMENTATION
#include "GMTK_Arguments.h"
#undef GMTK_ARGUMENTS_DOCUMENTATION

  Arg("dVeCptName", Arg::Req, dVeCptName, "Name of deep VE CPT to exercise"),

  // final one to signal the end of the list
  Arg()

};



/*
 * definition of needed global arguments
 */
RAND rnd(seedme);
GMParms GM_Parms;
#if 0
ObservationMatrix globalObservationMatrix;
#endif

FileSource *gomFS;
ObservationSource *globalObservationMatrix;

#include <iostream>

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
    "\nThis program trains deep neural networks for use with DeepVECPTs.\n");
  if(!parse_was_ok) {
    Arg::usage(); exit(-1);
  }


#define GMTK_ARGUMENTS_CHECK_ARGS
#include "GMTK_Arguments.h"
#undef GMTK_ARGUMENTS_CHECK_ARGS


  /////////////////////////////////////////////


  // Setup to read observation files
  infoMsg(IM::Max,"Opening Files ...\n");
  gomFS = instantiateFileSource();
  globalObservationMatrix = gomFS;
  infoMsg(IM::Max,"Finished opening files.\n");

  // Read in all GMTK object definitions (matrices, neural networks, etc.)
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
  
  printf("Finished reading in all parameters and structures\n");

  gomFS->openSegment(0);

  // Look up the Deep VE CPT to train 
  string dVeCptNameStr(dVeCptName);
  if (GM_Parms.deepVECptsMap.find(dVeCptNameStr) == GM_Parms.deepVECptsMap.end()) {
    error("Error: No Deep VE CPT named '%s' found\n", dVeCptName);
  }
  DeepVECPT *dvecpt = GM_Parms.deepVECpts[ GM_Parms.deepVECptsMap[dVeCptNameStr] ];

  struct rusage rus; /* starting time */
  struct rusage rue; /* ending time */
  getrusage(RUSAGE_SELF,&rus);

  // Setup TrainingSchedule to create training instances from observation files in desired order

  unsigned radius = dvecpt->windowRadius();
  unsigned offset = dvecpt->obsOffset();
  unsigned nfs    = dvecpt->numFeaturesPerFrame();
  unsigned cardinality = dvecpt->getDeepNN()->numOutputs();

  gomFS->setMinPastFrames( radius );
  gomFS->setMinFutureFrames( radius );

  if (gomFS->numSegments()==0)
    error("ERROR: no segments are available in observation file");
  Range* dcdrng = new Range(dcdrng_str,0,gomFS->numSegments());
  if (dcdrng->length() <= 0) {
    infoMsg(IM::Default,"Decoding range '%s' specifies empty set. Exiting...\n",
	  dcdrng_str);
    exit_program_with_status(0);
  }
  Range::iterator* dcdrng_it = new Range::iterator(dcdrng->begin());
  while (!dcdrng_it->at_end()) {
    const unsigned segment = (unsigned)(*(*dcdrng_it));
    if (gomFS->numSegments() < (segment+1)) 
      error("ERROR: only %d segments in file, segment must be in range [%d,%d]\n",
	    gomFS->numSegments(),
	    0,gomFS->numSegments()-1);

    infoMsg(IM::Max,"Loading segment %d ...\n",segment);
    const unsigned numFrames = GM_Parms.setSegment(segment);
    infoMsg(IM::Max,"Finished loading segment %d with %d frames.\n",segment,numFrames);
    
    RVInfo rvinfo;
    ObsDiscRV parent(rvinfo, 0, 2);
    assert(parent.frame() == 0);
    for (unsigned frame=0; frame < numFrames; frame+=1) {
      assert(parent.frame() == frame);
      printf("%u %u:", segment, frame);
      for (unsigned i=0; i < cardinality; i+=1) {
	printf(" %f", dvecpt->applyNN(i, &parent).unlog());
      }
      printf("\n");
      parent.adjustFrameBy(1);
    }
    (*dcdrng_it)++;
  }
  getrusage(RUSAGE_SELF,&rue);
  if (IM::messageGlb(IM::Default)) { 
    infoMsg(IM::Default,"### Final time (seconds) just for DMLP training stage: ");
    double userTime,sysTime;
    reportTiming(rus,rue,userTime,sysTime,stdout);
  }
  exit_program_with_status(0);
}

