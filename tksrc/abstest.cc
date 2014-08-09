/*
 * abstest.cc
 *
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (C) 2013 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */


#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>

#include "BatchSource.h"
#include "AsynchronousBatchSource.h"

#include "general.h"
#include "error.h"
#include "rand.h"
#include "arguments.h"
#include "debug.h"

#include "GMTK_RV.h"
#include "GMTK_GMParms.h"
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

#include "GMTK_RandomSampleSchedule.h"
#include "GMTK_LinearSchedule.h"
#include "GMTK_ShuffleSchedule.h"
#include "GMTK_PermutationSchedule.h"

#include "GMTK_WordOrganization.h"

VCID(HGID)

unsigned radius=0, labelOffset, obsOffset=0, numFeaturesPerFrame, outputSize;
bool oneHot;
static char const *trainingSchedule = "linear";
char *cppCommandOptions = NULL;
GMParms GM_Parms;
 

#define GMTK_ARG_OBS_FILES
/****************************      FILE RANGE OPTIONS             ***********************************************/
#define GMTK_ARG_FILE_RANGE_OPTIONS
#define GMTK_ARG_TRRNG
#define GMTK_ARG_START_END_SKIP
/************************  OBSERVATION MATRIX TRANSFORMATION OPTIONS   ******************************************/
#define GMTK_ARG_OBS_MATRIX_OPTIONS
#define GMTK_ARG_OBS_MATRIX_XFORMATION

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

Arg("radius", Arg::Opt, radius, "radius"),
Arg("obsOffset", Arg::Opt, obsOffset, "feature position"),
Arg("featuresPerFrame", Arg::Req, numFeaturesPerFrame, "features per frame"),
Arg("labelOffset", Arg::Req, labelOffset, "Position in observation file where output labels start"),
Arg("labelsPerFrame", Arg::Req, outputSize, "labels per frame"),
Arg("oneHot", Arg::Opt, oneHot, "If true, labelOffset is the single discrete correct parent value, "
                                "else the parent distribution starts ate labelOffset"),
Arg("trainingSchedule", Arg::Opt, trainingSchedule, "Order to process training data (linear, random, permute, shuffle)"),

#undef GMTK_ARGUMENTS_DOCUMENTATION



  // final one to signal the end of the list
  Arg()

};



/*
 * definition of needed global arguments
 */
RAND rnd(false);

FileSource *gomFS;
ObservationSource *globalObservationMatrix;


int
main(int argc,char*argv[])
{
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


  // Setup to read observation files
  infoMsg(IM::Max,"Opening Files ...\n");
  gomFS = instantiateFileSource();
  globalObservationMatrix = gomFS;
  infoMsg(IM::Max,"Finished opening files.\n");


  gomFS->openSegment(0);


  // Setup TrainingSchedule to create training instances from observation files in desired order

  gomFS->setMinPastFrames( radius );
  gomFS->setMinFutureFrames( radius );



  if (oneHot) {
    if (labelOffset < gomFS->numContinuous()) {
      error("ERROR: labelOffset (%u) must refer to a discrete feature (the first %u are continuous)\n", 
	    labelOffset, gomFS->numContinuous());
    }
    if (labelOffset >= gomFS->numFeatures()) {
      error("ERROR: labelOffset (%u) is too large for the number of available features (%u)\n",
	    labelOffset, gomFS->numFeatures());
    }
  } else {
    if (labelOffset >= gomFS->numContinuous()) {
      error("ERROR: labelOffset (%u) is too large for the number of continuous features (%u)\n",
	    labelOffset, gomFS->numContinuous());
    }
    if (labelOffset + outputSize > gomFS->numContinuous()) {
      error("ERROR: labelOffset (%u) + number of outputs (%u) is too large for the number of continuous features (%u)\n", 
	    labelOffset, outputSize, gomFS->numContinuous());
    }
  }

  TrainingSchedule *trainSched;
  if (strcasecmp(trainingSchedule, "linear") == 0) {
    trainSched = new LinearSchedule(obsOffset, numFeaturesPerFrame, labelOffset, outputSize, oneHot, radius, 1, gomFS, trrng_str);
  } else if (strcasecmp(trainingSchedule, "random") == 0) {
    trainSched = new RandomSampleSchedule(obsOffset, numFeaturesPerFrame, labelOffset, outputSize, oneHot, radius, 1, gomFS, trrng_str);
  } else if (strcasecmp(trainingSchedule, "permute") == 0) {
    trainSched = new PermutationSchedule(obsOffset, numFeaturesPerFrame, labelOffset, outputSize, oneHot, radius, 1, gomFS, trrng_str);
  } else if (strcasecmp(trainingSchedule, "shuffle") == 0) {
    trainSched = new ShuffleSchedule(obsOffset, numFeaturesPerFrame, labelOffset, outputSize, oneHot, radius, 1, gomFS, trrng_str);
  } else {
    error("ERROR: unknown training schedule '%s', must be one of linear, random, shuffle, or permute\n", trainingSchedule);
  }

#define BATCHSIZE 10
#define QUEUESIZE (10*BATCHSIZE)
  BatchSource *synchSrc = new ScheduleBatchSource(trainSched);
  BatchSource *batchSrc = new AsynchronousBatchSource(synchSrc, QUEUESIZE);

  unsigned batchNum = 0;
  do {
    printf("batch %u {\n", batchNum++);
    Matrix data, labels;
    batchSrc->getBatch(BATCHSIZE, data, labels);
    assert(data.NumC() == labels.NumC());
    assert(data.NumC() > 0);
    assert(data.NumR() == synchSrc->numDataRows());
    assert(labels.NumR() == synchSrc->numLabelRows());
    for (unsigned i=0; i < data.NumC(); i += 1) {
      printf(" ");
      for (unsigned j=0; j < data.NumR(); j += 1) {
	printf(" %7.3f", data.At(j,i));
      }
      printf(" |");
      for (unsigned j=0; j < labels.NumR(); j+=1) {
	printf(" %7.3f", labels.At(j,i));
      }
      printf("\n");
    }
    printf("}\n");
  } while (1);

  exit_program_with_status(0);
}

