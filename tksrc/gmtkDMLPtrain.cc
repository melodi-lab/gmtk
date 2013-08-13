/*
 * gmtkDMLPtrain.cc
 *
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (C) 2013 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 *
 */


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <float.h>
#include <assert.h>


#include "DBN.h"

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

#include "MNIST.h"

#define GMTK_ARG_OBS_FILES
/****************************      FILE RANGE OPTIONS             ***********************************************/
#define GMTK_ARG_FILE_RANGE_OPTIONS
#define GMTK_ARG_TRRNG
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


/*************************   DEEP MLP TRAINING OPTIONS                *******************************************/

#define GMTK_ARG_DMLP_TRAINING_OPTIONS
#define GMTK_ARG_DMLP_TRAINING_PARAMS


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


  infoMsg(IM::Max,"Opening Files ...\n");
  gomFS = instantiateFileSource();
  globalObservationMatrix = gomFS;
  infoMsg(IM::Max,"Finished opening files.\n");

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

#if 0
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
#endif
  
  
  printf("Finished reading in all parameters and structures\n");

  gomFS->openSegment(0);

  DeepVECPT *cpt = NULL;
  for (unsigned i=0; i < GM_Parms.deepVECpts.size(); i+=1) {
    cpt = GM_Parms.deepVECpts[i];
    if (cpt->name().compare(DVECPTName) == 0) break;
  }

  if (!cpt || cpt->name().compare(DVECPTName) != 0) {
    error("Error: No Deep VE CPT named '%s' found\n", DVECPTName);
  }

  printf("Total number of trainable parameters in input files = %u\n",
	 cpt->totalNumberDMLPParameters());

  struct rusage rus; /* starting time */
  struct rusage rue; /* ending time */
  getrusage(RUSAGE_SELF,&rus);

  int inputSize = (int)cpt->numInputs();
  int numLayers = (int)cpt->numLayers();
  int outputSize = (int)cpt->numOutputs();
  
  vector<int> hiddenSize(numLayers);
  for (unsigned i=0; i < numLayers; i+=1)
    hiddenSize[i] = (int)cpt->layerOutputs(i);

  vector<Layer::ActFunc> hActFunc(numLayers);
  for (unsigned i=0; i < numLayers; i+=1) {
    switch (cpt->getSquashFn(i)) {
    case DeepVECPT::SOFTMAX: 
      if (i != numLayers - 1) {
	error("ERROR: gmtkDMLPtrain only supports softmax for the output layer\n");
      }
      hActFunc[i] = Layer::ActFunc(Layer::ActFunc::LINEAR);
      break;
    case DeepVECPT::LOGISTIC: 
      if (i == numLayers - 1) {
	error("ERROR: gmtkDMLPtrain only supports linear or softmax for the output layer\n");
      }
      hActFunc[i] = Layer::ActFunc(Layer::ActFunc::LOG_SIG, (float)cpt->getBeta(i)); 
      assert((float)cpt->getBeta(i) == 1.0f);
      break;
    case DeepVECPT::TANH: 
      if (i == numLayers - 1) {
	error("ERROR: gmtkDMLPtrain only supports linear or softmax for the output layer\n");
      }
      hActFunc[i] = Layer::ActFunc(Layer::ActFunc::TANH); 
      break;
    case DeepVECPT::ODDROOT: 
      if (i == numLayers - 1) {
	error("ERROR: gmtkDMLPtrain only supports linear or softmax for the output layer\n");
      }
      hActFunc[i] = Layer::ActFunc(Layer::ActFunc::CUBIC); 
      break;
    case DeepVECPT::LINEAR:
      hActFunc[i] = Layer::ActFunc(Layer::ActFunc::LINEAR);
      break;
    case DeepVECPT::RECTLIN:
      if (i == numLayers - 1) {
	error("ERROR: gmtkDMLPtrain only supports linear or softmax for the output layer\n");
      }
      if (pretrainMode != DBN::NONE) {
	error("ERROR: gmtkDMLPtrain only supports rectified linear activation functions with -pretrainType none\n");
      }
      hActFunc[i] = Layer::ActFunc(Layer::ActFunc::RECT_LIN);
      break;
    default: 
      error("Error: unknown activation function\n");
    }
  }
  DBN dbn(numLayers, inputSize, hiddenSize, outputSize, iActFunc, hActFunc);

  vector<DBN::HyperParams> pretrainHyperParams(numLayers);
  for (int j = 0; j < numLayers; j+=1) {
    pretrainHyperParams[j].initStepSize     = ptInitStepSize;
    pretrainHyperParams[j].maxMomentum      = ptMaxMomentum;
    pretrainHyperParams[j].maxUpdate        = ptMaxUpdate;
    pretrainHyperParams[j].l2               = ptL2;
    pretrainHyperParams[j].numUpdates       = ptNumUpdates;
    pretrainHyperParams[j].numAnnealUpdates = ptNumAnnealUpdates;
    pretrainHyperParams[j].miniBatchSize    = ptMiniBatchSize;
    pretrainHyperParams[j].checkInterval    = ptCheckInterval;
    pretrainHyperParams[j].dropout          = ptDropout;
    pretrainHyperParams[j].pretrainType     = pretrainMode;
  }

  DBN::HyperParams bpHyperParams;
  bpHyperParams.initStepSize     = bpInitStepSize;
  bpHyperParams.maxMomentum      = bpMaxMomentum;
  bpHyperParams.maxUpdate        = bpMaxUpdate;
  bpHyperParams.l2               = bpL2;
  bpHyperParams.numUpdates       = bpNumUpdates;
  bpHyperParams.numAnnealUpdates = bpNumAnnealUpdates;
  bpHyperParams.miniBatchSize    = bpMiniBatchSize;
  bpHyperParams.checkInterval    = bpCheckInterval;
  bpHyperParams.dropout          = bpDropout;

  unsigned radius = cpt->windowRadius();
  gomFS->setMinPastFrames( radius );
  gomFS->setMinFutureFrames( radius );
  
  Range* trrng = new Range(trrng_str,0,gomFS->numSegments());
  if (trrng->length() <= 0) {
    error("Error: training range '%s' specifies empty set. Exiting...\n", trrng_str);
  }

  unsigned stride = gomFS->stride();
  unsigned numInstances = 0;
  Range::iterator* trrng_it = new Range::iterator(trrng->begin());
  while (!trrng_it->at_end()) {
    const unsigned segment = (unsigned)(*(*trrng_it));
    if (gomFS->numSegments() < (segment+1)) 
      error("ERROR: only %u segments in file, segment must be in range [%u,%u]\n",
	    gomFS->numSegments(),
	    0,gomFS->numSegments()-1);
    if (!gomFS->openSegment(segment))
      error("ERROR: unable to open segment %u\n", segment);

    numInstances += gomFS->numFrames();
    (*trrng_it)++;
  }
  delete trrng_it;

  double *doubleObsData = new double[inputSize * numInstances];
  double *p = doubleObsData;
  double *doubleObsLabel = new double[outputSize * numInstances];
  double *q = doubleObsLabel;
  unsigned obsOffset = cpt->obsOffset();

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
  trrng_it = new Range::iterator(trrng->begin());
  while (!trrng_it->at_end()) {
    const unsigned segment = (unsigned)(*(*trrng_it));
    infoMsg(IM::Max,"Loading segment %u ...\n",segment);
    if (!gomFS->openSegment(segment))
      error("ERROR: unable to open segment %u\n", segment);
    unsigned numFrames = gomFS->numFrames();

    for (unsigned i = 0; i < numFrames; i+=1) {
      Data32 const *obsData = gomFS->loadFrames(i, 1);
      for (int w = -radius; w < radius + 1; w+=1) {
	for (unsigned j=0; j < inputSize; j+=1) {
	  *(p++) = (double)( *((float *)(obsData + w * stride) + obsOffset + j) );
	}
      }
      if (oneHot) {
	for (unsigned j=0; j < outputSize; j+=1) {
	  unsigned label = *((unsigned *)obsData + labelOffset);
	  if ( label >= outputSize ) {
	    error("ERROR: oneHot label %u is too large for output size %u at frame %u in segment %u\n",
		  label, outputSize, i, segment);
	  }
	  *(q++) = (j == label)  ?  1.0 : 0.0;
	}
      } else {
	for (unsigned j=0; j < outputSize; j+=1) {
	  *(q++) = (double)( *((float *)obsData + labelOffset + j) );
	}
      }
      infoMsg(IM::Max,"Finished loading segment %u with %u frames.\n",segment,numFrames);
    }
    (*trrng_it)++;
  }

  Matrix   trainData(doubleObsData,  inputSize,  numInstances, inputSize,  false);
  Matrix trainLabels(doubleObsLabel, outputSize, numInstances, outputSize, false);
  
  DBN::ObjectiveType objType = ( cpt->getSquashFn(numLayers-1) == DeepVECPT::SOFTMAX ) ? DBN::SOFT_MAX : DBN::SQ_ERR;
  dbn.Train(trainData, trainLabels, objType, false, pretrainHyperParams, bpHyperParams);
  delete[] doubleObsLabel;
  delete[] doubleObsData;
  
  vector<DoubleMatrix *> layerMatrix = cpt->getMatrices();
  assert(layerMatrix.size() == numLayers);
  for (unsigned layer=0; layer < numLayers; layer+=1) {
    Matrix const W = dbn.getWeights(layer);
    cpt->setParams(layer, W.Start(), W.Ld(), dbn.getBias(layer).Start());
  }
  
  if (outputTrainableParameters != NULL) {
    char buff[2048];
    copyStringWithTag(buff,outputTrainableParameters,
		      CSWT_EMPTY_TAG,2048);
    oDataStreamFile of(buff,binOutputTrainableParameters);
    GM_Parms.writeTrainable(of);
  }

  // also write according to output master
  GM_Parms.write(outputMasterFile,cppCommandOptions);  

  getrusage(RUSAGE_SELF,&rue);
  if (IM::messageGlb(IM::Default)) { 
    infoMsg(IM::Default,"### Final time (seconds) just for DMLP training stage: ");
    double userTime,sysTime;
    reportTiming(rus,rue,userTime,sysTime,stdout);
  }

  exit_program_with_status(0);
}

