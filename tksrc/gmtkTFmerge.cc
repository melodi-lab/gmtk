/*
 * gmtkTFmerge.cc
 *    Merges optionally the P, C, and E triangulated partitions from up to three
 *    trifiles all of which must have the same boundary and same underlying .str file.
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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <float.h>
#include <assert.h>

#include "general.h"
#include "error.h"
#include "debug.h"
#include "rand.h"
#include "arguments.h"
#include "ieeeFPsetup.h"
#include "spi.h"
#include "version.h"

VCID("$Header$");

#include "GMTK_FileParser.h"
#include "GMTK_RandomVariable.h"
#include "GMTK_DiscreteRandomVariable.h"
#include "GMTK_ContinuousRandomVariable.h"
#include "GMTK_GMTemplate.h"
#include "GMTK_GMParms.h"
#include "GMTK_ObservationMatrix.h"
#include "GMTK_MixtureCommon.h"
#include "GMTK_GaussianComponent.h"
#include "GMTK_MeanVector.h"
#include "GMTK_DiagCovarVector.h"
#include "GMTK_DlinkMatrix.h"
#include "GMTK_ProgramDefaultParms.h"
#include "GMTK_BoundaryTriangulate.h"
#include "GMTK_JunctionTree.h"

/*
 * command line arguments
 */
static char *strFileName=NULL;
static double varFloor = GMTK_DEFAULT_VARIANCE_FLOOR;
static char *cppCommandOptions = NULL;
static bool loadParameters = false;

static char *inputMasterFile=NULL;
static char *inputTrainableParameters=NULL;
static char *outputTriangulatedFile=NULL;
static bool binInputTrainableParameters=false;

static char* Ptrifile = false;
static char* Ctrifile = false;
static char* Etrifile = false;

static unsigned numBackupFiles = 7;
static unsigned verbosity = IM::Default;
// static bool printResults = false;
static int allocateDenseCpts=0;

Arg Arg::Args[] = {

  /////////////////////////////////////////////////////////////
  // input parameter/structure file handling

  Arg("cppCommandOptions",Arg::Opt,cppCommandOptions,"Additional CPP command line"),
  Arg("strFile",Arg::Req,strFileName,"Graphical Model Structure File"),
  Arg("inputMasterFile",Arg::Opt,inputMasterFile,"Input file of multi-level master CPP processed GM input parameters"),
  Arg("inputTrainableParameters",Arg::Opt,inputTrainableParameters,"File of only and all trainable parameters"),

  /////////////////////////////////////////////////////////////
  // Triangulation file Options

  Arg("outputTriangulatedFile",Arg::Opt,outputTriangulatedFile,"File name to write resulting triangulation to"),

  Arg("Ptrifile",
      Arg::Opt,Ptrifile,
      "Tri-file to obtain triangulation for P partition"),
  Arg("Ctrifile",
      Arg::Opt,Ctrifile,
      "Tri-file to obtain triangulation for C partition"),
  Arg("Etrifile",
      Arg::Opt,Etrifile,
      "Tri-file to obtain triangulation for E partition"),

  Arg("numBackupFiles",Arg::Opt,numBackupFiles,"Number of backup output .trifiles (_bak0,_bak1,etc.) to keep."),
  // Arg("printResults",Arg::Opt,printResults,"Print information about result of final triangulation."),
  Arg("loadParameters",Arg::Opt,loadParameters,"Also load in all trainable parameters."),

  /////////////////////////////////////////////////////////////
  // General Options

  Arg("allocateDenseCpts",Arg::Opt,allocateDenseCpts,"Automatically allocate any undefined CPTs. arg = -1, no read params, arg = 0 noallocate, arg = 1 means use random initial CPT values. arg = 2, use uniform values"),
  Arg("verbosity",Arg::Opt,verbosity,"Verbosity (0 <= v <= 100) of informational/debugging msgs"),
  // final one to signal the end of the list
  Arg()

};

#define MYBS(x) ((x)?"T":"F")



/*
 *
 * definition of needed global arguments
 *
 */
RAND rnd(false);
GMParms GM_Parms;
ObservationMatrix globalObservationMatrix;

/*
 *
 * backupTriFile:
 *    Make a backup copys by renaming the file triFile since it might
 *    be a mistake to delete it and since these file scan take
 *    a while to generate. TODO: put this in general.cc
 *
 */
void
backupTriFile(const string &triFile) 
{
  if (numBackupFiles == 0)
    return;
  for (unsigned bk_num=(numBackupFiles-1);bk_num>0;bk_num--) {
    char buff[1024];
    sprintf(buff,"%d",bk_num-1);
    string curFile =  triFile + "_bak" + buff;
    if (fsize(curFile.c_str()) == 0)
      continue;

    sprintf(buff,"%d",bk_num);
    string backupFile = triFile + "_bak" + buff;
    if (rename(curFile.c_str(),backupFile.c_str()) != 0)
      infoMsg(IM::Warning,"Warning: could not create backup triangulation file %s.\n",
	      backupFile.c_str());
  }
  if (fsize(triFile.c_str()) == 0)
    return;
  string backupFile = triFile + "_bak0";
  if (rename(triFile.c_str(),backupFile.c_str()) != 0)
    infoMsg(IM::Warning,"Warning: could not create backup triangulation file %s.\n",
	    backupFile.c_str());
}



int
main(int argc,char*argv[])
{
  ////////////////////////////////////////////
  // set things up so that if an FP exception
  // occurs such as an "invalid" (NaN), overflow
  // or divide by zero, we actually get a FPE
  ieeeFPsetup();
  set_new_handler(memory_error);

  ////////////////////////////////////////////
  // parse arguments
  Arg::parse(argc,argv);
  (void) IM::setGlbMsgLevel(verbosity);


  if (Ptrifile == NULL) {
    if (Etrifile == NULL) {
      if (Ctrifile == NULL) {
	error("Error: All three trifiles have been unspecified on command line.\n");
      } else {
	Ptrifile = Ctrifile = Etrifile;
	infoMsg(IM::Warning,"Warning: Using Etrifile for both Ctrifile and Ptrifile.\n");
      }
    } else {
      Ptrifile = Etrifile;
      infoMsg(IM::Warning,"Warning: Using Etrifile for Ptrifile\n");
    }
  }
  if (Ctrifile == NULL) {
    Ctrifile = Ptrifile;
    infoMsg(IM::Warning,"Warning: Using Ptrifile for Ctrifile.\n");
  }
  if (Etrifile == NULL) {
    Etrifile = Ptrifile;    
    infoMsg(IM::Warning,"Warning: Using Ptrifile for Etrifile.\n");
  }

  MixtureCommon::checkForValidRatioValues();
  MeanVector::checkForValidValues();
  DiagCovarVector::checkForValidValues();
  DlinkMatrix::checkForValidValues();

  ////////////////////////////////////////////
  // set global variables/change global state from args
  GaussianComponent::setVarianceFloor(varFloor);

  // if (chunkSkip > maxNumChunksInBoundary)
  //  error("ERROR: Must have S<=M at this time.\n");

  /////////////////////////////////////////////
  if (loadParameters) {
    // read in all the parameters
    if (inputMasterFile) {
      // flat, where everything is contained in one file, always ASCII
      iDataStreamFile pf(inputMasterFile,false,true,cppCommandOptions);
      GM_Parms.read(pf);
    }
    if (inputTrainableParameters) {
      // flat, where everything is contained in one file
      iDataStreamFile pf(inputTrainableParameters,binInputTrainableParameters,true,cppCommandOptions);
      GM_Parms.readTrainable(pf);
    }
  }
  GM_Parms.finalizeParameters();

  /////////////////////////////
  // read in the structure of the GM, this will
  // die if the file does not exist.
  FileParser fp(strFileName,cppCommandOptions);
  infoMsg(IM::Tiny,"Finished reading in structure\n");

  // parse the file
  fp.parseGraphicalModel();
  // create the rv variable objects
  fp.createRandomVariableGraph();
  // Make sure that there are no directed loops in the graph.
  fp.ensureValidTemplate();

  if (loadParameters) {
    // link the RVs with the parameters that are contained in
    // the bn1_gm.dt file.
    if (allocateDenseCpts >= 0) {
      if (allocateDenseCpts == 0)
	fp.associateWithDataParams(FileParser::noAllocate);
      else if (allocateDenseCpts == 1)
	fp.associateWithDataParams(FileParser::allocateRandom);
      else if (allocateDenseCpts == 2)
	fp.associateWithDataParams(FileParser::allocateUniform);
      else
	error("Error: command line argument '-allocateDenseCpts d', must have d = {0,1,2}\n");
    }
  }


  // create with bogus M and S values for now, to be overwritten
  // once we read in the partitions.
  GMTemplate gm_template_for_P(fp,GMTEMPLATE_UNINITIALIZED_MS,GMTEMPLATE_UNINITIALIZED_MS);
  GMTemplate gm_template_for_C(fp,GMTEMPLATE_UNINITIALIZED_MS,GMTEMPLATE_UNINITIALIZED_MS);
  GMTemplate gm_template_for_E(fp,GMTEMPLATE_UNINITIALIZED_MS,GMTEMPLATE_UNINITIALIZED_MS);


  {
    iDataStreamFile is(Ptrifile,false,false);
    if (!fp.readAndVerifyGMId(is))
      error("ERROR: triangulation P tri-file '%s' does not match graph given in structure file '%s'\n",Ptrifile,strFileName);
    gm_template_for_P.readPartitions(is);
    gm_template_for_P.readMaxCliques(is);
    gm_template_for_P.triangulatePartitionsByCliqueCompletion();
  }
  {
    iDataStreamFile is(Ctrifile,false,false);
    if (!fp.readAndVerifyGMId(is))
      error("ERROR: triangulation C tri-file '%s' does not match graph given in structure file '%s'\n",Ctrifile,strFileName);
    gm_template_for_C.readPartitions(is);
    gm_template_for_C.readMaxCliques(is);
    gm_template_for_C.triangulatePartitionsByCliqueCompletion();
  }
  {
    iDataStreamFile is(Etrifile,false,false);
    if (!fp.readAndVerifyGMId(is))
      error("ERROR: triangulation E tri-file '%s' does not match graph given in structure file '%s'\n",Etrifile,strFileName);
    gm_template_for_E.readPartitions(is);
    gm_template_for_E.readMaxCliques(is);
    gm_template_for_E.triangulatePartitionsByCliqueCompletion();
  }


  // so they all match the str file, but we need to make sure they all
  // have the same S, M, and boundary.

  if (gm_template_for_E.maxNumChunksInBoundary() != gm_template_for_C.maxNumChunksInBoundary())
    error("Error: E trifile has M = %d != %d = M for C trifile\n");
  if (gm_template_for_C.maxNumChunksInBoundary() != gm_template_for_P.maxNumChunksInBoundary())
    error("Error: C trifile has M = %d != %d = M for P trifile\n");
  if (gm_template_for_E.chunkSkip() != gm_template_for_C.chunkSkip())
    error("Error: E trifile has S = %d != %d = S for C trifile\n");
  if (gm_template_for_C.chunkSkip() != gm_template_for_P.chunkSkip())
    error("Error: C trifile has S = %d != %d = S for P trifile\n");

  
  // chose C for the final destination version.
  gm_template_for_C.P.setCliquesFromAnotherPartition(gm_template_for_P.P);
  gm_template_for_C.E.setCliquesFromAnotherPartition(gm_template_for_E.E);

  // make sure triangulated
  BoundaryTriangulate triangulator(fp,gm_template_for_C.maxNumChunksInBoundary(),gm_template_for_C.chunkSkip(),1.0);
  triangulator.ensurePartitionsAreChordal(gm_template_for_C);

  string output_tri_file;
  if (outputTriangulatedFile == NULL) {
    output_tri_file = string("merged-output-trifile") + GMTemplate::fileExtension;
  } else {
    output_tri_file = string(outputTriangulatedFile);
  }

  char buff[2*4096];
  sprintf(buff,"Created from gmtkTFmerge with input files for P (%s), C (%s), and E(%s)\n",
	  Ptrifile,Ctrifile,Etrifile);
  string str = buff;

  backupTriFile(output_tri_file);
  oDataStreamFile os(output_tri_file.c_str());
  fp.writeGMId(os);
  gm_template_for_C.writePartitions(os,str);
  gm_template_for_C.writeMaxCliques(os);

  exit_program_with_status(0);
}
