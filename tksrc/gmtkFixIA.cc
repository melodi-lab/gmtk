/*
 * gmtkJT.cc
 *   compute probability of evidence or clique posteriors
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *            Richard Rogers <rprogers@uw.edu>
 *
 * Copyright (C) 2001, 2015 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
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
#include <time.h>

#include "general.h"
#include "error.h"
#include "debug.h"
#include "rand.h"
#include "arguments.h"
#include "ieeeFPsetup.h"

#include "GMTK_WordOrganization.h"

#include "GMTK_FileParser.h"
#include "GMTK_RV.h"
#include "GMTK_DiscRV.h"
#include "GMTK_ContRV.h"
#include "GMTK_GMTemplate.h"
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
#include "GMTK_Stream.h"
#include "GMTK_ZeroCliqueException.h"

// Supported inference tasks
#include "GMTK_ProbEvidenceTask.h"
#include "GMTK_ForwardBackwardTask.h"

// Section scheduling
#include "GMTK_LinearSectionScheduler.h"
#include "GMTK_IslandSectionScheduler.h"
#include "GMTK_ArchipelagosSectionScheduler.h"

// Supported within-sectin inference algorithms
#include "GMTK_SparseJoinInference.h"
#include "GMTK_PedagogicalInference.h"

#include "GMTK_MixtureCommon.h"
#include "GMTK_GaussianComponent.h"
#include "GMTK_MeanVector.h"
#include "GMTK_DiagCovarVector.h"
#include "GMTK_Dlinks.h"
#include "GMTK_DlinkMatrix.h"
#include "GMTK_ProgramDefaultParms.h"
#include "GMTK_BoundaryTriangulate.h"
#include "GMTK_MaxClique.h"


VCID(HGID)



/*************************   INPUT TRAINABLE PARAMETER FILE HANDLING  *******************************************/
#define GMTK_ARG_INPUT_TRAINABLE_FILE_HANDLING
#define GMTK_ARG_CPP_CMD_OPTS
#define GMTK_ARG_INPUT_MASTER_FILE
#define GMTK_ARG_DLOPEN_MAPPERS
#define GMTK_ARG_INPUT_TRAINABLE_PARAMS
#define GMTK_ARG_ALLOC_DENSE_CPTS

/*************************   INPUT STRUCTURE PARAMETER FILE HANDLING  *******************************************/
#define GMTK_ARG_INPUT_MODEL_FILE_HANDLING
#define GMTK_ARG_STR_FILE
#define GMTK_ARG_INPUT_TRI_FILE
#define GMTK_ARG_CHECK_TRI_FILE_CARD
#define GMTK_ARG_DLOPEN_MAPPERS
#define GMTK_ARG_INPUT_TRAINABLE_PARAMS
#define GMTK_ARG_OUTPUT_TRI_FILE

/****************************         GENERAL OPTIONS             ***********************************************/
#define GMTK_ARG_GENERAL_OPTIONS
#define GMTK_ARG_SEED
#define GMTK_ARG_SKIP_STARTUP_CHECKS
#define GMTK_ARG_VERB
#define GMTK_ARG_HELP
#define GMTK_ARG_VERSION

/****************************         INFERENCE OPTIONS           ***********************************************/
#define GMTK_ARG_INFERENCE_OPTIONS
#define GMTK_ARG_CLIQUE_VAR_ITER_ORDERS
#define GMTK_ARG_JT_OPTIONS


/************************            TRIANGULATION OPTIONS             ******************************************/
#define GMTK_ARG_TRIANGULATION_OPTIONS
#define GMTK_ARG_LOAD_PARAMETERS
#define GMTK_ARG_NUM_BACKUP_FILES
#define GMTK_ARG_JTW_UB
#define GMTK_ARG_JT_OPTIONS
#define GMTK_ARG_CROSSOVER_OPTIONS

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
 *
 * definition of needed global arguments
 *
 */
RAND rnd(seedme);
GMParms GM_Parms;
FileSource fileSource;
FileSource *gomFS = &fileSource;
ObservationSource *globalObservationMatrix = &fileSource;


/*
 *
 * backupTriFile:
 *    Make a backup copys by renaming the file triFile since it might
 *    be a mistake to delete it and since these file scan take
 *    a while to generate.
 *
 *  TODO: move this to general.{cc,h}
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
main(int argc,char*argv[]) {
  string input_tri_file, output_tri_file;
{ // use double so that we can destruct objects at end.

  ////////////////////////////////////////////
  // set things up so that if an FP exception
  // occurs such as an "invalid" (NaN), overflow
  // or divide by zero, we actually get a FPE
  ieeeFPsetup();
  set_new_handler(memory_error);

  //  CODE_TO_COMPUTE_ENDIAN;

  ////////////////////////////////////////////
  // parse arguments
  bool parse_was_ok = Arg::parse(argc,(char**)argv,
"\nThis program updates the IA file based on manual edits to\n"
"the section interface factorization\n");
  if(!parse_was_ok) {
    // Arg::usage(); 
    exit(-1);
  }

  infoMsg(IM::Max,"Finished parsing arguments\n");

#define GMTK_ARGUMENTS_CHECK_ARGS
#include "GMTK_Arguments.h"
#undef GMTK_ARGUMENTS_CHECK_ARGS


  /////////////////////////////////////////////
  // read in all the parameters

  dlopenDeterministicMaps(dlopenFilenames, MAX_NUM_DLOPENED_FILES);
  if (inputMasterFile) {
    // flat, where everything is contained in one file, always ASCII
    infoMsg(IM::Max,"Reading master file...\n");
    iDataStreamFile pf(inputMasterFile,false,true,cppCommandOptions);
    GM_Parms.read(pf);
    infoMsg(IM::Max,"Finished reading master file.\n");
  }
  if (inputTrainableParameters) {
    // flat, where everything is contained in one file
    infoMsg(IM::Max,"Reading trainable file...\n");
    iDataStreamFile pf(inputTrainableParameters,binInputTrainableParameters,true,cppCommandOptions);
    GM_Parms.readTrainable(pf);
    infoMsg(IM::Max,"Finished reading trainable file.\n");
  }
  GM_Parms.finalizeParameters();


  /////////////////////////////
  // read in the structure of the GM, this will
  // die if the file does not exist.
  infoMsg(IM::Max,"Reading structure file...\n");
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

  if (loadParameters) {
    // link the RVs with the parameters.
    infoMsg(IM::Max,"Allocating cpts...\n");
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

  //////////////////////////////////////////////////////////////////////
  // Get name of input triangulation
  //////////////////////////////////////////////////////////////////////
  if (inputTriangulatedFile == NULL) {
    input_tri_file = string(strFileName) + GMTemplate::fileExtension;
  }
  else {
    input_tri_file = string(inputTriangulatedFile);
    if (fsize(input_tri_file.c_str()) == 0) {
      error("ERROR: inference architecture file '%s' does not exist or is empty\n",
        input_tri_file.c_str() );
    }
  }

  //////////////////////////////////////////////////////////////////////
  // Name the output trifiles
  //////////////////////////////////////////////////////////////////////
  if (outputTriangulatedFile == NULL) {
    output_tri_file = string(strFileName) + GMTemplate::fileExtension;
  } else {
    output_tri_file = string(outputTriangulatedFile);
  }

  infoMsg(IM::Max,"Creating template...\n");
  GMTemplate gm_template(fp);


  // Instantiate the requested inference algorithm for the time series as a whole

  SectionScheduler *myjt = new SectionScheduler(gm_template, fp, gomFS);


  string tri_file;
  if (inputTriangulatedFile == NULL) 
    tri_file = string(strFileName) + GMTemplate::fileExtension;
  else 
    tri_file = string(inputTriangulatedFile);

  {
    // do this in scope so that is gets deleted now rather than later.
    iDataStreamFile is(tri_file.c_str());
    myjt->updateDataStructures(is, varSectionAssignmentPrior,varCliqueAssignmentPrior, checkTriFileCards);
  }

  printf("writing IA to %s\n", output_tri_file.c_str());
  backupTriFile(output_tri_file);
  myjt->printAllIAInfo(output_tri_file.c_str(), writeComments);

} // close brace to cause a destruct on valid end of program.
 exit_program_with_status(0); 
}

