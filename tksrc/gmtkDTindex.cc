/*
 * gmtkJT.cc
 * produce a junction tree
 *
 * Written by Chris Bartels & Jeff Bilmes <bilmes@ee.washington.edu>
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

VCID("$Header$");

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
#include "GMTK_RngDecisionTree.h"

/*
 * command line arguments
 */

static char *DTFiles           = NULL;
static char *inputMasterFile   = NULL;
static char *cppCommandOptions = NULL;
static bool print_version_and_exit = false;

Arg Arg::Args[] = {

  /////////////////////////////////////////////////////////////
  // input parameter/structure file handling

  Arg("decisionTreeFiles", Arg::Opt, DTFiles, "List of decision tree files"),

  Arg("inputMasterFile", Arg::Opt, inputMasterFile,
    "Input file of multi-level master CPP processed GM input parameters"),

  Arg("cppCommandOptions", Arg::Opt, cppCommandOptions,
    "Command line options to give to cpp"),

  Arg("version", Arg::Opt, print_version_and_exit,
    "Print GMTK version number and exit."),

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
  
  if (print_version_and_exit) {
    printf("%s\n",gmtk_version_id);
    exit(0);
  }
  
  if(!parse_was_ok) {
    Arg::usage(); exit(-1);
  }


  ////////////////////////////////////////////
  // Write index files for all decision trees in the master file 
  if (inputMasterFile != NULL) {
    iDataStreamFile pf(inputMasterFile,false,true,cppCommandOptions);

    GM_Parms.read(pf);
    GM_Parms.writeDecisionTreeIndexFiles();
  }

  ////////////////////////////////////////////
  // Write index files for all decision trees named in the list 
  if (DTFiles != NULL) {

    string DT_file_name;
    unsigned int i;

    i = 0;
    while( DTFiles[i]!='\0' ) { 

      RngDecisionTree decision_tree;

      ////////////////////////////////////////////
      // Skip white space 
      while( (DTFiles[i]==' ') || (DTFiles[i]=='\t') ) { 
        ++i;
      }
 
      ////////////////////////////////////////////
      // Get the file name, create a tree object, and write the index file 
      DT_file_name = "";
      while( (DTFiles[i]!=' ') && (DTFiles[i]!='\t') && (DTFiles[i]!='\0') ) {
        DT_file_name += DTFiles[i];
        ++i;
      }

      decision_tree.initializeIterableDT(DT_file_name);
      decision_tree.writeIndexFile();
    }
  }

  exit_program_with_status(0);
}

