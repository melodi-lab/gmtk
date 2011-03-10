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
#include "GMTK_ObservationMatrix.h"
#include "GMTK_MixtureCommon.h"
#include "GMTK_GaussianComponent.h"
#include "GMTK_MeanVector.h"
#include "GMTK_DiagCovarVector.h"
#include "GMTK_DlinkMatrix.h"
#include "GMTK_RngDecisionTree.h"


#define GMTK_ARG_INPUT_MASTER_FILE_OPT_ARG
#define GMTK_ARG_CPP_CMD_OPTS
#define GMTK_ARG_VERB
#define GMTK_ARG_VERSION

static char *DTFiles           = NULL;


#define GMTK_ARGUMENTS_DEFINITION
#include "GMTK_Arguments.h"
#undef GMTK_ARGUMENTS_DEFINITION


Arg Arg::Args[] = {

#define GMTK_ARGUMENTS_DOCUMENTATION
#include "GMTK_Arguments.h"
#undef GMTK_ARGUMENTS_DOCUMENTATION

  Arg("decisionTreeFiles", Arg::Opt, DTFiles, "List of decision tree files"),

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

