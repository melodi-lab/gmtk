/*-
 * gmtkEMtrain.cc
 *     Train up a GM using EM
 *
 * Written by:
 *   Jeff Bilmes <bilmes@ee.washington.edu>
 *   Geoffrey Zweig <gzweig@us.ibm.com>
 *
 * Copyright (c) 2001, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle, and Jeff Bilmes make no representations about
 * the suitability of this software for any purpose.  It is provided
 * "as is" without express or implied warranty.
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

#include "spi.h"

#include "GMTK_FileParser.h"
#include "GMTK_RandomVariable.h"
#include "GMTK_DiscreteRandomVariable.h"
#include "GMTK_ContinuousRandomVariable.h"
#include "GMTK_GM.h"
#include "GMTK_GMParms.h"
#include "GMTK_ObservationMatrix.h"

VCID("$Header$");

// the file of observation

bool seedme = false;
float pruneRatio=0.0;
char *obsFileName;
char *strFileName;
char *parmsFileName;

unsigned maxEMIterations=3;

bool randomizeParams = true;
bool enem = false;

ARGS ARGS::Args[] = {

 ARGS("obsFile",ARGS::Req,obsFileName,"File containing observations"),
 ARGS("parmsFile",ARGS::Req,parmsFileName,"GM Parms File"), 
 ARGS("strFile",ARGS::Req,strFileName,"GM Structure File"),

 ARGS("seed",ARGS::Opt,seedme,"Seed the RN generator"),

 ARGS("maxEmIters",ARGS::Opt,maxEMIterations,"Max number of EM iterations to do"),

 ARGS("pruneRatio",ARGS::Opt,pruneRatio,"Pruning Ratio, values less than this*max are pruned"),
 ARGS("random",ARGS::Opt,randomizeParams,"Randomize the parameters"),
 ARGS("enem",ARGS::Opt,enem,"Run enumerative EM"),

 ARGS()

};


RAND rnd(seedme);
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


  ARGS::parse(argc,argv);


  if (seedme)
    rnd.seed();

  /////////
  // read in all parameters
  iDataStreamFile parmsFile(parmsFileName);
  GM_Parms.readAll(parmsFile);

  /////////////////////////////
  // read in the structure

  /////////////////////////////
  // read in the structure of the GM, this will
  // die if the file does not exist.
  FileParser fp(strFileName);
  // parse the file
  fp.parseGraphicalModel();
  // create the rv variable objects
  fp.createRandomVariableGraph();
  // make sure that there are no directed loops in the graph
  // by imposing the S,SE,E,NE constrains
  fp.ensureS_SE_E_NE();
  // link the RVs with the parameters that are contained in
  // the bn1_gm.dt file.
  fp.associateWithDataParams();

  // now associate the RVs with a GM
  GMTK_GM gm;
  fp.addVariablesToGM(gm);

  gm.setExampleStream(obsFileName);

  gm.verifyTopologicalOrder();

  gm.GM2CliqueChain();
  // gm.showCliques();

  if (randomizeParams) {
    printf("randomizing parameters\n");
    GM_Parms.makeRandom();
    printf("writing the random parameters to random.gmp\n");
    oDataStreamFile of("random.gmp");
    GM_Parms.writeAll(of);
  }

  gm.setupForVariableLengthUnrolling(fp.firstChunkFrame(),fp.lastChunkFrame());
  if (enem) {
    cout << "\n\n\n WARNING: Doing enumerative EM!!! \n\n";
    gm.enumerativeEM(maxEMIterations);
  } else {
    gm.cliqueChainEM(maxEMIterations, pruneRatio);
  }

  return 0;  

}
