/*-
 * gmtkEMtrain.cc
 *     Train up a GM using EM
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
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
unsigned maxEmIterations=30;

ARGS ARGS::Args[] = {

 ARGS("obsFile",ARGS::Req,obsFileName,"File containing observations"),
 ARGS("parmsFile",ARGS::Req,parmsFileName,"GM Parms File"), 
 ARGS("strFile",ARGS::Req,strFileName,"GM Structure File"),

 ARGS("seed",ARGS::Opt,seedme,"Seed the RN generator"),
 ARGS("maxEmIters",ARGS::Opt,maxEmIterations,"Number of EM iterations to do"),
 ARGS("pruneRatio",ARGS::Opt,pruneRatio,"Pruning Ratio"),
 ARGS()

};


RAND rnd(seedme);
GMParms GM_Parms;
ObservationMatrix globalObservationMatrix;


int
main(int argc,char*argv[])
{
  ARGS::parse(argc,argv);

  if (seedme)
    rnd.seed();

  logpr beam=pruneRatio;

  /////
  // open the observation file
  globalObservationMatrix.openFile(obsFileName);

  /////////
  // read in all parameters
  iDataStreamFile parmsFile(parmsFileName);
  GM_Parms.read(parmsFile);

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

  gm.verifyTopologicalOrder();

  gm.GM2CliqueChain();
  // gm.showCliques();

  // printf("randomizing parameters\n");
  // gm.makeRandom();

  logpr previous_likelihood = 0.0;
  for (unsigned emIter=0;emIter<maxEmIterations;emIter++) {
    gm.emStartIteration();
    logpr likelihood = 0.0;
    for (unsigned segment=0;
	 segment<globalObservationMatrix.numSegments();
	 segment++) {
      globalObservationMatrix.loadSegment(segment);

      //////////////////////////////
      // TODO: right now, need to 're-unroll' the network in each case,
      // This is terribly inefficient and it needs to be fixed.
      
      if (gm.chain)
	delete gm.chain;

      if (((globalObservationMatrix.numFrames 
	    - fp.firstChunkframe()
	    - (fp.maxFrame() - fp.lastChunkframe()))
	   % 
	   (fp.lastChunkframe()-fp.firstChunkframe()+1)) != 0) {
	error("ERROR: segment %d, number of frames %d minus %d must be a multiple of %d\n",segment,
	      globalObservationMatrix.numFrames,
	      fp.firstChunkframe()+(fp.maxFrame()-fp.lastChunkframe()),
	      (fp.lastChunkframe()-fp.firstChunkframe()+1));
      }

      printf("Unrolling %d times\n",
	     (globalObservationMatrix.numFrames 
	      - fp.firstChunkframe()
	      - (fp.maxFrame() - fp.lastChunkframe()))
	     /
	     (fp.lastChunkframe()-fp.firstChunkframe()+1));

      gm.unroll(fp.firstChunkframe(), fp.lastChunkframe(), 
		(globalObservationMatrix.numFrames 
		 - fp.firstChunkframe()
		 - (fp.maxFrame() - fp.lastChunkframe()))
		/
		(fp.lastChunkframe()-fp.firstChunkframe()+1)
		);
      if (!gm.chain->computePosteriors(beam)) {
	printf("Skipping segment %d due to 0 probability\n",segment);
	continue;
      }
      likelihood *= gm.chain->dataProb;
    }
    printf("Total data likelihood is %e\n",likelihood.val());
    gm.emEndIteration();
    if (likelihood > previous_likelihood) {
      gm.emSwapCurAndNew();
    }
    previous_likelihood = likelihood;
  }

  return 0;  

}
