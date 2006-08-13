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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <float.h>
#include <assert.h>

// TODO: remove next 2 eventually 
#include <iostream.h>
#include <fstream.h>

#include "general.h"
#include "error.h"
#include "debug.h"
#include "rand.h"
#include "arguments.h"
#include "ieeeFPsetup.h"
#include "version.h"

#include "GMTK_WordOrganization.h"

VCID("$Header$")

#include "GMTK_FileParser.h"
#include "GMTK_RV.h"
#include "GMTK_DiscRV.h"
#include "GMTK_ContRV.h"
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
#include "GMTK_MaxClique.h"


/*****************************   OBSERVATION INPUT FILE HANDLING   **********************************************/
#define GMTK_ARG_OBS_FILES

/*************************   INPUT TRAINABLE PARAMETER FILE HANDLING  *******************************************/
#define GMTK_ARG_CPP_CMD_OPTS
#define GMTK_ARG_INPUT_MASTER_FILE
#define GMTK_ARG_INPUT_TRAINABLE_PARAMS
#define GMTK_ARG_ALLOC_DENSE_CPTS
#define GMTK_ARG_CPT_NORM_THRES

/*************************   INPUT STRUCTURE PARAMETER FILE HANDLING  *******************************************/
#define GMTK_ARG_STR_FILE
#define GMTK_ARG_TRI_FILE
#define GMTK_ARG_CHECK_TRI_FILE_CARD
#define GMTK_ARG_JT_INFO_FILE
#define GMTK_ARG_JTW_UB


/*************************   CONTINUOUS RANDOM VARIABLE OPTIONS       *******************************************/
#define GMTK_ARG_VAR_FLOOR
#define GMTK_ARG_VAR_FLOOR_ON_READ


/*************************          BEAM PRUNING OPTIONS              *******************************************/
#define GMTK_ARG_CBEAM
#define GMTK_ARG_CPBEAM
#define GMTK_ARG_CKBEAM
#define GMTK_ARG_CRBEAM
#define GMTK_ARG_CMBEAM
#define GMTK_ARG_SBEAM

/*************************          MEMORY MANAGEMENT OPTIONS         *******************************************/
#define GMTK_ARG_HASH_LOAD_FACTOR
#define GMTK_ARG_STORE_DETERMINISTIC_CHILDREN
#define GMTK_ARG_CLEAR_CLIQUE_VAL_MEM


/****************************      FILE RANGE OPTIONS             ***********************************************/
#define GMTK_ARG_DCDRNG
#define GMTK_ARG_START_END_SKIP

/****************************         GENERAL OPTIONS             ***********************************************/
#define GMTK_ARG_SEED
#define GMTK_ARG_VERB
#define GMTK_ARG_HELP
#define GMTK_ARG_VERSION


/****************************         INFERENCE OPTIONS           ***********************************************/
#define GMTK_ARG_ISLAND
#define GMTK_ARG_CE_SEP_DRIVEN
#define GMTK_ARG_MIXTURE_CACHE
#define GMTK_ARG_CLIQUE_VAR_ITER_ORDERS
#define GMTK_ARG_JT_OPTIONS
#define GMTK_ARG_VE_SEPS

/************************  OBSERVATION MATRIX TRANSFORMATION OPTIONS   ******************************************/
#define GMTK_ARG_OBS_MATRIX_XFORMATION

/************************            DECODING OPTIONS                  ******************************************/
#define GMTK_ARG_DECODING_OPTIONS


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
ObservationMatrix globalObservationMatrix;


int
main(int argc,char*argv[])
{

  ////////////////////////////////////////////
  // set things up so that if an FP exception
  // occurs such as an "invalid" (NaN), overflow
  // or divide by zero, we actually get a FPE
  ieeeFPsetup();
  set_new_handler(memory_error);

  CODE_TO_COMPUTE_ENDIAN;
   
  ////////////////////////////////////////////
  // parse arguments
  bool parse_was_ok = Arg::parse(argc,(char**)argv);
  if(!parse_was_ok) {
    Arg::usage(); exit(-1);
  }


#define GMTK_ARGUMENTS_CHECK_ARGS
#include "GMTK_Arguments.h"
#undef GMTK_ARGUMENTS_CHECK_ARGS

  globalObservationMatrix.openFiles(nfiles,
				    (const char**)&ofs,
				    (const char**)&frs,
				    (const char**)&irs,
				    (unsigned*)&nfs,
				    (unsigned*)&nis,
				    (unsigned*)&ifmts,
				    (bool*)&iswp,
				    startSkip,
				    endSkip,
				    Cpp_If_Ascii,
				    cppCommandOptions,
				    (const char**)&postpr,  //Frame_Range_Str,
				    Action_If_Diff_Num_Frames,
				    Action_If_Diff_Num_Sents,
				    Per_Stream_Transforms,
				    Post_Transforms,
				    Ftr_Combo,
				    (const char**)&sr,
				    (const char**)&prepr,
				    gpr_str
				    );


  /////////////////////////////////////////////
  // read in all the parameters

  if (inputTrainableParameters) {
    // flat, where everything is contained in one file
    iDataStreamFile pf(inputTrainableParameters,binInputTrainableParameters,true,cppCommandOptions);
    GM_Parms.readTrainable(pf);
  }

  if (inputMasterFile) {
    // flat, where everything is contained in one file, always ASCII
    iDataStreamFile pf(inputMasterFile,false,true,cppCommandOptions);
    GM_Parms.read(pf);
  }

  // comment for now Sun Jan 11 09:47:23 2004
  GM_Parms.finalizeParameters();

  /////////////////////////////
  // read in the structure of the GM, this will
  // die if the file does not exist.
  FileParser fp(strFileName,cppCommandOptions);
  infoMsg(IM::Tiny,"Finished reading in all parameters and structures\n");

  // parse the file
  fp.parseGraphicalModel();
  // create the rv variable objects
  fp.createRandomVariableGraph();
  // Make sure that there are no directed loops in the graph.
  fp.ensureValidTemplate();

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



  // make sure that all observation variables work
  // with the global observation stream.
  fp.checkConsistentWithGlobalObservationStream();
  GM_Parms.checkConsistentWithGlobalObservationStream();

  GM_Parms.setStride(globalObservationMatrix.stride());

  // Utilize both the partition information and elimination order
  // information already computed and contained in the file. This
  // enables the program to use external triangulation programs,
  // where this program ensures that the result is triangulated
  // and where it reports the quality of the triangulation.
  
  string tri_file;
  if (triFileName == NULL) 
    tri_file = string(strFileName) + GMTemplate::fileExtension;
  else 
    tri_file = string(triFileName);
  GMTemplate gm_template(fp);
  {
    // do this in scope so that is gets deleted now rather than later.
    iDataStreamFile is(tri_file.c_str());
    if (!fp.readAndVerifyGMId(is,checkTriFileCards))
      error("ERROR: triangulation file '%s' does not match graph given in structure file '%s'\n",tri_file.c_str(),strFileName);
    
    gm_template.readPartitions(is);
    gm_template.readMaxCliques(is);
  }
  gm_template.triangulatePartitionsByCliqueCompletion();
  if (1) { 
    // check that graph is indeed triangulated.
    // TODO: perhaps take this check out so that inference code does
    // not need to link to the triangulation code (either that, or put
    // the triangulation check in a different file, so that we only
    // link to tri check code).
    BoundaryTriangulate triangulator(fp,
				     gm_template.maxNumChunksInBoundary(),
				     gm_template.chunkSkip(),1.0);
    triangulator.ensurePartitionsAreChordal(gm_template);
  }


  ////////////////////////////////////////////////////////////////////
  // CREATE JUNCTION TREE DATA STRUCTURES
  infoMsg(IM::Default,"Creating Junction Tree\n"); fflush(stdout);
  JunctionTree myjt(gm_template);
  myjt.setUpDataStructures(varPartitionAssignmentPrior,varCliqueAssignmentPrior);
  myjt.prepareForUnrolling();
  if (jtFileName != NULL)
    myjt.printAllJTInfo(jtFileName);
  infoMsg(IM::Default,"DONE creating Junction Tree\n"); fflush(stdout);
  ////////////////////////////////////////////////////////////////////


  if (globalObservationMatrix.numSegments()==0) {
    infoMsg(IM::Default,"ERROR: no segments are available in observation file. Exiting...");
    exit_program_with_status(0);
  }

  Range* dcdrng = new Range(dcdrng_str,0,globalObservationMatrix.numSegments());
#if 0
  if (dcdrng->length() <= 0) {
    infoMsg(IM::Default,"Training range '%s' specifies empty set. Exiting...\n",
	  dcdrng_str);
    exit_program_with_status(0);
  }
#endif

  if (dcdrng->length() == 0) { 
    error("Decoding range must specify a non-zero length range. Range given is %s\n",
	  dcdrng_str);
  }

  logpr total_data_prob = 1.0;

  FILE* vitValsFile;
  if (showVitVals) {
    if (strcmp(showVitVals,"-") == 0)
      vitValsFile = stdout;
    else {
      if ((vitValsFile = fopen(showVitVals, "w")) == NULL)
	error("Can't open file '%s' for writing\n",vitValsFile);
    }
  } else {
    vitValsFile = NULL;
  }

  map<int, string> word_map;
  if (wordVar != NULL)
  {
    if (varMapFile==NULL)
      error("File to map from word index to string not specified.");
    if(transitionLabel == NULL)
      error("The label of the transition variable was not specified.");
    ifstream in(varMapFile);
    if (!in) { cout << "Unable to open " << varMapFile << endl; exit(1); }
    string name;
    int val;
    while (!in.eof())
    {
      in >> val >> name >> ws;
      if ( in.eof() )
	break;
      // printf("reading %d word = (%s)\n",val,name.c_str());
      word_map[val] = name;
    }
    in.close();
  }

  set<string> dumpVars;
  map<string,int> posFor;
  int ndv=0;
  if (dumpNames)
  {
      ifstream din(dumpNames);
      if (!din) {cout << "Unable to open " << dumpNames << endl; exit(1);}
      string s;
      while (!din.eof())
      {
          din >> s >> ws;
          dumpVars.insert(s);
          posFor[s] = ndv++;  // which position within a frame to put it
      }
      din.close();
  }

  vector<string> ofiles;
  if (ofilelist)
  {
      ifstream oin(ofilelist);
      if (!oin) {cout << "Unable to open " << ofilelist << endl; exit(1);} 
      string s;
      while (!oin.eof())
      {
          oin >> s >> ws;
          ofiles.push_back(s);
      }
      oin.close();
  }
  // We always do viterbi scoring/option in this program.
  JunctionTree::viterbiScore = true;

  struct rusage rus; /* starting time */
  struct rusage rue; /* ending time */
  getrusage(RUSAGE_SELF,&rus);

  total_data_prob = 1.0;
  Range::iterator* dcdrng_it = new Range::iterator(dcdrng->begin());

  if (island) {
    if (MixtureCommon::cacheMixtureProbabilities == true) {
      infoMsg(IM::Default,"NOTE: with island algorithm, might want to also try turning off Gaussian component caching with '-componentCache F'\n"); 
      fflush(stdout);
    }
  }
    
  while (!dcdrng_it->at_end()) {
    const unsigned segment = (unsigned)(*(*dcdrng_it));
    if (globalObservationMatrix.numSegments() < (segment+1)) 
      error("ERROR: only %d segments in file, decode range must be in range [%d,%d] inclusive\n",
	    globalObservationMatrix.numSegments(),
	    0,globalObservationMatrix.numSegments()-1);

    const unsigned numFrames = GM_Parms.setSegment(segment);

    logpr probe;
    if (island) {
      // warning("WARNING:: Island algorithm for decoding not yet debugged (but almost). Use at your own risk!!\n");
      unsigned numUsableFrames;
      myjt.collectDistributeIsland(numFrames,
				   numUsableFrames,
				   base,
				   lst,
				   false, // run EM algorithm
				   true,  // run viterbi algorithm
				   false  // localCliqueNormalization, unused here.
				   );
      probe = myjt.curProbEvidenceIsland();
      printf("Segment %d, after Island, viterbi log(prob(evidence)) = %f, per frame =%f, per numUFrams = %f\n",
	     segment,
	     probe.val(),
	     probe.val()/numFrames,
	     probe.val()/numUsableFrames);
      if (probe.not_essentially_zero()) {
	total_data_prob *= probe;
      }
    } else {
      // linear space inference
      unsigned numUsableFrames = myjt.unroll(numFrames);
      infoMsg(IM::Med,"Collecting Evidence\n");
      myjt.collectEvidence();
      infoMsg(IM::Med,"Done Collecting Evidence\n");
      probe = myjt.probEvidence();
      infoMsg(IM::Default,"Segment %d, after CE, viterbi log(prob(evidence)) = %f, per frame =%f, per numUFrams = %f\n",
	     segment,
	     probe.val(),
	     probe.val()/numFrames,
	     probe.val()/numUsableFrames);
      if (probe.essentially_zero()) {
	infoMsg(IM::Default,"Skipping segment since probability is essentially zero\n");
      } else {
	myjt.setRootToMaxCliqueValue();
	total_data_prob *= probe;
	infoMsg(IM::Low,"Distributing Evidence\n");
	myjt.distributeEvidence();
	infoMsg(IM::Low,"Done Distributing Evidence\n");
      }
    }
    if (vitValsFile) {
      if (probe.essentially_zero())
	printf("Not printing Viterbi values since segment has zero probability\n");
      else 
	myjt.printCurrentRVValues(vitValsFile);
    }

    if (dumpNames)
    {
        if (segment >= ofiles.size()) 
	  error("More utterances than output files");
        FILE *fp = fopen(ofiles[segment].c_str(), "wb");
        if (!fp) {cout << "Unable to open " << ofiles[segment] << endl; exit(1);}
        int *vals = new int[myjt.curNodes().size()];
        // just in case a variable isn't present in a slice, write a -1
        for (int i=0; i<int(myjt.curNodes().size()); i++) vals[i] = -1;
        int nv=0;
        for (int i=0; i<int(myjt.curNodes().size()); i++)
            if (dumpVars.count(myjt.curNodes()[i]->name())) 
            {
                if (!myjt.curNodes()[i]->discrete()) 
                    error("variables to dump must be discrete and hidden");
                int p = myjt.curNodes()[i]->frame()*dumpVars.size() 
                        + posFor[myjt.curNodes()[i]->name()];
                assert(p<int(myjt.curNodes().size()));
                vals[p] = RV2DRV(myjt.curNodes()[i])->val; 
                nv++;
            }
        fwrite((void *) vals, sizeof(int), nv, fp); 
        delete [] vals;
        fclose(fp);
    }
    if (wordVar && probe.not_essentially_zero())
    {
      // print the sequence of values for this variable
      // compress consecutive values into a single instance
      // the times are right if a word transition at time t means there is
      // a new word at t+1
      string pvn = string(wordVar);
      string tl = string(transitionLabel);
      for (int i=0, lv=-1, lf=0; i<int(myjt.curNodes().size()); i++)
      {
        if (myjt.curNodes()[i]->name() == pvn)
        {
          if (!myjt.curNodes()[i]->discrete()) 
            error("Can only print Viterbi values for discrete variables");
          if (RV2DRV(myjt.curNodes()[i])->cardinality != word_map.size())
            error("Word-val to string map file size %d does not match the number of words %d.",int(word_map.size()),RV2DRV(myjt.curNodes()[i])->cardinality);
          lv = RV2DRV(myjt.curNodes()[i])->val;
        }
        else if (myjt.curNodes()[i]->name()==tl)
        {
          if (!myjt.curNodes()[i]->discrete()) 
            error("Word transition variable should be discrete");
          if (RV2DRV(myjt.curNodes()[i])->cardinality != 2) 
            error("Word transition variable should have two values");
          if (RV2DRV(myjt.curNodes()[i])->val==1)  // a word transition
          {
            cout << word_map[lv] << " (" << lf << "-" 
                 << myjt.curNodes()[i]->frame() << ")\n" << flush;  
            lf = myjt.curNodes()[i]->frame()+1;
          }
        }
      }
    }
    (*dcdrng_it)++;
  }

  infoMsg(IM::Default,"Total data log prob is: %1.9e\n",
	  total_data_prob.val());

  if (vitValsFile && vitValsFile != stdout)
    fclose(vitValsFile);

  getrusage(RUSAGE_SELF,&rue);
  if (IM::messageGlb(IM::Default)) { 
    infoMsg(IM::Default,"### Final time (seconds) just for decoding stage: ");
    double userTime,sysTime;
    reportTiming(rus,rue,userTime,sysTime,stdout);
  }

  exit_program_with_status(0);
}
