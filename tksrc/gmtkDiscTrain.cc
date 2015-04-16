/*
 * gmtkKernel.cc
 *  Use a DGM as a Kernel (e.g., Fisher Kernel, or something else).
 *  In other words, for each observation file, we will write out 
 *  a vector that is on the order of the number of current system parameters.
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (C) 2001 Jeff Bilmes
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

#include <algorithm>
#include <cstdlib>
#include <ctime>

#include "general.h"
#include "error.h"
#include "debug.h"
#include "rand.h"
#include "arguments.h"
#include "ieeeFPsetup.h"

#include "GMTK_WordOrganization.h"

VCID(HGID)


#include "GMTK_FileParser.h"
#include "GMTK_RV.h"
#include "GMTK_DiscRV.h"
#include "GMTK_ContRV.h"
#include "GMTK_GMTemplate.h"
#include "GMTK_GMParms.h"
#if 0
#  include "GMTK_ObservationMatrix.h"
#else
#  include "GMTK_ObservationSource.h"
#  include "GMTK_FileSource.h"
#  include "GMTK_CreateFileSource.h"
#  include "GMTK_ASCIIFile.h"
#  include "GMTK_FlatASCIIFile.h"
#  include "GMTK_PFileFile.h"
#  include "GMTK_HTKFile.h"
#  include "GMTK_HDF5File.h"
#  include "GMTK_BinaryFile.h"
#  include "GMTK_Filter.h"
#  include "GMTK_Stream.h"
#endif
#include "GMTK_MixtureCommon.h"
#include "GMTK_GaussianComponent.h"
#include "GMTK_MeanVector.h"
#include "GMTK_DiagCovarVector.h"
#include "GMTK_DlinkMatrix.h"
#include "GMTK_ProgramDefaultParms.h"
#include "GMTK_BoundaryTriangulate.h"
#include "GMTK_JunctionTree.h"
#include "GMTK_MaxClique.h"

#include "GMTK_DiagGaussian.h"
#include "GMTK_MDCPT.h"
#include "GMTK_Dense1DPMF.h"


/*****************************   DISCRIMINITIVE TRAINING   **********************************************/

#define GMTK_ARG_MMI_TRAINING

/*****************************   OBSERVATION INPUT FILE HANDLING   **********************************************/
#define GMTK_ARG_OBS_FILES

/*************************   INPUT TRAINABLE PARAMETER FILE HANDLING  *******************************************/
#define GMTK_ARG_INPUT_TRAINABLE_FILE_HANDLING
#define GMTK_ARG_CPP_CMD_OPTS
#define GMTK_ARG_INPUT_MASTER_FILE
#define GMTK_ARG_OUTPUT_MASTER_FILE
#define GMTK_ARG_DLOPEN_MAPPERS
#define GMTK_ARG_INPUT_TRAINABLE_PARAMS
#define GMTK_ARG_OUTPUT_TRAINABLE_PARAMS
#define GMTK_ARG_WPAEEI
#define GMTK_ARG_ALLOC_DENSE_CPTS
#define GMTK_ARG_CPT_NORM_THRES

/*************************   INPUT STRUCTURE PARAMETER FILE HANDLING  *******************************************/
#define GMTK_ARG_INPUT_MODEL_FILE_HANDLING
#define GMTK_ARG_STR_FILE
#define GMTK_ARG_TRI_FILE
#define GMTK_ARG_CHECK_TRI_FILE_CARD
#define GMTK_ARG_JT_INFO_FILE
#define GMTK_ARG_JTW_UB
#define GMTK_ARG_LATTICE_PARAMS

/*************************   CONTINUOUS RANDOM VARIABLE OPTIONS       *******************************************/
#define GMTK_ARG_CONTINUOUS_RANDOM_VAR_OPTIONS
#define GMTK_ARG_VAR_FLOOR
#define GMTK_ARG_VAR_FLOOR_ON_READ


/*************************          BEAM PRUNING OPTIONS              *******************************************/
#define GMTK_ARG_BEAM_PRUNING_OPTIONS
#define GMTK_ARG_CBEAM
#define GMTK_ARG_CPBEAM
#define GMTK_ARG_CKBEAM
#define GMTK_ARG_CCBEAM
#define GMTK_ARG_CRBEAM
#define GMTK_ARG_CMBEAM
#define GMTK_ARG_SBEAM
#define GMTK_ARG_EBEAM

/*************************          MEMORY MANAGEMENT OPTIONS         *******************************************/
#define GMTK_ARG_MEMORY_MANAGEMENT_OPTIONS
#define GMTK_ARG_HASH_LOAD_FACTOR
#define GMTK_ARG_STORE_DETERMINISTIC_CHILDREN
#define GMTK_ARG_CLEAR_CLIQUE_VAL_MEM
#define GMTK_ARG_MEM_GROWTH
#define GMTK_ARG_USE_MMAP

/****************************      FILE RANGE OPTIONS             ***********************************************/
#define GMTK_ARG_DCDRNG
#define GMTK_ARG_START_END_SKIP

/****************************         GENERAL OPTIONS             ***********************************************/
#define GMTK_ARG_GENERAL_OPTIONS
#define GMTK_ARG_FILE_RANGE_OPTIONS
#define GMTK_ARG_SEED
#define GMTK_ARG_SKIP_STARTUP_CHECKS
#define GMTK_ARG_VERB
#define GMTK_ARG_HELP
#define GMTK_ARG_VERSION

/****************************         INFERENCE OPTIONS           ***********************************************/
#define GMTK_ARG_INFERENCE_OPTIONS
#define GMTK_ARG_ONLY_KEEP_SEPS
#define GMTK_ARG_ISLAND
#define GMTK_ARG_CLIQUE_TABLE_NORMALIZE
#define GMTK_ARG_CE_SEP_DRIVEN
#define GMTK_ARG_COMPONENT_CACHE
#define GMTK_ARG_CLIQUE_VAR_ITER_ORDERS
#define GMTK_ARG_JT_OPTIONS
#define GMTK_ARG_VE_SEPS
#define GMTK_ARG_FAIL_ON_ZERO_CLIQUE

/****************************         EM TRAINING OPTIONS         ***********************************************/
#define GMTK_ARG_KERNEL_OPTIONS
#define GMTK_ARG_KERNEL_PARAMS

/************************  OBSERVATION MATRIX TRANSFORMATION OPTIONS   ******************************************/
#define GMTK_ARG_OBS_MATRIX_OPTIONS
#define GMTK_ARG_OBS_MATRIX_XFORMATION


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


typedef vector<double> vector_d;


/*
 *
 * definition of needed global arguments
 *
 */
RAND rnd(seedme);
GMParms GM_Parms;
GMParms GM_Parms_n;
GMParms GM_Parms_d;
#if 0
ObservationMatrix globalObservationMatrix;
#endif

FileSource *gomFS;
ObservationSource *globalObservationMatrix;


void setUpGM_Parms(FileParser& fp, GMTemplate& gm_template, GMParms & GM_Parms, bool prime) {

    /////////////////////////////////////////////
    // read in all the parameters
    dlopenDeterministicMaps(dlopenFilenames, MAX_NUM_DLOPENED_FILES);
    if (prime && inputMasterFile) {
        // flat, where everything is contained in one file, always ASCII
        iDataStreamFile pf(inputMasterFile,false,true,cppCommandOptions);
        GM_Parms.read(pf);
    }
    else if(!prime && denomInputMasterFile) {
        iDataStreamFile pf(denomInputMasterFile,false,true,cppCommandOptions);
        GM_Parms.read(pf);
    }

    if (prime && inputTrainableParameters) {
        // flat, where everything is contained in one file
        iDataStreamFile pf(inputTrainableParameters, binInputTrainableParameters,true,cppCommandOptions);
        GM_Parms.readTrainable(pf);
    }
    else if(!prime && denomInputTrainableParameters) {
        iDataStreamFile pf(denomInputTrainableParameters, denomBinInputTrainableParameters,true,cppCommandOptions);
        GM_Parms.readTrainable(pf);
    }

  // comment for now Sun Jan 11 09:47:23 2004
  GM_Parms.finalizeParameters();
  GM_Parms.markObjectsToNotTrain(objsToNotUtilizeFile,cppCommandOptions);

  /////////////////////////////
  // read in the structure of the GM, this will
  // die if the file does not exist.
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

  GM_Parms.setStride(gomFS->stride());

  // Utilize both the partition information and elimination order
  // information already computed and contained in the file. This
  // enables the program to use external triangulation programs,
  // where this program ensures that the result is triangulated
  // and where it reports the quality of the triangulation.
  
    string tri_file;
    if (prime) {
        if (triFileName == NULL) 
            tri_file = string(strFileName) + GMTemplate::fileExtension;
        else 
            tri_file = string(triFileName);
    }
    else {
        if (denomTriFileName == NULL) 
            tri_file = string(denomStrFileName) + GMTemplate::fileExtension;
        else 
            tri_file = string(denomTriFileName);
    }
  
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

  //  printf("Dlinks: min lag %d    max lag %d\n", Dlinks::globalMinLag(), Dlinks::globalMaxLag());
  // FIXME - min past = min(dlinkPast, VECPTPast), likewise for future
  int dlinkPast = Dlinks::globalMinLag();
  dlinkPast = (dlinkPast < 0) ? -dlinkPast : 0;
  gomFS->setMinPastFrames( dlinkPast );
  
  int dlinkFuture = Dlinks::globalMaxLag();
  dlinkFuture = (dlinkFuture > 0) ? dlinkFuture : 0;
  gomFS->setMinFutureFrames( dlinkFuture );

}



  ////////////////////////////////////////////////////////////////////
  // CREATE JUNCTION TREE DATA STRUCTURES
void createJunctionTree(JunctionTree & myjt, GMParms & GM_Parms) {
  infoMsg(IM::Default,"Creating Junction Tree\n"); fflush(stdout);
  myjt.setUpDataStructures(varPartitionAssignmentPrior,varCliqueAssignmentPrior);
  myjt.prepareForUnrolling();
  if (jtFileName != NULL)
    myjt.printAllJTInfo(jtFileName);
  infoMsg(IM::Default,"DONE creating Junction Tree\n"); fflush(stdout);
  ////////////////////////////////////////////////////////////////////

  if (randomizeParams) {
    infoMsg(IM::Default,"WARNING: GMTK is randomizing all trainable parameters and writing them to file 'random.gmp'\n");
    GM_Parms.makeRandom();
    oDataStreamFile of("random.gmp");
    GM_Parms.writeTrainable(of);
  }

  if (gomFS->numSegments()==0) {
    infoMsg(IM::Default,"ERROR: no segments are available in observation file. Exiting...");
    exit_program_with_status(0);
  }
}


void writeTrainedOutput(string & index) {

    string fname(storeFeatureFile);
    fname += "." + index;
    oDataStreamFile outf(fname.c_str(), transFileIsBinary);

    outf.nl();

    GM_Parms.writeTrainable(outf, false);

    outf.nl();
    //outf.close();
}


void train(JunctionTree & myjt, const char* range_str, vector<double> & data_probs, bool initEM) {


    Range* dcdrng = new Range(range_str,0,gomFS->numSegments());


  if (dcdrng->length() <= 0) {
    infoMsg(IM::Default,"Training range '%s' specifies empty set. Exiting...\n",
	  dcdrng_str);
    exit_program_with_status(0);
  }

  struct rusage rus; /* starting time */
  struct rusage rue; /* ending time */
  getrusage(RUSAGE_SELF,&rus);

  // Now, do CE/DE iterations, writing out the parameter "increments" one per segment.
  Range::iterator* dcdrng_it = new Range::iterator(dcdrng->begin());
  {
    oDataStreamFile outf(storeFeatureFile,transFileIsBinary);
    bool firstTime = true;
    while (!dcdrng_it->at_end()) {
      const unsigned segment = (unsigned)(*(*dcdrng_it));
      try {
	if (gomFS->numSegments() < (segment+1)) 
	  error("ERROR: only %d segments in file, segment must be in range [%d,%d]\n",
		gomFS->numSegments(),
		0,gomFS->numSegments()-1);

	const unsigned numFrames = GM_Parms.setSegment(segment);

	logpr data_prob = 1.0;
    if(initEM) GM_Parms.emInitAccumulators(firstTime);

	unsigned numUsableFrames;
	if (island) {
	  myjt.collectDistributeIsland(numFrames,
				       numUsableFrames,
				       base,
				       lst,
				       rootBase, islandRootPower, 
				       true, // run EM algorithm,
				       false, // run Viterbi algorithm
				       localCliqueNormalization);
	  printf("Segment %d, after Island, log(prob(evidence)) = %f, per frame =%f, per numUFrams = %f, ",
		 segment,
		 myjt.curProbEvidenceIsland().val(),
		 myjt.curProbEvidenceIsland().val()/numFrames,
		 myjt.curProbEvidenceIsland().val()/numUsableFrames);
	  data_prob = myjt.curProbEvidenceIsland();
	} else if (onlyKeepSeparators) {

	  infoMsg(IM::Low,"Collecting Evidence (linear space)\n");
	  data_prob = myjt.collectEvidenceOnlyKeepSeps(numFrames, &numUsableFrames);
	  infoMsg(IM::Low,"Done Collecting Evidence\n");

	  infoMsg(IM::Low,"Distributing Evidence\n");
	  myjt.distributeEvidenceOnlyKeepSeps();
	  infoMsg(IM::Low,"Done Distributing Evidence\n");

	  printf("Segment %d, after CE/DE, log(prob(evidence)) = %f, per frame =%f, per numUFrams = %f, ",
		 segment,
		 data_prob.val(),
		 data_prob.val()/numFrames,
		 data_prob.val()/numUsableFrames);
      
	  myjt.emIncrement(data_prob,localCliqueNormalization);
	} else {

	  numUsableFrames = myjt.unroll(numFrames);
	  gomFS->justifySegment(numUsableFrames);
	  infoMsg(IM::Low,"Collecting Evidence\n");
	  myjt.collectEvidence();
	  infoMsg(IM::Low,"Done Collecting Evidence\n");
	  data_prob = myjt.probEvidence();

	  infoMsg(IM::Low,"Distributing Evidence\n");
	  myjt.distributeEvidence();
	  infoMsg(IM::Low,"Done Distributing Evidence\n");

	  printf("Segment %d, after CE/DE, log(prob(evidence)) = %f, per frame =%f, per numUFrams = %f, ",
		 segment,
		 data_prob.val(),
		 data_prob.val()/numFrames,
		 data_prob.val()/numUsableFrames);
      
	  myjt.emIncrement(data_prob,localCliqueNormalization);

	}

    data_probs.push_back(data_prob.val());

	printf("writing %s-kernel feature space vector ...\n",(fisherKernelP?"Fisher":"accumulator"));
	if (annotateTransformationOutput) {
	  char buff[1024];
	  sprintf(buff,"Segment %d : %d frames, %d usable frames, log(PE) = %f",segment,numFrames,numUsableFrames,data_prob.val());
	  outf.write(buff);
	  outf.nl();
	};
	outf.writeComment("segment %d  log(PE) %f\n", segment, data_prob.val());
	outf.write(data_prob.val());
	outf.nl();
	//GM_Parms.emWriteUnencodedAccumulators(outf,writeLogVals);
	GM_Parms.writeTrainable(outf, false);
	outf.nl();

      } catch (ZeroCliqueException &e) {
	warning("Segment %d aborted due to zero clique\n", segment);
      }
      
      (*dcdrng_it)++;
      firstTime = false;

        
    }
  }

  getrusage(RUSAGE_SELF,&rue);
  if (IM::messageGlb(IM::Default)) { 
    infoMsg(IM::Default,"### Final time (seconds) just for kernel computing stage: ");
    double userTime,sysTime;
    reportTiming(rus,rue,userTime,sysTime,stdout);
  }

}


void initAdaGrad(map<string, vector_d>& adagrad) {
    for(unsigned i=0; i<GM_Parms.components.size(); i++) {
        if(GM_Parms_n.components[i]->typeName() == "Diag Gaussian") {
            DiagGaussian* dg = (DiagGaussian*)GM_Parms_n.components[i];
            string mean_name = dg->getMean()->name();
            string covar_name = dg->getCovar()->name();

            sArray<float>& means = dg -> getMean() -> getMeans();
            sArray<float>& covars = dg -> getCovar() ->getCovars();

            for(unsigned j=0; j<means.size(); j++) adagrad[mean_name].push_back(1.0);
            for(unsigned j=0; j<covars.size(); j++) adagrad[covar_name].push_back(1.0);

        }

    }
    
    //GM_Parms.mdCpts.size() should be 3 for TIMIT
    for(unsigned i=0; i<GM_Parms.mdCpts.size(); i++) {
        MDCPT* mdcpt = (MDCPT*)GM_Parms.mdCpts[i];

        string name = mdcpt->name();
        if(name == "internal:UnityScore") continue;

        sArray<logpr>& mdcpts = mdcpt->getMdcpt();
        for(unsigned j=0; j<mdcpts.size(); j++) adagrad[name].push_back(1.0);
        
    }

    for(unsigned i=0; i<GM_Parms.dPmfs.size(); ++ i) {
	Dense1DPMF* dpmf = (Dense1DPMF*)GM_Parms.dPmfs[i];

        string name = dpmf->name();

	sArray<logpr>& pmfs = dpmf->getPmf();
	for(unsigned j=0; j<pmfs.size(); ++ j) adagrad[name].push_back(1.0);
    }
}




void updateAdagrad(map<string, vector_d>& adagrad) {
    for(unsigned i=0; i<GM_Parms.components.size(); i++) {
        if(GM_Parms_n.components[i]->typeName() == "Diag Gaussian") {
            DiagGaussian* dg = (DiagGaussian*)GM_Parms_n.components[i];
            string mean_name = dg->getMean()->name();
            string covar_name = dg->getCovar()->name();

            sArray<float>& next_means = dg -> getNextMeans();
            sArray<float>& next_covars = dg -> getNextCovars();

            for(unsigned j=0; j<next_means.size(); j++) {
                double old_val = adagrad[mean_name][j];
                adagrad[mean_name][j] = sqrt(old_val * old_val + next_means[j] * next_means[j]);
            }

            for(unsigned j=0; j<next_covars.size(); j++) {
                double old_val = adagrad[covar_name][j];
                adagrad[covar_name][j] = sqrt(old_val * old_val + next_covars[j] * next_covars[j]);
            
            }
        }
    }
    
    //GM_Parms.mdCpts.size() should be 3 for TIMIT
    for(unsigned i=0; i<GM_Parms.mdCpts.size(); i++) {
        MDCPT* mdcpt = (MDCPT*)GM_Parms.mdCpts[i];

        string name = mdcpt->name();
        if(name == "internal:UnityScore") continue;

        sArray<logpr>& next_mdcpts = mdcpt->getNextMdcpt();

        for(unsigned j=0; j<next_mdcpts.size(); j++) {
            double old_val = adagrad[name][j];
            double acc_val = next_mdcpts[j].unlog();
            adagrad[name][j] = sqrt(old_val * old_val + acc_val * acc_val);
        }
    }
}



void normalizeCPT(sArray<logpr>& mdcpts) {
    bool has_parent = false;
    
    double sum = 0.0;
    if(!has_parent) {
        for(unsigned i=0; i<mdcpts.size(); i++) {
            sum += mdcpts[i].unlog();
        }
        if(sum == 0.0) return;
    
        logpr log_sum(sum);
        for(unsigned i=0; i<mdcpts.size(); i++) {
            mdcpts[i] /= log_sum;
        }
    }


}





void updateBoth(double lr, double mean_lr, double covar_lr, map<string, vector_d>& adagrad, bool use_adagrad) {
    

    double local_lr = lr;

    for(unsigned i=0; i<GM_Parms_n.components.size(); i++) {
        if(GM_Parms_n.components[i]->typeName() == "Diag Gaussian") {
            DiagGaussian* dg = (DiagGaussian*)GM_Parms_n.components[i];
            string mean_name = dg->getMean()->name();
            string covar_name = dg->getCovar()->name();


            DiagGaussian* dg_d = NULL;
            bool found_match = false;


            //first check if the ith component for denominator and numerator are the same, for the reason that they are sometimes the same model
            if(GM_Parms_d.components[i]->typeName().compare("Diag Gaussian") == 0) {
                dg_d = (DiagGaussian*)GM_Parms_d.components[i];
                if(dg_d->getMean()->name().compare(mean_name) == 0 && dg_d->getCovar()->name().compare(covar_name) == 0) {
                    found_match = true;
                }
            }
            if(!found_match) {
                for(unsigned j=0; j<GM_Parms_d.components.size(); ++ j) {
                    if(GM_Parms_d.components[j]->typeName().compare("Diag Gaussian") != 0) continue;
                    else {
                        dg_d = (DiagGaussian*)GM_Parms_d.components[j];
                        if(dg_d->getMean()->name().compare(mean_name) == 0 && dg_d->getCovar()->name().compare(covar_name) == 0) {
                            found_match = true;
                            break;
                        }
                    }
                }
            }

            if(!found_match) {
                fprintf(stderr, "WARNING: mean/covar name different for numerator not found in denominator %u : %s, %s\n", i, mean_name.c_str(), covar_name.c_str());
                continue;
            }
            


            sArray<float>& means = dg -> getMean() -> getMeans();
            sArray<float>& covars = dg -> getCovar() ->getCovars();

            sArray<float>& means_d = dg_d -> getMean() -> getMeans();
            sArray<float>& covars_d = dg_d -> getCovar() ->getCovars();

            sArray<float>& next_means = dg -> getNextMeans();
            sArray<float>& next_covars = dg -> getNextCovars();

            sArray<float>& next_means_d = dg_d->getNextMeans();
            sArray<float>& next_covars_d = dg_d->getNextCovars();


            double acc_val;

	        if(update_mean && mean_name.compare("intensity_mean") != 0) {
	            for(unsigned j=0; j<means.size(); j++) {
		            acc_val = next_means[j] - down_weight_mean * next_means_d[j];
		            if(use_adagrad) {
		                double old_val = adagrad[mean_name][j];
			            local_lr = mean_lr / old_val;
			            adagrad[mean_name][j] = sqrt(old_val * old_val + acc_val * acc_val);
		            }
		            else local_lr = mean_lr;

		            means[j] += local_lr * acc_val;
		            means_d[j] = means[j];
		
		        }
	        }


	        //if(update_covar && covar_name.compare("covar0") != 0) {
            if(update_covar) {

		    //double down_weight = 1.0;
		    //if(covar_name.compare("covar0") != 0) down_weight = 1.2;
		    //else down_weight = 1.0;

	            for(unsigned j=0; j<covars.size(); j++) {
		            acc_val = (next_covars[j] - down_weight_covar * next_covars_d[j]); //* covar_lr;
		            //if(i == 13) fprintf(stderr, "s: %s i: %d n: %f %f, d: %f %f, acc: %f\n", covar_name.c_str(), j, covars[j], next_covars[j], covars_d[j], next_covars_d[j], acc_val);

		            if(use_adagrad) {
		                double old_val = adagrad[covar_name][j];
			            local_lr = covar_lr / old_val;
			            adagrad[covar_name][j] = sqrt(old_val * old_val + acc_val * acc_val);
		            }
		            else local_lr = covar_lr;

		            covars[j] += local_lr * acc_val;
		            if(covars[j] < varFloor) covars[j] = varFloor;
		            covars_d[j] = covars[j];
                    //if(i == 13) fprintf(stderr, "s2: %s i: %d n: %f %f, d: %f %f, acc2: %f %f %f\n", covar_name.c_str(), j, covars[j], next_covars[j], covars_d[j], next_covars_d[j], local_lr * acc_val, local_lr, lr);
		        }
	        }

            for(unsigned j=0; j<next_means.size(); ++ j) {
                next_means[j] = 0.0;
                next_means_d[j] = 0.0;
            }

            for(unsigned j=0; j<next_covars.size(); ++ j) {
                next_covars[j] = 0.0;
                next_covars_d[j] = 0.0;
            }

            dg -> getCovar()->preCompute();
            dg_d -> getCovar()->preCompute();


        }
    }

    

    for(unsigned i=0; i<GM_Parms_n.mdCpts.size(); ++ i) {
        MDCPT* mdcpt = (MDCPT*)GM_Parms_n.mdCpts[i];

        string name = mdcpt->name();
        if(name == "internal:UnityScore") continue;


        MDCPT* mdcpt_d = (MDCPT*)GM_Parms_d.mdCpts[i];
        bool found_match = false;

        
        if(mdcpt_d -> name().compare(name) == 0) {
            found_match = true;
        }

        if(!found_match) {
            for(unsigned j=0; j<GM_Parms_d.mdCpts.size(); ++ j) {
                mdcpt_d = (MDCPT*)GM_Parms_d.mdCpts[j];
                if(mdcpt_d -> name().compare(name) == 0) {
                    found_match = true;
                    break;
                }
            }
        }

        if(!found_match) {
            fprintf(stderr, "WARNING: mdcpt name in numerator not found in denominator model %u : %s\n", i, name.c_str());
            continue;
        }


        sArray<logpr>& mdcpts = mdcpt->getMdcpt();
        sArray<logpr>& next_mdcpts = mdcpt->getNextMdcpt();

        sArray<logpr>& mdcpts_d = mdcpt_d -> getMdcpt();
        sArray<logpr>& next_mdcpts_d = mdcpt_d->getNextMdcpt();


        if(update_CPT) {
            //TODO: need normalization to simplex
            for(unsigned j=0; j<mdcpts.size(); j++) {
                double acc_val = next_mdcpts[j].unlog() - next_mdcpts_d[j].unlog();
                if(use_adagrad) {
                    double old_val = adagrad[name][j];
                    local_lr = lr / old_val;
                    adagrad[name][j] = sqrt(old_val * old_val + acc_val * acc_val);
                }
                else local_lr = lr;

                double p = mdcpts[j].unlog() + local_lr * acc_val;
                if(p < 0.0) p = 0.0;
                mdcpts[j] = logpr(p);
            }
        }//end of updateCPT

        //normalizeCPT(mdcpts);
        //fprintf(stdout, "====mdcpt====\n");
        for(unsigned j=0; j<mdcpts.size(); j++){
            mdcpts_d[j] = mdcpts[j];
        //    fprintf(stdout, "%f\n", mdcpts[j].val());
        }
        //fprintf(stdout, "\n");


        for(unsigned j=0; j<next_mdcpts.size(); j++) {
            next_mdcpts[j] = 0.0;
            next_mdcpts_d[j] = 0.0;
        }

    }

    
    //DPMFs
    if(update_DPMF) {
        for(unsigned i=0; i<GM_Parms_n.dPmfs.size(); ++ i) {
            Dense1DPMF* dpmf = (Dense1DPMF*)GM_Parms_n.dPmfs[i];

            string name = dpmf->name();

            Dense1DPMF* dpmf_d = (Dense1DPMF*)GM_Parms_d.dPmfs[i];
            bool found_match = false;

            if(dpmf_d->name().compare(name) == 0) {
                found_match = true;
            }
            if(!found_match) {
                for(unsigned j=0; j<GM_Parms_d.dPmfs.size(); ++ j) {
                    dpmf_d = (Dense1DPMF*)GM_Parms_d.dPmfs[j];
                    if(dpmf_d->name().compare(name) == 0) {
                        found_match = true;
                        break;
                    }
                }
            }

            if(!found_match) {
                fprintf(stderr, "WARNING: dpmf name for numerator not found in denominator %u : %s\n", i, name.c_str());
                continue;
            }


            sArray<logpr>& pmfs = dpmf->getPmf();
            sArray<logpr>& next_pmfs = dpmf->getNextPmf();

            sArray<logpr>& pmfs_d = dpmf_d->getPmf();
            sArray<logpr>& next_pmfs_d = dpmf_d->getNextPmf();


            for(unsigned j=0; j<pmfs.size(); j++) {
                double acc_val = next_pmfs[j].unlog() - next_pmfs_d[j].unlog();
                if(use_adagrad) {
                    double old_val = adagrad[name][j];
                    local_lr = lr / old_val;
                    adagrad[name][j] = sqrt(old_val * old_val + acc_val * acc_val);
                }
                else local_lr = lr;

                double p = pmfs[j].unlog() + local_lr * acc_val;
                if(p < 0.0) p = 0.0;
                pmfs[j] = logpr(p);
            }

            //normalize??
            normalizeCPT(pmfs);
            
            //fprintf(stdout, "====dpmfs====\n");
            for(unsigned j=0; j<pmfs.size(); j++){
                pmfs_d[j] = pmfs[j];
            }
            //fprintf(stdout, "\n");


            for(unsigned j=0; j<next_pmfs.size(); j++) {
                next_pmfs[j] = 0.0;
                next_pmfs_d[j] = 0.0;
            }
        }
    }

}


void update(double lr, map<string, vector_d>& adagrad, bool use_adagrad) {

    if(use_adagrad) updateAdagrad(adagrad);

    for(unsigned i=0; i<GM_Parms.components.size(); i++) {
        if(GM_Parms_n.components[i]->typeName() == "Diag Gaussian") {
            DiagGaussian* dg = (DiagGaussian*)GM_Parms_n.components[i];
            string mean_name = dg->getMean()->name();
            string covar_name = dg->getCovar()->name();

            sArray<float>& means = dg -> getMean() -> getMeans();
            sArray<float>& covars = dg -> getCovar() ->getCovars();

            sArray<float>& next_means = dg -> getNextMeans();
            sArray<float>& next_covars = dg -> getNextCovars();



            for(unsigned j=0; j<means.size(); j++) {
                if(use_adagrad) means[j] += lr * next_means[j] / adagrad[mean_name][j]; 
                else means[j] += lr * next_means[j]; 
            }
            for(unsigned j=0; j<covars.size(); j++) {
                if(use_adagrad) covars[j] += lr * next_covars[j] / adagrad[covar_name][j]; 
                else covars[j] += lr * next_covars[j]; 
            }
        }

    }
    
    //GM_Parms.mdCpts.size() should be 3 for TIMIT
    for(unsigned i=0; i<GM_Parms.mdCpts.size(); i++) {
        MDCPT* mdcpt = (MDCPT*)GM_Parms.mdCpts[i];

        string name = mdcpt->name();
        if(name == "internal:UnityScore") continue;

        sArray<logpr>& mdcpts = mdcpt->getMdcpt();
        sArray<logpr>& next_mdcpts = mdcpt->getNextMdcpt();


        for(unsigned j=0; j<mdcpts.size(); j++) {
            if(use_adagrad)
                mdcpts[j] = logpr(mdcpts[j].unlog() + lr * next_mdcpts[j].unlog() / adagrad[name][j]);
            else
                mdcpts[j] = logpr(mdcpts[j].unlog() + lr * next_mdcpts[j].unlog());
        }
    }
}



void i2str(int i, string& s) {
    std::vector<int> cs;

    if(i == 0) {
        s = "0";
        return;
    }

    while(i > 0) {
        int n = i % 10;
        i = i / 10;
        cs.push_back(n);
    }

    if(cs.size() > 0)
        for(int j=cs.size() - 1; j>=0; j--) {
            s += char(cs[j] + '0');
        }

}



void getGmParmsValue() {

    for(unsigned i=0; i<GM_Parms.components.size(); i++) {
        if(GM_Parms.components[i]->typeName() == "Diag Gaussian") {
            DiagGaussian* dg = (DiagGaussian*)GM_Parms.components[i];
            string covar_name = dg->getCovar()->name();
            
            sArray<float>& covars = dg -> getCovar() ->getCovars();
            sArray<float>& next_covars = dg -> getNextCovars();
            
            for(unsigned j=0; j<covars.size(); j++) {
                if(i == 13) {
                    if(j < next_covars.size()) fprintf(stderr, "WA s: %s i: %d n: %f %f\n", covar_name.c_str(), j, covars[j], next_covars[j]);
                    else fprintf(stderr, "WA s: %s i: %d n: %f\n", covar_name.c_str(), j, covars[j]);
                
                }
            }
            
        }
    }

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

    CODE_TO_COMPUTE_ENDIAN;

    ////////////////////////////////////////////
    // parse arguments
    bool parse_was_ok = Arg::parse(argc,(char**)argv,
        "\nThis program uses a DGM as a Kernel (e.g., Fisher Kernel\n"
        "or accumulator). In other words, for each observation file,\n"
        "it will write out a vector that is on the order of the number\n"
        "of current system parameters.\n");
    if(!parse_was_ok) {
        Arg::usage(); exit(-1);
    }

    #define GMTK_ARGUMENTS_CHECK_ARGS
    #include "GMTK_Arguments.h"
    #undef GMTK_ARGUMENTS_CHECK_ARGS


    // Check if we want to do the Fisher kernel rather than the
    // accumulator kernel.
    EMable::fisherKernelMode = fisherKernelP;


    gomFS = instantiateFileSource();
    globalObservationMatrix = gomFS;

    // NOTICE: all 2's are for denominator models; n for numerator, d for denominator

    FileParser fp_n(strFileName, cppCommandOptions);
    GMTemplate gm_template_n(fp_n);
    setUpGM_Parms(fp_n, gm_template_n, GM_Parms, true);
    GM_Parms_n = GM_Parms;
    

    GM_Parms.clearParms();

    FileParser fp_d(denomStrFileName, cppCommandOptions);
    GMTemplate gm_template_d(fp_d);
    setUpGM_Parms(fp_d, gm_template_d, GM_Parms, false);
    GM_Parms_d = GM_Parms;



    GM_Parms = GM_Parms_n;
    JunctionTree myjt_n(gm_template_n);
    createJunctionTree(myjt_n, GM_Parms);
    GM_Parms_n = GM_Parms;


    GM_Parms = GM_Parms_d;
    JunctionTree myjt_d(gm_template_d);
    createJunctionTree(myjt_d, GM_Parms);
    GM_Parms_d = GM_Parms;


    double lr = init_lr;
    double covarLr = covar_lr;
    double meanLr = mean_lr;

    
    map<string, vector_d> adagrad;
    initAdaGrad(adagrad);


    unsigned orig_num_instance = gomFS->numSegments();
    unsigned num_instance = orig_num_instance;

    vector<unsigned> random_index;

    //unsigned batch_size = 20;
    num_instance = num_instance / batch_size;
    if(num_instance * batch_size < orig_num_instance) ++ num_instance;


    fprintf(stderr, "Parameter Settings: update_CPT:%d, update_covar:%d, use_adagrad:%d, max_iter:%d, batch_size:%u, init_iter=%d\nLearning rate Settings: use_decay_lr:%d, init_lr:%e, decay_lr_rate:%e, use_covar_decay_lr:%d, init_covar_lr:%e, decay_covar_lr_rate:%e, use_mean_decay_lr:%d, init_mean_lr:%e, decay_mean_lr_rate:%e\n", update_CPT, update_covar, use_adagrad, max_iter, batch_size, init_iter, use_decay_lr, init_lr, decay_lr_rate, use_covar_decay_lr, covar_lr, decay_covar_lr_rate, use_mean_decay_lr, mean_lr, decay_mean_lr_rate);


    for(unsigned i=0; i<num_instance; i++) {
        random_index.push_back(i);
    }
    std::srand(random_seed);
    if(shuffle_data) std::random_shuffle(random_index.begin(), random_index.end());


    double cond_prob_iterPrev = 0.0;
    double lldp = 0.0;

    for(unsigned iter=0; iter<max_iter; iter++) {
        fprintf(stderr, "\n========= ITER %u =========\n", iter);

        double cond_prob_iter = 0.0;

        for(unsigned i=0; i<num_instance; i++) {

            if(use_decay_lr) {
	            lr = init_lr / pow((iter + init_iter) * num_instance + i + 1.0, decay_lr_rate);
		        // fprintf(stderr, "\t### lr: %f ###\n", lr);
            }

            if(use_mean_decay_lr) {
	            meanLr = mean_lr / pow((iter + init_iter) * num_instance + i + 1.0, decay_mean_lr_rate);
		        // fprintf(stderr, "\t### lr: %f ###\n", lr);
            }

            if(use_covar_decay_lr) {
	            covarLr = covar_lr / pow((iter + init_iter) * num_instance + i + 1.0, decay_covar_lr_rate);
		        // fprintf(stderr, "\t### lr: %f ###\n", lr);
        }

        string range_str = "";
	    unsigned start = random_index[i] * batch_size;
	    unsigned end = (random_index[i] + 1) * batch_size - 1;
	    

	    if(start >= orig_num_instance) continue;
	    if(start > end) continue;

	    string s1 = "";
	    string s2 = ""; 

	    if (end >= orig_num_instance) end = orig_num_instance - 1;
	    i2str(start, s1);
	    i2str(end, s2);

        range_str += (s1 + ":" + s2);
	    fprintf(stderr, "range: %u %u\n", start, end);


            vector<double> d_probs;
            vector<double> n_probs;
            

            //denominator train/update
            GM_Parms = GM_Parms_d;
            //getGmParmsValue();
            train(myjt_d, range_str.c_str(), d_probs, iter == 0);
            //getGmParmsValue();
            GM_Parms_d = GM_Parms;



            //numerator train/update
            GM_Parms = GM_Parms_n;
            //getGmParmsValue();
            train(myjt_n, range_str.c_str(), n_probs, iter == 0);
            //getGmParmsValue();
            GM_Parms_n = GM_Parms;

            double cond_prob = 0.0;
            for(unsigned j=0; j< d_probs.size(); ++ j) {
                cond_prob += (n_probs[j] - d_probs[j]);
            }

            cond_prob_iter += cond_prob;

            fprintf(stderr, "### Conditional Log Likelihood of Current Batch: %f\n", cond_prob);

            updateBoth(lr, meanLr, covarLr, adagrad, use_adagrad);
            //getGmParmsValue();
	    
        }

        fprintf(stderr, "### Conditional Log Likelihood of Current Iteration: %f Per Utterance: %f\n", cond_prob_iter, cond_prob_iter / orig_num_instance);

        if(iter > 0) {
	        lldp = cond_prob_iter - cond_prob_iterPrev;

	        if(lldp < 0) {
	            lldp = -lldp / cond_prob_iterPrev;
            }
	        else {
	            lldp = lldp / cond_prob_iterPrev;
            }

	        fprintf(stderr, "CLL: %f, Per Utterance CLL: %f, LLDP: %f\n", cond_prob_iter, cond_prob_iter / orig_num_instance, lldp);
	      
	    } else {
	        fprintf(stderr, "CLL: %f, Per Utterance CLL: %f\n", cond_prob_iter, cond_prob_iter / orig_num_instance);
	    }
	    cond_prob_iterPrev = cond_prob_iter;



	string iter_index = "";
	i2str(iter, iter_index);
	writeTrainedOutput(iter_index);

    }


    exit_program_with_status(0);
}
