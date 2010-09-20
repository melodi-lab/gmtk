/**
* \file Multivariate MI calculation program -- Calculates the mutual
* information (unconditional or conditional on a discrete "hidden"
* variable) for arbitrary feature positions in a sequence of feature
* vectors.
*

*
* Input format: Data is assumed to be a collections of "sentences,"
* each of which is a sequence of vectors of fixed dimension.


  Karim Filali <karim@cs.washington.edu>

*/


#if HAVE_CONFIG_H
#include <config.h>
#endif


#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <errno.h>
#include <cstring>
#if HAVE_VALUES_H
#include <values.h>
#endif
#include <cmath>
#include <cassert>
#include <ctime>
#include <signal.h>

#include "general.h"
#include "error.h"
#include "range.h"
#include "mixNormal.h"
#include "mixNormalCollection.h"
#include "readRange.h"
#include "GMTK_ObservationMatrix.h"
#include "tests.h"
#include "arguments.h"
#include "GMTK_WordOrganization.h"
#include "rand.h"

///////////////////  Defines ////////////////////////////////////////////
#define _KEEP_MG_ 1
/////////////////////////////////////////////////////////////////////////


///////////////////  Global variables ///////////////////////////////////
RAND rnd(false);

static const char* program_name;
ObservationMatrix  globalObservationMatrix;
ObservationMatrix  globalLabelMatrix;


#define MAX_NUM_OBS_FILES (5)

char       *Input_Fname[MAX_NUM_OBS_FILES] = { NULL,NULL,NULL,NULL,NULL };  // Input file name.
unsigned    Num_Floats[MAX_NUM_OBS_FILES]  = { 0, 0, 0, 0, 0 };
unsigned    Num_Ints[MAX_NUM_OBS_FILES] = { 0, 0, 0, 0, 0 };
char *Float_Range_Str[MAX_NUM_OBS_FILES] = { "all", "all", "all", "all", "all" };
char *Int_Range_Str[MAX_NUM_OBS_FILES]   = { "all", "all", "all", "all", "all" };
char *Fmt_Str[MAX_NUM_OBS_FILES]  = { "pfile", "pfile", "pfile", "pfile", "pfile" };

char       *Output_Fname                = NULL;

int         Min_Num_Consecutive_Labels  = 0;

int         Num_Iterations_Between_Saves = 100; // Number of iteration between parameter saves
float       Cov_Noise_Constant          = 1e-6; // adds a small amount of noise to the diagonal entries of covariance matrices to prevent numerical issues.
float       Cov_Noise_Max_Rand          = 5e-7; // we add a random number between covAddConst-covAddEpsilon and covAddConst+Epsilon
double      Clamp_Covariance            = 1e-10;

double      MCVR                        = 20.0;

#ifdef INTV_WORDS_BIGENDIAN
bool     Swap[MAX_NUM_OBS_FILES] = {true,true,true,true,true};
#else
bool     Swap[MAX_NUM_OBS_FILES] = {false,false,false,false,false};
#endif

char*    Action_If_Diff_Num_Frames_Str[MAX_NUM_OBS_FILES]={"er","er","er","er","er"};   //
unsigned Action_If_Diff_Num_Frames[MAX_NUM_OBS_FILES]={FRAMEMATCH_ERROR,FRAMEMATCH_ERROR,FRAMEMATCH_ERROR,FRAMEMATCH_ERROR,FRAMEMATCH_ERROR};
char*    Action_If_Diff_Num_Sents_Str[MAX_NUM_OBS_FILES]={"te","te","te","te","te"};
unsigned Action_If_Diff_Num_Sents[MAX_NUM_OBS_FILES]={SEGMATCH_TRUNCATE_FROM_END,SEGMATCH_TRUNCATE_FROM_END,SEGMATCH_TRUNCATE_FROM_END,SEGMATCH_TRUNCATE_FROM_END,SEGMATCH_TRUNCATE_FROM_END};   //

bool     Cpp_If_Ascii        = true;
char*    Cpp_Command_Options = NULL;

bool     Seed                = true;

bool     Verbose             = false; // Print a lot of status messages.
bool     Quiet               = false; // Don't print any status message.  Overrides verbose.
bool     Print_Help                = false;


char *MI_Tuple_File = NULL;          // File name listing all the tuples we want to compute the mi/entropy of.

char *Input_MG_Fname  = NULL;       // parameter mg input file name
char *Output_MG_Fname = NULL;       // parameter mg output file name
char *Kmeans_Input_MG_Fname  = NULL;      // parameter kmeans input file name
char *Kmeans_Output_MG_Fname = NULL;      // parameter kmeans output file name

char *Sentence_Range_Str        = "all";         // sentence range string
char *Kmeans_Sentence_Range_Str = "all";     // kmeans sentence range string
char *MI_Sentence_Range_Str     = "all";         // mi data sentence range string
//char *Frame_Range_Str[MAX_NUM_OBS_FILES] = {NULL,NULL,NULL,NULL,NULL};   // per stream per sentence range string
char *Frame_Range_Str[MAX_NUM_OBS_FILES] = {"all","all","all","all","all"};   // per stream per sentence range string

int   Label_Position = -1;                   // If computing conditional mi and the labels are in the same pfile as features, labpos is the index of the label.  -1 means we are using the last discrete feature in each frame of the pfile as the label.

unsigned Num_Mixtures                  = 1;       // By default use one mixture to estimate densities i.e. assumes linear dependencies.
double   Log_Likelihood_Perc_Diff      = 0.01;    // When the change in log likelihood is less than lldp, we assume EM has converged for this tuple.
unsigned Max_Num_Kmeans_Iterations     = 3;       // Number of k-means iterations to perform.
unsigned Max_Num_EM_Iterations         = 100;     // The maximum number of EM iterations per tuple.
unsigned Num_Samples_Law_large_Numbers = 10000;   // The number of samples to use when using the law of large numbers to calculate MI.


bool Dont_Backup_MG_Files = false;               // When true we overwrite existing mg files.
bool Activate_All_MGs     = false;     // When true and starting from an existing mg file, all tuples are assumed not to have converged.
bool Use_Data_For_MI_Estimation = false;                 // When computing MI, use the same (or a subset) of the data used for estimating the densities instead of generating samples from the densities.
int  Num_Active_To_Inactive_Changes_For_Save = 1;     // number of active->inactive changes for a save to occur  (not counting any other timeouts).


bool Skip_EM     = false;               // If true, we don't run EM.
bool Skip_Kmeans = false;           // If true we skip kmeans.  Not implemented yet.  The main issue is to have an alternate way to initialize EM and kmeans (even if for one iteration) seems teh smartest way to go about it.

unsigned  Start_Skip  = 0;
unsigned  End_Skip    = 0;


char    *Per_Stream_Transforms[MAX_NUM_OBS_FILES]={NULL,NULL,NULL,NULL,NULL};   //
char    *Post_Transforms=NULL;

char    *Ftr_Combo_Str="none";
unsigned Ftr_Combo=FTROP_NONE;

int    numSecondsPerPrint       = NUM_SECONDS_PER_PRINT;
int    numSecondsPerSentPrint   = NUM_SECONDS_PER_SENT_PRINT;
int    printFrequency           = PRINT_FREQUENCY;
int    sentPrintFrequency       = SENT_PRINT_FREQUENCY;
int    activePrintFrequency     = ACTIVE_PRINT_FREQUENCY;
int    minTimePerPrintNumActive = MINTIMEPERPRINTNUMACTIVE;


unsigned Max_Num_Kmeans_Rerands = 100; // not used yet

 // corresponding to the tuple with index distNumToDump is written out
#if DEBUG
int distNumToDump=0;
#else
int    distNumToDump            = -1; // For debugging: distribution data
#endif
bool   mgBinFormat              = true; // do we read/write mg parameters in binary?
int    debug_level = 0;

bool Dont_Run = false;

bool Detailed_Output=false;

char* Label_Range_Str=NULL;

bool Marginalize_First_Parent_Out=false;

//////  At some point will be deleted /////
const char *lofs[]                   = { NULL };
unsigned    lnfs[]                   = { 0 };
unsigned    lnis[]                   = { 0 };
const char *lfrs[]                   = { "all" };
const char *lirs[]                   = { "all" };
const char *lfmts[]                  = { "pfile" };
bool        liswps[]                 = { false };
bool        gotLabelFile             = false;
bool        lswap                    = false;       // If true we need to byte swap the label file.


const char *label_fname              = NULL;    // Label pfile name (if any).  // This is deprecated.  To specify lables, now, just make sure there is at least an int across all input streams.  The last int by default will be used as the label, unless overriden using the labpos option.

int Num_Active_To_Stop = 0;        // Not used.  Might be useful to implement if it appears we are spending much time estimating the densities for a few tuples.
////////////////////////////////////////////

/////////////////////////////////////////////////////



Arg Arg::Args[] = {
  Arg("i",        Arg::Req, Input_Fname,"Input file.  Replace the X with the observation file number.",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("fmt",      Arg::Opt, Fmt_Str ,"Format of input file",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("nf",       Arg::Opt, Num_Floats,"Number of floats in input file",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("ni",       Arg::Opt, Num_Ints,"Number of ints (labels) in input file",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("fr",       Arg::Opt, Float_Range_Str,"Float range to use",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("ir",       Arg::Opt, Int_Range_Str,"Int range to use",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("o",        Arg::Opt, Output_Fname,"Output file (- for sdtdout)"),
  Arg("miTupleFile", Arg::Req, MI_Tuple_File, "File specifying the tuples to compute the MI/entropy of"),

  Arg("sr",       Arg::Opt, Sentence_Range_Str, "Sentence range for EM estimation of the probability densities"),
  Arg("ksr",      Arg::Opt, Kmeans_Sentence_Range_Str, "Sentence range for kmeans initial estimation of the probability densities"),
  Arg("lr",       Arg::Opt, Label_Range_Str, "Label range to condition on."),

  Arg("m",        Arg::Opt, Num_Mixtures, "Number of mixtures"),
  Arg("lll",      Arg::Opt, Num_Samples_Law_large_Numbers, "Number of samples to use for the Law of Large Numbers MI estimation"),
  Arg("data",     Arg::Tog, Use_Data_For_MI_Estimation, "Re-use the input data for MI estimation instead of generating LLL samples"),
  Arg("mr",       Arg::Opt, MI_Sentence_Range_Str, "Sentence range for the MI calculation using the data method above"),
  Arg("lldp",     Arg::Opt, Log_Likelihood_Perc_Diff, "Log likelihood percent difference below which density estimation is assumed to have converged"),
  Arg("maxEMIters",     Arg::Opt, Max_Num_EM_Iterations, "Maximum number of EM iterations"),
  Arg("maxKmeansIters", Arg::Opt, Max_Num_Kmeans_Iterations, "Maximum number of kmeans iterations"),
  Arg("noEM",     Arg::Tog, Skip_EM, "Skip EM density estimation"),
  Arg("noKmeans", Arg::Tog, Skip_Kmeans, "Skip kmeans"),
  Arg("krerands", Arg::Opt, Max_Num_Kmeans_Rerands, "Maximum number of kmeans cluster random re-assignements when null clusters are found"),
  Arg("nacps",    Arg::Opt, Num_Active_To_Inactive_Changes_For_Save,"Number of active->inactive changes for a save"),
  Arg("nips",     Arg::Opt, Num_Iterations_Between_Saves, "Number of iteartions between parameter saves"),

  Arg("clampCov",      Arg::Opt, Clamp_Covariance, "Value to clamp the covariance entries to if they fall below it."),
  Arg("addCovConst",   Arg::Opt, Cov_Noise_Constant, "Non random value to add to the diagonal entries of covariance matrices to avoid numerical problems (notably, when the entries become too small)"),
  Arg("addCovMaxRand", Arg::Opt, Cov_Noise_Max_Rand, "Value of maximum randon number to add/substract from the additive constant above.  Useful to prevent degenerate covariance matrices (identical rows for example)"),
  Arg("mcvr",          Arg::Opt, MCVR, "Vanishing ratio.  Needs to be  > 1."),
  Arg("labPosition",   Arg::Opt, Label_Position, "Position of the int used as a label when computing conditional MI/entropy"),
  Arg("mlps",          Arg::Opt, Min_Num_Consecutive_Labels, "Minimum consecutive labels per segment to compute conditional MI with"),

  Arg("pi",       Arg::Opt, Input_MG_Fname,  "Input mg filename"),
  Arg("po",       Arg::Opt, Output_MG_Fname,  "Output mg filename"),
  Arg("kpi",      Arg::Opt, Kmeans_Input_MG_Fname, "Input kmeans mg filename"),
  Arg("kpo",      Arg::Opt, Kmeans_Output_MG_Fname, "Output kmeans mg filename"),
  Arg("mgBinFormat",  Arg::Tog, mgBinFormat, "MG files are written out in binary"),
  Arg("nobak",        Arg::Tog, Dont_Backup_MG_Files, "Do not backup mg files before overwriting them"),
  Arg("activateAll",  Arg::Tog, Activate_All_MGs, "Activate all MGs in input mg file regardless of file status"),

  Arg("iswap",    Arg::Tog, Swap,"Do byte swapping on the input file",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("fdiffact", Arg::Opt, Action_If_Diff_Num_Frames_Str ,"Action if different number of frames in streams: error (er), repeat last frame (rl), first frame (rf), segmentally expand (se), truncate from start (ts), truncate from end (te)",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("sdiffact", Arg::Opt, Action_If_Diff_Num_Sents_Str ,"Action if different number of sentences in streams: error (er), truncate from end (te), repeat last sent (rl), and wrap around (wa).",Arg::ARRAY,MAX_NUM_OBS_FILES),
  Arg("cppifascii",        Arg::Tog, Cpp_If_Ascii,"Pre-process ASCII files using CPP"),
  Arg("cppCommandOptions", Arg::Opt, Cpp_Command_Options,"Additional CPP command line"),
  Arg("seed",     Arg::Tog, Seed, "Seed the random number generator"),
  //  Arg("detailedOutput", Arg::Tog, Detailed_Output, "For each tuple specification output, in addition to MI, the following statistics the # of data samples used, # of EM iterations, and the entropies of the X, Y and (X,Y) sets respectively."),
  Arg("verbose",  Arg::Tog, Verbose, "Be very verbose"),
  Arg("q",        Arg::Tog, Quiet,"Do not print any diagnostic messages"),
  Arg("norun",    Arg::Tog,  Dont_Run, "Don't run; just print the values of the arguments"),
  Arg("marginalizeFirstParent",     Arg::Tog,  Marginalize_First_Parent_Out, "Marginalize the first parent out when computing MI."),
  Arg("help",     Arg::Tog, Print_Help,"print this message"),
  // The argumentless argument marks the end of the above list.
  Arg()
};





/////////////////////////////////////////////////////////////////////////


// not yet used ///////////////////////
/**
* read a sentence into memory
* @param sentence sentence to be read im memory
* @param numFrames output parameter that stores the number of frames in the sentence
* @param numSamples output parameter storing the number of samples
*/
int readSentenceUncond(Range::iterator sentence,
		       size_t &numFrames,
		       Range &lrrng);



int readSentenceCond(Range::iterator sentence,
		     size_t &numFrames,
		     size_t &numSamples,
		     Range &lrrng,
		     int labpos,
		     unsigned &frameStart,
		     unsigned &firstFrame);

///////////////////////////////////////////////////

/**
 * reads features and labels and performs error checking.
 * NOTE:  loadSegment() seems to be called more than necessary.
 */
int readFeatures(Range::iterator krit, size_t &n_frames,
		 size_t &n_samps,
		 Range &lrrng, int labpos,
		 unsigned &frameStart,unsigned &firstFrame) {

  bool segAlreadyLoaded = false;
  if(frameStart != 0) segAlreadyLoaded = true;

  size_t n_labs;
  n_labs = (gotLabelFile?globalLabelMatrix.numDiscrete():globalObservationMatrix.numDiscrete());
  const bool stateCondMI = !lrrng.full();

  if (stateCondMI && n_labs < 1) {
    error("For conditional MI, number of discrete features per frame in pfile must be at least one.");
  }

  ObservationMatrix* obsMat;
  if(gotLabelFile)
    obsMat = &globalLabelMatrix;
  else
    obsMat = &globalObservationMatrix;
  if(!segAlreadyLoaded)
    obsMat->loadSegment((const unsigned)(*krit));
  n_frames = obsMat->numFrames();

  if ( n_frames == SIZET_BAD )
    error("%s couldn't find number of frames at sentence %lu in input pfile.\n", program_name, (unsigned long) *krit);

  if (!stateCondMI) {
    firstFrame = 0;
    n_samps = n_frames;
    if(n_samps == 0) return NO_DATA;
    return DONE;
  }
  else {
    unsigned label;
    unsigned frameno;
    int pos = (unsigned) labpos;
    if(pos == -1) //use default: last discrete feature
      pos =  obsMat->numFeatures() - 1;

    int numFound = 0;
    for(frameno = frameStart; frameno < n_frames; ++frameno) {
      label =  obsMat->unsignedAtFrame(frameno,(const unsigned) pos);
      frameStart = frameno;
      while (lrrng.contains(label)) {
	++numFound;
	++frameno;
	if(frameno >= n_frames) break;
	label =  obsMat->unsignedAtFrame(frameno,(const unsigned) pos);
      }
      if (numFound == 0 || numFound < Min_Num_Consecutive_Labels)
	numFound = 0; //reset
      else break;
    }
    //At this point we have read a continuous chunk
    //Case 1:  We haven't reached the last frame
    //         Everything is good, we return DATA_LEFT
    //Case2:   We've reached the last frame:
    //   Case 2.a: numFound > 0 : return DONE
    //   Case 2.b: numFound ==0 : return NO_DATA
    //
    n_samps = numFound;
    firstFrame = frameStart;
    if(frameno != n_frames) {
      if(!segAlreadyLoaded) {
	globalObservationMatrix.loadSegment((const unsigned)(*krit));
	if(gotLabelFile) { //check that they have the same number of frames
	  unsigned obs_n_frames = globalObservationMatrix.numFrames();
	  if(obs_n_frames != n_frames)
	    error("The number of observation frames is different from that of label frames in sentence %d\n",(int)(*krit));
	}
      }
      frameStart = frameno;  //next time we start from here
      return DATA_LEFT;
    }
    else if(numFound > 0) {
      if(!segAlreadyLoaded) {
	globalObservationMatrix.loadSegment((const unsigned)(*krit));
	if(gotLabelFile) { //check that they have the same number of frames
	  unsigned obs_n_frames = globalObservationMatrix.numFrames();
	  if(obs_n_frames != n_frames)
	    error("The number of observation frames is different from that of label frames in sentence %d\n",(int)(*krit));
	}
      }
      return DONE;
    }
    else { //numFound == 0
      return NO_DATA;
    }
  } //end else (stateCondMI)
  return DONE;
}



bool usr1_terminate = false;
void sigusr1(int flag) {
  usr1_terminate = true;
  flag = 3; // keep compiler from complaining
}

// use this for pmake.
bool usr2_terminate = false;
void sigusr2(int flag) {
  usr2_terminate = true;
  flag = 3; // keep compiler from complaining
}
void sigexit(int flag) {
  exit(-1);
  flag = 3;
}



/**
* estimates the mixture of Gaussians and calculates MI
*
*/

static void multivariateMI(FILE *mi_ofp, // where to put output MI values
	       FILE *pi_fp,  // where to get input MG params
	       FILE *po_fp,  // where to place output MG params
	       FILE *kpi_fp,  // where to get input KMEANS params
	       FILE *kpo_fp,  // where to place output KMEANS params
	       Range &srrng,
	       Range &lrrng,
	       Range &kmeansrng,
	       Range &mirng,
	       const char* tuple_fname,
	       unsigned numMixtures,
   	       unsigned numIterKmeans,
	       unsigned numMaxEpochs,
	       double lldp,
	       unsigned lll,
	       bool data,
	       int labpos,
	       const bool force_all_active,
	       int nacps,
	       int num_active_to_stop,
	       int nips,
	       bool skipKmeans,
	       bool skipEM,
               bool marginalizeFirstParentOut,
	       const bool quiet
			   ) {

  size_t n_frames, n_samps;
  const size_t n_ftrs =  globalObservationMatrix.numFeatures();

  RangeSetCollection rngSetCol(tuple_fname);
  if(rngSetCol.getSize() == 0) {
    error("Empty MI range set\n");
  }

  FILE* rangeFileFP;
  if( (rangeFileFP = fopen(tuple_fname,"r")) == NULL) {
    error("Couldn't open range file in multivariate_mi.cc:mutivariateMI()\n");
  }

  DBGFPRINTF((stderr,"Starting multivariateMI: numIterKmeans = %d\n",numIterKmeans));

  bool fullCovar = true;
  MixNormalCollection mg(rngSetCol, numMixtures, numIterKmeans, fullCovar,Cov_Noise_Constant,Cov_Noise_Max_Rand,Clamp_Covariance);
  mg.setMixtureCoeffVanishNumber(MCVR);

  int n_mis = (int) rngSetCol.getSize();
  int numActive = (int) rngSetCol.getSize();
  int prevNumActive = numActive;

  // set signals so user can stop the EM iterations
  // but we still don't loose the work.
#ifdef HAVE_SIGSET
  if (sigset(SIGUSR1,sigusr1) == SIG_ERR)
    error("Can't set sigusr1 signal.");
  if (sigset(SIGUSR2,sigusr2) == SIG_ERR)
    error("Can't set sigusr2 signal.");
  if (sigset(SIGXCPU,sigexit) == SIG_ERR)
    error("Can't set SIGXCPU signal.");
  if (sigset(SIGTERM,sigexit) == SIG_ERR)
    error("Can't set SIGTERM signal.");
#else
  if (signal(SIGUSR1,sigusr1) == SIG_ERR)
    error("Can't set SIGUSR1 signal.");
  if (signal(SIGUSR2,sigusr2) == SIG_ERR)
    error("Can't set SIGUSR2 signal.");
  if (signal(SIGXCPU,sigexit) == SIG_ERR)
    error("Can't set SIGXCPU signal.");
  if (signal(SIGTERM,sigexit) == SIG_ERR)
    error("Can't set SIGTERM signal.");

#endif

  // For debugging purposes: use -dumpdist <tuple num> to write out the distribution for the specified tuple
  if(distNumToDump >=0) {
    char dataFile[25];
    for(unsigned mixNum=0;mixNum<(unsigned)rngSetCol.getSize();++mixNum) {
      sprintf(dataFile,"DUMPED_DATA_POINTS.OUT.%d",mixNum);
      FILE* ofp=fopen(dataFile,"w");
    if(ofp == NULL) {
      fprintf(stderr,"Could not open dump file for writing\n");
      exit(-1);
    }

    cout<<"Dumping distribution data for mixture # "<<mixNum<<endl;
    dumpDistribSampleData(ofp,
			  &globalObservationMatrix,
			  rngSetCol,
			  lrrng,
			  kmeansrng,
			  numMixtures,
			  numIterKmeans,
			  labpos,
			  mixNum,
			  quiet);
    }
  } // end if(distNumToDump >=0)
  ////////////////////////////////////////////////////////////////////////////

  if (pi_fp != NULL && fsize(pi_fp) > 0) {
    prevNumActive = numActive =
      mg.readCurParams(pi_fp,force_all_active,mgBinFormat);
  }
  else if(!skipKmeans) {
    if(kpi_fp != NULL && fsize(kpi_fp) > 0) {
      cout<<"Reading input kmeans file...\n";
      mg.readCurKMeansParams(kpi_fp,mgBinFormat);
    }
    else {
      if(!quiet) cout<<"Running kMeans.\n";
      DBGFPRINTF((stderr,"Parameters passed to kmeans: numMixtures = %d, numIterKmeans = %d, labpos = %d\n",numMixtures, numIterKmeans,labpos));
      mg.kmeans(&globalObservationMatrix, rngSetCol, lrrng, kmeansrng, numMixtures, numIterKmeans,labpos,quiet);
      if(!quiet) cout << "Finished kMeans.\n";
      if(kpo_fp != NULL) {
	if(!quiet) cout<<"Writing KMeans parameters...\n";
	mg.writeCurKMeansParams(kpo_fp,mgBinFormat);
      }
    }
    if(!quiet) cout<<"Converting kmeans parameters to mg ones...\n";
    mg.calcB();
    if(!quiet) cout<<"Writing converted kmeans parameters to mg file...\n";
    mg.writeCurParams(po_fp,mgBinFormat);
  }
  else {
	if(!quiet) printf("Skipping kmeans.\n");
  }



  time_t timeOfLastPrint = time(0) - NUM_SECONDS_PER_PRINT - 1;
  time_t timeOfLastSaveParams = time(0);
  time_t timeOfLastPrintNumActive = time(0) - MINTIMEPERPRINTNUMACTIVE -1;
  unsigned prevSaveIter=0;
  unsigned em_iter=0;
  double maxDist=0.0, aveDist=0.0, minDist=0.0;
  int readStatus;
  unsigned frameStart,firstFrame;

  if(!skipEM) {  // Perform EM
    // do ... while the number of EM iterations is less than some maximum and convergence has not been achieved on all tuples.
    if(!quiet) {
      printf("Starting EM.\n"); fflush(stdout);
    }
    // replaced do loop whith while to avoid case in which we read an mg file with no active tuple and still
    // iterate needlessly.
    while (numActive != 0 && em_iter < numMaxEpochs) {
      em_iter++;

      if (!quiet)
	if  ( (time(0)-timeOfLastPrint) > numSecondsPerPrint || (em_iter % printFrequency == 0) ) {
	  printf("Iter %d: Starting Iter\n",em_iter);
	  fflush(stdout);
	  timeOfLastPrint = time(0);
	}

      mg.startEpoch();  // Initialize EM data structures
      // Iterate over sentences and accumulate EM statistics
      for ( Range::iterator srit = srrng.begin(); !srit.at_end(); srit++ ) {
	if ( ! quiet ) {
	  if ((time(0)-timeOfLastPrint) > numSecondsPerSentPrint || *srit % sentPrintFrequency == 0) {
	    printf("Iter %d, sentence %d\n",em_iter,(*srit));
	    fflush(stdout);
	    timeOfLastPrint = time(0);
	  }
	}

	// Read sentence in
	frameStart = 0;
	do{
	  readStatus =
	    readFeatures(srit, n_frames, n_samps,
			 lrrng,labpos, frameStart,firstFrame);
	  if(readStatus == NO_DATA) break;  //no frames were read
	  mg.addToEpoch(&globalObservationMatrix, n_ftrs, n_frames, n_samps, firstFrame,rngSetCol);
	} while(readStatus == DATA_LEFT);

      } // end of for loop that iterates overs sentences

      int rangeSpecNum;
      if( mg.noSamplesFound(rangeSpecNum) ) {
      error("ERROR:  There were no samples for at least one range spec (the %d th one).  Possible causes:  the label provided does not exist in the label file or there are too few frames with that label.\n",rangeSpecNum);
      }

      mg.endEpoch();  // Finish accumulating statistics and update EM parameters

      prevNumActive = numActive;
      numActive=mg.reComputeNumActive(maxDist,aveDist,minDist,lldp,em_iter);

      if (!quiet &&
	  ( ( numActive < prevNumActive ) ||
	    ( (time(0) - timeOfLastPrintNumActive ) > minTimePerPrintNumActive ) ||
	    ( em_iter % activePrintFrequency == 0 )
	    )
	  ){
	printf("Iter %d: NA=%d/%d(%.0f%%), PNA=%d, dist Max(%e) Avg(%e) Min(%e)\n",
	       em_iter,
	       numActive,n_mis,100*numActive/(double)n_mis,
	       prevNumActive,
	       maxDist,
	       aveDist,
	       minDist);
	fflush(stdout);
	timeOfLastPrint = timeOfLastPrintNumActive = time(0);
      }

      if (usr2_terminate) {
	// If we got a sigusr2 recently, then we expect to
	// soon be killed (by pmake) but we don't want to be killed in
	// the middle of saving the parameters. So, instead of saving, we forfeit
	// the work done during this em_iter for safety's sake.
	printf("Iter %d: Not Saving Mixture Parameters Since Received SIGUSR2.\n",em_iter); fflush(stdout);
      }
      else if (po_fp != NULL &&
	       ((prevSaveIter+nips <= em_iter) ||
		(numActive+nacps <= prevNumActive) ||
		(numActive == num_active_to_stop) ||
		((time(0) - timeOfLastSaveParams) > MINTIMEPERPARMSAVE))) {
	// save the current parameters.
	prevSaveIter=em_iter;
	mg.writeCurParams(po_fp,mgBinFormat);
	timeOfLastSaveParams = time(0);
	if(!quiet) {
	  printf("Iter %d: Saving Mixture Parameters.\n",em_iter); fflush(stdout);
	}
#if DEBUG
	mg.dumpCurIterParams(em_iter);
#endif
      }
    }

    if (usr2_terminate) {
      printf("Iter %d: Exiting early with failure due to received SIGUSR2\n",em_iter);
      exit (EXIT_FAILURE);
    }

    if(!quiet) printf("Finished computing mixtures.\n");

  }  // end if(!skipEM)
  else {
    cout<<"Skipping EM.\n";
  }

  // Compute MI quantities from the learned densities
  if( data ) // Use the original data (or a subset specfied by mirng) to compute MI
    mg.computeMIUsingData(globalObservationMatrix,rngSetCol, mirng, quiet, mi_ofp,lrrng,labpos,rangeFileFP);
  else       // Sample from the learned densities (lll is the number of samples to use)
    mg.computeMI(mi_ofp, lll, rngSetCol,rangeFileFP,marginalizeFirstParentOut);

  // restore signals.
  if (signal(SIGUSR1,SIG_DFL) == SIG_ERR)
    error("Can't unset SIGUSR1 signal.");
  if (signal(SIGUSR2,SIG_DFL) == SIG_ERR)
    error("Can't unset SIGUSR2 signal.");
  if (signal(SIGTERM,SIG_DFL) == SIG_ERR)
    error("Can't unset SIGTERM signal.");
  if (signal(SIGXCPU,SIG_DFL) == SIG_ERR)
    error("Can't unset SIGXCPU signal.");
}



/**
 * the main routine of the multivariate-mi program program
 *
 * @param argc the number of arguments
 * @param argv the string array of arguments
 * @return 0 if everything goes well.
 */
int main(int argc, const char *argv[]) {

  Range *sr_rng;                     // sentence range
  Range *mi_rng;                     // mi data range
  Range *kmeans_rng;                 // kmeans sentence range
  unsigned    ifmt[MAX_NUM_OBS_FILES];
  int num_files=0;

  Range *lr_rng;                     // label range

  ////// Figure out the Endian of the machine this is running on and set the swap defaults accordingly /////
  bool doWeSwap;

  ByteEndian byteEndian = getWordOrganization();
  switch(byteEndian) {
  case BYTE_BIG_ENDIAN:
    doWeSwap=false;
    break;
  case BYTE_LITTLE_ENDIAN:
    doWeSwap=true;
    break;
  default:
    // We weren't able to figure the Endian out.  Leave the swap defaults as they are.
#ifdef INTV_WORDS_BIGENDIAN
    doWeSwap=true;
#else
    doWeSwap=false;
#endif
  }

  for(int i=0; i<MAX_NUM_OBS_FILES; ++i) {
    Swap[i]=doWeSwap;
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////

  ///////  Parse Arguments //////
  bool successful_parse = Arg::parse(argc,(char**)argv);

  if(Dont_Run || Print_Help) {
    Arg::usage();
    exit(0);
  }

  if(!successful_parse) {
    Arg::usage();
    exit(-1);
  }


  if(Seed) {
    rnd.seed();
    srand((unsigned)(time(NULL)));
  }

  // TODO: put all the checks below in a new function

  if(MCVR < 1) {
    fprintf(stderr,"Vanishing ration (-mcvr) has to be greater than 1.0\n");
    Arg::usage();
    exit(-1);
  }

  for (int i=0;i<MAX_NUM_OBS_FILES;i++) {
    num_files += (Input_Fname[i] != NULL);
    if (strcmp(Fmt_Str[i],"htk") == 0)
      ifmt[i] = HTK;
    else if (strcmp(Fmt_Str[i],"binary") == 0)
      ifmt[i] = RAWBIN;
    else if (strcmp(Fmt_Str[i],"ascii") == 0)
      ifmt[i] = RAWASC;
    else if (strcmp(Fmt_Str[i],"pfile") == 0)
      ifmt[i] = PFILE;
    else
      error("ERROR: Unknown observation file format type: '%s'\n",Fmt_Str[i]);
  }

  for(int i=0; i < MAX_NUM_OBS_FILES; ++i) {
    if(Output_Fname!=NULL && Input_Fname[i] !=NULL && strcmp(Input_Fname[i],Output_Fname)==0) {
      error("Input and output filenames cannot be the same.");
    }
  }

  for(int i=0; i < MAX_NUM_OBS_FILES; ++i) {
   if(Input_Fname[i]!=NULL) {
     if (strcmp(Action_If_Diff_Num_Frames_Str[i],"er") == 0)
       Action_If_Diff_Num_Frames[i] = FRAMEMATCH_ERROR;
     else if (strcmp(Action_If_Diff_Num_Frames_Str[i],"rl") == 0)
       Action_If_Diff_Num_Frames[i] = FRAMEMATCH_REPEAT_LAST;
     else if (strcmp(Action_If_Diff_Num_Frames_Str[i],"rf") == 0)
       Action_If_Diff_Num_Frames[i] = FRAMEMATCH_REPEAT_FIRST;
     else if (strcmp(Action_If_Diff_Num_Frames_Str[i],"se") == 0)
       Action_If_Diff_Num_Frames[i] = FRAMEMATCH_EXPAND_SEGMENTALLY;
     else if (strcmp(Action_If_Diff_Num_Frames_Str[i],"ts") == 0)
       Action_If_Diff_Num_Frames[i] = FRAMEMATCH_TRUNCATE_FROM_START;
     else if (strcmp(Action_If_Diff_Num_Frames_Str[i],"te") == 0)
       Action_If_Diff_Num_Frames[i] = FRAMEMATCH_TRUNCATE_FROM_END;
     else
       error("ERROR: Unknown action when diff num of frames: '%s'\n", Action_If_Diff_Num_Frames_Str[i]);
   }
 }

for(int i=0; i < MAX_NUM_OBS_FILES; ++i) {
   if(Input_Fname[i]!=NULL) {
     if (strcmp(Action_If_Diff_Num_Sents_Str[i],"er") == 0)
       Action_If_Diff_Num_Sents[i] = SEGMATCH_ERROR;
     else if (strcmp(Action_If_Diff_Num_Sents_Str[i],"rl") == 0)
       Action_If_Diff_Num_Sents[i] = SEGMATCH_REPEAT_LAST;
     else if (strcmp(Action_If_Diff_Num_Sents_Str[i],"wa") == 0)
       Action_If_Diff_Num_Sents[i] = SEGMATCH_WRAP_AROUND;
     else if (strcmp(Action_If_Diff_Num_Sents_Str[i],"te") == 0)
       Action_If_Diff_Num_Sents[i] = SEGMATCH_TRUNCATE_FROM_END;
     else
       error("ERROR: Unknown action when diff num of sentences: '%s'\n",Action_If_Diff_Num_Sents_Str[i]);
   }
 }

 FILE *out_fp=NULL;
 if (Output_Fname==0 || !strcmp(Output_Fname,"-")) {
   out_fp = stdout;
 }
 else {
     if ((out_fp = fopen(Output_Fname, "w")) == NULL) {
       error("Couldn't open output file for writing.\n");
     }
 }

 // If we have a pfile, we can extract the number if features from the file directly
 for(int i=0; i < MAX_NUM_OBS_FILES; ++i) {
   if(Input_Fname[i]!=NULL) {
     if(ifmt[i]==PFILE) {
       FILE *in_fp = fopen(Input_Fname[i], "r");
       if (in_fp==NULL) error("Couldn't open input pfile for reading.");
       InFtrLabStream_PFile* in_streamp = new InFtrLabStream_PFile(debug_level,"",in_fp,1,Swap[i]);
       Num_Ints[i]=in_streamp->num_labs();
       Num_Floats[i]=in_streamp->num_ftrs();
       if (fclose(in_fp)) error("Couldn't close input pfile.");
       delete in_streamp;
     }

     if(Num_Ints[i]==0 && Num_Floats[i]==0) {
       error("The number of floats and the number of ints cannot be both zero.");
     }
   }
 }

  if (Num_Active_To_Inactive_Changes_For_Save < 1) {
      error("nacps (number of active->inactive changes per save must be >= 1");
  }
  if (Num_Active_To_Stop < 0) {
      error("-sac argument must be > 0");
  }

  // Create objects

  lr_rng = new Range(Label_Range_Str,0,MAX_LABEL_VAL);
  if(lr_rng->full() && label_fname != NULL)
    error("Cannot specify a label file when no label range is given\n");
  if (!lr_rng->full()) { //or better yet if lr_str != NULL
    // only bother to open this if the label range isn't full.
    if (label_fname!=NULL) {
      gotLabelFile = true;

      lofs[0] = label_fname;
      //lnfs[0] = 0; lnis[0] = 1;
      if(lswap == true) liswps[0] = true;
      if (lofs[0] != NULL && lnis[0] == 0)
	error("ERROR: command line must specify lni not zero");
      if(Label_Position < 0 && lnis[0] > 1) {
	printf("Using the last discrete feature of the label file as the label\n");
      }

      unsigned lifmts[1];
      if (strcmp(lfmts[0],"htk") == 0)
	lifmts[0] = HTK;
      else if (strcmp(lfmts[0],"binary") == 0)
	lifmts[0] = RAWBIN;
      else if (strcmp(lfmts[0],"ascii") == 0)
	lifmts[0] = RAWASC;
      else if (strcmp(lfmts[0],"pfile") == 0)
	lifmts[0] = PFILE;
      else
	error("ERROR: Unknown observation file format type: '%s'\n",lfmts[0]);

      globalLabelMatrix.openFiles(1,
				  (const char**)&lofs,
				  (const char**)&lfrs,
				  (const char**)&lirs,
				  (unsigned*)&lnfs,
				  (unsigned*)&lnis,
				  (unsigned*)&lifmts,
				  (bool*)&liswps);

      unsigned numLabelFeatures = globalLabelMatrix.numFeatures();
      unsigned numLabelContinuous = globalLabelMatrix.numContinuous();
      if(Label_Position != -1 &&
	 (Label_Position < (int)numLabelContinuous || Label_Position >= (int)numLabelFeatures) )
      error("Label position (%d) out of range (%d - %d): must be within the range of discrete observations of the label file\n",Label_Position,numLabelContinuous,numLabelFeatures-1);
      }
      else {
	gotLabelFile = false;
	if(Label_Position < 0) { //No label position has been given
	  printf("Using the last discrete feature of the last input file as the label\n");
	}
      }
    }

     globalObservationMatrix.openFiles(num_files,  // number of files.
                                   (const char**)&Input_Fname,
                                   (const char**)&Float_Range_Str,
                                   (const char**)&Int_Range_Str,
                                   (unsigned*)&Num_Floats,
                                   (unsigned*)&Num_Ints,
                                   (unsigned*)&ifmt,
                                   (bool*)&Swap,
                                   Start_Skip,  // startSkip
                                   End_Skip,  // endSkip  pr_rng takes care of these two
                                   Cpp_If_Ascii,
                                   Cpp_Command_Options,
                                   (const char**)&Frame_Range_Str,
                                   Action_If_Diff_Num_Frames,
                                   Action_If_Diff_Num_Sents,
                                   Per_Stream_Transforms,
                                   Post_Transforms,
                                   Ftr_Combo);

    unsigned numFeatures = globalObservationMatrix.numFeatures();
    unsigned numContinuous = globalObservationMatrix.numContinuous();
    if(!gotLabelFile && Label_Position != -1 &&
       (Label_Position < (int)numContinuous || Label_Position >= (int)numFeatures) )
      error("labpos (%d) out of range (%d - %d): must be within the range of discrete obsevations\n",Label_Position,numContinuous,numFeatures-1);

  sr_rng = new Range(Sentence_Range_Str,0,globalObservationMatrix.numSegments());
  kmeans_rng = sr_rng;
  if(Kmeans_Sentence_Range_Str != NULL)
    kmeans_rng = new Range(Kmeans_Sentence_Range_Str,0,globalObservationMatrix.numSegments());
  mi_rng = sr_rng;
  if(MI_Sentence_Range_Str != NULL)
    mi_rng = new Range(MI_Sentence_Range_Str,0,globalObservationMatrix.numSegments());


  // Open the input/output mixture of Gaussians files
  FILE *pi_fp = NULL;
  FILE *po_fp = NULL;
  if (Input_MG_Fname!=NULL && Output_MG_Fname!=NULL && !strcmp(Input_MG_Fname,Output_MG_Fname)) {

    if (!Dont_Backup_MG_Files && ((pi_fp = fopen(Input_MG_Fname, "r")) != NULL)) {
      // backup the file to : pi_fname + ".bak"
      const int bufsiz = strlen(Input_MG_Fname)+8096;
      size_t sz;
      char *buf = new char[bufsiz];
      sprintf(buf,"%s.bak",Input_MG_Fname);
      if ((po_fp = fopen(buf, "w")) == NULL)
	error("Can't open backup pi file.");
      while ((sz = fread(buf,1,8096,pi_fp))) {
	if (ferror(pi_fp))
	  error("Couldn't read buffer for backup file.");
	if (fwrite(buf,1,sz,po_fp) != sz)
	  error("Couldn't write full buffer to backup file.");
      }
      fclose(pi_fp); fclose(po_fp);
      delete [] buf;
    }
    else {
      // parameter file does not yet exist or we don't want to backup.
    }

    pi_fp = po_fp = fopen(Input_MG_Fname, "r+");
    if (pi_fp==NULL) // assume file doesn't exist.
      pi_fp = po_fp = fopen(Input_MG_Fname, "w+");
    if (pi_fp==NULL)
      error("Couldn't open i/o file for reading/writing.");
  }
  else {
    if (Input_MG_Fname!=0) {
      pi_fp = fopen(Input_MG_Fname, "r");
      if (pi_fp==NULL)
	error("Couldn't open input mg (pi) file for reading.");
    }
    if (Output_MG_Fname != 0) {
      po_fp = fopen(Output_MG_Fname, "w");
      if (po_fp==NULL)
	error("Couldn't open output mg (po) file for writing.");
    }
  }

  // Figure out which kmeans parameters input/output file names are passed.
  FILE *kpi_fp = NULL;
  FILE *kpo_fp = NULL;
  if ( Kmeans_Input_MG_Fname!=NULL ) {  // input kmeans specified
    if ( Kmeans_Output_MG_Fname!=NULL ) {  // output kmeans specified.  Cannot allow that.
      error("You cannot specify an output kmeans file along with an input one");
    }
    if ( Kmeans_Input_MG_Fname!=NULL ) {  // input mg specified
      fprintf(stderr,"mg input file specfied along with kmeans input file.  Kmeans file will be overidden if mg input file exists.\n");
    }
    kpi_fp = fopen(Kmeans_Input_MG_Fname, "r");
    if (kpi_fp==NULL)
      error("Couldn't open input kmeans file for reading.");
  }
  else if(Kmeans_Output_MG_Fname!=NULL) {
    if ( Input_MG_Fname!=0 ) {  // input mg specified
      fprintf(stderr,"Kmeans output file has no effect when mg input file is specfied i.e. kmeans is skipped altogether.\n");
    }
    if(Output_MG_Fname==NULL || strcmp(Kmeans_Output_MG_Fname,Output_MG_Fname)!=0) { // if the kmeans and mg output parameter filenames are different
      DBGFPRINTF((stderr,"Opening kpo file\n"));
      kpo_fp = fopen(Kmeans_Output_MG_Fname, "w");
      if (kpo_fp==NULL)
	error("Couldn't open input kmeans file for writing.");
    }
    else
      kpo_fp=po_fp;  // else we write to the mg output param file
  }


  // Control the amount of verbosity
  if(Verbose) {
    numSecondsPerPrint = VERBOSE_NUM_SECONDS_PER_PRINT;
    numSecondsPerSentPrint = VERBOSE_NUM_SECONDS_PER_SENT_PRINT;
    printFrequency = VERBOSE_PRINT_FREQUENCY;
    sentPrintFrequency = VERBOSE_SENT_PRINT_FREQUENCY;
    activePrintFrequency =  VERBOSE_ACTIVE_PRINT_FREQUENCY;
    minTimePerPrintNumActive = VERBOSE_MINTIMEPERPRINTNUMACTIVE;
  }


  // Do the work.
  multivariateMI(out_fp,
		 pi_fp,
		 po_fp,
		 kpi_fp,
		 kpo_fp,
		 *sr_rng,
		 *lr_rng,
		 *kmeans_rng,
		 *mi_rng,
		 MI_Tuple_File,
		 Num_Mixtures,
		 Max_Num_Kmeans_Iterations,
		 Max_Num_EM_Iterations,
		 Log_Likelihood_Perc_Diff,
		 Num_Samples_Law_large_Numbers,
		 Use_Data_For_MI_Estimation,
		 Label_Position,
		 Activate_All_MGs,
		 Num_Active_To_Inactive_Changes_For_Save,
		 Num_Active_To_Stop,
		 Num_Iterations_Between_Saves,
		 Skip_Kmeans,
		 Skip_EM,
		 Marginalize_First_Parent_Out,
		 Quiet);

  // Clean up and exit.

  if (pi_fp != NULL)
    fclose(pi_fp);
  if (po_fp != NULL && po_fp != pi_fp)
    fclose(po_fp);

  if (out_fp && fclose(out_fp))
    error("Couldn't close output file.");

  delete sr_rng;
  delete lr_rng;
  if (Kmeans_Sentence_Range_Str != NULL) delete kmeans_rng;
  if (MI_Sentence_Range_Str != NULL)     delete mi_rng;

  return 0;
}
