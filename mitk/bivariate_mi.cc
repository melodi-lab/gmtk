/*
**    $Header$
**  
**    This program performs a mutual-information (MI) analysis on the input 
**    pfile, by computing the MI between pairs of feature locations, between
**    all elements of current frame to all elements of all frames in the
**    past up to some threshold.
** 
**    Written by:
**       Jeff Bilmes <bilmes@ee.washington.edu>
**
**
**    Slight modifications to the interface to read input files using several types of formats besides pfile.  
**       Karim Filali <karim@cs.washington.edu> 
*/

// The maximum label value in a pfile label stream.
// The lr argument can not specify a label greater than
// this or the program will halt. It is safe to 
// increase this value as necessary.
#define MAX_LABEL_VAL 16384
#define MAX_CFR 20 

#define ABS_OF_INT(i) ((i)<0?(-(i)):(i))
#define MAX_OF_INT(i,j) ((i)>(j)?(i):(j))


#ifdef HAVE_NONSTANDARD_ARITHMETIC
extern "C" void nonstandard_arithmetic();
#endif


#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <values.h>
#include <math.h>
#include <signal.h>
#include <sys/types.h>
#include <time.h>

#include "sArray.h"
#include "general.h"

#  ifdef HAVE_IEEEFP_H
#include "ieeeFPsetup.h"
#  endif

#include "arguments.h"

//
// Turn on for debugging.
#define SET_IEEE_TRAPS 1

//
// for the mixture gaussian case, this is the
// minimum amount of time between parameter saves, in seconds.
//#define MINTIMEPERPARMSAVE (60*10)
#define MINTIMEPERPARMSAVE (60*20)

// for the mixture gaussian case, this is the
// minimum amount of time between printing the number of active
//#define MINTIMEPERPRINTNUMACTIVE (2)
#define MINTIMEPERPRINTNUMACTIVE (20)

// for the mixture gaussian case, this is the
// minimum amount of time between printing the iter I/II status msg.
//#define NUM_SECONDS_PER_PRINT (2)
#define NUM_SECONDS_PER_PRINT (20)


#include <assert.h>

#include "pfile.h"
#include "error.h"
#include "range.h"
#include "bp_range.h"
#include "MixBiNormal.h"
#include "rand.h"
//#include "spi.h"
#include "GMTK_ObservationMatrix.h"

#define MIN(a,b) ((a)<(b)?(a):(b))


//GLOBAL VARIBLES
static const char* program_name;
ObservationMatrix globalObservationMatrix;
ObservationMatrix globalLabelMatrix;
bool gotLabelFile = false;

#define MAX_NUM_OBS_FILES (3)
const char *ofs[MAX_NUM_OBS_FILES] = { NULL, NULL, NULL }; 
unsigned nfs[MAX_NUM_OBS_FILES] = { 0, 0, 0 };
unsigned nis[MAX_NUM_OBS_FILES] = { 0, 0, 0 };
const char *frs[MAX_NUM_OBS_FILES] = { "all", "all", "all" };
const char *irs[MAX_NUM_OBS_FILES] = { "all", "all", "all" };
const char *fmts[MAX_NUM_OBS_FILES] = { "pfile", "pfile", "pfile" };
bool iswps[MAX_NUM_OBS_FILES] = { false, false, false };

const char *lofs[] = { NULL}; 
unsigned lnfs[] = { 0};
unsigned lnis[] = { 0};
const char *lfrs[] = { "all"};
const char *lirs[] = { "all"};
const char *lfmts[] = { "pfile"};
bool liswps[] = {false };

//ARGS ARGS::Args[] = {
//  ARGS()
//};


static void
usage(const char* message = 0)
{
    if (message)
        fprintf(stderr, "%s: %s\n", program_name, message);
    fprintf(stderr, "Usage: %s <options>\n", program_name);
    fprintf(stderr,
	    "Where <options> include:\n"
	    "-help           Print this message.\n"
	    "-i <file-name>  Input file.\n"
	    "-ni1 <number of ints> for file 1.\n"
	    "-nf1 <number of floats> for file 1.\n"
	    "[-ir1 <int range>] for file 1.\n"
	    "[-fr1 <float range>] for file 1.\n"
	    "[-fmt1 <format(htk,bin,asc,pfile)>] for file 1.  Default= pfile\n"
	    "[-i2 <file-name>]   Opt. 2nd input file to be merged with 1st.\n"
	    "[-ni2 <number of ints>] for file 2.\n"
	    "[-nf2 <number of floats>] for file 2.\n"
	    "[-ir2 <int range>] for file 2.\n"
	    "[-fr2 <float range>] for file 2.\n"
	    "[-fmt2 <format(htk,bin,asc,pfile)>] for file 2.\n"
	    "[-iswap1]       Byte-swap input file 1\n"
	    "[-iswap2]       Byte-swap input file 2\n"
	    "[-lswap]        Byte-swap label file\n"
	    "-o <file-name>  Output file ('-' for stdout)\n"
	    "[-sr range]     Sentence range.  Default= all\n"
	    "-cfr bp-range   Temporal context range of frames to use.\n"
	    "[-lr range]     Subset of labels to condition on. Default= all.\n"
	    "[-labpos int]   Posit. of the label.  Default= last discrete feat.\n"
	    "[-lf file]      File containing labels. Use only with -lr option.\n"
	    "[-lni <number of ints>] for label file.\n"
	    "[-lnf <number of floats>] for label file.\n"
	    "[-lir <int range>] for label file.\n"
	    "[-labfr <float range>] for label file.\n"
	    "[-lfmt <format(htk,bin,asc,pfile)>] for the label file\n"
	    "[-entr]         Also print the zero lag self-info or entropies.\n"
	    "[-cc|-mg n t]   -cc means use corrcoeff to compute MI\n"
	    "                -mg=means use mixture of gaussians where:\n"
	    "                    n = number of mixture components\n"
	    "                    t = termination log-likelihood diff threshold\n"
	    "The following options only apply when the '-mg' option is active\n"
	    "     [-pi file]   File to obtain initial mixture parameters.\n"
	    "     [-po file]   File to place final mixture parameters.\n"
	    "     [-bfr bp-range] Base feature range of feature elements to use (def all).\n"
	    "     [-lfr bp-range] Range of lagged feature elements to use (def all).\n"
	    "     [-grid b n]  Use grid MI computation with b bins and n stds\n"
	    "     [-lll n]     Use law-of-large-nums MI comp with n samps\n"
	    "     [-sac n]     Stop when numactive=n (don't compute MI if n>0)\n"
	    "     [-nocmi]     Only finish computing mixtures, don't compute MI\n"
	    "     [-maxi n]    Stop after 'n' EM iterations, don't computeMI\n"
	    "     [-varfl f]   Variance floor, if var < f, re-randomize\n"
	    "     [-detfl f]   Determinant floor, if det <= f >= 0, re-randomize\n"
	    "     [-mcvr  f]   Mixture Coefficient Vanishing Ratio\n"
	    "     [-rerands d] Re-randomizations before mixture component drop\n"
	    "     [-rerandone] Only re-randomize the bad component (not all of them) when it goes awry\n"
	    "     [-nodrprrnd] Do *not* re-randomize all mixtures when a component drop occurs\n"
	    "     [-nacps d]   Number of active->inactive changes for a save\n"
	    "     [-nobak]     Don't create .bak files\n"
	    "     [-frcact]    Activate all MGs in input .mg file regardless of file status\n"
	    "     [-mlps d]    Minimum labels per segment to compute CMI with.\n"
	    // -ocaar is a hack option to make a copy of a .mg file activating range
	    "     [-ocaar rng]  Only Copy And Activate Range of the input .mg file\n"
	    "     [-seed]       Seed the random number generator (default = false)\n"
	    "-debug <level>  number giving level of debugging output to produce 0=none.\n"
    );
    exit(EXIT_FAILURE);
}



// ======================================================================
// ======================================================================
// -------------- Code for Linear cor-coef based MI computation ---------
// ======================================================================
// ======================================================================

class Cor_Col {

  const int n_ftrs;
  const int n_lags;
  const int n_means;
  const int n_mis;



  // mean vector E[X]
  double *ftr_sums;
  // mean squared vector E[X^2]  
  double *ftr_sumssq;
  // squared means vector E[X]^2
  double *ftr_sqsums;
  // mean cross correlation vector E[XY]
  double *ftr_sumsxy;
  // The actual MI values.
  double *ftr_mi;

  // Position where zero lag exists. 
  // E.g., ftr_sums + zero_lag_offset*n_ftrs is the
  // zero lag position.
  int zero_lag_offset;

  // number of samples used
  int counts;

public:

  Cor_Col(const int n_ftrs_a,
	  Range & cfr_rng);
  ~Cor_Col();

  void
  accumulateFrame(float *const ftr_buf_cur,
		  const int n_before,
		  const int n_after,
		  Range& cfr_rng);
  void normalize(Range & cfr_rng);
  void computeMI(Range& cfr_rng);
  void printMI(Range& cfr_rng,FILE *out_fp,const bool printEntropies);

};



Cor_Col::Cor_Col(const int n_ftrs_a,
		 Range & cfr_rng)
 : n_ftrs(n_ftrs_a),
   n_lags(cfr_rng.length()),
   n_means(n_ftrs*cfr_rng.length()),
   n_mis( n_ftrs*n_ftrs*(cfr_rng.length())
      - (cfr_rng.contains(0) ? n_ftrs*(n_ftrs-1)/2 : 0))
{
  if (!cfr_rng.contains(0)) {
    error("Context range must contain zero lag position.");
  }

  // Not very efficient, but this is only done once.
  zero_lag_offset = 0;
  for (Range::iterator cnit=cfr_rng.begin();
       !cnit.at_end();cnit++) {
    if ((*cnit) == 0)
      break;
    zero_lag_offset++;
  }

  ftr_sums = new double [n_means];
  ftr_sumssq = new double [n_means];
  ftr_sqsums = new double [n_means];
  ftr_sumsxy = new double [n_mis];
  ftr_mi = new double [n_mis];
  counts = 0;
  ::memset(ftr_sums,0,n_means*sizeof(double));
  ::memset(ftr_sumssq,0,n_means*sizeof(double));
  ::memset(ftr_sqsums,0,n_means*sizeof(double));
  ::memset(ftr_sumsxy,0,n_mis*sizeof(double));

}

Cor_Col::~Cor_Col()
{
  delete ftr_sums;
  delete ftr_sumssq;
  delete ftr_sqsums;
  delete ftr_sumsxy;
  delete ftr_mi;
}

void
Cor_Col::accumulateFrame(float *const ftr_buf_cur,
			 const int n_before,
			 const int n_after,
			 Range& cfr_rng)
{
  double* ftr_sums_p = ftr_sums;
  double* ftr_sumssq_p = ftr_sumssq;
  double* ftr_sumsxy_p = ftr_sumsxy;


  if (cfr_rng.first() < 0 && -cfr_rng.first() > n_before)
    return;
  if (cfr_rng.last() > 0 && cfr_rng.last() > n_after)
    return;
  counts++;

  for (Range::iterator cnit=cfr_rng.begin();
       !cnit.at_end();cnit++) {
    const int l = (*cnit);
    
    float *ftr_buf_lag_p;

    ftr_buf_lag_p = ftr_buf_cur + l*n_ftrs;
    
    for (int j=0;j<n_ftrs;j++) {
      // j-loop is over position in the lagged feature vector

      const double tmp = *ftr_buf_lag_p;
      *ftr_sums_p  += tmp;
      *ftr_sumssq_p += tmp*tmp;

      // k-loop is over position in the current feature vector	    
      if (l != 0) {
	float *ftr_buf_cur_p = ftr_buf_cur;
	for (int k=0;k<n_ftrs;k++) {
	  *ftr_sumsxy_p++ += (tmp)*(*ftr_buf_cur_p++);
	}
      } else {
	float *ftr_buf_cur_p = ftr_buf_cur+j;
	for (int k=(j);k<n_ftrs;k++) {
	  *ftr_sumsxy_p++ += (tmp)*(*ftr_buf_cur_p++);
	}
      }
      ftr_buf_lag_p++;

      ftr_sums_p++;
      ftr_sumssq_p++;
    }
  }
}

void 
Cor_Col::normalize(Range & cfr_rng)
{

  double *ftr_sumsxy_p = ftr_sumsxy;
  double *ftr_sums_p = ftr_sums;
  double *ftr_sumssq_p = ftr_sumssq;
  double *ftr_sqsums_p = ftr_sqsums;

  if (counts == 0)
    error("No samples used, so can't normalize.");
  const double inv_norm = 1.0/counts;
    
  for (Range::iterator cnit=cfr_rng.begin();
       !cnit.at_end();cnit++) {
    const int l = (*cnit);
    for (int j=0;j<n_ftrs;j++) {
      const double tmp = (*ftr_sums_p)*inv_norm;
      *ftr_sqsums_p++  = (tmp*tmp);
      *ftr_sums_p++ = tmp;
      *ftr_sumssq_p++ *= inv_norm;
      if (l != 0) {
	for (int k=0;k<n_ftrs;k++) {
	  *ftr_sumsxy_p++ *= inv_norm;
	}
      } else {
	for (int k=(j);k<n_ftrs;k++) {
	  *ftr_sumsxy_p++ *= inv_norm;
	}
      }
    }
  }
}

void
Cor_Col::computeMI(Range& cfr_rng)
{

  const double half_inv_log2 = 0.5/log(0.5);

  // pointers to the lag 0 case.
  double* const ftr_sums_l0 = ftr_sums + n_ftrs*zero_lag_offset;
  double* const ftr_sumssq_l0 = ftr_sumssq + n_ftrs*zero_lag_offset;
  double* const ftr_sqsums_l0 = ftr_sqsums + n_ftrs*zero_lag_offset;

  double *ftr_mi_p = ftr_mi;
  double *ftr_sumsxy_p = ftr_sumsxy;
  double *ftr_sums_p = ftr_sums;
  double *ftr_sumssq_p = ftr_sumssq;
  double *ftr_sqsums_p = ftr_sqsums;
  for (Range::iterator cnit=cfr_rng.begin();
       !cnit.at_end();cnit++) {
    const int l = (*cnit);
    for (int i=0;i<n_ftrs;i++) {
      // i-loop is over the laged values.
      
      const double Ex = *ftr_sums_p;
      const double Exx = *ftr_sumssq_p;
      const double ExEx = *ftr_sqsums_p;
      const double Exx_ExEx = Exx-ExEx;

      const size_t offset = (l==0?(i):0);
      const double *ftr_sums_l0_p = ftr_sums_l0+offset;
      const double *ftr_sumssq_l0_p = ftr_sumssq_l0+offset;
      const double *ftr_sqsums_l0_p = ftr_sqsums_l0+offset;
      for (int j=offset;j<n_ftrs;j++) {

	// j-loop is over the zero-lag values.

	const double Exy = *ftr_sumsxy_p;
	const double Ey = *ftr_sums_l0_p;
	const double Eyy = *ftr_sumssq_l0_p;
	const double EyEy = *ftr_sqsums_l0_p;


	if (l == 0 && (i==j)) {
	  // Compute the feature self entropy.
	  assert ( Ex == Ey );
	  const double var = Exy - Ex*Ey;
	  // printf("%d: var = %f\n",i,var);
	  *ftr_mi_p = 0.5*log(2.0*M_PI*M_E*var)/log(2.0);
	} else {
	  // Instead of computing the correlation coefficient as:
	  //    double corcoef = (Exy - Ex*Ey)/sqrt( (Exx_ExEx)*(Eyy-EyEy) );
	  // we compute the square of the above since that's what we need.
	  const double num = (Exy - Ex*Ey);
	  const double den = ((Exx_ExEx)*(Eyy-EyEy));
	  // If the denimonator is too small, we assume this is
	  // from a degenerate distribution for which I(X;Y) == 0.
	  if (fabs(den) < DBL_MIN) {
	    *ftr_mi_p = 0.0;
	  } else {
	    const double sq_corcoef = (num*num)/den;
	    if (sq_corcoef >= 1.0)
	      *ftr_mi_p = HUGE_VAL;
	    else
	      *ftr_mi_p = half_inv_log2*log(1-sq_corcoef);
	  }
	}

	ftr_mi_p++;
	ftr_sumsxy_p++;
	ftr_sums_l0_p++;
	ftr_sumssq_l0_p++;
	ftr_sqsums_l0_p++;
      }
      
      ftr_sums_p++;
      ftr_sumssq_p++;
      ftr_sqsums_p++;
    }
  }
}


void
Cor_Col::printMI(Range& cfr_rng,
		 FILE *out_fp,
		 const bool printEntropies)
{

  // MI array is as follows:
  // Assume X[i,l] is the i'th feature at lag l (zero offset for both)
  // and assume that L is the maximum lag, then
  // and assume that N is the number of features (i.e., N = n_ftrs), then
  // mi[0] = I(X[0,0],X[0,-L])
  // mi[1] = I(X[1,0],X[0,-L])    
  // ...
  // mi[N] = I(X[0,0],X[1,-L])
  // mi[N+1] = I(X[1,0],X[1,-L])
  // ...
  // mi[N*N-1] = I(X[N-1,0],X[N-1,-L])
  // mi[N*N] = I(X[0,0],X[0,-(L-1)])
  // ...
  // ...
  // mi[N*N*L-1] = I(X[N-1,0],X[N-1,0])
  // symmetric case for positive lags.

  double *ftr_mi_p = ftr_mi;
  for (Range::iterator cnit=cfr_rng.begin();
       !cnit.at_end();cnit++) {
    const int l = (*cnit);
    // l determines the lag
    for (int j=0;j<n_ftrs;j++) {
      // j determines the feature index
      // of X(t-l) where l=lag
      if (l!=0) {
	for (int k=0;k<n_ftrs;k++) {
	  // k determines the feature index of X(t)
	  /*fprintf(out_fp,
		  "%d:%d(%d) %.10e\n",
		  k,j,l,*ftr_mi_p);*/
	  fprintf(out_fp,
		  "[%d@0] [%d@%d] %.10e\n",
		  //"[%d@0]\t[%d@%d]\t%.10e\n",
		  k,j,l,*ftr_mi_p);
	  ftr_mi_p++;
	}
      } else {
	for (int k=(j);k<n_ftrs;k++) {
	  // k determines the feature index of X(t)
	  if (printEntropies || (k>j)) {
	    /*fprintf(out_fp,
		    "%d:%d(%d) %.10e\n",
		    k,j,l,*ftr_mi_p); */
	    fprintf(out_fp,
		    "[%d@0] [%d@%d] %.10e\n",
		    //"[%d@0]\t[%d@%d]\t%.10e\n",
		    k,j,l,*ftr_mi_p);
	  }
	  ftr_mi_p++;
	}
      }
    }
  }
}


// ======================================================================
// ======================================================================
// ----------  Main code for Linear  MI computation ---------------------
// ======================================================================
// ======================================================================

static void
pfile_mi_cc(FILE *out_fp,
	    Range& srrng,
	    Range& lrrng,
	    Range& cfr_rng,
	    const bool printEntropies,
	    const int min_labels_per_sentence,
	    int labpos)

{

    // Feature and label buffers are dynamically grown as needed.
    size_t buf_size = 1000;      // Start with storage for 300 frames.
    size_t n_labs;
    n_labs = (gotLabelFile?globalLabelMatrix.numDiscrete():globalObservationMatrix.numDiscrete());
    const size_t n_ftrs = globalObservationMatrix.numContinuous();

    // true if we do not ignore the labels to produce MI
    // that is conditional on the state.
    const bool stateCondMI = !lrrng.full();

    Cor_Col ccc(n_ftrs,cfr_rng);

    if (stateCondMI && n_labs < 1) { 
      error("For conditional MI, number of discrete features per frame in pfile must be at least one.");
    }

    size_t total_windows = 0;
    printf("Computing correlations\n");
    
    unsigned* lab_buf;
    ObservationMatrix* obsMat;
    if(gotLabelFile) {
      obsMat = &globalLabelMatrix;
      if(n_labs == 1) {
	lab_buf = (unsigned*) obsMat->features.ptr;
      }
      else
	lab_buf = new unsigned[buf_size * n_labs];
    }
    else {
      obsMat = &globalObservationMatrix;
      lab_buf = new unsigned[buf_size * n_labs];
    }
    
    // Go through input pfile to get the initial statistics,
    for (Range::iterator srit=srrng.begin();!srit.at_end();srit++) {

      obsMat->loadSegment((const unsigned)(*srit));
      const size_t n_frames = obsMat->numFrames();

      if ((*srit) % 100 == 0)
	printf("Processing sentence %d\n",(*srit));

      if (n_frames == SIZET_BAD)
	{
	  fprintf(stderr,
		  "%s: Couldn't find number of frames "
		  "at sentence %lu in input pfile.\n",
		  program_name, (unsigned long) (*srit));
	  error("Aborting.");
	}

      
	if (stateCondMI) {
	  //if we have a label file with one label per frame
	  //we just use obsMat->features.ptr as the label buffer
	  //Otherwise, we read copy the relevant label into
	  //a buffer
	  if(gotLabelFile && n_labs==1) {
	    lab_buf = (unsigned*) obsMat->features.ptr;
	  }
	  else { 
	    //printf("Building label buffer.\n");
	    if(n_frames > buf_size) {
	      printf("Buffer size is too small.  Increasing its size...\n");
	      delete [] lab_buf;
	      buf_size = 2*n_frames;
	      lab_buf = new unsigned[buf_size*n_labs];
	    }

	    unsigned label;
	    unsigned* lab_buf_ptr = lab_buf;
	    int pos = (unsigned) labpos;
	    if(pos == -1) //use default: last discrete feature
	      pos =  obsMat->numFeatures() - 1;
	    for(unsigned f = 0; f < n_frames; ++f) {
	      label =  obsMat->unsignedAtFrame(f,(const unsigned) pos);
	      *lab_buf_ptr++ = label;
	    }
	  }
	  const size_t n_read = n_frames;
	  
	  // check if any of the labels in lrrng are in
	  // the lab_buf. If there aren't any, don't bother to
	  // read the features, and just continue.
	  
	  if (min_labels_per_sentence > 0) {
	    int num_found=0;
	    for (size_t li=0;li<n_read;li++) {
	      if (lrrng.contains(lab_buf[li]))
		num_found++;
	      if (num_found >= min_labels_per_sentence)
		break;
	    }
	    if (num_found < min_labels_per_sentence)
	      goto continue_with_next_sentence;
	  }

	  globalObservationMatrix.loadSegment((const unsigned)(*srit));
      }

#if 0
	{
	  printf("Printing first 10\n");
	  for (int i =0;i<10;i++) {
	    printf("%f ",ftr_buf[i]);
	  }
	  printf("\n");
	}
#endif
      
      { // make the above 'goto' legal.
	float *ftr_buf_cur = (float*) globalObservationMatrix.features.ptr;
	unsigned *lab_buf_cur = lab_buf;
	unsigned *lab_buf_cur_endp = lab_buf + n_frames;
	size_t num_windows=0;
	for (size_t i=0;i<n_frames;i++) {
	  if (!stateCondMI || lrrng.contains(*lab_buf_cur)) {
	    num_windows++;
	    ccc.accumulateFrame(ftr_buf_cur,
				lab_buf_cur-lab_buf,
				lab_buf_cur_endp-lab_buf_cur-1,
				cfr_rng);
	  }
	  ftr_buf_cur += n_ftrs;
	  lab_buf_cur ++;
	}
	total_windows += num_windows;
      }

    continue_with_next_sentence:
      ;
    }

    if (total_windows == 0)
      error("ERROR: No frames have a label that matches lable range (-lr) criterion.");

    printf("Computing Mutual Information\n");
    ccc.normalize(cfr_rng);
    ccc.computeMI(cfr_rng);
    printf("..done\n");
    ccc.printMI(cfr_rng,out_fp,printEntropies);
    printf("Exiting...\n");
    if(!(gotLabelFile && n_labs == 1))
      delete [] lab_buf;
}



// ======================================================================
// ======================================================================
// -------------- Miscelaneous Support Routines -------------------------
// ======================================================================
// ======================================================================


#if 0
//
// Return the file size without moving the file
// pointer.
long unsigned
fsize(FILE*stream) 
{
  long curpos = ftell(stream);
  // rewind
  (void) fseek (stream, 0L, SEEK_END);
  long filesize = ftell(stream);
  (void) fseek (stream, curpos, SEEK_SET);
  return filesize;
}
#endif

//
// transposeAndCovertToDouble:
//   assumes m > 0 and n > 0
// transposes 'in' (a mXn float matrix) converting
// to doubles, to 'out' (a mXn double matrix).
void
transposeAndConvertToDouble(const float *const in,
			   const size_t m, // rows of 'in'
			   const size_t n, // cols of 'in'
			   double *const out)
{
  const float *inp = in;
  double *outp = out;
  size_t i;
  i = 0;
  do {
    double *outpp = outp;
    double *const outp_endp = outp + m;
    const float *inpp = inp;
    do {
      *outpp = *inpp;
      inpp += n;
      outpp++;
    } while (outpp != outp_endp);
    inp++;
    outp = outp_endp;
    i++;
  } while (i<n);
}


// ======================================================================
// ======================================================================
// -------------- Code for Gaussian Mixture  MI computation -------------
// ======================================================================
// ======================================================================


// ======================================================================
// ----------  Mixture Bi-variate Normal Collection Class ---------------
// ======================================================================


// A Collection of MixBiNormal objects.
class MBN_Col {

  const int n_mis;
  MixBiNormal *const ftr_mi;
  MixBiNormal *const ftr_mi_endp;

public:
  MBN_Col(const int n_mis_a)
    : n_mis(n_mis_a), ftr_mi(new MixBiNormal[n_mis]),
      ftr_mi_endp(ftr_mi+n_mis) {}
  ~MBN_Col() { delete [] ftr_mi; }

  static void setNumComps(int n_comps)
  { MixBiNormal::GlobalNumComps = n_comps; }
  static void setVarianceFloor(double vf)
  { MixBiNormal::varianceFloor = vf; }
  static void setDetFloor(double vf)
  { MixBiNormal::detFloor = vf; }
  static void setMixtureCoeffVanishNumerator(double vt)
    { MixBiNormal::mixtureCoeffVanishNumerator = vt; }
  static void setReRandsPerMixCompRedux(int re_rands)
  { MixBiNormal::reRandsPerMixCompRedux = re_rands; }
  static void setReRandomizeOnlyOneComp(bool b)
    { MixBiNormal::reRandomizeOnlyOneComp = b; }
  static void setNoReRandomizeOnDrop(bool b)
    { MixBiNormal::noReRandomizeOnDrop = b; }



  typedef void (MixBiNormal::*ATPFT)(double *,double *, const int);


  int readCurParams(FILE *const fp,const bool force_all_active=false);
  int readCurParams(FILE *const fp,Range& force_active_rng);

  void setAllToDirty();

  bool prepareCurrent();

  void writeCurParams(FILE *const fp);

  void startEpoch();
  bool addtoEpoch(double *const dftr_buf,
		  const size_t n_ftrs, // = num rows of dftr_buf
		  const size_t n_samps, // = num cols of dftr_buf
		  const size_t col_stride, // column stride
		  const size_t n_before, // additional cols before b[0]
		  const size_t n_after,  // additional cols after b[n_samps-1]
		  Range& cfr_rng,
		  Range& bfr_rng,
		  Range& lfr_rng,
		  size_t *const counts);
  void endEpoch();


  int recomputeNumActive(double &max_dist, double &avg_dist, 
			 double &min_dist,
			 const double term_dist);
  void computeMI(FILE *out_fp,const size_t n_ftrs,
		 Range& cfr_rng,
		 Range& bfr_rng,
		 Range& lfr_rng,
		 const bool grid,
		 const int n_bins,const double n_stds,
		 const bool printEntropies);

};


int
MBN_Col::readCurParams(FILE *const fp,const bool force_all_active)
{
  MixBiNormal *ftr_mi_p;    
  int numActive = n_mis;
  
  if (fp == NULL)
    return 0;

  ftr_mi_p = ftr_mi;
  while (ftr_mi_p != ftr_mi_endp) {
    ftr_mi_p->readCurParamsBin(fp);
    // since we just read it, it's not dirty since it lives on disk.
    ftr_mi_p->reSetDirty();
    ftr_mi_p++;
  }

  // At this point, we assume everything is already active.

  // read the active status as well, if it exists.
  char active;
  if (!force_all_active && (fread(&active,sizeof(char),1,fp) == 1)) {
    // then we presume there is status data.
    numActive = 0;
    ftr_mi_p = ftr_mi;
    if (active) {
      ftr_mi_p->setActive();
      numActive++;
    } else 
      ftr_mi_p->reSetActive();
    ftr_mi_p++;
    while (ftr_mi_p != ftr_mi_endp) {
      if (fread(&active,sizeof(char),1,fp) != 1)
	error("EOF encountered in parameter input file while reading active bits.");
      if (active) {
	ftr_mi_p->setActive();
	numActive++;
      } else 
	ftr_mi_p->reSetActive();
      ftr_mi_p++;
    }
  }
  return numActive;
}


//
// TODO: This routine (and the -ocaar option) is getting pretty hacky. i.e., I needed 
// a way to make a copy of a .mg file setting some (or all) of the components
// to be active. Ideally, I should write a different program to
// deal with this directly and more cleanly.
int
MBN_Col::readCurParams(FILE *const fp,Range& force_active_rng)
{
  MixBiNormal *ftr_mi_p;    
  int numActive = n_mis;
  
  if (fp == NULL)
    return 0;

  ftr_mi_p = ftr_mi;
  while (ftr_mi_p != ftr_mi_endp) {
    ftr_mi_p->readCurParamsBin(fp);
    // since we just read it, it's not dirty since it lives on disk.
    ftr_mi_p->reSetDirty();
    ftr_mi_p++;
  }

  // read the active status as well, if it exists.
  char active;
  if ((fread(&active,sizeof(char),1,fp) == 1)) {
    // then we presume there is status data.
    numActive = 0;
    int mgno=0;
    ftr_mi_p = ftr_mi;
    if (active || force_active_rng.contains(mgno)) {
      ftr_mi_p->setActive();
      numActive++;
    } else 
      ftr_mi_p->reSetActive();
    ftr_mi_p++;
    mgno++;
    while (ftr_mi_p != ftr_mi_endp) {
      if (fread(&active,sizeof(char),1,fp) != 1)
	error("EOF encountered in parameter input file while reading active bits.");
      if (active || force_active_rng.contains(mgno)) {
	ftr_mi_p->setActive();
	numActive++;
      } else 
	ftr_mi_p->reSetActive();
      ftr_mi_p++;
      mgno++;
    }
  }
  return numActive;
}


void
MBN_Col::setAllToDirty()
{
  MixBiNormal* ftr_mi_p = ftr_mi;
  while (ftr_mi_p != ftr_mi_endp) {
    ftr_mi_p->setDirty();
    ftr_mi_p++;
  }
}

bool
MBN_Col::prepareCurrent()
{
  bool rc = false;
  MixBiNormal* ftr_mi_p = ftr_mi;
  while (ftr_mi_p != ftr_mi_endp) {
    rc |= ftr_mi_p->prepareCurrent();
    ftr_mi_p++;
  }
  return rc;
}


void
MBN_Col::writeCurParams(FILE *const fp)
{
  MixBiNormal *ftr_mi_p;    

  if (fseek (fp, 0L, SEEK_SET) != 0)
    error("Error seeking to beginning of output parameter file.");
  ftr_mi_p = ftr_mi;
  while (ftr_mi_p != ftr_mi_endp) {
    if (ftr_mi_p->dirty()) {
      ftr_mi_p->printCurParamsBin(fp);
      ftr_mi_p->reSetDirty();
    } else { // assume already saved.
      ftr_mi_p->seekOverCurParamsBin(fp);
    }
    ftr_mi_p++;
  }
  // save the active status as well
  ftr_mi_p = ftr_mi;
  while (ftr_mi_p != ftr_mi_endp) {
    const char act=1; const char inact=0;
    if (ftr_mi_p->active()) {
      if (fwrite(&act,sizeof(char),1,fp) != 1)
	error("Error writing active status.");
    } else {
      if (fwrite(&inact,sizeof(char),1,fp) != 1)
	error("Error writing active status.");
    }
    ftr_mi_p++;
  }
  if (fflush(fp) != 0)
    error("Error flushing mg parameters.");
}




void
MBN_Col::startEpoch()
{
  MixBiNormal *ftr_mi_p = ftr_mi;
  while (ftr_mi_p != ftr_mi_endp) {
    if (ftr_mi_p->active())
      ftr_mi_p->startEpoch();
    ftr_mi_p++;
  }
}


bool
MBN_Col::addtoEpoch(double *const dftr_buf,
		    const size_t n_ftrs, // = num rows of dftr_buf
		    const size_t  n_samps, // = num cols of dftr_buf
		    const size_t col_stride, // column stride
		    const size_t n_before, // additional cols before b[0]
		    const size_t n_after,  // additional cols after b[n_samps-1]
		    Range& cfr_rng,
		    Range& bfr_rng,
		    Range& lfr_rng,
		    size_t *const counts)
{
  bool data_used = false;
  MixBiNormal *ftr_mi_p = ftr_mi;

  size_t *counts_p = counts;
  for (Range::iterator cnit=cfr_rng.begin();
       !cnit.at_end();cnit++,counts_p++) {
    const int offset = (*cnit);

    double*dftr_buf_cur_p;
    double*dftr_buf_lag_p;
    int num_samps;
    
    if (offset <= 0) {
      // negative lag case.
      if (-offset <= (int)n_before) {
	dftr_buf_cur_p = dftr_buf;
	dftr_buf_lag_p = dftr_buf + offset;
	num_samps = n_samps;
      } else {
	const int diff = -offset - n_before;
	dftr_buf_cur_p = dftr_buf + diff;	
	dftr_buf_lag_p = dftr_buf + offset + diff;
	num_samps = n_samps - diff;
      }
    } else {
      // positive lag case.
      if (offset <= (int)n_after) {
	dftr_buf_cur_p = dftr_buf;
	dftr_buf_lag_p = dftr_buf + offset;
	num_samps = n_samps;	
      } else {
	const int diff = offset - n_after;
	dftr_buf_cur_p = dftr_buf;
	dftr_buf_lag_p = dftr_buf + offset;
	num_samps = n_samps-diff;
      }
    }
    if (num_samps <= 0)
      continue;

    if (counts != NULL) 
      *counts_p += num_samps;
    data_used = true;

    if (bfr_rng.full() && lfr_rng.full()) {
      // then we need not worry about bfr_rng and lfr_rng

      if (offset == 0) {
	// this is the zero-lag case, only do upper off-diagonal.
	size_t j,k;
	for (j=0;j<n_ftrs;j++) {
	  double *dftr_buf_cur_pp = dftr_buf_cur_p + (j+1)*col_stride;
	  for (k=(j+1);k<n_ftrs;k++) {
	    if (ftr_mi_p->active())
	      ftr_mi_p->addtoEpoch(dftr_buf_lag_p,
				   dftr_buf_cur_pp,
				   num_samps);
	    dftr_buf_cur_pp += col_stride;
	    ftr_mi_p++;		  
	  }
	  dftr_buf_lag_p += col_stride;
	}
      } else {
	// non-zero lag case, do full matrix
	size_t j,k;
	for (j=0;j<n_ftrs;j++) {
	  double *dftr_buf_cur_pp = dftr_buf_cur_p;
	  for (k=0;k<n_ftrs;k++) {
	    if (ftr_mi_p->active())
	      ftr_mi_p->addtoEpoch(dftr_buf_lag_p,
				   dftr_buf_cur_pp,
				   num_samps);
	    dftr_buf_cur_pp += col_stride;
	    ftr_mi_p++;		  
	  }
	  dftr_buf_lag_p += col_stride;
	}
      }
    } else {
      // do both zero and non-zero lag case together.
      for (Range::iterator lfit=lfr_rng.begin();
	   !lfit.at_end();lfit++) {    
	const int j = (*lfit);
	// j determines the feature index
	// of X(t-l) where l=lag

	double *dftr_buf_lag_pp = 
	  dftr_buf_lag_p + j*col_stride;
	
	if (offset != 0) {
	  for (Range::iterator bfit=bfr_rng.begin();
	       !bfit.at_end();
	       bfit++) {
	    if (ftr_mi_p->active()) {
	      const int k = (*bfit);
	      double *dftr_buf_cur_pp = 
		dftr_buf_cur_p + k*col_stride;
	      ftr_mi_p->addtoEpoch(dftr_buf_lag_pp,
				   dftr_buf_cur_pp,
				   num_samps);
	    }
	    ftr_mi_p++;		  
	  }
	} else {
	  for (Range::iterator bfit=bfr_rng.begin();
	       !bfit.at_end();
	       bfit++) {
	    if (((*bfit) == (*lfit)) ||
		 (((*bfit) < (*lfit))
		  &&
		  (lfr_rng.contains(*bfit) && bfr_rng.contains(*lfit))))
	      continue;
	    if (ftr_mi_p->active()) {
	      const int k = (*bfit);
	      double *dftr_buf_cur_pp = 
		dftr_buf_cur_p + k*col_stride;
	      ftr_mi_p->addtoEpoch(dftr_buf_lag_pp,
				   dftr_buf_cur_pp,
				   num_samps);
	    }
	    ftr_mi_p++;		  
	  }
	}
      }
    }
  }
  return data_used;
}


void
MBN_Col::endEpoch()
{
  MixBiNormal *ftr_mi_p = ftr_mi;    
  int mc_drop=0;
  while (ftr_mi_p != ftr_mi_endp) {
    if (ftr_mi_p->active()) {
      if (ftr_mi_p->endEpoch())
	mc_drop++;
      // current parameters have not been saved.
      ftr_mi_p->setDirty();
    }
    ftr_mi_p++;
  }
  if (mc_drop > 0) {
    printf("Warning: %u mixture component drop(s)\n",mc_drop);
    fflush(stdout);      
  }
}



int 
MBN_Col::recomputeNumActive(double &max_dist, double &avg_dist, 
			    double &min_dist,
			    const double term_dist)
{
  MixBiNormal *ftr_mi_p = ftr_mi;
  max_dist = 0.0;
  avg_dist = 0.0;
  min_dist = HUGE;
  int numSum = 0;
  int numActive=0;
  while (ftr_mi_p != ftr_mi_endp) {
    if (ftr_mi_p->active()) {
      const double tmp = ftr_mi_p->llPercDiff();
      if (tmp < term_dist)
	ftr_mi_p->reSetActive();
      else
	numActive++;
      if (tmp > max_dist)
	max_dist = tmp;
      if (tmp < min_dist)
	min_dist = tmp;
      avg_dist += tmp;
      numSum++;
    }
    ftr_mi_p++;
  }
  avg_dist /= numSum;
  return numActive;
}

void	  
MBN_Col::computeMI(FILE *out_fp,
		   const size_t n_ftrs,
		   Range& cfr_rng,
		   Range& bfr_rng,
		   Range& lfr_rng,
		   const bool grid,
		   const int n_bins,
		   const double n_stds,
		   const bool printEntropies)
{
  // MI array is arranged as follows:
  // Assume X[i,l] is the i'th feature at lag l (zero offset for both)
  // and assume that L is the maximum lag in cfr_rng, then
  // and assume that N is the number of features (i.e., N = n_ftrs), then
  // mi[0] = I(X[0,0],X[0,-L])
  // mi[1] = I(X[1,0],X[0,-L])    
  // ...
  // mi[N] = I(X[0,0],X[1,-L])
  // mi[N+1] = I(X[1,0],X[1,-L])
  // ...
  // mi[N*N-1] = I(X[N-1,0],X[N-1,-L])
  // mi[N*N] = I(X[0,0],X[0,-(L-1)])
  // ...
  // ...
  // mi[N*N*(L-1)-1] = I(X[N-1,0],X[N-1,-1])
  // # here come the off-diagonal terms for the lag=0 case.
  // mi[N*N*(L-1)] = I(X[1,0],X[0,0])
  // mi[N*N*(L-1)+1] = I(X[2,0],X[0,0])
  // ...
  // mi[N*N*(L-1)+N-2] = I(X[N-1,0],X[0,0])
  // mi[N*N*(L-1)+N-1] = I(X[2,0],X[1,0])
  // mi[N*N*(L-1)+N] =   I(X[3,0],X[1,0])
  // ...
  // mi[N*N*(L-1)+N(N-1)/2-1] = I(X[N-1,0],X[N-2,0])
  // ... etc. symmetrically for positive lags.

  double *margx = new double[n_bins];
  double *margy = new double[n_bins];
  double *distxy = NULL;
  if (grid)
    distxy = new double[n_bins*n_bins];
  MixBiNormal *ftr_mi_p = ftr_mi;

  for (Range::iterator cnit=cfr_rng.begin();
       !cnit.at_end();cnit++) {
    const int l = (*cnit);
    // l is  the lag

    for (Range::iterator lfit=lfr_rng.begin();
	 !lfit.at_end();lfit++) {    
      const int j = (*lfit);

      // j determines the feature index
      // of X(t-l) where l=lag

      // for the zero-lag case (i.e., l == 0), 
      // there is a symmetry in the matrix.
      // There are three parts of the matrix
      // where bfr = base feature index, 
      //       lfr = lagged feature index.
      // 1) bfr > lfr
      // 2) bfr == lfr
      // 3) lfr > bfr
      // Part 1 and 3 might have redundancies and
      // part 2 is the entropy which we
      // don't need. Therefore, we skip
      // part 2, and skip any
      // of lfr > bfr (part 3) that
      // is redundant with (contained in) part 1.
      //
      //       +-----------------+
      //       |             / / |
      //       | lfr>bfr    / /  |
      //       |region 3   / /   |
      //       |          / /    |
      //       |         / /     |
      //       |        / /      |
      //    lfr|       / /       |
      //       |      / /        |
      //       |     / /         |
      //       |    / /          |
      //       |   / /           |
      //       |  / /  bfr>lfr   |
      //       | / /   region 1  |
      //       |/ /              |
      //       | /               |
      //       +-----------------+
      //               bfr
      //
      // I.e., for non-zero lag case, do all regions.
      // For zero-lag case, only do the portion of
      // region 1 determined by lfr_rng and bfr_rng
      // but also do any of region 3 that is specified
      // by lfr_rng and bfr_rng but not contained 
      // in the symmetrically opposite point in region 1.
      //


      for (Range::iterator bfit=bfr_rng.begin();
	   !bfit.at_end();bfit++) {
	if ( (l == 0)
	      &&
	      (((*bfit) == (*lfit)) ||
	       (((*bfit) < (*lfit))
		&&
		(lfr_rng.contains(*bfit) && bfr_rng.contains(*lfit))))
	     )
	  continue;

	const int k = (*bfit);
	// k determines the feature index of X(t)

	// We assume that ftr_mi_p->prepareCurrent() has
	// already been called.
	double mi_val,Hx,Hy,Hxy;
	if (grid)
	  mi_val = ftr_mi_p->mi(n_bins,n_stds,margx,margy,distxy,
				&Hx,&Hy,&Hxy);
	else 
	  mi_val = ftr_mi_p->mi(n_bins,margx,margy,
				&Hx,&Hy,&Hxy);

	// Print ouput  format:
	// "base_feature_index:lagged_feature_index(lag_index) MI-value"
	/*fprintf(out_fp,
		"%d:%d(%d) %.10e\n",
		k,j,l,mi_val); */
	fprintf(out_fp,
		"[%d@0] [%d@%d] %.10e\n",
		k,j,l,mi_val);

	if (printEntropies && l==0) {
	  // print out the entropy also.
	  // Hx contains the entropy at position at X(t-l)(j)
	  // Hy contains the entropy at position X(t)(k)

	  // **** BUG *****: this doesn't work right when n_ftrs == 1
	  // because this bit of code will never execute and 
	  // we'll never print out the entropy at 0:0(0).
	  // ignore for now.
	  if (k==1) {  
	    // use X(t-l)(j) == X(t)(0) since
	    //    l == 0 and when k == 1, j == 0
	    fprintf(out_fp,
		    "%d:%d(%d) %.10e\n",
		    0,0,0,Hx);
	  }
	  if (k == (j+1)) {
	    // use X(t)(k) 
	    fprintf(out_fp,
		    "%d:%d(%d) %.10e\n",
		    k,k,0,Hy);
	  }
	}
	ftr_mi_p++;
      }
    }
    if (fflush(out_fp) != 0)
      error("Error flushing MI output, lag %d.",l);
  }
  if (fflush(out_fp) != 0)
    error("Error flushing MI output");

  delete margx;
  delete margy;
  if (grid)
    delete distxy;
}

// ======================================================================
// ======================================================================
// ----------  Main code Gaussian Mixture MI computation ----------------
// ======================================================================
// ======================================================================


int compute_n_mis(Range& cfr_rng, // context frame range
		  Range& bfr_rng, // base feature element range
		  Range& lfr_rng) // lag feature element range
{
  int n_mis = cfr_rng.length()*bfr_rng.length()*lfr_rng.length();
  int subtract_off = 0;
  if (cfr_rng.contains(0)) {
    for (Range::iterator lfit=lfr_rng.begin();
	 !lfit.at_end();lfit++) {    
      for (Range::iterator bfit=bfr_rng.begin();
	   !bfit.at_end();bfit++) {    
	if (
	    ((*bfit) == (*lfit)) ||
	    (((*bfit) < (*lfit))
	     &&
	     (lfr_rng.contains(*bfit) && bfr_rng.contains(*lfit))))
	  subtract_off ++;
      }
    }
  }
  // subtract_off (== n_ftrs*(n_ftrs+1)/2 if bfr and lfr == all)
  // removed since zero-lag the zero-lag
  // case is symetric and its diagonal values aren't used.
  n_mis -= subtract_off;
  return n_mis;
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


//
//  pfile_mi_mg: Compute the mutual information (MI) using
//  EM according essentially to the following strategy:
// 
// 0. - If previous parameters exist,
//         read them as the current parameters, set clean.
//      else
//         randomize current parameters, set dirty.
//    - Precompute current
//    - ll = log(~0) (force previous parameters to be "very unlikely")
// 1. Start epoch (we assume the current parameters are legal and prepared)
//     a. zero the new accumulators
///    b. prev_ll = ll; ll = 0;
// 2. accumulate data
// 3. Update new means and variances.
//     a. update
//     b. check for any new alpha < thres
//          and if so, set up to cause a re-rand later on.
// 4. Swap new with current
// 5. - prepare current
//       check if invalid, if so
//          a. re-rand and potentially drop a mix. component
//          b. If we re-rand or drop, we do
//                 prev_ll = 2*log(~0), ll = log(~0)
//                 (this makes sure diff is big but that ll != 0)
//             but if we try to drop with only one comp left,
//             wet set to a special MI=0 and set prev_ll=ll<0 to end.
//    - set dirty
// 6. check param diff, set to inactive if diff < threshold
//     diff is defined as 100(ll-prev_ll)/abs(ll)
// 7. if number active changed or timeout
//       save all dirty, and set to undirty
// 8. If any active, goto 1.
// 9. Compute MI, with current valid and already "prepared" parameters.
// 
static void
pfile_mi_mg(FILE *out_fp, // where to put output MI values
	    FILE *pi_fp,  // where to get input MG params
	    FILE *po_fp,  // where to place output MG params
	    Range& srrng, // sentence range
	    Range& lrrng, // label range for conditional MI
	    Range& cfr_rng, // context frame range
	    Range& bfr_rng, // base feature element range
	    Range& lfr_rng, // lag feature element range
	    const int n_comps, // number of mixture components to use
	    const double term_dist, // terminal log likelihood diff
	    const bool grid, // use grid rather than law-of-large-nums
	    const int n_bins, // marginal grid size or number of samps
	    const double n_stds, // for grid, number of standard deviations.
	    const int num_active_to_stop, // stop when numActive reaches this.
	    bool compute_mi,
	    const int max_em_iterations,
	    const double varianceFloor,
	    const double detFloor,
	    const double mcvr,
	    const int re_rands,
	    const bool printEntropies,
	    const bool force_all_active,
	    const int nacps,
	    const int min_labels_per_sentence,
	    const bool rerandone,
	    const bool nodrprrnd,
	    int labpos)
{


    // Feature and label buffers are dynamically grown as needed.
    size_t buf_size = 300;      // Start with storage for 300 frames.

    size_t n_labs;
    n_labs = (gotLabelFile?globalLabelMatrix.numDiscrete():globalObservationMatrix.numDiscrete());
    const size_t n_ftrs = globalObservationMatrix.numContinuous();
    const size_t total_n_ftrs = globalObservationMatrix.numFeatures();

    // The number of MI values computed,
    const size_t n_mis =  compute_n_mis(cfr_rng,bfr_rng,lfr_rng);

    // true if we do not ignore the labels to produce MI
    // that is conditional on the state.
    const bool stateCondMI = !lrrng.full();

    // We only care about the number of labls if we're computing
    // conditional mutual information.
    if (stateCondMI && n_labs <1) { 
      error("Number of discrete features per frame must be at least one.");
    }

    int numActive = n_mis; // number of active pairs. 
    int prevNumActive = numActive;

    
    unsigned* lab_buf;
    ObservationMatrix* obsMat;
    if(gotLabelFile) {
      obsMat = &globalLabelMatrix;
      if(n_labs ==1) {
	lab_buf = (unsigned*) obsMat->features.ptr;
      }
      else {
	lab_buf = new unsigned[buf_size * n_labs];
	printf("Lab buffer copy.  Using -lni <pos of label>  option might avoid it\n");
      }
    }
    else {
      obsMat = &globalObservationMatrix;
      lab_buf = new unsigned[buf_size * n_labs];
    }

    sArray<double> dftr_buf;

    // The actual MI values.
    MBN_Col::setNumComps(n_comps);
    MBN_Col::setVarianceFloor(varianceFloor);
    MBN_Col::setDetFloor(detFloor);
    MBN_Col::setMixtureCoeffVanishNumerator(1.0/mcvr);
    MBN_Col::setReRandsPerMixCompRedux(re_rands);
    MBN_Col::setReRandomizeOnlyOneComp(rerandone);
    MBN_Col::setNoReRandomizeOnDrop(nodrprrnd);

    MBN_Col mbn(n_mis);


    if (pi_fp != NULL && fsize(pi_fp) > 0) {
      prevNumActive = numActive = 
	mbn.readCurParams(pi_fp,force_all_active);
      mbn.prepareCurrent();
    }

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


    // make sure they all do it the first time.
    time_t timeOfLastPrint = time(0) - NUM_SECONDS_PER_PRINT - 1;
    time_t timeOfLastSaveParams = time(0);
    time_t timeOfLastPrintNumActive = time(0) - MINTIMEPERPRINTNUMACTIVE -1;

    size_t em_iter=0;
    size_t n_obs_frames;
    size_t * const counts = new size_t [cfr_rng.length()];
     

    while (numActive > num_active_to_stop && !usr1_terminate && !usr2_terminate 
	   && em_iter < (unsigned)max_em_iterations) {

      em_iter ++;

      if ((time(0)-timeOfLastPrint) > NUM_SECONDS_PER_PRINT) {
	printf("Iter %d: Starting Iter\n",em_iter);
	fflush(stdout);
	timeOfLastPrint = time(0);
      }

      // do an iteration
      mbn.startEpoch();

      // @@@@@@@@@@
      // if requested, 
      //   load and/or accumulate accumulators.
     
      ::memset(counts,0,cfr_rng.length()*sizeof(size_t));
      // Go through input pfile to get the initial statistics,
      for (Range::iterator srit=srrng.begin();!srit.at_end();srit++) {

	obsMat->loadSegment((const unsigned)(*srit));
	const size_t n_frames = obsMat->numFrames();
	//	globalLabelMatrix.printSegmentInfo();

	if ((time(0)-timeOfLastPrint) > NUM_SECONDS_PER_PRINT) {
	    printf("Iter %d, sentence %d\n",em_iter,(*srit));
	    fflush(stdout);
	    timeOfLastPrint = time(0);
	}
	  
	if (n_frames == SIZET_BAD)
	  {
	    fprintf(stderr,
		    "%s: Couldn't find number of frames "
		    "at sentence %lu in input pfile.\n",
		    program_name, (unsigned long) (*srit));
	    error("Aborting.");
	  }

	if (stateCondMI) {
	  //if we have a label file with one label per frame
	  //we just use obsMat->features.ptr as the label buffer
	  //Otherwise, we read copy the relevant label into
	  //a buffer
	  if(gotLabelFile && n_labs==1) {
	    lab_buf = (unsigned*) obsMat->features.ptr;
	  }
	  else { 
	    //printf("Building label buffer.\n");
	    if(n_frames > buf_size) {
	      printf("Buffer size is too small.  Increasing its size...\n");
	      delete [] lab_buf;
	      buf_size = 2*n_frames;
	      lab_buf = new unsigned[buf_size*n_labs];
	    }

	    unsigned label;
	    unsigned* lab_buf_ptr = lab_buf;
	    int pos = (unsigned) labpos;
	    if(pos == -1) //use default: last discrete feature
	      pos =  obsMat->numFeatures() - 1;
	    for(unsigned f = 0; f < n_frames; ++f) {
	      label =  obsMat->unsignedAtFrame(f,(const unsigned) pos);
	      *lab_buf_ptr++ = label;
	    }
	  }
	  const size_t n_read = n_frames;
	  
	  // check if any of the labels in lrrng are in
	  // the lab_buf. If there aren't any, don't bother to
	  // read the features, and just continue.

#if 0
	  {
	    printf("Printing first 100\n");
	    for (int i =0;i<100;i++) {
	    //printf("%f ",ftr_buf[i]);
	      //printf("%f ",(float*) globalObservationMatrix.features.ptr);
	      //printf("%d ",(int) globalLabelMatrix.features.ptr[i]);
	      //printf("%d ",(int) lab_buf[i]);
	      printf("%d ",(int) globalLabelMatrix.features.ptr[i]);
	    }
	  printf("\n");
	  }
#endif
  
	  if (min_labels_per_sentence > 0) {
	    int num_found=0;
	    for (size_t li=0;li<n_read;li++) {
	      if (lrrng.contains(lab_buf[li]))
		num_found++;
	      if (num_found >= min_labels_per_sentence)
		break;
	    }
	    if (num_found < min_labels_per_sentence) {
#if 0
	      printf("No label in segment %d. Skipping...\n",(*srit));
#endif	     
	      goto continue_with_next_sentence;
	    }
	  }
	  //There are enough samples:  load the observation features
	  globalObservationMatrix.loadSegment((const unsigned)(*srit));
	  n_obs_frames = globalObservationMatrix.numFrames();
	  if(n_obs_frames != n_frames) {
	    //globalObservationMatrix.printSegmentInfo();
	    //globalLabelMatrix.printSegmentInfo();
	    error("ERROR: The number of frames in the label file (%d) does not match the number of frames in the observation file (%d) at sentence %d.  If you are using ascii files, make sure they are terminated with a new line.\n",n_frames,n_obs_frames,(int)(*srit));
	  }
	}
#if 0
	  {
	    printf("Printing first 78\n");
	    for (int i =0;i<78;i++) {
	    //printf("%f ",ftr_buf[i]);
	      //printf("%f ",(float*) globalObservationMatrix.features.ptr);
	      //printf("%d ",(int) globalLabelMatrix.features.ptr[i]);
	      printf("%d ",(int) lab_buf[i]);
	    }
	  printf("\n");
	  }
#endif
	  
	  // convert sentence to double precision and transpose matrix.
	  dftr_buf.growByNIfNeeded(2,n_frames*n_ftrs);
	  transposeAndConvertToDouble((float*) globalObservationMatrix.features.ptr,
				    n_frames,total_n_ftrs,dftr_buf.ptr);


	MixBiNormal::setScratchSize(n_frames,n_comps);

	if (!stateCondMI) {
	  mbn.addtoEpoch(dftr_buf.ptr,
			 n_ftrs,
			 n_frames,
			 n_frames,
			 0,0,
			 cfr_rng,
			 bfr_rng,
			 lfr_rng,
			 counts);
	} else {
	  unsigned *lab_buf_cur = lab_buf;
	  const unsigned *lab_buf_cur_endp = lab_buf + n_frames;
	  double *dftr_buf_cur = dftr_buf.ptr;

	  while (lab_buf_cur < lab_buf_cur_endp) {
	    size_t num_samps;
	    if (!lrrng.contains(*lab_buf_cur)) {
	      lab_buf_cur++;
	      dftr_buf_cur++;
	      continue;
	    } else {
	      unsigned *lab_buf_curp = lab_buf_cur+1;
	      while (lab_buf_curp < lab_buf_cur_endp
		     && lrrng.contains(*lab_buf_curp))
		lab_buf_curp++;
	      num_samps = lab_buf_curp-lab_buf_cur;
	    }
#if 0
	      printf("Sentence %d :  num_samps = %d\n",(int)(*srit),num_samps);
#endif
	    //A single sample is not worth adding to the epoch
	    if(num_samps > 1) {
#if 0
	      printf("Adding to epoch...\n");
#endif
	      mbn.addtoEpoch(dftr_buf_cur,
			     n_ftrs,
			     num_samps,
			     n_frames,
			     lab_buf_cur-lab_buf,
			     lab_buf_cur_endp-lab_buf_cur-num_samps,
			     cfr_rng,
			     bfr_rng,
			     lfr_rng,
			     counts);
	    
	    }
	    lab_buf_cur += num_samps;
	    dftr_buf_cur += num_samps;

	  }
	}

      continue_with_next_sentence:
	;
	// Only count these samps if at least some of the
	// data was used.
      }

      // @@@@ if requested,
      //    store accumulators
      //   goto cleanup and exit program successfully.


      size_t *counts_p = counts;
      for (Range::iterator cnit=cfr_rng.begin();
	   !cnit.at_end();
	   cnit++,counts_p++) {
	if (*counts_p == 0) {
	  error("ERROR: No matching frames at context offset %d.  One possible cause for the problem is that the wrong byte swapping is used for the label file.  Or simply the label(s) requested are not found in the sentence range.",(*cnit));
	}
      }

      mbn.endEpoch();
      if ((time(0)-timeOfLastPrint) > NUM_SECONDS_PER_PRINT) {
	printf("Iter %d: Finished Iter\n",em_iter);
	fflush(stdout);
	timeOfLastPrint = time(0);
      }

      double max_dist = 0.0;
      double avg_dist = 0.0;
      double min_dist = 0.0;

      prevNumActive = numActive;
      numActive = mbn.recomputeNumActive(max_dist,avg_dist,min_dist,term_dist);

      if ((numActive < prevNumActive) || 
	  ((time(0) - timeOfLastPrintNumActive) > MINTIMEPERPRINTNUMACTIVE)) {
	printf("Iter %d: NA=%d/%d(%.0f%%), PNA=%d, dist Max(%e) Avg(%e) Min(%e)\n",
	       em_iter,
	       numActive,n_mis,100*numActive/(double)n_mis,
	       prevNumActive,
	       max_dist,
	       avg_dist,
	       min_dist);
	fflush(stdout);
	timeOfLastPrint =
	  timeOfLastPrintNumActive = time(0);
      }

      if (usr2_terminate) {
	// If we got a sigusr2 recently, then we expect to
	// soon be killed (by pmake) but we don't want to be killed in 
        // the middle
	// of saving the parameters. So, instead of saving, we forfeit
	// the work done during this em_iter for safety's sake.
	printf("Iter %d: Not Saving Mixture Parameters Since Received SIGUSR2.\n",em_iter); fflush(stdout);
      } else if (po_fp != NULL && 
		 ((numActive+nacps <= prevNumActive) ||
		  (numActive == num_active_to_stop) ||
		  ((time(0) - timeOfLastSaveParams) > MINTIMEPERPARMSAVE)))
	{
	  // save the current parameters.
	  // printf("Iter %d: Saving Mixture Parameters.\n",em_iter); fflush(stdout);
	  mbn.writeCurParams(po_fp);
	  timeOfLastSaveParams = time(0);
	}
    }  // while (numActive > ...
    delete counts;

    if (usr2_terminate) {
      printf("Iter %d: Exiting early with failure due to received SIGUSR2\n",em_iter);
      exit (EXIT_FAILURE);
    }

    printf("Finished computing mixtures.\n");

    if (numActive == 0 && compute_mi) {
      printf("Computing Mutual Information\n");
      fflush(stdout);
      mbn.computeMI(out_fp,n_ftrs,cfr_rng,bfr_rng,lfr_rng,
		    grid,n_bins,n_stds,printEntropies);
      printf("Finished computing Mutual Information\n");
      fflush(stdout);
    } else {
      printf("Early stop at EM iter %d, Num Active=%d, not computing MI\n",
	     em_iter,numActive);
    }

    // restore signals.
    if (signal(SIGUSR1,SIG_DFL) == SIG_ERR)
      error("Can't unset SIGUSR1 signal.");
    if (signal(SIGUSR2,SIG_DFL) == SIG_ERR)
      error("Can't unset SIGUSR2 signal.");
    if (signal(SIGTERM,SIG_DFL) == SIG_ERR)
      error("Can't unset SIGTERM signal.");
    if (signal(SIGXCPU,SIG_DFL) == SIG_ERR)
      error("Can't unset SIGXCPU signal.");


    printf("Freeing memory...\n");
    // dftr_buf is an an sArrya object which is deleted in the destuctor
    if(!(gotLabelFile && n_labs == 1))
      delete [] lab_buf;
}


// ======================================================================
// ======================================================================
// ----------  Code for Gaussian Mixture parameter file copying ---------
// ======================================================================
// ======================================================================


static void
pfile_mi_mg_ocaar(FILE *out_fp, // where to put output MI values
		  FILE *pi_fp,  // where to get input MG params
		  FILE *po_fp,  // where to place output MG params
		  Range& cfr_rng, // context frame range
		  const int n_comps, // number of mixture components to use
		  const char *const ocaar)
{

  const size_t n_ftrs = globalObservationMatrix.numContinuous();
  //const size_t n_ftrs = in_streamp->n_ftrs();
  const size_t n_mis =  n_ftrs*n_ftrs*(cfr_rng.length())
    - (cfr_rng.contains(0) ? n_ftrs*(n_ftrs+1)/2 : 0);

  MBN_Col::setNumComps(n_comps);
  MBN_Col mbn(n_mis);

  Range ocaar_rng(ocaar,0,n_mis);

  if (pi_fp != NULL && fsize(pi_fp) > 0) {
    mbn.readCurParams(pi_fp,ocaar_rng);
    // coerce all to dirty so they'll write.
    mbn.setAllToDirty();
    mbn.writeCurParams(po_fp);
  }
}

// ======================================================================
// ======================================================================
// ----------  main() and command-line argument support -----------------
// ======================================================================
// ======================================================================

static long
parse_long(const char*const s)
{
    size_t len = strlen(s);
    char *ptr;
    long val;

    val = strtol(s, &ptr, 0);

    if (ptr != (s+len))
        error("Not an integer argument.");

    return val;
}

static float
parse_float(const char*const s)
{
    size_t len = strlen(s);
    char *ptr;
    double val;
    val = strtod(s, &ptr);
    if (ptr != (s+len))
        error("Not an floating point argument.");
    return val;
}


static double
parse_double(const char*const s)
{
    size_t len = strlen(s);
    char *ptr;
    double val;
    val = strtod(s, &ptr);
    if (ptr != (s+len))
        error("Not an floating point argument.");
    return val;
}


int main(int argc, const char *argv[])
{

    const char *input_fname = 0;   // Input pfile name.
    const char *input_fname2 = 0;   // 2nd input pfile name, if provided
    const char *output_fname = 0;   // Output file name.
    const char *label_fname = 0;   // Label pfile name (if any).

    const char *pi_fname = 0;   // parameter input file name
    const char *po_fname = 0;   // parameter output file name

    const char *sr_str = 0;   // sentence range string
    Range *sr_rng;

    const char *lr_str = 0;   // label range string
    Range *lr_rng;

    int labpos = -1;

    bool lswap = false;

    const char *cfr_str = 0;   // context frame range string
    Range *cfr_rng = NULL;

    const char *bfr_str = 0;   // base feature range string
    Range *bfr_rng = NULL;

    const char *lfr_str = 0;   // lagged feature range string
    Range *lfr_rng = NULL;

    int debug_level = 0;
    bool mg_mi = false; // true if we use the gaussian mixture MI measure.
    int n_comps=1;
    double term_dist=0.1;

    bool grid = true;

    int n_bins=100; // = nsamps if grid = false (i.e., law-of-large-nums = true)
    double n_stds=2;

    int num_active_to_stop = 0;

    bool compute_mi = true;
    int max_em_iterations = MAXINT;

    double varianceFloor = 0.0;
    int re_rands = 5;

    double detFloor = DBL_MIN;

    int min_labels_per_sentence = 1;

    // mixture coefficient vanishing ratio
    // i.e., we consider dropping a component
    // if it gets very improbable, i.e., 
    // if alpha[i] < (1.0/n_components)/mcvr
    double mcvr = 1e100;

    // number of active->inactive changes for a save to occur
    // (not counting any other timeouts).
    int nacps = 1;

    bool printEntropies = false;

    bool no_bak = false;

    bool force_all_active = false;
    
    // Only Copy And Activate Range
    const char *ocaar = NULL;

    // If True, only re-randomize the bad mixture component when
    // it goes awry. If false, all mixture components are re-randomized.
    bool rerandone=false;

    // If True, then do *not* re-randomize all mixtues when a component
    // drop occurs. If false (default behavior), all mixture components
    // get re-randomized when a drop occurs. This option could
    // speed convergence when there really are too many mixtures. If
    // this option is true, then the work spent training the non-dropped
    // components will not be wasted.
    bool nodrprrnd=false;

    program_name = *argv++;
    argc--;

    while (argc--)
    {
        char buf[BUFSIZ];
        const char *argp = *argv++;

        if (strcmp(argp, "-help")==0)
        {
            usage();
        }
        else if (strcmp(argp, "-i")==0)
        {
            // Input file name.
            if (argc>0)
            {
                // Next argument is input file name.
                input_fname = *argv++;
                argc--;
            }
            else
                usage("No input filename given.");
        }
        else if (strcmp(argp, "-i2")==0)
        {
            // Input file name.
            if (argc>0)
            {
                // Next argument is input file name.
                input_fname2 = *argv++;
                argc--;
            }
            else
                usage("No (2nd) input filename given.");
        } 
	else if (strcmp(argp, "-ni1")==0)
        {
	  if (argc>0) {
	    nis[0] = (int) parse_long(*argv++);
	    argc--;
	  } else
	    usage("No number of ints specified.");
        }

	else if (strcmp(argp, "-nf1")==0)
        {
	  if (argc>0) {
	    nfs[0] = (int) parse_long(*argv++);
	    argc--;
	  } else
	    usage("No number of floats specified.");
        }
	else if (strcmp(argp, "-fr1")==0)
        {
            if (argc>0)
            {
                frs[0] = *argv++;
                argc--;
            }
            else
                usage("No float range given.");
        } else if (strcmp(argp, "-ir1")==0)
        {
            if (argc>0)
            {
                irs[0] = *argv++;
                argc--;
            }
            else
                usage("No int range given for file 1.");
        }
	else if (strcmp(argp, "-fmt1")==0)
        {
            if (argc>0)
            {
                fmts[0] = *argv++;
                argc--;
            }
            else
                usage("No format given for file 1.");
        }
	else if (strcmp(argp, "-ni2")==0)
        {
	  if (argc>0) {
	    nis[1] = (int) parse_long(*argv++);
	    argc--;
	  } else
	    usage("No number of ints specified.");
        }
	
	else if (strcmp(argp, "-nf2")==0)
        {
	  if (argc>0) {
	    nfs[1] = (int) parse_long(*argv++);
	    argc--;
	  } else
	    usage("No number of floats specified.");
        } else if (strcmp(argp, "-ir2")==0)
        {
            if (argc>0)
            {
                irs[1] = *argv++;
                argc--;
            }
            else
                usage("No int range given for file 2.");
        }
	 else if (strcmp(argp, "-fr2")==0)
        {
            if (argc>0)
            {
                frs[1] = *argv++;
                argc--;
            }
            else
                usage("No float range given.");
        }
	else if (strcmp(argp, "-fmt2")==0)
        {
	  if (argc>0)
            {
	      fmts[1] = *argv++;
	      argc--;
            }
            else
	      usage("No format given for file 2.");
        }
        else if (strcmp(argp,"-iswap1") == 0)
	  {
	    iswps[0] = true;
	  }
        else if (strcmp(argp,"-iswap2") == 0)
          {
            iswps[1] = true;
	  }
        else if (strcmp(argp,"-lswap") == 0)
	  {
	    lswap = true;
	  }
        else if (strcmp(argp, "-o")==0)
        {
            // Output file name.
            if (argc>0)
            {
                // Next argument is input file name.
                output_fname = *argv++;
                argc--;
            }
            else
                usage("No output filename given.");
        }
        else if (strcmp(argp, "-pi")==0)
        {
            // Output file name.
            if (argc>0)
            {
                // Next argument is input file name.
                pi_fname = *argv++;
                argc--;
            }
            else
                usage("No pi filename given.");
        }
        else if (strcmp(argp, "-po")==0)
        {
            // Output file name.
            if (argc>0)
            {
                // Next argument is input file name.
                po_fname = *argv++;
                argc--;
            }
            else
                usage("No po filename given.");
        }
        else if (strcmp(argp, "-debug")==0)
        {
            if (argc>0)
            {
                // Next argument is debug level.
                debug_level = (int) parse_long(*argv++);
                argc--;
            }
            else
                usage("No debug level given.");
        }
        else if (strcmp(argp, "-entr")==0)
        {
	  printEntropies = true;
        }
        else if (strcmp(argp, "-cc")==0)
        {
	  mg_mi = false;
        }
        else if (strcmp(argp, "-nocmi")==0)
        {
	  compute_mi = false;
        }
        else if (strcmp(argp, "-nobak")==0)
        {
	  no_bak = true;
        }
        else if (strcmp(argp, "-rerandone")==0)
        {
	  rerandone = true;
        }
        else if (strcmp(argp, "-nodrprrnd")==0)
        {
	  nodrprrnd = true;
        }
        else if (strcmp(argp, "-frcact")==0)
        {
	  force_all_active = true;
        }
        else if (strcmp(argp, "-ocaar")==0)
        {
            if (argc>0)
            {
                ocaar = *argv++;
                argc--;
            }
            else
                usage("No -ocaar range given.");
        }
        else if (strcmp(argp, "-mg")==0)
        {
	  mg_mi = true;
	  if (argc>0) {
	      n_comps = (int) parse_long(*argv++);
	      argc--;
	  } else
	    usage("No -mg *n* t value given.");
	  if (argc>0) {
	      term_dist = parse_double(*argv++);
	      argc--;
	  } else
	    usage("No -mg n *t* value given.");
        }
        else if (strcmp(argp, "-grid")==0)
        {
	  grid = true;
	  if (argc>0) {
	      n_bins = (int) parse_long(*argv++);
	      argc--;
	  } else
	    usage("No -grid *b* n value given.");
	  if (argc>0) {
	      n_stds = parse_double(*argv++);
	      argc--;
	  } else
	    usage("No -grid b *n* value given.");
        }
        else if (strcmp(argp, "-lll")==0)
        {
	  grid=false;
	  if (argc>0) {
	      n_bins = (int) parse_long(*argv++);
	      argc--;
	  } else
	    usage("No -lll *n* value given.");
        }
        else if (strcmp(argp, "-sac")==0)
        {
	  if (argc>0) {
	    num_active_to_stop = (int) parse_long(*argv++);
	    argc--;
	  } else
	    usage("No -sac *n* value given.");
        }
        else if (strcmp(argp, "-maxi")==0)
        {
	  if (argc>0) {
	    max_em_iterations = (int) parse_long(*argv++);
	    argc--;
	  } else
	    usage("No -maxi *n* value given.");
        }
        else if (strcmp(argp, "-varfl")==0)
        {
	  if (argc>0) {
	    varianceFloor = parse_double(*argv++);
	    argc--;
	  } else
	    usage("No -varfl *f* value given.");
        }
        else if (strcmp(argp, "-detfl")==0)
        {
	  if (argc>0) {
	    detFloor = parse_double(*argv++);
	    argc--;
	  } else
	    usage("No -detfl *f* value given.");
        }
        else if (strcmp(argp, "-mcvr")==0)
        {
	  if (argc>0) {
	    mcvr = parse_double(*argv++);
	    argc--;
	  } else
	    usage("No -mcvr *f* value given.");
        }
        else if (strcmp(argp, "-rerands")==0)
        {
	  if (argc>0) {
	    re_rands = parse_long(*argv++);
	    argc--;
	  } else
	    usage("No -rerands *d* value given.");
        }
        else if (strcmp(argp, "-mlps")==0)
        {
	  if (argc>0) {
	    min_labels_per_sentence = parse_long(*argv++);
	    argc--;
	  } else
	    usage("No -mlps *d* value given.");
	  if (min_labels_per_sentence < 0)
	    error("-mlps argument must be positive");
        }
        else if (strcmp(argp, "-nacps")==0)
        {
	  if (argc>0) {
	    nacps = parse_long(*argv++);
	    argc--;
	  } else
	    usage("No -nacps *d* value given.");
        }
        else if (strcmp(argp, "-seed")==0)
        {
	  RAND rnd(true); 
	  rnd.coin(0.5); // keep compiler from complaining.
        }
        else if (strcmp(argp, "-cfr")==0)
        {
            if (argc>0)
            {
                cfr_str = *argv++;
                argc--;
            }
            else
                usage("No -cfr range given.");
        }
        else if (strcmp(argp, "-bfr")==0)
        {
            if (argc>0)
            {
                bfr_str = *argv++;
                argc--;
            }
            else
                usage("No -bfr range given.");
        }
        else if (strcmp(argp, "-lfr")==0)
        {
            if (argc>0)
            {
                lfr_str = *argv++;
                argc--;
            }
            else
                usage("No -lfr range given.");
        }
        else if (strcmp(argp, "-sr")==0)
        {
            if (argc>0)
            {
	      sr_str = *argv++;
	      argc--;
            }
            else
                usage("No range given.");
        }
        else if (strcmp(argp, "-lr")==0)
        {
            if (argc>0)
            {
	      lr_str = *argv++;
	      argc--;
            }
            else
                usage("No range given.");
        } 
        else if (strcmp(argp, "-labpos")==0)
        {
	  if (argc>0) {
	    labpos = parse_long(*argv++);
	    argc--;
	  } else
	    usage("No -labpos *d* value given.");
        }
        else if (strcmp(argp, "-lf")==0)
        {
            if (argc>0)
            {
	      label_fname = *argv++;
	      argc--;
            }
            else
                usage("No label file name given.");
        }
	else if (strcmp(argp, "-lni")==0)
        {
	  if (argc>0) {
	    lnis[0] = (int) parse_long(*argv++);
	    argc--;
	  } else
	    usage("No number of ints specified for label file (-lni).");
        }

	else if (strcmp(argp, "-lnf")==0)
        {
	  if (argc>0) {
	    lnfs[0] = (int) parse_long(*argv++);
	    argc--;
	  } else
	    usage("No number of floats specified for label file (-nlf).");
        }
	else if (strcmp(argp, "-labfr")==0)
        {
            if (argc>0)
            {
                lfrs[0] = *argv++;
                argc--;
            }
            else
                usage("No float range given for label file.");
        } else if (strcmp(argp, "-lir")==0)
        {
            if (argc>0)
            {
                lirs[0] = *argv++;
                argc--;
            }
            else
                usage("No int range given for label file.\n");
        }
	else if (strcmp(argp, "-lfmt")==0)
        {
            if (argc>0)
            {
                lfmts[0] = *argv++;
                argc--;
            }
            else
                usage("No format given for label file.");
        }
        else {
	  sprintf(buf,"Unrecognized argument (%s).",argp);
	  usage(buf);
	}
    }

    //////////////////////////////////////////////////////////////////////
    // Check all necessary arguments provided before creating objects.
    //////////////////////////////////////////////////////////////////////

    if (num_active_to_stop < 0) {
      error("-sac argument must be > 0");
    }
    if (varianceFloor < 0) {
      error("variance floor must be non-negative.");
    }
    if (detFloor < 0) {
      error("variance floor must be non-negative.");
    }
    if (mcvr <= DBL_MIN) {
      error("mcvr must be greater than %e.",DBL_MIN);
    }
    if (mcvr < 1.0) {
      fprintf(stderr,"WARNING: mcvr is set to %f < 1.0.",mcvr);
    }


    // If an output pfile name is not supplied, we just
    // compute the statistics.
    FILE *out_fp=NULL;
    if (output_fname==0 || !strcmp(output_fname,"-"))
      out_fp = stdout; // no output pfile desired.
    else {
      if ((out_fp = fopen(output_fname, "w")) == NULL) {
	error("Couldn't open output file for writing.");
      }
    }

    if (re_rands < 1) {
      error("Number of rerands must be >= 1");
    }
    if (nacps < 1) {
      error("nacps (number of active->inactive changes per save must be >= 1");
    }


    //////////////////////////////////////////////////////////////////////
    // Set IEEE FP exception masks
    //////////////////////////////////////////////////////////////////////

#if SET_IEEE_TRAPS
#ifdef HAVE_IEEE_FLAGS
    char *out;
    ieee_flags("set",
	       "exception",
	       "invalid",
	       &out);
    ieee_flags("set",
	       "exception",
	       "division",
	       &out);
#else
#  ifdef HAVE_FPSETMASK
    fpsetmask(
	      FP_X_INV      /* invalid operation exception */
	      /* | FP_X_OFL */     /* overflow exception */
	      /* | FP_X_UFL */     /* underflow exception */
	      | FP_X_DZ       /* divide-by-zero exception */
	      /* | FP_X_IMP */      /* imprecise (loss of precision) */
	      );
#  else // create a syntax error.
    //#    error No way known to trap FP exceptions
#  endif
#endif
#endif


#ifdef HAVE_NONSTANDARD_ARITHMETIC
    // this presumably sets a bit in the FPU that keeps
    // denormals from traping and being handled in sofware.
    // Instead, (again presumably), denormals are truncated to zero.
    nonstandard_arithmetic();
    // ieee_retrospective(stdout); 
#endif

    //////////////////////////////////////////////////////////////////////
    // Create objects.
    //////////////////////////////////////////////////////////////////////



    int fileno = 0;
    if (input_fname == NULL) {
      usage("No input pfile name supplied");
    }
    ofs[fileno++] = input_fname;
    if (input_fname2 != NULL) {
      ofs[fileno++] = input_fname2;
      printf("NOTE: Merging multiple pfiles frame-by-frame\n");
    }
        
    lr_rng = new Range(lr_str,0,MAX_LABEL_VAL);
    if (!lr_rng->full()) { //or better yet if lr_str != NULL
      // only bother to open this if the label range isn't full.
      if (label_fname!=NULL) {
	gotLabelFile = true;
	lofs[0] = label_fname;
	//	lnfs[0] = 0; lnis[0] = 1;
	if(lswap == true) liswps[0] = true;
	if (lofs[0] != NULL && lnis[0] == 0)
	  error("ERROR: command line must specify lni not zero");
	if(labpos < 0 && lnis[0] > 1) {
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
      if(labpos != -1 && 
	 (labpos < (int)numLabelContinuous || labpos >= (int)numLabelFeatures) )
      error("labpos (%d) out of range (%d - %d): must be within the range of discrete observations of the label file\n",labpos,numLabelContinuous,numLabelFeatures-1);
      }
     
      else {
	gotLabelFile = false;
	if(labpos < 0) { //No label position has been given
	  printf("Using the last discrete feature of the last input file as the label\n");
	}
      }
    }

    ////////////////////////////////////////////
    // check for valid argument values.
    int nfiles = 0;
    unsigned ifmts[MAX_NUM_OBS_FILES];
    for (int i=0;i<MAX_NUM_OBS_FILES;i++) {
      if (ofs[i] != NULL && nfs[i] == 0 && nis[i] == 0)
	error("ERROR: command line must specify one of nf%d and ni%d not zero",i+1,i+1);
      nfiles += (ofs[i] != NULL);
      if (strcmp(fmts[i],"htk") == 0)
	ifmts[i] = HTK;
      else if (strcmp(fmts[i],"binary") == 0)
	ifmts[i] = RAWBIN;
      else if (strcmp(fmts[i],"ascii") == 0)
	ifmts[i] = RAWASC;
      else if (strcmp(fmts[i],"pfile") == 0)
	ifmts[i] = PFILE;
      else
	error("ERROR: Unknown observation file format type: '%s'\n",fmts[i]);
    }
    
    globalObservationMatrix.openFiles(nfiles,
				      (const char**)&ofs,
				      (const char**)&frs,
				      (const char**)&irs,
				      (unsigned*)&nfs,
				      (unsigned*)&nis,
				      (unsigned*)&ifmts,
				      (bool*)&iswps);
 
    unsigned numFeatures = globalObservationMatrix.numFeatures();
    unsigned numContinuous = globalObservationMatrix.numContinuous();
    if(!gotLabelFile && labpos != -1 && 
       (labpos < (int)numContinuous || labpos >= (int)numFeatures) )
      error("labpos (%d) out of range (%d - %d): must be within the range of discrete obsevations\n",labpos,numContinuous,numFeatures-1);

    sr_rng = new Range(sr_str,0,globalObservationMatrix.numSegments());

    if (cfr_str == NULL) {
      error("Must specify a context frame range.");
    } else {
      //numFrames should be the the totalnumber of frames
      //not just for the current segment
      //cout<<"globalObservationMatrix.numFrames = "
      //	  <<globalObservationMatrix.numFrames()<<endl;
      //cfr_rng = new Range(cfr_str,
      //		  -(int)globalObservationMatrix.numFrames(),
      //		  (int)globalObservationMatrix.numFrames());
      cfr_rng = new Range(cfr_str, -MAX_CFR, MAX_CFR);
      if (cfr_rng->length()==0) {
	  error("Must specify a non-empty context frame range.");
      }
    }

    bfr_rng = new Range(bfr_str,0,numContinuous);
    lfr_rng = new Range(lfr_str,0,numContinuous);

    //////////////////////////////////////////////////////////////////////
    // Do the work.
    //////////////////////////////////////////////////////////////////////

    if (!mg_mi) {
      pfile_mi_cc(out_fp,
		  *sr_rng,
		  *lr_rng,
		  *cfr_rng,
		  printEntropies,
		  min_labels_per_sentence,
		  labpos);
    } else {
      FILE *pi_fp = NULL;
      FILE *po_fp = NULL;
      if (pi_fname!=0 && po_fname!=0 && !strcmp(pi_fname,po_fname)) {

	if (!no_bak && ((pi_fp = fopen(pi_fname, "r")) != NULL)) {
	  // backup the file to : pi_fname + ".bak"
	  const int bufsiz = strlen(pi_fname)+8096;
	  size_t sz;
	  char *buf = new char[bufsiz];
	  sprintf(buf,"%s.bak",pi_fname);
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
	} else {
	  // parameter file does not yet exist.
	}

	pi_fp = po_fp = fopen(pi_fname, "r+");
	if (pi_fp==NULL) // assume file doesn't exist.
	  pi_fp = po_fp = fopen(pi_fname, "w+");
	if (pi_fp==NULL)
	  error("Couldn't open i/o file for reading/writing.");
      } else {
	if (pi_fname!=0) {
	  pi_fp = fopen(pi_fname, "r");
	  if (pi_fp==NULL)
	    error("Couldn't open input pi file for reading.");
	}
	if (po_fname != 0) {
	  po_fp = fopen(po_fname, "w");
	  if (po_fp==NULL)
	    error("Couldn't open output po file for writing.");
	}
      }
      if (ocaar != NULL) {
	if (pi_fp == NULL || po_fp == NULL)
	  error("Error, -ocaar option given without corresponding -po and -pi options.");
	pfile_mi_mg_ocaar(out_fp,
			  pi_fp,
			  po_fp,
			  *cfr_rng,
			  n_comps,
			  ocaar);
      } else {
	pfile_mi_mg(out_fp,
		    pi_fp,
		    po_fp,
		    *sr_rng,
		    *lr_rng,
		    *cfr_rng,
		    *bfr_rng,
		    *lfr_rng,
		    n_comps,term_dist,grid,n_bins,n_stds,
		    num_active_to_stop,compute_mi,max_em_iterations,
		    varianceFloor,detFloor,mcvr,
		    re_rands,
		    printEntropies,
		    force_all_active,
		    nacps,
		    min_labels_per_sentence,
		    rerandone,nodrprrnd,
		    labpos);
      }
      if (pi_fp != NULL)
	fclose(pi_fp);
      if (po_fp != NULL && po_fp != pi_fp)
	fclose(po_fp);

    }

    //////////////////////////////////////////////////////////////////////
    // Clean up and exit.
    //////////////////////////////////////////////////////////////////////

    delete sr_rng;
    delete cfr_rng;
    if (lr_rng != NULL)
      delete lr_rng;
    if (out_fp && fclose(out_fp))
        error("Couldn't close output file.");
    //    if (lab_fp && fclose(lab_fp))
    //  error("Couldn't label input file.");
    return EXIT_SUCCESS;
}

// ======================================================================
// ---------------------------- END -------------------------------------
// ======================================================================
