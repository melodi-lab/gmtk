/**
 *: mixNormal.h
 *  This class contains the mixture of N normal distributions
 *    f(x) = \sum_{k=1}^N{alpha_k e^{-(x-u_k)^2/2/sigma^2}/(sqrt(2 pi) sigma)}
 *
 *  When using full convariance matrix, the inverse of the matrix
 *  M^(-1) is written in (1-B)'D(1-B) where * D is a diagonal matrix
 *  and B is a upper triangle matrix with diagonal zero.  * Actually
 *  in the implementation, I used B to denote the matrix 1 - B.
 */


#ifndef MIX_NORMAL_H
#define MIX_NORMAL_H


#include <cstdio>
#include <fstream>
#include <string>

#include "error.h"
#include "data-points.h"

//////////////////// class MixNormal ////////////////////

// at the first stage, we just deal with diagonal covariance

/**
 * a gaussian mixture compoent
 */
class MixNormal {

public:

/**
   * when we use the full covariance matrix,
   * we need to train one more matrix
   */
  bool _fullCoVar;

  // constructors and destructor
  MixNormal() {
    bitmask = 0; 
  }

  //////////////////// MixNormal ////////////////////

  /**
   * default constructor
   *
   * @param numVariables the number of variables
   * @param numMixures the origin number of mixtures
   *                   This might change after the learning
   * @exception OutOfMemoryError out of memory
   */
  MixNormal(unsigned index,
	    unsigned int numVariables,
	    unsigned int numMixtures,
	    unsigned maxKMeansIter,
	    bool fullCoVar=false) {
    bitmask = 0;
    setInit(index,numVariables, numMixtures, maxKMeansIter, fullCoVar,0,0,1e-10);
  } // end MixNormal

  virtual ~MixNormal();

  //_numSamples holds the number of samples associated with this pair
  unsigned getNumSamples() { return _numAccum;}
  unsigned getNumEMIterToFinish() { return _numEMIter; }
  void setNumEMIters(unsigned numEMIter) { _numEMIter = numEMIter; }

  // "active" is used to tag mixtures that have not converged yet.
  void setActive();
  void reSetActive();
  bool active() const;

  // "dirty" is used to specify whether this mixture has been saved to disk.
  void setDirty();
  void reSetDirty();
  bool dirty() const;

  // prepare some parameters for em training
  // initializes data structures
  void setInit(unsigned index, unsigned numVariables, unsigned numMixtures, unsigned maxKMeansIter, bool fullCoVar,float covAddConst,float covAddEpsilon,double clampCov);
  // set log-likelihood percent difference below which the mixture is said to have converged
  void setLLDP(double lldp);
  // return the percent difference between the current and the previous log likelihoods
  double llPercDiff() const;

  // Interface to the EM algorithm.
  void startEpoch();
  bool addToEpoch(const PointerSetToDataPoints &pointerSet);  // for the full covariance case
  bool addToEpochDiag(const PointerSetToDataPoints &pointerSet); // for the diagonal case
  bool endEpoch();

  // kmeans interface
  void startKMeansEpoch(bool randLabel, bool estCov);
  void addToKMeansEpoch(const PointerSetToDataPoints &pointerSet);
  unsigned endKMeansEpoch();
  

  // io operations
  void readParams(ifstream &ifs, bool bin = false);
  void readCurParams(FILE *fp);
  void readCurParamsBin(FILE *fp);
  void writeParams(ofstream &ofs, bool bin = false) const;
  void print();
  void printCurParams(FILE *fp) const;
  void printMeans(FILE *fp) const;
  void printCurParamsBin(FILE *fp) const;
  void seekOverCurParamsBin(FILE *fp) const;

  void notPosDef(unsigned index, unsigned compNum, unsigned dim, PARAM_DATA_TYPE* cov);
  
  // pdf calculation using the current paramters
  PARAM_DATA_TYPE prob_x(const PointerSetToDataPoints &pointerSet,unsigned sample) const;
  PARAM_DATA_TYPE prob_x(const PARAM_DATA_TYPE* vec, unsigned vecLen) const;
  PARAM_DATA_TYPE prob_x(const PARAM_DATA_TYPE* vec, unsigned vecLen,  PARAM_DATA_TYPE* meanVecs, PARAM_DATA_TYPE* B, PARAM_DATA_TYPE* varVecs, PARAM_DATA_TYPE* dets, PARAM_DATA_TYPE inv_pow_sqrt_2pi) const;

  // --------- sampling --------------------
  unsigned sampleComponent() const;
  unsigned sampleUniform(unsigned a, unsigned b) const;
  void sample(PARAM_DATA_TYPE *vecOutPtr) const;
  void sampleUsingCov(PARAM_DATA_TYPE *vecOutPtr) const;
  
  void generateData(FILE* ofp,unsigned numSamples);
  void generateDataUsingCov(FILE* ofp,unsigned numSamples);
  // -------------------------------------------

  // randomize the lth component with random dimension and parameters 
  // really shouldn't use it alone for mixture l
  void randomize_l(unsigned l, unsigned d);
  
  
  // invert a positive definite matrix
  PARAM_DATA_TYPE invertCov(PARAM_DATA_TYPE* cov, PARAM_DATA_TYPE* invCov, unsigned d) const;
  
  // partition a square matrix into blocks of matrices, split at dX
  void partitionCov(unsigned dX,
		    PARAM_DATA_TYPE* cov,
		    PARAM_DATA_TYPE* covX,
		    //PARAM_DATA_TYPE* covXY,
		    PARAM_DATA_TYPE* covY,
		    unsigned d) const;
  
  // marginalizes the first dX components of a GM
  void marginalize(unsigned dX);
  void marginalizeWithFirstParentOut(unsigned dX);
  
  // finds the probability of a sample using natural parameters
  PARAM_DATA_TYPE prob_x_Gaussian(const PointerSetToDataPoints& ps,
				  const unsigned sampleNum,
				  const PARAM_DATA_TYPE* mean, 
				  const PARAM_DATA_TYPE* invCov, 
				  const PARAM_DATA_TYPE invDet,
				  const double invSqrtTwoPiD) const;
  PARAM_DATA_TYPE prob_x_Gaussian(const PointerSetToDataPoints& ps,
				  const unsigned sampleNum,
				  const unsigned start, // start of subVector
				  const unsigned end,   // end of subVector 
				  const PARAM_DATA_TYPE* mean, 
				  const PARAM_DATA_TYPE* invCov, 
				  const PARAM_DATA_TYPE invDet,
				  const double invSqrtTwoPiD) const;
  PARAM_DATA_TYPE prob_x_Gaussian(const PARAM_DATA_TYPE *vec,
				  const unsigned vecSize,
				  const PARAM_DATA_TYPE* mean, 
				  const PARAM_DATA_TYPE* invCov, 
				  const PARAM_DATA_TYPE invDet,
				  const double invSqrtTwoPiD) const;

  // finds prob of a GM
 PARAM_DATA_TYPE prob_x_GM(const PointerSetToDataPoints& ps,
			   const unsigned sampleNum,
			   const PARAM_DATA_TYPE *alphas,
			   const PARAM_DATA_TYPE *means,  
			   const PARAM_DATA_TYPE *invCovs, 
			   const PARAM_DATA_TYPE *invDets,
			   const double invSqrtTwoPiD) const;
 PARAM_DATA_TYPE prob_x_GM(const PointerSetToDataPoints& ps,
			   const unsigned sampleNum,
			   const unsigned start, // start of subVector
			   const unsigned end,   // end of subVector 
			   const PARAM_DATA_TYPE *alphas,
			   const PARAM_DATA_TYPE *means,  
			   const PARAM_DATA_TYPE *invCovs, 
			   const PARAM_DATA_TYPE *invDets,
			   const double invSqrtTwoPiD) const;
 PARAM_DATA_TYPE prob_x_GM(const PARAM_DATA_TYPE* vec,
			   const unsigned vecSize,
			   const PARAM_DATA_TYPE *alphas,
			   const PARAM_DATA_TYPE *means,  
			   const PARAM_DATA_TYPE *invCovs, 
			   const PARAM_DATA_TYPE *invDets,
			   const double invSqrtTwoPiD) const;

 // convert moments to LDU for covariance
 bool moment2LDU(PARAM_DATA_TYPE*  cov, PARAM_DATA_TYPE* b,PARAM_DATA_TYPE* invVars, unsigned d) const;

 // convert LDU to covariance
 void LDU2moment(PARAM_DATA_TYPE* outM, PARAM_DATA_TYPE* b, PARAM_DATA_TYPE* invVars, unsigned d) const;

 
 // compute MI using all of mixtures by law of large numbers (LLN)
 PARAM_DATA_TYPE computeMI(unsigned  nSamples, unsigned  dX, PARAM_DATA_TYPE& Hx, PARAM_DATA_TYPE& Hy, PARAM_DATA_TYPE& Hxy);

 PARAM_DATA_TYPE computeMIWithFirstParentOut(unsigned nSamples, unsigned dX, PARAM_DATA_TYPE& Hx, PARAM_DATA_TYPE& Hy, PARAM_DATA_TYPE& Hxy);


 //-- compute MI using the actual data instead of sampling by LLN ---
 void startEpochMI(unsigned dX);
 void addToEpochMI(PointerSetToDataPoints& ps);
 PARAM_DATA_TYPE endEpochMI(PARAM_DATA_TYPE& Hx, PARAM_DATA_TYPE& Hy, PARAM_DATA_TYPE& Hxy);
 //-----------------------------------------------------------------

 // ----- calculate ENTROPY using sampling -----
 PARAM_DATA_TYPE computeEntropy(unsigned  nSamples, PARAM_DATA_TYPE& H);
 // ----- using data ------
 void startEpochEntropy();
 void addToEpochEntropy(PointerSetToDataPoints& ps);
 PARAM_DATA_TYPE endEpochEntropy(PARAM_DATA_TYPE& H);
 // ---------------------------------------------

 unsigned getNumberofMixtures() const;
 
 // the following should be in protected, just for debugging
 void randomize();
 
 // the following should be removed after debugging
 void printParams() const;
 
 // static parameters
 
 /**
  * if any variance is below this value, we do a re-randomization
  */
 //static ProcType varianceFloor;
 static PARAM_DATA_TYPE varianceFloor;
 
 /**
  * if the determinant of covariance matrix is below this value,
  * we do a re-randomization
  */
 //  static ProcType detFloor;
 static PARAM_DATA_TYPE detFloor;

  /**
   * if the coefficient of a component is below 1/numMixtures/mcvr,
   * we drop this component
   */
  //  static ProcType mcvr;
  static PARAM_DATA_TYPE mcvr;
  
  /**
   * the max mumber of re-randomization that occur before the
   * number of mixture is reduced
   */
  static unsigned maxNumReRands;

  /**
   * if true, re-randomize the bad mixture component when it goes
   * away. otherwise all mixture components are re-randomized.
   */
  static bool reRandOnlyOneComp;

  /**
   * if true, do not re-randomize all mixtures when a component
   * drop occurs.
   */
  static bool noReRandOnDrop;
  
  /**
   * number of re-randomizations occur before the number of mixture
   * component is reduced
   */
  static unsigned reRandsPerMixCompRedux;


  // KMEANS STUFF

  // calculate the covariance
  void calcCov(double *cov, double * mean, double* alpha);
  void calcSampleMeanCov(double* mean, double *cov);
  void kMeansNormalize();   //> normalize alphas
  unsigned nearestCluster(const PointerSetToDataPoints &pointerSet, unsigned offset) const;
  void setParameter(unsigned index,unsigned d, unsigned k, unsigned maxIter);

  void calcB();
  void checkB();

  void printCurKMeansParamsBin(FILE *fp) const;
  void readCurKMeansParamsBin(FILE *fp);
  void printCurKMeansParams(FILE *fp) const;
  void readCurKMeansParams(FILE *fp);
  
  protected:
  // randomize the parameters
  void randomize(unsigned l);

  void normalize();
  void normalizeWithNoise();
  
  // copy the paramters
  void copyNextToCur();
  bool prepareNext();

  void eliminateComponent(unsigned l);
  
  // return  Prob(l|x,theta);
  PARAM_DATA_TYPE prob_l_x_theta(unsigned l, const PointerSetToDataPoints &pointerSet,unsigned sample) const;
  // CHANGE THIS!!! ONLY DIFFERENCE IN THE CASE OF THE FIRST LETTER.
  void Prob_l_x_theta(PARAM_DATA_TYPE &Px, 
		      const PointerSetToDataPoints &pointerSet, 
		      PARAM_DATA_TYPE * postProbOutPtr,
		      unsigned sample) const;
  // return the probability from component l
  // P_l(x, theta_l)
  PARAM_DATA_TYPE prob_x_theta_l(const PointerSetToDataPoints &pointerSet, unsigned l, unsigned sample) const;
  PARAM_DATA_TYPE prob_x_theta_l(const PARAM_DATA_TYPE *vec, unsigned l, unsigned vecLen) const;
  PARAM_DATA_TYPE prob_x_theta_l(const PARAM_DATA_TYPE* vec, unsigned l, unsigned vecLen,  PARAM_DATA_TYPE* meanVecs, PARAM_DATA_TYPE* B, PARAM_DATA_TYPE* varVecs, PARAM_DATA_TYPE* dets, PARAM_DATA_TYPE inv_pow_sqrt_2pi) const;
  // the native parameters for the variables
  /** component weights */
    PARAM_DATA_TYPE *_alphas;

  /** means */
    PARAM_DATA_TYPE *_means;

  /** diagonal entries of the covariance matrix. */
    PARAM_DATA_TYPE *_invVars;

  /** inverse of determinents of covariance matrices */
  PARAM_DATA_TYPE *_invDets;

  /**
   * the matrix 1 - B in the model
   * in this case, we are using full covariance
   * matrix.  this is used only when fullCoVariance
   * is true. which means the S = (1-B)'D(1-B).
   */
  PARAM_DATA_TYPE *_b;
  
  /**
   * the minimum percentage difference
   * in log likelihood.  The em training
   * will stop at this point
   */
  double _lldp;


private:
  
  virtual inline void cleanup();
  
  /* The number of EM iterations it took this mixture to converge*/
  unsigned  _numEMIter;
  
  /** the number of variables we are talking about */
  unsigned _numVariables;

  /** the pdf has this # of mixture gaussian */
  unsigned _numMixtures;

  /**
   * the original value because we have chances
   * to reduce this number
   */
  unsigned _orgNumMixtures;

  /** the difference of log likelihood */
  double _diffLogLikelihood;

  /** the number of epoches in the training */
  unsigned _numEpoches;
  
  /** the total number of re-randomization */
  unsigned _totalCurrentReRands;

  /** the scale of variance for re-randomization */
  //ProcType _reRandScale;
  PARAM_DATA_TYPE _reRandScale;

  /** a bitmask giving the state of the object */
  enum {
    bm_active = 0x01,			// true if we have not yet converged.
    bm_dirty  = 0x02,			// true if current parameters have not been saved

    bm_initState = 0x0
  };
  
  /** save more space for flags */
  unsigned bitmask;

  // these variables are created in order to elminate the new and delete operations
  PARAM_DATA_TYPE *_nextAlphas;
  PARAM_DATA_TYPE *_nextMeans;  
  PARAM_DATA_TYPE *_nextInvVars;
  PARAM_DATA_TYPE *_nextB;
  PARAM_DATA_TYPE *_nextCov;
  
  /** the log likelihood of the convergence */
  double _logLikelihood;
  
  /** the previous log likelihood */
  double _preLogLikelihood;

  /** the counter for accumulation in em training */
  unsigned _numAccum;
  
  /** a number that is used very often */
  double _inv_pow_sqrt_2pi;
  
  
  /** MI VARIABLES */
  /** size of X, the length of first split */
  unsigned _dX;
  /** information theoretic quantities */
  PARAM_DATA_TYPE _Hx, _Hy, _Hxy;
  
  /** means for partial vectors X and Y, full vector [X' Y']' */
  PARAM_DATA_TYPE *_meansX, *_meansY;
  
  /** inverse covariances for X and Y */
  PARAM_DATA_TYPE *_invCovsX, *_invCovsY;
  // added temporarly.
  PARAM_DATA_TYPE *_invCovsXY;
  
  /** inverse determinant for X and Y */
  PARAM_DATA_TYPE *_invDetsX, *_invDetsY;  

  PARAM_DATA_TYPE *_invVarsX,* _invVarsY, *_bX,*_bY;

  // temporary matrix which stores the inverse of a given matrix
  PARAM_DATA_TYPE * _tmpM;
  PARAM_DATA_TYPE * _tmpM2;
  
  /** number of samples used in LLN calculation by using "real" data */
  unsigned _nSamplesMI;
  
  /** sqrt(2\pi)^(dX),(dY) */
  double _inv_pow_sqrt_2piX, _inv_pow_sqrt_2piY;

  // kMeans variables
  
  int _index; // holds the index number of the n-tuplet (same as the
  // order of the n-tuple in the range file)

  unsigned _d;       // dimensionality of the n-tuple (i.e. _d == n in n-tuple)
  unsigned _k;       // number of clusters to use
  unsigned _maxKMeansIter; // maximum number of kmeans iterations
  
  double _eqWeight;  // intialized to 1 / _numAccum
  
  bool _randLabel;
  bool _estCov;
  
  double *_cov;
  double * _accumMean;  // an accumulator for the vectors
  //  ----- end of kmeans variables

  float covAddConst; // adds a small amount of noise to the diagonal entries of covariance matrices to prevent numerical issues.
  float covAddEpsilon; // we add a random number between covAddConst-covAddEpsilon and covAddConst+Epsilon

  double clampCov; // Clamps the covariance to avoid sigularities when it gets too small 

}; // end class MixNormal


#endif
