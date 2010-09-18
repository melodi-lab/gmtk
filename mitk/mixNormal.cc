#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cfloat>

#if HAVE_CONFIG_H
#include <config.h>
#if HAVE_MATH_H
#include <math.h>
#endif
#if !HAVE_FINITE
#error "I need the finite function"
#endif
#endif

#include "rand.h"
#include "mixNormal.h"
#include "matrix-ops.h"

#if defined(WIN32)
#pragma warning( disable : 4800 )	      // Disable warning messages 4800
#endif

#ifdef __sparc
#include "ieeeFPsetup.h"
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_2PI
#define M_2PI 6.28318530717958647692
#endif

//extern bool Seed;  // declared and initialized in multivariate_mi.cc
//RAND rnd(Seed);  
//RAND rnd(true);  


/**
 * set up some static parameters
 */

PARAM_DATA_TYPE MixNormal::varianceFloor = 10e-20;
PARAM_DATA_TYPE MixNormal::mcvr = 40;
unsigned MixNormal::maxNumReRands = 5;
bool MixNormal::reRandOnlyOneComp = true;
bool MixNormal::noReRandOnDrop = true;
unsigned MixNormal::reRandsPerMixCompRedux = 5;


#ifdef DOUBLE_PROCESSING_DEFINED
double MixNormal::detFloor = DBL_MIN;
#else
double MixNormal::detFloor = FLT_MIN;
#endif

#ifndef LOG_ZERO
#ifdef DOUBLE_PROCESSING_DEFINED
const double LOG_ZERO = -1e250;
#else
const double LOG_ZERO = -1e-30F;
#endif
#endif

//////////////////// ~MixNormal ////////////////////
/**
 * destructor
 * cleanup the memory
 */
MixNormal::~MixNormal() {
  cleanup();
} // end MixNormal


//////////////////// setActive ////////////////////
/**
 * "active" is used to say if this mixture has converged
 */
void MixNormal::setActive() {
  bitmask |= bm_active;
} // end setActive


//////////////////// reSetActive ////////////////////
/**
 * "active" is used to say if this mixture has converged
 */
void MixNormal::reSetActive() {
  bitmask &= ~bm_active;
} // end reSetActive


//////////////////// active ////////////////////
/**
 * get the active status
 *
 * @return the active state
 */
bool MixNormal::active() const {
  return bitmask & bm_active;
} // end active


//////////////////// setDirty ////////////////////
/**
 * "dirty" is used to say if this mixture has been saved to disk
 */
void MixNormal::setDirty() {
  bitmask |= bm_dirty;
} // end setDirty


//////////////////// reSetDirty ////////////////////
/**
 * "dirty" is used to say if this mixture has been saved to disk
 */
void MixNormal::reSetDirty() {
  bitmask &= ~bm_dirty;
} // end reSetDirty


//////////////////// dirty ////////////////////
/**
 * get the dirty status
 *
 * @return the dirty status
 */
bool MixNormal::dirty() const {
  return bitmask & bm_dirty;
} // end dirty


//////////////////// setInit ////////////////////
/**
 * set the initial parameters
 *
 * @param numVariables the number of variables
 * @param numMixtures the number of mixtures
 * @param fullCoVar whether use full covariance matrix
 */
void MixNormal::setInit(unsigned index,
			unsigned numVariables,
			unsigned numMixtures,
			unsigned maxKMeansIter,
			bool fullCoVar,
			float covAddConst,
			float covAddEpsilon,
			double clampCov) {
  _index =  index;
  _numVariables = numVariables;
  _numMixtures = _orgNumMixtures = numMixtures;
  _fullCoVar = fullCoVar;
  _numEMIter = 0;

  if ( _numVariables <= 1 )
    _fullCoVar = false;		// since we only have 1-d

  _alphas      = new PARAM_DATA_TYPE [_numMixtures];
  _means       = new PARAM_DATA_TYPE [_numMixtures*_numVariables];
  _invVars     = new PARAM_DATA_TYPE [_numMixtures*_numVariables];
  _invDets     = new PARAM_DATA_TYPE [_numMixtures*_numVariables];

  _nextAlphas  = new PARAM_DATA_TYPE [_numMixtures*_numVariables];
  _nextMeans   = new PARAM_DATA_TYPE [_numMixtures*_numVariables];
  _nextInvVars = new PARAM_DATA_TYPE [_numMixtures*_numVariables];

  _meansX = new PARAM_DATA_TYPE [_numMixtures*_numVariables]; // do we
  // know dX at this point?  If so replace _numVariables with dX
  _invDetsX = new PARAM_DATA_TYPE [_numMixtures*_numVariables];
  _invCovsX = new PARAM_DATA_TYPE [_numMixtures*_numVariables*_numVariables];

  _invVarsX     = new PARAM_DATA_TYPE [_numMixtures*_numVariables];

  _meansY   = new PARAM_DATA_TYPE [_numMixtures*_numVariables];  // same comment as for meansX
  _invDetsY = new PARAM_DATA_TYPE  [_numMixtures*_numVariables];
  _invCovsY = new PARAM_DATA_TYPE [_numMixtures*_numVariables*_numVariables];

  _invVarsY     = new PARAM_DATA_TYPE [_numMixtures*_numVariables];

  _invCovsXY = new PARAM_DATA_TYPE [_numMixtures*_numVariables*_numVariables];

  _tmpM =  new PARAM_DATA_TYPE [_numVariables*_numVariables];
  _tmpM2 =  new PARAM_DATA_TYPE [_numVariables*_numVariables];


  if( _alphas == NULL || _means == NULL || _invVars == NULL 
      || _invDets == NULL || _nextAlphas == NULL 
      || _nextMeans == NULL || _nextInvVars == NULL 
      || _meansX == NULL || _invDetsX == NULL || _invCovsX == NULL 
      || _meansY == NULL || _invDetsY == NULL || _invCovsY == NULL )
    error("cannot create MixNormal");

  // for full covariance matrices
  
  if ( fullCoVar ) {
    _b = new PARAM_DATA_TYPE[_numMixtures * _numVariables * _numVariables];
    _bX = new PARAM_DATA_TYPE[_numMixtures * _numVariables * _numVariables];
    _bY = new PARAM_DATA_TYPE[_numMixtures * _numVariables * _numVariables];
    _nextB = new PARAM_DATA_TYPE [_numMixtures * _numVariables * _numVariables];
    _nextCov = new PARAM_DATA_TYPE [_numMixtures * _numVariables * _numVariables];
    if ( _b == NULL || _nextB == NULL || _nextCov == NULL ) error("cannot create MixNormal object");
  } else {
    _b = NULL;
    _nextB = NULL;
  }

  // some protected/private vaiables
  _logLikelihood = 0;
  _preLogLikelihood = 0;
  _diffLogLikelihood = DBL_MAX;
  _lldp = DBL_MIN;
  _numEpoches = 50;
  _reRandScale = 1;
  _totalCurrentReRands = 0;
  _inv_pow_sqrt_2pi = 1.0 / pow(sqrt(M_2PI), (int)_numVariables);
  //maxNumReRands = 50;

  // kMeans stuff
  _maxKMeansIter = maxKMeansIter;
  _k = numMixtures;  // redundant, have to get rid of them at some point
  _d = numVariables; // same
  _cov   = new double [_k*_d*_d];
  _accumMean = new double [_k*_d];
  _numAccum = 0;
  if(_cov == NULL) error("Could not allocate memory for MixNormal::_cov");

  this->covAddConst = covAddConst;
  this->covAddEpsilon = covAddEpsilon;

  this->clampCov = clampCov;

} // end setInit


//////////////////// setLLDP ////////////////////

/**
 * set the percentage difference in log-likelihood
 *
 * @param lldp the log likelihood difference in percentage
 */
void MixNormal::setLLDP(double lldp) {
  _lldp = lldp;
} // end setLLDP


//////////////////// llPercDiff ////////////////////

/**
 * get the percentage difference in log-likelihood
 *
 * @return the value of _lldp
 */
double MixNormal::llPercDiff() const {
  return _diffLogLikelihood;
} // end llPercDiff


//////////////////// startEpoch ////////////////////

/**
 * add new epoch to the em training
 */
void MixNormal::startEpoch() {
  unsigned d = _numVariables;

  // we reset the next parameters to zeros
  for ( unsigned i = 0; i < _numMixtures; i++ ) {
    _nextAlphas[i] = 0;
    for(unsigned j = 0; j < d; ++j){
      *(_nextMeans + i*d + j) = 0.0;
    }
    for(unsigned j = 0; j < d; ++j){
      *(_nextInvVars + i*d + j) = 0.0;
    }
  }
  
  // for full covariance matrices
  if ( _fullCoVar ) {
    for ( unsigned i = 0; i < _numMixtures; i++ ) {
      for(unsigned j = 0; j < d*d; ++j){
	*(_nextB + i*d*d + j) = 0.0; 
 	*(_nextCov + i*d*d + j) = 0.0;  
      }
    }
  }

  _preLogLikelihood = _logLikelihood;
  _logLikelihood = 0;
  _numAccum = 0;				// prepare for starting accumumlation
} // end startEpoch


//////////////////// addToEpoch ////////////////////

/**
 * accumulates the statistics 
 * parameters update are done in
 * the endEpoch method.
 *
 * @param pointerSet the pointer array to input data
 * @return the learning status
 */
bool MixNormal::addToEpoch(const PointerSetToDataPoints &pointerSet) {
  // update the parameters using the old value using:
  // alpha_l = \sum_i{p(l|x_i,theta)}/N
  // mue_l = \frac{\sum_i{x_i p(l| x_i, theta)}}{\sum_i{p(l|x_i,theta)}}
  // lambda_l = \frac{\sum_i{p(l|x_i,theta)}}{\sum_i{p(l|x_i,theta)(x_i-mue_i)^2}}
  // B_{ij} = (\sum_i{(x_i-(1-B)u)x_i' p(l|x,theta)})(\sum_i{x_i x_i' p(l|x_i,theta)})^(-1)

  PARAM_DATA_TYPE *pOut;
  PARAM_DATA_TYPE Px, prob_l;				// p(l|x,theta)
  PARAM_DATA_TYPE postProb[MAX_NUM_MIXTURES];              // vector of  p(l|x,theta)

  unsigned dim= _numVariables;

  for(unsigned i = 0; i < pointerSet.numSamples; i++){
    Prob_l_x_theta(Px, pointerSet,postProb,i);  //    postProb = Prob_l_x_theta(Px, x[i]);
    _logLikelihood -= log(Px);
    
    for(unsigned l = 0; l < _numMixtures; l++) {
      prob_l = postProb[l];
      _nextAlphas[l] += prob_l;

      for(unsigned j =0; j<dim;++j)        //_nextMeans[l] += x[i] * prob_l;
	*(_nextMeans + l*dim + j) += *(pointerSet.start[j] + i*pointerSet.skip) * prob_l;
      
      pOut = _nextCov+l*dim*dim;
      //  _nextCov[l] += (x[i]).toColVector() * (x[i]).toRowVector() * prob_l ;
      for(unsigned ii = 0; ii < dim; ++ii) // used to be vecProdAdd(_nextCov+l*dim*dim, pointerSet, i, prob_l, dim); 
	for(unsigned jj = 0; jj < dim; ++jj)
	    *(pOut + ii*dim + jj)  +=   *(pointerSet.start[ii] + i *pointerSet.skip) * *(pointerSet.start[jj] + i *pointerSet.skip) * prob_l;
    }
  }
  
  _numAccum += pointerSet.numSamples;
  
  return true;
} // end addToEpoch



//////////////////// addToEpochDiag ////////////////////

/**
 * same as addEpoch() except it is only for the case we have a
 * diagonal covariance.  I split these functions because they are the
 * inner loops of the program and need to be fast.  
 * accumulates the
 * statistics parameters update are done in the endEpoch method.
 *
 * @param x the Vector array of input data
 * @param numData the number of data
 * @return the learning status */
bool MixNormal::addToEpochDiag(const PointerSetToDataPoints &pointerSet) {
  // update the parameters using the old value using:
  // alpha_l = \sum_i{p(l|x_i,theta)}/N
  // mue_l = \frac{\sum_i{x_i p(l| x_i, theta)}}{\sum_i{p(l|x_i,theta)}}
  // lambda_l = \frac{\sum_i{p(l|x_i,theta)}}{\sum_i{p(l|x_i,theta)(x_i-mue_i)^2}}
  // B_{ij} = (\sum_i{(x_i-(1-B)u)x_i' p(l|x,theta)})(\sum_i{x_i x_i' p(l|x_i,theta)})^(-1)

  PARAM_DATA_TYPE *pOut;
  PARAM_DATA_TYPE Px, prob_l;				// p(l|x,theta)
  PARAM_DATA_TYPE postProb[MAX_NUM_MIXTURES];              // vector of  p(l|x,theta)

  unsigned dim= _numVariables;

  for(unsigned i = 0; i < pointerSet.numSamples; i++){
    Prob_l_x_theta(Px, pointerSet,postProb,i);  //    postProb = Prob_l_x_theta(Px, x[i]);
    _logLikelihood -= log(Px);
    
    for(unsigned l = 0; l < _numMixtures; l++) {
      prob_l = postProb[l];
      _nextAlphas[l] += prob_l;

      for(unsigned j =0; j<dim;++j)        //_nextMeans[l] += x[i] * prob_l;
	*(_nextMeans + l*dim + j) += *(pointerSet.start[j] + i*pointerSet.skip) * prob_l;

      pOut = _nextInvVars+l*dim;
      // we don't store the inverses here as the name _nextInvVars implies.  We will invert it in endEpoch();
      for ( unsigned k = 0; k < dim; ++k ) // used to be pointProdAdd(_nextInvVars+l*dim, pointerSet, i, prob_l,dim); 
	*pOut++ += *(pointerSet.start[k]+ i*pointerSet.skip) * *(pointerSet.start[k] + i * pointerSet.skip) * prob_l;
    }
  }

  _numAccum += pointerSet.numSamples;

  return true;
} // end addToEpochDiag



//////////////////// endEpoch ////////////////////

/**
 * end session of an epoch. this will update all the parameters
 * normalize  alphas, copy the value to previous
 * and prepare for next epoch
 *
 * @return the learning status
 */
bool MixNormal::endEpoch() {
  bool singular;
  unsigned dim = _numVariables;

  for ( unsigned l = 0; l < _numMixtures; l++ ) {
    // update alphas
    _alphas[l]  =  _nextAlphas[l] /_numAccum;
    // update means
    for ( unsigned i = 0; i < dim; i++ ) *(_means + l*dim+i) =  *(_nextMeans + l*dim+i) / _nextAlphas[l]; 

    // update variances and, in the case of a full covariance, covariances
    if ( _fullCoVar ) {
      for ( unsigned i = 0; i < dim*dim; i++ ) {
	*(_nextCov + l*dim*dim+i) /= _nextAlphas[l];
      }

      vecProdSub( _nextCov+l*dim*dim, _means+l*dim,_means+l*dim, dim); //_nextCov[l] -= _nextMeans[l].toColVector()*_nextMeans[l].toRowVector();

      if(covAddConst > 0)
	for ( unsigned i = 0; i < dim; i++ ) 
	  *(_nextCov + l*dim*dim+i*dim+i) += (covAddConst + ( rnd.drand48() - 0.5)*2*covAddEpsilon );
      
      if(clampCov > 0)
	for ( unsigned i = 0; i < dim; i++ ) {
	  if(*(_nextCov + l*dim*dim+i*dim+i) < clampCov && *(_nextCov + l*dim*dim+i*dim+i) >= 0) {
	    *(_nextCov + l*dim*dim+i*dim+i) = clampCov;
	    warning("Covariance of tuple # %d hit bottomf(%f). Clamping it.",_index,clampCov);
	  }
	}

      // Invert the covariance matrix
      singular = inverse(_nextCov+l*dim*dim, _tmpM, dim); //tmpM = _nextCov[l].inverse();
      if(singular) {
	fprintf(stderr,"While taking the inverse,\n");
	notPosDef(_index,l,dim,_nextCov+l*dim*dim);
      }
      
      // Perform Cholesky decomposition on the inverted matrix
      singular = diagCholeskyDecomp(_tmpM, _b+l*dim*dim, _invVars+l*dim, dim); 
      if(singular) {
	fprintf(stderr,"While doing Cholesky decomposition,\n");
	notPosDef(_index,l,dim,_tmpM);
      }
    } 
    else {
      for ( unsigned i = 0; i < dim; i++ ) *(_nextInvVars + l*dim+i) /= _nextAlphas[l]; 
      pointProdSub(_means+l*dim, _means+l*dim,_nextInvVars+l*dim, dim);  //_nextInvVars[l] -= _nextMeans[l].pointProd( _nextMeans[l]);
      invertVec(_nextInvVars+l*dim, _invVars+l*dim, dim); //_invVars[l] = _nextInvVars[l].inverse();
    }
  }

  normalize();

  _diffLogLikelihood = (_logLikelihood - _preLogLikelihood);
  _diffLogLikelihood = fabs(_diffLogLikelihood * 100.0 / _logLikelihood);

  return prepareNext();

} // end endEpoch



//////////////////// prob_x ////////////////////

/**
 * pdf using the current paramters P(x|theta)
 * = \sum_l{P(x|theta,l)P(l|theta)}
 *
 * @param ps the pointer set to the observation vector
 * @return the pdf using the current paramters
 */
PARAM_DATA_TYPE MixNormal::prob_x(const PointerSetToDataPoints &pointerSet,unsigned sample) const {
  PARAM_DATA_TYPE tmp = 0.0;

  for ( unsigned l = 0; l < _numMixtures; l++ ) {
    tmp += _alphas[l] * prob_x_theta_l(pointerSet, l,sample);
  }
  return tmp;
} // end prob_x

//////////////////// prob_x ////////////////////

/**
 * pdf using the current paramters P(x|theta)
 * = \sum_l{P(x|theta,l)P(l|theta)}
 *
 * @param vec the observation vector
 * @return the pdf using the current paramters
 */
PARAM_DATA_TYPE MixNormal::prob_x(const PARAM_DATA_TYPE* vec, unsigned vecLen) const {
  PARAM_DATA_TYPE sum = 0.0;

  assert(vecLen == _numVariables);
  
  for ( unsigned l = 0; l < _numMixtures; l++ ) {
    sum += _alphas[l] * prob_x_theta_l(vec, l,vecLen);
  }
  return sum;
} // end prob_x

//////////////////// prob_x ////////////////////

/**
 * pdf using the current paramters P(x|theta)
 * = \sum_l{P(x|theta,l)P(l|theta)}
 *
 * @param vec the observation vector
 * @param meanVecs -- mean vectors 
 * @param B -- B matrices 
 * @param varVecs -- variance vectors
 * @param dets -- determinants
 * @param inv_pow_sqrt_2pi
 * @return the pdf using the current paramters
 */
PARAM_DATA_TYPE MixNormal::prob_x(const PARAM_DATA_TYPE* vec, unsigned vecLen,  PARAM_DATA_TYPE* meanVecs, PARAM_DATA_TYPE* B, PARAM_DATA_TYPE* varVecs, PARAM_DATA_TYPE* dets, PARAM_DATA_TYPE inv_pow_sqrt_2pi) const {
  PARAM_DATA_TYPE sum = 0.0;
  
  for ( unsigned l = 0; l < _numMixtures; l++ ) {
    sum += _alphas[l] * prob_x_theta_l(vec, l,vecLen,meanVecs,B,varVecs,dets,inv_pow_sqrt_2pi);
  }
  return sum;
} // end prob_x

//////////////////// sampleComponent ////////////////////

/**
 * sample one component from the mixture
 *
 * @return a component label according to pdf
 */
unsigned MixNormal::sampleComponent() const {
  double tmp = rnd.drand48();

  unsigned l = 0;

  do {
    tmp -= _alphas[l];
    l++;
  } while ( tmp >0 && l < _numMixtures);

  return l - 1;
} // end sampleComponent


//////////////////// sample ////////////////////

/**
 * sample one set of value from the mixture
 *
 * @return a set of sample
 */
void MixNormal::sample(PARAM_DATA_TYPE *vecOutPtr) const {
  unsigned dim = _numVariables;
  unsigned l = sampleComponent();

  //  Another method to generate a sample ------------------
  //LDU2moment(_tmpM,_b+l*dim*dim,_invVars+l*dim,dim);
  //choleskyDecomp(_tmpM,_tmpM2,dim);
  //for ( unsigned i = 0; i < _numVariables; i++ ) {
  //vecOutPtr[i] = rnd.normal();                                                          
  //}
  //vecMatProduct(vecOutPtr,_tmpM2,dim);
  // for(unsigned j=0; j<dim; ++j) vecOutPtr[j] += *(_means+l*dim+j); 
  //--------------------------------------------------------

  for ( unsigned i = 0; i < _numVariables; i++ ) {
    vecOutPtr[i] = rnd.normal() / sqrt(*(_invVars+l*dim+i)); 
  }

  if (_fullCoVar) {
    // vecOutPtr = _b[l]^-1 * vecOutPtr;
    inverse(_b+l*dim*dim,_tmpM,dim);
    matVecProduct(_tmpM, vecOutPtr, dim);
  }
 
  for(unsigned j=0; j<dim; ++j) vecOutPtr[j] += *(_means+l*dim+j); 
} // end sample



//////////////////// randomize ////////////////////

/**
 * randomize the parameters
 */
void MixNormal::randomize() {
  DBGFPRINTF((stderr,"Randomizing parameters.\n")); 
  for ( unsigned l = 0; l < _numMixtures; l++ ) {
    randomize(l);
  }
  normalize();
} // end randomize


/**
 * randomize the prameters coresponding to mixture l
 *
 * @param l the label of mixture component
 */
void MixNormal::randomize(unsigned l) {
  DBGFPRINTF((stderr,"Randomizing parameters corresponding to mixture component %d.\n",l));
  unsigned dim = _numVariables;
  assert(l<_numMixtures);

  _alphas[l] = 1;
  for ( unsigned i = 0; i < _numVariables; i++ ) {
    *(_means + l*_numVariables +i) = 2.0 * rand() / RAND_MAX - 1.0;
    *(_invVars + l*_numVariables+i) = RAND_MAX * _reRandScale / (PARAM_DATA_TYPE)rand();
  }

  if ( _fullCoVar )
    unify(_b + l*dim*dim,dim);
} // end randomize


//////////////////// normalize ////////////////////

/**
 * normalize those alphas
 */
void MixNormal::normalize() {
  PARAM_DATA_TYPE tmp = 0.0;

  for ( unsigned l = 0; l < _numMixtures; l++ )
    tmp += _alphas[l];

  PARAM_DATA_TYPE invTmp = 1.0 / tmp;
  for ( unsigned l = 0; l < _numMixtures; l++ ) {
    _alphas[l] *= invTmp;
    // invDets is assigned the product of the elements of invVars
    for (unsigned j = 0; j < _numVariables; j++ )    *(_invDets +l) = 1.0;
    for (unsigned j = 0; j < _numVariables; j++ )
          *(_invDets + l) *= *(_invVars + l*_numVariables + j);
  }
} // end normalize


//////////////////// normalizeWithNoise ////////////////////

/**
 * normalize the alphas with noise
 */
void MixNormal::normalizeWithNoise() {
  PARAM_DATA_TYPE tmp = 0;

  unsigned l;
  for ( l = 0; l < _numMixtures; l++ )
    tmp += _alphas[l];

  // added to each one.
  PARAM_DATA_TYPE noise = tmp / _numMixtures;
  tmp += tmp;

  const PARAM_DATA_TYPE inv_tmp = 1.0 / tmp;
  for ( l=0; l < _numMixtures; l++ ) {
    _alphas[l] += noise;
    _alphas[l] *= inv_tmp;
  }
} // end normalizeWithNoise



//////////////////// prepareNext ////////////////////

/**
 * prepare the next epoch run
 *
 * @return returns true if there is a component drop, false otherwise
 */
bool MixNormal::prepareNext() {
  unsigned l = 0;
  bool mixCompDrop = false;
  unsigned dim = _numVariables;

  DBGFPRINTF((stderr,"At the start of prepareNext().\n"));
  DBGFPRINTF((stderr,"Relevant parameters are: mcvr(%.2f), varianceFloor(%.2f),  detFloor(%.2f).\n",mcvr,varianceFloor, detFloor));

start:
  if ( _totalCurrentReRands >= reRandsPerMixCompRedux ) {
    if ( _numMixtures > 1 ) {
      _totalCurrentReRands = 0;
      mixCompDrop = true;
      warning("Maximum number of rerandomizations exceeded. Dropping mixture component %d.\n",l);
      if ( noReRandOnDrop ){
	eliminateComponent(l);
      }
      else {
	_numMixtures--;  // since we are randomizing all components anyway, we just remove the last one.
	randomize();
      }

      _preLogLikelihood = 2 * LOG_ZERO;
      _logLikelihood = LOG_ZERO;
    } else {
      // what about the covariance matrix??
      _preLogLikelihood = _logLikelihood = -100;
    }
  }

  for ( l = 0; l < _numMixtures; l++ ) {
    if ( maxVal( (_invVars+l*dim), dim) * varianceFloor > 1) {	// the variance is too small
      warning("The variance of component %d hit the floor value %f\n",l,varianceFloor);
      if ( reRandOnlyOneComp ) {
	randomize(l);
	normalizeWithNoise();
	warning("Randomizing one component\n");
      } else {
	warning("Randomizing all components\n");
	randomize();
      }
      _totalCurrentReRands++;
      _preLogLikelihood = 2 * LOG_ZERO;
      _logLikelihood = LOG_ZERO;
      goto start;
    }

    const PARAM_DATA_TYPE det = 1.0 / *(_invDets + l*dim);
    if ( det < detFloor ) {				// det is too small
      if ( _numMixtures > 1 ) {
	if ( reRandOnlyOneComp ) {
	  randomize(l);
	  normalizeWithNoise();
	} else
	  randomize();
      }
      _totalCurrentReRands++;
      _preLogLikelihood = 2 * LOG_ZERO;
      _logLikelihood = LOG_ZERO;
      goto start;
    }

    const PARAM_DATA_TYPE invDet = 1.0 / det;

#if defined(WIN32)
    if ( ! _finite(invDet) ) {
#else
    if ( ! isfinite(invDet) ) {
#endif
      if ( reRandOnlyOneComp ) {
	randomize(l);
	normalizeWithNoise();
      } else
	randomize();

      _totalCurrentReRands++;
      _preLogLikelihood = 2 * LOG_ZERO;
      _logLikelihood = LOG_ZERO;
      goto start;
    }

    if ( _numMixtures * _alphas[l] * mcvr < 1.0 || _alphas[l] < DBL_MIN ) {
      warning("Component %d is below vanishing ratio.(%f)\n",l,mcvr);
      if ( reRandOnlyOneComp ) {
	randomize(l);
	normalizeWithNoise();
      } else
	randomize();

      _totalCurrentReRands++;
      _preLogLikelihood = 2 * LOG_ZERO;
      _logLikelihood = LOG_ZERO;
      goto start;
    }
  }

  return mixCompDrop;
} // end prepareNext


//////////////////// eliminateComponent ////////////////////

/**
 * eliminate one component from the mixture
 *
 * @param l the label of mixture component
 * @throws OutOfRangeException when l >= number of mixtures
 */
void MixNormal::eliminateComponent(unsigned l) {
  DBGFPRINTF((stderr,"Eliminating a mixture component in  MixNormal::eliminateComponent(unsigned l).\n"));
  assert(l < _numMixtures);

  // we don't bother to create a new buffer
  for ( unsigned j = l + 1; j < _numMixtures; j++ ) {
    _alphas[j] = _alphas[j-1];
    for ( unsigned k = 0; k < _numVariables; k++ ) {
      *(_means + j*_numVariables + k)  =  *(_means + (j-1)*_numVariables + k);  ///_means[j] = _means[j-1];
      *(_invVars + j*_numVariables + k)  =  *(_invVars + (j-1)*_numVariables + k); //_invVars[j] = _invVars[j-1];
      *(_b + j*_numVariables*_numVariables + k)  =  *(_b + (j-1)*_numVariables*_numVariables + k); // _b[j] = _b[j-1];
      *(_invDets + j*_numVariables + k)  =  *(_invDets + (j-1)*_numVariables + k); //_invDets[j] = _invDets[j-1];
    }
  }
  _numMixtures--;
  normalize();
} // end eliminateComponent



//////////////////// Prob_l_x_theta ////////////////////

/**
 * calculate the conditional probability under the current parameters
 * by Bayes rule, we get P(l|x,theta) = P(x|l,theta)P(l|theta)/P(x|theta)
 * where the P(l|theta) = alpha[l]
 *
 * @param x the observation vector
 * @param  P(l|x,theta) vector of probabilities for each mixture
 */
void MixNormal::Prob_l_x_theta(PARAM_DATA_TYPE &Px, 
			       const PointerSetToDataPoints &pointerSet, 
			       PARAM_DATA_TYPE * postProbOutPtr,
			       unsigned sample) const {
  Px = 0;
  for(unsigned l = 0; l < _numMixtures; l++){
    postProbOutPtr[l] = prob_x_theta_l(pointerSet, l,sample) * _alphas[l];
    Px = Px + postProbOutPtr[l];
  }
  if( Px <= FLT_MIN ){
    for(unsigned l = 0; l < _numMixtures; l++) postProbOutPtr[l] = FLT_MIN;
  }
  else{
    for(unsigned l = 0; l < _numMixtures; l++) postProbOutPtr[l] = postProbOutPtr[l] / Px;
  }
} // end prob_l_x_theta


//////////////////// prob_x_theta_l ////////////////////

/**
 * calculate the conditional probability under the current parameters
 *
 * @param pointerSet pointers to the observation vector
 * @param l the label of mixture component
 * @param sample the smaple we want to access via pointerSet
 * @return P(x|theta_l)
 */
PARAM_DATA_TYPE MixNormal::prob_x_theta_l(const PointerSetToDataPoints &pointerSet, unsigned l, unsigned sample) const {
  assert( l < _numMixtures );
  unsigned dim = _numVariables;
  PARAM_DATA_TYPE d[MAX_POINTER_SET_DIM];
  for(unsigned i=0; i<dim;++i) {
    d[i] = *(pointerSet.start[i] + sample*pointerSet.skip) - *(_means+l*dim +i);
  }

  if ( _fullCoVar ) matVecProduct( (_b+l*dim*dim), d,dim);  // d = _b[l] * d;
  double tmp = exp(-0.5 * dot( (_invVars+l*dim), d,d,dim)) * sqrt( *(_invDets+ l) ) * _inv_pow_sqrt_2pi;
  return (PARAM_DATA_TYPE) tmp;
} // end prob_x_theta_l

//////////////////// prob_x_theta_l ////////////////////

/**
 * calculate the conditional probability under the current parameters
 *
 * @param vec the observation vector
 * @param vecLen the length of the vector
 * @param l the label of mixture component
 * @return P(x|theta_l)
 */
PARAM_DATA_TYPE MixNormal::prob_x_theta_l(const PARAM_DATA_TYPE* vec, unsigned l, unsigned vecLen) const {
  assert(l < _numMixtures);
  //if ( l >= _numMixtures ) return 0;

  unsigned dim = _numVariables;

  assert(vecLen == dim);

  PARAM_DATA_TYPE d[MAX_POINTER_SET_DIM];

  for(unsigned i=0; i<dim;++i) {
    d[i] = *(vec + i) - *(_means+l*dim +i);
  }
  if ( _fullCoVar ) matVecProduct( (_b+l*dim*dim), d,dim);  // d = _b[l] * d;

  double tmp = exp(-0.5 * dot( (_invVars+l*dim), d,d,dim)) * sqrt( *(_invDets+ l) ) * _inv_pow_sqrt_2pi;
  return (PARAM_DATA_TYPE) tmp;
} // end prob_x_theta_l


//////////////////// prob_x_theta_l ////////////////////

/**
 * calculate the conditional probability under the current parameters
 *
 * @param vec the observation vector
 * @param l the label of mixture component
 * @param vecLen the length of the vector
 * @param meanVecs -- mean vectors 
 * @param B -- B matrices 
 * @param varVecs -- variance vectors
 * @param dets -- determinants
 * @param inv_pow_sqrt_2pi
 * @return P(x|theta_l)
 */
PARAM_DATA_TYPE MixNormal::prob_x_theta_l(const PARAM_DATA_TYPE* vec, unsigned l, unsigned vecLen,  PARAM_DATA_TYPE* meanVecs, PARAM_DATA_TYPE* B, PARAM_DATA_TYPE* varVecs, PARAM_DATA_TYPE* dets, PARAM_DATA_TYPE inv_pow_sqrt_2pi) const {
  assert(l < _numMixtures);
  //if ( l >= _numMixtures ) return 0;

  unsigned dim = vecLen;

  PARAM_DATA_TYPE d[MAX_POINTER_SET_DIM];

  for(unsigned i=0; i<dim;++i) {
    d[i] = *(vec + i) - *(meanVecs+l*dim +i);
  }
  if ( _fullCoVar ) matVecProduct( (B+l*dim*dim), d,dim);  // d = _b[l] * d;

  double tmp = exp(-0.5 * dot( (varVecs+l*dim), d,d,dim)) * sqrt( *(dets+ l) ) * inv_pow_sqrt_2pi;
  return (PARAM_DATA_TYPE) tmp;
} // end prob_x_theta_l


/**
 * for each mixture, calculate the B matrix and the inverse variances
 * from the covariance matrix using LDU decomposition 
 *
 *Preconditions:
 * _cov has to be initialized through kmeans and storage allocated to
 * _b, _invVars and _invDets 
 *
 *Postconditions: 
 * Matrices _b, _invVars and and invDets are updated.
 */
void MixNormal::calcB(){
  unsigned dim = _numVariables;
  bool malFormed = false;
  for( unsigned l = 0; l < _numMixtures; l++ ) {
    malFormed =  moment2LDU(_cov+l*dim*dim, _b+l*dim*dim, _invVars+l*dim,dim);
    // error if the covariance is not positive definite.  No need to
    // check the invVars because the covariance is positive definite
    // iff the diagonal entries in the LDU decomposition are positive
    if(malFormed) {
      notPosDef(_index, l, dim, _cov+l*dim*dim);
    }
    *(_invDets+l) = prod(_invVars+l*dim,dim); //*(_invDets+l) = _invVars[l].prod();
  }
}
 



//////////////////// cleanup ////////////////////

/**
 * cleanup the memory
 */
void MixNormal::cleanup() {
  delete [] _alphas;      _alphas = NULL;
  delete [] _means;       _means = NULL;
  delete [] _invVars;     _invVars = NULL;
  delete [] _invDets;	  _invDets = NULL;

  delete [] _nextAlphas;  _nextAlphas = NULL;
  delete [] _nextMeans;   _nextMeans = NULL;
  delete [] _nextInvVars; _nextInvVars = NULL;

  delete [] _meansX;   _meansX   = NULL;
  delete [] _invDetsX; _invDetsX = NULL;
  delete [] _invCovsX; _invCovsX = NULL;

  delete [] _invVarsX;     _invVarsX = NULL;

  delete [] _meansY;   _meansY   = NULL;
  delete [] _invDetsY; _invDetsY = NULL;
  delete [] _invCovsY; _invCovsY = NULL;

  delete [] _invVarsY;     _invVarsY = NULL;

  delete [] _tmpM; _tmpM = NULL;
  delete [] _tmpM2; _tmpM2 = NULL;

  // for full covariance matrices
  if ( _fullCoVar ) {
    delete [] _b;           _b = NULL;
    delete [] _bX;           _bX = NULL;
    delete [] _bY;           _bY = NULL;
    delete [] _nextB;       _nextB = NULL;
    delete [] _nextCov;     _nextCov = NULL;
  }

  _numVariables = 0;
  _numMixtures  = 0;
} // end cleanup


//// For debugging ////


/**
 * for each mixture, do the inverse LDU to check if calcB() is correct
 *
 */ 
void MixNormal::checkB(){
  fprintf(stderr,"In checkB()");
  PARAM_DATA_TYPE* newCov = new PARAM_DATA_TYPE[_numVariables*_numVariables];
  unsigned dim = _numVariables;
  for( unsigned l = 0; l < _numMixtures; l++ ) {
    fprintf(stderr,"%d",l);
    LDU2moment(newCov, _b+l*dim*dim, _invVars+l*dim,dim);
    fprintf(stderr,"InvVars: ");
    for( unsigned i = 0; i < _numVariables; i++ )
      fprintf(stderr,"%f ",*(_invVars+l*dim+i));
    fprintf(stderr,"\n");
    fprintf(stderr,"cov,  b: \n");
    for( unsigned i = 0; i < _numVariables*_numVariables; i++ ) { 
      fprintf(stderr,"%f %f %f\n",*(newCov+i),*(_cov+l*dim*dim+i),*( _b+l*dim*dim + i));
    }
  }
}
