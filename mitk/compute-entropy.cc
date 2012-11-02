#include <iostream>
#include "mixNormalCollection.h"
#include "string.h"

//// USING LLL ///////////////

//////////////////// MixNormal::computeEntropy ////////////////////

/** 
 * Compute the entropy of the reandom variable X using the law of large
 * numbers. First generate nSamples samples from the joint
 * distribution of X and Y and then calculate the entropy as:
 *  H(X) = sum log( 1/p(x) ) / nSamples 
 */
PARAM_DATA_TYPE MixNormal::computeEntropy(unsigned  nSamples, PARAM_DATA_TYPE& H) {
  PARAM_DATA_TYPE probX; 
  PARAM_DATA_TYPE sampleX[MAX_POINTER_SET_DIM];

  H = 0; 

  for( unsigned i = 0; i < nSamples; i++){
    sample(sampleX);   //    sampleX = sample();
#if 0
#if DEBUG
    for (unsigned i=0; i<_numVariables; ++i) {
      DBGFPRINTF((stderr,"%f ",*(sampleX+i)));
    }
    DBGFPRINTF((stderr,"\n"));
#endif
#endif   
    probX  = prob_x(sampleX,_numVariables);
    assert( (probX >= 0) );
    H  -= log(probX);
  }

  PARAM_DATA_TYPE invNlog2 = 1.0 / (PARAM_DATA_TYPE)(nSamples*log(2.0));
  H *= invNlog2;

  return H;
}


//// USING DATA ///////////

//////////////////// MixNormal::startEpochEntropy ////////////////////

/**
 * Initialize the necessary variables used in MI calculation by LLN
 * using "real" data.
 */
void MixNormal::startEpochEntropy() {
  _Hx  = 0;
  _Hy = 0;  //
  _Hxy = 0; // these two initializationa are needed because the same endEpoch routine used to calculate MI and entropy 
  _nSamplesMI = 0;
}

//////////////////// MixNormal::addToEpochEntropy ////////////////////

/**
 * Add nSamples data points to compute the mutual information 
 * between X and Y by using the actual samples. The input is an 
 * array of pointers to the data locations.
 *
 *  H(X) = sum log( 1/p(x) ) / N 
 */
void MixNormal::addToEpochEntropy(PointerSetToDataPoints& ps) {
  PARAM_DATA_TYPE probX; 

  for( unsigned i = 0; i < ps.numSamples; i++){
    probX  = prob_x(ps,i);
    // Have to get to the bottom of the differenc between prob_x_GM and prob_x
#if 0
    unsigned dim = _numVariables;
    for(unsigned l = 0; l < _numMixtures; l++) {
       for(unsigned i = 0; i < dim; i++) 
              *(_invCovsXY+l*dim*dim+i) = *(_invVars+l*dim + i);
     }
    probX  = prob_x_GM(ps,i, _alphas, _means, _invCovsXY, _invDets, _inv_pow_sqrt_2pi);
#endif
    DBGFPRINTF((stderr,"In MixNormal::addToEpochEntropy probX = %f\n",probX));
    assert( (probX >= 0) );
    _Hx  = _Hx  - log(probX);
  }
  _nSamplesMI = _nSamplesMI + ps.numSamples;
}

//////////////////// MixNormal::endEpochEntropy ////////////////////

/**
 * Finish MI calculation by normalizing with N and converting it to bits 
 * from nats.
 *
 *  I = sum log( p(x,y)/p(x)p(y) / N 
 */
PARAM_DATA_TYPE MixNormal::endEpochEntropy(PARAM_DATA_TYPE& H) {

  PARAM_DATA_TYPE invNlog2 = 1.0 / (PARAM_DATA_TYPE)(log(2.0)*_nSamplesMI);

  H *= invNlog2;

  return H;
}
