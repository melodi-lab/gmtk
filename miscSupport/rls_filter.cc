
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
#include "rls_filter.h"

VCID("$Header$")

RLSFilter::RLSFilter(unsigned _order, double _forgetting_coef)
: AdaptiveFilter(_order), forgetting_coef(_forgetting_coef)
{

  weights.resize(order);
  history.resize(order);
  concMat.resize(order*order);

  tmpVec.resize(order);
  double sum = 0;
  for (unsigned i = 0; i < order ; i++ ) {
    sum += (weights.ptr[i] = rnd.uniform(1.0));
    for (unsigned j = 0; j < order ; j++ ) {
      // initialize to a diagonal
      if (i == j) {
	concMat.ptr[order* i + j] = 1;
      } else
	concMat.ptr[order* i + j] = 0.0;
    }
  }
  sum = 1.0/sum;
  for (unsigned i = 0; i < order ; i++ ) {
    weights.ptr[i] *= sum;
  }

  
  // start off saying that the most recent sample is at the end.
  mostRecentSamplePosition = (order - 1);
  
  // we've loaded no samples
  numSamplesLoaded = 0;

}


void RLSFilter::init()
{
  // we've loaded no samples
  numSamplesLoaded = 0;

  // initialize to random values normalized to unity.
  bool randomInit = false;

  // position 0 of weights is most recent
  // position 1 of weights is next most recent
  // etc.

  if (randomInit) {
    double sum = 0;
    for (unsigned i = 0; i < order ; i++ ) {
      sum += (weights.ptr[i] = rnd.uniform(1.0));
    }
    if (sum == 0)
      sum = 1.0;
    else
      sum = 1.0/sum;
    for (unsigned i = 0; i < order ; i++ ) {
      weights.ptr[i] *= sum;
    }
  } else {
    // init to simple 2nd order filter if possible.
    if (order == 1) 
      weights.ptr[0] = 1.0;
    else {
      for (unsigned i = 0; i < order ; i++ ) {
	weights.ptr[i] = 0;
      }
      // simple fixed filter weights rest are zero.
      weights.ptr[0] = 2.0;
      weights.ptr[1] = - 1.0;
    }
  }
}

double RLSFilter::makePrediction()
{
  if (!readyToMakePrediction())
    return 0; // this should really be an error, or we could fake a prediction.

  double res = 0;
  for (unsigned i= 0; i < order; i++ ) {
    res += weights.ptr [i] * history.ptr [ (mostRecentSamplePosition + i ) % order ];
  }

  return res;

}


void RLSFilter::addNextSampleAndUpdate(double val) 
{

  if (readyToMakePrediction()) {
    // first, update the concentration matrix


    // first do the mat-vec operation tmpvec = concMat * history which
    // results in a quantity that we use multiple times.
    for (unsigned i=0; i<order; i++) {
      tmpVec.ptr[i] = 0.0;
    }
    for (unsigned i=0; i<order; i++) {
      for (unsigned j=0;j<order;j++) {
	tmpVec.ptr[i] += 
	  concMat.ptr[i*order + j]*
	  history.ptr[ (mostRecentSamplePosition + j ) % order ];
      }
    }

    // next update the concentration matrix
    double tmp1 = forgetting_coef;
    for (unsigned i=0; i<order; i++) {
      tmp1 += tmpVec.ptr[i]*
	history[ (mostRecentSamplePosition + i ) % order ];
    }
    // printf("tmp1 = %f\n",tmp1);
    tmp1 = 1.0/tmp1;

    double inv_fc = 1.0/forgetting_coef;
    for (unsigned i=0; i<order; i++) {
      for (unsigned j=0;j<order;j++) {
	concMat.ptr[i*order + j] =
	  inv_fc*(concMat.ptr[i*order + j]
		- tmpVec.ptr[i]*tmp1*tmpVec.ptr[j]);
      }
    }

    // update the weights

    // get the current prediction and error
    double prediction = makePrediction();
    double error = val - prediction;


#define JORDAN 0
#if JORDAN
    // recompute the new mat-vec operation tmpvec = concMat * history
    // with the new confMat
    for (unsigned i=0; i<order; i++) {
      tmpVec.ptr[i] = 0.0;
    }
    for (unsigned i=0; i<order; i++) {
      for (unsigned j=0;j<order;j++) {
	tmpVec.ptr[i] += 
	  history[ (mostRecentSamplePosition + j ) % order ]
	  *concMat[i*order + j];
      }
    }
#endif

    for (unsigned i = 0; i < order ; i ++) {
#if JORDAN
      weights.ptr[i] += tmpVec.ptr[i]*error;
#else
      weights.ptr[i] += tmpVec.ptr[i]*error*tmp1;
#endif

    }
  }

  // update history
  unsigned newPos = (mostRecentSamplePosition + order - 1) % order;
  numSamplesLoaded ++;
  history[newPos] = val;
  mostRecentSamplePosition = newPos;

}




////////////////////////////////////////////////////////////////////
//        Test Driver
////////////////////////////////////////////////////////////////////

#ifdef MAIN

#include "fileParser.h"

#include "fixed_filter.h"

RAND rnd(true);

int
main(int argc, char*argv[])
{
  // first create a simple slowly varying sinusoidal signal.

  const unsigned len = 1000;

  double signal[len];

  // have the signal go through 4 periods.
  for (unsigned i = 0; i < len; i++ ) {
    signal[i]
      = cos(i * i* 2 * M_PI * 2 / (double)len ) + rnd.normal()*.125/2.0;
  }

  unsigned order = 3;
  if (argc > 1)
    order = atoi(argv[1]);
  double fc = 1.0;
  if (argc > 2)  
    fc = atof(argv[2]);


  RLSFilter * l1 = new RLSFilter(order, fc);
  FixedFilter * ff = new FixedFilter();

  double err1 = 0;
  double err2 = 0;
  double sumsq = 0;
  unsigned count = 0;
  for (unsigned i = 0; i < len; i++ ) {
    l1->addNextSampleAndUpdate(signal[i]);
    ff->addNextSampleAndUpdate(signal[i]);
    if (l1->readyToMakePrediction() && (i+1) < len) {
      count ++;
      double pred1 = l1->makePrediction();
      double pred2 = ff->makePrediction();
      printf("i=%d: signal[i+1] = %f, pred1 = %f, pred2 = %f\n",i,
	     signal[i+1],pred1,pred2);
      sumsq += signal[i+1]*signal[i+1];
      err1 += (signal[i+1] - pred1)*(signal[i+1] - pred1);
      err2 += (signal[i+1] - pred2)*(signal[i+1] - pred2);
    }
  }
  printf("avg err1 = %f, avg err2 = %f, err1/err2 = %f\n",
	 err1/count,err2/count,err1/err2);

  delete l1;
  delete ff;

}


#endif
