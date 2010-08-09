
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
#include "lms_filter.h"

VCID("$Header$")

LMSFilter::LMSFilter(unsigned _order, double _learning_rate)
: AdaptiveFilter(_order), learning_rate(_learning_rate)
{
  // printf("order = %d\n",_order);

  weights.resize(order);
  history.resize(order);


  // start off saying that the most recent sample is at the end.
  mostRecentSamplePosition = (order - 1);

  init();

}

void LMSFilter::init()
{
  // we've loaded no samples
  numSamplesLoaded = 0;

  // initialize to random values normalized to unity.
  bool randomInit = false;


  if (randomInit) {
    double sum = 0;
    for (unsigned i = 0; i < order ; i++ ) {
      sum += (weights.ptr[i] = rnd.drand48pe());
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


double LMSFilter::makePrediction()
{
  if (!readyToMakePrediction())
    return 0; // this should really be an error, or we could fake a prediction.

  double res = 0;
  for (unsigned i= 0; i < order; i++ ) {
    res += weights.ptr [i] * history.ptr [ (mostRecentSamplePosition + i ) % order ];
  }

  return res;

}


void LMSFilter::addNextSampleAndUpdate(double val) 
{

  if (readyToMakePrediction()) {
    // get the current prediction and error
    double prediction = makePrediction();
    double error = val - prediction;
  
    // now update the filter  
    double norm = 0;
    // TODO: optimize this recursively, subtract off old, and add on new
    for (unsigned i = 0; i < order ; i ++) {
      norm += history.ptr[i]*history.ptr[i];
    }
    // printf("norm = %f\n",norm);
    if (norm <= 1e-15)
      norm = 1;
    else
      norm = 1.0/norm;
    for (unsigned i = 0; i < order ; i ++) {
      weights.ptr[i] +=
	norm * learning_rate * error * history[ (mostRecentSamplePosition + i ) % order ];
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
      = cos(i * i* 2 * M_PI * 2 / len ) + rnd.normal()*.125/2.0;
  }

  unsigned order = 3;
  if (argc > 1)
    order = atoi(argv[1]);
  double lr = 0.5;
  if (argc > 2)
    lr = atof(argv[2]);

  AdaptiveFilter * l1 = new LMSFilter(order, lr);
  AdaptiveFilter * ff = new FixedFilter();


  double err1 = 0;
  double err2 = 0;
  double sumsq = 0;
  double sum_err_ratio = 0;
  unsigned count = 0;
  for (unsigned i = 0; i < len; i++ ) {
    l1->addNextSampleAndUpdate(signal[i]);
    ff->addNextSampleAndUpdate(signal[i]);
    if (l1->readyToMakePrediction() && (i+1) < len) {
      count ++;
      double pred1 = l1->makePrediction();
      double pred2 = ff->makePrediction();
      sumsq += signal[i+1]*signal[i+1];
      err1 += (signal[i+1] - pred1)*(signal[i+1] - pred1);
      err2 += (signal[i+1] - pred2)*(signal[i+1] - pred2);
      double rat = (signal[i+1] - pred1)*(signal[i+1] - pred1)
	/ ( (signal[i+1] - pred2)*(signal[i+1] - pred2) ) ;
      printf("i=%d: signal[i+1] = %f, pred1 = %f, pred2 = %f, rat = %f\n",i,
	     signal[i+1],pred1,pred2,rat);
      sum_err_ratio += rat;

    }
  }
  printf("avg err1 = %f, avg err2 = %f, err1/err2 = %f, sum_err_rat = %f\n",
	 err1/count,err2/count,err1/err2,sum_err_ratio/count);

  delete l1;
  delete ff;

}


#endif
