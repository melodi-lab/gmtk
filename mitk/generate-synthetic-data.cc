#include "mixNormalCollection.h"
#include "mixNormal.h"
#include "matrix-ops.h"
#include "rand.h"

//extern bool Seed;
//RAND rnd(Seed);
//RAND rnd(true);

//////////////////// sampleUsingCov ////////////////////

/**
 * sample one set of value from the mixture
 *
 * @return a set of sample
 */
void MixNormal::sampleUsingCov(double *sampleVec) const {
  unsigned l = sampleComponent();
  unsigned n = _numVariables;
  double* b = new double [n*n];
  double* tmpM = new double [n*n];
  double* var = new double [n];
  bool malFormed =false;

   //moment2LDU(cov,b+l*n*n,vars+l*n,n);
  malFormed = diagCholeskyDecomp(_cov+l*n*n,b,var,n); 
  if(malFormed) {
    error("ERROR: in MixNormal::sampleUsingCov(), covariance is not positive definite.");
  }
  //inverse(b,tmpM,n);
  //transpose(b,tmpM,n);
  transpose(tmpM,b,n);
  //-------------------------------------

  for ( unsigned i = 0; i < n; i++ ) {
    //sampleVec[i] = inverse_normal_func(drand48()) *  sqrt(*(var+i));
    sampleVec[i] = rnd.normal() *  sqrt(*(var+i));
  }
  //inverse(b+l*n*n,tmpM,n);
  matVecProduct(tmpM, sampleVec, n);
  //matVecProduct(b, sampleVec, n);
 
  for(unsigned j=0; j<n; ++j) sampleVec[j] += *(_means+l*n+j); 
} // end sample



void MixNormal::generateData(FILE* ofp,unsigned numSamples) {
  unsigned n = _numVariables;
  double* sampleVec = new double [n];

  for(unsigned i = 0; i< numSamples; ++i) {
    sample(sampleVec);
    for(unsigned j=0; j<n; ++j) {
      fprintf(ofp,"%.3f ",sampleVec[j]);
    }
    fprintf(ofp,"\n");
  }
}



void MixNormal::generateDataUsingCov(FILE* ofp,unsigned numSamples) {
  unsigned n = _numVariables;
  double* sampleVec = new double [n];

  for(unsigned i = 0; i< numSamples; ++i) {
    sampleUsingCov(sampleVec);
    for(unsigned j=0; j<n; ++j) {
      fprintf(ofp,"%.3f ",sampleVec[j]);
    }
    fprintf(ofp,"\n");
  }
}



void MixNormalCollection::generateData(FILE* ofp,unsigned numSamples,unsigned mixtureNum) {
  MixNormal * p = _ftrMI + mixtureNum;
  p->generateData(ofp,numSamples);
}


void MixNormalCollection::generateDataUsingCov(FILE* ofp,unsigned numSamples, unsigned mixtureNum) {
  MixNormal * p = _ftrMI + mixtureNum;
  p->generateDataUsingCov(ofp,numSamples);
}
