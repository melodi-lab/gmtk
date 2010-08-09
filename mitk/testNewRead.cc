#include <iostream>
#include "readRange.h"
#include "ObservationMatrix.h"


void main() {




}


//////////////////// bufferToVector ////////////////////

/**
 *  extracts features corresponding to range and puts them 
 *  into vectors
 */
void MixNormalCollection::bufferToVector(ProcVector *&x,
					 float *dataBuffer, 
					 size_t featureSize,
					 size_t numberOfSample,
					 RangeSet range){ 
  unsigned d, d1, d2, i, j;
  float **featureP;

  d1 = range.getSize(0);
  d2 = range.getSize(1);
  d = d1 + d2;

  featureP = new float * [d];
  for( i = 0; i < d1; i++ ) 
    featureP[i] = dataBuffer + featureSize * range.set[0][i].lag 
                  + range.set[0][i].feat;
  for( i = 0; i < d2; i++ ) 
    featureP[i+d1] = dataBuffer + featureSize * range.set[1][i].lag 
                     + range.set[1][i].feat;    
  for( i = 0; i < (unsigned) numberOfSample - range.max_lag; i++ ){
    x[i].resize(d);
    for( j = 0; j < d; j++ ){
      x[i][j] = *featureP[j];
      featureP[j] += featureSize;
    }     
  }
  delete [] featureP; featureP = NULL;
}





