/**
 *: mixNormalCollection.cc
 */

#include <iostream>

#include "mixNormalCollection.h"
#include "mixNormal.h"
#include "range.h"
#include "readRange.h"
#include "error.h"

/* Constructor */
MixNormalCollection::MixNormalCollection(RangeSetCollection rangeSetCol) :
  _numMis(rangeSetCol.getSize()), 
  _ftrMI(new MixNormal [_numMis]), 
  _ftrMI_endp(_ftrMI + _numMis) 
{
  int k=0;
  for ( MixNormal * p = _ftrMI; p != _ftrMI_endp; p++ ) {
    int dimX =  rangeSetCol.rs[k].getSize(0);
    int setSize = dimX +  rangeSetCol.rs[k].getSize(1);
    p->setInit(setSize, dimX);
    ++k;
  }
}
/**
 * prepare to start an EM epoch
 */
void MixNormalCollection::startEpoch() {
  //Initialize all components
  //for ( MixNormal * p = _ftrMI; p != _ftrMI_endp; p++ )
  //p->startEpoch();
}

/**
 * *
 * @param obsMat ObservationMatrix
 * @param n_ftrs number of features
 * @param n_samps number of samples
 * @param rng the pointer to set of nodes
 * @return true if no componenet dropped
 */

bool MixNormalCollection::addToEpoch(ObservationMatrix* obsMat,
				     const size_t n_ftrs,
				     const size_t n_frames,
				     const size_t n_samps,
				     unsigned firstFrame,
				     const RangeSetCollection &rngCol) {
  bool dataUsed = true;
  RangeSet rng;
  int setSize1, setSize2, setSize, j, k = 0;
  PtrArray<unsigned> ptrArray(obsMat);

  for ( MixNormal *p = _ftrMI; p != _ftrMI_endp; p++ ) {
      rng = rngCol.rs[k];
      setSize1 = rng.getSize(0);
      setSize2 = rng.getSize(1);
      setSize = setSize1 + setSize2;
      ptrArray.setSize(setSize);

      unsigned* dftr_buf = (unsigned*) obsMat->features.ptr;
      
      int startFrame = firstFrame;
      if(rng.min_lag < 0 && (-rng.min_lag > (int) firstFrame)) {
	startFrame = - rng.min_lag;
      }

      dftr_buf = (unsigned*) (dftr_buf + startFrame * n_ftrs);

      int endFrame = firstFrame + n_samps - 1;
      if(rng.max_lag > 0 && 
	 ((int)n_frames - rng.max_lag - 1) < (int) endFrame) {
	endFrame = n_frames - rng.max_lag - 1;
      }

      if(endFrame < startFrame) {
	dataUsed = false;
	continue;
      }
 	
      // initialize the buff pointers
      for ( j = 0; j < setSize1; j++ ) {
	ptrArray.ptrVec[j] = dftr_buf + n_ftrs * rng.set[0][j].lag
	  + rng.set[0][j].feat;
      }
      for ( j = 0; j < setSize2; j++ ) {
	ptrArray.ptrVec[j+setSize1] = dftr_buf + n_ftrs * rng.set[1][j].lag
	  + rng.set[1][j].feat;
      }

      // now feed the data into the EM learning
      p->addToEpoch(&ptrArray, endFrame - startFrame + 1);
  
      k++;
  }
  dataUsed = true;

  return dataUsed;
}


/**
 *   For each range spec, compute the mutual information
 *
 *  */
void MixNormalCollection::endEpoch(const size_t n_samps,
				   const RangeSetCollection &rngCol,
				   FILE* ofp) 
{
  RangeSet rng;
  double I;
  int k=0, numSamples;
  for ( MixNormal * p = _ftrMI; p != _ftrMI_endp; p++ ) {
    rng = rngCol.rs[k];
    numSamples = n_samps - rng.max_lag; 
    I = p->endEpoch(numSamples);
    fprintf(ofp,"%e\n", I);
    ++k;
  }
}






