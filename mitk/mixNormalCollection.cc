#include <iostream>



#include "mixNormalCollection.h"
#include "mixNormal.h"
#include "range.h"
#include "readRange.h"
#include "error.h"

using namespace std;

bool MixNormalCollection::noSamplesFound(int& rangeSpecNum) {

  rangeSpecNum = -1;
  int k = 0;
  for ( MixNormal * p = _ftrMI; p != _ftrMI_endp; p++ ) {
    if( p->active() && p->getNumSamples() == 0) {
      rangeSpecNum = k;
      return true;
    }
    k++;
  }
  return false;
}


void MixNormalCollection::printNumSamples() {

  int k = 0;
  for ( MixNormal * p = _ftrMI; p != _ftrMI_endp; p++ ) {
    cout<<"Mixture no "<<k<<" : "<<p->getNumSamples()<<" samples\n";
    k++;
  }
}


/**
 * Second constructor
 *
 * @param numMis number of items in the collection
 */
MixNormalCollection::MixNormalCollection(RangeSetCollection rangeSetCol, unsigned numComps,unsigned maxKMeansIter,bool fullCovar,float covAddConst,float covAddEpsilon, double clampCov):  _numComps(numComps), _numMis(rangeSetCol.getSize()), _ftrMI(new MixNormal [_numMis]), _ftrMI_endp(_ftrMI + _numMis)
{
  
  int k=0;
  for ( MixNormal * p = _ftrMI; p != _ftrMI_endp; p++ ) {
    int setSize = rangeSetCol.rs[k].getSize(0) +  rangeSetCol.rs[k].getSize(1);
    p->setActive();
    p->setInit(k,setSize, numComps, maxKMeansIter,fullCovar,covAddConst,covAddEpsilon,clampCov);
    ++k;
  }
} // end MixNormalCollection


// NOT USED.
/**
 * set the initial parameters for MixNormal
 *
 * @param dimension the dimension of the gaussian mixture
 * @param numComps the number of mixture componenets
 * @param fullCoVar true if use full covariace matrix
 */
/*void MixNormalCollection::setInitParams(unsigned dimension,
  unsigned numComps,
  unsigned maxKMeansIter,
  bool fullCoVar) {
  unsigned i=0;
  for ( MixNormal * p = _ftrMI; p != _ftrMI_endp; p++ ) {
  p->setInit(i,dimension, numComps, maxKMeansIter,fullCoVar,covAddConst,covAddEpsilon);
  ++i;
  }
}
*/


/**
 * set the variance floor for gaussian mixture
 *
 * @param vf the variance floor
 */
void MixNormalCollection::setVarianceFloor(double vf) {
  MixNormal::varianceFloor = vf;
}


/**
 * set the determinent floor for gaussian mixture
 *
 * @param df the determinent floor
 */
void MixNormalCollection::setDetFloor(double df) {
  MixNormal::detFloor = df;
}


/**
 * set the mixture coefficient vanish number
 *
 * @param mcvf the ratio that eleminate a component
 */
void MixNormalCollection::setMixtureCoeffVanishNumber(float mcvf) {
  MixNormal::mcvr = mcvf;
}


/**
 * set whether re-randomize only one component
 *
 * @param b the boolean value
 */
void MixNormalCollection::setReRandOnlyOneComp(bool b) {
  MixNormal::reRandOnlyOneComp = b;
}


/**
 * set whether re-randomize on component dropping
 *
 * @param b the boolean value
 */
void MixNormalCollection::setNoReRandOnDrop(bool b) {
  MixNormal::noReRandOnDrop = b;
}


/**
 * set number of re-randomization per mixture component
 *
 * @param r the number of re-randomzation
 */
void MixNormalCollection::setReRandsPerMixCompRedux(unsigned r) {
  MixNormal::reRandsPerMixCompRedux = r;
}



/**
 * prepare to start an EM epoch
 */
void MixNormalCollection::startEpoch() {
  for ( MixNormal * p = _ftrMI; p != _ftrMI_endp; p++ )
    if ( p->active() )
      p->startEpoch();
}


/**
 * add more data set into the EM learning accumulation
 *
 * @param obsMat ObservationMatrix 
 * @param n_ftrs number of features
 * @param n_samps number of samples
 * @param rng the pointer to set of nodes
 * @return true if no componenet dropped
 */
void MixNormalCollection::addToEpoch(ObservationMatrix* obsMat,
				  size_t featureVecDim,
				  size_t totalNumFramesInSentence,
				  size_t numFramesToProcess,
				  unsigned firstFrameToProcess,
				  RangeSetCollection &tupleCol){
  
  assert(numFramesToProcess > 0);

  unsigned numFramesProcessed;
  RangeSet tuple;
  PointerSetToDataPoints pointerSet(featureVecDim);

  unsigned i = 0;
  for( MixNormal *p = _ftrMI; p != _ftrMI_endp; p++ ){
    if ( p->active() ) {
      tuple = tupleCol.rs[i];
      pointerSet.setDim(tuple.getSize());
      BUFFER_DATA_TYPE* obsMatPtr = (BUFFER_DATA_TYPE*) obsMat->baseAtFrame(0);
      numFramesProcessed = pointerSet.initialize(obsMatPtr, totalNumFramesInSentence, numFramesToProcess, firstFrameToProcess, tuple);
      if (numFramesProcessed != 0)
	if(p->_fullCoVar)
	  p->addToEpoch(pointerSet);
	else
	  p->addToEpochDiag(pointerSet);
    }
    ++i;
  }
}


/**
 * end an epoch of EM learning
 *
 * @return true if success
 */
void MixNormalCollection::endEpoch() {
  unsigned mcDrop = 0;

  for ( MixNormal * p = _ftrMI; p != _ftrMI_endp; p++ ) {
    if ( p->active() ) {
      if ( p->endEpoch() )
	mcDrop++;

      p->setDirty();
    }
  }
  if ( mcDrop > 0 )
    cout << "warning: " << mcDrop << " mixtures components dropped" << endl;
  fflush(stdout);
}



/**
 * recalculate the number of active learnings
 *
 * @param maxDist the maximum distance
 * @param aveDist the average distance
 * @param minDist the minimum distance
 * @param termDist the terminate distance
 */
int MixNormalCollection::reComputeNumActive(double &maxDist,
					    double &aveDist,
					    double &minDist,
					    const double termDist,
					    unsigned numEMIter) {
  maxDist = 0.0;
  aveDist = 0.0;
  minDist = HUGE;
  double tmp;

  int numSum = 0;
  int numActive = 0;

  for ( MixNormal *p = _ftrMI; p != _ftrMI_endp; p++ ) {
    if ( p->active() ) {
      tmp = p->llPercDiff();
      if ( tmp < termDist ) {
	p->reSetActive();
	p->setNumEMIters(numEMIter);
      }
      else
	numActive++;

      if ( tmp > maxDist )
	maxDist = tmp;
      if ( tmp < minDist )
	minDist = tmp;
      aveDist += tmp;
      numSum++;
    }
  }

  aveDist /= numSum;

  return numActive;
}





//////////////////// calcB //////////

/**
 * calculate the B matrix from the covariance obtained in the kmeans step
 * 
 */
void MixNormalCollection::calcB(){
   for ( MixNormal * p = _ftrMI; p != _ftrMI_endp; p++ ) {
     p->calcB();
     p->setDirty();
   }
}


//////////////////// checkB //////////

/**
 * calculate the B matrix from the covariance obtained in the kmeans step
 * 
 */
void MixNormalCollection::checkB(){
   for ( MixNormal * p = _ftrMI; p != _ftrMI_endp; p++ ) {
     p->checkB();
   }
}


