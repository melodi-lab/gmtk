#include <math.h>
#include <sys/time.h>

#include "mixNormalCollection.h"
#include "matrix-ops.h"


int readFeatures(Range::iterator krit, size_t &n_frames,
		 size_t &n_samps,
		 Range &lrrng, int labpos,
		 unsigned &frameStart,unsigned &firstFrame);


extern int sentPrintFrequency;

/////////////////////////// MixNormalCollection::kmeans /////////
/** 
 * runs kMeans: iterates over sentences and accumulates sufficient
 * statistics.  At the end of each epoch, kMeans parameters are updated.
 * 
 * */
void MixNormalCollection::kmeans(ObservationMatrix* obs_mat,
				 RangeSetCollection rangeSetCol,
				 Range &lrrng,
				 Range &kMeansRange,
				 unsigned numMixtures,
				 unsigned maxIter,
				 int labpos,
				 const bool quiet){
  unsigned i, nullCluster;
  const size_t n_ftrs = obs_mat->numFeatures();
  size_t n_frames, n_samps;
  bool estCov, randLabel;
  int readStatus;
  unsigned frameStart,firstFrame;
  
  DBGFPRINTF((stderr,"In MixNormalCollection::kmeans. Will perform %d iterations of kmeans.\n",maxIter));

  i = 0;
  while( i < maxIter ){
    if( i == 0 ) randLabel = true; else randLabel = false;
    if( i == maxIter-1) estCov = true; else estCov = false;
    if(!quiet) cout << "k-means iteration " << i << endl;
    DBGFPRINTF((stderr,"In MixNormalCollection::kmeans.  Starting iteration %d of %d\n",i,maxIter));
    startKMeansEpoch(randLabel, estCov);
    for( Range::iterator sent = kMeansRange.begin(); !sent.at_end(); sent++ ){
      if( ! quiet )
	if( *sent % sentPrintFrequency == 0 )
	  std::cout << "Processing sentence no " << *sent << std::endl;   

      frameStart = 0;
      do {
	readStatus = readFeatures(sent, n_frames,n_samps,lrrng,labpos,frameStart,firstFrame);
	if(readStatus == NO_DATA) break;

	addToKMeansEpoch(obs_mat, n_ftrs, n_frames, n_samps, firstFrame, rangeSetCol);
      } while(readStatus == DATA_LEFT);
    }
    int rangeSpecNum;
    if( noSamplesFound(rangeSpecNum) ) {
      error("ERROR: There were no samples for at least one range spec (the %d th one).  Possible causes:  the label provided does not exist in the label file or there are too few frames with that label.\n",rangeSpecNum);
    }
    
    nullCluster = endKMeansEpoch();
    i++;
    if( nullCluster > 0) {
      cout <<nullCluster <<" clusters were null.  They were assigned the sample mean and covariance"<< endl;
    }
  } // end 

}



/////////////////////////// MixNormalCollection::startKMeansEpoch /////////
/** 
 * calls MixNormal:startKMeansEpoch for each n-tuplet kmeans 
 * 
 *
 */
void MixNormalCollection::startKMeansEpoch(bool randLabel, bool estCov){

  for( MixNormal * p = _ftrMI; p != _ftrMI_endp; p++ ) 
    p->startKMeansEpoch(randLabel, estCov);
}


/////////////////////////// MixNormalCollection::addToKMeansEpoch /////////
/** 
 * calls MixNormal::addEpoch for each n-tuplet kmeans 
 * 
 *
 */
void MixNormalCollection::addToKMeansEpoch(ObservationMatrix* obs_mat,
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
    tuple = tupleCol.rs[i++];
    pointerSet.setDim(tuple.getSize());
    BUFFER_DATA_TYPE* obs_mat_ptr = (BUFFER_DATA_TYPE*) obs_mat->baseAtFrame(0);
    numFramesProcessed = pointerSet.initialize(obs_mat_ptr, totalNumFramesInSentence, numFramesToProcess, firstFrameToProcess, tuple);
    if (numFramesProcessed > 0)  
      p->addToKMeansEpoch(pointerSet);
  }

}

/////////////////////////// MixNormalCollection::endKMeansEpoch /////////
/** 
 * calls endEpoch for each n-tuplet kmeans 
 * 
 * @return the number of empty clusters
 */
unsigned MixNormalCollection::endKMeansEpoch(){
  unsigned nullCluster = 0;
  for( MixNormal * p = _ftrMI; p != _ftrMI_endp; p++ ) { 
    nullCluster += p->endKMeansEpoch();
    p->setDirty();
  }
  fflush(stdout);

  return nullCluster;
}



/////////////////////// MixNormal::nearestCluster //////////////

/**
 * find the closest clusters to a given feature vector
 *
 */
unsigned MixNormal::nearestCluster(const PointerSetToDataPoints &pointerSet, unsigned offset) const {

  unsigned iLabel;
  double diff, dist,distNew;
  unsigned dim = pointerSet.dim;

  iLabel = 0;
  dist=0.0;
  for(unsigned j = 0; j < dim; ++j) {
    diff = *(pointerSet.start[j] + offset*pointerSet.skip ) - *(_means + iLabel*dim + j);
    dist += (diff * diff);
  }
  for( unsigned i = 1; i < _k; i++ ){
    distNew=0.0;
    for(unsigned j = 0; j < dim; ++j) {
      diff = *(pointerSet.start[j] + offset*pointerSet.skip ) - *(_means + i*dim + j);
      distNew += (diff * diff);
      if(distNew >= dist) goto nextLabel; 
    }
    iLabel = i;
    dist = distNew;
  nextLabel:;
  }

  return iLabel;
}

 /////////////////////////////// MixNormal::startKMeansEpoch ////////////////////////

/**
 * initializes data structures before performing one epoch of kMeans
 *  
 */
void MixNormal::startKMeansEpoch(bool randLabel, bool estCov){ 
  unsigned i;

  _eqWeight  = 1.0 / _k;
  _randLabel = randLabel;
  _estCov    = estCov;
  _numAccum = 0; //  _accumData = 0;
  for( i = 0; i < _k; i++){
    *(_alphas+i) = 0;
    for(unsigned j = 0; j < _d; ++j)  *(_accumMean + i*_d + j) = 0.0;
    if( _estCov )
      for(unsigned j = 0; j < _d*_d; ++j)  *(_cov + i*_d*_d + j) = 0.0;
  }

}

/////////////////////////////// MixNormal::addToKMeansEpoch ////////////////////////

/**
 * accumulate data into the means and the covariances (if it is the
 * last epoch and we want to estimate the cov)
 */
void MixNormal::addToKMeansEpoch(const PointerSetToDataPoints &pointerSet){

  unsigned i, iLabel;
  unsigned dim = pointerSet.dim;

  assert(_d == dim);
  
  for( i = 0; i < pointerSet.numSamples; i++){
    if( _randLabel ) {
      iLabel = rand() % _k;
    }
    else {
      iLabel = nearestCluster(pointerSet, i);
    }
    (*(_alphas+iLabel))++;//_alpha[iLabel]++;
    for(unsigned j = 0; j < dim; ++j) {
      *(_accumMean + iLabel*dim +j) +=  *(pointerSet.start[j] + i*pointerSet.skip );
    }
    if( _estCov ) {
      for(unsigned j = 0; j < dim; ++j) {  
	for(unsigned k = 0; k < dim; ++k) {  
	  *(_cov + iLabel*dim*dim + j*dim + k) += *(pointerSet.start[j] + i*pointerSet.skip ) * *(pointerSet.start[k] + i*pointerSet.skip );
	}
      }
    }
  }
  // Accumulate the number of samples so far
  _numAccum += pointerSet.numSamples;    //    _accumData += pointerSet.numSamples;

  
}
  
/////////////////////////////// MixNormal::endKMeansEpoch ////////////////////////

/**
 * finish computing the mean and cov from the accumulated quantities using kMeans 
 *
 */
unsigned MixNormal::endKMeansEpoch(){

  unsigned i, numNullClusters = 0;
  //bool performedSampleMeanCovCalc = false;
  double* sampleMean;
  double* sampleCov;
  bool doNotDropComp = true; // true if we do not drop null clusters (we use the sample mean and cov instead)
  double deter;
  // To make to this sample Mean and Cov more efficient, only compute
  // them the first time they are needed.  This way, we save
  // computation when there are no drops.

  if(doNotDropComp) {
    sampleMean = new double [_d];
    sampleCov =  new double [_d*_d];
    // zero mean and cov
    for( unsigned j = 0; j < _d;    ++j) *(sampleMean + j)  = 0.0;
    if( _estCov ) for( unsigned j = 0; j < _d*_d; ++j) *(sampleCov  + j)  = 0.0;

    for( unsigned clusterCounter = 0; clusterCounter < _k; clusterCounter++ ) {
      for( unsigned j = 0; j < _d; ++j) *(sampleMean + j)  += *(_accumMean + clusterCounter*_d + j);
      if( _estCov ) {
	for(unsigned j = 0; j < _d; ++j)
	  for(unsigned k = 0; k < _d; ++k)
	    *(sampleCov + j*_d + k) += *(_cov + clusterCounter*_d*_d + j*_d + k);

      }
    }

    for( unsigned j = 0; j < _d; ++j) *(sampleMean + j) /= _numAccum;
    if(_estCov)  {
      double tmpAccumData = (double) _numAccum;

      calcCov(sampleCov , sampleMean, &tmpAccumData);

      if( (deter=determinant(sampleCov,_d)) <= 0 ) {
	cout<<"Determinant of sample covariance matrix of the component of n-tuple # "<<_index<<" is negative\n";
	cout<<"determinant = "<<deter<<". Cov = "<<endl;
	for( unsigned j = 0; j < _d*_d; ++j) printf("%f ",*(sampleCov  + j));
	cout<<endl;
	exit(-1);
      }
    }
  }

  for( i = 0; i < _k; i++){
    if( *(_alphas+i) > 0 ){
      // find mean: divide each entry of the accumulated mean vector by
      // alpha (which hold the number of vectors in the cluster)
      for(unsigned j = 0; j < _d; ++j)
	*(_means +i*_d + j)  = *(_accumMean +i*_d +j) / *(_alphas+i);

      // calculate the covariance
      if( _estCov ) {
	calcCov(_cov + i*_d*_d, _means + i*_d, _alphas + i);
      }
      *(_alphas+i) /= _numAccum;
    }
    else{
      if(doNotDropComp) {
	// assign sample mean and cov to the current mean and cov (i-th index) 
	*(_alphas+i) = _eqWeight;
	
	for( unsigned j = 0; j < _d;    ++j) *(_means + i*_d    + j) = *(sampleMean + j);
	if(_estCov) for( unsigned j = 0; j < _d*_d; ++j) *(_cov  + i*_d*_d + j) = *(sampleCov  + j);
      }
      numNullClusters++;
    }
    if( _estCov ) if( (deter = determinant(_cov,_d) ) <= 0 ) {
      cout<<"Determinant ("<<deter<<") of covariance matrix of "<<i<<"th component of n-tuple # "<<_index<<" is negative\n";
      cout<<"determinant = "<<deter<<". Cov = "<<endl;
      for( unsigned j = 0; j < _d*_d; ++j) printf("%f ",*(_cov +i*_d*_d + j));
      cout<<endl;
      cout<<"Replacing the covariance with the sample covariance\n";
      for( unsigned j = 0; j < _d*_d; ++j) *(_cov  + i*_d*_d + j) = *(sampleCov  + j);
    }
  }
  // normalize the alphas in case there were component drops
  kMeansNormalize();
  
  return numNullClusters;
}

/**
 * calculate the covariance
 * COV: _d x _d matrix, MEAN: _d column vector  
 * COV = (1/_alpha)*COV - MEAN * transpose(MEAN)
 */
void MixNormal::calcCov(double *cov, double * mean, double* alpha) {
  for(unsigned j = 0; j < _d; ++j)
    for(unsigned k = 0; k < _d; ++k)
      *(cov + j*_d + k)  = ( *(cov + j*_d + k) / (*alpha) ) - ( *(mean + j) * *(mean + k) );
}

/** 
 * Normalize alphas.  Needed when there is a component drop.
 */
void MixNormal::kMeansNormalize(){
  unsigned i;
  double sum = 0;

for( i = 0; i < _k; i++ ) sum += *(_alphas+i);
  for( i = 0; i < _k; i++ ) *(_alphas+i) /= sum;
}

/**
 * auxiliary function to compute the sample mean and covariance.
 * Called by endKMeansEpoch. 
 *
 */
void MixNormal::calcSampleMeanCov(double* sampleMean, double *sampleCov) {
  double deter;
  // zero mean and cov
  for( unsigned j = 0; j < _d;    ++j) *(sampleMean + j)  = 0.0;
  if( _estCov ) for( unsigned j = 0; j < _d*_d; ++j) *(sampleCov  + j)  = 0.0;
  
  for( unsigned clusterCounter = 0; clusterCounter < _k; clusterCounter++ ) {
    for( unsigned j = 0; j < _d; ++j) *(sampleMean + j)  += *(_accumMean + clusterCounter*_d + j);
    if( _estCov ) {
      for(unsigned j = 0; j < _d; ++j)
	for(unsigned k = 0; k < _d; ++k)
	  *(sampleCov + j*_d + k) += *(_cov + clusterCounter*_d*_d + j*_d + k);
      
    }
  }

  for( unsigned j = 0; j < _d; ++j) *(sampleMean + j) /= _numAccum;
  if(_estCov)  {
    double tmpAccumData = (double) _numAccum;
    
    calcCov(sampleCov , sampleMean, &tmpAccumData);
    
    if( (deter=determinant(sampleCov,_d)) <= 0 ) {
      cout<<"Determinant of sample covariance matrix of the component of n-tuple # "<<_index<<" is negative\n";
      cout<<"determinant = "<<deter<<". Cov = "<<endl;
      for( unsigned j = 0; j < _d*_d; ++j) printf("%f ",*(sampleCov  + j));
      cout<<endl;
      exit(-1);
    }
  }
  
}
