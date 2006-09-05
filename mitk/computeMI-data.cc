#include <iostream>
#include "mixNormalCollection.h"

extern int sentPrintFrequency;

/**
 * compute the MI using the data.  Loops over sentence and uses tuples
 * as samples ( [x y] ) in the following computation:
 * I = 1/N log[ p(x,y) / (p(x)p(y)) ] 
 *
 * Preconditions: The paramaters (_means, _invVars, _b) have to be
 * properly initialized either by running kmeans, EM, both, or by
 * reading them from a file.
 *
 */
void MixNormalCollection::computeMIUsingData(ObservationMatrix& obsMat,
					     RangeSetCollection rangeSetCol,
					     Range &sentenceRange, 
					     const bool quiet,
					     FILE *outFileMI, Range &lrrng,int labpos,
					     FILE* rangeFileFP){
  size_t n_frames, n_samps;
  size_t n_ftrs = obsMat.numContinuous();
  int readStatus;
  unsigned frameStart,firstFrame;
  startEpochMI(rangeSetCol);
  for( Range::iterator sent = sentenceRange.begin(); !sent.at_end(); sent++ ){
    if( ! quiet )
      if( *sent % sentPrintFrequency == 0 )
	std::cout << "Processing sentence " << *sent << std::endl;      
    frameStart = 0;
    do {
      readStatus = readFeatures(sent, n_frames,
				n_samps,lrrng,labpos,
				frameStart,firstFrame);
      if(readStatus == NO_DATA) break;
      addToEpochMI(&obsMat, n_ftrs, 
				n_frames, n_samps,firstFrame,
				rangeSetCol);	
    } while(readStatus == DATA_LEFT);

  }
  endEpochMI(outFileMI,rangeFileFP);
}



//////////////////// MixNormalCollection::startEpochMI ////////////////////

/**
 *  starts the initialization of MI calculation by LLN using the
 *  actual data for each mixture 
 */
void MixNormalCollection::startEpochMI(const RangeSetCollection &tupleCol) {
  unsigned i, dX,scndSetSize;
  RangeSet tuple;

  i = 0;
  for( MixNormal *p = _ftrMI; p != _ftrMI_endp; p++ ){
    tuple = tupleCol.rs[i];
    dX = tuple.getSize(0);
    scndSetSize = tuple.getSize(1);
    if(scndSetSize == 0) //entropy
      p->startEpochEntropy();
    else // MI
      p->startEpochMI(dX);
    i++;
  }
}

//////////////////// MixNormalCollection::addToEpochMI ////////////////////

/**
 *  adds the data in the input buffer to MI calculation 
 *  for each mixture
 */
void MixNormalCollection::addToEpochMI(ObservationMatrix* obsMat,
				     size_t featureVecDim,
				     size_t totalNumFramesInSentence,
				     size_t numFramesToProcess,
				     unsigned firstFrameToProcess,
				     const RangeSetCollection &tupleCol){
  
  assert(numFramesToProcess > 0);
  unsigned numFramesProcessed;
  RangeSet tuple;
  PointerSetToDataPoints pointerSet(featureVecDim);

  unsigned i = 0;
  for( MixNormal *p = _ftrMI; p != _ftrMI_endp; p++ ){
    tuple = tupleCol.rs[i++];
    pointerSet.setDim(tuple.getSize());
    BUFFER_DATA_TYPE* obsMatPtr = (BUFFER_DATA_TYPE*) obsMat->features.ptr;
    numFramesProcessed = pointerSet.initialize(obsMatPtr, totalNumFramesInSentence, numFramesToProcess, firstFrameToProcess, tuple);
    if (numFramesProcessed != 0) {
      if(tuple.getSize(1) == 0) //entropy
	p->addToEpochEntropy(pointerSet);
      else //MI
	p->addToEpochMI(pointerSet);
    }
  }
}

//////////////////// MixNormalCollection::endEpochMI ////////////////////

/**
 *  finishes MI calculation for each mixture
 */
void MixNormalCollection::endEpochMI(FILE *outFileMI, FILE* rangeFileFP){
  unsigned i;
  PARAM_DATA_TYPE I, Hx, Hy, Hxy;
  char line[MAX_LINE_LEN];
  char *pos = line;

  i = 0;
  for( MixNormal *p = _ftrMI; p != _ftrMI_endp; p++ ){
    I = p->endEpochMI(Hx, Hy, Hxy); // When we calculate entropy Hy and Hxy will be zero.
    fgets(line,MAX_LINE_LEN,rangeFileFP); 
    // remove comments
    do {
      fgets(line,MAX_LINE_LEN,rangeFileFP); 
      if( ( pos = strchr(line,'#') ) != NULL ) *pos = '\0';  //remove comments
      else pos = line;
      while(isspace(*pos)) ++pos;  //eat spaces
    } while(*pos=='\0' || *pos =='\n');
    // --------------------
    if(line[strlen(line)-1] == '\n') line[strlen(line)-1] = '\0';
    fprintf(outFileMI,"%s %e\n", line, I);
    i++;
  }
  fflush(outFileMI);
}

//// end of MixNormalCollection MI routines ////////////


/// start of MixNormal MI routines ////////////////////

//////////////////// MixNormal::startEpochMI ////////////////////

/**
 * Initialize the necessary variables used in MI calculation by LLN
 * using "real" data.
 */
void MixNormal::startEpochMI(unsigned dX) {
  _dX  = dX;
  _Hx  = 0;
  _Hy  = 0;
  _Hxy = 0;
  _nSamplesMI = 0;
  marginalize(dX);
}

//////////////////// MixNormal::addToEpochMI ////////////////////

/**
 * Add nSamples data points to compute the mutual information 
 * between X and Y by using the actual samples. The input is an 
 * array of pointers to the data locations.
 *
 *  I = sum log( p(x,y)/p(x)p(y) / N 
 */
void MixNormal::addToEpochMI(PointerSetToDataPoints& ps) {
  PARAM_DATA_TYPE probX, probY, probXY; 

  if( ( _dX < 1 ) || ( _dX >= _numVariables ) || ps.numSamples < 1)
    error("ERROR:Inappropriate dX\n");
  cout<<"MixNormal::addEpochMI()\n";
  for( unsigned i = 0; i < ps.numSamples; i++){
    probX  = prob_x_GM(ps,i,0,_dX-1, _alphas, _meansX, _invCovsX, _invDetsX, _inv_pow_sqrt_2piX);
    probY  = prob_x_GM(ps,i,_dX,_numVariables-1, _alphas, _meansY, _invCovsY, _invDetsY, _inv_pow_sqrt_2piY);
    // There is redundancy here, get rid of on eofthe functions.
    //probXY = prob_x(ps, i);
    probXY  = prob_x_GM(ps,i, _alphas, _means, _invCovsXY, _invDets, _inv_pow_sqrt_2pi);
    //probXY  = prob_x_GM(ps,i, 0,_numVariables-1,_alphas, _means, _invCovsXY, _invDets, _inv_pow_sqrt_2pi);
    assert( (probX >= 0) && (probY >= 0) && (probXY >= 0) );
    //    assert(probXY >= probX*probY);
    if( (probX > 0) && (probY > 0) && (probXY > 0) ){
      _Hx  = _Hx  - log(probX);
      _Hy  = _Hy  - log(probY);
      _Hxy = _Hxy - log(probXY);
    }
  }
  _nSamplesMI = _nSamplesMI + ps.numSamples;
}

//////////////////// MixNormal::endEpochMI ////////////////////

/**
 * Finish MI calculation by normalizing with N and converting it to bits 
 * from nats.
 *
 *  I = sum log( p(x,y)/p(x)p(y) / N 
 */
PARAM_DATA_TYPE MixNormal::endEpochMI(PARAM_DATA_TYPE& Hx, PARAM_DATA_TYPE& Hy, PARAM_DATA_TYPE& Hxy) {

  PARAM_DATA_TYPE invNlog2 = 1.0 / (PARAM_DATA_TYPE)(log(2.0)*_nSamplesMI);

  Hx  = _Hx  * invNlog2;
  Hy  = _Hy  * invNlog2;
  Hxy = _Hxy * invNlog2;

  return Hx+Hy-Hxy;
}
