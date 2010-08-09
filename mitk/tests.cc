#include "global-parameters.h"
#include "general.h"
#include "error.h"
//#include "Range.H"
#include "mixNormal.h"
#include "mixNormalCollection.h"
#include "readRange.h"
#include "GMTK_ObservationMatrix.h"

void generateSyntheticDataFromLearnedParams(MixNormalCollection* mg, int numSamples, FILE* pi_fp) {
  
  FILE* ofp2=fopen("SYNTHETIC_DATA_WITH_COV.OUT","w");
  FILE* ofp3=fopen("SYNTHETIC_DATA_WITH_B.OUT","w");
  if(ofp2 == NULL) {
    fprintf(stderr,"Could not open SYNTHETIC_DATA.OUT file for writing\n");
    exit(-1);
  }
  if (pi_fp == NULL || fsize(pi_fp) <= 0) 
    mg->generateDataUsingCov(ofp2,10000, 0);
  mg->generateData(ofp3,10000, 0);
  
}

void dumpDistribSampleData(FILE* ofp, 
			   ObservationMatrix * obsMat,
			   RangeSetCollection &tupleCol,
			   Range &lrrng,
			   Range &sentRange,
			   unsigned numMixtures,
			   unsigned maxIter,
			   int labpos,
			   unsigned tupleNum,
			   const bool quiet) {
  
  const size_t featureVecDim = obsMat->numFeatures();
  RangeSet tuple = tupleCol.rs[tupleNum];
  PointerSetToDataPoints pointerSet(featureVecDim,tuple.getSize());
  size_t totalNumFramesInSentence, numFramesToProcess;
  int readStatus;
  unsigned frameStart,firstFrameToProcess;
  unsigned numFramesProcessed;

  
  for( Range::iterator sent = sentRange.begin(); !sent.at_end(); sent++ ){
    //if( ! quiet )
    //if( *sent % FREQUENCY == 0 )
    //	cout << "Processing sentence no " << *sent <<endl;   
    
    frameStart = 0;
    do {
      readStatus = readFeatures(sent, totalNumFramesInSentence,numFramesToProcess,lrrng,labpos,frameStart,firstFrameToProcess);
      //cout << "Processing sentence no " << *sent <<endl;
      if(readStatus == NO_DATA) break;
  
      
      //BUFFER_DATA_TYPE* obsMatPtr = (BUFFER_DATA_TYPE*) obsMat->features.ptr;
      BUFFER_DATA_TYPE* obsMatPtr = (BUFFER_DATA_TYPE*) obsMat->baseAtFrame(0);
      numFramesProcessed = pointerSet.initialize(obsMatPtr, totalNumFramesInSentence, numFramesToProcess, firstFrameToProcess, tuple);
      //cout << "Processing sentence no " << *sent <<endl;
      pointerSet.print(ofp);
      //cout << "Processing sentence no " << *sent <<endl;
    } while(readStatus == DATA_LEFT);
#if DEBUG
    fprintf(stderr,"First 2 frames of input file:\n");
    obsMat->printFrame(stderr,0);
    obsMat->printFrame(stderr,1);
#endif    
  }

}


void MixNormalCollection::dumpParameters(FILE* ofp, 
					 unsigned tupleNum) {
  MixNormal * p;
  unsigned i;
  for (  p = _ftrMI, i = 0; (p != _ftrMI_endp && i < tupleNum); p++ );
  if(p == _ftrMI_endp) {
    fprintf(stderr,"Tuple number larger than number of tuples\n");
    return;
  }
  //p->printCurParams(ofp);
  p->printMeans(ofp);
}


void MixNormalCollection::dumpCurIterParams(int iter){
  MixNormal *ftr_mi_p;    

  
  char mgFilename[50];
 
  int mixNum=0;
  ftr_mi_p = _ftrMI;
  while (ftr_mi_p != _ftrMI_endp) {
    FILE* fp;
    sprintf(mgFilename,"mg.%d.%d",mixNum,iter);
    if( (fp = fopen(mgFilename,"w"))==NULL) {
      error("Could not open output mg file in MixNormalCollection::dumpCurIterParams (tests.cc)\n");
    }
    ftr_mi_p->printCurParams(fp);
    ftr_mi_p++;
    mixNum++;
    fclose(fp);
  }
}


//////////////////// printMeans ////////////////////

/**
 * write to the c-like files
 *
 * @param fp the output file pointer
 * @exception ErrorWritingException error writing to file
 */
void MixNormal::printMeans(FILE *fp) const {
  unsigned i, l;

  for ( l = 0; l < _numMixtures; l++ ) {

    for ( i = 0; i < _numVariables; i++ )
      fprintf(fp, "%f ", *(_means+l*_numVariables+i));
    fprintf(fp, "\n");
  }
} // end printMeans
