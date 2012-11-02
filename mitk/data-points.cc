#include "data-points.h"

/**
   * initializes the set of pointers to the start and end of the data points
   * @param obsMatrPtr matrix of dimensions totalNumFramesInSentence x skip
   * @param totalNumFramesInSentence 
   * @param numFramesToProcess 
   * @param firstFrameToProcess first frame to process in the matrix obsMatrPtr 
   * @param tuple tuple specification
   * @return number of processed frames == (end - start) / skip
   */
unsigned PointerSetToDataPoints::initialize(BUFFER_DATA_TYPE* obsMatPtr, 
				   size_t totalNumFramesInSentence,
				   size_t numFramesToProcess,
				   unsigned firstFrameToProcess,
				   RangeSet tuple) {
  //  unsigned startFrame,endFrame;

  // 1. Figure out the start and end frames given the tuple

  // lastFrame has to be iniatialized before firstFrameToProcess is modified
  unsigned lastFrame = firstFrameToProcess + numFramesToProcess - 1;
  int lastAllowableFrame = (int)totalNumFramesInSentence - tuple.max_lag - 1; // see below for meaning

  if(lastAllowableFrame < 0) return 0; // no frames are processed for this tuple for this sentence segment

  // if given the position of the first frame in the matrix, the
  // smallest lag (negative lag; if positive, there is no problem)
  // would point outside the matrix (negative index), make the first
  // frame be equal to the absolute value of the lag.
  if(tuple.min_lag < 0 && (-tuple.min_lag > (int) firstFrameToProcess) ) {
      firstFrameToProcess = - tuple.min_lag;
    }

  // same thing for the maximum lag (positive): adjust the last frame
  // to process to make sure the lag does not point to a position
  // outside the matrix bounds.
  if(tuple.max_lag > 0 && lastAllowableFrame < (int) lastFrame) {
    lastFrame = lastAllowableFrame;
  }
  
  //if(endFrame < startFrame) {
  if(lastFrame < firstFrameToProcess) {
    numSamples = 0;
    return numSamples;
  }

  
  // 2. Initialize pointers

  BUFFER_DATA_TYPE * pointerToStart = (BUFFER_DATA_TYPE*) (obsMatPtr +  firstFrameToProcess * skip);
  BUFFER_DATA_TYPE * pointerToEnd   = (BUFFER_DATA_TYPE*) (obsMatPtr +  lastFrame * skip);

  unsigned d1, d2;// d1: dimension of first part of the tuple, d2: dim of second part
  d1 = tuple.getSize(0);
  d2 = tuple.getSize(1);
  dim  = d1 + d2;  // initialize the dimension of the pointer vectors start and end

  // initialize first part of the tuple [0,d1]
  for( unsigned i = 0; i < d1; i++ ) {
    start[i] = pointerToStart + skip * tuple.set[0][i].lag + tuple.set[0][i].feat;
    end[i] = pointerToEnd + skip * tuple.set[0][i].lag + tuple.set[0][i].feat;
  }

  // second part of tuple (d1,d]
  for( unsigned i = 0; i < d2; i++ ) {
    start[i+d1] = pointerToStart + skip * tuple.set[1][i].lag + tuple.set[1][i].feat;
    end[i+d1] = pointerToEnd + skip * tuple.set[1][i].lag + tuple.set[1][i].feat;
  }

  numSamples = lastFrame - firstFrameToProcess + 1;
  return numSamples;

}

  /**
   * prints the datapoints
   * @param ofp the file handle to which data is printed out
   */
void PointerSetToDataPoints::print(FILE* ofp) {
  for(unsigned i = 0; i < numSamples; ++i) {
    for(unsigned j=0; j < dim; ++j) {
      fprintf(ofp,"%f ",*(start[j]));
      start[j] += skip;
    }
      fprintf(ofp,"\n");      
  }  
  
}
