#ifndef _DATA_POINTS_
#define _DATA_POINTS_

#include "global-parameters.h"
#include "GMTK_ObservationMatrix.h"
#include "readRange.h"

/**
 * set of pointers to the data points for which the joint prob distribution is to be estimated

*/

class PointerSetToDataPoints {

 public:

  unsigned dim;  // dimension of the vectors start and end 

  BUFFER_DATA_TYPE * start[MAX_POINTER_SET_DIM];  // array of pointers to the first set of data points
  BUFFER_DATA_TYPE * end[MAX_POINTER_SET_DIM];    // array of pointers to the last set of data points
  unsigned skip; // by how much to advance to get the next set of data
                   // points.  Is equal to the dimension of the
                   // *feature* vectors , whcih is not equl to the
                   // dimension of the vectors start and end.

  unsigned numSamples;  // can be infered from (end - start) / (sizeof(BUFFER_DAT_TYPE) but his is "cleaner"

  PointerSetToDataPoints(unsigned featureVecDim) {
    skip = featureVecDim;
  }

  PointerSetToDataPoints(unsigned featureVecDim, unsigned dim) {
    skip = featureVecDim;
    this->dim = dim;
  }
  
  void setDim(unsigned dim) { this->dim = dim;}

  /**
   * initializes the set of pointers to the start and end of the data points
   * @param obsMatrPtr matrix of dimensions totalNumFramesInSentence x skip
   * @param totalNumFramesInSentence 
   * @param numFramesToProcess 
   * @param firstFrameToProcess first frame to process in the matrix obsMatrPtr 
   * @param tuple tuple specification
   * @return number of processed frames == (end - start) / skip
   */
  
  unsigned initialize(BUFFER_DATA_TYPE* obsMatPtr, 
		      size_t totalNumFramesInSentence,
		      size_t numFramesToProcess,
		      unsigned firstFrameToProcess,
		      RangeSet tuple);
  /**
   * prints the datapoints
   * @param ofp the file handle to which data is printed out
   */
  void print(FILE* ofp);

};


#endif
