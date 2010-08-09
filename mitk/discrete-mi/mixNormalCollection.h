/**
 *: mixnormalcollection.h
 */


#ifndef MIX_NORMAL_COLLECTION_H
#define MIX_NORMAL_COLLECTION_H

#include <cstdio>


#include "mixNormal.h"
#include "range.h"
#include "readRange.h"
#include "GMTK_ObservationMatrix.h"

/**
 * collection of MixNormal objects
 */
class MixNormalCollection {

public:
  // constructors and destructor
  MixNormalCollection(RangeSetCollection rangeSetCol);
  ~MixNormalCollection(){delete [] _ftrMI;}

  // interface to EM learning
  void startEpoch();
  bool addToEpoch(ObservationMatrix* obsMat, const size_t n_ftrs, 
		  const size_t n_frames,
		  const size_t n_samps, 
		  unsigned firstFrame,
		  const RangeSetCollection &rngCol);
  void endEpoch(const size_t n_samps, const RangeSetCollection &rngCol,
		FILE*ofp);

private:
  const unsigned _numMis;  
  const unsigned _numVars;
  const unsigned _dimX;
  MixNormal * const _ftrMI;
  MixNormal * const _ftrMI_endp;
};

#endif





