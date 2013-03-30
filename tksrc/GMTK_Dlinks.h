/*-
 * GMTK_DLINKS
 *      .h file the .cc file.
 *
 *  Written by Jeff Bilmes <bilmes@ee.washington.edu>
 * 
 *  $Header$
 * 
 * Copyright (C) 2001 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 *
 *
 */


#ifndef GMTK_DLINKS_H
#define GMTK_DLINKS_H

#include <vector>


#include "fileParser.h"
#include "logp.h"
#include "sArray.h"

#include "GMTK_NamedObject.h"

class DlinkMatrix;
class LinMeanCondDiagGaussian;

class Dlinks : public NamedObject {
  friend class DlinkMatrix;
  friend class LinMeanCondDiagGaussian;
  friend class GMParms;

  ///////////////////////////////////////////////////////
  // Structure for representing the dependency links to
  // other observation variables.
  struct Dlink {
    int lag;
    int offset;
  };

  ///////////////////////////////////////
  // an array of arrays of Dlinks. Each
  // member of the outer array corresponds to the dependencies
  // of one feature. The inner arrays point to the positions
  // in a feature vector of where the dependencies come from.
  vector< vector<Dlink> > dIndices;

  // use an sArray for speed
  sArray<int> preComputedOffsets;

  // object specific min/max lags
  int _minLag;
  int _maxLag;

  // global min/max lags and min/max link offsets in feature vectors.
  static int _globalMinLag;
  static int _globalMaxLag;
  static int _globalMinOffset;
  static int _globalMaxOffset;
  
  /////////////////////////////////////////////////////////
  // Support for a DlinkMatrix objects who share
  // this dlink structure. Rather than needing to compute
  // the xz and zz arrays multiple times, we compute
  // it only once here, and cache it according to the current
  // feature base.
  // --- cache support
  const Data32* arrayCacheTag;
  sArray <float> zzArrayCache;
  sArray <float> xzArrayCache;
  sArray <float> zArrayCache;
  void cacheArrays(const Data32* const base,
		   const float*const f);
  // --- precomputed length support
  unsigned _zzAccumulatorLength;
  /////////////////////////////////////////////////////////
  

public:

  ///////////////////////////////////////////////////////////  
  // General constructor
  Dlinks(); 

  ///////////////////////////////////////////////////////////  
  // dim: return the number of features for this 
  // collection of links corresponds to.
  int dim() { return dIndices.size(); }

  ///////////////////////////////////////////////////////////  
  // return total number of dlinks refered to by this.
  unsigned totalNumberLinks() { return (unsigned) preComputedOffsets.len(); }

  ///////////////////////////////////////////////////////////  
  // For EM, return the length of the z'z accumulators
  unsigned zzAccumulatorLength() { return _zzAccumulatorLength; }

  ///////////////////////////////////////////////////////////  
  // numLinks: return the number of links for the ith
  // feature.
  int numLinks(const int i) { 
    assert ( i >= 0 && i < dim() );
    return dIndices[i].size(); 
  }


  ////////////////////////////////////////////////
  // precomputes the offset array 
  void preCompute(const unsigned stride);

  ////////////////////////////////////////////////
  // clear the internal cache
  void clearArrayCache();

  ///////////////////////////////////////////////////////////    
  // read in the basic parameters, assuming file pointer 
  // is located at the correct position.
  void read(iDataStreamFile& is);

  ///////////////////////////////////////////////////////////    
  // write out the basic parameters, starting at the current
  // file position.
  void write(oDataStreamFile& os);

  // return the min and maximum lags respectively.
  int minLag() { return _minLag; }
  int maxLag() { return _maxLag; }

  static int globalMinLag() { return _globalMinLag; }
  static int globalMaxLag() { return _globalMaxLag; }

};


#endif

