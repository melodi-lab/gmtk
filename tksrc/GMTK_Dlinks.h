/*-
 * GMTK_DLINKS
 *      .h file the .cc file.
 *
 *  Written by Jeff Bilmes <bilmes@ee.washington.edu>
 * 
 *  $Header$
 * 
 * Copyright (c) 2001, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about
 * the suitability of this software for any purpose.  It is provided
 * "as is" without express or implied warranty.
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

  unsigned totalNumberLinks() { return (unsigned) preComputedOffsets.len(); }

  int _minLag;
  int _maxLag;

  
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
  sArray <float> zArray;
  void cacheArrays(const Data32* const base,
		   const float*const f);
  // --- precomputed length support
  unsigned zzAccumulatorLength;
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
  // numLinks: return the number of links for the ith
  // feature.
  int numLinks(const int i) { 
    assert ( i >= 0 && i < dim() );
    return dIndices[i].size(); 
  }

  /////////////////////////////////////////////////////
  // returns true if objects are compatible
  bool compatibleWith(DlinkMatrix&d);

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

};


#endif

