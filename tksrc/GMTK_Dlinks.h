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

class Dlinks : public NamedObject {


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

public:

  ///////////////////////////////////////////////////////////  
  // General constructor
  Dlinks(); 

  ///////////////////////////////////////////////////////////  
  // numFeats: return the number of features for this 
  // collection of links corresponds to.
  int numFeats() { return dIndices.size(); }

  ///////////////////////////////////////////////////////////  
  // numLinks: return the number of links for the ith
  // feature.
  int numLinks(const int i) { 
    assert ( i >= 0 && i < numFeats() );
    return dIndices[i].size(); 
  }


  ///////////////////////////////////////////////////////////    
  // read in the basic parameters, assuming file pointer 
  // is located at the correct position.
  void read(iDataStreamFile& is);

  ///////////////////////////////////////////////////////////    
  // write out the basic parameters, starting at the current
  // file position.
  void write(oDataStreamFile& os);


};


#endif
