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


#ifndef GMTK_DLINKS
#define GMTK_DLINKS


#include "fileParser.h"
#include "logp.h"
#include "sArray.h"

class Dlinks {


  ///////////////////////////////////////////////////////
  // Structure for representing the dependency links to
  // other observation variables.
  struct Dlink {
    int lag;
    int offset;
  };
  
  sArray< sArray<Dlink> > dIndices;

public:

  ///////////////////////////////////////////////////////////  
  // General constructor
  Dlinks() 

  ///////////////////////////////////////////////////////////  
  // numFeats: return the number of features for this 
  // collection of links corresponds to.
  int numFeats() { return dIndices.len(); }

  ///////////////////////////////////////////////////////////  
  // numLinks: return the number of links for the ith
  // feature.
  int numLinks(const int i) { 
    assert ( i >= 0 && i < numFeats() );
    return dIndices[i].len(); 
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


#endif // defined GMTK_CPT
