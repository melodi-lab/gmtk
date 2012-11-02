/*-
 * GMTK_Weightmatrix.h
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


#ifndef GMTK_DIRICHLETTABLE_H
#define GMTK_DIRICHLETTABLE_H

#include "fileParser.h"
#include "logp.h"
#include "sArray.h"

#include "GMTK_NamedObject.h"

class RV;

class DirichletTable : public NamedObject  {

  ///////////////////////////////////////////////////////////  
  // The data values
  sArray < logpr > table;
  ///////////////////////////////////////////////////////////  


  // This table of counts corresponds to a multi-dimensional array,
  // where the table is (cardinalities.size()+1)-D, and where the
  // first cardinalities.size() dimension lengths are contained in
  // cardinalities[j], for dimension , and where the final dimension
  // length is of value _card.

  unsigned _numDimensions;
  // the dimension lengths of all but the last dimension
  vector < unsigned > dimensionLengths;

  // the table must corresponds to a Dense CPT, so to relate
  // this to that, we have for an MDCPT:
  // DirichletTable::_numDimensions == MDCPT::_numParents + 1
  // DirichletTable::dimensionLengths.size() = MDCPT::cardinalities.size()+1
  // For i = 1 ... MDCPT::cardinalities.size()
  //     DirichletTable::dimensionLengths[i] = MDCPT::cardinalities[i]
  // For i = MDCPT::cardinalities.size()+1
  //     DirichletTable::dimensionLengths[i] = MDCPT::_card
 
  

public:

  ///////////////////////////////////////////////////////////  
  // General constructor
  DirichletTable();
  ~DirichletTable() {}


  unsigned numDimensions() {
    return _numDimensions;
  }
  unsigned dimension(const unsigned int par) {
    assert ( par < _numDimensions );
    return dimensionLengths[par];
  }
  unsigned lastDimension() {
    return dimensionLengths[dimensionLengths.size()-1];
  }
  
  unsigned tableSize() {
    return table.size();
  }

  logpr tableValue(unsigned i) 
  { 
    return table.ptr[i]; 
  }

  ///////////////////////////////////////////////////////////    
  void read(iDataStreamFile& is);

  void write(oDataStreamFile& os);

  const string typeName() { return "Dirichlet Table"; }
  //////////////////////////////////

};



#endif 
