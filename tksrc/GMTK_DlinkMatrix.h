/*-
 * GMTK_DlinkMatrix.h
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


#ifndef GMTK_DLINKMATRIX_H
#define GMTK_DLINKMATRIX_H

#include "fileParser.h"
#include "logp.h"
#include "sArray.h"

#include "GMTK_EMable.h"
#include "GMTK_NamedObject.h"
#include "GMTK_RandomVariable.h"

class Dlinks;
class LinMeanCondDiagGaussian;


class DlinkMatrix : public EMable, public NamedObject  {

  friend class Dlinks;
  friend class LinMeanCondDiagGaussian;

  ///////////////////////////////////////////////////////////  
  // the number of links per feature
  sArray < int > _numLinks;

  //////////////////////////////////
  // The acutal matrix data values, packed
  // into one 1D array
  sArray< float > arr;
  ///////////////////////////////////////////////////////////  

  //////////////////////////////////
  // Data structures support for EM
  sArray< float > nextArr;

public:

  ///////////////////////////////////////////////////////////  
  // General constructor
  DlinkMatrix();
  ~DlinkMatrix() { }

  //////////////////////////////////
  // set all current parameters to random values
  void makeRandom();

  void makeUniform();

  //////////////////////////////////////////////
  // read/write basic parameters
  void read(iDataStreamFile& is);
  void write(oDataStreamFile& os);

  ///////////////////////////////////////////////////////////  
  // num number of features (dimensionality) for this 
  int dim() { return _numLinks.len(); }

  ///////////////////////////////////////////////////////////  
  // numLinks: return the number of links for the ith
  // feature.
  int numLinks(const int i) { 
    assert ( i >=0 && i < dim() );
    return _numLinks[i];
  }

  bool compatibleWith(Dlinks&d);

  //////////////////////////////////
  // Public interface support for EM
  //////////////////////////////////
  void emStartIteration() {}
  void emIncrement(RandomVariable*,logpr prob) {}
  void emEndIteration() {}
  void emSwapCurAndNew() {}
  void emStoreAccumulators(oDataStreamFile& ofile) {}
  void emLoadAccumulators(iDataStreamFile& ifile) {}
  void emAccumulateAccumulators(iDataStreamFile& ifile) {}
  //////////////////////////////////


};



#endif
