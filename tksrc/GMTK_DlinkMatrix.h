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
#include "GMTK_Dlinks.h"

class LinMeanCondDiagGaussian;

class DlinkMatrix : public EMable, public NamedObject  {
  friend class Dlinks;
  friend class LinMeanCondDiagGaussian;

  ///////////////////////////////////////////////////////
  // The actual dlink structure
  Dlinks* dLinks;

  ///////////////////////////////////////////////////////////
  // The acutal matrix data values, packed
  // into one 1D array
  sArray< float > arr;
  ///////////////////////////////////////////////////////////  

  ///////////////////////////////////////////////////////////
  // Data structures support for EM
  sArray< float > nextArr;

  /////////////////////////////////////////////////
  // counts the number of gaussian components
  // that are sharing this mean.
  unsigned refCount;

public:

  ///////////////////////////////////////////////////////////
  // General constructor
  DlinkMatrix();
  ~DlinkMatrix() { }


  // When noisy cloning an object, this gives
  // the fraction to multiply to get the STD of the noise.
  static double cloneSTDfrac;
  static void checkForValidValues();

  //////////////////////////////////
  // set all current parameters to random/uniform values
  void makeRandom();
  void makeUniform();

  //////////////////////////////////////////////
  // read/write basic parameters
  void read(iDataStreamFile& is);
  void write(oDataStreamFile& os);

  // create a copy of self, but with slightly perturbed
  // means values.
  DlinkMatrix* noisyClone();

  ///////////////////////////////////////////////////////////  
  // num number of features (dimensionality) for this 
  int dim() { return dLinks->dim(); }

  ///////////////////////////////////////////////////////////  
  // numLinks: return the number of links for the ith
  // feature.
  int numLinks(const int i) { 
    assert ( i >=0 && i < dim() );
    return dLinks->numLinks(i);
  }

  //////////////////////////////////
  // Public interface support for EM
  //////////////////////////////////
  void emStartIteration(sArray<float>& xzAccumulators,
			sArray<float>& zzAccumulators,
			sArray<float>& zAccumulators);
  void emIncrement(const logpr prob,
		   const float fprob,
		   const float* const f,
		   const Data32* const base,
		   const int stride,
		   float* xzAccumulators,
		   float* zzAccumulators,
		   float* zAccumulators);
  void emEndIteration(const float*const xzAccumulators);
  void emSwapCurAndNew();
  void emStoreAccumulators(oDataStreamFile& ofile);
  void emLoadAccumulators(iDataStreamFile& ifile);
  void emAccumulateAccumulators(iDataStreamFile& ifile);
  //////////////////////////////////


};



#endif
