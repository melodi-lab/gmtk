/*-
 * GMTK_VariableParams.h
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


#ifndef GMTK_VARIABLEPARAMS_H
#define GMTK_VARIABLEPARAMS_H


#include "fileParser.h"
#include "logp.h"
#include "sArray.h"

class VariableParams {

  /////////////////////////////  
  // Collection of prob. mass functions
  sArray< Discrete1DPDF* > pmfs;

  /////////////////////////////
  // Collection of means
  sArray< MeanVector* > means;

  ////////////////////////////////
  // Collection of diag. covariances
  sArray< DiagCovarVector* > covars;

  ////////////////////////////////
  // Collection of DLINKS
  sArray< Dlinks* > dlinks;

  ////////////////////////////////
  // Collection of packed sparse matrices
  sArray< PackedSparseRealMatrix* > sparseMats;

  ////////////////////////////////
  // Collection of 2D Dense matrices
  sArray< RealMatrix* > denseMats;

  //////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////
  // Collection of diag. covariance Gaussians
  sArray< DiagGaussian* > diagGaussians;

  ////////////////////////////////
  // Collection of DLINKS
  sArray< Dlinks* > dlinks;

  ////////////////////////////////
  // Collection of diagonal covariance Gaussians with linear mean 
  // dependency links (these cover the case full Covariance 
  // Gaussians, plus other forms such as banded or block diagonal, 
  // factored sparse inverse covariances, and so on.
  sArray< LinMeanCondDiagGaussian* > linMeanCondGaussians;

  ////////////////////////////////
  // Collection of diagonal covariance Gaussians with linear and/or 
  // non-linear mean dependency links (these cover the case 
  // of "non-linear" Gaussians, and so on).
  sArray< NLinMeanCondDiagGaussian* > nLinMeanCondGaussians;

  //////////////////////////////////////////////////////////////////


public:

  ///////////////////////////////////////////////////////////  
  // General constructor
  VariableParams(); 

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
