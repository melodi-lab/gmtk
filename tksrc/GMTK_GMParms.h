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


  //////////////////////////////////////////////////////////////////
  // Basic shared low-level parameters: These are the objects that 
  //  higher level objects (see below) might share together.
  //  All of these objects are "EMable" in the sense that
  //  they may be trained using EM (plus possibly some other gradient
  //  based training method).
  //////////////////////////////////////////////////////////////////

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
  // Collection of objects
  // used for linear dependencies via
  // a dlink topology structure.
  sArray< DlinkMatrix* > dLinkMats;

  ////////////////////////////////
  // Collection of 2D Dense matrices, used
  // for weight matrices of MLPs, or 
  // for logistic regression.
  sArray< WeightMatrix* > weightMats;

  ////////////////////////////////
  // Collection of multi-dimensional dense CPTs
  sArray< MDCPT* > mdCpts;

  ///////////////////////////////////
  // Collection of multi-dimensional sparse CPTs (transition matrices, etc.)
  sArray< SMDCPT* > sMdCpts;


  ///////////////////////////////////
  // Collection of multi-dimensional decision-tree based sparse CPTs 
  sArray< DTMDCPT* > dtMdCpts;
  

  //////////////////////////////////////////////////////////////////
  // Basic Gaussian Components
  //////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////
  // Collection of diag. covariance Gaussians
  sArray< DiagGaussian* > diagGaussians;

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
  //////////////////////////////////////////////////////////////////
  // These are only the Possible OBSERVATION DISTRIBUTIONS 
  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////

  ////////////////////////////////
  // Mixtures of Gaussians (could be a heterogeneous mixutre of
  // different types above)
  sArray < GMTK_MixGaussians* > mixGaussians;

  ////////////////////////////////
  // Switching mixtures of Gaussians. The switching is
  // implemented with Gaussians.
  sArray < GMTK_GausSwitchingMixGaussians* > gausSwitchMixGaussians;

  ////////////////////////////////
  // Switching mixtures of Gaussians. The switching is
  // implemented with logistic regression (i.e., 1 layer MLP)
  sArray< GMTK_LogitSwitchingMixGaussians* > logitSwitchMixGaussians;


  ////////////////////////////////
  // Switching mixtures of Gaussians. The switching is
  // implemented with 2 layer (2 weight matrix) MLP
  sArray< GMTK_LogitSwitchingMixGaussians* > mlpSwitchMixGaussians;

  // The following are for discrete observations

  ...

  ...



public:

  ///////////////////////////////////////////////////////////  
  // General constructor
  VariableParams(); 

  ///////////////////////////////////////////////////////////    
  // read in all the basic parameters, assuming file pointer 
  // is located at the correct position.
  void read(iDataStreamFile& is);

  ///////////////////////////////////////////////////////////    
  // write out the basic parameters, starting at the current
  // file position.
  void write(oDataStreamFile& os);

};


#endif // defined GMTK_CPT
